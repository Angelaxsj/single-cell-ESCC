## NMF program
library(tidyr)
library(data.table)
library(ComplexHeatmap)
library(dplyr)
library(ggplot2)
require(ggpubr)
library(Seurat)
library(NMF)
library(ggradar)
library(corrr)
library(circlize)
library(tibble)
library(GSVA)

######==========================Figure5B===========================#####
scESCC <- readRDS("~/Project/scESCC/CellType/Maincell/Rdata/scESCC_annotation.rds")
sampleinfo <- read.table(file = "~/Project/scESCC/ESCC_metadata_20240705.txt",header = T,sep = "\t")
metadata_scESCC <- scESCC@meta.data
tmp <- dplyr::select(metadata_scESCC, c("orig.ident", "celltype"))
tmp$celltype <- as.character(tmp$celltype)
#tmp <- tmp[!tmp$celltype%in%exclded_celltype,]

df <- data.frame()
for(i in unique(tmp$orig.ident)){
  df_i <- subset(tmp, tmp$orig.ident==i) %>% pull(celltype) %>% table()
  df_i <- data.frame(sample=rep(i, length(df_i)), value=df_i/nrow(metadata_scESCC[which(metadata_scESCC$orig.ident==i),])*100)
  names(df_i) <- c("SampleID", "celltype", "value")
  df <- rbind(df, df_i)
}

df <- merge(df,sampleinfo,by = "SampleID")
# df_baseline <- df[df$Timepoint=="Baseline",]
# df_ESCC <- df_baseline[,c(1:3)]
exclded_celltype <- c("MP1","MP2","MP3","MP4","MP5","MP6")
df_ESCC <- df[!df$celltype%in%exclded_celltype,]
df_heatmap <- reshape2::dcast(df_ESCC,SampleID~celltype)
df_heatmap[is.na(df_heatmap)] <- 0
df_heatmap <- column_to_rownames(df_heatmap,var = "SampleID") %>% t( ) %>% as.data.frame()





# NMF analysis

mtx = t(df_heatmap)
nmf.res = nmf(t(mtx), rank=3:6, nrun=20, seed=123456)
plot(nmf.res)
nmf.res = nmf(t(mtx), rank=5, nrun=20, seed=123456) 
w = basis(nmf.res)
rownames(w) = colnames(mtx)
colnames(w) = paste0('NMF', seq(1,5))

w_ = w %>% t %>% scale
col = colorRamp2(breaks = c(-2, 0, 2), colors = c('#798ED3', 'white', '#EFA3A3'))
p=ComplexHeatmap::Heatmap(w_, width = 11, height = 4, column_km = 5, 
                          cluster_rows = F, clustering_method_columns = 'ward.D',col = col) 

heatmap = draw(p)





######==========================Figure5C===========================#####
nrank=5
h = coef(nmf.res) %>% t
colnames(h) = paste0('NMF', seq(1,nrank))

h_df = merge(h,sampleinfo,by.x=0,by.y="SampleID") 
names(h_df)[1] <- "SampleID"
patientID_comb <- c("P14","P17","P19","P20","P21","P22","P24")
h_df$Timepoint[h_df$Timepoint=="Baseline"] <- "Pre"
h_df$Timepoint[h_df$Timepoint=="Treatment"] <- "Post"
h_df_comb <- h_df[h_df$PatientID%in%patientID_comb,]
input = h_df_comb %>% 
  pivot_longer(names_to = 'NMF', values_to = 'value', cols = contains('NMF')) %>%
  pivot_wider(names_from = 'Timepoint', values_from = 'value',id_cols =c( 'PatientID','NMF','Response','Treatment') )%>%
  mutate(dnmc=Post-Pre)
pdf(file = "final_figure/Figure5/Figure5C.pdf",width = 8,height = 3)
input$Response <- factor(input$Response ,levels = c("responder","non-responder"))
input %>%
  ggplot(aes(x=Response, y=dnmc)) +
  geom_boxplot(aes(col=Response,fill=Response),alpha=.5, width=.6, lwd=.8, outlier.color = NA) +
  geom_point(aes(col=Response), position = position_jitter(width = .15)) +
  geom_hline(yintercept = 0, linetype='dashed') +
  facet_wrap(~NMF, nrow=1) +
  scale_color_manual(values = c('#DD9442','#3C7DAA','#366E82'), name='Response') +
  scale_fill_manual(values = c('#DD9442','#3C7DAA','#366E82'), name='Response')+
  theme_pubr() + ylab('NMF score dynamic') + xlab('') +
  theme(plot.title = element_text(face = "bold", size=16),
        strip.text = element_text(size = 15, color='black'),
        strip.background =element_rect(fill=NA, color = NA, size = .5), 
        axis.text.x = element_text(size=15, angle = 0, hjust = .5, vjust=.5),
        axis.text.y = element_text(size=15),
        axis.title = element_text(size = 18),
        axis.line = element_line(size = .5),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.position = 'right',
        legend.key.size = unit(1.4,'line'),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 15, face = 'bold'),
        plot.margin = margin(t = 0, l=0))
  dev.off()


######==========================单细胞水平相关性===========================#####
sc_Tumor_clean <- readRDS("~/Project/scESCC/CellType/Epi/Rdata/sc_Tumor_clean.rds")
metadata <- scESCC@meta.data
tmp <- dplyr::select(metadata, c("orig.ident", "celltype"))
df <- data.frame()
for(i in unique(tmp$orig.ident)){
  df_i <- subset(tmp, tmp$orig.ident==i) %>% pull(celltype) %>% table()
  df_i <- data.frame(sample=rep(i, length(df_i)), value=df_i/nrow(tmp[which(tmp$orig.ident==i),])*100)
  names(df_i) <- c("SampleID", "celltype", "value")
  df <- rbind(df, df_i)
}  
df_MP <- reshape2::dcast(df,SampleID~celltype)


df_cor <- merge(df_MP,h_df,by="SampleID")


cor.test(df_cor$MP6,df_cor$NMF4,method = "spearman")
df_cor_sub <- df_cor[,colnames(df_cor)%in%c("MP6","NMF1","NMF2","NMF3","NMF4","NMF5")]
pdf(file = "final_figure/Figure5/Figrue5e_correlation_MP6_EC3.pdf",width = 4,height = 5)
df_cor_MP6 <- reshape2::melt(df_cor_sub,id.vars = "MP6",measure.vars = c("NMF1","NMF2","NMF3","NMF4","NMF5"))

pdf(file="final_figure/Figure5/Figure5_EC_MP6_correlation.pdf",width =4.5 ,height = 5)
ggplot(df_cor_MP6, aes(x = MP6, y = value)) +
 # geom_rug(aes(color =variable)) +
  geom_smooth(aes(color = variable), method = lm, 
              se = FALSE, fullrange = F)+
  scale_color_manual(values = c("#8DB2C6", "#97D1B3", "#CC9999","#FFCC99","#F49FD5"))+
  #ylim(0,0.06)+
  theme_classic()+
  ggpubr::stat_cor(aes(color = variable),method = "spearman")

dev.off()






######==========================Nichenet分析具体的调控机制===========================#####
###########nichenet 分析###############
sc_Tumor_clean <- readRDS("~/Project/scESCC/CellType/Epi/Rdata/sc_Tumor_clean.rds")
sc_Tumor_clean <- metadata_combine(sc_Tumor_clean,sampleinfo)
sc_TME4 <- subset(sc_TME,subset=celltype%in%EC4_celltype)
sc_TME4 <- metadata_combine(sc_TME4,sampleinfo)
sc_TME4_sub <- subset(sc_TME4,subset=Response=="non-responder"&Treatment=="Chemo+PD1")

Idents(sc_TME4_sub) <- "Timepoint"

MP6 <- subset(sc_Tumor_clean,subset=MP_module=="MP6"&Response=="non-responder"&Treatment=="Chemo+PD1")

sc_TME4_sub@meta.data <- sc_TME4_sub@meta.data[,c(1,14:21)]
MP6@meta.data <- MP6@meta.data [,c(1,17,24:30)]
names(sc_TME4_sub@meta.data)[2] <- "stype"
sc_TME4_sub$stype <-
names(MP6@meta.data)[2] <- "stype"
sce <- merge(sc_TME4_sub,MP6)
Idents(sce) <- "stype"
# sce <- NormalizeData(sce) %>% FindVariableFeatures() %>% ScaleData()
# sce <- RunPCA(sce,npcs = 50,verbose = F)
# sce <- RunHarmony(sce,max.iter.harmony = 15,group.by.vars = "orig.ident")
# metadata <- scRNA_singlet@meta.data
# colnames(metadata)
# Idents(scRNA_singlet_clean) <- "annotation"
weighted_networks <- readRDS("/home/xiangsj/Project/scRNA_ES/ESCC_Rproject/weighted_networks.rds")
ligand_target_matrix <- readRDS("/home/xiangsj/Project/scRNA_ES/ESCC_Rproject/ligand_target_matrix.rds")
lr_network <- readRDS("/home/xiangsj/Project/scRNA_ES/ESCC_Rproject/lr_network.rds")

nichenet_output = nichenet_seuratobj_aggregate(seurat_obj = sce, 
                                               top_n_ligands = 20,
                                               receiver = c("EC4"),  ####受体细胞
                                               sender  = c("MP6"),  ####配体细胞
                                               condition_colname = "Timepoint", ####分组
                                               condition_oi = "Treatment", ###实验组
                                               condition_reference = "Baseline", ####对照组
                                               ligand_target_matrix = ligand_target_matrix, 
                                               lr_network = lr_network, 
                                               weighted_networks = weighted_networks, 
                                               organism = "human",
                                               assay_oi = "RNA") 
nichenet_output$ligand_activities####首先查看最高活性的配体是哪些
nichenet_output$top_ligands

#######热图可视化####
nichenet_output$ligand_differential_expression_heatmap

pdf("final_figure/Figure5/MP6_EC4_interaction_in_nonresponder.pdf",width = 9,height = 5)
nichenet_output$ligand_target_heatmap + 
  scale_fill_gradient2(low = "whitesmoke",  high = "royalblue")
dev.off()




######==========================在bulk队列中进行进行验证===========================#####
clustermarkers_noripo <- read.csv(file = "final_figure/Figure5/subcell_marker_genes_noripo.csv",sep=",")
top_marker <- clustermarkers_noripo[clustermarkers_noripo$avg_log2FC>0.8,]
top_marker <-  top_n(group_by(top_marker,cluster), n = 100, wt = avg_log2FC)
#View(top_marker)
EC1_signature <- top_marker[top_marker$cluster=="EC1",]

EC2_signature <- top_marker[top_marker$cluster=="EC2",]

EC3_signature <- top_marker[top_marker$cluster=="EC3",]

EC4_signature <- top_marker[top_marker$cluster=="EC4",]

EC5_signature <- top_marker[top_marker$cluster=="EC5",]

geneset <- list("EC1"=EC1_signature$gene,"EC2"=EC2_signature$gene,"EC3"=EC3_signature$gene,"EC4"=EC4_signature$gene,"EC5"=EC5_signature$gene)



NICE1_study <- readRDS("Bulk/Cancer_Cell_IOC_NICE1.RDS")
SCH_Clinical_data <- NICE1_study$pheno
SCH_expr <- NICE1_study$TPM$count.combat.symbol
gsva_data <- gsva(as.matrix(SCH_expr),geneset,method = "zscore",kcdf="Gaussian")
gsva_data[1:5,1:5]
gsva_data <- as.data.frame(t(gsva_data))
gsva_data <- apply(gsva_data,2,normalize) %>% as.data.frame()
df <- merge(SCH_Clinical_data,gsva_data,by.x="Sample",by.y=0)
ggbetweenstats(df,x="ResponseGroup2",y="EC4")+stat_compare_means()
df$feature <- df$EC3-df$EC4
df$Response <- ifelse(df$ResponseGroup2=="poor-no response","NR","R")



###======================生存分析======================####
library(survival)
library(survminer)
res.cut=surv_cutpoint(data =df,time = "Diseasefreesurvival_time",event = "Diseasefreesurvival_status",
                      variables = "EC4" ,minprop = 0.1)#注意调整minprop这个参数
summary(res.cut)
res.cut=surv_categorize(res.cut)
head(res.cut)
fit=survfit(Surv(Diseasefreesurvival_time,Diseasefreesurvival_status)~EC4,data = res.cut)

pdf("final_figure/Figure5/NICE1_study_EC3_divided_EC4_score_DFS.pdf",width = 3.5 ,height=5,onefile = FALSE)
ggsurvplot(fit = fit,data = res.cut,
           #legend.titile="TCR_Ex",
           #legend.labs=c("TCR_Ex_high","TCR_Ex_low"),
           conf.int = F,##添加置信区间
           surv.median.line = "hv",
           pval = T,
           pval.method = T,
           risk.table = T,
           tables.height = 0.2,
           tables.theme = theme_cleantable(),
           palette =  c("#FFCC99","#B6B5CC"),
           size=2,
           ggtheme = theme_bw())+labs(title = "EC3/EC4")

dev.off()


######==========================在验证队列的单细胞水平进行验证===========================#####
library(AUCell)
sc_ESCC_cc2023 <- readRDS("~/Project/scESCC/CellType/Maincell/Rdata/sc_ESCC_cc2023.rds")
sampleinfo_2023 <- read.table(file="ESCC_validation_sccohort/Sampleinfo_validation.txt",header = T,sep = "\t")
sc_ESCC_cc2023_TME <- subset(sc_ESCC_cc2023,subset=main.celltype=="Epithelial",invert=T)
cells_rankings_validation <- AUCell_buildRankings(sc_ESCC_cc2023_TME@assays$RNA@counts,nCores=10, plotStats=F)##基因排序
cells_AUC <- AUCell_calcAUC(geneset, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05) ##计算AUC值。
#cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=1,assign=TRUE)##挑选阈值
aucs <- getAUC(cells_AUC)
aucs <- as.data.frame(t(aucs))
sce <- AddMetaData(sc_Tumor_clean,aucs)


value_df <- h_df[,c(2:6,15)]
data <- data.frame()
for(i in 1:5){
  data_i<- apply(value_df, 1, function(x){ x[i]/x[6]} )
  data_i <- as.data.frame(t(data_i))
 # names(data_i) <- paste0("NMF",i)
 data <- rbind(data,data_i)

}

rownames(data) <- c("NMF1","NMF2","NMF3","NMF4","NMF5")
colnames(data) <- h_df$SampleID
data <- as.data.frame(t(data))%>% merge(sampleinfo,by.x=0,by.y="SampleID")
data_combpaired <- data[data$PatientID%in%patientID_comb,]

# data_base <- data[data$Timepoint=="Treatment"&data$Treatment=="Chemo+PD1"&data$Response%in%c("responder","non-responder"),]
# ggbetweenstats(data_base,x="Response",y="NMF3")+stat_compare_means()




data$Timepoint[data$Timepoint=="Baseline"] <- "Pre"
data$Timepoint[data$Timepoint=="Treatment"] <- "Post"
data <- data[data$Treatment=="Chemo+PD1",]
input = data %>% 
  pivot_longer(names_to = 'NMF', values_to = 'value', cols = contains('NMF')) %>%
  pivot_wider(names_from = 'Timepoint', values_from = 'value',id_cols =c( 'PatientID','NMF','Response','Treatment') )%>%
  mutate(dnmc=Post-Pre)
input$group <- "Monotreat"
input$group[input$Treatment=="Chemo+PD1"&input$Response=="responder"] <- "Combresponder"
input$group[input$Treatment=="Chemo+PD1"&input$Response=="non-responder"] <- "Combnonresponder"
input$group <- factor(input$group,levels = c("Combresponder","Combnonresponder","Monotreat"))
pdf(file = "final_figure/Figure5/Figure5B_dynamic_boxplot.pdf",width = 10,height = 3)
input %>%
  ggplot(aes(x=Response, y=dnmc)) +
  geom_boxplot(aes(col=Response,fill=Response),alpha=.5, width=.6, lwd=.8, outlier.color = NA) +
  geom_point(aes(col=Response), position = position_jitter(width = .15)) +
  geom_hline(yintercept = 0, linetype='dashed') +
  facet_wrap(~NMF, nrow=1) +
  scale_color_manual(values = c('#DD9442','#3C7DAA',), name='Response') +
  scale_fill_manual(values = c('#DD9442','#3C7DAA',), name='Response')+
  theme_pubr() + ylab('NMF score dynamic') + xlab('') +
  theme(plot.title = element_text(face = "bold", size=16),
        strip.text = element_text(size = 15, color='black'),
        strip.background =element_rect(fill=NA, color = NA, size = .5), 
        axis.text.x = element_text(size=15, angle = 0, hjust = .5, vjust=.5),
        axis.text.y = element_text(size=15),
        axis.title = element_text(size = 18),
        axis.line = element_line(size = .5),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.position = 'right',
        legend.key.size = unit(1.4,'line'),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 15, face = 'bold'),
        plot.margin = margin(t = 0, l=0))+stat_compare_means(method = "wilcox.test",comparisons = list(c("responder","non-responder")
                                                                                                      ))
dev.off()

input_nmf_NMF3 <- input[input$NMF=="NMF3",]




