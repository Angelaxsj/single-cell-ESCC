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
