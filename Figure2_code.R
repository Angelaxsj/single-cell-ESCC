library(patchwork)
library(tibble)
library(data.table)
library(stringi)
library(Seurat)
library(ggplot2)
library(ggpubr)
library(ggstatsplot)
library(tidydr)
library(scales)
library(pheatmap)
library(dplyr)
library(ggraph)
library(ggprism)
######==========================Function===========================#####
theme_black <- function(base_size = 12, base_family = "") {
  theme_grey(base_size = base_size, base_family = base_family) %+replace%
    theme(
      axis.line = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      legend.background = element_rect(color = NA, fill = "black"),
      legend.key = element_rect(color = "white",  fill = "black"),
      legend.key.size = unit(1.2, "lines"),
      legend.key.height = NULL,
      legend.key.width = NULL,
      legend.text = element_text(size = base_size*0.8, color = "white"),
      legend.title = element_text(size = base_size*0.8, face = "bold", hjust = 0, color = "white"),
      legend.position = "none",
      legend.text.align = NULL,
      legend.title.align = NULL,
      legend.direction = "vertical",
      legend.box = NULL,
      panel.background = element_rect(fill = "black", color  =  NA),
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing = unit(0, "lines"),
      strip.background = element_rect(fill = "grey30", color = "grey10"),
      strip.text.x = element_text(size = base_size*0.8, color = "white"),
      strip.text.y = element_text(size = base_size*0.8, color = "white",angle = -90),
      plot.background = element_rect(color = "black", fill = "black"),
      plot.title = element_text(size = base_size*1.2, color = "white"),
      plot.margin = unit(rep(0, 4), "lines")
    )
}
get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}
metadata_combine <- function(obj,sampleinfo){
metadata <- obj@meta.data
metadata <- rownames_to_column(metadata,var = "umi")
metadata <- merge(metadata,sampleinfo,by.x="orig.ident",by.y="SampleID") %>% column_to_rownames(var = "umi")
obj@meta.data <- metadata
return(obj)
}

######==========================Figure2A===========================#####
Colors_for_CD8_umap <- c("#E9CFC6","#DAC699","#99568B","#B475A4","#C9A4C5","#9D8EB7","#EFE385","#CFE3C6","#343D74","#545896","#BF6C70")
p1=DimPlot(scRNA,label = F,raster=FALSE,group.by = "celltype",cols = Colors_for_CD8_umap)
# ggsave(p1,file="final_figure/Figures1/ESCC_celltype_umap.tiff",height = 20,width = 20)

p_umap2<-p1+theme_dr(xlength = 0.22, ylength = 0.22,
                     arrow = grid::arrow(length = unit(0.1, "inches"), type = "closed"))+theme(panel.grid = element_blank(),legend.position = "bottom");p_umap2

ggsave(p_umap2,file="final_figure/Figure2/Figure2A.pdf",height = 7,width = 6)


######==========================Figure2B===========================#####
Differentiation <- c("Naive", "Effector", "Exhaustion","Senescence")
Function <- c("TCR_signaling", "Cytotoxicity", "Cytokine",
              "Chemokine", 
              "Costimulatory", "Stress",  "Adhesion",
              "IFNG")
Metabolism <- c("OXPHOS", "Glycolysis", "Lipid_metabolism")
Apoptosis <- c("Pro_apoptosis", "Anti_apoptosis")
MarkerNameVector <- c(Differentiation, Function, Metabolism, Apoptosis)
FunctionScoreMatrix <- matrix(0,
                              ncol = length(unique(sce_CD8T$celltype)),
                              nrow = length(MarkerNameVector))
colnames(FunctionScoreMatrix) <- paste0("CD8_c", 0:10)
rownames(FunctionScoreMatrix) <- MarkerNameVector
for(ci in 1:ncol(FunctionScoreMatrix)){
  for(ri in 1:nrow(FunctionScoreMatrix)){
    FunctionVec <- as_tibble(sce_CD8T@meta.data) %>% pull(MarkerNameVector[ri])
    fv <- mean(FunctionVec[sce_CD8T$celltype == levels(sce_CD8T$celltype)[ci]])
    FunctionScoreMatrix[ri, ci] <- fv
  }
}
FunctionScoreMatrix <- t(apply(FunctionScoreMatrix, 1, rescale, to=c(-1, 1)))
#orderC = c("CD8_c13", "CD8_c3", "CD8_c6", "CD8_c0", "CD8_c11", "CD8_c9", "CD8_c10", "CD8_c12", "CD8_c8", "CD8_c2", "CD8_c7", "CD8_c4", "CD8_c5", "CD8_c1")
#FunctionScoreMatrix <- FunctionScoreMatrix[,orderC]
my.breaks <- c(seq(-1, 0, by=0.01), seq(0.1, 1, by=0.01))
my.colors <- c(
  colorRampPalette(colors = c("#6DCCFD", "white"))(length(my.breaks)/2),
  colorRampPalette(colors = c("white", "#FD9AA0"))(length(my.breaks)/2))

## cellType_col <- data.frame(cell.type = CD8_Obj_CellType)
## rownames(cellType_col) <- colnames(FunctionScoreMatrix)
signatureType_row <- data.frame(Signature.type = c(
  rep("Differentiation", length(Differentiation)),
  rep("Function", length(Function)),
  rep("Metabolism", length(Metabolism)),
  rep("Apoptosis", length(Apoptosis))))
rownames(signatureType_row) <- MarkerNameVector
pdf("final_figure/Figure2/Figure2B.pdf",width = 6.5,height = 7)
pheatmap(FunctionScoreMatrix,
         show_colnames = T,
         show_rownames = T,
         ## annotation_col = cellType_col,
         annotation_row = signatureType_row,
         gaps_row = c(4, 12,15),
         cluster_rows = F,
         cluster_cols = F,
         breaks = my.breaks,
         color = my.colors,
         border_color = "NA",
         fontsize = 8,
         width = 5,
         height = 3.8,
         #filename = file.path(figurePath, paste0("CD8_FunctionScore_heatmap.pdf"))
)
dev.off()





######==========================Figure2C===========================#####
p1=ggplot(df_sub_Trm,aes(group,value,color=group,fill=group))+
  geom_bar(stat="summary",fun=mean,position="dodge")+ #绘制柱状图
  stat_summary(geom = "errorbar",fun.data = 'mean_se', width = 0.3, color="black")+
  labs(x="Samples",y=NULL)+#标题
  theme(panel.grid = element_blank(),
        axis.line = element_line(colour = 'black', size = 1),
        panel.background = element_blank(),
        plot.title = element_text(size = 15, hjust = 0.5),
        plot.subtitle = element_text(size = 15, hjust = 0.5),
        axis.text = element_text(size = 15, color = 'black'),
        axis.title = element_text(size = 15, color = 'black')) +
  scale_color_manual(values = c("#F1D5B3","#DD9542","#B1CBDD"))+ # 自定义边框颜色
  scale_fill_manual(values = c("#F1D5B3","#DD9542","#B1CBDD"))+
  #
  stat_compare_means(comparisons = list(c("Pre","Post_R"),
                                        c("Pre","Post_NR"),
                                        c("Post_R","Post_NR")
                                        
  ),
  # 设置需要比较的组
  map_signif_level = F, #是否使用星号显示
  test = "wilcox.test", ##计算方法
  y_position = c(50,55,60),#图中横线位置设置
  tip_length = c(c(0.01,0.01),
                 c(0.01,0.01),
                 c(0.01,0.01)),#横线下方的竖线设置
  size=6,color="black")
p1

p2=ggplot(df_sub_Tex,aes(group,value,color=group,fill=group))+
  geom_bar(stat="summary",fun=mean,position="dodge")+ #绘制柱状图
  stat_summary(geom = "errorbar",fun.data = 'mean_se', width = 0.3, color="black")+
  labs(x="Samples",y=NULL)+#标题
  theme(panel.grid = element_blank(),
        axis.line = element_line(colour = 'black', size = 1),
        panel.background = element_blank(),
        plot.title = element_text(size = 15, hjust = 0.5),
        plot.subtitle = element_text(size = 15, hjust = 0.5),
        axis.text = element_text(size = 15, color = 'black'),
        axis.title = element_text(size = 15, color = 'black')) +
  scale_color_manual(values = c("#F1D5B3","#DD9542","#B1CBDD"))+ # 自定义边框颜色
  scale_fill_manual(values = c("#F1D5B3","#DD9542","#B1CBDD"))+
  #
  stat_compare_means(comparisons = list(c("Pre","Post_R"),
                                        c("Pre","Post_NR"),
                                        c("Post_R","Post_NR")
                                        
  ),
  # 设置需要比较的组
  map_signif_level = F, #是否使用星号显示
  test = "wilcox.test", ##计算方法
  y_position = c(50,55,60),#图中横线位置设置
  tip_length = c(c(0.01,0.01),
                 c(0.01,0.01),
                 c(0.01,0.01)),#横线下方的竖线设置
  size=6,color="black")
p2
p=p1+p2
ggsave(p,file="final_figure/Figure2/Figure2C.pdf",height = 5,width = 5)


######==========================Figure2D&E===========================#####
library(ggplot2) ##绘图使用
library(ggprism) ##设置主题私用
library(ggrepel) ##给火山图加标签使用

p1=ggplot(Timepoint_deg_comb, aes(x =avg_log2FC, y=-log10(p_val_adj), colour=DEG)) +
  geom_point(alpha=0.6, size=5) +  #点的透明度和大小
  scale_color_manual(values=c('steelblue','brown')) + #调整点的颜色
  xlim(c(-3, 3)) +  ##调整x轴的取值范围，max(abs(BRCA_Match_DEG$logFC))最大值是多少，再进行取舍
  geom_vline(xintercept=c(0,0),lty=2,col="black",lwd=0.8) + #添加x轴辅助线,lty函数调整线的类型
  geom_hline(yintercept = -log10(0.05), lty=4,col="black",lwd=0.8) +  #添加y轴辅助线
  labs(x="log2FC", y="-log10FDR") +  #x、y轴标签
  ggtitle("Tem DEG within Chemo+Anti-PD1") + #标题
  #theme(plot.title = element_text(hjust = 0.5),legend.position="right",legend.title = element_blank())+
  geom_text_repel(data = up,
                  aes(x = avg_log2FC, y = -log10(p_val_adj), label = Gene)) +theme_prism(border = T) +
  geom_text_repel(data = down,
                  aes(x = avg_log2FC, y = -log10(p_val_adj), label = Gene))
#调整主题
p1


ggsave(p1,file="final_figure/Figure2/Figure2D.pdf",width = 7.5,height = 5)

p2=ggplot(Timepoint_deg_mono, aes(x =avg_log2FC, y=-log10(p_val_adj), colour=DEG)) +
  geom_point(alpha=0.6, size=5) +  #点的透明度和大小
  scale_color_manual(values=c('steelblue','brown')) + #调整点的颜色
  xlim(c(-3, 3)) +  ##调整x轴的取值范围，max(abs(BRCA_Match_DEG$logFC))最大值是多少，再进行取舍
  geom_vline(xintercept=c(0,0),lty=2,col="black",lwd=0.8) + #添加x轴辅助线,lty函数调整线的类型
  geom_hline(yintercept = -log10(0.05), lty=4,col="black",lwd=0.8) +  #添加y轴辅助线
  labs(x="log2FC", y="-log10FDR") +  #x、y轴标签
  ggtitle("Tem DEG within Chemo+Anti-PD1") + #标题
  #theme(plot.title = element_text(hjust = 0.5),legend.position="right",legend.title = element_blank())+
  geom_text_repel(data = up,
                  aes(x = avg_log2FC, y = -log10(p_val_adj), label = Gene)) +theme_prism(border = T) +
  geom_text_repel(data = down,
                  aes(x = avg_log2FC, y = -log10(p_val_adj), label = Gene))
#调整主题
p2


ggsave(p2,file="final_figure/Figure2/Figure2E.pdf",width = 7.5,height = 5)



######==========================Figure2G===========================#####
ibrary(ComplexHeatmap)
library(Startrac)
TCR_ex <- subT.meta[subT.meta$Clone_NUM>1,]
names(TCR_ex)[3] <- "clone.id"
TCR_ex_post <- TCR_ex[TCR_ex$Timepoint=="Treatment",]
tmp <- dplyr::select(TCR_ex_post, c("celltype", "orig.ident")) 
out3 <- Startrac.run(TCR_ex_post,verbose=F)
dat.plot <- as.matrix(subset(out3@pIndex.tran,aid==out3@proj)[,c(-1,-2,-3)])
rownames(dat.plot) <- subset(out3@pIndex.tran,aid==out3@proj)[,3]
dat.plot[is.na(dat.plot)] <- 0
yrange <- pretty(dat.plot)
col.heat <- colorRamp2(seq(0,max(yrange),length=15),
                       colorRampPalette(rev(brewer.pal(n=7,name="RdBu")))(15),
                       space = "LAB")
col.heat<- colorRampPalette(c( "white","#B2D0E7","#3E7BB7","#E7B731","#F2B129","#C84529"))(200)

pdf("final_figure/Figure2/Figure2G.pdf",width = 6,height = 5)
Heatmap(dat.plot,name="pIndex.tran",col = col.heat,border = T)+
  labs(title = "CD8_TCR_clonetype_sharing_across_clusters")
dev.off()
######==========================Figure2H===========================#####
sampleinfo <- read.table(file = "~/Project/scESCC/ESCC_metadata_20240705.txt",header = T,sep = "\t")
CD8T <- scRNA
CD4T <- readRDS("~/Project/scESCC/CellType/Tcell/Rdata/CD4T.annotation.rds")

library(ggpubr)
metadata <- CD8T@meta.data
tmp <- dplyr::select(metadata, c("orig.ident", "celltype"))
df <- data.frame()
for(i in unique(tmp$orig.ident)){
  df_i <- subset(tmp, tmp$orig.ident==i) %>% pull(celltype) %>% table()
  df_i <- data.frame(sample=rep(i, length(df_i)), value=df_i/nrow(metadata[which(metadata$orig.ident==i),])*100)
  names(df_i) <- c("SampleID", "celltype", "value")
  df <- rbind(df, df_i)
}  

df <- merge(df,sampleinfo,by = "SampleID") 
df_CD8T <- df
df_CD8T_cor <- reshape2::dcast(df_CD8T,SampleID~celltype)


metadata <- CD4T@meta.data
tmp <- dplyr::select(metadata, c("orig.ident", "celltype"))
df <- data.frame()
for(i in unique(tmp$orig.ident)){
  df_i <- subset(tmp, tmp$orig.ident==i) %>% pull(celltype) %>% table()
  df_i <- data.frame(sample=rep(i, length(df_i)), value=df_i/nrow(metadata[which(metadata$orig.ident==i),])*100)
  names(df_i) <- c("SampleID", "celltype", "value")
  df <- rbind(df, df_i)
}  

df <- merge(df,sampleinfo,by = "SampleID") 
df_CD4T <- df
df_CD4T_cor <- reshape2::dcast(df_CD4T,SampleID~celltype)

df_cor <- merge(df_CD4T_cor,df_CD8T_cor,by="SampleID")
cor.test(df_cor$`CD4_c2_Treg_TNFRSF9+`,df_cor$CD8_c4_Tex_ITGAE)
pdf(file = "final_figure/Figure2/Figure2G.pdf",height = 5,width = 3)
ggscatter(df_cor, x = "CD4_c2_Treg_TNFRSF9", y = "CD8_c4_Tex_ITGAE",
          add = "reg.line", conf.int = TRUE,size=2,
          add.params = list(fill = "#6DCCFD")
          
)+
  stat_cor(method = "spearman"
  )
dev.off()
######==========================Figure2H  ===========================#####
ds <- list(PD1 = density(PT_1[names(PT_1) %in% PD1_Cells]),
           Chemo_PD1 = density(PT_1[names(PT_1) %in% Chemo_PD1_Cells]))


# Scale the pseudotime from 0 to 100
ds$PD1$x <- rescale(ds$PD1$x, to = c(0, 100))
ds$Chemo_PD1$x <- rescale(ds$Chemo_PD1$x, to = c(0, 100))

xlim <- range(c(ds$PD1$x, ds$Chemo_PD1$x))
ylim <- range(c(ds$PD1$y, ds$Chemo_PD1$y))

pdf(file = "final_figure/Figure2/Figure2H_lineage1_combtreat_vs_monotreat.pdf",width = 4.5,height = 5)
plot(xlim, ylim, col = "white", xlab = "Pseudotime", ylab = "", cex.lab=0.5, cex.axis=0.5)

polygon(c(min(ds$PD1$x), ds$PD1$x, max(ds$PD1$x)), c(0, ds$PD1$y, 0),
        col = alpha(Density_colors[1], alpha = .5))

polygon(c(min(ds$Chemo_PD1$x), ds$PD1$x, max(ds$PD1$x)), c(0, ds$Chemo_PD1$y, 0), # 2 midlle ones wrong?
        col = alpha(Density_colors[2], alpha = .5))
legend("topleft", legend = c("anti-PD1", "TP+anti-PD1"), 
       fill = alpha(Density_colors, alpha = .5), bty = "n", cex = 0.5, x.intersp = 1, y.intersp = 1.5)
title("Lineage1", cex.main = 0.75, adj=0, line = 0.25)
dev.off()

######==========================Figure2I  ===========================#####
subT.meta <- readRDS(file = "final_figure/Figure2/subT.meta.rds")
TCR_ex <- subT.meta[subT.meta$clone.status=="Clonal",]
#sc_CD8TCR <- readRDS(file = "final_figure/Figure2/sc_CD8TCR.rds")

patientID_comb <- c("P14","P17","P19","P20","P21","P22","P24")
TCR_ex_paired <- TCR_ex[TCR_ex$PatientID%in%patientID_comb,]
Pre <- TCR_ex_paired %>% filter( Timepoint== "Baseline")
Post <- TCR_ex_paired %>% filter(Timepoint == "Treatment")

#Pre
ann <- Pre %>% group_by(PatientID, stype, celltype) %>% 
  summarize(nr_Pre = n_distinct(Clone_AA)) 
tot_pre <- Pre %>%  group_by(PatientID, stype) %>% summarize(tot_pre = n_distinct(Clone_AA)) %>%
  dplyr::rename(stype_pre = stype)
nr_TCR_pre <- ann %>% left_join(tot_pre) %>% dplyr::rename(Pre = celltype)

#Post
ann <- Post %>% group_by(PatientID, stype, celltype) %>% 
  summarize(nr_Post = n_distinct(Clone_AA)) 
tot_post <- Post %>%  group_by(PatientID, stype) %>% summarize(tot_post = n_distinct(Clone_AA)) %>%
  dplyr::rename(stype_post = stype)
nr_TCR_post <- ann %>% left_join(tot_post) %>% dplyr::rename(Post = celltype)


Tcellsharing3$combination <- paste(Tcellsharing3$Pre, Tcellsharing3$Post, sep = "_")

tmp <- Tcellsharing3 %>% full_join(tot_pre) %>% full_join(tot_post) %>% ungroup
patient_combination <- tmp %>% tidyr::expand(combination, PatientID) %>% filter(!is.na(combination))
combination_info <- Tcellsharing3 %>% dplyr::select(Pre, Post, stype_Pre, stype_Post, combination) %>% distinct()     
patient_combination <- patient_combination %>% left_join(combination_info)

prop <- Tcellsharing3 %>% full_join(patient_combination) %>% 
  tidyr::replace_na(list(Freq_Pre=0, Freq_Post=0)) %>% 
  left_join(tot_pre) %>% left_join(tot_post) %>% 
  mutate(Prop_pre = Freq_Pre / tot_pre, Prop_post = Freq_Post / tot_post)      

meta <- subT.meta%>% dplyr::select(PatientID, Treatment) %>% distinct()
prop <- left_join(prop,meta )
#write.csv(prop, paste0(outdir, "TCRs_shared_", name, "_proportion_corr_CD4_CD8_Tcells.csv"), row.names = F)


# Calculate average per group to plot
prop_group <- prop %>% dplyr::group_by(combination, Treatment) %>% 
  mutate(mean_prop_pre = mean(Prop_pre, na.rm = T), mean_prop_post = mean(Prop_post, na.rm = T)) %>% ungroup() %>% 
  dplyr::select(!c(PatientID,  Freq_Pre, Freq_Post, tot_pre, tot_post, Prop_pre, Prop_post)) %>% distinct() %>% arrange(combination)
#write.csv(prop_group, paste0(outdir, "TCRs_shared_", name, "_prop_group_corr_CD4_CD8_Tcells.csv"), row.names = F)

# Long format
long <- to_lodes_form(prop_group, key = "loc", axes = c("Pre", "Post"))
df <- long %>% dplyr::rename(Pre = mean_prop_pre, Post = mean_prop_post) %>% 
  tidyr::pivot_longer(cols = c(Pre, Post)) %>%  filter(loc == name)%>% dplyr::rename(mean_prop = value) %>% dplyr::select(!name)
#write.csv(df, paste0(outdir, "TCRs_shared_", name, "_prop_group_corr_CD4_CD8_Tcells_long_alluvial_format.csv"), row.names = F)


### plot per group


df$loc <- factor(df$loc, levels = c("Pre", "Post"))


pdf(paste0(outdir, name, "_alluvial_relative_values.pdf"), height = 6, width = 12)

df <- df %>% mutate(mean_prop_0 = tidyr::replace_na(value, 0))


plot <- ggplot(data = df, aes(x = loc, y = mean_prop, stratum = stratum, alluvium = alluvium, fill = stratum, label = stratum)) +
  geom_flow(stat = "flow", knot.pos = 1/4, aes.flow = "forward") + 
  geom_stratum() +
  scale_fill_manual(values = Colors_for_CD8_umap) + 
  xlab("") +
  ylab("Proportion shared TCRs") +
  #ggtitle(name) +
  theme(text=element_text(size=23),
        legend.title=element_blank(),
        axis.text.x = element_text(size=23, colour = "black"),
        axis.text = element_text(colour = "black"),
        strip.text = element_text(size = 26, colour = "black"), 
        axis.ticks.x = element_blank(), 
        strip.background = element_blank(), 
        panel.border = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black")) +
  scale_y_continuous(expand = expansion(mult = c(0, .1)))
plot






