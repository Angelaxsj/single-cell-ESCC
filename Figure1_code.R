library(patchwork)
library(tibble)
library(data.table)
library(stringi)
library(Seurat)
library(ggplot2)
library(ggpubr)
library(ggstatsplot)
library(tidydr)
library(Nebulosa)
library(dplyr)
######==========================Figure1C===========================#####

Colors_for_umap <- c("#A0C974","#AD98C3","#2F73A2","#CFA5C8","#A7C1CC","#A25836","#6C6D9E","#B95094","#CEDBA8","#EB8D8F",
                     "#95BED1","#F0D1C8","#E2C99A","#C86F74","#A65795","#C379AF","#D3A8CC","#A392BE","#F7EA86","#D2E8CA",
                     "#383F7A","#595C9F","#F3B576","#EE8975","#C0D975","#C3C0D9","#9FD3CA","#FAF7BF","#A1BED1","#33A0D1",
                     "#A6A6D2","#70B5A3","#955187","#D57A2C","#4B7393","#3C83B4",'#E59840',"#366E82","#F89EC0","#853E0A",
                     "#C6AFB1","#9ED0C8","#BD9A7D","#D42A28","#4B9044","#64BCAA","#58944F","#B4264B","#ED864F","#D14F40",
                     "#CBC7E1","#9D3761","#2D3465","#C6AFB1","#996FA8","#E7E4A1","#74A1BD","#918138","#82AB8E","#E35530",
                     "#707F45","#D1D9A6","#93743B","#E4C069","#E6D1A2","#844640","#CB6C73","#DD9DA0","#C89FC0","#88C383",
                     "#5EB29B","#EB8365","#CC7EAF","#D37E30","#93C565","#EBC44E")
p1=DimPlot(scESCC_annotation,label = F,raster=FALSE,group.by = "celltype",cols = Colors_for_umap)
# ggsave(p1,file="final_figure/Figures1/ESCC_celltype_umap.tiff",height = 20,width = 20)

p_umap2<-p1+theme_dr(xlength = 0.22, ylength = 0.22,
                     arrow = grid::arrow(length = unit(0.3, "inches"), type = "closed"))+theme(panel.grid = element_blank())+NoLegend();p_umap2

ggsave(p_umap2,file="final_figure/Figure1/ESCC_celltype_umap.tiff",height = 10,width = 10)


######==========================Figure1D===========================#####
scESCC_annotation$maincell <- scESCC_annotation$main.celltype
Colors_for_mainumap <- c("#8DC284","#DCA3A3","#E4856A","#487CA8","#C4C0D8","#EE9AB9","#E0C16F","#6E6D9D")

metadata <- scESCC_annotation@meta.data
tmp <- dplyr::select(metadata, c("orig.ident", "maincell"))

df <- data.frame()
for(i in unique(tmp$orig.ident)){
  df_i <- subset(tmp, tmp$orig.ident==i) %>% pull(maincell) %>% table()
  df_i <- data.frame(sample=rep(i, length(df_i)), value=df_i/nrow(metadata[which(metadata$orig.ident==i),])*100)
  names(df_i) <- c("SampleID", "celltype", "value")
  df <- rbind(df, df_i)
}  
df <- merge(df,sampleinfo,by = "SampleID")

df <- df[order(df$Timepoint,df$Treatment),]
sampleID <-unique(df$SampleID)
df$SampleID <- factor(df$SampleID,levels = sampleID)
p <- ggplot(df, aes(x=SampleID, y=value, fill=celltype)) + 
  geom_bar(stat= "identity", position = "fill")  +
  labs(x = 'Samples', y = 'Relative Abundance', title = 'Samples composition') +
  scale_fill_manual(values = Colors_for_mainumap)+
  theme(panel.grid = element_blank(), panel.background = element_blank(), 
        legend.position = "bottom",
        axis.text.x = element_text(angle=45))+coord_flip()
p  
ggsave(p,file="final_figure/Figure1/Figure1D.pdf",height = 8,width=6)

######==========================Figure1E===========================#####
p1=plot_density(scESCC_annotation,features = c("CD3D","EPCAM","DCN","CD68"),pal = "plasma")
p2=plot_density(scESCC_annotation,features = c("CD19","IGHG1","KIT","PECAM1"),pal = "plasma")
p=p1|p2
ggsave(p,file="final_figure/Figure1/Figure1F.pdf",height = 6,width=3)
p
dev.off()

######==========================Figure1F===========================#####
df <- merge(df,sampleinfo,by = "SampleID") 
patientID_mono <- c("P2","P3","P4","P12")
patientD_comb <- c("P14","P17","P19","P20","P21","P22","P24")
df_sub <- df[df$PatientID%in%patientD_comb,]
p=ggplot(df_sub, aes(x = Timepoint, y = value)) +
  geom_boxplot(aes(fill = Timepoint), show.legend = F, width = 0.6,alpha=0.8) +  #箱线图
  geom_point(aes(fill = Timepoint),size = 3,shape=21) +  #绘制散点
  scale_fill_manual(values = c("#9F89B2","#564478")) +  #设置颜色
  geom_line(aes(group = PatientID), color = '#BBB4C9', lwd = 0.5) +  #配对样本间连线
  theme(panel.grid = element_blank(),
        axis.line = element_line(colour = 'black', size = 1),
        panel.background = element_blank(),
        plot.title = element_text(size = 15, hjust = 0.5),
        plot.subtitle = element_text(size = 15, hjust = 0.5),
        axis.text = element_text(size = 15, color = 'black'),
        axis.title = element_text(size = 15, color = 'black')) +
  labs(x = '', y = 'gene.expression', title = 'Baseline vs Treatment')+
  facet_wrap(.~celltype,                     # 同时做几个图，这几个图分别对应数据转换前的第1列到第4列
             scale="free")+
  stat_compare_means(method = "wilcox.test",paired = TRUE,
                     comparisons=list(c("Baseline", "Treatment")))

p  
ggsave(p,file="final_figure/Figure1/comb_Maincelltype_paired_comparison.pdf",height = 8,width=6)



df_sub <- df[df$PatientID%in%patientID_mono,]
p=ggplot(df_sub, aes(x = Timepoint, y = value)) +
  geom_boxplot(aes(fill = Timepoint), show.legend = F, width = 0.6,alpha=0.8) +  #箱线图
  geom_point(aes(fill = Timepoint),size = 3,shape=21) +  #绘制散点
  scale_fill_manual(values = c("#9F89B2","#564478")) +  #设置颜色
  geom_line(aes(group = PatientID), color = '#BBB4C9', lwd = 0.5) +  #配对样本间连线
  theme(panel.grid = element_blank(),
        axis.line = element_line(colour = 'black', size = 1),
        panel.background = element_blank(),
        plot.title = element_text(size = 15, hjust = 0.5),
        plot.subtitle = element_text(size = 15, hjust = 0.5),
        axis.text = element_text(size = 15, color = 'black'),
        axis.title = element_text(size = 15, color = 'black')) +
  labs(x = '', y = 'gene.expression', title = 'Baseline vs Treatment')+
  facet_wrap(.~celltype,                     # 同时做几个图，这几个图分别对应数据转换前的第1列到第4列
             scale="free")+
  stat_compare_means(method = "wilcox.test",paired = TRUE,
                     comparisons=list(c("Baseline", "Treatment")))

p  
ggsave(p,file="final_figure/Figure1/mono_Maincelltype_paired_comparison.pdf",height = 8,width=6)




