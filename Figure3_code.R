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
library(ggprism)
library(ComplexHeatmap)
library(dplyr)
library(harmony)
library(RColorBrewer)
library(AUCell)

######==========================Figure3A===========================#####
ggplot(data = nmf_intersect_meltI_NEW, aes(x=Var1, y=Var2, fill=100*value/(100-value), color=100*value/(100-value))) + 
  geom_tile() + 
  scale_color_gradient2(limits=c(2,25), low=custom_magma[1:111],  mid =custom_magma[112:222], high = custom_magma[223:333], midpoint = 13.5, oob=squish, name="Similarity\n(Jaccard index)") +                                
  scale_fill_gradient2(limits=c(2,25), low=custom_magma[1:111],  mid =custom_magma[112:222], high = custom_magma[223:333], midpoint = 13.5, oob=squish, name="Similarity\n(Jaccard index)")  +
  theme( axis.ticks = element_blank(), panel.border = element_rect(fill=F), panel.background = element_blank(),  axis.line = element_blank(), axis.text = element_text(size = 11), axis.title = element_text(size = 12), legend.title = element_text(size=11), legend.text = element_text(size = 10), legend.text.align = 0.5, legend.justification = "bottom") + 
  
  guides(fill = guide_colourbar(barheight = 4, barwidth = 1))

######==========================Figure3B===========================#####

data_scale<-as.matrix(apply(expr,1,scale))
data_scale <- t(data_scale)
heat_colors2 <- colorRamp2(c(-2.5, 0, 2.5), c("blue", "white", "red"))
Heatmap(data_scale,
                  name="Expression",
                  clustering_distance_columns="pearson",
                  cluster_rows = FALSE,
                  col=heat_colors2,
                  color_space = "RGB",
                  cluster_columns = T,
                  row_order=NULL,
                  column_order=NULL,
                  show_column_names = FALSE,
                  #show_row_names = FALSE,
                  #row_names_gp = gpar(fontsize = 5),
                  #split=matSplit,
                  gap = unit(1, "mm"),
                  column_title = "",
                  column_title_gp = gpar(fontsize = 5),
                  width=unit(8, "cm"),
                  show_heatmap_legend = T,
                  # heatmap_legend_param=list(labels_gp = gpar(fontsize = 9),
                  #                           title_gp = gpar(fontsize = 9, fontface = "bold"))
) 
ht1org



######==========================Figure3C===========================#####
pdf(file = "final_figure/Figure3/Figure3C.pdf",height = 5,width = 4)
ggscatter(df_cor, x = "CD8T", y = "Epithelial",
          add = "reg.line", conf.int = TRUE,size=2,
          add.params = list(fill = "#6DCCFD")
          
)+
  stat_cor(method = "spearman"
  )
dev.off()

######==========================Figure3D===========================#####
pvalue_df <- data.frame(combined=c(-1.07881,-3.015935,-1.386294,-1.237874,-0.4307829,5.278515),mono=c(0.01005034,-0.1984509,-1.427116,1.771957,0.02020271,-0.01005034))
heat_colors2 <- colorRamp2(c(-5,-2.5, 0,2.5, 5), c("#282A62", "#2B6BA5", "white","#F6C63C","#D96558"))
pdf(file = "final_figure/Figure3/Figure3C.pdf",width = 2,height = 4.5)
Heatmap(pvalue_df,
        name="Expression",
        clustering_distance_columns="pearson",
        cluster_rows = FALSE,
        col=heat_colors2,
        color_space = "RGB",cell_fun=function(j,i,x,y,width,height,fill){
          grid.text(sprintf("%.1f",pvalue_df[i,j]),x,y,gp=gpar(fontsize=10))})
dev.off()



#######==========================Figure3F===========================#####
MP2_cytokine_correlation_df <- readRDS(file = "final_figure/Figure3/MP2_cytokine_correlation_df.rds")
MP2_cytokine_correlation_df$cohort <- rownames(MP2_cytokine_correlation_df)
MP2_cytokine_correlation_df$log10p <- -log10(MP2_cytokine_correlation_df$pvalue) 
MP2_cytokine_correlation_df <- MP2_cytokine_correlation_df[order(MP2_cytokine_correlation_df$log10p,decreasing = T),]
MP2_cytokine_correlation_df$cohort <- factor(MP2_cytokine_correlation_df$cohort,levels = c("NC","NICE-1","TCGA","li_etal"))
pdf(file = "final_figure/Figure3/Figure3J_correlation.pdf",width = 6,height = 4)
ggplot(MP2_cytokine_correlation_df,
       aes(x = cohort, y = log10p)) +
geom_segment(aes(x = cohort, xend = cohort, y = 0, yend = log10p),
               linetype = "solid", 
               size = 1,
               color = "black" 
  ) +
  # y轴0刻度的水平线
  geom_hline(
    yintercept = -log10(0.05),  # 水平线位置
    linetype = "dashed",  # 虚线
    size = 1, # 连线的粗细
    colour="#E4007F" # 连线的颜色
  ) +
  geom_point(aes(color = cohort,size = abs(correlation)),
             
  )+
  scale_size_continuous(range = c(10, 16))+ ###控制点的大小范围
  theme_classic() +
  theme(
    plot.title = element_text(size = 14, hjust = 0.5),
    axis.text.x = element_text(size = 12, angle = 45, hjust=1, face = 'bold'),
    axis.text.y = element_text(size = 12, face = 'bold'),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 12, face = 'bold'),
    plot.margin = unit(c(1, 1, 1, 1), "cm")
  )+
  scale_color_manual(values = c("#CB6C73","#996FA8","#9ED0C8","#E6D1A2"))
dev.off()
#######==========================Figure3G===========================#####
cellchat <- createCellChat(object = data,
                           meta = metadata,
                           group.by = "celltype")
cellchat@DB <- CellChatDB.human
cellchat <- subsetData(cellchat) %>% identifyOverExpressedGenes() %>% identifyOverExpressedInteractions() %>%
  computeCommunProb(raw.use = T,population.size = T) %>%
  filterCommunication( min.cells = 10) %>%
  computeCommunProbPathway() %>%
  aggregateNet() %>% netAnalysis_computeCentrality(slot.name = "netP")
groupSize <- as.numeric(table(cellchat@idents))
levels(cellchat@idents)
pdf(file = "final_figure/Figure3/Figure3G_MP_toTcell.pdf",width = 5,height = 6)
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength",sources.use = c(1:11),targets.use = c(12:17))
dev.off()
pdf(file = "final_figure/Figure3/Figure3G_Tcell_toMP.pdf",width = 5,height = 6)
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength",targets.use = c(1:11),sources.use = c(12:17))
dev.off()


#######==========================Figure3H===========================#####
cellchat.list <- list(Baseline=cellchat1,Treatment=cellchat2)
cellchatcomb <- mergeCellChat(cellchat.list,add.names = names(cellchat.list),cell.prefix = T)
levels(cellchatcomb@idents$joint)
unique(CellChatDB.human[["interaction"]][["pathway_name"]])
pdf("final_figure/Figure3/Figure3H_LR_CD8_to_Tumor.pdf",width = 9,height = 5)
netVisual_bubble(cellchatcomb,sources.use = c(1:11),targets.use = c(13),comparison = c(1,2),angle.x = 45,
                       signaling = c("TIGIT","CD6","CD46","CD226","CD96","TNF","SEMA4","MIF","CD39","SEMA4"
                                     ),

                       #remove.isolate=T,#max.dataset = 1,
                       #return.data = T,
                       title.name = "Decreased signaling after TP+anti-PD1")
dev.off()
pdf("final_figure/Figure3/Figure3H_LR_Tumor_to_CD8T.pdf",width = 9,height = 8)
p2=netVisual_bubble(cellchatcomb,targets.use = c(1:11),sources.use = c(13),comparison = c(1,2),angle.x = 45,
                 signaling = c("MHC-I","TNF","TIGIT","NECTIN","ALCAM","WNT","ITGB2","JAM","L1CAM","ICAM","LAMININ","XCR"
                               ,"ICAM","IFN-I"),
                 
                 #remove.isolate=T,#max.dataset = 1,
                 return.data = T,
                 title.name = "Decreased signaling in Responder")
p2
dev.off()

pc <- p1+p2
pc
ggsave(pc,file="final_figure/Figure2/Figure3G.pdf",width = 14,height=7)



#######==========================Figure3K===========================#####
exp <- GetAssayData(sc_Tumor_clean,assay = "RNA",slot = "data")
exp <- as.matrix(exp)
exp <- as.data.frame(t(exp))
exp <- exp[exp$NECTIN2!=0,]

df_cor <- merge(exp,sce@meta.data,by=0)
df_cor[1:5,1:5]
cor.test(df_cor$NECTIN2,df_cor$stemness,method = "spearman")


df_cor$NECTIN2 <- normalize(df_cor$NECTIN2)
pdf(file = "final_figure/Figure3/Figure3K.pdf",width = 5,height = 5)

ggscatter(df_cor, x = "NECTIN2", y = "stemness",
add = "reg.line", conf.int = TRUE,size=0.5,
add.params = list(fill = "#6DCCFD")
)+
stat_cor(method = "spearman"
)
dev.off()








