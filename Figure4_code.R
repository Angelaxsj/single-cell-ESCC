#####===================计算各细胞组分================####
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
library(cowplot)
library(reshape2)
library(pheatmap)
library(readxl)
library(readr)
library(xlsx)
library(coin)


######==========================Data import===========================#####
sc_Tumor_NMF <- readRDS("~/Project/scESCC/CellType/Epi/Rdata/sc_Tumor_NMF.rds")
sc_Tumor_clean <- readRDS("~/Project/scESCC/CellType/Epi/Rdata/sc_Tumor_clean.rds")
sampleinfo <- read.table(file = "~/Project/scESCC/ESCC_metadata_20240705.txt",header = T,sep = "\t")
signature_df <- read.table(file = "CellType/Epi/tables/signature_gene_df.txt",header = T,sep = "\t")
sc_Fibroblast_clean <- readRDS("~/Project/scESCC/CellType/Fibroblast/sc_Fibroblast_clean.rds")
######==========================Figure 4A===========================#####
library(AUCell)
sce <- readRDS(file="final_figure/Figure4/sc_Tumor_score.rds")
sce <- metadata_combine(sc_Tumor_NMF,sampleinfo)


dotplot <- DotPlot(sce_sub_base,features = "MP6",group.by = "orig.ident")
data <- dotplot$data

data <- merge(data,sampleinfo,by.x="id",by.y="SampleID")

data$Response <- as.factor(data$Response)

pdf("/home/xiangsj/Project/scESCC/final_figure/Figure4/Figure4A.pdf",width = 4.5,height = 5)
ggplot(data, 
          aes(Response,
              avg.exp.scaled,  
              color = Response, 
              fill = Response))+ 
  geom_jitter(width=.15,size=4,alpha=0.6)+ 
  geom_boxplot(width=.4, alpha=.4,aes(fill=Response),color="black")+ #
  scale_color_manual(values = c("#3C7DA9","#DC9442"))+ 
  scale_fill_manual(values = c("#3C7DA9","#DC9442"))+ 
  theme_test()+ 
  labs(y="Length of petal", x="")+ 
 
  scale_y_continuous(expand = c(0.1,0.1))+# 横纵坐标不要统一成一样的，各个图根据自己情况来自动设置横纵坐标的范围
  theme(strip.background = element_blank(), #把头上的灰框去掉
        axis.text = element_text(color="black", size=10),
        axis.line = element_blank())
dev.off()






geneset <- as.list(signature_df)

geneset <- lapply(geneset,function(x){x=c(x,"EPCAM", "SFN","CDH1")})

#
######==========================Figure 4C===========================#####
library(survival)
library(survminer)
res.cut=surv_cutpoint(data =df,time = "Diseasefreesurvival_time",event = "Diseasefreesurvival_status",
                      variables = "MP6" ,minprop = 0.1)
summary(res.cut)
res.cut=surv_categorize(res.cut)
head(res.cut)
fit=survfit(Surv(Diseasefreesurvival_time,Diseasefreesurvival_status)~MP6,data = res.cut)
pdf("final_figure/Figure4/Figure4C_li.etal.pdf",width = 3.5 ,height=5,onefile = FALSE)
ggsurvplot(fit = fit,data = res.cut,
           
           conf.int = F,##添加置信区间
           surv.median.line = "hv",
           pval = T,
           pval.method = T,
           risk.table = T,
           tables.height = 0.2,
           tables.theme = theme_cleantable(),
           palette =  c("#E4305D","#4B6DB1"),
           size=2,
           ggtheme = theme_bw())+labs(title = "MP6_score")

dev.off()

######==========================Figure 4D===========================#####
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
cellchat_count <- cellchat_all@net$count
pdf(file = "final_figure/Figure4/Figure4D_MP_to_TME.pdf",width = 5,height = 6)
netVisual_circle(cellchat_all@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= T, title.name = "Interaction weights/strength",sources.use = c(1),targets.use = c(2:10))
dev.off()

######==========================Figure 4E===========================#####
cellchat <- readRDS("final_figure/Figure4/MP6_fibrolast_cellchat.rds")

#####===================关注MP6和Fibroblast受体配体对的情况=======================#####
responder_cells <- rownames(scRNA_cellchat@meta.data)[scRNA_cellchat$Treatment%in%"Chemo+PD1"&
                                                        scRNA_cellchat$Timepoint%in%"Baseline"&
                                                        scRNA_cellchat$Response%in%"responder"]
nonresponder_cells <- rownames(scRNA_cellchat@meta.data)[scRNA_cellchat$Treatment%in%"Chemo+PD1"&
                                                           scRNA_cellchat$Timepoint%in%"Baseline"&
                                                           scRNA_cellchat$Response%in%"non-responder"]
data <- as.matrix(scRNA_cellchat@assays$RNA@data)
data.use1 <- data[,responder_cells]

data.use2 <- data[,nonresponder_cells]
# cell.use3 <- rownames(scRNA_cellchat@meta.data)[scRNA_cellchat$subgroup=="Monotreat"]
# data.use3 <- data[,cell.use3]
metadata <- scRNA_cellchat@meta.data
metadata.use1 <- metadata[responder_cells ,]
metadata.use2 <- metadata[nonresponder_cells ,]
#metadata.use3 <- metadata[cell.use3 ,]


data <- data.use1
metadata <- metadata.use1
cellchat1 <- createCellChat(object = data,
                            meta = metadata,
                            group.by = "celltype")
cellchat1@DB <- CellChatDB.human
cellchat1 <- subsetData(cellchat1) %>% identifyOverExpressedGenes() %>% identifyOverExpressedInteractions() %>%
  computeCommunProb(raw.use = T,population.size = T) %>%
  filterCommunication( min.cells = 10) %>%
  computeCommunProbPathway() %>%
  aggregateNet() %>% netAnalysis_computeCentrality(slot.name = "netP")
### 2.2  Combtreatment cellchat construction
data <- data.use2
metadata <- metadata.use2
cellchat2 <- createCellChat(object = data,
                            meta = metadata,
                            group.by = "celltype")
cellchat2@DB <- CellChatDB.human
cellchat2 <- subsetData(cellchat2) %>% identifyOverExpressedGenes() %>% identifyOverExpressedInteractions() %>%
  computeCommunProb(raw.use = T,population.size = T) %>%
  filterCommunication( min.cells = 10) %>%
  computeCommunProbPathway() %>%
  aggregateNet() %>% netAnalysis_computeCentrality(slot.name = "netP")


###==================================细胞通讯结果可视化=============================###
###合并cellchat对象
cellchat.list <- list(Responder=cellchat1,Non_Responder=cellchat2)
cellchatcomb <- mergeCellChat(cellchat.list,add.names = names(cellchat.list),cell.prefix = T)
###可以选择需要可视化的细胞
levels(cellchatcomb@idents$joint)
unique(CellChatDB.human[["interaction"]][["pathway_name"]])
pdf("final_figure/Figure4/Figure4E.pdf",width = 5.5,height = 5.5)
netVisual_bubble(cellchatcomb,
                 sources.use = c(1),targets.use = c(2:7),
                 comparison = c(1,2),angle.x = 45,
                 signaling = c("COLLAGEN"),
                 remove.isolate=F,
                 #max.dataset = 1,
                 #return.data = T,
                 title.name = "Decreased signaling after TP+anti-PD1")
dev.off()









