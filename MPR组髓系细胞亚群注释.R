library(Seurat)
library(tidyverse)
library(harmony)
library(dplyr)
library(patchwork)
library(SeuratWrappers)

MPRNMPR_object_miMPRobe<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/MPRNMPR_object_miMPRobe_v1_harmony.rds")
MPRNMPR_object_miMPRobe = UpdateSeuratObject(MPRNMPR_object_miMPRobe)
MPRNMPR_object_miMPRobe
View(MPRNMPR_object_miMPRobe)

MPRNMPR_object_miMPRobe_remove<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/MPRNMPR_object_miMPRobe_remove.rds")
MPRNMPR_object_miMPRobe_remove = UpdateSeuratObject(MPRNMPR_object_miMPRobe_remove)
MPRNMPR_object_miMPRobe_remove

setwd("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/髓系细胞亚型/")
myeloid_cells_object<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/myeloid_cells_object.rds")

myeloid_cells_object=readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/myeloid_cells_object_filter使用.rds")
myeloid_cells_object = UpdateSeuratObject(myeloid_cells_object)
myeloid_cells_object

####提取myeloid cells subject----
Idents(MPRNMPR_object_miMPRobe)<-"cell_type_new"

myeloid_cells_object=subset(MPRNMPR_object_miMPRobe,idents = c('Myeloid cells',"Mast cells"))

table(Idents(myeloid_cells_object))

colnames(myeloid_cells_object@meta.data)

pbmc<- myeloid_cells_object

pbmc <- NormalizeData(pbmc,normalization.method = "LogNormalize",
                      scale.factor = 1e4)
GetAssay(pbmc,assay = "RNA")
pbmc <- FindVariableFeatures(pbmc,
                             selection.method = "vst",nfeatures = 2000)
pbmc<- ScaleData(pbmc)
pbmc <- RunPCA(object = pbmc,pc.genes = VariableFeatures(pbmc))
pbmc <- FindNeighbors(pbmc,dims = 1:20)
for (res in c(0.01, 0.05, 0.1, 0.2 ,0.3,0.4,0.5, 0.6,0.8, 1)) {
  pbmc = FindClusters(pbmc,
                      resolution = res, algorithm = 1)
}
table(pbmc$RNA_snn_res.0.4)
sel.clust <- "RNA_snn_res.0.4"

pbmc <- SetIdent(pbmc,value = sel.clust)
metadata <- pbmc@meta.data

pbmc <- RunUMAP(pbmc,dims = 1:20)#UMAP二次降维，取前20个维度
pbmc <- RunTSNE(pbmc, dims = 1:20, reduction.key = "tsne_")
plot3 = DimPlot(pbmc, reduction = "umap", label=T) #查看clusters在UMAP降维图中的分布
plot3
dev.new()
plot4 = DimPlot(pbmc, reduction = "umap", group.by='orig.ident') #查看每个样本在UMAP降维图中的分布

plot5 = DimPlot(pbmc, reduction = "umap", split.by='orig.ident') #查看每个样本在UMAP降维图中的分面图
plotc <- plot4+plot5
plotc
#harmony分析***
pbmc <- RunHarmony(object =pbmc, reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")#harmony需要pca降维的结果，而scRNA对象此前已经做过PCA降维了，所以此处不再重复
pbmc <- RunUMAP(pbmc, reduction = "harmony", dims = 1:20,reduction.name = "umap.harmony")
pbmc <- RunTSNE(object = pbmc,reduction = "harmony", reduction.key = "tsneharmony_",dims.use = 1:20,reduction.name = "tsne.harmony")
pbmc <- FindNeighbors(pbmc, reduction = "umap.harmony", dims = 1:2)
pbmc <- FindClusters(pbmc)#标准聚类

#绘图和保存
p1 <- DimPlot(pbmc, group.by = "orig.ident",reduction = "umap.harmony", pt.size=0.1) +
  ggtitle("Integrated by harmony")#去批次后的可视化
p1
p2 <- DimPlot(pbmc, group.by = "RNA_snn_res.0.4",reduction = "umap.harmony", pt.size=0.1,label = T) +
  ggtitle("RNA_snn_res.0.4")#去批次后的可视化
p2

p3 <- DimPlot(pbmc, group.by = "RNA_snn_res.0.4",reduction = "tsne.harmony", pt.size=0.1,label = T) +
  ggtitle("RNA_snn_res.0.4")#去批次后的可视化
p3

p4 <- DimPlot(pbmc, group.by = "RNA_snn_res.0.4",reduction = "tsne.harmony", pt.size=0.1,label = T) +
  ggtitle("Integrated by harmony")#去批次后的可视化
p4

p5 <- DimPlot(pbmc, group.by = "orig.ident",reduction = "tsne.harmony", pt.size=0.1)
p5

p6 <- DimPlot(pbmc, group.by = "orig.ident",reduction = "tsne_", pt.size=0.1)
p6

p7 <- DimPlot(pbmc, group.by = "RNA_snn_res.0.4",reduction = "umap", pt.size=0.1,label = T)

p7

p8 <- DimPlot(pbmc, group.by = "RNA_snn_res.0.4",reduction = "umap", pt.size=0.1)

p8

DefaultAssay(pbmc) <- "RNA"

markers <- FindAllMarkers(myeloid_cells_object, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
all.markers =markers %>% dplyr::select(gene, everything()) %>% subset(p_val<0.05)
top20 =all.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
Idents(myeloid_cells_object) <- myeloid_cells_object$RNA_snn_res.0.4
p_heatmap<-DoHeatmap(myeloid_cells_object,features = top20$gene)+ #热图复现
  theme(text = element_text(size = 5))

write.csv(all.markers,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/髓系细胞亚型/myeloid_cells_marker_all.markers.csv")
ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/髓系细胞亚型/myeloid_cells_topmarkergene_heatmap.jpg",p_heatmap, width = 10,height = 15)

myeloid_cells_object<-pbmc
#保存RDS
saveRDS(myeloid_cells_object,"D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/myeloid_cells_object.rds")


####marker基因小提琴图----
setwd("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/髓系细胞亚型/")
myeloid_cells_object<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/myeloid_cells_object.rds")

cell_marker<- read.csv("CXCR2.csv")
######删除找不到的基因
nofinder<- c("TPSAB2", "CD209A", "ILRB4", "LY6I", "MS4A4C", "CX3CRL", "LY6C", "CD43SPN", "CD62L","SELL", 
             "CD43", "CD11B", "CD115", "CD62L", "CD11C", "LY6G", "RETNLG", "FCNB", "CHIL3", "NGP", "CCL6", "STFA2L1",
             "CD169", "C10A", "C1QG", "IFN", "IDO")

cell_marker_filtered <- cell_marker[!cell_marker$gene_name %in% nofinder, ]
Idents(myeloid_cells_object)<-"RNA_snn_res.0.4"
unique(cell_marker_filtered$gene_name)
library(randomcoloR)

my36colors = distinctColorPalette(84)


my36colors <-c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
               '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
               '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
               '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
               '#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
               '#968175'
)


####ggplot2画小提琴图----

vln.dat=FetchData(myeloid_cells_object,c(unique(cell_marker_filtered$gene_name),"RNA_snn_res.0.4"))
colnames(vln.dat)
#宽转长
vln.dat.melt <- reshape2::melt(vln.dat, id.vars = c("RNA_snn_res.0.4"), 
                               measure.vars = unique(cell_marker_filtered$gene_name),
                               variable.name = "gene", 
                               value.name = "Expr") 

vln.dat.melt$gene<- factor(vln.dat.melt$gene,levels = unique(vln.dat.melt$gene))

View(vln.dat.melt)

library(cowplot)
vln.dat.melt$RNA_snn_res.0.4<- factor(vln.dat.melt$RNA_snn_res.0.4,levels = c("0","1","2","3","4","5","6",
                                                                              "7","8","9","10","11","12","13","14","15",
                                                                              "16","17","18","19","20","21"))
vln.dat.melt$cells_type <- factor(vln.dat.melt$cells_type,levels=c("Lateral root cap", "Cortex","Endodermis","Pericycle","Xylem",
                                                                   "Atrichoblast","Root hairs","Unknow" ))


colnames(vln.dat.melt)
library(cowplot)

p1 <- ggplot(vln.dat.melt, aes(gene, Expr, fill = gene)) +
  geom_violin(scale = "width", adjust = 1, trim = TRUE) +
  scale_y_continuous(expand = c(0, 0), position="right", labels = function(x)
    c(rep(x = "", times = length(x)-2), x[length(x) - 1], "")) +
  facet_grid(rows = vars(RNA_snn_res.0.4), scales = "free", switch = "y") +
  scale_fill_manual(values = my36colors) + 
  theme_cowplot(font_size = 10) +
  theme(legend.position = "none", panel.spacing = unit(0, "lines"),
        plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = NA, color = "black"),
        plot.margin = margin(7, 7, 0, 7, "pt"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        strip.text.y.left = element_text(angle = 0),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, color = "black")
  ) +ylab("Expression Level")
p1

ggsave("myeloid_marker小提琴图RNA_snn_res.0.4_9.jpg",p1,width = 7,height = 7)
ggsave("myeloid_marker小提琴图RNA_snn_res.0.4.pdf",p1,width = 7,height = 7)


####featureplot----


FeaturePlot(object=myeloid_cells_object,reduction = "umap.harmony",raster=FALSE,features =c("TPSB2",
                                                                                       "TPSAB1",
                                                                                       "TPSAB2",
                                                                                       "CPA3",
                                                                                       "MS4A2",
                                                                                       "KIT")) 

FeaturePlot(object=myeloid_cells_object,reduction = "umap.harmony",raster=FALSE,features =c("CD14",
                                                                                            "CD300E",
                                                                                            "CD244",
                                                                                            "HLA-DRA",
                                                                                            "S100A4",
                                                                                            "S100A12"
                                                                                            
                                                                                            )) 
VlnPlot(myeloid_cells_object, features = "LY9",raster=FALSE)

VlnPlot(myeloid_cells_object, 
        features = c("CXCR2"),
        pt.size = 0,
        ncol = 2) 

cluster20<- subset(myeloid_cells_object,idents = c("20"))
FeaturePlot(object=cluster20,reduction = "umap.harmony",raster=FALSE,features =c("CD14",
                                                                                            "CD300E",
                                                                                            "CD244",
                                                                                            "HLA-DRA",
                                                                                            "S100A4",
                                                                                            "S100A12"
                                                                                            
)) 


####命名cluster----
###去除cluster20
Idents(myeloid_cells_object)<-"RNA_snn_res.0.4"
myeloid_cells_object<- subset(myeloid_cells_object,idents=c(20),invert=TRUE)
levels(Idents(myeloid_cells_object))

old_to_new_groups <- c("19" = "Mast cells",
                       "16" = "Mast cells",
                       "4" ="Mast cells",
                       "21" = "Mast cells",
                       "0"="Macrophages",
                       "1"="Macrophages",
                       "2"="Macrophages",
                       "3"="Macrophages",
                       "5"="Macrophages",
                       "10"="Macrophages",
                       "15"="Macrophages",
                       "8"="Macrophages",
                       "9"="Macrophages",
                       "11"="Macrophages",
                       "12"="Macrophages",
                       "13"="Macrophages",
                       "17"="Macrophages",
                       "14"="Migratory cDCs",#migratory cDCs conventional dendritic cells
                       "6"="cDC2",#Dendritic cells
                       "7"="Monocytes",
                       "18"="Monocytes")


myeloid_cells_object@meta.data$myeloid_cells_type <- plyr::mapvalues(myeloid_cells_object@meta.data$RNA_snn_res.0.4,
                                                         from = names(old_to_new_groups), 
                                                         to = old_to_new_groups)

table(myeloid_cells_object@meta.data$myeloid_cells_type)

DimPlot(myeloid_cells_object, group.by = "myeloid_cells_type",reduction = "umap.harmony", pt.size=0.1)

saveRDS(myeloid_cells_object,"D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/myeloid_cells_object.rds")



Idents(myeloid_cells_object)<-"myeloid_cells_type" 

macrophage_object<- subset(myeloid_cells_object,idents = c("Macrophages"))


####细分macrophage----

####marker基因小提琴图----
setwd("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/髓系细胞亚型/")
myeloid_cells_object<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/myeloid_cells_object.rds")

cell_marker<- read.csv("marker10.csv")
######删除找不到的基因
nofinder<- c("TPSAB2", "CD209A", "ILRB4", "LY6I", "MS4A4C", "CX3CRL", "LY6C", "CD43SPN", "CD62L","SELL", 
             "CD43", "CD11B", "CD115", "CD62L", "CD11C", "LY6G", "RETNLG", "FCNB", "CHIL3", "NGP", "CCL6", "STFA2L1",
             "CD169", "C10A", "C1QG", "IFN", "IDO")

cell_marker_filtered <- cell_marker[!cell_marker$gene_name %in% nofinder, ]
macrophage_object<-subset(myeloid_cells_object_filter,idents = c("0","1","2","3","5","10_macro","15",
"8", "11_macro")) 

cluster8<-subset(myeloid_cells_object_filter,idents = c("8")) 

Idents(macrophage_object)<-"RNA_snn_res.0.4"
levels(Idents(macrophage_object))
unique(cell_marker_filtered$gene_name)
library(randomcoloR)

my36colors = distinctColorPalette(84)

vln.dat=FetchData(macrophage_object,c(unique(cell_marker_filtered$gene_name),"RNA_snn_res.0.4"))
colnames(vln.dat)
#宽转长
vln.dat.melt <- reshape2::melt(vln.dat, id.vars = c("RNA_snn_res.0.4"), 
                               measure.vars = unique(cell_marker_filtered$gene_name),
                               variable.name = "gene", 
                               value.name = "Expr") 

vln.dat.melt$gene<- factor(vln.dat.melt$gene,levels = unique(vln.dat.melt$gene))

View(vln.dat.melt)

library(cowplot)
vln.dat.melt$RNA_snn_res.0.4<- factor(vln.dat.melt$RNA_snn_res.0.4,levels = c("0","1","2","3","5","10_macro","15",
                                                                              "8", "11_macro"))

colnames(vln.dat.melt)
library(cowplot)

p1 <- ggplot(vln.dat.melt, aes(gene, Expr, fill = gene)) +
  geom_violin(scale = "width", adjust = 1, trim = TRUE) +
  scale_y_continuous(expand = c(0, 0), position="right", labels = function(x)
    c(rep(x = "", times = length(x)-2), x[length(x) - 1], "")) +
  facet_grid(rows = vars(RNA_snn_res.0.4), scales = "free", switch = "y") +
  scale_fill_manual(values = my36colors) + 
  theme_cowplot(font_size = 10) +
  theme(legend.position = "none", panel.spacing = unit(0, "lines"),
        plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = NA, color = "black"),
        plot.margin = margin(7, 7, 0, 7, "pt"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        strip.text.y.left = element_text(angle = 0),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, color = "black")
  ) +ylab("Expression Level")
p1

ggsave("myeloid_marker小提琴图RNA_snn_res.0.4_macrophage1.jpg",p1,width = 7,height = 7)


FeaturePlot(object=myeloid_cells_object,reduction = "umap.harmony",raster=FALSE,features =c("IFN1")) 

FeaturePlot(object=cluster8,reduction = "umap.harmony",raster=FALSE,features =c("ZG16B")) 

old_to_new_groups <- c("19" = "Mast cells",
                       "16" = "Mast cells",
                       "4" ="Mast cells",
                       "21" = "Mast cells",
                       "0"="Macro-IL1B",#Macrophages
                       "1"="Macro-SPP1+",#Macrophages
                       "2"="Macro-IL1B",#Macrophages
                       "3"="Macro-SPP1+",#Macrophages
                       "5"="Macro-C1QB+",#Macrophages
                       "10"="Macro-SPP1+",#Macrophages
                       "15"="Macro-C1QB+",#Macrophages
                       "8"="Macro-IL1B",#Macrophages
                       "9"="Macro-SPP1+",#Macrophages
                       "11"="Macro-ISG15+",#Macrophages
                       "12"="Macro-LA",#Macrophages
                       "13"="Macro-LA",#Macrophages
                       "17"="Macro-SPP1+",#Macrophages
                       "14"="Migratory cDCs",#migratory cDCs conventional dendritic cells
                       "6"="cDC2",#Dendritic cells
                       "7"="Monocytes",
                       "18"="Monocytes")


myeloid_cells_object@meta.data$myeloid_cells_type <- plyr::mapvalues(myeloid_cells_object@meta.data$RNA_snn_res.0.4,
                                                                     from = names(old_to_new_groups), 
                                                                     to = old_to_new_groups)

table(myeloid_cells_object@meta.data$myeloid_cells_type)
Idents(myeloid_cells_object)<-"myeloid_cells_type"

saveRDS(myeloid_cells_object,"D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/myeloid_cells_object.rds")

DoHeatmap(macrophage_object, features = unique(cell_marker_filtered$gene_name), size = 3)

DotPlot(macrophage_object, features = unique(cell_marker_filtered$gene_name)) + RotatedAxis()

FeaturePlot(object=myeloid_cells_object,reduction = "umap.harmony",raster=FALSE,features =c("CD163"
)) 


####拆分cluster----

Idents(myeloid_cells_object)<-"RNA_snn_res.0.4"
cluster11<- subset(myeloid_cells_object,idents = c("11"))

cluster11$RNA_snn_res.0.4 = ifelse(cluster11@assays$RNA@counts['CD14',]>0&cluster11@assays$RNA@counts['CD68',]>0 ,'11_macro','11_neg')
myeloid_cells_object$RNA_snn_res.0.4 <- as.character(myeloid_cells_object$RNA_snn_res.0.4)
cluster11$RNA_snn_res.0.4 <- as.character(cluster11$RNA_snn_res.0.4)

myeloid_cells_object$RNA_snn_res.0.4[match(colnames(cluster11),colnames(myeloid_cells_object))]=cluster11$RNA_snn_res.0.4

table(myeloid_cells_object$RNA_snn_res.0.4)

Idents(myeloid_cells_object)<-"RNA_snn_res.0.4"


cluster10<- subset(myeloid_cells_object,idents = c("10"))

cluster10$RNA_snn_res.0.4 = ifelse(cluster10@assays$RNA@counts['CD14',]>0&cluster10@assays$RNA@counts['CD68',]>0 ,'10_macro','10_neg')

myeloid_cells_object$RNA_snn_res.0.4 <- as.character(myeloid_cells_object$RNA_snn_res.0.4)

cluster10$RNA_snn_res.0.4 <- as.character(cluster10$RNA_snn_res.0.4)

myeloid_cells_object$RNA_snn_res.0.4[match(colnames(cluster10),colnames(myeloid_cells_object))]=cluster10$RNA_snn_res.0.4

table(myeloid_cells_object$RNA_snn_res.0.4)

####删除其他cluster----

saveRDS(myeloid_cells_object,"D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/myeloid_cells_object.rds")

Idents(myeloid_cells_object_filter)<-"RNA_snn_res.0.4"
levels(Idents(myeloid_cells_object_filter))

myeloid_cells_object_filter <- subset(myeloid_cells_object,idents = c("10","11","20","9","12","13","17","8",),invert=TRUE)
table(Idents(myeloid_cells_object_filter))

markers <- FindAllMarkers(myeloid_cells_object_filter, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
all.markers =markers %>% dplyr::select(gene, everything()) %>% subset(p_val<0.05)
top20 =all.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
Idents(myeloid_cells_object_filter) <- myeloid_cells_object_filter$RNA_snn_res.0.4

write.csv(all.markers,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/髓系细胞亚型/myeloid_cells_marker_all.markers_过滤.csv")
write.csv(top20,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/髓系细胞亚型/myeloid_cells_marker_top20_过滤.csv")


saveRDS(myeloid_cells_object_filter,"D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/myeloid_cells_object_filter.rds")

####重新聚类----
pbmc<- myeloid_cells_object_filter

pbmc <- NormalizeData(pbmc,normalization.method = "LogNormalize",
                      scale.factor = 1e4)
GetAssay(pbmc,assay = "RNA")
pbmc <- FindVariableFeatures(pbmc,
                             selection.method = "vst",nfeatures = 2000)
pbmc<- ScaleData(pbmc)
pbmc <- RunPCA(object = pbmc,pc.genes = VariableFeatures(pbmc))
pbmc <- FindNeighbors(pbmc,dims = 1:20)
for (res in c(0.01, 0.05, 0.1, 0.2 ,0.3,0.4,0.5, 0.6,0.8, 1)) {
  pbmc = FindClusters(pbmc,
                      resolution = res, algorithm = 1)
}
table(pbmc$RNA_snn_res.0.4)
sel.clust <- "RNA_snn_res.0.4"

pbmc <- SetIdent(pbmc,value = sel.clust)
metadata <- pbmc@meta.data
pbmc <- RunUMAP(pbmc,dims = 1:20)
pbmc <- RunUMAP(pbmc,dims = 1:20)#UMAP二次降维，取前20个维度
pbmc <- RunTSNE(pbmc, dims = 1:20, reduction.key = "tsne_")
plot3 = DimPlot(pbmc, reduction = "umap", label=T) #查看clusters在UMAP降维图中的分布
plot3
dev.new()
plot4 = DimPlot(pbmc, reduction = "umap", group.by='orig.ident') #查看每个样本在UMAP降维图中的分布

plot5 = DimPlot(pbmc, reduction = "umap", split.by='orig.ident') #查看每个样本在UMAP降维图中的分面图
plotc <- plot4+plot5
plotc
#harmony分析***
pbmc <- RunHarmony(object =pbmc, reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")#harmony需要pca降维的结果，而scRNA对象此前已经做过PCA降维了，所以此处不再重复
pbmc <- RunUMAP(pbmc, reduction = "harmony", dims = 1:20,reduction.name = "umap.harmony")
pbmc <- RunTSNE(object = pbmc,reduction = "harmony", reduction.key = "tsneharmony_",dims.use = 1:20,reduction.name = "tsne.harmony")
pbmc <- FindNeighbors(pbmc, reduction = "umap.harmony", dims = 1:2)
pbmc <- FindClusters(pbmc)#标准聚类

#绘图和保存
p1 <- DimPlot(pbmc, group.by = "orig.ident",reduction = "umap.harmony", pt.size=0.1) +
  ggtitle("Integrated by harmony")#去批次后的可视化
p1
p2 <- DimPlot(pbmc, group.by = "RNA_snn_res.0.4",reduction = "umap.harmony", pt.size=0.1,label = T) +
  ggtitle("RNA_snn_res.0.4")#去批次后的可视化
p2
p3 <- DimPlot(pbmc, group.by = "RNA_snn_res.0.4",reduction = "tsne.harmony", pt.size=0.1,label = T) +
  ggtitle("RNA_snn_res.0.4")#去批次后的可视化
p3
p4 <- DimPlot(pbmc, group.by = "RNA_snn_res.0.4",reduction = "tsne.harmony", pt.size=0.1,label = T) +
  ggtitle("Integrated by harmony")#去批次后的可视化
p4
p5 <- DimPlot(pbmc, group.by = "orig.ident",reduction = "tsne.harmony", pt.size=0.1)
p5
p6 <- DimPlot(pbmc, group.by = "orig.ident",reduction = "tsne_", pt.size=0.1)
p6
p7 <- DimPlot(pbmc, group.by = "RNA_snn_res.0.4",reduction = "umap", pt.size=0.1,label = T)
p7
p8 <- DimPlot(pbmc, group.by = "RNA_snn_res.0.4",reduction = "umap", pt.size=0.1)
p8
DefaultAssay(pbmc) <- "RNA"

markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
all.markers =markers %>% dplyr::select(gene, everything()) %>% subset(p_val<0.05)
top20 =all.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)

write.csv(all.markers,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/髓系细胞亚型/myeloid_cells_marker_all.markers_过滤.csv")

myeloid_cells_object_filter <- pbmc
saveRDS(myeloid_cells_object_filter,"D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/myeloid_cells_object_filter.rds")

####marker小提琴图----
setwd("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/髓系细胞亚型")
cell_marker<- read.csv("marker10_TAMs.csv")
######删除找不到的基因
nofinder<- c("TPSAB2", "CD209A", "ILRB4", "LY6I", "MS4A4C", "CX3CRL", "LY6C", "CD43SPN", "CD62L","SELL", 
             "CD43", "CD11B", "CD115", "CD62L", "CD11C", "LY6G", "RETNLG", "FCNB", "CHIL3", "NGP", "CCL6", "STFA2L1",
             "CD169", "C10A", "C1QG", "IFN", "IDO")

cell_marker_filtered <- cell_marker[!cell_marker$gene_name %in% nofinder, ]

Idents(pbmc) <- "RNA_snn_res.0.4"
DoHeatmap(pbmc, features = unique(cell_marker_filtered$gene_name), size = 3)

VlnPlot(pbmc, features = unique(cell_marker_filtered$gene_name),raster=FALSE)

DotPlot(pbmc, features = unique(cell_marker_filtered$gene_name)) + RotatedAxis()

FeaturePlot(object=pbmc,reduction = "umap.harmony",raster=FALSE,features =unique(cell_marker_filtered$gene_name)) 

vln.dat=FetchData(pbmc,c(unique(cell_marker_filtered$gene_name),"RNA_snn_res.0.4"))
colnames(vln.dat)

vln.dat.melt <- reshape2::melt(vln.dat, id.vars = c("RNA_snn_res.0.4"), 
                               measure.vars = unique(cell_marker_filtered$gene_name),
                               variable.name = "gene", 
                               value.name = "Expr") 

vln.dat.melt$gene<- factor(vln.dat.melt$gene,levels = unique(vln.dat.melt$gene))

View(vln.dat.melt)

library(cowplot)

colnames(vln.dat.melt)
my36colors<- distinctColorPalette(118)
p1 <- ggplot(vln.dat.melt, aes(gene, Expr, fill = gene)) +
  geom_violin(scale = "width", adjust = 1, trim = TRUE) +
  scale_y_continuous(expand = c(0, 0), position="right", labels = function(x)
    c(rep(x = "", times = length(x)-2), x[length(x) - 1], "")) +
  facet_grid(rows = vars(RNA_snn_res.0.4), scales = "free", switch = "y") +
  scale_fill_manual(values = my36colors) + 
  theme_cowplot(font_size = 10) +
  theme(legend.position = "none", panel.spacing = unit(0, "lines"),
        plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = NA, color = "black"),
        plot.margin = margin(7, 7, 0, 7, "pt"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        strip.text.y.left = element_text(angle = 0),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, color = "black")
  ) +ylab("Expression Level")
p1
ggsave("myeloid_marker小提琴图RNA_snn_res.0.4_删除之后marker.jpg",p1,width =15,height = 7)

myeloid_cells_object_filte <- pbmc

saveRDS(myeloid_cells_object_filter,"D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/myeloid_cells_object_filter.rds")
myeloid_cells_object_filter<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/myeloid_cells_object_filter.rds")
####重新命名cluster----

FeaturePlot(object=myeloid_cells_object_filter,reduction = "umap.harmony",raster=FALSE,features =c("VEGFA",
                                                                                                   "SPP1"
                                                                                                   
                                                                                                   
)) 

Idents(myeloid_cells_object_filter)<-"RNA_snn_res.0.4"

levels(Idents(myeloid_cells_object_filter))

old_to_new_groups <- c("5" = "Mast cells",
                       "11" = "Mast cells",
                       "14" ="Mast cells",
                       "15" = "Mast cells",
                       "0"="Macrophages",
                       "2"="Macrophages",
                       "3"="Macrophages",
                       "4"="Macrophages",
                       "6"="Macrophages",
                       "7"="Macrophages",
                       "10"="Macrophages",
                       "1"="Monocyte-derived dendritic cell", 
                       "9"="Migratory cDCs",#migratory cDCs conventional dendritic cells
                       "13"="cDC1",
                       "8"="cDC2",#Dendritic cells
                       "3"="Monocytes",
                       "12"="Monocytes")

old_to_new_groups <- c("5" = "Mast cells",
                       "11" = "Mast cells",
                       "14" ="Mast cells",
                       "15" = "Mast cells",
                       "0"="ZBTB16+TAMs",
                       "2"="SPP1+TAMs",
                       "3"="SPP1+TAMs",
                       "4"="CXCL10+TAMs",
                       "6"="FOLR2+TAMs",
                       "7"="SPP1+TAMs",
                       "10"="FOLR2+TAMs",
                       "1"="Monocyte-derived dendritic cells", 
                       "9"="Migratory cDCs",#migratory cDCs conventional dendritic cells
                       "13"="cDC1",
                       "8"="cDC2",#Dendritic cells
                       "12"="Neutrophil")



myeloid_cells_object_filter@meta.data$myeloid_cells_type <- plyr::mapvalues(myeloid_cells_object_filter@meta.data$RNA_snn_res.0.4,
                                                                     from = names(old_to_new_groups), 
                                                                     to = old_to_new_groups)

table(myeloid_cells_object_filter@meta.data$myeloid_cells_type)
unique(myeloid_cells_object_filter@meta.data$myeloid_cells_type)

DimPlot(myeloid_cells_object_filter, group.by = "myeloid_cells_type",reduction = "umap.harmony", pt.size=0.1)
Idents(myeloid_cells_object_filter)<-"myeloid_cells_type"

DoHeatmap(myeloid_cells_object_filter, features = unique(cell_marker_filtered$gene_name), size = 3)

DotPlot(myeloid_cells_object_filter, features = unique(cell_marker_filtered$gene_name)) + RotatedAxis()

FeaturePlot(object=myeloid_cells_object,reduction = "umap.harmony",raster=FALSE,features =c("CD14"
)) 

saveRDS(myeloid_cells_object_filter,"D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/myeloid_cells_object_filter.rds")

markers <- FindAllMarkers(myeloid_cells_object_filter, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
all.markers =markers %>% dplyr::select(gene, everything()) %>% subset(p_val<0.05)
top20 =all.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)

write.csv(all.markers,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/髓系细胞亚型/myeloid_cells_marker_all.markers.cluster重命名.csv")

write.csv(top20,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/髓系细胞亚型/myeloid_cells_marker_top20.cluster重命名.csv")

####细胞注释正式画图----
setwd("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/髓系细胞亚型/")
myeloid_cells_object_filter<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/myeloid_cells_object_filter.rds")
myeloid_cells_object_filter = UpdateSeuratObject(myeloid_cells_object_filter)
library(tidyverse)
myeloid_cells.harmonyredu<- myeloid_cells_object_filter@reductions$umap.harmony@cell.embeddings%>%
  as.data.frame() %>% 
  cbind(cell_type = myeloid_cells_object_filter@meta.data$myeloid_cells_type)%>% 
  cbind(Patients = myeloid_cells_object_filter@meta.data$Patients) %>% 
  cbind(Treatments = myeloid_cells_object_filter@meta.data$Treatments) %>% 
  cbind(MPR_Response = myeloid_cells_object_filter@meta.data$MPR_Response) %>%
  cbind(MPR_Response2 = myeloid_cells_object_filter@meta.data$MPR_Response2) %>%
  cbind(cell_id = myeloid_cells_object_filter@meta.data$cell_id) %>%
  cbind(orig.ident = myeloid_cells_object_filter@meta.data$orig.ident)

names(myeloid_cells.harmonyredu)

colnames(myeloid_cells.harmonyredu) <- c("umapharmony_1","umapharmony_2","cell_type","Patients",
                                          "Treatments","MPR_Response", "MPR_Response2","cell_id",
                                          "orig.ident")

write.csv(myeloid_cells.harmonyredu,"myeloid_cells_umap降维数据_过滤.csv")


###颜色设置###画小坐标轴####https://cloud.tencent.com/developer/article/1924260
library(randomcoloR)
cell_type_color = distinctColorPalette(10)
cell_type_color<- c("#0000FF","#FFa07a","#00bfff","#7fff00","#8a2be2","#00ffff","#9bcd9b","#FF00ff","#ffd700","#8b795e")

unique(myeloid_cells.harmonyredu$cell_type)

myeloid_cells.harmonyredu$cell_type <- factor(myeloid_cells.harmonyredu$cell_type,levels= c("Mast cells","SPP1+TAMs","FOLR2+TAMs",
                                                                                            "ZBTB16+TAMs","CXCL10+TAMs","Neutrophil",
                                                                                            "Monocyte-derived dendritic cells",
                                                                                            "cDC1","cDC2","Migratory cDCs"))
library(ggplot2)

myeloid_cell_type <- ggplot(myeloid_cells.harmonyredu,aes(x= umapharmony_1, y =umapharmony_2,color = cell_type)) +  
  geom_point(size = 0.5, alpha =1)+  
  scale_color_manual(values = cell_type_color)+
  xlab("UMAP1")+
  ylab("UMAP2")+
  theme(panel.grid.major = element_blank(), #主网格线
        panel.grid.minor = element_blank(), #次网格线
        panel.border = element_blank(), #边框
        axis.title = element_blank(),  #轴标题
        axis.text = element_blank(), # 文本
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = 'white'), #背景色
        plot.background=element_rect(fill="white"))+
  theme(legend.title = element_blank(), #去掉legend.title 
        legend.key=element_rect(fill='white'), #
        legend.text = element_text(size=20), #设置legend标签的大小
        legend.key.size=unit(1,'cm'),
        legend.position = "None") +  # 设置legend标签之间的大小
  guides(color = guide_legend(override.aes = list(size=5)))+ #设置legend中 点的大小 
  geom_segment(aes(x = min(myeloid_cells.harmonyredu$umapharmony_1) , y = min(myeloid_cells.harmonyredu$umapharmony_2) ,
                   xend = min(myeloid_cells.harmonyredu$umapharmony_1) +3, yend = min(myeloid_cells.harmonyredu$umapharmony_2) ),
               colour = "black", size=1,arrow = arrow(length = unit(0.5,"cm")))+ 
  geom_segment(aes(x = min(myeloid_cells.harmonyredu$umapharmony_1)  , y = min(myeloid_cells.harmonyredu$umapharmony_2)  ,
                   xend = min(myeloid_cells.harmonyredu$umapharmony_1) , yend = min(myeloid_cells.harmonyredu$umapharmony_2) +3),
               colour = "black", size=1,arrow = arrow(length = unit(0.5,"cm"))) +
  annotate("text", x = min(myeloid_cells.harmonyredu$umapharmony_1) +1.5, y = min(myeloid_cells.harmonyredu$umapharmony_2)-1, label = "UMPA_1",
           color="black",size = 4, fontface="bold" ) + 
  annotate("text", x = min(myeloid_cells.harmonyredu$umapharmony_1)-1, y = min(myeloid_cells.harmonyredu$umapharmony_2) + 1.5, label = "UMAP_2",
           color="black",size = 4, fontface="bold" ,angle=90)

myeloid_cell_type

myeloid_cell_type_med <- myeloid_cells.harmonyredu %>%
  group_by(cell_type) %>%
  summarise(
    umapharmony_1 = median(umapharmony_1),
    umapharmony_2 = median(umapharmony_2)
  )

library(ggrepel)
myeloid_cell_type_med

myeloid_cell_type2<-myeloid_cell_type +geom_label_repel(aes(label=cell_type),size=4,color="black",fontface="bold",data = myeloid_cell_type_med,
                                                      point.padding=unit(0.5, "lines"),fill = alpha(c("white"),0.5),
                                                      segment.size=0.5, nudge_x=1, nudge_y = 0.5,direction="y",max.overlaps=50)+
  theme(legend.position = "right")


ggsave('./myeloid_cell_type.jpg',myeloid_cell_type2,width =13,height = 6)
ggsave('./myeloid_cell_type.pdf',myeloid_cell_type2,width =13,height = 6)

Idents(myeloid_cells_object_filter)<-"myeloid_cells_type"

#### 细胞注释marker 热图----
cell_marker<- read.csv("marker画图.csv")
markers <-cell_marker$gene_name
markers <- as.data.frame(markers)

vln.dat=FetchData(myeloid_cells_object_filter,c(unique(cell_marker$gene_name),"myeloid_cells_type"))
colnames(vln.dat)
#宽转长
vln.dat.melt <- reshape2::melt(vln.dat, id.vars = c("myeloid_cells_type"), 
                               measure.vars = unique(cell_marker$gene_name),
                               variable.name = "gene", 
                               value.name = "Expr") 
vln.dat.melt_mean <- aggregate(vln.dat.melt$Expr, by = list(vln.dat.melt$myeloid_cells_type, vln.dat.melt$gene), FUN = mean)

colnames(vln.dat.melt_mean) <- c("myeloid_cells_type","gene","value")
vln.dat.melt<-vln.dat.melt_mean
aver_dt_df_long<- vln.dat.melt
colnames(aver_dt_df_long)
View(aver_dt_df_long)
write.csv(aver_dt_df_long,"myeloid_cells_object_filter_marker热图数据.csv")
####热图顶部色块
library(tidyverse)
library(ggnewscale)
library(MetBrewer)
library(patchwork)
library(ggtext)
marker_group_list<-c("Mast cells","SPP1+TAMs","FOLR2+TAMs",
                     "ZBTB16+TAMs","CXCL10+TAMs","Neutrophil",
                     "Monocyte-derived dendritic cells",
                     "cDC1","cDC2","Migratory cDCs")

aver_dt_df_long<- read.csv("myeloid_cells_object_filter_marker热图数据.csv",header = T,row.names = 1)
aver_dt_df_long$gene <- factor(aver_dt_df_long$gene,levels =rev(unique(aver_dt_df_long$gene)))
aver_dt_df_long$myeloid_cells_type <- factor(aver_dt_df_long$myeloid_cells_type,levels =marker_group_list)


p1_heat<-ggplot(aver_dt_df_long) +
  geom_tile(aes(myeloid_cells_type,gene,fill=value), colour = "white", linewidth = 0.5)+
  scale_fill_gradientn(name='Mean scaled expression',
                       #colours=c("white","red","black"))+
                       colours=c("blue","yellow","red"))+
  theme_minimal() +
  theme(axis.line=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = 10,color = 'black',
        ),
        panel.border=element_rect(fill=NA,color="black",size=0.5,linetype="solid"),
        axis.text.x = element_text(size = 10,color = 'black',angle =40,hjust =1,vjust=1),
        legend.title = element_text(size=10, color = "black"),
        legend.text = element_text(size=10,color = "black",angle =0),
        legend.position = "right") +
  scale_y_discrete(position = "left")
p1_heat


sub_markers_count$id<- factor(sub_markers_count$id, levels = marker_group_list)

pdot <-
  ggplot(sub_markers_count,aes(x = features.plot,y = id)) +
  geom_point(aes(fill =avg.exp.scaled,size = pct.exp),
             color = 'black',
             shape = 21) +
  theme_bw(base_size = 14) +
  xlab('') + ylab('') +
  scale_fill_gradient2(low = 'white',mid = '#EB1D36',high = '#990000',
                       midpoint = 5,
                       limits = c(0,10),
                       #breaks = seq(0,2,10),
                       #labels = seq(0,2,10),
                       name = 'Average Expression') +
  scale_size(range = c(1,5),
             limits = c(0, 100),
             breaks = seq(20,100,20),
             labels = seq(20,100,20)
  ) +
  theme(panel.grid = element_blank(),
        legend.title = element_text(size=12,face="bold"),
        legend.text=element_text(size=9,face="bold"),
        axis.text = element_text(color = 'black'),
        aspect.ratio = 0.5,
        plot.margin = margin(t = 0.5,r = 0.5,b = 0.5,l = 0.5,unit = 'cm'),
        axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5, face = 'bold.italic',size=15),
        axis.text.y = element_text(color = "black",size=15,face="bold")
  ) +
  coord_cartesian(clip = 'off') +
  guides(
    fill = guide_colorbar(
      direction = "horizontal",
      title.position = "top",
      barwidth = unit(5, "cm")
    ),
    size = guide_legend(
      title = "Percent Expressed",
      direction = "horizontal",
      title.position = "top",
      label.position = "bottom",
      override.aes = list(
        color = "black",
        fill = "grey"
      )
    )
  )

pdot

####dotplot数据画热图----
cell_marker<- read.csv("marker画图.csv")
p1 <- DotPlot(myeloid_cells_object_filter, features =unique(cell_marker$gene_name),
              assay='RNA' )
p1
sub_markers_count<- p1$data

names(sub_markers_count)
colnames(sub_markers_count)
sub_markers_count<- sub_markers_count %>% data.frame()
sub_markers_count$features.plot <- factor(sub_markers_count$features.plot,levels =rev(unique(sub_markers_count$features.plot)))
sub_markers_count$id <- factor(sub_markers_count$id,levels =marker_group_list)
p1_heat<-ggplot(sub_markers_count) +
  geom_tile(aes(id,features.plot,fill=avg.exp.scaled), colour = "white", linewidth = 0.5)+
  scale_fill_gradientn(name='Mean scaled expression',
                       #colours=c("white","red","black"))+
                       colours=c("blue","yellow","red"))+
  theme_minimal() +
  theme(axis.line=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = 15,color = 'black',
        ),
        panel.border=element_rect(fill=NA,color="black",size=0.5,linetype="solid"),
        axis.text.x = element_text(size = 15,color = 'black',angle =40,hjust =1,vjust=1),
        legend.title = element_text(size=15, color = "black"),
        legend.text = element_text(size=15,color = "black",angle =0),
        legend.position = "right") +
  scale_y_discrete(position = "left")
p1_heat

df_right <- sub_markers_count %>% group_by(features.plot) %>% filter (! duplicated(features.plot))
View(df_right)
df_right$id<- factor(df_right$id,levels =marker_group_list)
df_right$features.plot <- factor(df_right$features.plot,levels =rev(unique(df_right$features.plot)))
df_right$group <- cell_marker[match(cell_marker$gene_name,df_right$features.plot),1]
df_right$group <- factor(df_right$group,levels =c("Mast cell","Macrophages","TAMs",
                                                  "Neutrophil","Monocytes","DC",
                                                  "cDC1","cDC2" ,"Migratory cDCs"))

top_color<- c("#0000FF","#cd96cd","#eee5de","#00ffff","#9bcd9b","#CD853F","#FF00ff","#ffd700","#8b795e")


View(df_right)

p1_right <- df_right %>% ggplot(aes(x=id,y=features.plot))+
  geom_tile(data=df_right,color="black",aes(fill=group))+
  scale_fill_manual(values=top_color)+
  new_scale_fill()+
  labs(x=NULL,y=NULL)+
  theme(axis.text.x = element_blank(),
        axis.text.y=element_blank(),
        plot.margin=unit(c(0,1,-20,0.5),unit="cm"),
        plot.background = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank(),
        legend.text = element_text(color="black"),
        legend.position = "none")
p1_right

p1_right<- p1_right + theme(plot.margin = unit(c(1, 1, 0.5, 0.5), "cm"))
p1_heat <- p1_heat + theme(plot.margin = unit(c(0.5, 0.5, 1, 1), "cm"))

library(cowplot)
p_data_plot<-plot_grid(
                       p1_heat,
                       p1_right,
                       align = 'h',
                       nrow=1,
                       rel_heights= c(1,1),
                       rel_widths = c(2,0.3))

p_data_plot


ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/髓系细胞亚型/myeloid_cells_细胞亚型分类热图.jpg",p_data_plot,height =13,width = 8)
ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/髓系细胞亚型/myeloid_cells_细胞亚型分类热图.pdf",p_data_plot,height =13,width = 8)


####髓系细胞不过滤分类----
setwd("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/髓系细胞亚型/髓系细胞不过滤分类")

myeloid_cells_object<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/myeloid_cells_object.rds")

p1 <- DimPlot(myeloid_cells_object, group.by = "orig.ident",reduction = "umap.harmony", pt.size=0.1) +
  ggtitle("Integrated by harmony")#去批次后的可视化
p1
p2 <- DimPlot(myeloid_cells_object, group.by = "RNA_snn_res.0.4",reduction = "umap.harmony", pt.size=0.1,label = T) +
  ggtitle("RNA_snn_res.0.4")#去批次后的可视化
p2

Idents(myeloid_cells_object)<-"RNA_snn_res.0.4"

cell_marker<- read.csv("marker画图.csv")
p1 <- DotPlot(myeloid_cells_object, features =unique(cell_marker$gene_name),
              assay='RNA' )
p1
sub_markers_count<- p1$data

names(sub_markers_count)
colnames(sub_markers_count)
sub_markers_count<- sub_markers_count %>% data.frame()
sub_markers_count$features.plot <- factor(sub_markers_count$features.plot,levels =rev(unique(sub_markers_count$features.plot)))
sub_markers_count$id <- factor(sub_markers_count$id,levels =marker_group_list)
p1_heat<-ggplot(sub_markers_count) +
  geom_tile(aes(id,features.plot,fill=avg.exp.scaled), colour = "white", linewidth = 0.5)+
  scale_fill_gradientn(name='Mean scaled expression',
                       #colours=c("white","red","black"))+
                       colours=c("blue","yellow","red"))+
  theme_minimal() +
  theme(axis.line=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = 12,color = 'black',
        ),
        panel.border=element_rect(fill=NA,color="black",size=0.5,linetype="solid"),
        axis.text.x = element_text(size = 15,color = 'black',angle =40,hjust =1,vjust=1),
        legend.title = element_text(size=15, color = "black"),
        legend.text = element_text(size=15,color = "black",angle =0),
        legend.position = "right") +
  scale_y_discrete(position = "left")
p1_heat

df_right <- sub_markers_count %>% group_by(features.plot) %>% filter (! duplicated(features.plot))
View(df_right)
df_right$id<- factor(df_right$id,levels =marker_group_list)
df_right$features.plot <- factor(df_right$features.plot,levels =rev(unique(df_right$features.plot)))
df_right$group <- cell_marker[match(cell_marker$gene_name,df_right$features.plot),1]
df_right$group <- factor(df_right$group,levels =c("Mast cell","Macrophages","TAMs",
                                                  "Neutrophil","Monocytes","DC",
                                                  "cDC1","cDC2" ,"Migratory cDCs"))

top_color<- c("#0000FF","#cd96cd","#eee5de","#00ffff","#9bcd9b","#CD853F","#FF00ff","#ffd700","#8b795e")


View(df_right)

p1_right <- df_right %>% ggplot(aes(x=id,y=features.plot))+
  geom_tile(data=df_right,color="black",aes(fill=group))+
  scale_fill_manual(values=top_color)+
  new_scale_fill()+
  labs(x=NULL,y=NULL)+
  theme(axis.text.x = element_blank(),
        axis.text.y=element_blank(),
        plot.margin=unit(c(0,1,-20,0.5),unit="cm"),
        plot.background = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank(),
        legend.text = element_text(color="black"),
        legend.position = "none")
p1_right

p1_right<- p1_right + theme(plot.margin = unit(c(1, 1, 0.5, 0.5), "cm"))
p1_heat <- p1_heat + theme(plot.margin = unit(c(0.5, 0.5, 1, 1), "cm"))

library(cowplot)
p_data_plot<-plot_grid(
  p1_heat,
  p1_right,
  align = 'h',
  nrow=1,
  rel_heights= c(1,1),
  rel_widths = c(2,0.3))

p_data_plot


ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/髓系细胞亚型/髓系细胞不过滤分类/myeloid_cells_细胞亚型分类热图.jpg",p_data_plot,height =13,width = 8)
ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/髓系细胞亚型/髓系细胞不过滤分类/myeloid_cells_细胞亚型分类热图.pdf",p_data_plot,height =13,width = 8)


####拆分cluster----

Idents(myeloid_cells_object)<-"RNA_snn_res.0.4"
cluster11<- subset(myeloid_cells_object,idents = c("11"))

cluster11$RNA_snn_res.0.4 = ifelse(cluster11@assays$RNA@counts['CD14',]>0&cluster11@assays$RNA@counts['CD68',]>0 ,'11_macro','11_neg')
myeloid_cells_object$RNA_snn_res.0.4 <- as.character(myeloid_cells_object$RNA_snn_res.0.4)
cluster11$RNA_snn_res.0.4 <- as.character(cluster11$RNA_snn_res.0.4)

myeloid_cells_object$RNA_snn_res.0.4[match(colnames(cluster11),colnames(myeloid_cells_object))]=cluster11$RNA_snn_res.0.4

table(myeloid_cells_object$RNA_snn_res.0.4)

Idents(myeloid_cells_object)<-"RNA_snn_res.0.4"


cluster10<- subset(myeloid_cells_object,idents = c("10"))

cluster10$RNA_snn_res.0.4 = ifelse(cluster10@assays$RNA@counts['CD14',]>0&cluster10@assays$RNA@counts['CD68',]>0 ,'10_macro','10_neg')

myeloid_cells_object$RNA_snn_res.0.4 <- as.character(myeloid_cells_object$RNA_snn_res.0.4)

cluster10$RNA_snn_res.0.4 <- as.character(cluster10$RNA_snn_res.0.4)

myeloid_cells_object$RNA_snn_res.0.4[match(colnames(cluster10),colnames(myeloid_cells_object))]=cluster10$RNA_snn_res.0.4

table(myeloid_cells_object$RNA_snn_res.0.4)


saveRDS(myeloid_cells_object,"D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/myeloid_cells_object.rds")

####有无微生物细胞差异分析----
MPRNMPR_object_miMPRobe_remove<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/MPRNMPR_object_miMPRobe_remove.rds")
MPRNMPR_object_miMPRobe_remove = UpdateSeuratObject(MPRNMPR_object_miMPRobe_remove)
MPRNMPR_object_miMPRobe_remove

Idents(MPRNMPR_object_miMPRobe_remove) <- "cell_type_new"
MPRNMPR_object_miMPRobe_remove@meta.data$MPR_Response2
Myeloid_cells_object=subset(MPRNMPR_object_miMPRobe_remove,idents =  c("Myeloid cells","Mast cells"))

Myeloid_cells_object@meta.data$group_microbe2
Idents(Myeloid_cells_object) <- "group_microbe2"
Myeloid_cells_object_deg_all=FindMarkers(Myeloid_cells_object, ident.1 = c("Bacteria+"), ident.2 = "Bacteria-")
dim(Myeloid_cells_object_deg_all)
names(Myeloid_cells_object_deg_all)
write.csv(Myeloid_cells_object_deg_all,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/免疫细胞/Myeloid_cells_object_deg_all.csv")
Myeloid_cells_object_deg_all<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/免疫细胞/Myeloid_cells_object_deg_all.csv",header = T)

####根据过滤掉pt 小于0.1 的行

names(Myeloid_cells_object_deg_all)

Myeloid_cells_object_deg_all_filter<- Myeloid_cells_object_deg_all %>% filter(pct.1>0.1 &pct.2>0.1)

####读取免疫细胞相关基因，从差异基因集中筛选

gene_name<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/免疫细胞/免疫细胞gene_name.csv")

intersection <- intersect(gene_name$Gene.symbol, Myeloid_cells_object_deg_all_filter$X)

Myeloid_cells_object_deg_all_filter_inter<- Myeloid_cells_object_deg_all_filter %>% filter(Myeloid_cells_object_deg_all_filter$X %in%intersection )

df <- Myeloid_cells_object_deg_all_filter_inter
View(df)
df$group2<-as.factor(ifelse(df$p_val_adj < 0.05 & abs(df$avg_log2FC) >= 0.25, 
                            ifelse(df$avg_log2FC>= 0.25 ,'Up','Down'),'NS'))
df$group2 <- factor(df$group2, levels = c("Up", "Down", "NS"))

##确定添加标签的数据

df$label<-ifelse(df$p_val_adj<0.05&abs(df$avg_log2FC)>=0.25,"Y","N")
df$label<-ifelse(df$label == 'Y', as.character(df$X), '')

library(tidyverse)
library(ggrepel)

library(ggplot2) ##绘图使用
library(ggprism) ##设置主题私用
library(ggrepel)

degp<- ggplot(df, aes(x =avg_log2FC, y=-log10(p_val_adj), colour=group2)) +
  geom_point(alpha=0.85, size=3) +  #点的透明度和大小
  scale_color_manual(values=c('brown','steelblue','gray')) + #调整点的颜色
  xlim(c(-1, 1)) +  ##调整x轴的取值范围，max(abs(BRCA_Match_DEG$logFC))最大值是多少，再进行取舍
  geom_vline(xintercept=c(-0.25,0.25),lty=4,col="black",lwd=0.8) + #添加x轴辅助线,lty函数调整线的类型
  geom_hline(yintercept = -log10(0.05), lty=4,col="black",lwd=0.8) +  #添加y轴辅助线
  labs(x="avg_log2FC", y="-log10(p_val_adj)") +  #x、y轴标签
  ggtitle("Myeloid cells DEG (Bacteria+ vs Bacteria-)") + #标题
  theme(plot.title = element_text(hjust = 0.5),legend.position="right",legend.title = element_blank())+
  guides(color = guide_legend(override.aes = list(size = 6)))+
  geom_label_repel(data = df, aes(label = label),##添加标签
                   size = 5,box.padding = unit(0.5, "lines"),
                   point.padding = unit(0.8, "lines"),
                   segment.color = "black",
                   show.legend = FALSE, max.overlaps = 10000,
                   fill = NA,        # 去除背景
                   color = "black",  # 标签文字颜色
                   label.size = NA )+# 去除边框+
  theme_prism(border = T) +
  theme(legend.text = element_text(size = 12))

degp
ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/免疫细胞/Myeloid细胞火山图.jpg",degp,width = 10,height=7)
ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/免疫细胞/Myeloid细胞火山图.pdf",degp,width = 10,height=7)


####髓系细胞OR指数----

library("sscVis")
library("data.table")
library("grid")
library("cowplot")
library("ggrepel")
library("readr")
library("plyr")
library("ggpubr")
library("ggplot2")
#设置图片输出目录
myeloid_cells_object
out.prefix <- "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/髓系细胞亚型/"
myeloid_cells_object<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/myeloid_cells_object_filter使用.rds")

Idents(myeloid_cells_object)<-"myeloid_cells_type"

TAMs_object=subset(myeloid_cells_object,idents = c("SPP1+TAMs","FOLR2+TAMs",
                                                            "ZBTB16+TAMs","CXCL10+TAMs"))

DCs_object=subset(myeloid_cells_object,idents = c("Monocyte-derived dendritic cells",
                                                  "cDC1","cDC2","Migratory cDCs"))
  
library(dplyr)

# 转换元数据为数据框
metadata <- as.data.frame(myeloid_cells_object@meta.data)
colnames(metadata)
unique(metadata$myeloid_cells_type)


do.tissueDist <- function(cellInfo.tb = meta.tb,
                          meta.cluster = meta.tb$cell_type_new,
                          colname.patient = "patient",
                          loc = meta.tb$MPR_Response,
                          out.prefix,
                          pdf.width=3,
                          pdf.height=5,
                          verbose=1){
  ##input data 
  library(data.table)
  dir.create(dirname(out.prefix),F,T)
  
  cellInfo.tb = data.table(cellInfo.tb)
  cellInfo.tb$meta.cluster = as.character(meta.cluster)
  
  if(is.factor(loc)){
    cellInfo.tb$loc = loc
  }else{cellInfo.tb$loc = as.factor(loc)}
  
  loc.avai.vec <- levels(cellInfo.tb[["loc"]])
  count.dist <- unclass(cellInfo.tb[,table(meta.cluster,loc)])[,loc.avai.vec]
  freq.dist <- sweep(count.dist,1,rowSums(count.dist),"/")
  freq.dist.bin <- floor(freq.dist * 100 / 10)
  print(freq.dist.bin)
  
  {
    count.dist.melt.ext.tb <- test.dist.table(count.dist)
    p.dist.tb <- dcast(count.dist.melt.ext.tb,rid~cid,value.var="p.value")
    OR.dist.tb <- dcast(count.dist.melt.ext.tb,rid~cid,value.var="OR")
    OR.dist.mtx <- as.matrix(OR.dist.tb[,-1])
    rownames(OR.dist.mtx) <- OR.dist.tb[[1]]
  }
  
  sscVis::plotMatrix.simple(OR.dist.mtx,
                            out.prefix=sprintf("%s.OR.dist",out.prefix),
                            show.number=F,
                            waterfall.row=T,par.warterfall = list(score.alpha = 2,do.norm=T),
                            exp.name=expression(italic(OR)),
                            z.hi=4,
                            palatte=viridis::viridis(7),
                            pdf.width = 4, pdf.height = pdf.height)
  if(verbose==1){
    return(list("count.dist.melt.ext.tb"=count.dist.melt.ext.tb,
                "p.dist.tb"=p.dist.tb,
                "OR.dist.tb"=OR.dist.tb,
                "OR.dist.mtx"=OR.dist.mtx))
  }else{
    return(OR.dist.mtx)
  }
}

test.dist.table <- function(count.dist,min.rowSum=0)
{
  count.dist <- count.dist[rowSums(count.dist)>=min.rowSum,,drop=F]
  sum.col <- colSums(count.dist)
  sum.row <- rowSums(count.dist)
  count.dist.tb <- as.data.frame(count.dist)
  setDT(count.dist.tb,keep.rownames=T)
  count.dist.melt.tb <- melt(count.dist.tb,id.vars="rn")
  colnames(count.dist.melt.tb) <- c("rid","cid","count")
  count.dist.melt.ext.tb <- as.data.table(ldply(seq_len(nrow(count.dist.melt.tb)), function(i){
    this.row <- count.dist.melt.tb$rid[i]
    this.col <- count.dist.melt.tb$cid[i]
    this.c <- count.dist.melt.tb$count[i]
    other.col.c <- sum.col[this.col]-this.c
    this.m <- matrix(c(this.c,
                       sum.row[this.row]-this.c,
                       other.col.c,
                       sum(sum.col)-sum.row[this.row]-other.col.c),
                     ncol=2)
    res.test <- fisher.test(this.m)
    data.frame(rid=this.row,
               cid=this.col,
               p.value=res.test$p.value,
               OR=res.test$estimate)
  }))
  count.dist.melt.ext.tb <- merge(count.dist.melt.tb,count.dist.melt.ext.tb,
                                  by=c("rid","cid"))
  #count.dist.melt.ext.tb[,adj.p.value:=p.adjust(p.value,"BH")]
  return(count.dist.melt.ext.tb)
}

# importFrom(data.table,":=")

myeloid_cells_object$myeloid_cells_type

meta.tb <- myeloid_cells_object@meta.data

TAMs_object$myeloid_cells_type
meta.tb <- TAMs_object@meta.data

DCs_object$myeloid_cells_type
meta.tb <- DCs_object@meta.data


names(meta.tb)
View(meta.tb)

A <- do.tissueDist(cellInfo.tb = meta.tb,
                   meta.cluster = meta.tb$myeloid_cells_type,
                   colname.patient = "patient",
                   loc = meta.tb$MPR_Response2,
                   out.prefix,
                   pdf.width=3,
                   pdf.height=5,
                   verbose=1)

#查看并保存文件
A$OR.dist.mtx #做热图数据，OR值
A$p.dist.tb #p值
A$OR.dist.tb #OR值
A$count.dist.melt.ext.tb#组合表，adjust-pvalue等

#自己做图
data <- A$count.dist.melt.ext.tb
unique(data$rid)
data$rid <- factor(data$rid, levels=c("Mast cells","SPP1+TAMs","FOLR2+TAMs",
                                      "ZBTB16+TAMs","CXCL10+TAMs","Neutrophil",
                                      "Monocyte-derived dendritic cells",
                                      "cDC1","cDC2","Migratory cDCs"))

data$rid <- factor(data$rid, levels=c("SPP1+TAMs","FOLR2+TAMs",
                                      "ZBTB16+TAMs","CXCL10+TAMs"))

data$rid <- factor(data$rid, levels=c("Monocyte-derived dendritic cells",
                                      "cDC1","cDC2","Migratory cDCs"))

cell_type_color<- c("#0000FF","#FFa07a","#00bfff","#7fff00","#8a2be2","#00ffff","#9bcd9b","#FF00ff","#ffd700","#8b795e")


bl <- colorRampPalette(c("navy","royalblue","lightskyblue"))(100)                      
re <- colorRampPalette(c("mistyrose", "red2","darkred"))(100)
library(RColorBrewer)
OR_plot<-ggplot(data, aes(rid,cid)) + 
  geom_tile(aes(fill = OR), colour = "black", size = 0.6)+
  scale_fill_gradientn(name='OR',
                       colours=colorRampPalette(rev(brewer.pal(n = 5, name ="RdYlBu")))(100))+
  theme_minimal() + 
  theme(axis.line=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(), 
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = 12,color = 'black'),
        axis.text.x = element_text(size = 12,color = 'black',angle =90,hjust =0.4,vjust=0.2),
        legend.title = element_text(size=12, color = "black"), 
        legend.text = element_text(size=12,color = "black")) + 
  scale_y_discrete(position = "right")

OR_plot

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/髓系细胞亚型/OR值/髓系_OR_plot_MPR.pdf",OR_plot,width=9,height=5)
ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/髓系细胞亚型/OR值/髓系_OR_plot_MPR_TAMs.pdf",OR_plot,width=9,height=5)

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/髓系细胞亚型/OR值/髓系_OR_plot_MPR_DCs.pdf",OR_plot,width=9,height=5)


####髓系细胞比例计算----
myeloid_cells_object$myeloid_cells_type

Cellratio <- prop.table(table(myeloid_cells_object@meta.data$myeloid_cells_type,myeloid_cells_object@meta.data$orig.ident), margin = 2)

Cellratio <- as.data.frame(Cellratio)
colnames(Cellratio) <- c("Cluster","Group","Freq")
library(reshape2)
cellper <- dcast(Cellratio,Group~Cluster, value.var = "Freq")#长数据转为宽数据
rownames(cellper) <- cellper[,1]
cellper <- cellper[,-1]

View(cellper)
rownames(cellper)

###t添加分组信息###
sample <- c("HTX_SC_pre","LM_SC_pre","LXD_SC_post","LXD_SC_pre","QHH_SC_pre",
            "WRY_SC_pre","XL_SC_post","XL_SC_pre","XYH_SC_post","XYH_SC_pre",
            "YMS_SC_post","YMS_SC_pre","YXC_SC_pre","YXJ_SC_pre")

MPR_Response <- c("NMPR","NMPR","MPR","MPR","MPR","MPR",
                  "MPR","MPR","MPR","MPR","MPR","MPR","NMPR","NMPR")
Treatments <- c("Pre-treat","Pre-treat","Post-treat","Pre-treat","Pre-treat","Pre-treat",
                "Post-treat","Pre-treat","Post-treat","Pre-treat","Post-treat","Pre-treat"
                ,"Pre-treat","Pre-treat")

MPR_Response2 <- c("Pre_NMPR","Pre_NMPR","Post_MPR","Pre_MPR","Pre_MPR",
                   "Pre_MPR","Post_MPR","Pre_MPR","Post_MPR","Pre_MPR",
                   "Post_MPR","Pre_MPR","Pre_NMPR","Pre_NMPR")
groups<- data.frame(sample,MPR_Response, Treatments,MPR_Response2)#创建数据框
rownames(groups)=groups$sample
cellper$MPR_Response <- groups[,"MPR_Response"]#R添加列
cellper$Treatments <- groups[,'Treatments']#R添加列
cellper$MPR_Response2 <- groups[,'MPR_Response2']
View(cellper)

names(cellper)


cell_type_groups = c("Mast cells","SPP1+TAMs","FOLR2+TAMs",
                     "ZBTB16+TAMs","CXCL10+TAMs","Neutrophil",
                     "Monocyte-derived dendritic cells",
                     "cDC1","cDC2","Migratory cDCs")
library(ggplot2)
library(dplyr)
library(ggpubr)
library(cowplot)

treatment_color <- c("#4974a4","#4dae47","#f29600")

pplist = list()
for(cell_type in cell_type_groups){
  cellper_cell_type = cellper %>% select(one_of(c('Treatments','MPR_Response',"MPR_Response2",cell_type)))#选择一组数据
  colnames(cellper_cell_type) = c('Treatments','MPR_Response','MPR_Response2','percent')#对选择数据列命名
  cellper_cell_type$percent = as.numeric(cellper_cell_type$percent)#数值型数据
  cellper_cell_type <- cellper_cell_type%>% group_by(MPR_Response) %>% mutate(upper =  quantile(percent, 0.75), 
                                                                              lower = quantile(percent, 0.25),
                                                                              mean = mean(percent),
                                                                              median = median(percent))#上下分位数
  print(cell_type)
  print(cellper_cell_type$median)
  
  cellper_cell_type$MPR_Response2 <- factor(cellper_cell_type$MPR_Response2,levels=c("Pre_NMPR","Pre_MPR","Post_MPR"))
  
  pp1 = ggplot(cellper_cell_type,aes(x=MPR_Response2,y=percent)) + #ggplot作图
    geom_boxplot(aes(fill=MPR_Response2),outlier.shape = NA,lwd= 0.7)+
    #geom_boxplot(outlier.colour="red", outlier.shape=7,outlier.size=1) +
    scale_fill_manual(values = treatment_color)+
    #geom_jitter(shape = 21,aes(fill=Treatments$Treatments),width = 0.25) + 
    #stat_summary(fun=mean, geom="point", color="grey60") +
    theme_cowplot() +
    theme(axis.text = element_text(size = 12),
          axis.line.x=element_line(linetype=1,color="black",size=1),
          axis.ticks.x=element_line(color="black",size=2,lineend = 1),
          axis.line.y=element_line(linetype=1,color="black",size=1),
          axis.ticks.y=element_line(color="black",size=2,lineend = 1),
          axis.title = element_text(size = 10),
          legend.text = element_text(size = 12),
          axis.text.x = element_text(vjust=1,hjust=1,angle=30,size=12),
          axis.text.y = element_text(vjust=1,hjust=1,angle=0,size=12),
          legend.title = element_text(size = 10),
          plot.title = element_text(size = 13),
          legend.position = 'none') + 
    labs(title = cell_type,y='Percentage',x="") +
    guides(fill = guide_legend(title = NULL))
  #geom_errorbar(aes(ymin = lower, ymax = upper),col = "grey60",width =  1)
  
  ###组间t检验分析
  my_comparisons <- list(c("Pre_MPR", "Pre_NMPR"),c("Pre_MPR", "Post_MPR"),c("Pre_NMPR", "Post_MPR"))
  # my_comparisons<-list(c("MPR", "NMPR"))
  # hide.ns = TRUE
  pp2<-pp1+stat_compare_means(comparisons=my_comparisons,
                              #                             #label.y = c(28, 32, 36),
                              method="wilcox.test",size=3.5,bracket.size = 1)
  # pp3<- pp2+stat_compare_means(aes(group = MPR_Response2),
  #                              method="wilcox.test",
  #                              label="p.signif",
  #                              show.legend = F,size=10)
  #label="p.signif",
  #show.legend = F)
  
  #label="p.signif",
  #show.legend = F)
  pplist[[cell_type]] = pp2
}

###"Pre-B cells", "Primarily mature naive B cells", "Naive B cells", "Active B cells", "Naive plasma cell", "Plasma cells"
# library(cowplot)


cells_<-plot_grid(pplist[['Mast cells']],
                  pplist[['SPP1+TAMs']],
                  pplist[['FOLR2+TAMs']],
                  pplist[['ZBTB16+TAMs']],
                  pplist[['CXCL10+TAMs']],
                  pplist[['Neutrophil']],
                  pplist[['Monocyte-derived dendritic cells']],
                  pplist[['cDC1']],
                  pplist[['cDC2']],
                  pplist[['Migratory cDCs']],
                  align = "h",  #axis = 'l',
                  nrow =2)

ggsave('D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/髓系细胞亚型/细胞百分比/myeloid_cells_MPR_percentage.pdf',cells_,width =12,height =8,limitsize = FALSE)



####
zbtb16_tams_metadata <- TAMs_object@meta.data[TAMs_object$myeloid_cells_type == "ZBTB16+TAMs",]
zbtb16_tams_metadata$group_microbe2
table(zbtb16_tams_metadata$group_microbe2)

####GSEA----
library(gggsea)
library(ggplot2)
library(GSVA)
library(enrichplot)
library(GSEABase)
library(fgsea)
library(clusterProfiler)
unique(myeloid_cells_object$MPR_Response2)
Idents(myeloid_cells_object)<-"MPR_Response2"
deg_all=FindMarkers(myeloid_cells_object, ident.1 = "Pre_MPR", ident.2 = "Pre_NMPR")
write.csv(deg_all, "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/髓系细胞亚型/M1M2特征得分和GESA分析/Pre_MPR_Pre_NMPR_deg_all.csv", row.names = T)
deg_all2=FindMarkers(myeloid_cells_object, ident.1 = "Post_MPR", ident.2 = "Pre_MPR")
write.csv(deg_all2, "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/髓系细胞亚型/M1M2特征得分和GESA分析/Post_MPR_Pre_MPR_deg_all.csv", row.names = T)
deg_all3=FindMarkers(myeloid_cells_object, ident.1 = "Post_MPR", ident.2 = "Pre_NMPR")
write.csv(deg_all3, "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/髓系细胞亚型/M1M2特征得分和GESA分析/Post_MPR_Pre_NMPR_deg_all.csv", row.names = T)

####pre_nmpr_pre_mpr_deg_gsea----

deg_all<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/髓系细胞亚型/M1M2特征得分和GESA分析/Pre_MPR_Pre_NMPR_deg_all.csv")
colnames(deg_all)[1]<- "gene_name"
deg<- deg_all
head(deg)

gmtfile <- "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/h.all.v2024.1.Hs.symbols.gmt"
txtfile<- "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/髓系细胞亚型/M1M2特征得分和GESA分析/特征得分marker.txt"

hallmark <- read.gmt(gmtfile)
hallmark <- read.gmt(txtfile)
head(hallmark)
View(hallmark)
alldiff <- deg[order(deg$avg_log2FC,decreasing = T),]
id <- alldiff$avg_log2FC
names(id) <- alldiff$gene_name
## 这里去掉了基因集前缀
hallmark$term <- gsub('HALLMARK_','',hallmark$term)
hallmark.list <- hallmark %>% split(.$term) %>% lapply( "[[", 2)

gsea.re1<-clusterProfiler::GSEA(id,TERM2GENE=hallmark,verbose = T)
g1<-as.data.frame(gsea.re1)
g1<-subset(g1,p.adjust<0.05)
g1<-g1[order(g1$NES,decreasing = T),]
gsea.re2 <- fgseaMultilevel(pathways = hallmark.list,
                            stats = id,
                            minSize=1,
                            maxSize=10000)
colnames(gsea.re2)

gsea.re2$pathway
gsea.re2$padj
gsea.re2$ES
gsea.re2$NES
sum(gsea.re2$padj < 0.05, na.rm = TRUE)
g2 <- gsea.re2[gsea.re2$padj< 0.05,]
g2 <- g2[order(g2$NES,decreasing = T),]
g2$padj
g2$pathway

####ggsea 可视化
names(hallmark.list)
se_hall<-c(head(g2$pathway,10),tail(g1$pathway,10))
sig1<-g2[g2$pathway%in%se_hall,]

hallmark.se<-hallmark.list[sig1$pathway]
df.new <- gseaCurve(id, hallmark.se, sig1)
pal_line<-pal_lancet()(9)#各曲线颜色

p_nes<-ggplot() + 
  geom_gsea(df.new,
            prettyGSEA=T,
            tickcolor='grey30',#竖线颜色
            ticksize=0.1,#线条粗细
            #colour='grey',#色块边框
            linecolor=pal_line,
            linesize=1,lty=1,
            #多个图形时使用，设定每行列图形个数
            nrow = 1,
            ncol = 3
  ) + 
  theme_bw()+
  theme(strip.text = element_text(size =10,face = 'bold.italic'),#增大分面标题字体
        strip.background = element_rect(fill = 'white'))+
  xlab(bquote(italic('Rank')))+ylab(bquote(italic('Enrichment Score')))+
  theme(axis.text.x = element_text(size=10,angle = -30,face = 'plain',hjust = 0.5),
        axis.text.y = element_text(size=10,angle = 0,face = 'plain',vjust = 0.5),
        axis.title.x = element_text(size=12,face = 'plain',hjust = 0.5),
        axis.title.y = element_text(size=12,face = 'plain',vjust = 0.5))
p_nes

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/髓系细胞亚型/M1M2特征得分和GESA分析/Pre_MPR_Pre_NMPR_gsea.pdf",p_nes,width=7,height=3)

library(ggsci)
col_gsea1<-pal_simpsons()(16)
num2=3
gseaplot2(gsea.re1,geneSetID = rownames(g1)[1:num2],
          title = "",#标题
          color = col_gsea1[1:num2],#颜色
          base_size = 14,#基础大小
          rel_heights = c(1, 0.2, 0.4),#小图相对高度
          subplots = 1:3,#展示小图
          pvalue_table = T,#p值表格
          ES_geom = "line"#line or dot
)

######post_mpr_pre_mpr----

# deg_all2=FindMarkers(myeloid_cells_object, ident.1 = "Post_MPR", ident.2 = "Pre_MPR")
# write.csv(deg_all2, "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/髓系细胞亚型/M1M2特征得分和GESA分析/Post_MPR_Pre_MPR_deg_all.csv", row.names = T)
# deg_all3=FindMarkers(myeloid_cells_object, ident.1 = "Post_MPR", ident.2 = "Pre_NMPR")
# write.csv(deg_all3, "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/髓系细胞亚型/M1M2特征得分和GESA分析/Post_MPR_Pre_NMPR_deg_all.csv", row.names = T)

deg_all2<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/髓系细胞亚型/M1M2特征得分和GESA分析/Post_MPR_Pre_MPR_deg_all.csv")
colnames(deg_all2)[1]<- "gene_name"
deg<- deg_all2
head(deg)

txtfile<- "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/髓系细胞亚型/M1M2特征得分和GESA分析/特征得分marker2.txt"

hallmark <- read.gmt(txtfile)
head(hallmark)
View(hallmark)

alldiff <- deg[order(deg$avg_log2FC,decreasing = T),]
id <- alldiff$avg_log2FC
names(id) <- alldiff$gene_name
## 这里去掉了基因集前缀
hallmark$term <- gsub('HALLMARK_','',hallmark$term)
hallmark.list <- hallmark %>% split(.$term) %>% lapply( "[[", 2)

gsea.re1<-clusterProfiler::GSEA(id,TERM2GENE=hallmark,verbose = T)
g1<-as.data.frame(gsea.re1)
g1<-subset(g1,p.adjust<0.05)
g1<-g1[order(g1$NES,decreasing = T),]
gsea.re2 <- fgseaMultilevel(pathways = hallmark.list,
                            stats = id,
                            minSize=1,
                            maxSize=10000)
colnames(gsea.re2)

gsea.re2$pathway
gsea.re2$padj
gsea.re2$ES
gsea.re2$NES
sum(gsea.re2$padj < 0.05, na.rm = TRUE)
g2 <- gsea.re2[gsea.re2$padj< 0.05,]
g2 <- g2[order(g2$NES,decreasing = T),]
g2$padj
g2$pathway

####ggsea 可视化
names(hallmark.list)
se_hall<-c(head(g2$pathway,10),tail(g1$pathway,10))
sig1<-g2[g2$pathway%in%se_hall,]

hallmark.se<-hallmark.list[sig1$pathway]
df.new <- gseaCurve(id, hallmark.se, sig1)
pal_line<-pal_lancet()(9)#各曲线颜色

p_nes<-ggplot() + 
  geom_gsea(df.new,
            prettyGSEA=T,
            tickcolor='grey30',#竖线颜色
            ticksize=0.1,#线条粗细
            #colour='grey',#色块边框
            linecolor=pal_line,
            linesize=1,lty=1,
            #多个图形时使用，设定每行列图形个数
            nrow = 3,
            ncol = 3
  ) + 
  theme_bw()+
  theme(strip.text = element_text(size =10,face = 'bold.italic'),#增大分面标题字体
        strip.background = element_rect(fill = 'white'))+
  xlab(bquote(italic('Rank')))+ylab(bquote(italic('Enrichment Score')))+
  theme(axis.text.x = element_text(size=10,angle = -30,face = 'plain',hjust = 0.5),
        axis.text.y = element_text(size=10,angle = 0,face = 'plain',vjust = 0.5),
        axis.title.x = element_text(size=12,face = 'plain',hjust = 0.5),
        axis.title.y = element_text(size=12,face = 'plain',vjust = 0.5))
p_nes

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/髓系细胞亚型/M1M2特征得分和GESA分析/Post_MPR_Pre_MPR_gsea.pdf",p_nes,width=7,height=5)

library(ggsci)
col_gsea1<-pal_simpsons()(16)
num2=5
gseaplot2(gsea.re1,geneSetID = rownames(g1)[1:num2],
          title = "",#标题
          color = col_gsea1[1:num2],#颜色
          base_size = 14,#基础大小
          rel_heights = c(1, 0.2, 0.4),#小图相对高度
          subplots = 1:3,#展示小图
          pvalue_table = F,#p值表格
          ES_geom = "line"#line or dot
)


######post_mpr_pre_nmpr----

# deg_all2=FindMarkers(myeloid_cells_object, ident.1 = "Post_MPR", ident.2 = "Pre_MPR")
# write.csv(deg_all2, "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/髓系细胞亚型/M1M2特征得分和GESA分析/Post_MPR_Pre_MPR_deg_all.csv", row.names = T)
# deg_all3=FindMarkers(myeloid_cells_object, ident.1 = "Post_MPR", ident.2 = "Pre_NMPR")
# write.csv(deg_all3, "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/髓系细胞亚型/M1M2特征得分和GESA分析/Post_MPR_Pre_NMPR_deg_all.csv", row.names = T)

deg_all3<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/髓系细胞亚型/M1M2特征得分和GESA分析/Post_MPR_Pre_NMPR_deg_all.csv")
colnames(deg_all3)[1]<- "gene_name"
deg<- deg_all3
head(deg)

txtfile<- "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/髓系细胞亚型/M1M2特征得分和GESA分析/特征得分marker2.txt"

hallmark <- read.gmt(txtfile)
head(hallmark)
View(hallmark)

alldiff <- deg[order(deg$avg_log2FC,decreasing = T),]
id <- alldiff$avg_log2FC
names(id) <- alldiff$gene_name
## 这里去掉了基因集前缀
hallmark$term <- gsub('HALLMARK_','',hallmark$term)
hallmark.list <- hallmark %>% split(.$term) %>% lapply( "[[", 2)

gsea.re1<-clusterProfiler::GSEA(id,TERM2GENE=hallmark,verbose = T)
g1<-as.data.frame(gsea.re1)
g1<-subset(g1,p.adjust<0.05)
g1<-g1[order(g1$NES,decreasing = T),]
gsea.re2 <- fgseaMultilevel(pathways = hallmark.list,
                            stats = id,
                            minSize=1,
                            maxSize=10000)
colnames(gsea.re2)

gsea.re2$pathway
gsea.re2$padj
gsea.re2$ES
gsea.re2$NES
sum(gsea.re2$padj < 0.05, na.rm = TRUE)
g2 <- gsea.re2[gsea.re2$padj< 0.05,]
g2 <- g2[order(g2$NES,decreasing = T),]
g2$padj
g2$pathway

####ggsea 可视化
names(hallmark.list)
se_hall<-c(head(g2$pathway,10),tail(g1$pathway,10))
sig1<-g2[g2$pathway%in%se_hall,]

hallmark.se<-hallmark.list[sig1$pathway]
df.new <- gseaCurve(id, hallmark.se, sig1)
pal_line<-pal_lancet()(9)#各曲线颜色

p_nes<-ggplot() + 
  geom_gsea(df.new,
            prettyGSEA=T,
            tickcolor='grey30',#竖线颜色
            ticksize=0.1,#线条粗细
            #colour='grey',#色块边框
            linecolor=pal_line,
            linesize=1,lty=1,
            #多个图形时使用，设定每行列图形个数
            nrow = 3,
            ncol = 3
  ) + 
  theme_bw()+
  theme(strip.text = element_text(size =10,face = 'bold.italic'),#增大分面标题字体
        strip.background = element_rect(fill = 'white'))+
  xlab(bquote(italic('Rank')))+ylab(bquote(italic('Enrichment Score')))+
  theme(axis.text.x = element_text(size=10,angle = -30,face = 'plain',hjust = 0.5),
        axis.text.y = element_text(size=10,angle = 0,face = 'plain',vjust = 0.5),
        axis.title.x = element_text(size=12,face = 'plain',hjust = 0.5),
        axis.title.y = element_text(size=12,face = 'plain',vjust = 0.5))
p_nes

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/髓系细胞亚型/M1M2特征得分和GESA分析/Post_MPR_Pre_NMPR_gsea.pdf",p_nes,width=6,height=3)

library(ggsci)
col_gsea1<-pal_simpsons()(16)
num2=2
gseaplot2(gsea.re1,geneSetID = rownames(g1)[1:num2],
          title = "",#标题
          color = col_gsea1[1:num2],#颜色
          base_size = 14,#基础大小
          rel_heights = c(1, 0.2, 0.4),#小图相对高度
          subplots = 1:3,#展示小图
          pvalue_table = T,#p值表格
          ES_geom = "line"#line or dot
)




######每个细胞类型GESA----
cell_types <- c("Mast cells", "SPP1+TAMs", "FOLR2+TAMs", 
               "ZBTB16+TAMs", "CXCL10+TAMs", "Neutrophil", 
               "Monocyte-derived dendritic cells", "cDC1", 
               "cDC2", "Migratory cDCs")

# 定义文件保存路径
output_dir <- "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/髓系细胞亚型/M1M2特征得分和GESA分析/Pre_MPR_Pre_NMPR_每种细胞类型得分/"

# 遍历每个细胞类型，计算差异表达基因并保存为 CSV 文件
for (type in cell_type) {
  # 提取当前细胞类型的子集
  cell_type_object <- subset(myeloid_cells_object, subset = myeloid_cells_type == type)
  
  # 检查是否有数据进行分析
  if (nrow(cell_type_object@meta.data) > 0) {
    # 计算差异表达基因
    deg_all <- FindMarkers(cell_type_object, ident.1 = c("Pre_MPR"), ident.2 = "Pre_NMPR")
    
    # 定义输出文件路径
    output_file <- paste0(output_dir, type, "_deg_all.csv")
    
    # 保存结果为 CSV 文件
    write.csv(deg_all, output_file)
    
    cat("保存了", type, "的差异表达结果到文件:", output_file, "\n")
  } else {
    cat(type, "没有足够的细胞数据进行分析。\n")
  }
}

Migratory_cDCs_deg<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/髓系细胞亚型/M1M2特征得分和GESA分析/Pre_MPR_Pre_NMPR_每种细胞类型得分/Migratory cDCs_deg_all.csv")
colnames(Migratory_cDCs_deg)[1]<-"gene_name"
alldiff <- Migratory_cDCs_deg[order(Migratory_cDCs_deg$avg_log2FC,decreasing = T),]
filtered_genes <- alldiff  %>%
  filter(avg_log2FC > 0.25, p_val < 0.05)

id <- filtered_genes$avg_log2FC
names(id) <- filtered_genes$gene_name
Migratory_cDCs_gsea<- fgseaMultilevel(pathways = hallmark.list,
                stats = id,
                minSize=1,
                maxSize=10000)

Migratory_cDCs_gsea<-clusterProfiler::GSEA(id,TERM2GENE=hallmark,verbose = T)
Migratory_cDCs_gsea$pathway
Migratory_cDCs_gsea$pval
Migratory_cDCs_gsea$padj
Migratory_cDCs_gsea$log2err
Migratory_cDCs_gsea$ES
Migratory_cDCs_gsea$NES
Migratory_cDCs_gsea$size
Migratory_cDCs_gsea$leadingEdge
Migratory_cDCs_gsea_result <- data.frame(
  pathway = as.character(Migratory_cDCs_gsea$pathway),
  pval = as.numeric(Migratory_cDCs_gsea$pval),
  padj = as.numeric(Migratory_cDCs_gsea$padj),
  log2err = as.numeric(Migratory_cDCs_gsea$log2err),
  ES = as.numeric(Migratory_cDCs_gsea$ES),
  NES = as.numeric(Migratory_cDCs_gsea$NES),
  size = as.numeric(Migratory_cDCs_gsea$size),
  leadingEdge = sapply(Migratory_cDCs_gsea$leadingEdge, function(x) paste(x, collapse = ",")) # Convert list to comma-separated string
)

# Write to CSV
write.csv(Migratory_cDCs_gsea_result, "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/髓系细胞亚型/M1M2特征得分和GESA分析/Pre_MPR_Pre_NMPR_每种细胞gsea输出结果/Migratory_cDCs_gsea_result.csv", row.names = FALSE)


for (cell_type in cell_types) {
  # Load the corresponding data
  print(cell_type)
  filepath <- paste0("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/髓系细胞亚型/M1M2特征得分和GESA分析/Pre_MPR_Pre_NMPR_每种细胞类型得分/", cell_type, "_deg_all.csv")
  deg_data <- read.csv(filepath)
  colnames(deg_data)[1] <- "gene_name"
  
  # Process data
  alldiff <- deg_data[order(deg_data$avg_log2FC, decreasing = TRUE),]
  # filtered_genes <- alldiff  %>%
  #   filter(avg_log2FC > 0.25, p_val < 0.05)
  
  id <- alldiff$avg_log2FC
  names(id) <- alldiff$gene_name

  # Perform GSEA
  gsea_result <- fgseaMultilevel(pathways = hallmark.list,
                                 stats = id,
                                 minSize = 1,
                                 maxSize = 10000)
  # Save each result to CSV
  gsea_result3 <- data.frame(
    pathway = as.character(gsea_result$pathway),
    pval = as.numeric(gsea_result$pval),
    padj = as.numeric(gsea_result$padj),
    log2err = as.numeric(gsea_result$log2err),
    ES = as.numeric(gsea_result$ES),
    NES = as.numeric(gsea_result$NES),
    size = as.numeric(gsea_result$size),
    leadingEdge = sapply(gsea_result$leadingEdge, function(x) paste(x, collapse = ",")) # Convert list to comma-separated string
  )
  
  results_list[[cell_type]] <- gsea_result2
  
  write.csv(gsea_result2, 
            file = paste0("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/髓系细胞亚型/M1M2特征得分和GESA分析/Pre_MPR_Pre_NMPR_每种细胞gsea输出结果/", cell_type, "_gsea_result.csv"), 
            row.names = FALSE)
}




# Combine all results into one data frame
combined_results <- do.call(rbind, results_list)

combined_results$cell_type<- rownames(combined_results)
combined_results$cell_type <- gsub("\\.[0-9]+$", "", combined_results$cell_type)
# Save combined results to CSV
write.csv(combined_results, "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/髓系细胞亚型/M1M2特征得分和GESA分析/Pre_MPR_Pre_NMPR_每种细胞gsea输出结果/combined_gsea_results.csv", row.names = FALSE)

####绘制每个细胞对应的NES的heatmap----

library(tidyverse)
library(MetBrewer)
library(grid)


df <- read_tsv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/髓系细胞亚型/M1M2特征得分和GESA分析/ggplot2绘制热图细节拉满/data.xls") %>% select(1,where(is.numeric)) %>% select(1:20) %>% 
  pivot_longer(-sample.ID) %>% 
  mutate(value=log(value)) %>% 
  dplyr::rename("median centered<br>log-transformed FPKM"="value")

View(df)

combined_results <- na.omit(combined_results)
head(combined_results)
combined_results1<-combined_results[,c(1,6,9)]
View(combined_results1)
head(combined_results1)
combined_results_long<- combined_results1 %>%
  pivot_longer(-c(cell_type,pathway))


head(combined_results_long)
View(combined_results_long)

p1 <- ggplot() +
  geom_tile(data=combined_results_long,color="black",
            aes(x=cell_type,y=pathway,fill=value)) +
  scale_fill_gradientn(colors=met.brewer("Cassatt1"),na.value = NA)+
  labs(x=NULL,y=NULL)+
  coord_cartesian(clip = "off")+
  theme(axis.text.x =element_text(angle = 90,vjust=0.5,hjust=1,size=8,color="black"),
        axis.text.y=element_text(size=8,color="black"),
        legend.position ="right",
        plot.margin=unit(c(0.5,0.5,0.5,0.5),unit="cm"),
        plot.background = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank())
p1

#####M1M2基因集打分----

m1m2_pws<-read_lines("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/髓系细胞亚型/M1M2特征得分和GESA分析/GSE5099_CLASSICAL_M1_VS_ALTERNATIVE_M2_MACROPHAGE_UP.v2024.1.Hs.gmt")%>%
  lapply(str_split,"\\t")%>%
  unlist(recursive=F)%>%
  lapply(function(x)setNames(list(x[-c(1:2)]),x[1]))%>%
  unlist(recursive=F)

###导入M2
m1m2_pws<-append(m1m2_pws,read_lines("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/髓系细胞亚型/M1M2特征得分和GESA分析/GSE5099_CLASSICAL_M1_VS_ALTERNATIVE_M2_MACROPHAGE_DN.v2024.1.Hs.gmt")%>%
                     lapply(str_split,"\\t")%>%
                     unlist(recursive=F)%>%
                     lapply(function(x)setNames(list(x[-c(1:2)]),x[1]))%>%
                     unlist(recursive=F))

###运用addmodulescore函数，加入富集分数到metadata中
Idents(myeloid_cells_object)<-"myeloid_cells_type"
M1 <-c("IRF1","CXCL11","CXCL9","MARCO","IL1B","CD86","TNF","IL2RA","CXCL10","FCGR1A")
M2<-c("CLEC7A", "GAS7", "CCL18", "CD209", "LIPA", "F13A1", "CTSD", "MS4A4A", "MAF", "CSF1R", "CCL23", "CCL7", "HMOX1", "FN1", "CCL2", "CCL17", "IL27RA", "CXCR4", "PPARG")

######使用M1M2 marker----
M1<-c("CCL3","CCL4","CCL2","IFNG","IRF5","IRF1","IL23A","TNF","KYNU","IL6","CD40","CXCL9","CXCL10","CXCL11","CCL5","CCR7","IL1A","IL1B","CD86","CD80","CD68","NOS2","HLA-DPB1","MARCO","IL2RA","IL15","CD14","FCGR3A","FCGR1A","FCGR1B")									
M2<-c("CD163","ARG2","IL4R","CCL13","CCL17","CCL18","CCL22","CCL24","LYVE1","VEGFA","VEGFB","VEGFC","VEGFD","EGF","CTSA","CTSB","CTSC","CTSD","TGFB1","TGFB2","TGFB3","MMP14","MMP19","MMP9","CLEC7A","WNT7B","FASLG","TNFSF12","TNFSF8","CD276","FN1","IRF4","IDO1","FABP4","CCR2","CD1B","CD1A","ALOX15","CCL26","CHN2")


m1m2_pws<- list(M1,M2)

imm_anno<-AddModuleScore(object=myeloid_cells_object,features=m1m2_pws,name=c("m1up","m1dn"),nbin=12)

head(imm_anno)
###比较通路，当然这个没有显著性

use_colors<-c(
  "Mast cells"="#0000FF", 
  "SPP1+TAMs"="#FFa07a", 
  "FOLR2+TAMs"="#00bfff", 
  "ZBTB16+TAMs"="#7fff00",
  "CXCL10+TAMs"="#8a2be2",
  "Neutrophil"="#00ffff", 
  "Monocyte-derived dendritic cells"="#9bcd9b", 
  "cDC1"="#FF00ff", 
  "cDC2"="#ffd700", 
  "Migratory cDCs"="#8b795e")

cell_type_color<- c("#0000FF","#FFa07a","#00bfff","#7fff00","#8a2be2","#00ffff","#9bcd9b","#FF00ff","#ffd700","#8b795e")
imm_anno$myeloid_cells_type<- factor(imm_anno$myeloid_cells_type,levels=c("Mast cells", "SPP1+TAMs", "FOLR2+TAMs", 
                                                                          "ZBTB16+TAMs", "CXCL10+TAMs", "Neutrophil", 
                                                                          "Monocyte-derived dendritic cells", "cDC1", 
                                                                          "cDC2", "Migratory cDCs"))

VlnPlot(imm_anno,features=c("m1up1","m1dn2"),
        group.by="myeloid_cells_type",pt.size=0,idents=c("Mast cells", "SPP1+TAMs", "FOLR2+TAMs", 
                                                    "ZBTB16+TAMs", "CXCL10+TAMs", "Neutrophil", 
                                                    "Monocyte-derived dendritic cells", "cDC1", 
                                                    "cDC2", "Migratory cDCs"),cols=use_colors)

VlnPlot(imm_anno,features=c("m1up1","m1dn2"),
        group.by="myeloid_cells_type",pt.size=0,idents=c("Mast cells", "SPP1+TAMs", "FOLR2+TAMs", 
                                                         "ZBTB16+TAMs", "CXCL10+TAMs", "Neutrophil", 
                                                         "Monocyte-derived dendritic cells", "cDC1", 
                                                         "cDC2", "Migratory cDCs"),cols=use_colors)


imm_anno_df <- imm_anno@meta.data
names(imm_anno_df)
imm_anno_df <-imm_anno_df %>% mutate(m1_m2_ratio=m1up1/m1dn2 ) 
imm_anno_filter<-imm_anno_df  %>% select(myeloid_cells_type,m1up1,m1dn2,m1_m2_ratio,MPR_Response2)

library(tidyverse)
library(gghalves)
library(rstatix)
library(ggpubr)
library(ggsci)
df<-imm_anno_filter
df_p_val1 <- df %>%
  wilcox_test(m1up1~ myeloid_cells_type) %>%
  adjust_pvalue(p.col = "p", method = "bonferroni") %>%
  add_significance(p.col = "p.adj") %>% 
  add_xy_position(x = "myeloid_cells_type", dodge = 0.8) 

View(df_p_val1)

treatment_color <- c("#4974a4","#4dae47","#f29600")
cell_type_color<- c("#0000FF","#FFa07a","#00bfff","#7fff00","#8a2be2","#00ffff","#9bcd9b","#FF00ff","#ffd700","#8b795e")

df$myeloid_cells_type<- factor(df$myeloid_cells_type,levels=c("Mast cells", "SPP1+TAMs", "FOLR2+TAMs", 
                                                                          "ZBTB16+TAMs", "CXCL10+TAMs", "Neutrophil", 
                                                                          "Monocyte-derived dendritic cells", "cDC1", 
                                                                          "cDC2", "Migratory cDCs"))
df$MPR_Response2<- factor(df$MPR_Response2,levels = c("Pre_NMPR","Pre_MPR","Post_MPR"))
M1_sig<- df %>%
  ggplot(aes(myeloid_cells_type,m1up1)) +
  geom_half_boxplot(fill=cell_type_color,color="black",side="l",errorbar.draw = T,
                    outlier.shape = NA,width=0.8,lwd= 0.7) +
  geom_half_point(color=cell_type_color,side = "r",alpha=0.1,
                  transformation_params = list(height = 0,width = 0.001,seed = 2))+
  #stat_pvalue_manual(df_p_val1,label = "p.adj.signif",label.size=5,hide.ns = F,bracket.size = 0.5,y.position = c(0.5,0.6,0.7))+
  #facet_wrap(.~continent,nrow=1)+
  #scale_y_continuous(limits = c(0,5),breaks = seq(0,2,5))+
  #facet_wrap(~MPR_Response2, ncol = 1)+
  scale_fill_npg()+
  scale_color_npg()+
  labs(title = "M1 signature",y='Signature score',x="")+
  theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),units=,"cm"),
        plot.title = element_text(hjust = 0.5, size = 13),
        panel.background = element_blank(),
        plot.background = element_blank(),
        strip.text = element_text(size=12),
        axis.title.y = element_text(size = 14),
        axis.line = element_line(color = "black",size = 0.4),
        axis.text.y = element_text(color="black",size=12),
        axis.text.x = element_text(margin = margin(t = 2),color="black",size=12,angle = 30,vjust = 1,hjust = 1),
        legend.position = "none",
        panel.spacing = unit(0,"lines"),
        axis.ticks.x=element_line(color="black",size=0.5,lineend = 1),
        axis.line.y=element_line(linetype=1,color="black",size=0.7),
        axis.ticks.y=element_line(color="black",size=0.5,lineend = 1),
        axis.line.x=element_line(linetype=1,color="black",size=0.7))+
  
  coord_cartesian()
M1_sig

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/髓系细胞亚型/M1M2特征得分和GESA分析/M1_signature_boxplot.pdf",M1_sig,width = 6,height = 5)


df_p_val2 <- df %>%
  wilcox_test(m1dn2~ myeloid_cells_type) %>%
  adjust_pvalue(p.col = "p", method = "bonferroni") %>%
  add_significance(p.col = "p.adj") %>% 
  add_xy_position(x = "myeloid_cells_type", dodge = 0.8) 
treatment_color <- c("#4974a4","#4dae47","#f29600")
cell_type_color<- c("#0000FF","#FFa07a","#00bfff","#7fff00","#8a2be2","#00ffff","#9bcd9b","#FF00ff","#ffd700","#8b795e")

df$myeloid_cells_type<- factor(df$myeloid_cells_type,levels=c("Mast cells", "SPP1+TAMs", "FOLR2+TAMs", 
                                                              "ZBTB16+TAMs", "CXCL10+TAMs", "Neutrophil", 
                                                              "Monocyte-derived dendritic cells", "cDC1", 
                                                              "cDC2", "Migratory cDCs"))
df$MPR_Response2<- factor(df$MPR_Response2,levels = c("Pre_NMPR","Pre_MPR","Post_MPR"))
M2_sig<- df %>%
  ggplot(aes(myeloid_cells_type,m1dn2)) +
  geom_half_boxplot(fill=cell_type_color,color="black",side="l",errorbar.draw = T,
                    outlier.shape = NA,width=0.8,lwd= 0.7) +
  geom_half_point(color=cell_type_color,side = "r",alpha=0.1,
                  transformation_params = list(height = 0,width = 0.001,seed = 2))+
  #stat_pvalue_manual(df_p_val1,label = "p.adj.signif",label.size=5,hide.ns = F,bracket.size = 0.5,y.position = c(0.5,0.6,0.7))+
  #facet_wrap(.~continent,nrow=1)+
  #scale_y_continuous(limits = c(0,5),breaks = seq(0,2,5))+
  scale_fill_npg()+
  scale_color_npg()+
  labs(title = "M2 signature",y='Signature score',x="")+
  theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),units=,"cm"),
        plot.title = element_text(hjust = 0.5, size = 13),
        panel.background = element_blank(),
        plot.background = element_blank(),
        strip.text = element_text(size=12),
        axis.title.y = element_text(size = 14),
        axis.line = element_line(color = "black",size = 0.4),
        axis.text.y = element_text(color="black",size=12),
        axis.text.x = element_text(margin = margin(t = 2),color="black",size=12,angle = 30,vjust = 1,hjust = 1),
        legend.position = "none",
        panel.spacing = unit(0,"lines"),
        axis.ticks.x=element_line(color="black",size=0.5,lineend = 1),
        axis.line.y=element_line(linetype=1,color="black",size=0.7),
        axis.ticks.y=element_line(color="black",size=0.5,lineend = 1),
        axis.line.x=element_line(linetype=1,color="black",size=0.7))+
  
  coord_cartesian()
M2_sig

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/髓系细胞亚型/M1M2特征得分和GESA分析/M2_signature_boxplot.pdf",M2_sig,width = 6,height = 5)


####分组分面----
use_colors<-c(
  "Mast cells"="#0000FF", 
  "SPP1+TAMs"="#FFa07a", 
  "FOLR2+TAMs"="#00bfff", 
  "ZBTB16+TAMs"="#7fff00",
  "CXCL10+TAMs"="#8a2be2",
  "Neutrophil"="#00ffff", 
  "Monocyte-derived dendritic cells"="#9bcd9b", 
  "cDC1"="#FF00ff", 
  "cDC2"="#ffd700", 
  "Migratory cDCs"="#8b795e")

M1_sig_facet<- df %>%
  ggplot(aes(myeloid_cells_type,m1up1)) +
  geom_half_boxplot(aes(fill = myeloid_cells_type),color="black",side="l",errorbar.draw = T,
                    outlier.shape = NA,width=0.8,lwd= 0.5) +
  geom_half_point(aes(),side = "r",alpha=0.1,size=0.1,
                  transformation_params = list(height = 0,width = 0.001,seed = 2))+
  #stat_pvalue_manual(df_p_val6,label = "p.adj.signif",label.size=5,hide.ns = F)+
  facet_wrap(.~MPR_Response2,ncol=1)+
  scale_fill_manual(values= use_colors)+
  theme_minimal()+
  #scale_y_continuous(limits = c(0,95),breaks = seq(0,95,20))+
  #scale_fill_npg()+
  scale_color_npg()+
  labs(x=NULL,y=NULL)+
  #labs(title = "Progenitor exhausted CD8 T cell signature",y='Signature score',x="")+
  labs(title = "M1 signature",y='Signature score',x="")+
  theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),units=,"cm"),
        plot.title = element_text(hjust = 0.5, size = 13),
        strip.text = element_text(size = 11),   # 控制分面标题文字样式
        strip.background = element_rect(fill = "#D3D3D3", color = "white"),
        panel.background = element_blank(),
        plot.background = element_blank(),
        axis.title.y = element_text(size = 14),
        axis.line = element_line(color = "black",size = 0.4),
        axis.text.y = element_text(color="black",size=12),
        axis.text.x = element_text(margin = margin(t = 2),color="black",size=12,angle=30,vjust = 1,hjust = 1),
        legend.position = "none",
        panel.spacing = unit(0,"lines"),
        axis.ticks.x=element_line(color="black",size=0.5,lineend = 1),
        axis.line.y=element_line(linetype=1,color="black",size=0.5),
        axis.ticks.y=element_line(color="black",size=0.5,lineend = 1),
        axis.line.x=element_line(linetype=1,color="black",size=0.5))+
  coord_cartesian()

M1_sig_facet

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/髓系细胞亚型/M1M2特征得分和GESA分析/M1_sig_facet_boxplot.pdf",M1_sig_facet,width = 6,height = 8)

M2_sig_facet<- df %>%
  ggplot(aes(myeloid_cells_type,m1dn2)) +
  geom_half_boxplot(aes(fill = myeloid_cells_type),color="black",side="l",errorbar.draw = T,
                    outlier.shape = NA,width=0.8,lwd= 0.5) +
  geom_half_point(aes(),side = "r",alpha=0.1,size=0.1,
                  transformation_params = list(height = 0,width = 0.001,seed = 2))+
  #stat_pvalue_manual(df_p_val6,label = "p.adj.signif",label.size=5,hide.ns = F)+
  facet_wrap(.~MPR_Response2,ncol=1)+
  scale_fill_manual(values= use_colors)+
  theme_minimal()+
  #scale_y_continuous(limits = c(0,95),breaks = seq(0,95,20))+
  #scale_fill_npg()+
  scale_color_npg()+
  labs(x=NULL,y=NULL)+
  #labs(title = "Progenitor exhausted CD8 T cell signature",y='Signature score',x="")+
  labs(title = "M2 signature",y='Signature score',x="")+
  theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),units=,"cm"),
        plot.title = element_text(hjust = 0.5, size = 13),
        strip.text = element_text(size = 11),   # 控制分面标题文字样式
        strip.background = element_rect(fill = "#D3D3D3", color = "white"),
        panel.background = element_blank(),
        plot.background = element_blank(),
        axis.title.y = element_text(size = 14),
        axis.line = element_line(color = "black",size = 0.4),
        axis.text.y = element_text(color="black",size=12),
        axis.text.x = element_text(margin = margin(t = 2),color="black",size=12,angle=30,vjust = 1,hjust = 1),
        legend.position = "none",
        panel.spacing = unit(0,"lines"),
        axis.ticks.x=element_line(color="black",size=0.5,lineend = 1),
        axis.line.y=element_line(linetype=1,color="black",size=0.5),
        axis.ticks.y=element_line(color="black",size=0.5,lineend = 1),
        axis.line.x=element_line(linetype=1,color="black",size=0.5))+
  coord_cartesian()

M2_sig_facet

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/髓系细胞亚型/M1M2特征得分和GESA分析/M2_sig_facet_boxplot.pdf",M2_sig_facet,width = 6,height = 8)


library(MetBrewer)
met.brewer
combined_results_long<- imm_anno_filter %>% select(m1up1,m1dn2,myeloid_cells_type,MPR_Response2) %>%
  pivot_longer(-c(myeloid_cells_type,MPR_Response2))
names(combined_results_long)
head(combined_results_long)
combined_results_long$MPR_Response2<- factor(combined_results_long$MPR_Response2,levels = c("Pre_NMPR","Pre_MPR","Post_MPR"))

p1_heatmap <- ggplot() +
  geom_tile(data=combined_results_long,color="black",
            aes(x=myeloid_cells_type,y=name,fill=value)) +
  scale_fill_gradientn(colors=met.brewer("Cassatt1"),na.value = NA)+
  facet_wrap(~MPR_Response2,ncol=1,strip.position = "right")+
  labs(x=NULL,y="")+
  coord_cartesian(clip = "off")+
  theme(axis.text.x =element_text(size=12,angle = 90,vjust=1,hjust=1,color="black"),
        axis.text.y=element_text(size=12,color="black"),
        legend.position ="right",
        plot.margin=unit(c(0.5,0.5,0.5,0.5),unit="cm"),
        plot.background = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank(),
        strip.text = element_text(size = 10),
              strip.background = element_rect(color="black",fill="grey90"))
p1_heatmap

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/髓系细胞亚型/M1M2特征得分和GESA分析/M1M2_heatmap.pdf",p1_heatmap,width =7,height = 6)

####表达量相关性分析----

myeloid_cells_object=readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/myeloid_cells_object_filter使用.rds")
myeloid_cells_object = UpdateSeuratObject(myeloid_cells_object)
myeloid_cells_object

T_cells_object<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/T_cells_object.rds")
T_cells_object = UpdateSeuratObject(T_cells_object)
T_cells_object

merged_obj <- merge(myeloid_cells_object, y = T_cells_object, add.cell.ids = c("myeloid_cells_object", "T_cells_object"))

View(merged_obj)
merged_obj$cd4_cd8_group
merged_obj$myeloid_cells_type
merged_obj@meta.data$merged_group <- paste(merged_obj@meta.data$cd4_cd8_group, 
                                           merged_obj@meta.data$myeloid_cells_type, sep = "_")

# 查看更新后的 meta.data
head(merged_obj@meta.data)
merged_obj@meta.data$merged_group <- gsub("^NA_", "", merged_obj@meta.data$merged_group)
merged_obj@meta.data$merged_group <- gsub("_NA$", "", merged_obj@meta.data$merged_group)

unique(merged_obj@meta.data$merged_group)

unique(merged_obj$MPR_Response2)

merged_obj_filter<- subset(merged_obj,subset=merged_group %in% c("SPP1+TAMs","FOLR2+TAMs","cDC2",                           
                                                             "Migratory cDCs","cDC1","Monocyte-derived dendritic cells",
                                                             "ZBTB16+TAMs","Mast cells","CXCL10+TAMs",                    
                                                             "Neutrophil","CD8_C1-GNLY","CD8_C4-GZMK","CD8_C3-IL7R"))
merged_obj_filter$merged_group <- gsub("-", "_", merged_obj_filter$merged_group)
unique(merged_obj_filter@meta.data$merged_group)


merged_obj_filter$Treatments
merged_obj_filter_Post <- subset(merged_obj_filter,subset=Treatments %in% c("Post"))
merged_obj_filter$MPR_Response2
merged_obj_filter_Pre_MPR <- subset(merged_obj_filter,subset=MPR_Response2 %in% c("Pre_MPR"))

avg_expression <- AverageExpression(merged_obj_filter_Pre_MPR, group.by = "merged_group", return.seurat = FALSE)

head(avg_expression$RNA)

# 提取基因表达量矩阵
expr_matrix <- as.matrix(avg_expression$RNA)

# 检查矩阵的维度
print(dim(expr_matrix))
rownames(expr_matrix)
colnames(expr_matrix)
# 计算相关性矩阵
cor_matrix <- cor(expr_matrix, method = "pearson")
dim(cor_matrix)
# cor_matrix <- cor(t(expr_matrix), method = "pearson")
# 查看相关性矩阵
print(cor_matrix)

library(pheatmap)

pheatmap(cor_matrix, 
         clustering_distance_rows = "euclidean", 
         clustering_distance_cols = "euclidean", 
         color = colorRampPalette(c("blue", "white", "red"))(50))

####相关性棒棒图----

library(psych)
expr_matrix <- as.matrix(avg_expression$RNA)%>% as.data.frame()
names(expr_matrix)
expr_matrix_GNLY <- expr_matrix[,-c(2,3)] %>% as.data.frame()
colnames(expr_matrix_GNLY)<-c("CD8_C1_GNLY","cDC1",                            
                              "cDC2","CXCL10+TAMs",                     
                              "FOLR2+TAMs","Mast cells",                      
                              "Migratory cDCs","Monocyte-derived dendritic cells",
                              "Neutrophil","SPP1+TAMs",                       
                              "ZBTB16+TAMs")
View(expr_matrix_GNLY)
CD8_C1_GNLY <- expr_matrix_GNLY$CD8_C1_GNLY
myeloid_cells<- expr_matrix_GNLY[,-c(1)]
pearson <- corr.test(CD8_C1_GNLY, myeloid_cells, method = 'pearson', adjust = 'bonferroni')
r <- data.frame(pearson$r)  #pearson 相关系数矩阵
p <- data.frame(pearson$p.adj)  #p 值矩阵

#p <- data.frame(pearson$p)
#结果整理以便于作图
# r$CD8_C1_GNLY <- rownames(r)
# p$CD8_C1_GNLY <- rownames(p)
# r <- melt(r, id = 'CD8_C1_GNLY')
# p <- melt(p, id = 'CD8_C1_GNLY')
dim(r)
dim(p)
pearson_result <- rbind(r, p)
#pearson <- cbind(r, p$value)
View(pearson_result)
pearson_result<-t(pearson_result)%>% as.data.frame()

pearson_result$myeloid_cells<- rownames(pearson_result)

colnames(pearson_result) <- c('pearson_correlation', 'p.adjusted.value',"myeloid_cells")
pearson_filtered <- subset(pearson_result, p.adjusted.value< 0.05)

pearson_filtered$myeloid_cells <- factor(pearson_filtered$myeloid_cells, levels = c("Mast.cells", "SPP1.TAMs", "FOLR2.TAMs", 
                                                                                    "ZBTB16.TAMs", "CXCL10.TAMs", "Neutrophil", 
                                                                                    "Monocyte.derived.dendritic.cells", "cDC1", 
                                                                                    "cDC2", "Migratory.cDCs"))
head(pearson_filtered$myeloid_cells)  #整理好的环境变量和物种丰度的 pearson 相关性统计表


col_set<- c("Mast.cells"="#0000FF",
                    "SPP1.TAMs"="#FFa07a",
                    "FOLR2.TAMs"="#00bfff",
                    "ZBTB16.TAMs"="#7fff00",
                    "CXCL10.TAMs"="#8a2be2",
                    "Neutrophil"="#00ffff",
                    "Monocyte.derived.dendritic.cells"="#9bcd9b",
                    "cDC1" ="#FF00ff",
                    "cDC2"="#ffd700",
                    "Migratory.cDCs"="#8b795e")

p4 <- ggplot(pearson_filtered,
             aes(x = myeloid_cells, y = pearson_correlation)) +
  # 棒棒图的连线
  geom_segment(aes(x = myeloid_cells, xend = myeloid_cells, y = 0, yend = pearson_correlation),  
               linetype = "solid", # 实线
               size = 1, # 连线的粗细 
               color = "gray40" # 连线的颜色
  ) + 
  # y轴0刻度的水平线
  geom_hline( 
    yintercept = 0,  # 水平线位置
    linetype = "dashed",  # 虚线
    size = 1, # 连线的粗细 
    colour="gray40" # 连线的颜色
  ) +
  # 绘制点
  geom_point(aes(color = myeloid_cells),  
             color = col_set,
             size = 14) +   
  # 添加相关性R值标签
  geom_text(aes(label = round(pearson_correlation, 2)), 
            #color = ifelse(round(spearman_correlation, 2) != 0.96, "black", 'red'), 
            size = 4) +
  # 添加相关性P值标签,这一步最精华的一点是，根据r值正负调整水平移动的位置。
  geom_text(aes(label = paste("p=",round(p.adjusted.value, 2),sep="")),
            hjust =-1,
            #ifelse(spearman_filtered$spearman_correlation >= 0, 1.5, -0.5),
            vjust = -0.2,
            angle = 90,
            fontface = 'italic',
            #color = ifelse(data$p != 'p<0.001', "black", 'red'), 
            size = 4) +
  # y轴刻度设置
  scale_y_continuous( #设置y轴
    limits = c(0, 1.0),
    breaks = c(0, 0.4, 0.8),
    labels = c(0, 0.4, 0.8)) +
  # 坐标title和图形title设置
  labs(y = "Pearson correlation coefficient",
       title = paste0("Correlations between CD8-C1-GNLY expression \n with Myeloid subsets")) +
  theme_classic() +
  theme(axis.line = element_line(size = 0.7),
        plot.title = element_text(size = 13, hjust = 0.5),
        axis.text.x = element_text(size = 13, angle = 45, hjust=1,colour = "black"), 
        axis.text.y = element_text(size = 13),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 13),
        plot.margin = unit(c(1, 1, 1, 1), "cm")
  )

p4

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/髓系细胞亚型/相关性分析/GNLY和髓系亚型相关性棒棒图.pdf",p4,width = 6,height =6)


#####"CD8_C4-GZMK","CD8_C3-IL7R" 棒棒图----
library(psych)
expr_matrix <- as.matrix(avg_expression$RNA)%>% as.data.frame()
names(expr_matrix)
expr_matrix_GZMK <- expr_matrix[,-c(1,2)] %>% as.data.frame()
colnames(expr_matrix_GZMK)<-c("CD8_C1_GZMK","cDC1",                            
                              "cDC2","CXCL10+TAMs",                     
                              "FOLR2+TAMs","Mast cells",                      
                              "Migratory cDCs","Monocyte-derived dendritic cells",
                              "Neutrophil","SPP1+TAMs",                       
                              "ZBTB16+TAMs")
View(expr_matrix_GZMK)
CD8_C1_GZMK <- expr_matrix_GZMK$CD8_C1_GZMK
myeloid_cells<- expr_matrix_GZMK[,-c(1)]
pearson <- corr.test(CD8_C1_GZMK, myeloid_cells, method = 'pearson', adjust = 'bonferroni')
r <- data.frame(pearson$r)  #pearson 相关系数矩阵
p <- data.frame(pearson$p.adj)  #p 值矩阵

dim(r)
dim(p)
pearson_result <- rbind(r, p)
#pearson <- cbind(r, p$value)
View(pearson_result)
pearson_result<-t(pearson_result)%>% as.data.frame()

pearson_result$myeloid_cells<- rownames(pearson_result)

colnames(pearson_result) <- c('pearson_correlation', 'p.adjusted.value',"myeloid_cells")
pearson_filtered <- subset(pearson_result, p.adjusted.value< 0.05)

pearson_filtered$myeloid_cells <- factor(pearson_filtered$myeloid_cells, levels = c("Mast.cells", "SPP1.TAMs", "FOLR2.TAMs", 
                                                                                    "ZBTB16.TAMs", "CXCL10.TAMs", "Neutrophil", 
                                                                                    "Monocyte.derived.dendritic.cells", "cDC1", 
                                                                                    "cDC2", "Migratory.cDCs"))
head(pearson_filtered$myeloid_cells)  #整理好的环境变量和物种丰度的 pearson 相关性统计表


col_set<- c("Mast.cells"="#0000FF",
            "SPP1.TAMs"="#FFa07a",
            "FOLR2.TAMs"="#00bfff",
            "ZBTB16.TAMs"="#7fff00",
            "CXCL10.TAMs"="#8a2be2",
            "Neutrophil"="#00ffff",
            "Monocyte.derived.dendritic.cells"="#9bcd9b",
            "cDC1" ="#FF00ff",
            "cDC2"="#ffd700",
            "Migratory.cDCs"="#8b795e")

p4 <- ggplot(pearson_filtered,
             aes(x = myeloid_cells, y = pearson_correlation)) +
  # 棒棒图的连线
  geom_segment(aes(x = myeloid_cells, xend = myeloid_cells, y = 0, yend = pearson_correlation),  
               linetype = "solid", # 实线
               size = 1, # 连线的粗细 
               color = "gray40" # 连线的颜色
  ) + 
  # y轴0刻度的水平线
  geom_hline( 
    yintercept = 0,  # 水平线位置
    linetype = "dashed",  # 虚线
    size = 1, # 连线的粗细 
    colour="gray40" # 连线的颜色
  ) +
  # 绘制点
  geom_point(aes(color = myeloid_cells),  
             color = col_set,
             size = 14) +   
  # 添加相关性R值标签
  geom_text(aes(label = round(pearson_correlation, 2)), 
            #color = ifelse(round(spearman_correlation, 2) != 0.96, "black", 'red'), 
            size = 4) +
  # 添加相关性P值标签,这一步最精华的一点是，根据r值正负调整水平移动的位置。
  geom_text(aes(label = paste("p=",round(p.adjusted.value, 2),sep="")),
            hjust =-1,
            #ifelse(spearman_filtered$spearman_correlation >= 0, 1.5, -0.5),
            vjust = -0.2,
            angle = 90,
            fontface = 'italic',
            #color = ifelse(data$p != 'p<0.001', "black", 'red'), 
            size = 4) +
  # y轴刻度设置
  scale_y_continuous( #设置y轴
    limits = c(0, 1.0),
    breaks = c(0, 0.4, 0.8),
    labels = c(0, 0.4, 0.8)) +
  # 坐标title和图形title设置
  labs(y = "Pearson correlation coefficient",
       title = paste0("Correlations between CD8-C1-GZMK expression \n with Myeloid subsets")) +
  theme_classic() +
  theme(axis.line = element_line(size = 0.7),
        plot.title = element_text(size = 13, hjust = 0.5),
        axis.text.x = element_text(size = 13, angle = 45, hjust=1,colour = "black"), 
        axis.text.y = element_text(size = 13),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 13),
        plot.margin = unit(c(1, 1, 1, 1), "cm")
  )

p4

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/髓系细胞亚型/相关性分析/GZMK和髓系亚型相关性棒棒图.pdf",p4,width = 6,height =6)


#####CD8_C3-IL7R" 棒棒图----
library(psych)
expr_matrix <- as.matrix(avg_expression$RNA)%>% as.data.frame()
names(expr_matrix)
expr_matrix_IL7R <- expr_matrix[,-c(1,3)] %>% as.data.frame()
colnames(expr_matrix_IL7R)<-c("CD8_C1_IL7R","cDC1",                            
                              "cDC2","CXCL10+TAMs",                     
                              "FOLR2+TAMs","Mast cells",                      
                              "Migratory cDCs","Monocyte-derived dendritic cells",
                              "Neutrophil","SPP1+TAMs",                       
                              "ZBTB16+TAMs")
View(expr_matrix_IL7R)
CD8_C1_IL7R <- expr_matrix_IL7R$CD8_C1_IL7R
myeloid_cells<- expr_matrix_IL7R[,-c(1)]
pearson <- corr.test(CD8_C1_IL7R, myeloid_cells, method = 'pearson', adjust = 'bonferroni')
r <- data.frame(pearson$r)  #pearson 相关系数矩阵
p <- data.frame(pearson$p.adj)  #p 值矩阵

dim(r)
dim(p)
pearson_result <- rbind(r, p)
#pearson <- cbind(r, p$value)
View(pearson_result)
pearson_result<-t(pearson_result)%>% as.data.frame()

pearson_result$myeloid_cells<- rownames(pearson_result)

colnames(pearson_result) <- c('pearson_correlation', 'p.adjusted.value',"myeloid_cells")
pearson_filtered <- subset(pearson_result, p.adjusted.value< 0.05)

pearson_filtered$myeloid_cells <- factor(pearson_filtered$myeloid_cells, levels = c("Mast.cells", "SPP1.TAMs", "FOLR2.TAMs", 
                                                                                    "ZBTB16.TAMs", "CXCL10.TAMs", "Neutrophil", 
                                                                                    "Monocyte.derived.dendritic.cells", "cDC1", 
                                                                                    "cDC2", "Migratory.cDCs"))
head(pearson_filtered$myeloid_cells)  #整理好的环境变量和物种丰度的 pearson 相关性统计表


col_set<- c("Mast.cells"="#0000FF",
            "SPP1.TAMs"="#FFa07a",
            "FOLR2.TAMs"="#00bfff",
            "ZBTB16.TAMs"="#7fff00",
            "CXCL10.TAMs"="#8a2be2",
            "Neutrophil"="#00ffff",
            "Monocyte.derived.dendritic.cells"="#9bcd9b",
            "cDC1" ="#FF00ff",
            "cDC2"="#ffd700",
            "Migratory.cDCs"="#8b795e")

p4 <- ggplot(pearson_filtered,
             aes(x = myeloid_cells, y = pearson_correlation)) +
  # 棒棒图的连线
  geom_segment(aes(x = myeloid_cells, xend = myeloid_cells, y = 0, yend = pearson_correlation),  
               linetype = "solid", # 实线
               size = 1, # 连线的粗细 
               color = "gray40" # 连线的颜色
  ) + 
  # y轴0刻度的水平线
  geom_hline( 
    yintercept = 0,  # 水平线位置
    linetype = "dashed",  # 虚线
    size = 1, # 连线的粗细 
    colour="gray40" # 连线的颜色
  ) +
  # 绘制点
  geom_point(aes(color = myeloid_cells),  
             color = col_set,
             size = 14) +   
  # 添加相关性R值标签
  geom_text(aes(label = round(pearson_correlation, 2)), 
            #color = ifelse(round(spearman_correlation, 2) != 0.96, "black", 'red'), 
            size = 4) +
  # 添加相关性P值标签,这一步最精华的一点是，根据r值正负调整水平移动的位置。
  geom_text(aes(label = paste("p=",round(p.adjusted.value, 2),sep="")),
            hjust =-1,
            #ifelse(spearman_filtered$spearman_correlation >= 0, 1.5, -0.5),
            vjust = -0.2,
            angle = 90,
            fontface = 'italic',
            #color = ifelse(data$p != 'p<0.001', "black", 'red'), 
            size = 4) +
  # y轴刻度设置
  scale_y_continuous( #设置y轴
    limits = c(0, 1.0),
    breaks = c(0, 0.4, 0.8),
    labels = c(0, 0.4, 0.8)) +
  # 坐标title和图形title设置
  labs(y = "Pearson correlation coefficient",
       title = paste0("Correlations between CD8-C1-IL7R expression \n with Myeloid subsets")) +
  theme_classic() +
  theme(axis.line = element_line(size = 0.7),
        plot.title = element_text(size = 13, hjust = 0.5),
        axis.text.x = element_text(size = 13, angle = 45, hjust=1,colour = "black"), 
        axis.text.y = element_text(size = 13),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 13),
        plot.margin = unit(c(1, 1, 1, 1), "cm")
  )

p4

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/髓系细胞亚型/相关性分析/IL7R和髓系亚型相关性棒棒图.pdf",p4,width = 6,height =6)



####髓系和t细胞通讯分析----
library(CellChat)
library(tidyverse)
library(Seurat)

myeloid_cells_object=readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/myeloid_cells_object_filter使用.rds")
myeloid_cells_object = UpdateSeuratObject(myeloid_cells_object)
myeloid_cells_object

T_cells_object<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/T_cells_object.rds")
T_cells_object = UpdateSeuratObject(T_cells_object)
T_cells_object

merged_obj <- merge(myeloid_cells_object, y = T_cells_object, add.cell.ids = c("myeloid_cells_object", "T_cells_object"))

View(merged_obj)
merged_obj$cd4_cd8_group
merged_obj$myeloid_cells_type
merged_obj@meta.data$merged_group <- paste(merged_obj@meta.data$cd4_cd8_group, 
                                           merged_obj@meta.data$myeloid_cells_type, sep = "_")

# 查看更新后的 meta.data
head(merged_obj@meta.data)
merged_obj@meta.data$merged_group <- gsub("^NA_", "", merged_obj@meta.data$merged_group)
merged_obj@meta.data$merged_group <- gsub("_NA$", "", merged_obj@meta.data$merged_group)

unique(merged_obj@meta.data$merged_group)

unique(merged_obj$MPR_Response2)

merged_obj


####分组比较----


KS_cellchat <- function(input_obj,
                        assay= NULL,
                        group.by = NULL,
                        workers,
                        species=c('human','mouse'),
                        CellChatDB.use=NULL,#a character vector, which is a subset of c("Secreted Signaling","ECM-Receptor","Cell-Cell Contact","Non-protein Signaling")
                        PPIuse=F,
                        type="triMean",# c("triMean", "truncatedMean", "thresholdedMean", "median")
                        min.cells = 10
){
  
  
  
  cellchat.obj = createCellChat(input_obj, assay = assay, group.by = group.by)
  
  if(species=='human'){
    
    CellChatDB <- CellChatDB.human
    ppi = PPI.human
  }
  
  if(species =="mouse"){
    
    CellChatDB <- CellChatDB.mouse
    ppi = PPI.mouse
  }
  
  
  if(is.null(CellChatDB.use)){
    
    cellchat.obj@DB <- CellChatDB
    
  }else{
    
    CellChatDB <- subsetDB(CellChatDB, search = CellChatDB.use, key = "annotation")
    cellchat.obj@DB <- CellChatDB
  }
  
  cellchat.obj <- subsetData(cellchat.obj) 
  future::plan("multisession", workers = workers) # 并行计算
  cellchat.obj <- identifyOverExpressedGenes(cellchat.obj)
  cellchat.obj <- identifyOverExpressedInteractions(cellchat.obj)
  
  if(PPIuse==F){
    
    cellchat.obj <- computeCommunProb(cellchat.obj, type = type)
    
  }else{
    
    cellchat.obj <- projectData(cellchat.obj, ppi)
    cellchat.obj <- computeCommunProb(cellchat.obj, raw.use=F, type = type)
  }
  
  
  cellchat.obj <- filterCommunication(cellchat.obj, min.cells = min.cells)
  cellchat.obj <- computeCommunProbPathway(cellchat.obj)
  
  #互作网络整合,可以设置soure和target，我们这里默认了，就是全部的
  cellchat.obj <- aggregateNet(cellchat.obj)
  
  return(cellchat.obj)
  
}

merged_obj_filter<- subset(merged_obj,subset=merged_group %in% c("SPP1+TAMs","FOLR2+TAMs","cDC2",                           
                                                                 "Migratory cDCs","cDC1","Monocyte-derived dendritic cells",
                                                                 "ZBTB16+TAMs","Mast cells","CXCL10+TAMs",                    
                                                                 "Neutrophil","CD8_C1-GNLY","CD8_C4-GZMK","CD8_C3-IL7R"))

merged_obj_filter<- myeloid_cells_object
pre_nmpr_myeloid_t_chat<- subset(merged_obj_filter, subset=MPR_Response2=="Pre_NMPR")
unique(pre_nmpr_myeloid_t_chat$MPR_Response2)
pre_mpr_myeloid_t_chat<- subset(merged_obj_filter, subset=MPR_Response2=="Pre_MPR")
unique(pre_mpr_myeloid_t_chat$MPR_Response2)
post_mpr_myeloid_t_chat<- subset(merged_obj_filter, subset=MPR_Response2=="Post_MPR")
unique(post_mpr_myeloid_t_chat$MPR_Response2)
unique(post_mpr_myeloid_t_chat$CRNCR_Response2)
merged_obj_filter$group_microbe2

myeloid_t_bacteria_yes_chat<- subset(merged_obj_filter, subset=group_microbe2=="Bacteria+")

myeloid_t_bacteria_yes_pre_nmpr_myeloid_t_chat<- subset(myeloid_t_bacteria_yes_chat, subset=MPR_Response2=="Pre_NMPR")
unique(myeloid_t_bacteria_yes_pre_nmpr_myeloid_t_chat$MPR_Response2)
myeloid_t_bacteria_yes_pre_mpr_myeloid_t_chat<- subset(myeloid_t_bacteria_yes_chat, subset=MPR_Response2=="Pre_MPR")
unique(myeloid_t_bacteria_yes_pre_mpr_myeloid_t_chat$MPR_Response2)
myeloid_t_bacteria_yes_post_mpr_myeloid_t_chat<- subset(myeloid_t_bacteria_yes_chat, subset=MPR_Response2=="Post_MPR")
unique(myeloid_t_bacteria_yes_post_mpr_myeloid_t_chat$MPR_Response2)

myeloid_t_bacteria_yes_pre_nmpr_myeloid_t_chat<- KS_cellchat(myeloid_t_bacteria_yes_pre_nmpr_myeloid_t_chat, assay = 'RNA', group.by = "merged_group",
                                                            workers=2, species='human')

myeloid_t_bacteria_yes_pre_mpr_myeloid_t_chat<- KS_cellchat(myeloid_t_bacteria_yes_pre_mpr_myeloid_t_chat, assay = 'RNA', group.by = "merged_group",
                                                           workers=2, species='human')

myeloid_t_bacteria_yes_post_mpr_myeloid_t_chat<- KS_cellchat(myeloid_t_bacteria_yes_post_mpr_myeloid_t_chat, assay = 'RNA', group.by = "merged_group",
                                                            workers=2, species='human')


saveRDS(myeloid_t_bacteria_yes_pre_mpr_myeloid_t_chat,
        file = "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/髓系细胞亚型/细胞通讯分析/myeloid_t_bacteria_yes_pre_mpr_myeloid_t_chat.rds")

saveRDS(myeloid_t_bacteria_yes_pre_mpr_myeloid_t_chat,
        file = "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/髓系细胞亚型/细胞通讯分析/myeloid_t_bacteria_yes_pre_mpr_myeloid_t_chat.rds")

saveRDS(myeloid_t_bacteria_yes_post_mpr_myeloid_t_chat,
        file = "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/髓系细胞亚型/细胞通讯分析/myeloid_t_bacteria_yes_post_mpr_myeloid_t_chat.rds")



myeloid_t_bacteria_no_chat <- subset(merged_obj_filter, subset=group_microbe2=="Bacteria-")

myeloid_t_bacteria_no_pre_nmpr_myeloid_t_chat<- subset(myeloid_t_bacteria_no_chat, subset=MPR_Response2=="Pre_NMPR")
unique(myeloid_t_bacteria_no_pre_nmpr_myeloid_t_chat$MPR_Response2)
myeloid_t_bacteria_no_pre_mpr_myeloid_t_chat<- subset(myeloid_t_bacteria_no_chat, subset=MPR_Response2=="Pre_MPR")
unique(myeloid_t_bacteria_no_pre_mpr_myeloid_t_chat$MPR_Response2)
myeloid_t_bacteria_no_post_mpr_myeloid_t_chat<- subset(myeloid_t_bacteria_no_chat, subset=MPR_Response2=="Post_MPR")
unique(myeloid_t_bacteria_no_post_mpr_myeloid_t_chat$MPR_Response2)

myeloid_t_bacteria_no_pre_nmpr_myeloid_t_chat<- KS_cellchat(myeloid_t_bacteria_no_pre_nmpr_myeloid_t_chat, assay = 'RNA', group.by = "merged_group",
                                 workers=2, species='human')

myeloid_t_bacteria_no_pre_mpr_myeloid_t_chat<- KS_cellchat(myeloid_t_bacteria_no_pre_mpr_myeloid_t_chat, assay = 'RNA', group.by = "merged_group",
                                                            workers=2, species='human')

myeloid_t_bacteria_no_post_mpr_myeloid_t_chat<- KS_cellchat(myeloid_t_bacteria_no_post_mpr_myeloid_t_chat, assay = 'RNA', group.by = "merged_group",
                                                           workers=2, species='human')


saveRDS(myeloid_t_bacteria_no_pre_mpr_myeloid_t_chat,
        file = "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/髓系细胞亚型/细胞通讯分析/myeloid_t_bacteria_no_pre_mpr_myeloid_t_chat.rds")

saveRDS(myeloid_t_bacteria_no_pre_mpr_myeloid_t_chat,
        file = "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/髓系细胞亚型/细胞通讯分析/myeloid_t_bacteria_no_pre_mpr_myeloid_t_chat.rds")

saveRDS(myeloid_t_bacteria_no_post_mpr_myeloid_t_chat,
        file = "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/髓系细胞亚型/细胞通讯分析/myeloid_t_bacteria_no_post_mpr_myeloid_t_chat.rds")



post_mpr_myeloid_t_chat_cr<- subset(post_mpr_myeloid_t_chat, subset=CRNCR_Response2=="Post_CR")
post_mpr_myeloid_t_chat_ncr<- subset(post_mpr_myeloid_t_chat, subset=CRNCR_Response2=="Post_NCR")


pre_nmpr.cellchat <- KS_cellchat(pre_nmpr_myeloid_t_chat, assay = 'RNA', group.by = "merged_group",
                            workers=2, species='human')

pre_nmpr.cellchat <- KS_cellchat(pre_nmpr_myeloid_t_chat, assay = 'RNA', group.by = "myeloid_cells_type",
                                 workers=2, species='human')

saveRDS(pre_nmpr.cellchat, file = "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/髓系细胞亚型/细胞通讯分析/pre_nmpr.cellchat_IL7R.rds")

options(future.globals.maxSize = 1e+09) 
pre_mpr.cellchat <- KS_cellchat(pre_mpr_myeloid_t_chat, assay = 'RNA', group.by = "merged_group",
                                 workers=2, species='human')

pre_mpr.cellchat <- KS_cellchat(pre_mpr_myeloid_t_chat, assay = 'RNA', group.by = "myeloid_cells_type",
                                workers=2, species='human')

saveRDS(pre_mpr.cellchat , file = "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/髓系细胞亚型/细胞通讯分析/pre_mpr.cellchat_IL7R.rds")


post_mpr.cellchat <- KS_cellchat(post_mpr_myeloid_t_chat, assay = 'RNA', group.by = "merged_group",
                                 workers=2, species='human')

post_mpr.cellchat <- KS_cellchat(post_mpr_myeloid_t_chat, assay = 'RNA', group.by = "myeloid_cells_type",
                                 workers=2, species='human')

saveRDS(post_mpr.cellchat, file = "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/髓系细胞亚型/细胞通讯分析/post_mpr.cellchat_IL7R.rds")


##################################################

myeloid_t_bacteria_yes_chat.cellchat <- KS_cellchat(myeloid_t_bacteria_yes_chat, assay = 'RNA', group.by = "merged_group",
                                                    workers=2, species='human')

saveRDS(myeloid_t_bacteria_yes_chat.cellchat, file = "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/髓系细胞亚型/细胞通讯分析/myeloid_t_bacteria_yes.cellchat.rds")


options(future.globals.maxSize = 4e+09)
myeloid_t_bacteria_no_chat.cellchat <- KS_cellchat(myeloid_t_bacteria_no_chat, assay = 'RNA', group.by = "merged_group",
                                                   workers=2, species='human')

saveRDS(myeloid_t_bacteria_no_chat.cellchat, file = "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/髓系细胞亚型/细胞通讯分析/myeloid_t_bacteria_no.cellchat.rds")

unique(myeloid_t_bacteria_no_chat.cellchat@meta$merged_group)

unique(myeloid_t_bacteria_yes_chat.cellchat@meta$merged_group)



##############cr 和ncr
post_mpr_cr.cellchat <- KS_cellchat(post_mpr_myeloid_t_chat_cr, assay = 'RNA', group.by = "merged_group",
                                 workers=2, species='human')

saveRDS(post_mpr_cr.cellchat, file = "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/髓系细胞亚型/细胞通讯分析/post_mpr_cr.cellchat.rds")

post_mpr_ncr.cellchat_ncr <- KS_cellchat(post_mpr_myeloid_t_chat, assay = 'RNA', group.by = "merged_group",
                                 workers=2, species='human')

saveRDS(post_mpr_ncr.cellchat_ncr, file = "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/髓系细胞亚型/细胞通讯分析/post_mpr_ncr.cellchat_ncr.rds")


####cellchat多组比较分析----

pre_nmpr.cellchat<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/髓系细胞亚型/细胞通讯分析/pre_nmpr.cellchat_IL7R_GNLY_GZMK.rds")

pre_mpr.cellchat<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/髓系细胞亚型/细胞通讯分析/pre_mpr.cellchat_IL7R_GNLY_GZMK.rds")

post_mpr.cellchat<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/髓系细胞亚型/细胞通讯分析/post_mpr.cellchat_IL7R_GNLY_GZMK.rds")


myeloid_t_bacteria_yes.cellchat<-readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/髓系细胞亚型/细胞通讯分析/myeloid_t_bacteria_yes.cellchat.rds")


myeloid_t_bacteria_no.cellchat<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/髓系细胞亚型/细胞通讯分析/myeloid_t_bacteria_no.cellchat.rds")


object.list <- list(bacteria_yes= myeloid_t_bacteria_yes.cellchat, bacteria_no= myeloid_t_bacteria_no.cellchat)

object.list <- list(Pre_NMPR= pre_nmpr.cellchat, Pre_MPR= pre_mpr.cellchat,Post_MPR=post_mpr.cellchat)

object.list <- list(Pre_NMPR= myeloid_t_bacteria_yes_pre_nmpr_myeloid_t_chat, 
                    Pre_MPR= myeloid_t_bacteria_yes_pre_mpr_myeloid_t_chat,
                    Post_MPR=myeloid_t_bacteria_yes_post_mpr_myeloid_t_chat)

object.list <- list(Pre_NMPR= myeloid_t_bacteria_no_pre_nmpr_myeloid_t_chat, 
                    Pre_MPR= myeloid_t_bacteria_no_pre_mpr_myeloid_t_chat,
                    Post_MPR=myeloid_t_bacteria_no_post_mpr_myeloid_t_chat)

cellchat <- mergeCellChat(object.list, add.names = names(object.list))

# object.list2 <- list(Post_MPR_CR=post_mpr_cr.cellchat,Post_MPR_NCR=post_mpr_ncr.cellchat_ncr)
# cellchat2 <- mergeCellChat(object.list2, add.names = names(object.list2),cell.prefix = TRUE)

gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
p1<-gg1 + gg2


gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2,3))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2,3), measure = "weight")
p1<-gg1 + gg2

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/髓系细胞亚型/细胞通讯分析/分组比较柱状图_IL7R_GNLY_GZMK.pdf",p1,width=6,height=4)


gg1 <- netVisual_heatmap(cellchat)
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
pdf(file = "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/髓系细胞亚型/细胞通讯分析/分组比较热图_IL7R_GNLY_GZMK.pdf", width = 10, height = 6)
gg1 + gg2
dev.off()

mat<-cellchat@net$bacteria_yes$weight

bacteria_yes<-pheatmap(as.matrix(mat),
                   cluster_rows=F, 
                   cluster_cols=F, 
                   show_colnames=TRUE, 
                   main="Heatmap of Interaction Strength")

mat<-cellchat@net$bacteria_no$weight

bacteria_no<-pheatmap(as.matrix(mat),
                       cluster_rows=F, 
                       cluster_cols=F, 
                       show_colnames=TRUE, 
                       main="Heatmap of Interaction Strength")
bacteria_yes+bacteria_no

num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  # 计算网络中心性得分
  object.list[[i]] <- netAnalysis_computeCentrality(object.list[[i]], slot.name = "netP") 
  # 然后生成散点图
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}
patchwork::wrap_plots(plots = gg)

p_combined <- patchwork::wrap_plots(plots = gg)
ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/髓系细胞亚型/细胞通讯分析/combined_plots_bacteria_yes_no.pdf", p_combined, width = 10, height = 4) 

?rankNet
treatment_color <- c("#4974a4","#4dae47","#f29600")
gg1 <- rankNet(cellchat, mode = "comparison",comparison = c(1, 2), stacked = T,font.size = 10, sources.use = "CD8_C3-IL7R",targets.use = c("SPP1+TAMs","FOLR2+TAMs","cDC2",                           
                                                                                                                                           "Migratory cDCs","cDC1","Monocyte-derived dendritic cells",
                                                                                                                                           "ZBTB16+TAMs","Mast cells","CXCL10+TAMs",                    
                                                                                                                                           "Neutrophil"),do.stat = TRUE, color.use = c("#4974a4","#4dae47"))

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/髓系细胞亚型/细胞通讯分析/有无细菌信息柱状图.pdf",gg1,width=5,height=7)

gg2 <- rankNet(cellchat, mode = "comparison",comparison = c(1, 3), stacked = T, font.size = 10,sources.use = "CD8_C3-IL7R",targets.use = c("SPP1+TAMs","FOLR2+TAMs","cDC2",                           
                                                                                                                                           "Migratory cDCs","cDC1","Monocyte-derived dendritic cells",
                                                                                                                                           "ZBTB16+TAMs","Mast cells","CXCL10+TAMs",                    
                                                                                                                                           "Neutrophil"),do.stat = TRUE,color.use = c("#4974a4","#f29600"))
gg3 <- rankNet(cellchat, mode = "comparison",comparison = c(2, 3), stacked = T, font.size = 10,sources.use = "CD8_C3-IL7R",targets.use = c("SPP1+TAMs","FOLR2+TAMs","cDC2",                           
                                                                                                                                           "Migratory cDCs","cDC1","Monocyte-derived dendritic cells",
                                                                                                                                           "ZBTB16+TAMs","Mast cells","CXCL10+TAMs",                    
                                                                                                                                           "Neutrophil"),do.stat = TRUE,color.use = c("#4dae47","#f29600"))
gg1+gg2+gg3


gg4 <- rankNet(cellchat, mode = "comparison",comparison = c(1, 2), stacked = T,font.size = 10, sources.use = "CD8_C1-GNLY",targets.use = c("SPP1+TAMs","FOLR2+TAMs","cDC2",                           
                                                                                                                                           "Migratory cDCs","cDC1","Monocyte-derived dendritic cells",
                                                                                                                                           "ZBTB16+TAMs","Mast cells","CXCL10+TAMs",                    
                                                                                                                                           "Neutrophil"),do.stat = TRUE, color.use = c("#4974a4","#4dae47"))

gg5 <- rankNet(cellchat, mode = "comparison",comparison = c(1, 3), stacked = T, font.size = 10,sources.use = "CD8_C1-GNLY",targets.use = c("SPP1+TAMs","FOLR2+TAMs","cDC2",                           
                                                                                                                                           "Migratory cDCs","cDC1","Monocyte-derived dendritic cells",
                                                                                                                                           "ZBTB16+TAMs","Mast cells","CXCL10+TAMs",                    
                                                                                                                                           "Neutrophil"),do.stat = TRUE,color.use = c("#4974a4","#f29600"))
gg6 <- rankNet(cellchat, mode = "comparison",comparison = c(2, 3), stacked = T, font.size = 10,sources.use = "CD8_C1-GNLY",targets.use = c("SPP1+TAMs","FOLR2+TAMs","cDC2",                           
                                                                                                                                           "Migratory cDCs","cDC1","Monocyte-derived dendritic cells",
                                                                                                                                           "ZBTB16+TAMs","Mast cells","CXCL10+TAMs",                    
                                                                                                                                           "Neutrophil"),do.stat = TRUE,color.use = c("#4dae47","#f29600"))

gg4+gg5+gg6






gg7 <- rankNet(cellchat, mode = "comparison",comparison = c(1, 2), stacked = T,font.size = 10,do.stat = TRUE, color.use = c("#4974a4","#4dae47"))

gg8 <- rankNet(cellchat, mode = "comparison",comparison = c(1, 3), stacked = T, font.size = 10,do.stat = TRUE,color.use = c("#4974a4","#f29600"))
gg9 <- rankNet(cellchat, mode = "comparison",comparison = c(2, 3), stacked = T, font.size = 10,do.stat = TRUE,color.use = c("#4dae47","#f29600"))

gg7+gg8+gg9


library(ComplexHeatmap)
i = 1
# combining all the identified signaling pathways from different datasets 
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
object.list[[i]] <- netAnalysis_computeCentrality(object.list[[i]])
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 6, height = 12)
object.list[[i+1]] <- netAnalysis_computeCentrality(object.list[[i+1]])
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 6, height = 12)
object.list[[i+2]] <- netAnalysis_computeCentrality(object.list[[i+2]])
ht3 = netAnalysis_signalingRole_heatmap(object.list[[i+2]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+2], width = 6, height = 12)
draw(ht1 + ht2+ht3, ht_gap = unit(0.5, "cm"))


library(ComplexHeatmap)
i = 1
# combining all the identified signaling pathways from different datasets 
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
object.list[[i]] <- netAnalysis_computeCentrality(object.list[[i]])
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 6, height = 12)
object.list[[i+1]] <- netAnalysis_computeCentrality(object.list[[i+1]])
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+1], width = 6, height = 12)
object.list[[i+2]] <- netAnalysis_computeCentrality(object.list[[i+2]])
ht3 = netAnalysis_signalingRole_heatmap(object.list[[i+2]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+2], width = 6, height = 12)
draw(ht1 + ht2+ht3, ht_gap = unit(0.5, "cm"))

library(ComplexHeatmap)
i = 1
# combining all the identified signaling pathways from different datasets 
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
object.list[[i]] <- netAnalysis_computeCentrality(object.list[[i]])
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all", signaling = pathway.union, title = names(object.list)[i], width = 6, height = 12)
object.list[[i+1]] <- netAnalysis_computeCentrality(object.list[[i+1]])
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "all", signaling = pathway.union, title = names(object.list)[i+1], width = 6, height = 12)
object.list[[i+2]] <- netAnalysis_computeCentrality(object.list[[i+2]])
ht3 = netAnalysis_signalingRole_heatmap(object.list[[i+2]], pattern = "all", signaling = pathway.union, title = names(object.list)[i+2], width = 6, height = 12)
draw(ht1 + ht2+ht3, ht_gap = unit(0.5, "cm"))
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
pathway.union


pbubble<-netVisual_bubble(cellchat, sources.use = "CD8_C1-GNLY", targets.use =c("SPP1+TAMs","FOLR2+TAMs","cDC2",                           
                                                                       "Migratory cDCs","cDC1","Monocyte-derived dendritic cells",
                                                                       "ZBTB16+TAMs","Mast cells","CXCL10+TAMs",                    
                                                                       "Neutrophil"),comparison = c(1,2,3), angle.x = 45,font.size = 12,color.text =c("#4974a4","#4dae47","#f29600"))

pbubble<-netVisual_bubble(cellchat, sources.use = "ZBTB16+TAMs", targets.use =c("CD8_C1-GNLY"),comparison = c(1,2,3), angle.x = 45,font.size = 12,color.text =c("#4974a4","#4dae47","#f29600"))


pbubble<-netVisual_bubble(cellchat, sources.use = "CD8_C1-GNLY", targets.use =c("ZBTB16+TAMs"),comparison = c(1,2,3), angle.x = 45,font.size = 12,color.text =c("#4974a4","#4dae47","#f29600"))

?netVisual_bubble
 
ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/髓系细胞亚型/细胞通讯分析/GNLY_受配体_combined_pbubble.pdf", pbubble, width = 7, height = 7) 


pbubble<-netVisual_bubble(cellchat, sources.use = "cDC1", targets.use =c("CD8_C1-GNLY","CD8_C4-GZMK","CD8_C3-IL7R"),comparison = c(1,2,3), angle.x = 45,font.size = 12,color.text =c("#4974a4","#4dae47","#f29600"))


pbubble<-netVisual_bubble(cellchat, sources.use = c("CD8_C1-GNLY","CD8_C4-GZMK","CD8_C3-IL7R"), targets.use =c("cDC1"),comparison = c(1,2), angle.x = 45,font.size = 12,color.text =c("#4974a4","#4dae47","#f29600"))


?netVisual_bubble

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/髓系细胞亚型/细胞通讯分析/cDC1_受配体_combined_pbubble有细菌存在.pdf", pbubble, width = 9, height = 9) 


pbubble<-netVisual_bubble(cellchat, sources.use = "CD8_C4-GZMK", targets.use =c("SPP1+TAMs","FOLR2+TAMs","cDC2",                           
                                                                                "Migratory cDCs","cDC1","Monocyte-derived dendritic cells",
                                                                                "ZBTB16+TAMs","Mast cells","CXCL10+TAMs",                    
                                                                                "Neutrophil"),comparison = c(1,2,3), angle.x = 45,font.size = 12,color.text =c("#4974a4","#4dae47","#f29600"))
?netVisual_bubble

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/髓系细胞亚型/细胞通讯分析/GZMK_受配体_combined_pbubble.pdf", pbubble, width = 9, height = 9) 





 
pathways.show <- c("PD-L1")
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show)
par(mfrow = c(1,3), xpd=TRUE)
par(mfrow = c(1,2), xpd=TRUE)
pathways.show <- c("PD-L1")
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}


pathways.show <- c("MHC-I")
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max =10, signaling.name = paste(pathways.show, names(object.list)[i]))
}


pathways.show <- c("MHC-II")
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}


pathways.show <- c("MIF")
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}

pathways.show <- c("ITGB2")
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}


netVisual_chord_gene(object.list[1], lab.cex = 0.5,legend.pos.y =50)



gg1 <- netVisual_heatmap(cellchat)
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
gg1 + gg2

####nichenetr v2----

devtools::install_github("saeyslab/nichenetr")

library(nichenetr) 
library(Seurat)
library(SeuratObject)
library(tidyverse)


lr_network <- readRDS('D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/髓系细胞亚型/细胞通讯分析/nichenetr/nichent-human-network/lr_network_human_21122021.rds')
ligand_target_matrix <- readRDS('D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/髓系细胞亚型/细胞通讯分析/nichenetr/nichent-human-network/ligand_target_matrix_nsga2r_final.rds')# target genes in rows, ligands in columns
weighted_networks <- readRDS('D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/髓系细胞亚型/细胞通讯分析/nichenetr/nichent-human-network/weighted_networks_nsga2r_final.rds')

lr_network <- lr_network %>% distinct(from, to)

myeloid_cells_object=readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/myeloid_cells_object_filter使用.rds")
myeloid_cells_object = UpdateSeuratObject(myeloid_cells_object)
myeloid_cells_object

T_cells_object<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/T_cells_object.rds")
T_cells_object = UpdateSeuratObject(T_cells_object)
T_cells_object

merged_obj <- merge(myeloid_cells_object, y = T_cells_object, add.cell.ids = c("myeloid_cells_object", "T_cells_object"))
merged_obj$cd4_cd8_group
merged_obj$myeloid_cells_type
merged_obj@meta.data$merged_group <- paste(merged_obj@meta.data$cd4_cd8_group, 
                                           merged_obj@meta.data$myeloid_cells_type, sep = "_")

# 查看更新后的 meta.data
head(merged_obj@meta.data)
merged_obj@meta.data$merged_group <- gsub("^NA_", "", merged_obj@meta.data$merged_group)
merged_obj@meta.data$merged_group <- gsub("_NA$", "", merged_obj@meta.data$merged_group)

unique(merged_obj@meta.data$merged_group)

unique(merged_obj$MPR_Response2)

merged_obj
merged_obj_filter<- subset(merged_obj,subset=merged_group %in% c("SPP1+TAMs","FOLR2+TAMs","cDC2",                           
                                                                 "Migratory cDCs","cDC1","Monocyte-derived dendritic cells",
                                                                 "ZBTB16+TAMs","Mast cells","CXCL10+TAMs",                    
                                                                 "Neutrophil","CD8_C1-GNLY","CD8_C4-GZMK","CD8_C3-IL7R"))
Idents(merged_obj_filter)<- "merged_group"
table(Idents(merged_obj_filter))
receiver = "CD8_C1-GNLY"
expressed_genes_receiver <- get_expressed_genes(receiver, merged_obj, pct = 0.1)

#获取定义的receiver cell中表达基因中receprter基因
all_receptors <- unique(lr_network$to)  
expressed_receptors <- intersect(all_receptors, expressed_genes_receiver)
#将潜在配体定义为所有同源受体表达的配体。
potential_ligands <- lr_network %>% filter(to %in% expressed_receptors) %>% pull(from) %>% unique()

#define sender cell types，作用于receiver，也可以包括receiver cell自身
#还是一样，需要返回sender cell返回在给定细胞簇中表达的基因。
sender_celltypes <- c("SPP1+TAMs","FOLR2+TAMs","cDC2",                           
                      "Migratory cDCs","cDC1","Monocyte-derived dendritic cells",
                      "ZBTB16+TAMs","Mast cells","CXCL10+TAMs",                    
                      "Neutrophil")
list_expressed_genes_sender <- sender_celltypes %>% unique() %>% lapply(get_expressed_genes, merged_obj_filter, 0.1)
expressed_genes_sender <- list_expressed_genes_sender %>% unlist() %>% unique()
#sender cell表达的配体
potential_ligands_focused <- intersect(potential_ligands, expressed_genes_sender) 

seurat_obj_receiver <- subset(merged_obj_filter, idents = receiver)
DE_table_receiver <-  FindMarkers(object = seurat_obj_receiver,
                                  ident.1 = "Pre_MPR", 
                                  ident.2 = "Pre_NMPR",
                                  group.by = "MPR_Response2",
                                  min.pct = 0.1) %>% rownames_to_column("gene")

geneset_oi <- DE_table_receiver %>% filter(p_val_adj <= 0.05 & abs(avg_log2FC) >= 0.25) %>% pull(gene)
geneset_oi <- geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]

background_expressed_genes <- expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]


ligand_activities <- predict_ligand_activities(geneset = geneset_oi,
                                               background_expressed_genes = background_expressed_genes,
                                               ligand_target_matrix = ligand_target_matrix,
                                               potential_ligands = potential_ligands)
ligand_activities <- ligand_activities %>% arrange(-aupr_corrected) %>% mutate(rank = rank(desc(aupr_corrected)))

best_upstream_ligands <- ligand_activities %>% top_n(30, aupr_corrected) %>% arrange(-aupr_corrected) %>% pull(test_ligand)


active_ligand_target_links_df <- best_upstream_ligands %>%
  lapply(get_weighted_ligand_target_links,
         geneset = geneset_oi,
         ligand_target_matrix = ligand_target_matrix,
         n = 200) %>%
  bind_rows() %>% drop_na()

ligand_receptor_links_df <- get_weighted_ligand_receptor_links(
  best_upstream_ligands, expressed_receptors,
  lr_network, weighted_networks$lr_sig) 

ligand_activities_all <- ligand_activities 
best_upstream_ligands_all <- best_upstream_ligands

ligand_activities <- ligand_activities %>% filter(test_ligand %in% potential_ligands_focused)
best_upstream_ligands <- ligand_activities %>% top_n(30, aupr_corrected) %>% arrange(-aupr_corrected) %>%
  pull(test_ligand) %>% unique()

active_ligand_target_links_df <- best_upstream_ligands %>%
  lapply(get_weighted_ligand_target_links,
         geneset = geneset_oi,
         ligand_target_matrix = ligand_target_matrix,
         n = 200) %>%
  bind_rows() %>% drop_na()

ligand_receptor_links_df <- get_weighted_ligand_receptor_links(
  best_upstream_ligands, expressed_receptors,
  lr_network, weighted_networks$lr_sig) 

DotPlot(subset(merged_obj_filter, merged_group %in% sender_celltypes),
        features = rev(best_upstream_ligands), cols = "RdYlBu") + 
  coord_flip() +
  scale_y_discrete(position = "right")

celltype_order <- levels(Idents(merged_obj_filter)) 

DE_table_top_ligands <- lapply(
  celltype_order[celltype_order %in% sender_celltypes],
  get_lfc_celltype, 
  seurat_obj = merged_obj_filter,
  condition_colname = "MPR_Response2",
  condition_oi = "Pre_MPR",
  condition_reference = "Pre_NMPR",
  celltype_col = "merged_group",
  min.pct = 0, logfc.threshold = 0,
  features = best_upstream_ligands 
) 

DE_table_top_ligands <- DE_table_top_ligands %>%  reduce(., full_join) %>% column_to_rownames("gene") 
vis_ligand_lfc <- as.matrix(DE_table_top_ligands[rev(best_upstream_ligands), , drop = FALSE])

p_lfc <- make_threecolor_heatmap_ggplot(vis_ligand_lfc,
                                        "Prioritized ligands", "LFC in Sender",
                                        low_color = "midnightblue", mid_color = "white",
                                        mid = median(vis_ligand_lfc), high_color = "red",
                                        legend_title = "LFC")

p_lfc


# compare rankings between the sender-agnostic and sender-focused approach
(make_line_plot(ligand_activities = ligand_activities_all,
                potential_ligands = potential_ligands_focused) +
    theme(plot.title = element_text(size=11, hjust=0.1, margin=margin(0, 0, -5, 0))))

####nichenet高级点图绘制----
####"CD8_C1-GNLY","CD8_C4-GZMK","CD8_C3-IL7R"

receiver="CD8_C3-IL7R"  
sender=c("SPP1+TAMs","FOLR2+TAMs","cDC2",                           
         "Migratory cDCs","cDC1","Monocyte-derived dendritic cells",
         "ZBTB16+TAMs","Mast cells","CXCL10+TAMs",                    
         "Neutrophil")

sender=c("ZBTB16+TAMs")

receiver="ZBTB16+TAMs"

sender=c("CD8_C4-GZMK")


receiver="CD8_C3-IL7R"

sender=c("ZBTB16+TAMs")

receiver="cDC1"

sender=c("CD8_C3-IL7R","CD8_C4-GZMK","CD8_C1-GNLY")

myeloid_t_bacteria_no_chat <- subset(merged_obj_filter, subset=group_microbe2=="Bacteria-")
myeloid_t_bacteria_yes_chat <- subset(merged_obj_filter, subset=group_microbe2=="Bacteria+")

merged_obj_filter$cell_cluster = ifelse(myeloid_t_bacteria_no_chat$merged_group==receiver, receiver,"other_cells")
DEG_Fib<-FindMarkers(merged_obj_filter, ident.1=receiver, ident.2="other_cells", only.pos=TRUE, logfc.threshold=0.2, group.by = "cell_cluster")
expressed_genes_receiver = get_expressed_genes(receiver, merged_obj_filter, pct = 0.1) %>% .[. %in% rownames(DEG_Fib)] 


receiver="ZBTB16+TAMs"

sender=c("CD8_C3-IL7R")

expressed_genes_receiver = get_expressed_genes(receiver, merged_obj_filter, pct = 0.1)


#sender expressed_genes
list_expressed_genes_sender <- sender %>% unique() %>% lapply(get_expressed_genes,merged_obj_filter, 0.1)
expressed_genes_sender <- list_expressed_genes_sender %>% unlist() %>% unique()


#potential_ligands
expressed_receptors <- intersect(unique(lr_network$to) , expressed_genes_receiver)
expressed_ligands = intersect(unique(lr_network$from),expressed_genes_sender)
potential_ligands = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% pull(from) %>% unique()


#Define the gene set of interest
seurat_obj_receiver <- subset(merged_obj_filter, idents = receiver)
seurat_obj_receiver$group_cluster = ifelse(seurat_obj_receiver$MPR_Response2=="Post_MPR", "Post_MPR","other_group")
DE_table_receiver <-  FindMarkers(object = seurat_obj_receiver,ident.1 = "Post_MPR", ident.2 = "other_group",group.by = "group_cluster",min.pct = 0.1) %>% rownames_to_column("gene")
geneset_oi <- DE_table_receiver %>% filter(p_val_adj <= 0.05 & abs(avg_log2FC) >= 0.25) %>% pull(gene)
geneset_oi <- geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]

#Define the background genes
background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)] 


#nichenet analysis
ligand_activities = predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)
ligand_activities = ligand_activities %>% arrange(-pearson) %>% mutate(rank = rank(desc(pearson)))
best_upstream_ligands = ligand_activities %>% top_n(20, pearson) %>% arrange(-pearson) %>% pull(test_ligand) %>% unique()

# targets
active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 200) %>% bind_rows() %>% drop_na()

# receptors
ligand_receptor_links_df <- get_weighted_ligand_receptor_links(best_upstream_ligands, expressed_receptors,lr_network, weighted_networks$lr_sig) 
#分析结束

#-------------------------------------------------------------------------------------------
#可视化1：
#p1-ligands排序
vis_ligand_aupr <- ligand_activities %>% filter(test_ligand %in% best_upstream_ligands) %>%
  column_to_rownames("test_ligand") %>% select(aupr_corrected) %>% arrange(aupr_corrected) %>% as.matrix(ncol = 1)

p1 = make_heatmap_ggplot(vis_ligand_aupr,"Prioritized ligands","Ligand activity",legend_title = "AUPR", color ="#EFBD37")+
  theme_bw()+
  theme(axis.text.x.top = element_blank(),
        axis.text = element_text(colour = 'black'),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_text(size = 10,face ="bold",hjust = 1))+
  xlab("Ligand\nactivity")

#p2-expression sender
p2 <- DotPlot(subset(merged_obj_filter, merged_group %in% sender),
              features = rev(best_upstream_ligands), cols = "RdYlBu") + 
  coord_flip() +
  scale_y_discrete(position = "right")+
  theme_bw()+
  geom_point(shape=21, aes(size=pct.exp),stroke=1)+
  theme(axis.text.x = element_text(colour = 'black'), size=10)+
  theme(axis.ticks = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 10,face ="bold"),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 9,  angle = 90, hjust = 0)) +
  ylab("Ligand Expression\nin Sender")

#p3-ligands-target heatmap
# Target gene plot
active_ligand_target_links <- prepare_ligand_target_visualization(
  ligand_target_df = active_ligand_target_links_df,
  ligand_target_matrix = ligand_target_matrix,
  cutoff = 0) 

order_ligands <- intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
order_targets <- active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links))

vis_ligand_target <- t(active_ligand_target_links[order_targets,order_ligands])

p3 <- make_heatmap_ggplot(vis_ligand_target, "Prioritized ligands", "Predicted target genes",
                          color = '#810f7c', legend_title = "Regulatory potential") +
  theme(axis.text.x.top = element_text(colour = 'black', face = "italic", angle = 90, hjust = 0),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(size = 10,face ="bold"),
        panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))

#组合展示
library(patchwork)
#组合图，将注释的图形宽度设置小一点，添加参数guides = 'collect'，这样两幅图的legend就结合在一起了。
pp=p1+p2+p3+plot_layout(ncol = 3, widths  = c(2, length(sender)+7, ncol(vis_ligand_target)),guides = 'collect')
pp & theme(plot.margin = margin(1,2,2,0),
           legend.position = "bottom",
           legend.direction = 'horizontal',
           # legend.justification=c(0,1),
           legend.key.width=unit(0.5,"cm"),
           legend.key.height = unit(0.3, "cm"),
           legend.title.position = 'top')

#receptor
vis_ligand_receptor_network <- prepare_ligand_receptor_visualization(
  ligand_receptor_links_df,
  best_upstream_ligands,
  order_hclust = "both") 

p4  = make_heatmap_ggplot(t(vis_ligand_receptor_network), 
                          y_name = "Ligands", x_name = "Receptors",  
                          color ="springgreen4", legend_title = "Prior interaction potential")+
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(size = 10,face ="bold"),
        panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))


#组合
pp2=p1+p4+p3+plot_layout(ncol = 3, widths  = c(2, nrow(vis_ligand_receptor_network), ncol(vis_ligand_target)),guides = 'collect')
pp2 & theme(plot.margin = margin(1,2,2,0),
            legend.position = "bottom",
            legend.direction = 'horizontal',
            # legend.justification=c(0,1),
            legend.key.width=unit(0.5,"cm"),
            legend.key.height = unit(0.3, "cm"),
            legend.title.position = 'top')

#-------------------------------------------------------------------------------------------
#可视化2：弦图展示受配体
#Circos plots to visualize ligand-target and ligand-receptor interactions
#Assign ligands to sender cells
Assign_ligands_to_celltype <- function(object,ligands,celltype_col,sender){
  
  Average_sender_ligand_exp <- AverageExpression(object,features = ligands,group.by = celltype_col,slot = 'data') 
  Average_sender_ligand_exp <- as.data.frame(Average_sender_ligand_exp$RNA)
  Average_sender_ligand_exp <- Average_sender_ligand_exp[,sender]
  
  func.assign <- function(x) {mean(x)+sd(x)}
  sender_ligand_assignment <- Average_sender_ligand_exp %>% apply(1, function(ligand_expression){
    ligand_expression > func.assign(ligand_expression)}) %>% t()
  
  sender_ligand_assignment <- sender_ligand_assignment %>% apply(2, function(x){x[x == TRUE]}) %>% purrr::keep(function(x){length(x) > 0})
  
  all_assigned_ligands = sender_ligand_assignment %>% lapply(function(x){names(x)}) %>% unlist()
  unique_ligands = all_assigned_ligands %>% table() %>% .[. == 1] %>% names()
  general_ligands = ligands %>% setdiff(unique_ligands)
  
  ligand_type_indication_df <- lapply(names(sender_ligand_assignment), function(sender) {
    unique_ligands_sender <- names(sender_ligand_assignment[[sender]]) %>% setdiff(general_ligands)
    
    if (length(unique_ligands_sender) > 0) {
      return(data.frame(ligand_type = sender, ligand = unique_ligands_sender))
    }}) %>% bind_rows()
  
  ligand_type_indication_df <- bind_rows(ligand_type_indication_df,data.frame(ligand = general_ligands) %>% mutate(ligand_type = "General"))
  
  
  
}


ligand_type_indication_df <- Assign_ligands_to_celltype(merged_obj_filter,
                                                        best_upstream_ligands,
                                                        celltype_col = " merged_group",
                                                        sender = sender) 


#Define the ligand-target links of interest
active_ligand_target_links_df$target_type <- "Post_MPR" # needed for joining tables
circos_links <- get_ligand_target_links_oi(ligand_type_indication_df,
                                           active_ligand_target_links_df,
                                           cutoff = 0.25) 



ligand_colors <- c("General" = "#377EB8", "Tcell" = "#4DAF4A", "Endo" = "#984EA3",
                   "LY" = "#FF7F00", "Fib" = "#F781BF") 
target_colors <- c("Post_MPR" = "#999999") 

vis_circos_obj <- prepare_circos_visualization(circos_links,
                                               ligand_colors = ligand_colors,
                                               target_colors = target_colors,
                                               celltype_order = NULL) 


make_circos_plot(vis_circos_obj, transparency = FALSE,  args.circos.text = list(cex = 0.5)) 



#ligand-receptor
lr_network_top_df <- get_weighted_ligand_receptor_links(best_upstream_ligands,
                                                        expressed_receptors,
                                                        lr_network,
                                                        weighted_networks$lr_sig) %>%
  rename(ligand=from, target=to) %>% 
  mutate(target_type = "SDDDD") %>% 
  inner_join(ligand_type_indication_df)

receptor_colors <- c("SDDDD" = "#E41A1C")

vis_circos_receptor_obj <- prepare_circos_visualization(lr_network_top_df,
                                                        ligand_colors = ligand_colors,
                                                        target_colors = receptor_colors) 


make_circos_plot(vis_circos_receptor_obj, transparency = T,link.visible = TRUE,  args.circos.text = list(cex = 0.8)) 



####中央记忆T细胞通讯分析----

library(CellChat)
library(tidyverse)
library(Seurat)

myeloid_cells_object=readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/myeloid_cells_object_filter使用.rds")
myeloid_cells_object = UpdateSeuratObject(myeloid_cells_object)
myeloid_cells_object

T_cells_object<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/T_cells_object.rds")
T_cells_object = UpdateSeuratObject(T_cells_object)
T_cells_object

merged_obj <- merge(myeloid_cells_object, y = T_cells_object, add.cell.ids = c("myeloid_cells_object", "T_cells_object"))

View(merged_obj)
merged_obj$cd4_cd8_group
merged_obj$myeloid_cells_type
merged_obj@meta.data$merged_group <- paste(merged_obj@meta.data$cd4_cd8_group, 
                                           merged_obj@meta.data$myeloid_cells_type, sep = "_")

# 查看更新后的 meta.data
head(merged_obj@meta.data)
merged_obj@meta.data$merged_group <- gsub("^NA_", "", merged_obj@meta.data$merged_group)
merged_obj@meta.data$merged_group <- gsub("_NA$", "", merged_obj@meta.data$merged_group)

unique(merged_obj@meta.data$merged_group)

unique(merged_obj$MPR_Response2)

merged_obj


####分组比较----

KS_cellchat <- function(input_obj,
                        assay= NULL,
                        group.by = NULL,
                        workers,
                        species=c('human','mouse'),
                        CellChatDB.use=NULL,#a character vector, which is a subset of c("Secreted Signaling","ECM-Receptor","Cell-Cell Contact","Non-protein Signaling")
                        PPIuse=F,
                        type="triMean",# c("triMean", "truncatedMean", "thresholdedMean", "median")
                        min.cells = 10
){
  
  
  
  cellchat.obj = createCellChat(input_obj, assay = assay, group.by = group.by)
  
  if(species=='human'){
    
    CellChatDB <- CellChatDB.human
    ppi = PPI.human
  }
  
  if(species =="mouse"){
    
    CellChatDB <- CellChatDB.mouse
    ppi = PPI.mouse
  }
  
  
  if(is.null(CellChatDB.use)){
    
    cellchat.obj@DB <- CellChatDB
    
  }else{
    
    CellChatDB <- subsetDB(CellChatDB, search = CellChatDB.use, key = "annotation")
    cellchat.obj@DB <- CellChatDB
  }
  
  cellchat.obj <- subsetData(cellchat.obj) 
  future::plan("multisession", workers = workers) # 并行计算
  cellchat.obj <- identifyOverExpressedGenes(cellchat.obj)
  cellchat.obj <- identifyOverExpressedInteractions(cellchat.obj)
  
  if(PPIuse==F){
    
    cellchat.obj <- computeCommunProb(cellchat.obj, type = type)
    
  }else{
    
    cellchat.obj <- projectData(cellchat.obj, ppi)
    cellchat.obj <- computeCommunProb(cellchat.obj, raw.use=F, type = type)
  }
  
  
  cellchat.obj <- filterCommunication(cellchat.obj, min.cells = min.cells)
  cellchat.obj <- computeCommunProbPathway(cellchat.obj)
  
  #互作网络整合,可以设置soure和target，我们这里默认了，就是全部的
  cellchat.obj <- aggregateNet(cellchat.obj)
  
  return(cellchat.obj)
  
}

merged_obj_filter<- subset(merged_obj,subset=merged_group %in% c("SPP1+TAMs","FOLR2+TAMs","cDC2",                           
                                                                 "Migratory cDCs","cDC1","Monocyte-derived dendritic cells",
                                                                 "ZBTB16+TAMs","Mast cells","CXCL10+TAMs",                    
                                                                 "Neutrophil","CD8_C3-IL7R","CD8_C4-GZMK"))


pre_nmpr_myeloid_t_chat<- subset(merged_obj_filter, subset=MPR_Response2=="Pre_NMPR")
unique(pre_nmpr_myeloid_t_chat$MPR_Response2)
pre_mpr_myeloid_t_chat<- subset(merged_obj_filter, subset=MPR_Response2=="Pre_MPR")
unique(pre_mpr_myeloid_t_chat$MPR_Response2)
post_mpr_myeloid_t_chat<- subset(merged_obj_filter, subset=MPR_Response2=="Post_MPR")
unique(post_mpr_myeloid_t_chat$MPR_Response2)
unique(post_mpr_myeloid_t_chat$CRNCR_Response2)


pre_nmpr.cellchat <- KS_cellchat(pre_nmpr_myeloid_t_chat, assay = 'RNA', group.by = "merged_group",
                                 workers=2, species='human')

saveRDS(pre_nmpr.cellchat, file = "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/组织驻留T细胞/细胞通讯分析/和髓系细胞通讯分析/pre_nmpr.cellchat.rds")

options(future.globals.maxSize = 1e+09) 
pre_mpr.cellchat <- KS_cellchat(pre_mpr_myeloid_t_chat, assay = 'RNA', group.by = "merged_group",
                                workers=2, species='human')

saveRDS(pre_mpr.cellchat , file = "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/组织驻留T细胞/细胞通讯分析/和髓系细胞通讯分析/pre_mpr.cellchat.rds")


post_mpr.cellchat <- KS_cellchat(post_mpr_myeloid_t_chat, assay = 'RNA', group.by = "merged_group",
                                 workers=2, species='human')

saveRDS(post_mpr.cellchat, file = "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/组织驻留T细胞/细胞通讯分析/和髓系细胞通讯分析/post_mpr.cellchat.rds")


##############cr 和ncr


####cellchat多组比较分析----

pre_nmpr.cellchat<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/组织驻留T细胞/细胞通讯分析/和髓系细胞通讯分析/pre_nmpr.cellchat.rds")

pre_mpr.cellchat<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/组织驻留T细胞/细胞通讯分析/和髓系细胞通讯分析/pre_mpr.cellchat.rds")

post_mpr.cellchat<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/组织驻留T细胞/细胞通讯分析/和髓系细胞通讯分析/post_mpr.cellchat.rds")

object.list <- list(Pre_NMPR= pre_nmpr.cellchat, Pre_MPR= pre_mpr.cellchat,Post_MPR=post_mpr.cellchat)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))


gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2,3))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2,3), measure = "weight")
p1<-gg1 + gg2

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/组织驻留T细胞/细胞通讯分析/和髓系细胞通讯分析/分组比较柱状图.pdf",p1,width=6,height=4)


gg1 <- netVisual_heatmap(cellchat)
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
pdf(file = "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/组织驻留T细胞/细胞通讯分析/和髓系细胞通讯分析/分组比较热图.pdf", width = 10, height = 6)
gg1 + gg2
dev.off()

num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  # 计算网络中心性得分
  object.list[[i]] <- netAnalysis_computeCentrality(object.list[[i]], slot.name = "netP") 
  # 然后生成散点图
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}
patchwork::wrap_plots(plots = gg)

p_combined <- patchwork::wrap_plots(plots = gg)
ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/组织驻留T细胞/细胞通讯分析/和髓系细胞通讯分析/combined_plots.pdf", p_combined, width = 10, height = 4) 

?rankNet
treatment_color <- c("#4974a4","#4dae47","#f29600")

gg1 <- rankNet(cellchat, mode = "comparison",comparison = c(1, 2), stacked = T,font.size = 10, sources.use = "CD8_C3-IL7R",targets.use = c("SPP1+TAMs","FOLR2+TAMs","cDC2",                           
                                                                                                                                           "Migratory cDCs","cDC1","Monocyte-derived dendritic cells",
                                                                                                                                           "ZBTB16+TAMs","Mast cells","CXCL10+TAMs",                    
                                                                                                                                           "Neutrophil"),do.stat = TRUE, color.use = c("#4974a4","#4dae47"))

gg2 <- rankNet(cellchat, mode = "comparison",comparison = c(1, 3), stacked = T, font.size = 10,sources.use = "CD8_C3-IL7R",targets.use = c("SPP1+TAMs","FOLR2+TAMs","cDC2",                           
                                                                                                                                           "Migratory cDCs","cDC1","Monocyte-derived dendritic cells",
                                                                                                                                           "ZBTB16+TAMs","Mast cells","CXCL10+TAMs",                    
                                                                                                                                           "Neutrophil"),do.stat = TRUE,color.use = c("#4974a4","#f29600"))
gg3 <- rankNet(cellchat, mode = "comparison",comparison = c(2, 3), stacked = T, font.size = 10,sources.use = "CD8_C3-IL7R",targets.use = c("SPP1+TAMs","FOLR2+TAMs","cDC2",                           
                                                                                                                                           "Migratory cDCs","cDC1","Monocyte-derived dendritic cells",
                                                                                                                                           "ZBTB16+TAMs","Mast cells","CXCL10+TAMs",                    
                                                                                                                                           "Neutrophil"),do.stat = TRUE,color.use = c("#4dae47","#f29600"))

gg4 <- rankNet(cellchat, mode = "comparison",comparison = c(1, 2), stacked = T,font.size = 10, sources.use = "CD8_C4-GZMK",targets.use = c("SPP1+TAMs","FOLR2+TAMs","cDC2",                           
                                                                                                                                           "Migratory cDCs","cDC1","Monocyte-derived dendritic cells",
                                                                                                                                           "ZBTB16+TAMs","Mast cells","CXCL10+TAMs",                    
                                                                                                                                           "Neutrophil"),do.stat = TRUE, color.use = c("#4974a4","#4dae47"))

gg5 <- rankNet(cellchat, mode = "comparison",comparison = c(1, 3), stacked = T, font.size = 10,sources.use = "CD8_C4-GZMK",targets.use = c("SPP1+TAMs","FOLR2+TAMs","cDC2",                           
                                                                                                                                           "Migratory cDCs","cDC1","Monocyte-derived dendritic cells",
                                                                                                                                           "ZBTB16+TAMs","Mast cells","CXCL10+TAMs",                    
                                                                                                                                           "Neutrophil"),do.stat = TRUE,color.use = c("#4974a4","#f29600"))
gg6 <- rankNet(cellchat, mode = "comparison",comparison = c(2, 3), stacked = T, font.size = 10,sources.use = "CD8_C4-GZMK",targets.use = c("SPP1+TAMs","FOLR2+TAMs","cDC2",                           
                                                                                                                                           "Migratory cDCs","cDC1","Monocyte-derived dendritic cells",
                                                                                                                                           "ZBTB16+TAMs","Mast cells","CXCL10+TAMs",                    
                                                                                                                                           "Neutrophil"),do.stat = TRUE,color.use = c("#4dae47","#f29600"))

gg1+gg2+gg3

gg4+gg5+gg6


gg7 <- rankNet(cellchat, mode = "comparison",comparison = c(1, 2), stacked = T,font.size = 10,do.stat = TRUE, color.use = c("#4974a4","#4dae47"))

gg8 <- rankNet(cellchat, mode = "comparison",comparison = c(1, 3), stacked = T, font.size = 10,do.stat = TRUE,color.use = c("#4974a4","#f29600"))
gg9 <- rankNet(cellchat, mode = "comparison",comparison = c(2, 3), stacked = T, font.size = 10,do.stat = TRUE,color.use = c("#4dae47","#f29600"))

gg7+gg8+gg9


library(ComplexHeatmap)
i = 1
# combining all the identified signaling pathways from different datasets 
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
object.list[[i]] <- netAnalysis_computeCentrality(object.list[[i]])
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 6, height = 12)
object.list[[i+1]] <- netAnalysis_computeCentrality(object.list[[i+1]])
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 6, height = 12)
object.list[[i+2]] <- netAnalysis_computeCentrality(object.list[[i+2]])
ht3 = netAnalysis_signalingRole_heatmap(object.list[[i+2]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+2], width = 6, height = 12)
draw(ht1 + ht2+ht3, ht_gap = unit(0.5, "cm"))


library(ComplexHeatmap)
i = 1
# combining all the identified signaling pathways from different datasets 
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
object.list[[i]] <- netAnalysis_computeCentrality(object.list[[i]])
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 6, height = 12)
object.list[[i+1]] <- netAnalysis_computeCentrality(object.list[[i+1]])
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+1], width = 6, height = 12)
object.list[[i+2]] <- netAnalysis_computeCentrality(object.list[[i+2]])
ht3 = netAnalysis_signalingRole_heatmap(object.list[[i+2]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+2], width = 6, height = 12)
draw(ht1 + ht2+ht3, ht_gap = unit(0.5, "cm"))


pbubble<-netVisual_bubble(cellchat, sources.use = "CD8_C3-IL7R", targets.use =c("SPP1+TAMs","FOLR2+TAMs","cDC2",                           
                                                                                "Migratory cDCs","cDC1","Monocyte-derived dendritic cells",
                                                                                "ZBTB16+TAMs","Mast cells","CXCL10+TAMs",                    
                                                                                "Neutrophil"),comparison = c(1,2,3), angle.x = 45,font.size = 12,color.text =c("#4974a4","#4dae47","#f29600"))
?netVisual_bubble

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/组织驻留T细胞/细胞通讯分析/和髓系细胞通讯分析/IL7R_受配体_combined_pbubble.pdf", pbubble, width = 9, height = 9) 


pbubble<-netVisual_bubble(cellchat, sources.use = "CD8_C3-GZMK", targets.use =c("SPP1+TAMs","FOLR2+TAMs","cDC2",                           
                                                                                "Migratory cDCs","cDC1","Monocyte-derived dendritic cells",
                                                                                "ZBTB16+TAMs","Mast cells","CXCL10+TAMs",                    
                                                                                "Neutrophil"),comparison = c(1,2,3), angle.x = 45,font.size = 12,color.text =c("#4974a4","#4dae47","#f29600"))
?netVisual_bubble

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/组织驻留T细胞/细胞通讯分析/和髓系细胞通讯分析/GZMK_受配体_combined_pbubble.pdf", pbubble, width = 9, height = 9) 



pathways.show <- c("PD-L1")
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show)
par(mfrow = c(1,3), xpd=TRUE)


pathways.show <- c("MHC-I")
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}

pathways.show <- c("MHC-II")
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}

pathways.show <- c("MIF")
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}

pathways.show <- c("ITGB2")

for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}


netVisual_chord_gene(object.list[1], lab.cex = 0.5,legend.pos.y =50)



####髓系细胞相关性热图----

myeloid_cells_object=readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/myeloid_cells_object_filter使用.rds")
myeloid_cells_object = UpdateSeuratObject(myeloid_cells_object)
myeloid_cells_object

T_cells_object<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/T_cells_object.rds")
T_cells_object = UpdateSeuratObject(T_cells_object)
T_cells_object

merged_obj <- merge(myeloid_cells_object, y = T_cells_object, add.cell.ids = c("myeloid_cells_object", "T_cells_object"))
merged_obj$cd4_cd8_group
merged_obj$myeloid_cells_type
merged_obj@meta.data$merged_group <- paste(merged_obj@meta.data$cd4_cd8_group, 
                                           merged_obj@meta.data$myeloid_cells_type, sep = "_")

# 查看更新后的 meta.data
head(merged_obj@meta.data)
merged_obj@meta.data$merged_group <- gsub("^NA_", "", merged_obj@meta.data$merged_group)
merged_obj@meta.data$merged_group <- gsub("_NA$", "", merged_obj@meta.data$merged_group)

unique(merged_obj@meta.data$merged_group)

unique(merged_obj$MPR_Response2)

merged_obj
merged_obj_filter<- subset(merged_obj,subset=merged_group %in% c("SPP1+TAMs","FOLR2+TAMs","cDC2",                           
                                                                 "Migratory cDCs","cDC1","Monocyte-derived dendritic cells",
                                                                 "ZBTB16+TAMs","Mast cells","CXCL10+TAMs",                    
                                                                 "Neutrophil","CD8_C1-GNLY","CD8_C4-GZMK","CD8_C3-IL7R"))

avg_expression <- AverageExpression(merged_obj_filter, group.by = "merged_group", return.seurat = FALSE)

unique(merged_obj_filter@meta.data$merged_group)

head(avg_expression$RNA)

# 提取基因表达量矩阵
expr_matrix <- as.matrix(avg_expression$RNA)

# 检查矩阵的维度
print(dim(expr_matrix))
rownames(expr_matrix)
colnames(expr_matrix)
# 计算相关性矩阵
cor_matrix <- cor(expr_matrix, method = "pearson")
dim(cor_matrix)
# cor_matrix <- cor(t(expr_matrix), method = "pearson")
# 查看相关性矩阵
print(cor_matrix)

library(pheatmap)

pheatmap(cor_matrix, 
         clustering_distance_rows = "euclidean", 
         clustering_distance_cols = "euclidean", 
         color = colorRampPalette(c("blue", "white", "red"))(50))

library(tidyverse)
library(RColorBrewer)
library(ggtext)
library(magrittr)
library(reshape)
library(psych)
# install.packages("devtools")
# devtools::install_github("Hy4m/linkET", force = TRUE)
library(linkET)
library(ggnewscale)
library(corrplot)


# 计算表达量相关性
cor_matrix <- cor(expr_matrix, method = "pearson")
dim(cor_matrix)
rownames(cor_matrix)
res1 <- cor.mtest(expr_matrix, conf.level = 0.95)

####提取mantel 计算cd8 GNLY 和cd8 CTSW
View(cor_matrix)

cor_df <- as.data.frame(cor_matrix)

cor_df$Other_Group <- rownames(cor_df)

# 使用 pivot_longer 将两列转换为长格式
names()

long_df <- cor_df %>%
  pivot_longer(cols = c(`CD8-C1-GNLY`,"CD8-C3-IL7R","CD8-C4-GZMK"), 
               names_to = "Group", 
               values_to = "Value") %>% as.data.frame()


long_df <- long_df[, c("Group", "Other_Group", "Value")]


names(long_df)<-c("spec","env","r")


p_df<- as.data.frame(res1$p)
p_df$Other_Group <- rownames(p_df)




p_long_df <- p_df %>%
  pivot_longer(cols = c(`CD8-C1-GNLY`, "CD8-C3-IL7R","CD8-C4-GZMK"), 
               names_to = "Group", 
               values_to = "Value") %>% as.data.frame()


p_long_df <- p_long_df[, c("Group", "Other_Group", "Value")]

names(p_long_df)<-c("spec","env","p")
View(p_long_df)
View(long_df)


library(dplyr)

# 假设 df1 是包含 'p' 列的数据框，df2 是包含 'r' 列的数据框
# 通过 spec 和 env 两列进行合并
merged_df <- left_join(long_df, p_long_df, by = c("spec", "env"))

# 查看合并后的结果
head(merged_df)
mantel<-merged_df %>%
  mutate(rd = cut(r, breaks = c(-Inf, 0.5,0.8, 0.9, Inf),
                  labels = c("< 0.5", "0.5-0.8","0.8 - 0.9", ">= 0.9")),
         pd = cut(p, breaks = c(-Inf, 0.0001, 0.001, Inf),
                  labels = c("< 0.0001", "0.0001 - 0.001", ">= 0.001")))


####其他分组矩阵
remove_cells <- c("CD8-C1-GNLY", "CD8-C2-CTSW")

remove_cells <- c("CD8-C1-GNLY", "CD8-C3-IL7R","CD8-C4-GZMK")

remove_cells <- c("CD8-C3-IL7R","CD8-C4-GZMK")

# 去除行和列
cor_matrix_filtered1 <- cor_matrix[!rownames(cor_matrix) %in% remove_cells, 
                                   !colnames(cor_matrix) %in% remove_cells]

# 查看过滤后的矩阵
print(cor_matrix_filtered1)
cor_matrix_filtered<-cor_matrix_filtered1
# 使用 Pearson 方法计算表格 table1 的相关系数矩阵
cor_matrix_filtered[lower.tri(cor_matrix_filtered)] <- NA
# 将相关系数矩阵的下三角部分（不含对角线）设为 NA
cor_matrix_filtered[cor_matrix_filtered == 1] <- NA
# 将相关系数为 1（即完全相关）的元素设为 NA
cor_data <- melt(cor_matrix_filtered)  # 将相关系数矩阵转换为长格式数据框，方便后续处理
# 新增一列，标识是否是对角线元素（即 X1 是否等于 X2）
cor_data$is_diag <- cor_data$X1 == cor_data$X2

####计算P值
# 计算p值
res1 <- cor.mtest(expr_matrix, conf.level = 0.95)
View(res1$p)
res1_p<- res1$p
res1_p_filtered <- res1_p[!rownames(res1_p) %in% remove_cells, 
                          !colnames(res1_p) %in% remove_cells]
# 查看过滤后的矩阵
print(res1_p_filtered)

res1_p_filtered[lower.tri(res1_p_filtered)] <- NA

# 转化p值格式
p_value <- res1_p_filtered%>% as.data.frame() %>%
  rownames_to_column(var="id") %>% 
  pivot_longer(-id) %>% 
  drop_na() %>% 
  filter(id !=name) %>% 
  mutate(p_signif=symnum(value,cutpoints = c(0, 0.001, 0.01, 0.05,1),
                         symbols = c("***", "**", "*", "")))


bt_plot<- qcorrplot(correlate(cor_matrix_filtered1,method = "pearson"),diag=T,type="lower",grid_col = "white")+
  # 根据cor_data数据设置边框颜色，若是对角线数据则设置为白色，如此则多出一个空位用于文本放置
  geom_tile(data=cor_data %>% drop_na(),
            aes(x = X1,y = X2,color = is_diag),
            fill="white",linewidth = 0.4,width=1,show.legend = F,inherit.aes = F)+
  scale_color_manual(values = c("black","white"))+
  new_scale_color()+ # 定义新的颜色映射
  # 添加相关性
  geom_point(data=cor_data %>% drop_na(),color="white",pch=22,
             inherit.aes = F,
             aes(x = X1,y = X2,fill=value,size=abs(value)))+
  # 添加p值文本
  geom_text(data=p_value,aes(id,name,label=p_signif),inherit.aes = F,vjust=1)+
  guides(size="none")+
  scale_size_continuous(range = c(0,10)) +
  new_scale("size")+
  scale_fill_gradientn(colors=rev(COL2('RdBu',10)),na.value = "white")+
  # 添加链接线
  geom_couple(aes(colour=pd,size=rd),data=mantel,
              label.colour = "black",
              curvature=nice_curvature(),
              label.fontface=0,
              label.size =4,drop = T,
              # 定义点的边框颜色及填充颜色
              node.colour = c("white", "white"),
              node.fill = c("#984EA3","#5785C1"),
              node.size = c(4.5,4),   # 定义点的大小
              node.shape=c(23,21))+   # 定义点的形状
  # 添加对角线文本
  geom_diag_label(angle=0,parse = T,size=3,color="black",
                  nudge_y =-0.1,nudge_x =-0.4,hjust=0)+
  coord_cartesian(clip="off")+
  scale_size_manual(values = c(0.3,0.5,1,1.5))+
  scale_colour_manual(values = c("#984EA3", "#4DAF4A", "grey","red"))+
  guides(fill = guide_colorbar(title = "pearson's r",order = 1),
         color = guide_legend(title = "*P* value",order = 2,
                              theme = theme(legend.title = 
                                              element_markdown(color="black"))),
         size = guide_legend(title = "pearson's r",order = 3))+
  theme(plot.margin = unit(c(0.5,0,0.5,0.5),units="cm"),
        axis.text.y=element_blank(),
        axis.text.x=element_markdown(size=9,color=c(rep("black",13),"white")),
        axis.ticks =element_blank(),
        panel.background = element_blank(), 
        legend.key = element_blank(), 
        legend.background = element_blank())

bt_plot
ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/髓系细胞亚型/相关性分析/相关性热图.pdf",bt_plot,width=6,height =6)






