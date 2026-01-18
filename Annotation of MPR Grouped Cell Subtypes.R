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

T_cells_object<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/T_cells_object.rds")
T_cells_object = UpdateSeuratObject(T_cells_object)
T_cells_object
table(T_cells_object$RNA_snn_res.0.8)
View(T_cells_object)

saveRDS(T_cells_object,"D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/T_cells_object.rds")


##################################################################################
#细胞亚群微生物####################################################################################
######################################################################################
###读取数据
library(Seurat)
library(tidyverse)
library(harmony)
library(dplyr)
library(patchwork)
library(SeuratWrappers)
setwd("./OSCC_data/MPR_NMPR_group")

MPRNMPR_object_miMPRobe<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/MPRNMPR_object_miMPRobe_v1_harmony.rds")
MPRNMPR_object_miMPRobe = UpdateSeuratObject(MPRNMPR_object_miMPRobe)
View(MPRNMPR_object_miMPRobe)
Idents(MPRNMPR_object_miMPRobe)

DimPlot(MPRNMPR_object_miMPRobe, reduction = 'umap.harmony', group.by = 'cell_type_new',
        label = TRUE,raster=FALSE) + NoLegend()

####seurat个性化细胞注释并把细分亚群放回总群

setwd("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释")

# colnames(MPRNMPR_object_miMPRobe@meta.data)
#
# MPRNMPR_object_miMPRobe@meta.data<-MPRNMPR_object_miMPRobe@meta.data[,-c(58:570)]

# saveRDS(MPRNMPR_object_miMPRobe,"D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/MPRNMPR_object_miMPRobe_v1_harmony.rds")

####提取t cell 亚群#####
T_cells_object=subset(MPRNMPR_object_miMPRobe,idents = 'T cells')
colnames(T_cells_object@meta.data)
pbmc<- T_cells_object
pbmc <- NormalizeData(pbmc,normalization.method = "LogNormalize",
                      scale.factor = 1e4)
GetAssay(pbmc,assay = "RNA")
pbmc <- FindVariableFeatures(pbmc,
                             selection.method = "vst",nfeatures = 5000)
pbmc<- ScaleData(pbmc)
pbmc <- RunPCA(object = pbmc,pc.genes = VariableFeatures(pbmc))
pbmc <- FindNeighbors(pbmc,dims = 1:20)
for (res in c(0.01, 0.05, 0.1, 0.2 ,0.3, 0.5, 0.6,0.8, 1)) {
  pbmc = FindClusters(pbmc,
                      resolution = res, algorithm = 1)

}
table(pbmc$RNA_snn_res.0.8)
sel.clust <- "RNA_snn_res.0.8"
pbmc <- SetIdent(pbmc,value = sel.clust)
metadata <- pbmc@meta.data

# cell_cluster <- data.frame(cell_ID=rownames(metadata), cluster_ID=metadata$RNA_snn_res.0.8) #将聚类结果另存为data.frame
#
# write.csv(cell_cluster,'D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/tcell_cluster.csv',row.names = F, quote = F)

pbmc <- RunUMAP(pbmc,dims = 1:20)#UMAP二次降维，取前20个维度
plot3 = DimPlot(pbmc, reduction = "umap", label=T) #查看clusters在UMAP降维图中的分布
plot3
plot4 = DimPlot(pbmc, reduction = "umap", group.by='orig.ident') #查看每个样本在UMAP降维图中的分布
plot5 = DimPlot(pbmc, reduction = "umap", split.by='orig.ident') #查看每个样本在UMAP降维图中的分面图
plotc <- plot4+plot5
plotc
#harmony分析***
pbmc <- RunHarmony(object =pbmc, reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")#harmony需要pca降维的结果，而scRNA对象此前已经做过PCA降维了，所以此处不再重复
pbmc <- RunUMAP(pbmc, reduction = "harmony", dims = 1:20,reduction.name = "umap.harmony")
pbmc <- RunTSNE(object = pbmc,reduction = "harmony", reduction.key = "tsneharmony_",dims.use = 1:20,reduction.name = "tsne.harmony")
pbmc <- FindNeighbors(pbmc, reduction = "umap.harmony", dims = 1:20)
pbmc <- FindClusters(pbmc)#标准聚类

#绘图和保存
p1 <- DimPlot(pbmc, group.by = "orig.ident",reduction = "umap.harmony", pt.size=0.1) +
  ggtitle("Integrated by harmony")#去批次后的可视化
p1
p2 <- DimPlot(pbmc, group.by = "RNA_snn_res.0.8",reduction = "umap.harmony", pt.size=0.1,label = T) +
  ggtitle("Integrated by harmony")#去批次后的可视化
p2

p3 <- DimPlot(pbmc, group.by = "RNA_snn_res.0.8",reduction = "tsne.harmony", pt.size=0.1,label = T) +
  ggtitle("Integrated by harmony")#去批次后的可视化
p3


DefaultAssay(pbmc) <- "RNA"
markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
all.markers =markers %>% dplyr::select(gene, everything()) %>% subset(p_val<0.05)
top20 =all.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
Idents(pbmc) <- pbmc$RNA_snn_res.0.8
p_heatmap<-DoHeatmap(pbmc,features = top20$gene)+ #热图复现
  theme(text = element_text(size = 5))

write.csv(top20,"tcell_marker_top20.csv")
ggsave("tcell_topmarkergene_heatmap.jpg",p_heatmap, width = 10,height = 15)

T_cells_object<-pbmc
#保存RDS
saveRDS(T_cells_object,"D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/T_cells_object.rds")


setwd("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释")
tcell_marker_top20<- read.csv("./cd4_cd8_markers.csv")
View(tcell_marker_top20)
marker_all<- read.csv("Cell_marker_All.csv")
View(marker_all)
marker_all_tcell<- marker_all[which(marker_all$marker %in% tcell_marker_top20$gene),]
marker_all_tcell2<- marker_all_tcell[match(tcell_marker_top20$gene,marker_all_tcell$marker),]

View(marker_all_tcell2)
write.csv(marker_all_tcell2,"marker_all_tcell_cd4_cd8.csv")

Patients_color <- c("#40e0d0","#ee82ee","#7ccd7c","#551a8b","#ffc0cb","#ff6347",
                    "#f4a460","#b2dfee","#00ff7f","#63b8ff")
library(randomcoloR)
palette <- distinctColorPalette(10)
##调整legend 顺序
MPRNMPR.umap.harmonyredu$Patients <- factor(MPRNMPR.umap.harmonyredu$Patients, levels=c('P1', 'P2', 'P3', 'P4',
                                                                                        'P5',"P6","P7","P8","P9","P10"))
DimPlot(pbmc, group.by = "MPR_Response2",reduction = "umap.harmony", pt.size=0.1) + scale_color_manual(values = palette)+
  ggtitle("MPR_Response2")#去批次后的可视化+

pbmc@meta.data$Patients
#########
####T cell 亚群注释----
####修改高可变基因----
top10<- head(VariableFeatures(T_cells_object),10)
plot1<-VariableFeaturePlot(T_cells_object)
plot2<-LabelPoints(plot = plot1,points=top10,repel=TRUE)
pbmc<- T_cells_object
View(pbmc)
pbmc@assays[["RNA"]]@var.features

var.feature<- read.csv("T_cells_object_var.features.csv",header = T)
pbmc@assays[["RNA"]]@var.features<-var.feature$gene
View(pbmc)



Idents(T_cells_object)<-"RNA_snn_res.0.8"
T_cells_object@meta.data$CRNCR_Response2
###cd4+ t cell###----
####拆分cluster----
FeaturePlot(T_cells_object,features ="CD4",reduction = "umap.harmony",raster=FALSE)
table(T_cells_object$RNA_snn_res.0.8)
sub_object_17 = T_cells_object[,Idents(T_cells_object)== "17"]
cluster20<- subset(myeloid_cells_object,idents = c("20"))
T_cells_object$RNA_snn_res.0.8 = as.character(Idents(T_cells_object))
table(T_cells_object$RNA_snn_res.0.8)
sub_object_17$RNA_snn_res.0.8 = ifelse(sub_object_17@assays$RNA@counts['CD4',]>0,'17_cd4','17_neg')
T_cells_object$RNA_snn_res.0.8[match(colnames(sub_object_17),colnames(T_cells_object))]=sub_object_17$RNA_snn_res.0.8
table(T_cells_object$RNA_snn_res.0.8)
Idents(T_cells_object)<-"RNA_snn_res.0.8"

sub_object_6 = T_cells_object[,Idents(T_cells_object)== "6"]
T_cells_object$RNA_snn_res.0.8 = as.character(Idents(T_cells_object))
table(T_cells_object$RNA_snn_res.0.8)
sub_object_6$RNA_snn_res.0.8 = ifelse(sub_object_6@assays$RNA@counts['CD4',]>0,'6_cd4','6_neg')
T_cells_object$RNA_snn_res.0.8[match(colnames(sub_object_6),colnames(T_cells_object))]=sub_object_6$RNA_snn_res.0.8
table(T_cells_object$RNA_snn_res.0.8)
Idents(T_cells_object)<-"RNA_snn_res.0.8"

sub_object_19 = T_cells_object[,Idents(T_cells_object)== "19"]
T_cells_object$RNA_snn_res.0.8 = as.character(Idents(T_cells_object))
table(T_cells_object$RNA_snn_res.0.8)
sub_object_19$RNA_snn_res.0.8 = ifelse(sub_object_19@assays$RNA@counts['CD4',]>0,'19_cd4','19_neg')
T_cells_object$RNA_snn_res.0.8[match(colnames(sub_object_19),colnames(T_cells_object))]=sub_object_19$RNA_snn_res.0.8
table(T_cells_object$RNA_snn_res.0.8)
Idents(T_cells_object)<-"RNA_snn_res.0.8"

sub_object_10 = T_cells_object[,Idents(T_cells_object)== "10"]
T_cells_object$RNA_snn_res.0.8 = as.character(Idents(T_cells_object))
table(T_cells_object$RNA_snn_res.0.8)
sub_object_10$RNA_snn_res.0.8 = ifelse(sub_object_10@assays$RNA@counts['CD4',]>0,'10_cd4','10_neg')
T_cells_object$RNA_snn_res.0.8[match(colnames(sub_object_10),colnames(T_cells_object))]=sub_object_10$RNA_snn_res.0.8
table(T_cells_object$RNA_snn_res.0.8)
Idents(T_cells_object)<-"RNA_snn_res.0.8"


sub_object_2 = T_cells_object[,Idents(T_cells_object)== "2"]
T_cells_object$RNA_snn_res.0.8 = as.character(Idents(T_cells_object))
table(T_cells_object$RNA_snn_res.0.8)
sub_object_2$RNA_snn_res.0.8 = ifelse(sub_object_2@assays$RNA@counts['CD4',]>0,'2_cd4','2_neg')
T_cells_object$RNA_snn_res.0.8[match(colnames(sub_object_2),colnames(T_cells_object))]=sub_object_2$RNA_snn_res.0.8
table(T_cells_object$RNA_snn_res.0.8)
Idents(T_cells_object)<-"RNA_snn_res.0.8"

sub_object_3 = T_cells_object[,Idents(T_cells_object)== "3"]
T_cells_object$RNA_snn_res.0.8 = as.character(Idents(T_cells_object))
table(T_cells_object$RNA_snn_res.0.8)
sub_object_3$RNA_snn_res.0.8 = ifelse(sub_object_3@assays$RNA@counts['CD4',]>0,'3_cd4','3_neg')
T_cells_object$RNA_snn_res.0.8[match(colnames(sub_object_3),colnames(T_cells_object))]=sub_object_3$RNA_snn_res.0.8
table(T_cells_object$RNA_snn_res.0.8)
Idents(T_cells_object)<-"RNA_snn_res.0.8"


sub_object_8 = T_cells_object[,Idents(T_cells_object)== "8"]
T_cells_object$RNA_snn_res.0.8 = as.character(Idents(T_cells_object))
table(T_cells_object$RNA_snn_res.0.8)
sub_object_8$RNA_snn_res.0.8 = ifelse(sub_object_8@assays$RNA@counts['CD4',]>0,'8_cd4','8_neg')
T_cells_object$RNA_snn_res.0.8[match(colnames(sub_object_8),colnames(T_cells_object))]=sub_object_8$RNA_snn_res.0.8
table(T_cells_object$RNA_snn_res.0.8)
Idents(T_cells_object)<-"RNA_snn_res.0.8"

sub_object_9 = T_cells_object[,Idents(T_cells_object)== "9"]
T_cells_object$RNA_snn_res.0.8 = as.character(Idents(T_cells_object))
table(T_cells_object$RNA_snn_res.0.8)
sub_object_9$RNA_snn_res.0.8 = ifelse(sub_object_9@assays$RNA@counts['CD4',]>0,'9_cd4','9_neg')
T_cells_object$RNA_snn_res.0.8[match(colnames(sub_object_9),colnames(T_cells_object))]=sub_object_9$RNA_snn_res.0.8
table(T_cells_object$RNA_snn_res.0.8)
Idents(T_cells_object)<-"RNA_snn_res.0.8"

sub_object_11 = T_cells_object[,Idents(T_cells_object)== "11"]
T_cells_object$RNA_snn_res.0.8 = as.character(Idents(T_cells_object))
table(T_cells_object$RNA_snn_res.0.8)
sub_object_11$RNA_snn_res.0.8 = ifelse(sub_object_11@assays$RNA@counts['CD4',]>0,'11_cd4','11_neg')
T_cells_object$RNA_snn_res.0.8[match(colnames(sub_object_11),colnames(T_cells_object))]=sub_object_11$RNA_snn_res.0.8
table(T_cells_object$RNA_snn_res.0.8)
Idents(T_cells_object)<-"RNA_snn_res.0.8"


sub_object_12 = T_cells_object[,Idents(T_cells_object)== "12"]
T_cells_object$RNA_snn_res.0.8 = as.character(Idents(T_cells_object))
table(T_cells_object$RNA_snn_res.0.8)
sub_object_12$RNA_snn_res.0.8 = ifelse(sub_object_12@assays$RNA@counts['CD4',]>0,'12_cd4','12_neg')
T_cells_object$RNA_snn_res.0.8[match(colnames(sub_object_12),colnames(T_cells_object))]=sub_object_12$RNA_snn_res.0.8
table(T_cells_object$RNA_snn_res.0.8)
Idents(T_cells_object)<-"RNA_snn_res.0.8"


sub_object_18 = T_cells_object[,Idents(T_cells_object)== "18"]
T_cells_object$RNA_snn_res.0.8 = as.character(Idents(T_cells_object))
table(T_cells_object$RNA_snn_res.0.8)
sub_object_18$RNA_snn_res.0.8 = ifelse(sub_object_18@assays$RNA@counts['CD4',]>0,'18_cd4','18_neg')
T_cells_object$RNA_snn_res.0.8[match(colnames(sub_object_18),colnames(T_cells_object))]=sub_object_18$RNA_snn_res.0.8
table(T_cells_object$RNA_snn_res.0.8)
Idents(T_cells_object)<-"RNA_snn_res.0.8"


sub_object_20 = T_cells_object[,Idents(T_cells_object)== "20"]
T_cells_object$RNA_snn_res.0.8 = as.character(Idents(T_cells_object))
table(T_cells_object$RNA_snn_res.0.8)
sub_object_20$RNA_snn_res.0.8 = ifelse(sub_object_20@assays$RNA@counts['CD4',]>0,'20_cd4','20_neg')
T_cells_object$RNA_snn_res.0.8[match(colnames(sub_object_20),colnames(T_cells_object))]=sub_object_20$RNA_snn_res.0.8
table(T_cells_object$RNA_snn_res.0.8)
Idents(T_cells_object)<-"RNA_snn_res.0.8"


sub_object = T_cells_object[,Idents(T_cells_object)== "15"]

FeaturePlot(sub_object,features =c("CD8A"),reduction = "umap.harmony",raster=FALSE)

####CD8+全部----
#### 0, 1,10_cd4,10_neg,11_cd4,11_neg,12_cd4,12_neg,
#### 13,14,15,16,17_cd4,17_neg,18_cd4,18_neg,
#### 19_cd4,19_neg,2_cd4,2_neg,20_cd4,20_neg,3_cd4,
#### 3_neg,4 ,5,6_cd4,6_neg,7,8_cd4,8_neg,9_cd4,9_neg


####cluster 1,11_neg,13,17_neg,6_neg

###18_neg 不表达cd4 cd8  19_neg,20_neg

table(T_cells_object1$RNA_snn_res.0.8)
DimPlot(T_cells_object, group.by = "RNA_snn_res.0.8",reduction = "umap.harmony", pt.size=0.1,label = T) +
  ggtitle("Integrated by harmony")
sub_object<- subset(T_cells_object,idents=c("11_cd8"))

FeaturePlot(T_cells_object,features =c("CD4","CD8B","CD8A"),reduction = "umap.harmony",raster=FALSE)

T_cells_object$RNA_snn_res.0.8 = as.character(Idents(T_cells_object))
table(T_cells_object$RNA_snn_res.0.8)
sub_object_10<- subset(T_cells_object,idents=c("10_neg"))
sub_object_10$RNA_snn_res.0.8 = ifelse(sub_object_10@assays$RNA@counts['CD8A',]>0|sub_object_10@assays$RNA@counts['CD8B',]>0,'10_cd8','10_neg')
T_cells_object$RNA_snn_res.0.8[match(colnames(sub_object_10),colnames(T_cells_object))]=sub_object_10$RNA_snn_res.0.8
table(T_cells_object$RNA_snn_res.0.8)
Idents(T_cells_object)<-"RNA_snn_res.0.8"

T_cells_object$RNA_snn_res.0.8 = as.character(Idents(T_cells_object))
table(T_cells_object$RNA_snn_res.0.8)
sub_object_12<- subset(T_cells_object,idents=c("12_neg"))
sub_object_12$RNA_snn_res.0.8 = ifelse(sub_object_12@assays$RNA@counts['CD8A',]>0|sub_object_12@assays$RNA@counts['CD8B',]>0,'12_cd8','12_neg')
T_cells_object$RNA_snn_res.0.8[match(colnames(sub_object_12),colnames(T_cells_object))]=sub_object_12$RNA_snn_res.0.8
table(T_cells_object$RNA_snn_res.0.8)
Idents(T_cells_object)<-"RNA_snn_res.0.8"

T_cells_object$RNA_snn_res.0.8 = as.character(Idents(T_cells_object))
table(T_cells_object$RNA_snn_res.0.8)
sub_object_2<- subset(T_cells_object,idents=c("2_neg"))
sub_object_2$RNA_snn_res.0.8 = ifelse(sub_object_2@assays$RNA@counts['CD8A',]>0|sub_object_2@assays$RNA@counts['CD8B',]>0,'2_cd8','2_neg')
T_cells_object$RNA_snn_res.0.8[match(colnames(sub_object_2),colnames(T_cells_object))]=sub_object_2$RNA_snn_res.0.8
table(T_cells_object$RNA_snn_res.0.8)
Idents(T_cells_object)<-"RNA_snn_res.0.8"


T_cells_object$RNA_snn_res.0.8 = as.character(Idents(T_cells_object))
table(T_cells_object$RNA_snn_res.0.8)
sub_object_3<- subset(T_cells_object,idents=c("3_neg"))
sub_object_3$RNA_snn_res.0.8 = ifelse(sub_object_3@assays$RNA@counts['CD8A',]>0|sub_object_3@assays$RNA@counts['CD8B',]>0,'3_cd8','3_neg')
T_cells_object$RNA_snn_res.0.8[match(colnames(sub_object_3),colnames(T_cells_object))]=sub_object_3$RNA_snn_res.0.8
table(T_cells_object$RNA_snn_res.0.8)
Idents(T_cells_object)<-"RNA_snn_res.0.8"

T_cells_object$RNA_snn_res.0.8 = as.character(Idents(T_cells_object))
table(T_cells_object$RNA_snn_res.0.8)
sub_object_7<- subset(T_cells_object,idents=c("7"))
sub_object_7$RNA_snn_res.0.8 = ifelse(sub_object_7@assays$RNA@counts['CD8A',]>0|sub_object_7@assays$RNA@counts['CD8B',]>0,'7_cd8','7_neg')
T_cells_object$RNA_snn_res.0.8[match(colnames(sub_object_7),colnames(T_cells_object))]=sub_object_7$RNA_snn_res.0.8
table(T_cells_object$RNA_snn_res.0.8)
Idents(T_cells_object)<-"RNA_snn_res.0.8"

T_cells_object$RNA_snn_res.0.8 = as.character(Idents(T_cells_object))
table(T_cells_object$RNA_snn_res.0.8)
sub_object_8<- subset(T_cells_object,idents=c("8_neg"))
sub_object_8$RNA_snn_res.0.8 = ifelse(sub_object_8@assays$RNA@counts['CD8A',]>0|sub_object_8@assays$RNA@counts['CD8B',]>0,'8_cd8','8_neg')
T_cells_object$RNA_snn_res.0.8[match(colnames(sub_object_8),colnames(T_cells_object))]=sub_object_8$RNA_snn_res.0.8
table(T_cells_object$RNA_snn_res.0.8)
Idents(T_cells_object)<-"RNA_snn_res.0.8"

T_cells_object$RNA_snn_res.0.8 = as.character(Idents(T_cells_object))
table(T_cells_object$RNA_snn_res.0.8)
sub_object_9<- subset(T_cells_object,idents=c("9_neg"))
sub_object_9$RNA_snn_res.0.8 = ifelse(sub_object_9@assays$RNA@counts['CD8A',]>0|sub_object_9@assays$RNA@counts['CD8B',]>0,'9_cd8','9_neg')
T_cells_object$RNA_snn_res.0.8[match(colnames(sub_object_9),colnames(T_cells_object))]=sub_object_9$RNA_snn_res.0.8
table(T_cells_object$RNA_snn_res.0.8)
Idents(T_cells_object)<-"RNA_snn_res.0.8"

FeaturePlot(T_cells_object,features =c("CD3E",
                                   "CD3D"),reduction = "umap.harmony",raster=FALSE)

####NK 7_neg
# 0 1 10_cd4 10_cd8 10_neg 11_cd4 11_neg 12_cd4 12_cd8 12_neg  13 14 15
# 16 17_cd4 17_neg 18_cd4 18_neg 19_cd4 19_neg  2_cd4  2_cd8  2_neg 20_cd4 20_neg  3_cd4
# 3_cd8  3_neg  4  5  6_cd4  6_neg  7_cd8  7_neg  8_cd4  8_cd8  8_neg  9_cd4  9_cd8
# 9_neg
VlnPlot(object = T_cells_object, features =c("TNFRSF14",
                                             "CD28",
                                             "ICOS",
                                             "TNFRSF9"
                                             ))
sub_object<- subset(T_cells_object,idents=c("16"))
FeaturePlot(sub_object,features =c("CD4",
                                       "CD8A",
                                       "CD8B"
),reduction = "umap.harmony",raster=FALSE)



####CD8 naive CD8----####cluster 1,11_neg,13,17_neg,6_neg

#### CD8+ effector----

FeaturePlot(T_cells_object,features =c("CX3CR1",
                                       "KLRG1"),reduction = "umap.harmony",raster=FALSE)




####CD4全部##----
##"10_cd4","2_cd4", "3_cd4","12_cd4","18_cd4"

setwd("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释")

CD4<-subset(T_cells_object,idents=c(
                                    ))

sub_object<- subset(T_cells_object,idents=c("CD4"))

FeaturePlot(T_cells_object,features =c("CXCL13",
                                       "PDCD1",
                                       "MAF",
                                       "BCL6",
                                       "TOX2",
                                       "IL6ST",
                                       "FKBP5"
                                       ),reduction = "umap.harmony",raster=FALSE)



markers <- FindAllMarkers(T_cells_object, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
View(markers)
write.csv(markers,"cd4_cd8_markers.csv")

sub_object<-subset(T_cells_object,idents=c("6_cd4","19_cd4"))
FeaturePlot(sub_object,features =c("NKG7",
                                   "KLRB1",
                                   "KLRD1",
                                   "XCL",
                                   "CD4"),reduction = "umap.harmony",raster=FALSE)


####Th1(CD4）----

sub_object<- subset(T_cells_object,idents=c("CCL5"))

VlnPlot(object = T_cells_object, features = "CCL5")

FeaturePlot(sub_object,features =c("TBX21",
                                       "IFNG",
                                       "CCR5",
                                       "CXCR3",
                                       "STAT4",
                                       "IL2",
                                       "IL10",
                                       "TNF",
                                       "CD4"),reduction = "umap.harmony",raster=FALSE)



####Th 2(CD4）----

VlnPlot(object = T_cells_object, features =c("GATA3","CCR8"))

FeaturePlot(T_cells_object,features =c("CCR3",
                                   "GATA3",
                                   "IL4",
                                   "CCR4",
                                   "CCR8",
                                   "CXCR4",
                                   "CD4"),reduction = "umap.harmony",raster=FALSE)


####Th 9(CD4）----
VlnPlot(object = T_cells_object, features =c("IL4R","TGFBR2"))

FeaturePlot(T_cells_object,features =c("IL9","IL4R","TGFBR2",
                                       "CD4"),reduction = "umap.harmony",raster=FALSE)



####Th17(CD4）----
VlnPlot(object = T_cells_object, features =c("IL17A",
                                             "IL17F",
                                             "RORA",
                                             "IL26"
))

FeaturePlot(CD4_Th17,features =c("IL17A",
                                       "IL17F",
                                       "RORA",
                                       "IL26",
                                       "CD4"),reduction = "umap.harmony",raster=FALSE)



####Th22(CD4)----

VlnPlot(object = T_cells_object, features =c("IL22",
                                             "AHR",
                                             "STAT3"
))

FeaturePlot(T_cells_object,features =c("IL22",
                                 "AHR",
                                 "STAT3",
                                 "CD4"),reduction = "umap.harmony",raster=FALSE)




###CD4-treg "FOXP3""IL2RA"
FeaturePlot(T_cells_object,features =c("CTLA4",
                                       "FOXP3",
                                       "ICOS",
                                       "TNFRSF18",
                                       "GITR",
                                       "IKZF2",
                                       "CCM2",
                                       "CD5",
                                       "CD4"),reduction = "umap.harmony",raster=FALSE)
FeaturePlot(T_cells_object,features =c("IL2RA"
                                       ),reduction = "umap.harmony",raster=FALSE)


CD4_treg<-subset(T_cells_object,idents=c("0","16","17_cd4"))

FeaturePlot(CD4_treg,features =c("CTLA4",
                                       "FOXP3",
                                       "ICOS",
                                       "TNFRSF18",
                                       "GITR",
                                       "IKZF2",
                                       "CCM2",
                                       "CD5",
                                 "CD4"),reduction = "umap.harmony",raster=FALSE)

####Tfh----
VlnPlot(object = T_cells_object, features =c("CXCR5",
                                             "ICOS",
                                            "PDCD1"
))


####nkT----
VlnPlot(object = T_cells_object, features =c("NCR3",
                                             "NCR2",
                                             "NCR1"))

####γδT----
VlnPlot(object = T_cells_object, features =c("GZMA",
                                             "Trdv4",
                                             "Trdc",
                                             "TRV9"
))




####NK----
VlnPlot(object = T_cells_object, features =c("NKG7",
                                             "KLRB1",
                                             "KLRD1",
                                             "XCL"
))

####immature NK----

VlnPlot(object = T_cells_object, features =c("KLRC1"
))

####mature NK

VlnPlot(object = T_cells_object, features =c("FCGR3A",
                                             "PRF1",
                                             "FCGR3A",
                                             "KLRC2"
))

####NK CD16----

VlnPlot(object = T_cells_object, features =c("CX3CR1",
                                             "KLRG1"
))


####Resident----

VlnPlot(object = T_cells_object, features =c("CD69",
                                             "RUNX3",
                                             "NR4A1"

))



####CD4 naive CD8 #####
gene<- c("TCF7",
         "SELL",
         "LEF1",
         "CCR7",
         "CD4")

VlnPlot(object = T_cells_object, features = as.character(gene))

VlnPlot(object = T_cells_object, features = "TCF7")

FeaturePlot(T_cells_object,features =c("TCF7",
                                 "SELL",
                                 "LEF1",
                                 "CCR7",
                                 "CD4"),reduction = "umap.harmony",raster=FALSE)



CD4_naive<-subset(T_cells_object,idents=c("14","15","10_cd4","2_cd4",
                                          "3_cd4","12_cd4","18_cd4","9_cd4"))
CD4_naive<-subset(T_cells_object,idents=c("10_cd4","2_cd4",
                                          "3_cd4","12_cd4","18_cd4"))
CD4_naive<-subset(T_cells_object,idents=c("19_cd4"))
FeaturePlot(CD4_naive,features =c("CD2",
                                       "SELL",
                                       "LEF1",
                                       "CCR7",
                                       "CD4"),reduction = "umap.harmony",raster=FALSE)

####effector CD4 cytotoxic CD8/ ----

CD4_effector<-subset(T_cells_object,idents=c("14","15","10_cd4","2_cd4",
                                          "3_cd4","12_cd4","18_cd4","9_cd4"))
CD4_effector<-subset(T_cells_object,idents=c("10_cd4","2_cd4",
                                          "3_cd4","12_cd4","18_cd4"))
CD4_effector<-subset(T_cells_object,idents=c("17_cd4"))

VlnPlot(object = T_cells_object, features = "GNLY")

FeaturePlot(T_cells_object,features =c("IFNG",
                                     "CCL5",
                                     "GZMB",
                                     "PRF1",
                                     "GNLY",
                                     "CD4"),reduction = "umap.harmony",raster=FALSE)



####exhaustedCD4/CD8----

VlnPlot(object = T_cells_object, features = "FKBP5")

sub_object<-subset(T_cells_object,idents=c("6_cd4","19_cd4"))
FeaturePlot(T_cells_object,features =c("PDCD1",
                                       "LAG3",
                                       "CD244",
                                       "CD160",
                                       "HAVCR2",
                                       "CD4"),reduction = "umap.harmony",raster=FALSE)



####memory CD4/CD8----
###"6_cd4","9_cd4"

CD4_memory<-subset(T_cells_object,idents=c("16","0","17_cd4"))

VlnPlot(object = T_cells_object, features = "LTB")


FeaturePlot(T_cells_object,features =c("ARHGAP45",
                                       "MS4A6B",
                                       "ARHGAP15",
                                       "MIR142HG",
                                       "GIMAP1",
                                       "PASK",
                                       "GPR183",
                                       "LTB"),reduction = "umap.harmony",raster=FALSE)


####Treg（CD4）----

VlnPlot(object = T_cells_object, features = "PMAIP1")

FeaturePlot(T_cells_object,features =c("FOXP3",
                                       "IKZF2",
                                       "CCM2",
                                       "CD5",
                                       "IL2RA",
                                       "CTLA4",
                                       "TNFRSF4",
                                       "PMAIP1"),reduction = "umap.harmony",raster=FALSE)





DimPlot(T_cells_object,group.by="cd4_cd8",reduction = "umap.harmony",label = T,raster=FALSE)

ggsave("tcell_dimplot.jpg",tcell_simplot, width = 10,height = 10)

DotPlot(sce_T,features = genes_to_check)+coord_flip()
dev.off()

FeaturePlot(T_cells_object,features ="PD1",reduction = "umap.harmony",raster=FALSE)

FeaturePlot(T_cells_object,features ="CD4",reduction = "umap.harmony",raster=FALSE)

FeaturePlot(T_cells_object,features =c("ANXA1","LMNA","GPR183","CCR7"),reduction = "umap.harmony",raster=FALSE)


celltype=data.frame(ClusterID=0:4,
                    celltype=0:4)
celltype[celltype$ClusterID %in% c(0),2]="TCF1hi"
celltype[celltype$ClusterID %in% c(1,2),2]="T cyto"
celltype[celltype$ClusterID %in% c(3),2]="Cycling 1"
celltype[celltype$ClusterID %in% c(4),2]="Cycling 2"
for (i in 1:nrow(celltype)) {
  sce_T@meta.data[which(sce_T@meta.data$RNA_snn_res.0.8==celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]

}
Idents(sce_T) <- sce_T$celltype
DimPlot(sce_T,reduction = "umap",label = T)
ggsave("dimplot.pdf")

###########################################################
############################################################
###Tcell marker 热图----

setwd("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/marker热图分型")

# 检查基因名称是否存在
genes_to_check <- c("CD8A", "CD8B", "CD4", "CD3E", "CD3D", "CD3G")
missing_genes <- genes_to_check[!(genes_to_check %in% rownames(T_cells_object@assays$RNA@counts))]

if (length(missing_genes) > 0) {
  stop(paste("Missing genes in the dataset:", paste(missing_genes, collapse = ", ")))
} else {
  # 如果基因存在，执行ifelse操作
  T_cells_object$cd4_cd8_group <- ifelse(T_cells_object@assays$RNA@counts['CD8A', ] > 0 | T_cells_object@assays$RNA@counts['CD8B', ] > 0, "CD8+",
                                         ifelse(T_cells_object@assays$RNA@counts['CD4', ] > 0, "CD4+",
                                                ifelse((T_cells_object@assays$RNA@counts['CD3E', ] > 0 | T_cells_object@assays$RNA@counts['CD3D', ] > 0 |
                                                          T_cells_object@assays$RNA@counts['CD3G', ] > 0) & T_cells_object@assays$RNA@counts['CD4', ] == 0 &
                                                         T_cells_object@assays$RNA@counts['CD8A', ] == 0 & T_cells_object@assays$RNA@counts['CD8B', ] == 0, "CD3+CD4-CD8-",
                                                       ifelse(T_cells_object@assays$RNA@counts['CD3E', ] == 0 | T_cells_object@assays$RNA@counts['CD3D', ] == 0 |
                                                                T_cells_object@assays$RNA@counts['CD3G', ] == 0, "NK", NA))))
}


markerdata <- ScaleData(T_cells_object, features = as.character(unique(tcell_marker$gene_name)), assay = "RNA")

tcell_marker<- read.csv("T_NK_marker.csv")
markers <-tcell_marker$gene_name
markers <- as.data.frame(markers)

aver_dt<-AverageExpression(markerdata,
                           features = as.character(unique(tcell_marker$gene_name)),
                           group.by =c('cd4_cd8_group'))


aver_dt_df <- as.data.frame(aver_dt)
View(aver_dt_df)
aver_dt_df$marker<- rownames(aver_dt_df)
aver_dt_df<- aver_dt_df[,-5]
colnames(aver_dt_df)<- c("CD3+CD4-CD8-","CD4+","CD8+","NK","marker")
unique_markers <- unique(aver_dt_df$marker)
aver_dt_df$marker <- factor(aver_dt_df$marker, levels = unique_markers)

aver_dt_df_long <- pivot_longer(aver_dt_df, cols = c("CD3+CD4-CD8-","CD4+","CD8+","NK"), names_to = "Tcell_group", values_to = "value")

library(RColorBrewer)
library(MetBrewer)

aver_dt<-AverageExpression(T_cells_object,
                           features = as.character(unique(tcell_marker$gene_name)),
                           group.by =c('RNA_snn_res.0.8'))


aver_dt_df <- as.data.frame(aver_dt)
View(aver_dt_df)
aver_dt_df$marker<- rownames(aver_dt_df)
unique_markers <- unique(aver_dt_df$marker)

aver_dt_df_long <- pivot_longer(aver_dt_df, cols = starts_with("RNA"), names_to = "Tcell_group", values_to = "value")
aver_dt_df_long$marker <- factor(aver_dt_df$marker, levels = rev(unique_markers))

aver_dt_df_long$Tcell_group<- factor(aver_dt_df_long$Tcell_group, levels = c("RNA.g0", "RNA.g10.cd4",
                                                                        "RNA.g11.cd4", "RNA.g12.cd4", "RNA.g13","RNA.g14",
                                                                        "RNA.g15","RNA.g16",
                                                                        "RNA.g17.cd4","RNA.g18.cd4",
                                                                        "RNA.g19.cd4","RNA.g2.cd4",
                                                                        "RNA.g20.cd4","RNA.g3.cd4",
                                                                        "RNA.g4","RNA.g5","RNA.g6.cd4",
                                                                        "RNA.g8.cd4","RNA.g9.cd4","RNA.g1" ,"RNA.g10.cd8",
                                                                        "RNA.g11.neg","RNA.g12.cd8","RNA.g17.neg","RNA.g18.neg",
                                                                        "RNA.g19.neg","RNA.g2.cd8","RNA.g20.neg","RNA.g3.cd8",
                                                                        "RNA.g6.neg","RNA.g7.cd8", "RNA.g8.cd8", "RNA.g9.cd8",
                                                                        "RNA.g10.neg" ,
                                                                        "RNA.g12.neg",
                                                                        "RNA.g2.neg",
                                                                        "RNA.g3.neg",
                                                                        "RNA.g8.neg" ,"RNA.g9.neg",
                                                                        "RNA.g7.neg"))

library(RColorBrewer)
library(MetBrewer)

cd4_group<- c("RNA.g0", "RNA.g10.cd4","RNA.g11.cd4", "RNA.g12.cd4", "RNA.g13","RNA.g14",
                                           "RNA.g15","RNA.g16",
                                             "RNA.g17.cd4","RNA.g18.cd4",
                                              "RNA.g19.cd4","RNA.g2.cd4",
                                               "RNA.g20.cd4","RNA.g3.cd4",
                                                  "RNA.g4","RNA.g5","RNA.g6.cd4",
                                                 "RNA.g8.cd4","RNA.g9.cd4","RNA.g3.cd4")
aver_dt_df_long_cd4 <- aver_dt_df_long[aver_dt_df_long$Tcell_group %in% cd4_group, ]
tcell_marker_plot<-ggplot(aver_dt_df_long_cd4, aes(Tcell_group,marker)) +
  geom_tile(aes(fill =value), colour = "white", linewidth = 0.1)+
  scale_fill_gradientn(name='T cell marker genes mean\nscaled expression',
                       colours=c("white","blue","yellow","red"))+
  theme_minimal() +
  theme(axis.line=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = 10,color = 'black',
        ),
        axis.text.x = element_text(size = 10,color = 'black',angle =40,hjust =1,vjust=1),
        legend.title = element_text(size=12, color = "black"),
        legend.text = element_text(size=10,color = "black",angle =45),
        legend.position = "top") +
  scale_y_discrete(position = "right")

library(pheatmap)
aver_dt_df_wider_cd4 <- aver_dt_df_long_cd4 %>%
  pivot_wider(names_from =Tcell_group , values_from = value)

aver_dt_df_wider_cd4 <-t(aver_dt_df_wider_cd4)
aver_dt_df_wider_cd4<-aver_dt_df_wider_cd4[-1,]
aver_dt_df_wider_cd4<-aver_dt_df_wider_cd4 %>% as.data.frame()
numeric_cols <- aver_dt_df_wider_cd4 %>%
mutate_all(as.numeric)
content <- scale(numeric_cols)
pheatmap(content, cellwidth = 20, cellheight = 20,
         main = "heatmap",filename = "heatmap_cd4.png")

content <- scale(content)#缩放数据
content <- t(content)
gg <- hclust(dist(content))#对行聚类
zz <- hclust(dist(t(content))) #对列聚类
content <- content[gg$order,]#行，按照聚类结果排序


library("ggtree")

v <- ggtree(zz)+layout_dendrogram()# 绘制行聚类树
#####热图和聚类树拼在一起#####
library("aplot")
tcell_marker_plot%>% insert_top(v,height = 0.1)# 使用 aplot包里的函数进行拼图
ggsave(file = "heatmap_cd4.png",width = 15, height = 5,dpi=600)


cd8_group<- c("RNA.g1" ,"RNA.g10.cd8",
              "RNA.g11.neg","RNA.g12.cd8","RNA.g17.neg","RNA.g18.neg",
              "RNA.g19.neg","RNA.g2.cd8","RNA.g20.neg","RNA.g3.cd8",
              "RNA.g6.neg","RNA.g7.cd8", "RNA.g8.cd8", "RNA.g9.cd8")
aver_dt_df_long_cd8 <- aver_dt_df_long[aver_dt_df_long$Tcell_group %in% cd8_group, ]


library(pheatmap)
aver_dt_df_wider_cd8 <- aver_dt_df_long_cd8 %>%
  pivot_wider(names_from =Tcell_group , values_from = value)

aver_dt_df_wider_cd8 <-t(aver_dt_df_wider_cd8)
aver_dt_df_wider_cd8<-aver_dt_df_wider_cd8[-1,]
aver_dt_df_wider_cd8<-aver_dt_df_wider_cd8 %>% as.data.frame()
numeric_cols <- aver_dt_df_wider_cd8 %>%
  mutate_all(as.numeric)
content <- scale(numeric_cols)
pheatmap(content, cellwidth = 20, cellheight = 20,
         main = "heatmap",filename = "heatmap1.png")

####T cell marker正式画热图一起#####----


Idents(T_cells_object)<-"cd4_cd8_group"
unique(Idents(T_cells_object))
markers <- FindAllMarkers(T_cells_object, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

write.csv(markers,"T_cell_markers.csv")
markers<- read.csv("T_cell_markers.csv")

marker_all<- read.csv("Cell_marker_All.csv")
marker_all_tcell<- marker_all[match(markers$gene,marker_all$marker),]

write.csv(marker_all_tcell,"marker_all_tcell_match.csv")

FindMarkers
Idents(T_cells_object) <-"cd4_cd8_group"
unique(Idents(T_cells_object))
subcd8_object<- subset(T_cells_object,idents=c("CD8_C1-RGS1","CD8_C2-LYST","CD8_C3-IL7R","CD8_C4-GZMK","CD8_C5-MKI67"))
deg_all= FindAllMarkers(subcd8_object,only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(deg_all,"CD8_marker_deg_all.csv")
deg_all=FindMarkers(subcd8_object, ident.1 = c("CD8_C2-LYST"))
write.csv(deg_all,"CD8_C2-LYST2.csv")


subcd4_object<- subset(T_cells_object,idents=c("CD4_C1-PDCD1","CD4_C2-FOXP3","CD4_C3-ITGA2","CD4_C4-IL17F"))
deg_all= FindAllMarkers(subcd4_object,only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(deg_all,"CD4_marker_deg_all.csv")
deg_all=FindMarkers(subcd4_object, ident.1 = c("CD4_C3-ITGA2"))
write.csv(deg_all,"CD4_C3-ITGA2_deg.csv")


setwd("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/marker热图分型")

old_to_new_groups <- c("17_cd4" = "CD4_C1-TNFRSF4",
                       "0" = "CD4_C1-TNFRSF4",
                       "16" = "CD4_C1-TNFRSF4",
                      "20_cd4" = "CD4_C1-TNFRSF4",
                      "11_cd4" = "CD4_C2-GZMB",
                       "14" = "CD4_C3-IL17F",
                       "2_cd4" = "CD4_C3-IL17F",
                       "15" = "CD4_C3-IL17F",
                       "18_cd4" = "CD4_C3-IL17F",
                       "10_cd4" = "CD4_C3-IL17F",
                       "12_cd4" = "CD4_C3-IL17F",
                       "9_cd4" = "CD4_C3-IL17F",
                       "3_cd4" ="CD4_C3-IL17F",
                       "19_cd4" = "CD4_C4-Unknow",
                       "6_cd4" = "CD4_C4-Unknow",
                       "11_neg" = "CD8_C1-LAG3",
                       "17_neg" = "CD8_C2-MKI67",
                       "18_neg" = "CD8_C3-SELL",
                       "2_cd8" = "CD8_C3-SELL",
                       "19_neg" = "CD8_C4-Unknow",
                       "6_neg" = "CD8_C4-Unknow",
                       "7_cd8" = "CD8_C5-Unknow",
                       "12_cd8" = "CD8_C5-Unknow",
                       "3_cd8" = "CD8_C5-Unknow",
                       "1" = "CD8_C6-Unknow",
                       "8_cd8" = "CD8_C6-Unknow",
                       "10_cd8" = "CD8_C6-Unknow",
                       "9_cd8" = "CD8_C6-Unknow",
                      "8_cd4" = "CD8_C6-Unknow",
                      "5" = "CD8_C6-Unknow",
                      "13" = "CD8_C6-Unknow",
                      "4" = "CD8_C6-Unknow",
                       "7_neg" = "NK",
                       "2_neg" ="CD3+CD4-CD8-",
                       "12_neg"="CD3+CD4-CD8-",
                       "10_neg"="CD3+CD4-CD8-",
                       "8_neg" ="CD3+CD4-CD8-",
                       "9_neg" ="CD3+CD4-CD8-",
                       "3_neg" ="CD3+CD4-CD8-",
                      "20_neg" = "CD3+CD4-CD8-")

old_to_new_groups2 <- c("17_cd4" = "CD4_C3-FOXP3",
                       "0" = "CD4_C3-FOXP3",
                       "16" = "CD4_C3-FOXP3",
                       "20_cd4" = "CD4_C3-FOXP3",
                       "11_cd4" = "CD4_C1-PDCD1",
                       "14" = "CD4_C2-IL7R",
                       "2_cd4" = "CD4_C2-IL7R",
                       "15" = "CD4_C2-IL7R",
                       "18_cd4" = "CD4_C2-IL7R",
                       "10_cd4" = "CD4_C2-IL7R",
                       "12_cd4" = "CD4_C2-IL7R",
                       "9_cd4" = "CD4_C2-IL7R",
                       "3_cd4" ="CD4_C2-IL7R",
                       "19_cd4" = "CD4_C4-FOSB",
                       "6_cd4" = "CD4_C4-FOSB",
                       "11_neg" = "CD8_C1-GNLY",
                       "17_neg" = "CD8_C5-STMN1",
                       "18_neg" = "CD8_C3-IL7R",
                       "2_cd8" = "CD8_C3-IL7R",
                       "19_neg" = "CD8_C2-CTSW",
                       "6_neg" = "CD8_C2-CTSW",
                       "7_cd8" = "CD8_C4-GZMK",
                       "12_cd8" = "CD8_C4-GZMK",
                       "3_cd8" = "CD8_C4-GZMK",
                       "1" = "CD8_C4-GZMK",
                       "8_cd8" = "CD8_C4-GZMK",
                       "10_cd8" = "CD8_C4-GZMK",
                       "9_cd8" = "CD8_C4-GZMK",
                       "8_cd4" = "CD8_C4-GZMK",
                       "5" = "CD8_C4-GZMK",
                       "13" = "CD8_C4-GZMK",
                       "4" = "CD8_C4-GZMK",
                       "7_neg" = "NK",
                       "2_neg" ="CD3+CD4-CD8-",
                       "12_neg"="CD3+CD4-CD8-",
                       "10_neg"="CD3+CD4-CD8-",
                       "8_neg" ="CD3+CD4-CD8-",
                       "9_neg" ="CD3+CD4-CD8-",
                       "3_neg" ="CD3+CD4-CD8-",
                       "20_neg" = "CD3+CD4-CD8-")


T_cells_object@meta.data$cd4_cd8_group <- plyr::mapvalues(T_cells_object@meta.data$RNA_snn_res.0.8,
                                                          from = names(old_to_new_groups2), 
                                                          to = old_to_new_groups2)


#saveRDS(T_cells_object,"D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/T_cells_object.rds")

tcell_marker<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/marker热图分型/T_NK_marker.csv")
markers <-tcell_marker$gene_name
markers <- as.data.frame(markers)
t_vln<-VlnPlot(T_cells_object, features = unique(tcell_marker$gene_name), 
        stack = TRUE, 
        sort = TRUE, 
        cols = my36colors,
        split.by ="cd4_cd8_group" , #每种cluster 一个颜色
        flip = TRUE) +
  theme(legend.position = "none") + 
  ggtitle("Identity on x-axis")
ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/marker热图分型/marker基因小提琴图.jpg",t_vln,width = 10,height = 15)

vln.dat=FetchData(T_cells_object,c(unique(tcell_marker$gene_name),"cd4_cd8_group"))



#宽转长
vln.dat.melt <- reshape2::melt(vln.dat, id.vars = c("cd4_cd8_group"), 
                               measure.vars = unique(tcell_marker$gene_name),
                               variable.name = "gene", 
                               value.name = "Expr") 

vln.dat.melt_mean <- aggregate(vln.dat.melt$Expr, by = list(vln.dat.melt$cd4_cd8_group, vln.dat.melt$gene), FUN = mean)

colnames(vln.dat.melt_mean) <- c("Tcell_group","marker","value")
vln.dat.melt<-vln.dat.melt_mean

colnames(vln.dat.melt)<-c("Tcell_group","marker","value","fillcolor")

vln.dat.melt$marker_group <- tcell_marker[match(vln.dat.melt$marker,tcell_marker$gene_name),3]

vln.dat.melt$cell_group <- plyr::mapvalues(vln.dat.melt$Tcell_group,
                                              from = names(cell_type_list), 
                                              to = cell_type_list)




my36colors = distinctColorPalette(65)
p1 <- ggplot(vln.dat.melt, aes(gene, Expr, fill = gene)) +
  geom_violin(scale = "width", adjust = 1, trim = TRUE) +
  scale_y_continuous(expand = c(0, 0), position="right", labels = function(x)
    c(rep(x = "", times = length(x)-2), x[length(x) - 1], "")) +
  facet_grid(rows = vars(cd4_cd8_group), scales = "free", switch = "y") +
  scale_fill_manual(values = my36colors) + 
  theme_cowplot(font_size = 12) +
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
  ) +
  ggtitle("Feature on x-axis with annotation") + ylab("Expression Level")
p1

ggsave("T细胞marker小提琴图.jpg",p1,width = 15,height = 14)

aver_dt_df_long<-vln.dat.melt


aver_dt<-AverageExpression(T_cells_object,
                           features = as.character(unique(tcell_marker$gene_name)),
                           group.by =c('cd4_cd8_group'))


aver_dt_df <- as.data.frame(aver_dt)
View(aver_dt_df)
aver_dt_df$marker<- rownames(aver_dt_df)
colnames(aver_dt_df)<- c("CD3+CD4-CD8-","CD4_C1-PDCD1","CD4_C2-IL7R","CD4_C3-FOXP3", 
                         "CD4_C4-FOSB","CD8_C1-GNLY","CD8_C2-CTSW","CD8_C3-IL7R",  
                         "CD8_C4-GZMK","CD8_C5-STMN1","NK", "marker")



unique_markers <- unique(aver_dt_df$marker)
aver_dt_df$marker <- factor(aver_dt_df$marker, levels = unique_markers)

aver_dt_df_long <- pivot_longer(aver_dt_df, cols =c("CD3+CD4-CD8-","CD4_C1-PDCD1","CD4_C2-IL7R","CD4_C3-FOXP3", 
                                                    "CD4_C4-FOSB","CD8_C1-GNLY","CD8_C2-CTSW","CD8_C3-IL7R",  
                                                    "CD8_C4-GZMK","CD8_C5-STMN1","NK"), names_to = "Tcell_group", values_to = "value")
aver_dt_df_long <- aver_dt_df_long%>%
  mutate(value=scale(value))
View(aver_dt_df_long )

annotation_rows <- c("#E16181",
                     "#7E76DC",
                     "#ABDC8B",
                     "#BF45E1",
                     "#DDDC51",
                     "#DC8E50",
                     "#7BADDA",
                     "#68E3B0",
                     "#DCD69A",
                     "#D270C5",
                     "#D2DED9")

cell_type_list<- c("CD3+CD4-CD8-"="CD3+CD4-CD8-",
                   "CD4_C1-PDCD1"="CD4+",
                   "CD4_C2-IL7R"="CD4+",
                   "CD4_C3-FOXP3"="CD4+",
                   "CD4_C4-FOSB"="CD4+",
                   "CD8_C1-GNLY"="CD8+",
                   "CD8_C2-CTSW"="CD8+",
                   "CD8_C3-IL7R"="CD8+",
                   "CD8_C4-GZMK"="CD8+", 
                   "CD8_C5-STMN1"="CD8+",
                   "NK"="NK")


aver_dt_df_long$marker_group <- tcell_marker[match(aver_dt_df_long$marker,tcell_marker$gene_name),3]

aver_dt_df_long$cell_group <- plyr::mapvalues(aver_dt_df_long$Tcell_group,
                                              from = names(cell_type_list), 
                                              to = cell_type_list)

View(aver_dt_df_long)

colnames(aver_dt_df_long)
write.csv(aver_dt_df_long,"T细胞小类marker画图数据.csv")
####热图顶部色块
library(tidyverse)
library(ggnewscale)
library(MetBrewer)
library(patchwork)
library(ggtext)
marker_group_list<-c("CD4_C1-PDCD1","CD4_C2-IL7R","CD4_C3-FOXP3", 
                     "CD4_C4-FOSB","CD8_C1-GNLY","CD8_C2-CTSW","CD8_C3-IL7R",  
                     "CD8_C4-GZMK","CD8_C5-STMN1","CD3+CD4-CD8-","NK")

df_top<- aver_dt_df_long %>% filter(marker=="CD3D")
df_top$cell_group <- factor(df_top$cell_group,levels =c("CD4+","CD8+","CD3+CD4-CD8-","NK"))
df_top$Tcell_group <- factor(df_top$Tcell_group ,levels =marker_group_list)

p1_top <- df_top %>% ggplot(aes(x=Tcell_group,y=marker))+
  geom_tile(data=df_top,color="black",aes(fill=cell_group))+
  scale_fill_manual(values=c("#A6CEE3","#B2DF8A","#FDBF6F","#CAB2D6"))+
  new_scale_fill()+
  labs(x=NULL,y=NULL)+
  theme(axis.text.x = element_blank(),
        axis.text.y=element_blank(),
        legend.title = element_text(size=12, color = "black"),
        legend.position = "right",
        plot.margin=unit(c(0,1,-20,0.5),unit="cm"),
        plot.background = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank(),
        legend.text = element_text(color="black"))

p1_top

library(randomcoloR)
annotation_rows = distinctColorPalette(13)
#[1] "#D8A7D6" "#76E1C2" "#D9DD55" "#8879D6" "#D8DDD9" "#7CB8DF" "#D5DE9F" "#E09960" "#7BE36E"
#[10] "#DA5ABE" "#998E80" "#B349E7" "#E26E7F"


df_right <- aver_dt_df_long %>% group_by(marker) %>% filter (! duplicated(marker))
View(df_right)
df_right$marker_group <- factor(df_right$marker_group,levels =unique(df_right$marker_group))
df_right$marker <- factor(df_right$marker,levels =rev(unique(df_right$marker)))

p1_right <- df_right %>% ggplot(aes(x=Tcell_group,y=marker))+
  geom_tile(data=df_right,color="black",aes(fill=marker_group))+
  scale_fill_manual(values=annotation_rows)+
  new_scale_fill()+
  labs(x=NULL,y=NULL)+
  theme(axis.text.x = element_blank(),
        axis.text.y=element_blank(),
        plot.margin=unit(c(0,1,-20,0.5),unit="cm"),
        plot.background = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank(),
        legend.text = element_text(color="black"),
        legend.position = "right")
p1_right

df_data_cd4<- aver_dt_df_long%>% filter(grepl("CD4_",Tcell_group)) %>%as.data.frame()

df_data_cd4$marker <- factor(df_data_cd4$marker,levels =rev(unique(df_data_cd4$marker)))
df_data_cd4$Tcell_group <- factor(df_data_cd4$Tcell_group ,levels =marker_group_list)

p1_cd4<-df_data_cd4 %>% ggplot(aes(Tcell_group,marker)) +
  geom_tile(data=df_data_cd4,aes(fill=scale(value)), colour = "white", linewidth = 0.05)+
  scale_fill_gradientn(name='T cell marker genes mean\nscaled expression',
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
        legend.title = element_text(size=12, color = "black"),
        legend.text = element_text(size=10,color = "black",angle =45),
        legend.position = "none") +
  scale_y_discrete(position = "left")
p1_cd4


df_data_cd8<- aver_dt_df_long%>% filter(grepl("CD8_",Tcell_group)) %>%as.data.frame()
df_data_cd8$marker <- factor(df_data_cd8$marker,levels =rev(unique(df_data_cd8$marker)))
df_data_cd8$Tcell_group <- factor(df_data_cd8$Tcell_group ,levels =marker_group_list)

p1_cd8<-df_data_cd8 %>% ggplot(aes(Tcell_group,marker)) +
  geom_tile(data=df_data_cd8,aes(fill=value), colour = "white", linewidth = 0.05)+
  scale_fill_gradientn(name='T cell marker genes mean\nscaled expression',
                       #colours=c("white","red","black"))+
                       colours=c("blue","yellow","red"))+
  theme_minimal() +
  theme(axis.line=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        panel.border=element_rect(fill=NA,color="black",size=0.5,linetype="solid"),
        axis.text.x = element_text(size = 10,color = 'black',angle =40,hjust =1,vjust=1),
        legend.title = element_text(size=12, color = "black"),
        legend.text = element_text(size=10,color = "black",angle =45),
        legend.position = "none") +
  scale_y_discrete(position = "left")
p1_cd8

df_data_cd3<- aver_dt_df_long%>% filter(Tcell_group=="CD3+CD4-CD8-") %>%as.data.frame()
df_data_cd3$marker <- factor(df_data_cd3$marker,levels =rev(unique(df_data_cd3$marker)))
df_data_cd3$Tcell_group <- factor(df_data_cd3$Tcell_group ,levels =marker_group_list)


p1_cd3<-df_data_cd3 %>% ggplot(aes(Tcell_group,marker)) +
  geom_tile(data=df_data_cd3,aes(fill=value), colour = "white", linewidth = 0.05)+
  scale_fill_gradientn(name='T cell marker genes mean\nscaled expression',
                       #colours=c("white","red","black"))+
                       colours=c("blue","yellow","red"))+
  theme_minimal() +
  theme(axis.line=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y =  element_blank(),
        panel.border=element_rect(fill=NA,color="black",size=0.5,linetype="solid"),
        axis.text.x = element_text(size = 10,color = 'black',angle =40,hjust =1,vjust=1),
        legend.title = element_text(size=12, color = "black"),
        legend.text = element_text(size=10,color = "black",angle =45),
        legend.position = "none") +
  scale_y_discrete(position = "left")
p1_cd3

df_data_nk<- aver_dt_df_long%>% filter(Tcell_group=="NK") %>%as.data.frame()
df_data_nk$marker <- factor(df_data_nk$marker,levels =rev(unique(df_data_nk$marker)))
df_data_nk$Tcell_group <- factor(df_data_nk$Tcell_group ,levels =marker_group_list)


p1_nk<-df_data_nk %>% ggplot(aes(Tcell_group,marker)) +
  geom_tile(data=df_data_nk,aes(fill=value), colour = "white", linewidth = 0.05)+
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
        axis.text.y = element_blank(),
        panel.border=element_rect(fill=NA,color="black",size=0.5,linetype="solid"),
        axis.text.x = element_text(size = 10,color = 'black',angle =40,hjust =1,vjust=1),
        legend.title = element_text(size=12, color = "black"),
        legend.text = element_text(size=10,color = "black",angle =45),
        legend.position = "top") +
  scale_y_discrete(position = "left")
p1_nk

library(cowplot)
p_data_plot<-plot_grid(p1_cd4,
                       p1_cd8,
                       p1_cd3,
                       p1_nk,
                       p1_right,
                       align = 'h',
                       nrow=1,
                       rel_widths = c(5,5,1.2,1.2,5.5))
  
# plot<-plot_grid(p1_top,
#                  p_data_plot,
#                  nrow=2,
#                 align = 'h',
#                  rel_heights= c(0.2,5),
#                 rel_widths = c(0.1,5))
ggsave("tcell_细胞亚型分类3.jpg",p_data_plot,height = 11,width = 10)
ggsave("tcell_细胞亚型分类2.jpg",p_data_plot,height = 11,width = 10)
ggsave("tcell_细胞亚型分类2.pdf",p_data_plot,height = 11,width = 10)

dev.new()

####featureplot----
Idents(T_cells_object)<-"cd4_cd8_group"

tcell_cell_type<-DimPlot(T_cells_object, reduction = 'umap.harmony', group.by = 'cd4_cd8_group',
                           label = TRUE,raster=FALSE,cols=annotation_rows)

ggsave("细胞注释所有.pdf",cell_type_new_all,width =9, height =7)

####提取降维信息画featureplot----
library(tidyverse)

Tcell_microbe<-read.csv("clipboard")
View(Tcell_microbe)
microbe_all<- read.csv("filter.combined.genus.umi.matrix.csv",header = T,row.names = 1)

Tcell_microbe_top<- microbe_all[,which(colnames(microbe_all) %in% Tcell_microbe$name)]

Tcell_microbe_top$barcode <- rownames(Tcell_microbe_top)

metadata<- FetchData(T_cells_object,"cell_id")

metadata$cell_id <- rownames(metadata)

metadata<- left_join(x=metadata,y=Tcell_microbe_top,by = join_by(cell_id==barcode))

rownames(metadata)<-metadata$cell_id

T_cells_object<- AddMetaData(T_cells_object,metadata = metadata)

table(rownames(T_cells_object@meta.data) %in% Tcell_microbe_top$barcode)

Tcell_microbe$name
DimPlot(T_cells_object, reduction = 'umap.harmony', group.by = 'Massilia_group',
        label = TRUE,raster=FALSE)

x <- ifelse(T_cells_object@meta.data$Herbaspirillum>0,"Herbaspirillum+","Herbaspirillum-")
T_cells_object@meta.data$Herbaspirillum_group <- x

x <- ifelse(T_cells_object@meta.data$Massilia>0,"Massilia+","Massilia-")
T_cells_object@meta.data$Massilia_group <- x

####细胞注释umap图----
# x <- ifelse(rownames(T_cells_object@meta.data) %in% Tcell_microbe_top$barcode,"Bacteria+","Bacteria-")
# 
# T_cells_object@meta.data$group_microbe2 <- x

T_cells_object@meta.data$Massilia

Tcells.umap.harmonyredu<- T_cells_object@reductions$umap.harmony@cell.embeddings%>%
  as.data.frame() %>% 
  cbind(cell_type = T_cells_object@meta.data$cd4_cd8_group)%>% 
  cbind(Patients = T_cells_object@meta.data$Patients) %>% 
  cbind(Treatments = T_cells_object@meta.data$Treatments) %>% 
  cbind(MPR_Response2 = T_cells_object@meta.data$MPR_Response2)%>%
  cbind(Herbaspirillum = T_cells_object@meta.data$Herbaspirillum) %>%
  cbind(Massilia = T_cells_object@meta.data$Massilia) %>%
cbind(Herbaspirillum_group = T_cells_object@meta.data$Herbaspirillum_group) %>%
  cbind(Massilia_group = T_cells_object@meta.data$Massilia_group)%>%
  cbind(microbe_group = T_cells_object@meta.data$group_microbe2)

unique(Tcells.umap.harmonyredu$cell_type)
Tcells.umap.harmonyredu$cell_type <- factor(Tcells.umap.harmonyredu$cell_type,
                                            levels=c("CD4_C1-PDCD1","CD4_C2-IL7R",
                                                     "CD4_C3-FOXP3","CD4_C4-FOSB",
                                                     "CD8_C1-GNLY","CD8_C2-CTSW","CD8_C3-IL7R",
                                                     "CD8_C4-GZMK","CD8_C5-STMN1","CD3+CD4-CD8-",
                                                     "NK"))



library(randomcoloR)
Tcell_type_color = distinctColorPalette(11)


Tcell_type_color<-c("#DF66B0","#9888DB","#DC8062","#D8ACC9","#83B5D3","#D5D8C5", "#77DCCF", "#BA4EE4",
                    "#AEE94F","#89DB8E","#E0D175")

Tcells.umap.harmonyredu$microbe_group <- factor(Tcells.umap.harmonyredu$microbe_group,levels=c("Bacteria+","Bacteria-"))

cell_type_plot <- ggplot(Tcells.umap.harmonyredu,aes(x= umapharmony_1, y = umapharmony_2,color = cell_type)) +  
  geom_point(size = 0.7, alpha =1)  +  
  scale_color_manual(values = Tcell_type_color)+
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
        legend.text = element_text(size=12,face = "bold"), #设置legend标签的大小
        legend.key.size=unit(1,'cm') ) +  # 设置legend标签之间的大小
  guides(color = guide_legend(override.aes = list(size=5)))+ #设置legend中 点的大小 
  geom_segment(aes(x = min(Tcells.umap.harmonyredu$umapharmony_1) , y = min(Tcells.umap.harmonyredu$umapharmony_2),
                   xend = min(Tcells.umap.harmonyredu$umapharmony_1) +4, yend = min(Tcells.umap.harmonyredu$umapharmony_2) ),
               colour = "black", size=1,arrow = arrow(length = unit(0.5,"cm")))+ 
  geom_segment(aes(x = min(Tcells.umap.harmonyredu$umapharmony_1)  , y = min(Tcells.umap.harmonyredu$umapharmony_2)  ,
                   xend = min(Tcells.umap.harmonyredu$umapharmony_1) , yend = min(Tcells.umap.harmonyredu$umapharmony_2) + 4),
               colour = "black", size=1,arrow = arrow(length = unit(0.5,"cm"))) +
  annotate("text", x = min(Tcells.umap.harmonyredu$umapharmony_1) +2, y = min(Tcells.umap.harmonyredu$umapharmony_2) -1, label = "UMAP_1",
           color="black",size = 4, fontface="bold") + 
  annotate("text", x = min(Tcells.umap.harmonyredu$umapharmony_1)-1, y = min(Tcells.umap.harmonyredu$umapharmony_2) + 2, label = "UMAP_2",
           color="black",size = 4, fontface="bold" ,angle=90)

cell_type_plot

cell_type_med <- Tcells.umap.harmonyredu %>%
  group_by(cell_type) %>%
  summarise(
    umapharmony_1 = median(umapharmony_1),
    umapharmony_2 = median(umapharmony_2)
  )

library(ggrepel)
cell_type_med

library(scales)
cell_type_plot2<-cell_type_plot+geom_label_repel(aes(label=cell_type),size=4,color="black",fontface="bold",data = cell_type_med,
                                                        point.padding = NA,label.size = NA, fill = alpha("white",0.2),
                                                        segment.size=0.5,force = 1,nudge_x=0.5, nudge_y = 0,direction="y",max.overlaps=50)+
  theme(legend.position = "right")


ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/featureplot/T cell 亚群注释umap图.jpg",cell_type_plot2,width = 8,height=6)
ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/featureplot/T cell 亚群注释umap图.pdf",cell_type_plot2,width = 8,height=6)
T_cells_object@assays$RNA@counts['CD4',]

######featureplot----
library(Seurat)
library(ggplot2)
library(cowplot)
library(dplyr)
library(tidydr)
library(stringr)
library(viridis)
devtools::install_github("r-lib/rlang")
library(scCustomize)

gene<-c("CD3D","CD4","CD8A","GNLY","CTSW","PDCD1",
        "IL7R","GZMK","STMN1","FOXP3","FOSB")

pal <- viridis(n = 10, option = "C")
# "#0D0887FF" "#47039FFF" "#7301A8FF" "#9C179EFF" "#BD3786FF" "#D8576BFF" "#ED7953FF" "#FA9E3BFF" "#FDC926FF" "#F0F921FF"
pal <- viridis(n = 15, option = "D", direction = -1)

p<- FeaturePlot(object =T_cells_object,features = gene,reduction = "umap.harmony",raster=FALSE, cols = pal, order = T)

ggsave(file="T cell_marker_featureplot1.jpg",p,width=15,height=9)
ggsave(file="T cell_marker_featureplot1.pdf",p,width=15,height=9)

setwd("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/featureplot")

gene<-c("CD3E","CD4","CD8A","GNLY","CTSW","PDCD1",
        "IL7R","GZMK","STMN1","FOXP3","FOSB")
i=1
plots=list()
for (i in 1:length(gene)){
  plots[[i]]=FeaturePlot(object=T_cells_object,features = gene[i],reduction = "umap.harmony",raster=FALSE, cols = pal, order = T)+
    xlab("UMAP_1")+
    ylab("UMAP_1")+
    theme(legend.position = "none")
}

library(patchwork)
p<-wrap_plots(plots, ncol = 4)+plot_annotation(tag_levels = "A");p

p_legend<-FeaturePlot(object=T_cells_object,features = "FOSB",reduction = "umap.harmony",raster=FALSE, cols = pal, order = T)+
  xlab("UMAP_1")+
  ylab("UMAP_1")+
  theme(legend.position = "right")

legend <- get_legend(p_legend+
                       #guides(color = guide_legend(nrow = 1)) +
                       theme(legend.position = "right"))
plt_featureplot2<- plot_grid(p,legend, ncol=2,rel_widths = c(1,0.1))

ggsave(file="./Tcell_marker_featureplot.jpg",plt_featureplot2,width=15,height=9)
ggsave(file="./Tcell_marker_featureplot.pdf",plt_featureplot2,width=15,height=9)

####细胞比例计算----
setwd("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/细胞比例")

Idents(T_cells_object)<-"cd4_cd8_group"
table(T_cells_object@meta.data$MPR_Response2)
prop.table(table(Idents(T_cells_object)))
table(Idents(T_cells_object),T_cells_object@meta.data$MPR_Response2)
Cellratio <- prop.table(table(Idents(T_cells_object),T_cells_object@meta.data$MPR_Response2), margin = 2)
Cellratio <- as.data.frame(Cellratio)
colourCount = length(unique(Cellratio$Var1))
library(ggplot2)
ggplot(Cellratio) + 
  geom_bar(aes(x =Var2, y= Freq, fill = Var1),stat = "identity",width = 0.7,size = 0.5,colour = '#222222')+ 
  theme_classic() +
  labs(x='Sample',y = 'Ratio')+
  coord_flip()+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))

table(T_cells_object$orig.ident)#查看各组细胞数
prop.table(table(Idents(T_cells_object)))
table(Idents(T_cells_object), T_cells_object$orig.ident)#各组不同细胞群细胞数
Cellratio <- prop.table(table(Idents(T_cells_object), T_cells_object$orig.ident), margin = 2)#计算各组样本不同细胞群比例
Cellratio <- data.frame(Cellratio)

cellper <- dcast(Cellratio,Var2~Var1, value.var = "Freq")#长数据转为宽数据

rownames(cellper) <- cellper[,1]
cellper <- cellper[,-1]
rownames(cellper)
###添加分组信息
microbe_group<-cbind(T_cells_object@meta.data$orig.ident)%>%
  cbind(cell_type = T_cells_object@meta.data$cd4_cd8_group)%>% 
  cbind(Patients = T_cells_object@meta.data$Patients) %>% 
  cbind(Treatments = T_cells_object@meta.data$Treatments) %>% 
  cbind(MPR_Response2 = T_cells_object@meta.data$MPR_Response2)%>%
  cbind(Herbaspirillum = T_cells_object@meta.data$Herbaspirillum) %>%
  cbind(Massilia = T_cells_object@meta.data$Massilia) %>%
  cbind(Herbaspirillum_group = T_cells_object@meta.data$Herbaspirillum_group) %>%
  cbind(Massilia_group = T_cells_object@meta.data$Massilia_group)%>%
  cbind(microbe_group = T_cells_object@meta.data$group_microbe2) %>% as.data.frame()

colnames(microbe_group) <- c("samples","cell_type","Patients","Treatments",          
                            "MPR_Response2","Herbaspirillum","Massilia","Herbaspirillum_group",
                            "Massilia_group","microbe_group")
View(microbe_group)
microbe_group$samples
colnames(microbe_group)
colnames(cellper)
cellper$samples<- microbe_group[match(rownames(cellper),microbe_group$samples),1]
cellper$MPR_Response2<- microbe_group[match(rownames(cellper),microbe_group$samples),5]

write.csv(cellper,"T cell_细胞比例cellper.csv")
cellper<- read.csv("T cell_细胞比例cellper.csv",row.names = 1,header = T)
cell_type =c("CD4_C1-PDCD1","CD4_C2-IL7R","CD4_C3-FOXP3","CD4_C4-FOSB",
"CD8_C1-GNLY","CD8_C2-CTSW","CD8_C3-IL7R","CD8_C4-GZMK","CD8_C5-STMN1","CD3+CD4-CD8-","NK")

library(ggplot2)
library(dplyr)
library(ggpubr)
library(tidyverse)
library(ggh4x)
library(rstatix)
library(ggpubr)

####长数据添加分组信息
Cellratio$samples<- microbe_group[match(Cellratio$Var2,microbe_group$samples),1]

Cellratio$MPR_Response2<- microbe_group[match(Cellratio$Var2,microbe_group$samples),5]
colnames(Cellratio)<-c("cell_type","samples_id","Freq","samples","MPR_response2")

Cellratio$MPR_response2<- factor(Cellratio$MPR_response2,levels=c("Pre_NMPR","Pre_MPR","Post_MPR"))
Cellratio$cell_type<- factor(Cellratio$cell_type,levels=c("CD4_C1-PDCD1","CD4_C2-IL7R","CD4_C3-FOXP3", 
                                                            "CD4_C4-FOSB","CD8_C1-GNLY","CD8_C2-CTSW","CD8_C3-IL7R",  
                                                            "CD8_C4-GZMK","CD8_C5-STMN1","CD3+CD4-CD8-","NK"))

df_p_val1 <- Cellratio %>% group_by(cell_type)%>%
  wilcox_test(Freq~MPR_response2) %>%
  adjust_pvalue(p.col="p",method="bonferroni") %>%
  add_significance(p.col="p.adj")
  # add_xy_position(x="MPR_Response2",dodge=0.1)
Tcell_type_color<-c("#77DCCB","#D4D2CE","#D9ADD9","#D8DC52","#D16BC6","#D86A70","#DDA45C"
                    ,"#A851E2" ,"#77E367","#7A9BD2","#C7DC9A")
p_Cellratio <- Cellratio %>% ggplot(aes(MPR_response2,Freq))+
  # annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 1, 
  #          fill = "cornflowerblue", alpha = .3, color = NA)+
  # annotate("rect", xmin = -Inf, xmax = Inf, ymin =-1, ymax = Inf, 
  #          fill = "#FCAE12", alpha = .3, color = NA) +
  geom_violin(aes(fill=MPR_response2),trim = FALSE,show.legend = F)+
  geom_boxplot(width = 0.2,outliers = FALSE, staplewidth = 0.5) +
  scale_y_continuous(sec.axis = sec_axis(~ ., name = ""))+
  stat_summary(fun=mean,geom="point",col="black",fill="#F98400",
               shape=23,show.legend = F)+
  stat_summary(fun=mean, geom="line", aes(group=MPR_response2), col="black")+
  # stat_pvalue_manual(df_p_val1,label="p.adj.signif",hide.ns=F,y.position=0.1,
  #                    tip.length = 0,label.size = 5,color="black")+
  # geom_hline(yintercept = 1,linetype=2)+
  # geom_hline(yintercept = -1,linetype=2)+
  # facet_wrap(.~cell_type,scales="free_y")+
  facet_nested_wrap(. ~cell_type,scales="free_y",
                    strip = strip_nested(background_x =
                                           elem_list_rect(fill=Tcell_type_color))) +
  scale_fill_manual(values = c("#3B9AB2","#7294D4","#E6A0C4"))+
  labs(x=NULL,y=NULL)+
  theme(axis.text.x=element_text(angle = 0,vjust=0.5,hjust=0.5,color="black"),
        axis.text.y=element_text(color="black"),
        plot.background = element_rect(fill="white"), 
        panel.background = element_rect(fill="white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing.x = unit(0,"cm"),
        plot.margin=unit(c(0.5,0.5,0.5,0.5),unit="cm"),
        axis.line.x=element_line(color="black"),
        axis.line.y.left = element_line(color="grey30"),
        axis.line.y.right = element_line(color="grey30"),
        axis.line.x.bottom = element_line(color="grey30"),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank())+
  guides(y = guide_axis(minor.ticks = TRUE))

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/细胞比例/t细胞比例.jpg",p_Cellratio,width = 15,height = 10)
ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/细胞比例/t细胞比例.pdf",p_Cellratio,width = 15,height = 10)


######OR指数比较单细胞亚群的组织偏好----
BiocManager::install("impute")
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
out.prefix <- "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/细胞比例/T_cell_OR"

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

meta.tb <- T_cells_object@meta.data

names(meta.tb)
View(meta.tb)

out.prefix <- "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/细胞比例/T cell_RO"
A <- do.tissueDist(cellInfo.tb = meta.tb,
                   meta.cluster = meta.tb$cd4_cd8_group,
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

library(RColorBrewer)
#自己做图
data <- A$count.dist.melt.ext.tb

data$cid <- factor(data$cid,levels=c("Pre_NMPR","Pre_MPR","Post_MPR"))
data$rid <- factor(data$rid, levels=rev(c("CD4_C1-PDCD1","CD4_C2-IL7R","CD4_C3-FOXP3", 
                                      "CD4_C4-FOSB","CD8_C1-GNLY","CD8_C2-CTSW","CD8_C3-IL7R",  
                                      "CD8_C4-GZMK","CD8_C5-STMN1","CD3+CD4-CD8-","NK")))

bl <- colorRampPalette(c("navy","royalblue","lightskyblue"))(100)                      
re <- colorRampPalette(c("mistyrose", "red2","darkred"))(100)


OR_plot<-ggplot(data, aes(cid, rid)) + 
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

ggsave("OR_plot_MPR.jpg",OR_plot,width=5,height=7)
ggsave("OR_plot_MPR.pdf",OR_plot,width=5,height=7)

ggsave("OR_plot_MPR_treatments.jpg",OR_plot,width=6,height=7)

ggsave("OR_plot_MPR_treatments.pdf",OR_plot,width=6,height=7)




######CD3细胞比例----

Idents(T_cells_object)

levels(Idents(T_cells_object))

c("CD8_C4-GZMK","CD8_C2-CTSW","CD4_C3-FOXP3","NK","CD3+CD4-CD8-","CD4_C2-IL7R","CD8_C3-IL7R", 
"CD8_C1-GNLY","CD8_C5-STMN1","CD4_C4-FOSB","CD4_C1-PDCD1")


T_cells_CD3_object<- subset(T_cells_object,idents = c("CD8_C4-GZMK","CD8_C2-CTSW","CD4_C3-FOXP3",
                                                      "CD3+CD4-CD8-","CD4_C2-IL7R","CD8_C3-IL7R", 
                                                      "CD8_C1-GNLY","CD8_C5-STMN1","CD4_C4-FOSB","CD4_C1-PDCD1"))
  
  
levels(Idents(T_cells_CD3_object))

Cellratio1 <- prop.table(table(Idents(T_cells_CD3_object),T_cells_CD3_object@meta.data$orig.ident), margin = 2)
Cellratio1 <- as.data.frame(Cellratio1)
colnames(Cellratio1) <- c("cell_type","samples_id","Freq")

microbe_group<-cbind(T_cells_object@meta.data$orig.ident)%>%
  cbind(cell_type = T_cells_object@meta.data$cd4_cd8_group)%>% 
  cbind(Patients = T_cells_object@meta.data$Patients) %>% 
  cbind(Treatments = T_cells_object@meta.data$Treatments) %>% 
  cbind(MPR_Response2 = T_cells_object@meta.data$MPR_Response2)%>%
  cbind(Herbaspirillum = T_cells_object@meta.data$Herbaspirillum) %>%
  cbind(Massilia = T_cells_object@meta.data$Massilia) %>%
  cbind(Herbaspirillum_group = T_cells_object@meta.data$Herbaspirillum_group) %>%
  cbind(Massilia_group = T_cells_object@meta.data$Massilia_group)%>%
  cbind(microbe_group = T_cells_object@meta.data$group_microbe2) %>% as.data.frame()

colnames(microbe_group) <- c("samples","cell_type","Patients","Treatments",          
                             "MPR_Response2","Herbaspirillum","Massilia","Herbaspirillum_group",
                             "Massilia_group","microbe_group")

Cellratio1$samples<- microbe_group[match(Cellratio1$samples,microbe_group$samples),1]


Cellratio1$MPR_Response2<- microbe_group[match(Cellratio1$samples,microbe_group$samples),5]
colnames(Cellratio1)<-c("cell_type","samples_id","Freq","samples","MPR_response2")
Cellratio1$group<- c(rep("CD3+ T cells",140))

Cellratio1$MPR_response2<- factor(Cellratio1$MPR_response2,levels=c("Pre_NMPR","Pre_MPR","Post_MPR"))
Cellratio1$cell_type<- factor(Cellratio1$cell_type,levels=c("CD4_C1-PDCD1","CD4_C2-IL7R","CD4_C3-FOXP3", 
                                                          "CD4_C4-FOSB","CD8_C1-GNLY","CD8_C2-CTSW","CD8_C3-IL7R",  
                                                          "CD8_C4-GZMK","CD8_C5-STMN1","CD3+CD4-CD8-","NK"))



df_p_val1 <- Cellratio1 %>% group_by(group)%>%
  wilcox_test(Freq~MPR_response2) %>%
  adjust_pvalue(p.col="p",method="bonferroni") %>%
  add_significance(p.col="p.adj")%>% 
  add_xy_position(x="MPR_response2")

library(randomcoloR)
color<- distinctColorPalette(20)

df_p_val1 <- Cellratio %>% group_by(group)%>%
  wilcox_test(Freq~MPR_response2) %>%
  adjust_pvalue(p.col="p",method="bonferroni") %>%
  add_significance(p.col="p.adj")%>% 
  add_xy_position(x="MPR_response2",dodge=2)

p_Cellratio <- Cellratio %>% ggplot(aes(MPR_response2,Freq))+
  geom_violin(aes(fill=MPR_response2),trim = T,show.legend = F)+
  geom_boxplot(width = 0.2,outliers = FALSE, staplewidth = 0.5) +
  scale_y_continuous(sec.axis = sec_axis(~ ., name = ""))+
  stat_summary(fun=mean,geom="point",col="black",fill="#F98400",
               shape=23,show.legend = F)+
  #stat_summary(fun=mean, geom="line", aes(group=MPR_response2), col="black")+
  stat_pvalue_manual(df_p_val1,label="p.adj.signif",hide.ns=F,
                     tip.length = 0,label.size = 5,face="bold",color="black")+
  scale_fill_manual(values = c("#3B9AB2","#7294D4","#E6A0C4"))+
  #scale_y_continuous(limits = c(-0.5, 1.5))+
  facet_grid(~group,scales = 'free',space = "free")+
  labs(x=NULL,y="Cell percent ratio")+
  theme(strip.background = element_rect(fill ="#CCE5DA", color = "black",size = 1),  # 背景颜色和边框颜色
        strip.text.x = element_text(color = "black", size = 13, face = "bold"),  
        axis.text.x=element_text(size=13,angle = 0,vjust=0.5,hjust=0.5,color="black"),
        axis.text.y=element_text(size=15,color="black"),
        axis.title.y=element_text(size=15,color="black"),
        panel.spacing.x = unit(0,"cm"),
        panel.background = element_rect(fill = "white", colour = "black", size = 1),  # Set panel background to white
        panel.border = element_rect(colour = "black", fill = NA, size = 1),  # Set panel border
        plot.background = element_rect(fill = "white", colour = NA),  # Optional: set plot background to white
        axis.line = element_line(colour = "black", size = 1),
        plot.margin=unit(c(0.5,0.5,0.5,0.5),unit="cm"),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank())+
  guides(y = guide_axis(minor.ticks = TRUE))

p_Cellratio

######CD3+CD4+CD8-细胞比例----

CD3CD4_object<- subset(T_cells_object,idents = c("CD4_C3-FOXP3","CD4_C2-IL7R","CD4_C4-FOSB","CD4_C1-PDCD1"))

levels(Idents(CD3CD4_object))

Cellratio2 <- prop.table(table(Idents(CD3CD4_object),CD3CD4_object@meta.data$orig.ident), margin = 2)
Cellratio2 <- as.data.frame(Cellratio2)
colnames(Cellratio2) <- c("cell_type","samples_id","Freq")

Cellratio2$samples<- microbe_group[match(Cellratio2$samples,microbe_group$samples),1]
Cellratio2$MPR_Response2<- microbe_group[match(Cellratio2$samples,microbe_group$samples),5]
colnames(Cellratio2)<-c("cell_type","samples_id","Freq","samples","MPR_response2")
Cellratio2$group<- c(rep("CD3+CD4+CD8- T cells",56))

Cellratio2$MPR_response2<- factor(Cellratio2$MPR_response2,levels=c("Pre_NMPR","Pre_MPR","Post_MPR"))
Cellratio2$cell_type<- factor(Cellratio2$cell_type,levels=c("CD4_C1-PDCD1","CD4_C2-IL7R","CD4_C3-FOXP3", 
                                                            "CD4_C4-FOSB"))

Cellratio2 <- Cellratio2 %>%
  filter(Freq != 0)

Cellratio2$MPR_response2<- factor(Cellratio2$MPR_response2,levels=c("Pre_NMPR","Pre_MPR","Post_MPR"))
Cellratio2$cell_type<- factor(Cellratio2$cell_type,levels=c("CD4_C1-PDCD1","CD4_C2-IL7R","CD4_C3-FOXP3", 
                                                            "CD4_C4-FOSB"))



df_p_val2 <- Cellratio2 %>% group_by(group)%>%
  wilcox_test(Freq~MPR_response2) %>%
  adjust_pvalue(p.col="p",method="bonferroni") %>%
  add_significance(p.col="p.adj")%>% 
  add_xy_position(x="MPR_response2")

library(randomcoloR)
color<- distinctColorPalette(3)

df_p_val2 <- Cellratio2 %>% group_by(group)%>%
  wilcox_test(Freq~MPR_response2) %>%
  adjust_pvalue(p.col="p",method="bonferroni") %>%
  add_significance(p.col="p.adj")%>% 
  add_xy_position(x="MPR_response2",dodge=2)

p_Cellratio2 <- Cellratio2 %>% ggplot(aes(MPR_response2,Freq))+
  geom_violin(aes(fill=MPR_response2),trim = T,show.legend = F)+
  geom_boxplot(width = 0.2,outliers = FALSE, staplewidth = 0.5) +
  scale_y_continuous(sec.axis = sec_axis(~ ., name = ""))+
  stat_summary(fun=mean,geom="point",col="black",fill="#F98400",
               shape=23,show.legend = F)+
  #stat_summary(fun=mean, geom="line", aes(group=MPR_response2), col="black")+
  stat_pvalue_manual(df_p_val2,label="p.adj.signif",hide.ns=F,
                     tip.length = 0,label.size = 5,face="bold",color="black")+
  scale_fill_manual(values = c("#3B9AB2","#7294D4","#E6A0C4"))+
  #scale_y_continuous(limits = c(-0.5, 1.5))+
  facet_grid(~group,scales = 'free',space = "free")+
  labs(x=NULL,y=NULL)+
  theme(strip.background = element_rect(fill = "#B3DB7C", color = "black",size = 1),  # 背景颜色和边框颜色
        strip.text.x = element_text(color = "black", size = 13, face = "bold"),  
        axis.text.x=element_text(size=13,angle = 0,vjust=0.5,hjust=0.5,color="black"),
        axis.text.y=element_text(size=15,color="black"),
        panel.spacing.x = unit(0,"cm"),
        panel.background = element_rect(fill = "white", colour = "black", size = 1),  # Set panel background to white
        panel.border = element_rect(colour = "black", fill = NA, size = 1),  # Set panel border
        plot.background = element_rect(fill = "white", colour = NA),  # Optional: set plot background to white
        axis.line = element_line(colour = "black", size = 1),
        plot.margin=unit(c(0.5,0.5,0.5,0.5),unit="cm"),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank())+
  guides(y = guide_axis(minor.ticks = TRUE))

p_Cellratio2

####CD3+CD8+细胞比例----

CD3CD8_object<- subset(T_cells_object,idents = c("CD8_C1-GNLY","CD8_C2-CTSW","CD8_C3-IL7R",  
                                                 "CD8_C4-GZMK","CD8_C5-STMN1"))

levels(Idents(CD3CD8_object))

Cellratio3 <- prop.table(table(Idents(CD3CD8_object),CD3CD8_object@meta.data$orig.ident), margin = 2)
Cellratio3 <- as.data.frame(Cellratio3)
colnames(Cellratio3) <- c("cell_type","samples_id","Freq")

Cellratio3$samples<- microbe_group[match(Cellratio3$samples,microbe_group$samples),1]


Cellratio3$MPR_Response2<- microbe_group[match(Cellratio3$samples,microbe_group$samples),5]

colnames(Cellratio3)<-c("cell_type","samples_id","Freq","samples","MPR_response2")

Cellratio3$group<- c(rep("CD3+CD4-CD8+ T cells",70))

Cellratio3$MPR_response2<- factor(Cellratio3$MPR_response2,levels=c("Pre_NMPR","Pre_MPR","Post_MPR"))
Cellratio3$cell_type<- factor(Cellratio3$cell_type,levels=c("CD8_C1-GNLY","CD8_C2-CTSW","CD8_C3-IL7R",  
                                                            "CD8_C4-GZMK","CD8_C5-STMN1"))

Cellratio3 <- Cellratio3 %>%
  filter(Freq != 0)

Cellratio3$MPR_response2<- factor(Cellratio3$MPR_response2,levels=c("Pre_NMPR","Pre_MPR","Post_MPR"))
Cellratio3$cell_type<- factor(Cellratio3$cell_type,levels=c("CD8_C1-GNLY","CD8_C2-CTSW","CD8_C3-IL7R",  
                                                            "CD8_C4-GZMK","CD8_C5-STMN1"))

library(randomcoloR)
color<- distinctColorPalette(3)


df_p_val3 <- Cellratio3 %>% group_by(group)%>%
  wilcox_test(Freq~MPR_response2) %>%
  adjust_pvalue(p.col="p",method="bonferroni") %>%
  add_significance(p.col="p.adj")%>% 
  add_xy_position(x="MPR_response2",dodge=2)

p_Cellratio3 <- Cellratio3 %>% ggplot(aes(MPR_response2,Freq))+
  geom_violin(aes(fill=MPR_response2),trim =T,show.legend = F)+
  geom_boxplot(width = 0.2,outliers = FALSE, staplewidth = 0.5) +
  scale_y_continuous(sec.axis = sec_axis(~ ., name = ""))+
  stat_summary(fun=mean,geom="point",col="black",fill="#F98400",
               shape=23,show.legend = F)+
  #stat_summary(fun=mean, geom="line", aes(group=MPR_response2), col="black")+
  stat_pvalue_manual(df_p_val3,label="p.adj.signif",hide.ns=F,
                     tip.length = 0,label.size = 5,face="bold",color="black")+
  scale_fill_manual(values = c("#3B9AB2","#7294D4","#E6A0C4"))+
  #scale_y_continuous(limits = c(-0.5, 1.5))+
  facet_grid(~group,scales = 'free',space = "free")+
  labs(x=NULL,y=NULL)+
  theme(strip.background = element_rect(fill = "lightblue", color = "black",size = 1),  # 背景颜色和边框颜色
         strip.text.x = element_text(color = "black", size = 13, face = "bold"),  
        axis.text.x=element_text(size=13,angle = 0,vjust=0.5,hjust=0.5,color="black"),
        axis.text.y=element_text(size=12,color="black"),
        panel.spacing.x = unit(0,"cm"),
        panel.background = element_rect(fill = "white", colour = "black", size = 1),  # Set panel background to white
        panel.border = element_rect(colour = "black", fill = NA, size = 1),  # Set panel border
        plot.background = element_rect(fill = "white", colour = NA),  # Optional: set plot background to white
        axis.line = element_line(colour = "black", size = 1),
        plot.margin=unit(c(0.5,0.5,0.5,0.5),unit="cm"),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank())+
  guides(y = guide_axis(minor.ticks = TRUE))

p_Cellratio3


p_Cellratio+p_Cellratio2+p_Cellratio3

#组合
library(cowplot)
p<-cowplot::plot_grid(p_Cellratio,p_Cellratio2,p_Cellratio3,nrow=1, rel_widths = c(1, 1,1), labels=LETTERS[1:3])

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/细胞比例/CD3cd4cd8小提琴图.jpg",p,width = 15,height=6)

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/细胞比例/CD3cd4cd8小提琴图.pdf",p,width = 15,height=6)





####AddModuleScore----
#通过msigdb官网找到自己想要的通路集合后直接通过msigdbr包下载数据集，
#使用AddModuleScore对细胞进行通路打分，这里还可以使用其他的打分方法，
#比如之前介绍的AUCell（AUCell：对单个细胞的基因集活性打分 (qq.com)）等。

library(msigdbr)
library(RColorBrewer)
gene_sets <- msigdbr(species = "Homo sapiens", category = "C2") %>%
  dplyr::filter(gs_name=="REACTOME_TCR_SIGNALING")
TCR_l <-list(TCR_signaling=unique(gene_sets$gene_symbol))
sce_T <- AddModuleScore(sce_T,TCR_l)
FeaturePlot(sce_T,'Cluster1',cols=rev(brewer.pal(10, name = "RdBu")))
ggsave("TCR signaling.pdf")
#cell cycle
cell_cycling <- msigdbr(species = "Homo sapiens", category = "C5") %>%
  dplyr::filter(gs_name=="GOBP_CELL_CYCLE")
cycling_l <- list(cell_cycling=unique(cell_cycling$gene_symbol))
sce_T <- AddModuleScore(sce_T,cycling_l)
FeaturePlot(sce_T,'Cluster1',cols = rev(brewer.pal(10, name = "RdBu")))
ggsave("cell cycling.pdf")
saveRDS(sce_T,"sce_T.rds")

####把细分的亚群放回原来的总群----

pbmc=readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/T_cells_object.rds")
dim(bcell)
dim(pbmc)
Idents(pbmc, cells = colnames(bcell)) <- Idents(bcell)

DimPlot(pbmc,label = TRUE)

#### Tcell微生物----

library(tidyverse)
umap.harmonyredu<- T_cells_object@reductions$umap.harmony@cell.embeddings%>%
  as.data.frame() %>%
  cbind(microbe_detected2 = T_cells_object@meta.data$group_microbe2)

names(umap.harmonyredu)

write.csv(umap.harmonyredu,"T_cell_微生物降维数据.csv")

###颜色设置###画小坐标轴####https://cloud.tencent.com/developer/article/1924260

umap.harmonyredu$microbe_detected2 <- factor(umap.harmonyredu$microbe_detected2,levels=c("Bacteria+","Bacteria-"))

bacteria_plot <- ggplot(umap.harmonyredu,aes(x= umapharmony_1, y = umapharmony_2,color = microbe_detected2)) +
  geom_point(size = 0.7, alpha =0.5)  +
  scale_color_manual(values = c("#ff7f24","grey97"))+
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
        legend.text = element_text(size=12,face = "bold"), #设置legend标签的大小
        legend.key.size=unit(1,'cm') ) +  # 设置legend标签之间的大小
  guides(color = guide_legend(override.aes = list(size=5)))+ #设置legend中 点的大小
  geom_segment(aes(x = min(umap.harmonyredu$umapharmony_1) , y = min(umap.harmonyredu$umapharmony_2) ,
                   xend = min(umap.harmonyredu$umapharmony_1) +5, yend = min(umap.harmonyredu$umapharmony_2) ),
               colour = "black", size=1,arrow = arrow(length = unit(0.5,"cm")))+
  geom_segment(aes(x = min(umap.harmonyredu$umapharmony_1)  , y = min(umap.harmonyredu$umapharmony_2)  ,
                   xend = min(umap.harmonyredu$umapharmony_1) , yend = min(umap.harmonyredu$umapharmony_2) + 5),
               colour = "black", size=1,arrow = arrow(length = unit(0.5,"cm"))) +
  annotate("text", x = min(umap.harmonyredu$umapharmony_1) +1.5, y = min(umap.harmonyredu$umapharmony_2) -1, label = "UMAP_1",
           color="black",size = 4, fontface="bold" ) +
  annotate("text", x = min(umap.harmonyredu$umapharmony_1) -1, y = min(umap.harmonyredu$umapharmony_2) + 1.5, label = "UMAP_2",
           color="black",size = 4, fontface="bold" ,angle=90)


bacteria_plot


cell_type_med <- umap.harmonyredu %>%
  group_by(cell_type) %>%
  summarise(
    umapharmony_1 = median(umapharmony_1),
    umapharmony_2 = median(umapharmony_2)
  )

library(ggrepel)
cell_type_med

cell_type2<-cell_type +geom_label_repel(aes(label=cell_type),size=4,color="black",fontface="bold",data = MPRNMPR_cell_type_med,
                                        point.padding = NA,label.size = NA, fill = alpha(c("white"),0.2),
                                        segment.size=0.5,force = 1,nudge_x=0.5, nudge_y = 0,direction="y",max.overlaps=50)+
  theme(legend.position = "none")

cell_type2<-MPRNMPR_cell_type +geom_label_repel(aes(label=cell_type),size=4,color="black",fontface="bold",data = MPRNMPR_cell_type_med,
                                                point.padding = NA,label.size = NA, fill = alpha(c("white"),0.2),
                                                segment.size=0.5,force = 1,nudge_x=0.5, nudge_y = 0,direction="y",max.overlaps=50)+
  theme(legend.position = "top")


####T cell 微生物差异基因表达----



####T cell PDCD1有无微生物----

View(T_cells_object)
markers <- c("PDCD1")
markers <- as.data.frame(markers)

markerdata <- ScaleData(T_cells_object, features = as.character(unique(markers$markers)), assay = "RNA")

Idents(T_cells_object) <- "cell_id"

Idents(T_cells_object)

aver_dt<-AverageExpression(T_cells_object,
                           features = "PDCD1",
                           group.by = "cell_id")


aver_dt <- as.data.frame(aver_dt$RNA)
aver_dt<-as.data.frame(t(aver_dt))
colnames(aver_dt)<- "Average_Expression"

group<-cbind(T_cells_object@meta.data$cell_id)%>%
  cbind(T_cells_object@meta.data$orig.ident) %>%
  cbind(T_cells_object@meta.data$group_microbe2) %>%
  cbind(T_cells_object@meta.data$MPR_Response2) %>%
  as.data.frame()
colnames(group) <- c("cell_id","orig.ident","group_microbe2","MPR_Response2")

rownames(aver_dt)<- gsub('-', '_', rownames(aver_dt), fixed = TRUE)

aver_dt$microbe_group <- group[match(rownames(aver_dt),group$cell_id),3]
aver_dt$MPR_Response2 <- group[match(rownames(aver_dt),group$cell_id),4]
colnames(aver_dt)


setwd("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/PDCD1有无微生物表达对比")

library(tidyverse)
#devtools::install_github("psyteachr/introdataviz")
library(introdataviz)
library(ggpubr)
library(scales)
library(patchwork)

sessionInfo()
df<-aver_dt

df <- df %>%
  filter(Average_Expression != 0)

df$MPR_Response2 <- factor(df$MPR_Response2,levels = c("Pre_NMPR","Pre_MPR","Post_MPR"))


a <- ggplot(df,aes(x = MPR_Response2, y =Average_Expression,fill=microbe_group))+
  geom_split_violin(trim =F,color = NA,adjust = 1.5)+
  guides(fill=guide_legend(title="group"))+
  stat_summary(fun.data = "mean_sd",position=position_dodge(0.15),geom = "errorbar",width = .1) +
  stat_summary(fun = "mean", geom = "point", position=position_dodge(0.15),show.legend = F)+
  stat_compare_means(aes(group = microbe_group), label = "p.signif",label.y=22, method="wilcox.test",size=7)+
  scale_fill_manual(values=c("#788FCE","#E6956F"))+
  scale_y_continuous(limits = c(0, 25))+
  labs(x=NULL,y="PDCD1 \n average expression")+
  theme(
    axis.text.x=element_text(angle =0,hjust =0.5,vjust =0.5,color="black",size = 16,margin = margin(b =2)),
    axis.text.y=element_text(color="black",size = 15,margin = margin(r =1)),
    axis.title.y=element_text(size = 17),
    panel.background = element_rect(fill = NA,color = NA),
    panel.grid.minor= element_line(size=0.2,color="#e5e5e5"),
    panel.grid.major = element_line(size=0.2,color="#e5e5e5"),
    panel.border = element_rect(fill=NA,color="black",size=1,linetype="solid"),
    legend.key=element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(color="black",size=13),
    legend.spacing.x=unit(0.1,'cm'),
    legend.key.width=unit(0.6,'cm'),
    legend.key.height=unit(0.6,'cm'),
    legend.position = "top",
    legend.background=element_blank())

a

ggsave("PDCD1_response_bacteria.jpg",a,width=7,height=7)
ggsave("PDCD1_response_bacteria.jpg.pdf",a,width=5,height=3)


######
library(viridis)
library(scCustomize)
library(RColorBrewer)

aver_dt<-AverageExpression(markerdata,
                           features = "PDCD1",
                           group.by =c('group_microbe2','MPR_Response2'))


aver_dt <- as.data.frame(aver_dt$RNA)

aver_dt<-as.data.frame(t(aver_dt))
colnames(aver_dt)<- "Average_Expression"
aver_dt$cell_type <- rownames(aver_dt)

aver_dt <- aver_dt %>% separate(cell_type, c('cell_type', 'response'),sep="_P")

aver_dt$response <- paste("P",aver_dt$response,sep = "")
# aver_dt$response <- paste("PDCD1",aver_dt$response,sep = "")

df_p_val1 <- aver_dt %>% group_by(response)%>%
  wilcox_test(cell_type ~Average_Expression) %>%
  adjust_pvalue(p.col="p",method="bonferroni") %>%
  add_significance(p.col="p.adj") %>%
  add_xy_position(x="year",dodge=0.8)

PDCD1_plot<-ggplot(aver_dt, aes(cell_type,response)) +
  geom_tile(aes(fill =Average_Expression), colour = "white", size = 0.1)+
  scale_fill_gradientn(name='PDCD1 Mean\nscaled expression',
                       colours=colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100))+
  theme_minimal() +
  theme(axis.line=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = 8,color = 'black',
        ),
        axis.text.x = element_text(size = 7,color = 'black',angle =90,hjust =1,vjust=1),
        legend.title = element_text(size=7, color = "black"),
        legend.text = element_text(size=7,color = "black")) +
  scale_y_discrete(position = "right")



####Tcell细胞比例----

library(ggplot2)
library(dplyr)
library(ggalluvial)
Ratio <- Bcell4@meta.data %>% #Bcell4为你的单细胞seurat对象
  group_by(Sample, celltype) %>% # Sample为你的seurat对象meta信息中的样本编号信息，celltype为细胞类型信息
  summarise(n=n()) %>%
  mutate(relative_freq = n/sum(n))


cli=read.csv("D:\\leonard\\all_sampleid.csv",header = T)
cli_neu=cli[match(Ratio$Sample,cli$sample),]
Ratio_cli=cbind(Ratio,cli_neu)
head(Ratio_cli)

library(RColorBrewer)
library(ggpubr)
library(ggplot2)
library(cowplot)

Ratio_cli[1:4,1:4]
Exp=Ratio_cli[c(1,2,4,6)]

library(reshape2)
#dcast()函数将长格式转化为宽格式
data_wide_d<-dcast(Exp, Sample~Exp$celltype,
                   value.var = 'relative_freq')
group=Exp[match(data_wide_d$Sample,Exp$Sample),]
data_wide_d$group=group$HBV_state

Exp_plot <- data_wide_d
Exp_plot$group <- factor(Exp_plot$group, levels = c("HBV_n", "HBV_a", "HBV_i"))
rownames(Exp_plot)=Exp_plot$Sample
Exp_plot=Exp_plot[-1]



cellname=c(colnames(Exp_plot)[-length(Exp_plot)])
# 为每个不同的组设置颜色。
col <- c("#5CB85C","#D9534F", "#F0AD4E")
comparisons <- list(c("HBV_n", "HBV_a"),
                    c("HBV_n", "HBV_i"),
                    c("HBV_a", "HBV_i"))
plist <- list() # 创建一个空的列表 'plist'，用于存储后续循环中生成的图形
for (i in 1:length(cellname)) {
  bar_tmp <- Exp_plot[, c(cellname[i], "group")]  # 从 'Exp_plot' 中提取当前基因的表达信息和样本组
  colnames(bar_tmp) <- c("Frequency", "group")
  pb1 <- ggboxplot(bar_tmp, # Create a boxplot using ggboxplot.
                   x = "group", # X-axis is for groups.
                   y = "Frequency", # Y-axis is for expression levels.
                   color = "group", # Fill by sample group.
                   fill = NULL,
                   add = "jitter", # Add jitter points.
                   bxp.errorbar.width = 0.8,
                   width = 0.5,
                   size = 0.1,
                   font.label = list(size = 20),
                   palette = col) +
    theme(panel.background = element_blank())
  pb1 <- pb1 + theme(axis.line = element_line(colour = "black")) +
    theme(axis.title.x = element_blank()) # 调整坐标轴

  pb1 <- pb1 + theme(axis.title.y = element_blank()) +
    theme(axis.text.x = element_text(size = 15, angle = 45, vjust = 1, hjust = 1))

  pb1 <- pb1 + theme(axis.text.y = element_text(size = 15)) +
    ggtitle(gene[i]) +
    theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"))

  pb1 <- pb1 + theme(legend.position = "NA")  # 删除图例（因为样本类型已经显示在横轴上）

  pb1 <- pb1 + stat_compare_means(method = "t.test", hide.ns = FALSE,
                                  comparisons = comparisons,label = "p.signif",vjust=0.02,bracket.size=0.6) # 执行显著性测试，使用 t 检验，并添加不同组别之间的比较

  plist[[i]] <- pb1 # 将生成的图形存储在 'plist' 中
}
pdf("boxplot_Ex3.pdf", height = 12, width = 10)
# Align and arrange the plots into a grid.
plot_grid(plist[[1]], plist[[2]], plist[[3]],
          plist[[4]], plist[[5]], plist[[6]],
          plist[[7]], plist[[8]], plist[[9]],
          plist[[10]], plist[[11]], plist[[12]],
          plist[[13]], ncol = 4) # ncol = 4 indicates the number of columns in the grid.
dev.off()



