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

FeaturePlot(T_cells_object,features ="LEF1",reduction = "umap.harmony",raster=FALSE)

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
         "LTB")

gene<-c("PDCD1",
        "LAYN",
        "HAVCR2",
        "LAG3",
        "CTLA4",
        "TIGIT",
        "TOX",
        "RBPJ",
        "VCAM1")

gene<-c("GZMB",
        "MYO7A",
        "CD244",
        "VSIR",
        "BTLA",
        "ENTPD1",
        "CD160",
        "LAIR1")

gene<-c("GZMK",
        "TOP2A",
        "MKI67",
        "STMN1","HMGB2")

VlnPlot(object = T_cells_object, features = as.character(gene))

VlnPlot(object = T_cells_object, features = "TCF7")

FeaturePlot(T_cells_object,features =gene,reduction = "umap.harmony",raster=FALSE)



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

p1 <- DotPlot(T_cells_object, features =unique(tcell_marker$gene_name),
              assay='RNA' )
p1
sub_markers_count<- p1$data
names(sub_markers_count)
colnames(sub_markers_count)
sub_markers_count<- sub_markers_count %>% data.frame()
sub_markers_count$features.plot <- factor(sub_markers_count$features.plot,levels =rev(unique(sub_markers_count$features.plot)))
sub_markers_count$id <- factor(sub_markers_count$id,levels =marker_group_list)

colnames(sub_markers_count)<-c("avg.exp","pct.exp","marker","Tcell_group","value")
sub_markers_count$marker_group <- tcell_marker[match(sub_markers_count$marker,tcell_marker$gene_name),3]

sub_markers_count$cell_group <- plyr::mapvalues(sub_markers_count$Tcell_group,
                                              from = names(cell_type_list), 
                                              to = cell_type_list)



aver_dt_df_long <- sub_markers_count


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
FeaturePlot(object =T_cells_object,features = "LTB",reduction = "umap.harmony",raster=FALSE, cols = pal, order = T)


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
                                           elem_list_rect(fill="white",color = "black"),
                                         text_x = elem_list_text(size = 12))) +
  scale_fill_manual(values = c("#4974a4","#4dae47","#f29600"))+
  labs(x=NULL,y=NULL)+
  theme(axis.text.x=element_text(size=10,angle = 30,vjust=0.5,hjust=0.5,color="black"),
        axis.text.y=element_text(size=10,color="black"),
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

p_Cellratio

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/细胞比例/t细胞比例.jpg",p_Cellratio,width = 8,height = 5)
ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/细胞比例/t细胞比例.pdf",p_Cellratio,width = 8,height = 5)

####T细胞亚型比例计算----

Cellratio <- prop.table(table(T_cells_object@meta.data$cd4_cd8_group,T_cells_object@meta.data$orig.ident), margin = 2)

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

cell_type_groups =c("CD4_C1-PDCD1","CD4_C2-IL7R","CD4_C3-FOXP3", 
                    "CD4_C4-FOSB","CD8_C1-GNLY","CD8_C2-CTSW","CD8_C3-IL7R",  
                    "CD8_C4-GZMK","CD8_C5-STMN1","CD3+CD4-CD8-","NK")
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
cells_<-plot_grid(pplist[['CD4_C1-PDCD1']],
                  pplist[['CD4_C2-IL7R']],
                  pplist[['CD4_C3-FOXP3']],
                  pplist[['CD4_C4-FOSB']],
                  pplist[['CD8_C1-GNLY']],
                  pplist[['CD8_C2-CTSW']],
                  pplist[['CD8_C3-IL7R']],
                  pplist[['CD8_C4-GZMK']],
                  pplist[['CD8_C5-STMN1']],
                  pplist[['CD3+CD4-CD8-']],
                  pplist[['NK']],
                  align = "h",  #axis = 'l',
                  nrow =2)



ggsave('D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/细胞比例/T_cells_MPR_boxplot.pdf',cells_,width =12,height =8,limitsize = FALSE)

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
out.prefix <- "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/细胞比例/T_cell_OR2"

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

out.prefix <- "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/细胞比例/T cell_ROCR"
A <- do.tissueDist(cellInfo.tb = meta.tb,
                   meta.cluster = meta.tb$cd4_cd8_group,
                   colname.patient = "patient",
                   loc = meta.tb$CRNCR_Response2,
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

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/细胞比例/OR_plot_MPR.jpg",OR_plot,width=5,height=7)
ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/细胞比例/OR_plot_MPR.pdf",OR_plot,width=5,height=7)

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/细胞比例/OR_plot_MPR_treatments.jpg",OR_plot,width=6,height=7)

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/细胞比例/OR_plot_MPR_treatments.pdf",OR_plot,width=6,height=7)




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

####有无微生物T细胞差异分析----
MPRNMPR_object_miMPRobe_remove<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/MPRNMPR_object_miMPRobe_remove.rds")
MPRNMPR_object_miMPRobe_remove = UpdateSeuratObject(MPRNMPR_object_miMPRobe_remove)
MPRNMPR_object_miMPRobe_remove

Idents(MPRNMPR_object_miMPRobe_remove) <- "cell_type_new"
MPRNMPR_object_miMPRobe_remove@meta.data$MPR_Response2
T_cells_object=subset(MPRNMPR_object_miMPRobe_remove,idents = c("T cells"))

T_cells_object@meta.data$group_microbe2
Idents(T_cells_object) <- "group_microbe2"
T_cells_object_deg_all=FindMarkers(T_cells_object, ident.1 = c("Bacteria+"), ident.2 = "Bacteria-")
dim(T_cells_object_deg_all)
names(T_cells_object_deg_all)
write.csv(T_cells_object_deg_all,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/免疫细胞/T_cells_object_deg_all.csv")
T_cells_object_deg_all<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/免疫细胞/T_cells_object_deg_all.csv",header = T)

####根据过滤掉pt 小于0.1 的行

names(T_cells_object_deg_all)

T_cells_object_deg_all_filter<- T_cells_object_deg_all %>% filter(pct.1>0.1 &pct.2>0.1)

####读取免疫细胞相关基因，从差异基因集中筛选

gene_name<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/免疫细胞/免疫细胞gene_name.csv")

tex_gene_name<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/有无微生物T细胞功能/t_gene_name.csv")

intersection <- intersect(gene_name$Gene.symbol, T_cells_object_deg_all_filter$X)
intersection <- intersect(tex_gene_name$Chemokine.Chemokine.receptor, T_cells_object_deg_all_filter$X)

T_cells_object_deg_all_filter_inter<- T_cells_object_deg_all_filter %>% filter(T_cells_object_deg_all_filter$X %in%intersection )

df <- T_cells_object_deg_all_filter_inter
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
  ggtitle("T cells DEG (Bacteria+ vs Bacteria-)") + #标题
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
ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/免疫细胞/T细胞火山图.jpg",degp,width = 10,height=7)
ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/免疫细胞/T细胞火山图.pdf",degp,width = 10,height=7)


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

####T细胞有无微生物GESA分析----
library(gggsea)
library(ggplot2)
library(GSVA)
library(enrichplot)
library(GSEABase)
library(fgsea)
library(clusterProfiler)
######GSEA_h.all.v2024.1.Hs.symbols.gmt####

##ID转换（将gene_symbol转为ENTREZID）
deg<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/免疫细胞/T_cells_object_deg_all.csv",header = T)
names(deg)

colnames(deg)[1] <- "gene_name"
gmtfile <- "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/h.all.v2024.1.Hs.symbols.gmt"
hallmark <- read.gmt(gmtfile)
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


# gsea.re2 <- fgseae(pathways = hallmark.list,
#                   stats = id,
#                   minSize=1,
#                   maxSize=10000,
#                   nperm=1000)

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
####可视化
library(ggsci)
col_gsea1<-pal_simpsons()(16)

num1=1
gseaplot2(gsea.re1,geneSetID = rownames(g1)[1:num1],
          title = "",#标题
          color = col_gsea1[1:num1],#颜色
          base_size = 14,#基础大小
          rel_heights = c(1, 0.2, 0.4),#小图相对高度
          subplots = 1:3,#展示小图
          pvalue_table = FALSE,#p值表格
          ES_geom = "line"#line or dot
)

####设置多个名称
num2=3
gseaplot2(gsea.re1,geneSetID = rownames(g1)[1:num2],
          title = "",#标题
          color = col_gsea1[1:num2],#颜色
          base_size = 14,#基础大小
          rel_heights = c(1, 0.2, 0.4),#小图相对高度
          subplots = 1:3,#展示小图
          pvalue_table = FALSE,#p值表格
          ES_geom = "line"#line or dot
)

####plotsea 可视化

plotGseaTable(hallmark.list[g2$pathway],
              id, 
              gsea.re2,gseaParam = 0.5,
              colwidths = c(0.5,0.2,0.1,0.1,0.1)
)

####ggsea 可视化
names(hallmark.list)
se_hall<-c(head(g2$pathway,9),tail(g1$pathway,9))

#3 所选中的基因集
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
            linesize=2,lty=1,
            #多个图形时使用，设定每行列图形个数
            nrow = 6,
            ncol = 2
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
ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/免疫细胞/免疫细胞GSEA/免疫细胞hallmark_GSEA.pdf",p_nes,width = 8,height=10)

####Go数据库----

library(dplyr)
library(org.Hs.eg.db) #物种注释包(Homo sapiens)
library(clusterProfiler) #富集分析
library(enrichplot)
library(ComplexHeatmap)
library(circlize)
gene <- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/免疫细胞/T_cells_object_deg_all.csv",header = T)
colnames(gene)[1] <- "gene_name"
gene$regulate <-ifelse(gene$avg_log2FC>0,"up","down")
gene_ID=bitr(gene$gene_name,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db")
colnames(gene_ID)
View(gene_ID)
table(is.na(gene_ID))
gene_df <- base::merge(gene_ID,gene,by.x="SYMBOL",by.y="gene_name")

#按照LOG_FOIDCHANGE进行排序
colnames(gene_df)
gene_df_sort <- gene_df[order(gene_df$avg_log2FC, decreasing = T),]
gene_fc <- gene_df_sort$avg_log2FC
names(gene_fc) <- gene_df_sort$ENTREZID

#GSEA分析——以GO数据库为例

Go_gseresult <- gseGO(gene_fc, 'org.Hs.eg.db', keyType = "ENTREZID", ont="all", 
                      minGSSize = 10, maxGSSize = 1000, pvalueCutoff=0.05,
                      pAdjustMethod = "BH",
                      verbose = FALSE,
                      eps = 0)

Go_ges <- setReadable(Go_gseresult,
                      OrgDb = org.Hs.eg.db,
                      keyType = "ENTREZID")

Go_ges_result <- Go_ges@result
Go_ges_result$Description

gseaplot2(Go_gseresult,18,pvalue_table = TRUE)
View(Go_ges)
my_colors <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b")
####


####p8p9正式画图----
library(GseaVis)
#指定名称：图纸
selected_gene_sets <- c(5, 8, 17)
p1<-gseaNb(object = Go_ges,
           geneSetID =selected_gene_sets)

View(Go_ges)
my_colors <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b")
# clsaasic with pvalue
p8<-gseaNb(object = Go_ges,
           geneSetID = selected_gene_sets,
           curveCol=c("#1f77b4", "#ff7f0e", "#2ca02c"),
           addPval = F,
           pvalX = 0.75,pvalY = 0.8,
           pCol = 'black',
           pHjust = 0,
           subPlot = 2)
p8
selected_gene_sets=c(25, 57, 66)
p9<-gseaNb(object = Go_ges,
           geneSetID = selected_gene_sets,
           curveCol=c("#d62728", "#9467bd", "#d3a193"),
           addPval = F,
           pvalX = 0.75,pvalY = 0.8,
           pCol = 'black',
           pHjust = 0,
           subPlot = 2)
p9
##############柱状图############################3
library(stringr)
library(ggplot2)
Go_ges_result <- Go_ges@result
Go_ges_result$Description
Go_ges_result_df <- cbind(ONTOLOGY=Go_ges_result$ONTOLOGY,ID=Go_ges_result$ID,
                          Description=Go_ges_result$Description,setSize=Go_ges_result$setSize,
                          enrichmentScore=Go_ges_result$enrichmentScore,NES=Go_ges_result$NES,
                          p.adjust=Go_ges_result$pvalue,Go_ges_result$p.adjust,
                          rank=Go_ges_result$qvalue,Go_ges_result$rank) %>% as.data.frame()
Go_ges_result_df$NES <- as.numeric(as.character(Go_ges_result_df$NES))
Go_ges_result_df <- Go_ges_result_df[order(Go_ges_result_df$NES),]
Go_ges_result_df$yax <- ifelse(Go_ges_result_df$NES >0, -0.02, 0.02)
Go_ges_result_df$col <- ifelse(Go_ges_result_df$NES > 0, "blue","red")
Go_ges_result_df$NES
write.csv(Go_ges_result_df,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/有无微生物T细胞功能/Go_ges_result_df柱状图数据.csv")
View(Go_ges_result_df)
Go_ges_result_df_filter<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/有无微生物T细胞功能/Go_ges_result_df柱状图数据_filter.csv")
p <- 
  ggplot(Go_ges_result_df_filter,aes(y=NES,x=reorder(Description2,NES),label = Description2))+
  geom_text(aes(y = yax,hjust = NES > 0),size=4)+ 
  geom_bar(stat="identity",aes(fill=col),width = 0.8)+
  geom_hline(yintercept = 0)+
  coord_flip(ylim = c(-2.5,2.5))+
  labs(y = "Normalised enrichment score", x = "",title="")+
  scale_fill_manual(values = c("#FB6542","#375E97"),name="Expression",labels=c("Up","Down"))+
  scale_x_discrete(breaks = NULL)+
  theme_classic(base_size =  12)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.y = element_blank(),
        legend.position = 'bottom')
p
ggsave(p,file="D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/有无微生物T细胞功能/gsea_t_cell_柱状图.pdf",width = 5,height = 5)



####T细胞亚群组织偏好性----
View(CD8_tex_object)
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
out.prefix <- "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/不同分组T亚群拟时序分析/CD8_OR_NCR"

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

meta.tb <- CD8_tex_object@meta.data

names(meta.tb)
View(meta.tb)

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

#自己做图
data <- A$count.dist.melt.ext.tb
data$rid <- factor(data$rid, levels=c("CD4_C1-PDCD1","CD4_C2-IL7R","CD4_C3-FOXP3","CD4_C4-FOSB",
                                      "CD8_C1-GNLY","CD8_C2-CTSW","CD8_C3-IL7R","CD8_C4-GZMK","CD8_C5-STMN1",
                                      "CD3+CD4-CD8-","NK" ))

data$rid <- factor(data$rid, levels=c("CD8_C1-GNLY","CD8_C5-STMN1" ))
head(data)
bl <- colorRampPalette(c("navy","royalblue","lightskyblue"))(100)                      
re <- colorRampPalette(c("mistyrose", "red2","darkred"))(100)
ex_color=c("#ffffff","#6238ff","#ff220e")
library(viridis)

#colorRampPalette(rev(brewer.pal(n = 5, name ="RdYlBu")))(100))
library(RColorBrewer)
OR_plot<-ggplot(data, aes(rid,cid)) + 
  geom_tile(aes(fill = OR), colour = "gray80", size = 1)+
  scale_fill_gradientn(name='OR',
                       colours=colorRampPalette(ex_color)(100))+
                       # limits = c(0,100),  # 颜色映射到数据值范围 0.5 - 2
                       # values = scales::rescale(c(0, 5, 40, 80,100))) + # 控制颜色对应的数值位置)+
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
  scale_y_discrete(position = "right")+
  theme(
    legend.position = "top",             # 将图例放在顶部
    legend.direction = "horizontal"      # 将图例横向排列
  )
OR_plot

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/不同分组T亚群拟时序分析和特征得分/CD8_OR_plot_CRNCR.jpg",OR_plot,width=3,height=7)

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/不同分组T亚群拟时序分析和特征得分/CD8_OR_plot_CRNCR.pdf",OR_plot,width=3,height=7)


####T Cell Exhuastion signature and T Cell Mediated Immune Response to Tumor Cell T cell Cytotoxic signature----
T_cells_object<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/T_cells_object.rds")
T_cells_object = UpdateSeuratObject(T_cells_object)
T_cells_object

###读取耗竭Marker---
T_exhausted<-read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/不同分组T亚群拟时序分析和特征得分/T细胞耗竭毒性所用marker.csv")
###毒性marker
# T_exhausted <-c("PDCD1","LAYN","HAVCR2","LAG3","CTLA4","TIGIT","TOX","RBPJ","VCAM1","GZMB","MYO7A","CD244","VSIR","BTLA","ENTPD1","CD160","LAIR1")
# cytotoxicity <-c("GZMA","GZMB","GZMH","GZMK","GZMH","GNLY","PRF1","IFNG","TNF","SERPINB1",
# "SERPINB6","SERPINB9","CTSA","CTSB","CTSC","CTSD","CTSW","CST3","CST7",	
# "CSTB","LAMP1","LAMP3","CAPN2","KLRK1","KLRK1","KLRB1","NKG7")
# Mediated_Immune_Response<-c("HLA-A","HLA-DRB1","HLA-DRB3","MR1","HMGB1","HSPD1","FBXO38","SLC22A13")
# Tex_progenitor<-c("IL7R","GPR183","LMNA","NR4A3","TCF7","MGAT4A","CD55",	
#                   "AIM1","PER1","FOSL2","EGR1","TSPYL2","YPEL5","CSRNP1",
#                   "REL","SKIL","PIK3R1","FOXP1","RGCC","PFKFB3","MYADM",
#                   "ZFP36L2","USP36","TC2N","FAM177A1","BTG2","TSC22D2",
#                   "FAM65B","STAT4","RGPD5","NEU1","IFRD1","PDE4B","NR4A1")

T_cell_markers<-list(T_exhausted$T.Cell.Exhuastion.signature,T_exhausted$T.Cell.Mediated.Immune.Response.to.Tumor.Cell,
                     T_exhausted$Terminally.Exhausted.CD8..T.Cell.Signature,T_exhausted$Terminally.Exhausted.CD8..T.Cell.Signature_filter,T_exhausted$Progenitor.Exhausted.CD8..T.Cell.Signature,
                     T_exhausted$Cytotoxicity,T_exhausted$Progenitor.Exhausted.CD8..T.Cell.Signature_filter)

imm_T<-AddModuleScore(T_cells_object,features=T_cell_markers,name=c("T.Cell.Exhuastion.signature",
                                                                    "T.Cell.Mediated.Immune.Response.to.Tumor.Cell",
                                                                    "Terminally.Exhausted.CD8..T.Cell.Signature",
                                                                    "Terminally.Exhausted.CD8..T.Cell.Signature_filter",
                                                                    "Progenitor.Exhausted.CD8..T.Cell.Signature",
                                                                    "Cytotoxicity","Progenitor.Exhausted.CD8..T.Cell.Signature_filter"),nbin=12)
imm_T$T.Cell.Exhuastion.signature1
imm_T$T.Cell.Mediated.Immune.Response.to.Tumor.Cell2
imm_T$Terminally.Exhausted.CD8..T.Cell.Signature3
imm_T$Terminally.Exhausted.CD8..T.Cell.Signature_filter4
imm_T$Progenitor.Exhausted.CD8..T.Cell.Signature5
imm_T$Cytotoxicity6
imm_T$Progenitor.Exhausted.CD8..T.Cell.Signature_filter7


####提取数据集画boxplot 并加入显著性标志----
library(tidyverse)
library(gghalves)
library(rstatix)
library(ggpubr)
library(ggsci)

sig_score<- cbind(imm_T$cd4_cd8_group,imm_T$MPR_Response2,imm_T$T.Cell.Exhuastion.signature1,
                  imm_T$T.Cell.Mediated.Immune.Response.to.Tumor.Cell2,
                  imm_T$Terminally.Exhausted.CD8..T.Cell.Signature3,
                  imm_T$Terminally.Exhausted.CD8..T.Cell.Signature_filter4,
                  imm_T$Progenitor.Exhausted.CD8..T.Cell.Signature5,
                  imm_T$Cytotoxicity6,
                  imm_T$Progenitor.Exhausted.CD8..T.Cell.Signature_filter7) %>% as.data.frame()

names(sig_score)<- c("cd4_cd8_group","MPR_Response2","T.Cell.Exhuastion.signature1",
                     "T.Cell.Mediated.Immune.Response.to.Tumor.Cell2",
                     "Terminally.Exhausted.CD8..T.Cell.Signature3",
                     "Terminally.Exhausted.CD8..T.Cell.Signature_filter4",
                     "Progenitor.Exhausted.CD8..T.Cell.Signature5",
                     "Cytotoxicity6",
                     "Progenitor.Exhausted.CD8..T.Cell.Signature_filter7")

df<-sig_score
df$T.Cell.Exhuastion.signature1 <- as.numeric(df$T.Cell.Exhuastion.signature1)
df$T.Cell.Mediated.Immune.Response.to.Tumor.Cell2 <- as.numeric(df$T.Cell.Mediated.Immune.Response.to.Tumor.Cell2)
df$Terminally.Exhausted.CD8..T.Cell.Signature3 <- as.numeric(df$Terminally.Exhausted.CD8..T.Cell.Signature3)
df$Terminally.Exhausted.CD8..T.Cell.Signature_filter4 <- as.numeric(df$Terminally.Exhausted.CD8..T.Cell.Signature_filter4)
df$Progenitor.Exhausted.CD8..T.Cell.Signature5 <- as.numeric(df$Progenitor.Exhausted.CD8..T.Cell.Signature5)
df$Cytotoxicity6 <- as.numeric(df$Cytotoxicity6)
df$Progenitor.Exhausted.CD8..T.Cell.Signature_filter7 <- as.numeric(df$Progenitor.Exhausted.CD8..T.Cell.Signature_filter7)

gene_list <-as.vector(T_exhausted$Terminally.Exhausted.CD8..T.Cell.Signature[1:5])
FeaturePlot(object=CD8_tex_object,features =c("TCF7"),reduction = "umap.harmony",raster=FALSE, order = T)+
  xlab("UMAP_1")+
  ylab("UMAP_1")+
  theme(legend.position = "right")

Idents(CD8_tex_object)<-"cd4_cd8_group"

FeaturePlot(object=CD8_tex_object,features =c("ZNF683",
                                              "HOPX",
                                              "ITGAE"),reduction = "umap.harmony",raster=FALSE, order = T)+
  xlab("UMAP_1")+
  ylab("UMAP_1")+
  theme(legend.position = "right")

FeaturePlot(object=CD8_tex_object,features =c("PDCD1","HAVCR2","SLAMF6","SLAMF1"),reduction = "umap.harmony",raster=FALSE, order = T)+
  xlab("UMAP_1")+
  ylab("UMAP_1")+
  theme(legend.position = "right")


VlnPlot(object = CD8_tex_object, features =as.vector(T_exhausted$Terminally.Exhausted.CD8..T.Cell.Signature[1:5]))
DotPlot(object = CD8_tex_object, features =as.vector(T_exhausted$Terminally.Exhausted.CD8..T.Cell.Signature[20:25]))

DotPlot(object = CD8_tex_object, features =c("TBX21","PDCD1","HAVCR2","LEF1","TCF7","CD69","SLAMF6","SLAMF1","CD244"))


# df <- read_tsv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/不同分组T亚群拟时序分析/ggplot2优雅绘制半边箱线图20240427/data.xls") %>% 
#   mutate(year=as.character(year))

# stat_compare_means(aes(group = microbe_group), label = "p.signif",label.y=22, method="wilcox.test",size=7)
####耗竭boxplot----
df_p_val1 <- df %>%
  wilcox_test(T.Cell.Exhuastion.signature1~ MPR_Response2) %>%
  adjust_pvalue(p.col = "p", method = "bonferroni") %>%
  add_significance(p.col = "p.adj") %>% 
  add_xy_position(x = "MPR_Response2", dodge = 0.8) 
treatment_color <- c("#4974a4","#4dae47","#f29600")

df$MPR_Response2<- factor(df$MPR_Response2,levels = c("Pre_NMPR","Pre_MPR","Post_MPR"))
t_exhausted<- df %>%
  ggplot(aes(MPR_Response2,T.Cell.Exhuastion.signature1)) +
  geom_half_boxplot(fill=treatment_color,color="black",side="l",errorbar.draw = T,
                    outlier.shape = NA,width=0.8,lwd= 0.7) +
  geom_half_point(color=treatment_color,side = "r",alpha=0.1,
                  transformation_params = list(height = 0,width = 0.001,seed = 2))+
  stat_pvalue_manual(df_p_val1,label = "p.adj.signif",label.size=5,hide.ns = F,bracket.size = 0.5,y.position = c(1.3,1.5,1.4))+
  #facet_wrap(.~continent,nrow=1)+
  #scale_y_continuous(limits = c(0,5),breaks = seq(0,2,5))+
  scale_fill_npg()+
  scale_color_npg()+
  labs(title = "T cells exhausted signature",y='Signature score',x="")+
  theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),units=,"cm"),
        plot.title = element_text(hjust = 0.5, size = 13),
        panel.background = element_blank(),
        plot.background = element_blank(),
        strip.text = element_text(size=12),
        axis.title.y = element_text(size = 14),
        axis.line = element_line(color = "black",size = 0.4),
        axis.text.y = element_text(color="black",size=12),
        axis.text.x = element_text(margin = margin(t = 2),color="black",size=12),
        legend.position = "none",
        panel.spacing = unit(0,"lines"),
        axis.ticks.x=element_line(color="black",size=0.5,lineend = 1),
        axis.line.y=element_line(linetype=1,color="black",size=0.7),
        axis.ticks.y=element_line(color="black",size=0.5,lineend = 1),
        axis.line.x=element_line(linetype=1,color="black",size=0.7))+
  
  coord_cartesian()

t_exhausted

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/不同分组T亚群拟时序分析/T cells exhausted signature.pdf",t_exhausted,width = 6,height = 5)

####毒力boxplot----
df_p_val2 <- df %>%
  wilcox_test(Cytotoxicity6~ MPR_Response2) %>%
  adjust_pvalue(p.col = "p", method = "bonferroni") %>%
  add_significance(p.col = "p.adj") %>% 
  add_xy_position(x = "MPR_Response2", dodge = 0.8) 

t_cytotoxicity<- df %>%
  ggplot(aes(MPR_Response2,Cytotoxicity6)) +
  geom_half_boxplot(fill=treatment_color,color="black",side="l",errorbar.draw = T,
                    outlier.shape = NA,width=0.8,lwd= 0.7) +
  geom_half_point(color=treatment_color,side = "r",alpha=0.1,
                  transformation_params = list(height = 0,width = 0.001,seed = 2))+
  stat_pvalue_manual(df_p_val2,label = "p.adj.signif",label.size=5,hide.ns = F,bracket.size = 0.5,y.position = c(1.3,1.5,1.4))+
  #facet_wrap(.~continent,nrow=1)+
  #scale_y_continuous(limits = c(0,5),breaks = seq(0,2,5))+
  scale_fill_npg()+
  scale_color_npg()+
  labs(title = "T cells cytotoxic signature",y='Signature score',x="")+
  theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),units=,"cm"),
        plot.title = element_text(hjust = 0.5, size = 13),
        panel.background = element_blank(),
        plot.background = element_blank(),
        strip.text = element_text(size=12),
        axis.title.y = element_text(size = 14),
        axis.line = element_line(color = "black",size = 0.4),
        axis.text.y = element_text(color="black",size=12),
        axis.text.x = element_text(margin = margin(t = 2),color="black",size=12),
        legend.position = "none",
        panel.spacing = unit(0,"lines"),
        axis.ticks.x=element_line(color="black",size=0.5,lineend = 1),
        axis.line.y=element_line(linetype=1,color="black",size=0.7),
        axis.ticks.y=element_line(color="black",size=0.5,lineend = 1),
        axis.line.x=element_line(linetype=1,color="black",size=0.7))+
  
  coord_cartesian()
t_cytotoxicity

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/不同分组T亚群拟时序分析/T cells cytotoxicity signature.pdf",t_cytotoxicity,width = 6,height = 5)

####免疫boxplot----
df_p_val3 <- df %>%
  wilcox_test(T.Cell.Mediated.Immune.Response.to.Tumor.Cell2~ MPR_Response2) %>%
  adjust_pvalue(p.col = "p", method = "bonferroni") %>%
  add_significance(p.col = "p.adj") %>% 
  add_xy_position(x = "MPR_Response2", dodge = 0.8) 

t_Immune<- df %>%
  ggplot(aes(MPR_Response2,T.Cell.Mediated.Immune.Response.to.Tumor.Cell2)) +
  geom_half_boxplot(fill=treatment_color,color="black",side="l",errorbar.draw = T,
                    outlier.shape = NA,width=0.8,lwd= 0.7) +
  geom_half_point(color=treatment_color,side = "r",alpha=0.1,
                  transformation_params = list(height = 0,width = 0.001,seed = 2))+
  stat_pvalue_manual(df_p_val3,label = "p.adj.signif",label.size=5,hide.ns = F,bracket.size = 0.5,y.position = c(1.7,1.9,1.8))+
  #facet_wrap(.~continent,nrow=1)+
  #scale_y_continuous(limits = c(0,5),breaks = seq(0,2,5))+
  scale_fill_npg()+
  scale_color_npg()+
  labs(title = "T cells mediated immune \n response to tumor cell",y='Signature score',x="")+
  theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),units=,"cm"),
        plot.title = element_text(hjust = 0.5, size = 13),
        panel.background = element_blank(),
        plot.background = element_blank(),
        strip.text = element_text(size=12),
        axis.title.y = element_text(size = 14),
        axis.line = element_line(color = "black",size = 0.4),
        axis.text.y = element_text(color="black",size=12),
        axis.text.x = element_text(margin = margin(t = 2),color="black",size=12),
        legend.position = "none",
        panel.spacing = unit(0,"lines"),
        axis.ticks.x=element_line(color="black",size=0.5,lineend = 1),
        axis.line.y=element_line(linetype=1,color="black",size=0.7),
        axis.ticks.y=element_line(color="black",size=0.5,lineend = 1),
        axis.line.x=element_line(linetype=1,color="black",size=0.7))+
  
  coord_cartesian()
t_Immune
ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/不同分组T亚群拟时序分析/T cells t_Immune signature.pdf",t_Immune,width = 6,height = 5)

Tcell_type_color<-c("#DF66B0","#9888DB","#DC8062","#D8ACC9","#83B5D3","#D5D8C5", "#77DCCF", "#BA4EE4",
                    "#AEE94F","#89DB8E","#E0D175")

####耗竭细胞类型boxplot----

use_colors<-c("CD4_C1-PDCD1"="#DF66B0", 
  "CD4_C2-IL7R"="#9888DB", 
  "CD4_C3-FOXP3"="#DC8062", 
  "CD4_C4-FOSB"="#D8ACC9",
  "CD8_C1-GNLY"="#83B5D3",
  "CD8_C2-CTSW"="#D5D8C5", 
  "CD8_C3-IL7R"="#77DCCF", 
  "CD8_C4-GZMK"="#BA4EE4", 
  "CD8_C5-STMN1"="#AEE94F", 
  "CD3+CD4-CD8-"="#89DB8E",
  "NK"="#E0D175")

df_p_val4 <- df %>%
  wilcox_test(T.Cell.Exhuastion.signature1~ cd4_cd8_group) %>%
  adjust_pvalue(p.col = "p", method = "bonferroni") %>%
  add_significance(p.col = "p.adj")%>% 
  add_xy_position(x = "cd4_cd8_group", dodge = 0.8)

df$cd4_cd8_group<- factor(df$cd4_cd8_group,levels=c("CD4_C1-PDCD1","CD4_C2-IL7R","CD4_C3-FOXP3","CD4_C4-FOSB",
                                                     "CD8_C1-GNLY","CD8_C2-CTSW","CD8_C3-IL7R","CD8_C4-GZMK","CD8_C5-STMN1",
                                                     "CD3+CD4-CD8-","NK" ))
unique(df$cd4_cd8_group)

t_exhausted_cell_type<- df %>%
  ggplot(aes(cd4_cd8_group,T.Cell.Exhuastion.signature1)) +
  # geom_half_boxplot(fill=Tcell_type_color,color="black",side="l",errorbar.draw = T,
  #                   outlier.shape = NA,width=0.8,lwd= 0.7) +
  # geom_half_point(color=Tcell_type_color,side = "r",alpha=0.1,
  #                 transformation_params = list(height = 0,width = 0.001,seed = 2))+
  #stat_pvalue_manual(df_p_val4,label = "p.adj.signif",label.size=5,hide.ns = T,bracket.size = 0.5)+
  #facet_wrap(.~continent,nrow=1)+
  #scale_y_continuous(limits = c(0,5),breaks = seq(0,2,5))+
  geom_half_boxplot(aes(fill=cd4_cd8_group),color="black",side="l",errorbar.draw = T,
                    outlier.shape = NA,width=0.8,lwd= 0.5) +
  geom_half_point(aes(),side = "r",alpha=0.1,size=0.1,
                  transformation_params = list(height = 0,width = 0.001,seed = 2))+
  scale_fill_manual(values= use_colors)+
  facet_wrap(.~MPR_Response2,ncol=1)+

  scale_fill_npg()+
  scale_color_npg()+
  labs(title = "T cells exhausted signature",y='Signature score',x="")+
  theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),units=,"cm"),
        plot.title = element_text(hjust = 0.5, size = 13),
        panel.background = element_blank(),
        plot.background = element_blank(),
        strip.text = element_text(size=12),
        axis.title.y = element_text(size = 14),
        axis.line = element_line(color = "black",size = 0.4),
        axis.text.y = element_text(color="black",size=12),
        axis.text.x = element_text(color="black",size=12,angle=30,vjust = 1,hjust = 1),
        legend.position = "none",
        panel.spacing = unit(0,"lines"),
        axis.ticks.x=element_line(color="black",size=0.5,lineend = 1),
        axis.line.y=element_line(linetype=1,color="black",size=0.7),
        axis.ticks.y=element_line(color="black",size=0.5,lineend = 1),
        axis.line.x=element_line(linetype=1,color="black",size=0.7))+
  
  coord_cartesian()

t_exhausted_cell_type

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/耗竭特征得分/T_cells_exhausted_cell_type.pdf",t_exhausted_cell_type,width = 7,height = 7)


####毒力细胞类型----
df_p_val5 <- df %>%
  wilcox_test(Cytotoxicity6~ cd4_cd8_group) %>%
  adjust_pvalue(p.col = "p", method = "bonferroni") %>%
  add_significance(p.col = "p.adj") %>% 
  add_xy_position(x = "cd4_cd8_group", dodge = 0.8) 

t_cytotoxicity_cell_type<- df %>%
  ggplot(aes(cd4_cd8_group,Cytotoxicity6)) +
  geom_half_boxplot(fill=Tcell_type_color,color="black",side="l",errorbar.draw = T,
                    outlier.shape = NA,width=0.8,lwd= 0.7) +
  geom_half_point(color=Tcell_type_color,side = "r",alpha=0.1,
                  transformation_params = list(height = 0,width = 0.001,seed = 2))+
  #stat_pvalue_manual(df_p_val5,label = "p.adj.signif",label.size=5,hide.ns = F,bracket.size = 0.5,y.position = c(1.3,1.5,1.4))+
  #facet_wrap(.~continent,nrow=1)+
  #scale_y_continuous(limits = c(0,5),breaks = seq(0,2,5))+
  scale_fill_npg()+
  scale_color_npg()+
  labs(title = "T cells cytotoxic signature",y='Signature score',x="")+
  theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),units=,"cm"),
        plot.title = element_text(hjust = 0.5, size = 13),
        panel.background = element_blank(),
        plot.background = element_blank(),
        strip.text = element_text(size=12),
        axis.title.y = element_text(size = 14),
        axis.line = element_line(color = "black",size = 0.4),
        axis.text.y = element_text(color="black",size=12),
        axis.text.x = element_text(margin = margin(t = 2),color="black",size=12,angle=30,vjust = 1,hjust = 1),
        legend.position = "none",
        panel.spacing = unit(0,"lines"),
        axis.ticks.x=element_line(color="black",size=0.5,lineend = 1),
        axis.line.y=element_line(linetype=1,color="black",size=0.7),
        axis.ticks.y=element_line(color="black",size=0.5,lineend = 1),
        axis.line.x=element_line(linetype=1,color="black",size=0.7))+
  
  coord_cartesian()
t_cytotoxicity_cell_type
ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/不同分组T亚群拟时序分析/T cells cytotoxicity cell_type.pdf",t_cytotoxicity_cell_type,width = 7,height = 5)

####免疫细胞类型boxplot----
df_p_val6 <- df %>%
  wilcox_test(T.Cell.Mediated.Immune.Response.to.Tumor.Cell2~ cd4_cd8_group) %>%
  adjust_pvalue(p.col = "p", method = "bonferroni") %>%
  add_significance(p.col = "p.adj") %>% 
  add_xy_position(x = "MPR_Response2", dodge = 0.8) 

t_Immune_cell_type<- df %>%
  ggplot(aes(cd4_cd8_group,T.Cell.Mediated.Immune.Response.to.Tumor.Cell2)) +
  # geom_half_boxplot(fill=Tcell_type_color,color="black",side="l",errorbar.draw = T,
  #                   outlier.shape = NA,width=0.8,lwd= 0.7) +
  # geom_half_point(color=Tcell_type_color,side = "r",alpha=0.1,
  #                 transformation_params = list(height = 0,width = 0.001,seed = 2))+
  geom_half_boxplot(aes(fill = cd4_cd8_group),color="black",side="l",errorbar.draw = T,
                    outlier.shape = NA,width=0.8,lwd= 0.5) +
  geom_half_point(aes(),side = "r",alpha=0.1,size=0.1,
                  transformation_params = list(height = 0,width = 0.001,seed = 2))+
  #stat_pvalue_manual(df_p_val5,label = "p.adj.signif",label.size=5,hide.ns = F,bracket.size = 0.5,y.position = c(1.7,1.9,1.8))+
  facet_wrap(.~MPR_Response2,ncol=1)+
  #scale_y_continuous(limits = c(0,5),breaks = seq(0,2,5))+
  scale_fill_manual(values= use_colors)+
  scale_fill_npg()+
  scale_color_npg()+
  labs(title = "T cells mediated immune \n response to tumor cell",y='Signature score',x="")+
  theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),units=,"cm"),
        plot.title = element_text(hjust = 0.5, size = 13),
        panel.background = element_blank(),
        plot.background = element_blank(),
        strip.text = element_text(size=12),
        axis.title.y = element_text(size = 14),
        axis.line = element_line(color = "black",size = 0.4),
        axis.text.y = element_text(color="black",size=12),
        axis.text.x = element_text(margin = margin(t = 2),color="black",size=12,angle=30,vjust = 1,hjust = 1),
        legend.position = "none",
        panel.spacing = unit(0,"lines"),
        axis.ticks.x=element_line(color="black",size=0.5,lineend = 1),
        axis.line.y=element_line(linetype=1,color="black",size=0.7),
        axis.ticks.y=element_line(color="black",size=0.5,lineend = 1),
        axis.line.x=element_line(linetype=1,color="black",size=0.7))+
  
  coord_cartesian()
t_Immune_cell_type
ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/不同分组T亚群拟时序分析/T cells t_Immune cell type.pdf",t_Immune_cell_type,width = 7,height = 5)
ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/不同分组T亚群拟时序分析/T cells t_Immune cell type分面.pdf",t_Immune_cell_type,width = 7,height =7)

####早期耗竭和晚期耗竭boxplot----
df$PT_ratio<-df$Progenitor.Exhausted.CD8..T.Cell.Signature4/df$Terminally.Exhausted.CD8..T.Cell.Signature3
names(sig_score)<- c("cd4_cd8_group","MPR_Response2","T.Cell.Exhuastion.signature1",
                     "T.Cell.Mediated.Immune.Response.to.Tumor.Cell2",
                     "Terminally.Exhausted.CD8..T.Cell.Signature3",
                     "Progenitor.Exhausted.CD8..T.Cell.Signature4",
                     "Cytotoxicity5")
View(df)
names(df)
treatment_color <- c("#4974a4","#4dae47","#f29600")
df1<- df %>% filter(cd4_cd8_group %in% c("CD8_C1-GNLY", "CD8_C5-STMN1"))

df1<- df %>% filter(cd4_cd8_group %in% c("CD4_C1-PDCD1"))
names(df1)
df_p_val6 <- df1 %>% group_by(cd4_cd8_group) %>%
  wilcox_test(Terminally.Exhausted.CD8..T.Cell.Signature3 ~ MPR_Response2) %>%
  adjust_pvalue(p.col = "p", method = "bonferroni") %>%
  add_significance(p.col = "p.adj") %>% 
  add_xy_position(x = "MPR_Response2", dodge = 0.8)

df_p_val6_1 <- df1 %>% 
  wilcox_test(Terminally.Exhausted.CD8..T.Cell.Signature3 ~ cd4_cd8_group) %>%
  adjust_pvalue(p.col = "p", method = "bonferroni") %>%
  add_significance(p.col = "p.adj") %>% 
  add_xy_position(x = "MPR_Response2", dodge = 0.8)
#df1$cd4_cd8_group<- factor(df1$cd4_cd8_group,levels=c("CD8_C5-STMN1","CD8_C1-GNLY"))
df1$MPR_Response2<- factor(df1$MPR_Response2,levels=c("Pre_NMPR","Pre_MPR","Post_MPR"))
TEM_plot<- df1 %>%
  ggplot(aes(MPR_Response2,Terminally.Exhausted.CD8..T.Cell.Signature3)) +
  geom_half_boxplot(aes(fill = MPR_Response2),color="black",side="l",errorbar.draw = T,
                    outlier.shape = NA,width=0.8,lwd= 0.5) +
  geom_half_point(aes(),side = "r",
                  transformation_params = list(height = 0,width = 0.001,seed = 2))+
  stat_pvalue_manual(df_p_val6,label = "p.adj.signif",label.size=5,hide.ns = F)+
  facet_wrap(.~cd4_cd8_group,nrow=1)+
  scale_fill_manual(values= c("Pre_NMPR" = "#4974a4",
                              "Pre_MPR" = "#4dae47",
                              "Post_MPR"="#f29600"))+
  theme_minimal()+
  #scale_y_continuous(limits = c(0,95),breaks = seq(0,95,20))+
  #scale_fill_npg()+
  scale_color_npg()+
  labs(x=NULL,y=NULL)+
  #labs(title = "Progenitor exhausted CD8 T cell signature",y='Signature score',x="")+
  labs(title = "Terminally exhausted CD8 T cell signature",y='Signature score',x="")+
  theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),units=,"cm"),
        plot.title = element_text(hjust = 0.5, size = 13),
        strip.text = element_text(size = 11),   # 控制分面标题文字样式
        strip.background = element_rect(fill = "#D3D3D3", color = "white"),
        panel.background = element_blank(),
        plot.background = element_blank(),
        axis.title.y = element_text(size = 14),
        axis.line = element_line(color = "black",size = 0.4),
        axis.text.y = element_text(color="black",size=12),
        axis.text.x = element_text(margin = margin(t = 2),color="black",size=12,angle=0,vjust = 0,hjust = 0.5),
        legend.position = "none",
        panel.spacing = unit(0,"lines"),
        axis.ticks.x=element_line(color="black",size=0.5,lineend = 1),
        axis.line.y=element_line(linetype=1,color="black",size=0.5),
        axis.ticks.y=element_line(color="black",size=0.5,lineend = 1),
        axis.line.x=element_line(linetype=1,color="black",size=0.5))+
  coord_cartesian()


TEM_plot

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/不同分组T亚群拟时序分析和特征得分/Terminally exhausted CD4-PDCD1 T cell signatures食管鳞癌的marker.pdf",TEM_plot,width = 7,height = 5)

df_p_val7 <- df1 %>% group_by(cd4_cd8_group) %>%
  wilcox_test(Progenitor.Exhausted.CD8..T.Cell.Signature5  ~ MPR_Response2) %>%
  adjust_pvalue(p.col = "p", method = "bonferroni") %>%
  add_significance(p.col = "p.adj") %>% 
  add_xy_position(x = "MPR_Response2", dodge = 0.8)

df_p_val7_1 <- df1  %>%
  wilcox_test(Progenitor.Exhausted.CD8..T.Cell.Signature5  ~cd4_cd8_group) %>%
  adjust_pvalue(p.col = "p", method = "bonferroni") %>%
  add_significance(p.col = "p.adj") %>% 
  add_xy_position(x = "MPR_Response2", dodge = 0.8)

df1$MPR_Response2<- factor(df1$MPR_Response2,levels=c("Pre_NMPR","Pre_MPR","Post_MPR"))
Pro_plot<- df1 %>%
  ggplot(aes(MPR_Response2,Progenitor.Exhausted.CD8..T.Cell.Signature5)) +
  geom_half_boxplot(aes(fill = MPR_Response2),color="black",side="l",errorbar.draw = T,
                    outlier.shape = NA,width=0.8,lwd= 0.5) +
  geom_half_point(aes(),side = "r",
                  transformation_params = list(height = 0,width = 0.001,seed = 2))+
  stat_pvalue_manual(df_p_val7,label = "p.adj.signif",label.size=5,hide.ns = F)+
  facet_wrap(.~cd4_cd8_group,nrow=1)+
  scale_fill_manual(values= c("Pre_NMPR" = "#4974a4",
                              "Pre_MPR" = "#4dae47",
                              "Post_MPR"="#f29600"))+
  theme_minimal()+
  #scale_y_continuous(limits = c(0,95),breaks = seq(0,95,20))+
  #scale_fill_npg()+
  scale_color_npg()+
  labs(x=NULL,y=NULL)+
  labs(title = "Progenitor exhausted CD8 T cell signature",y='Signature score',x="")+
  #labs(title = "Terminally exhausted CD8 T cell signature",y='Signature score',x="")+
  theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),units=,"cm"),
        plot.title = element_text(hjust = 0.5, size = 13),
        strip.text = element_text(size = 11),   # 控制分面标题文字样式
        strip.background = element_rect(fill = "#D3D3D3", color = "white"),
        panel.background = element_blank(),
        plot.background = element_blank(),
        axis.title.y = element_text(size = 14),
        axis.line = element_line(color = "black",size = 0.4),
        axis.text.y = element_text(color="black",size=12),
        axis.text.x = element_text(margin = margin(t = 2),color="black",size=12,angle=0,vjust = 0,hjust = 0.5),
        legend.position = "none",
        panel.spacing = unit(0,"lines"),
        axis.ticks.x=element_line(color="black",size=0.5,lineend = 1),
        axis.line.y=element_line(linetype=1,color="black",size=0.5),
        axis.ticks.y=element_line(color="black",size=0.5,lineend = 1),
        axis.line.x=element_line(linetype=1,color="black",size=0.5))+
  coord_cartesian()

Pro_plot

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/不同分组T亚群拟时序分析和特征得分/Progenitor exhausted CD4-PDCD1 T cell signature食管鳞癌marker.pdf",Pro_plot,width = 7,height = 5)



####耗竭性标志物featureplot----

ppp<- FeaturePlot(object=CD8_tex_object,features =c("PDCD1","HAVCR2","ENTPD1"),ncol=1,reduction = "umap.harmony",raster=FALSE, order = T)+
  xlab("UMAP_1")+
  ylab("UMAP_1")+
  theme(legend.position = "right")
 
ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/不同分组T亚群拟时序分析和特征得分/CD8耗竭featureplot.pdf",width=5,height = 7)


##CD8细胞耗竭得分----
TEX1<- c("CXCL13","HAVCR2","PDCD1","TIGIT","LAG3","CTLA4","LAYN", "RBPJ","VCAM1","GZMB","TOX","MYO7A")
Tex2<- c("CD69","PDCD1","HAVCR2","ENTPD1","CTLA4","LAG3","TIGIT","TBX21")

m1m2_pws<- list(TEX1,Tex2)

imm_anno<-AddModuleScore(object=T_cells_object,features=m1m2_pws,name=c("Tex1","Tex2"),nbin=12)

head(imm_anno)

use_colors<-c(
  "CD4_C1-PDCD1"="#DF66B0", 
  "CD4_C2-IL7R"="#9888DB", 
  "CD4_C3-FOXP3"="#DC8062", 
  "CD4_C4-FOSB"="#D8ACC9",
  "CD8_C1-GNLY"="#83B5D3",
  "CD8_C2-CTSW"="#D5D8C5", 
  "CD8_C3-IL7R"="#77DCCF", 
  "CD8_C4-GZMK"="#BA4EE4", 
  "CD8_C5-STMN1"="#AEE94F", 
  "CD3+CD4-CD8-"="#89DB8E",
  "NK"="#E0D175")

library(randomcoloR)
Tcell_type_color = distinctColorPalette(11)
Tcell_type_color<-c("#DF66B0","#9888DB","#DC8062","#D8ACC9","#83B5D3","#D5D8C5", "#77DCCF", "#BA4EE4",
                    "#AEE94F","#89DB8E","#E0D175")

imm_anno$cd4_cd8_group<- factor(imm_anno$cd4_cd8_group,levels=c("CD4_C1-PDCD1", "CD4_C2-IL7R","CD4_C3-FOXP3","CD4_C4-FOSB","CD8_C1-GNLY",
                                                                "CD8_C2-CTSW","CD8_C3-IL7R","CD8_C4-GZMK","CD8_C5-STMN1","CD3+CD4-CD8-","NK"))
imm_anno_df <- imm_anno@meta.data
names(imm_anno_df)
imm_anno_filter<-imm_anno_df  %>% select(cd4_cd8_group,Tex11,Tex22,MPR_Response2)

library(tidyverse)
library(gghalves)
library(rstatix)
library(ggpubr)
library(ggsci)
df<-imm_anno_filter
df_p_val1 <- df %>%
  wilcox_test(Tex11~ cd4_cd8_group) %>%
  adjust_pvalue(p.col = "p", method = "bonferroni") %>%
  add_significance(p.col = "p.adj") %>% 
  add_xy_position(x = "Tex11", dodge = 0.8) 

View(df_p_val1)

treatment_color <- c("#4974a4","#4dae47","#f29600")
cell_type_color<- c("#DF66B0","#9888DB","#DC8062","#D8ACC9","#83B5D3","#D5D8C5", "#77DCCF", "#BA4EE4",
                    "#AEE94F","#89DB8E","#E0D175")

df$cd4_cd8_group<- factor(df$cd4_cd8_group,levels=c("CD4_C1-PDCD1", "CD4_C2-IL7R","CD4_C3-FOXP3","CD4_C4-FOSB","CD8_C1-GNLY",
                                                          "CD8_C2-CTSW","CD8_C3-IL7R","CD8_C4-GZMK","CD8_C5-STMN1","CD3+CD4-CD8-","NK"))
df$MPR_Response2<- factor(df$MPR_Response2,levels = c("Pre_NMPR","Pre_MPR","Post_MPR"))
Tex_sig<- df %>%
  ggplot(aes(cd4_cd8_group,Tex22)) +
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
  labs(title = "T cell exhuastion signature",y='Signature score',x="")+
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
Tex_sig

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/髓系细胞亚型/M1M2特征得分和GESA分析/M1_signature_boxplot.pdf",M1_sig,width = 6,height = 5)





####CD8STMN1到CD8GNLY转化的分子机制----
# CD8_tex_object
library(gggsea)
library(ggplot2)
library(GSVA)
library(enrichplot)
library(GSEABase)
library(fgsea)
library(clusterProfiler)
library(org.Hs.eg.db)

Idents(CD8_tex_object)<-"cd4_cd8_group"
diff_expr_genes <- FindMarkers(CD8_tex_object, ident.1 = "CD8_C1-GNLY", ident.2 = "CD8_C5-STMN1")
write.csv(diff_expr_genes,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/GNLY_STMN1_DEG.csv")
deg<-diff_expr_genes
deg$gene_name<-rownames(deg)
gmtfile <- "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/h.all.v2024.1.Hs.symbols.gmt"
hallmark <- read.gmt(gmtfile)
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
####可视化
library(ggsci)
col_gsea1<-pal_simpsons()(16)

num1=1
gseaplot2(gsea.re1,geneSetID = rownames(g1)[1:num1],
          title = "",#标题
          color = col_gsea1[1:num1],#颜色
          base_size = 14,#基础大小
          rel_heights = c(1, 0.2, 0.4),#小图相对高度
          subplots = 1:3,#展示小图
          pvalue_table = FALSE,#p值表格
          ES_geom = "line"#line or dot
)

####设置多个名称
num2=3
gseaplot2(gsea.re1,geneSetID = rownames(g1)[1:num2],
          title = "",#标题
          color = col_gsea1[1:num2],#颜色
          base_size = 14,#基础大小
          rel_heights = c(1, 0.2, 0.4),#小图相对高度
          subplots = 1:3,#展示小图
          pvalue_table = FALSE,#p值表格
          ES_geom = "line"#line or dot
)

####plotsea 可视化

plotGseaTable(hallmark.list[g2$pathway],
              id, 
              gsea.re2,gseaParam = 0.5,
              colwidths = c(0.5,0.2,0.1,0.1,0.1)
)

####ggsea 可视化
names(hallmark.list)

se_hall<-c(head(g2$pathway,9),tail(g1$pathway,9))

#3 所选中的基因集
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
            linesize=2,lty=1,
            #多个图形时使用，设定每行列图形个数
            nrow = 6,
            ncol = 2
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
ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/免疫细胞/免疫细胞GSEA/免疫细胞hallmark_GSEA.pdf",p_nes,width = 8,height=10)

######Go数据库----

library(dplyr)
library(org.Hs.eg.db) #物种注释包(Homo sapiens)
library(clusterProfiler) #富集分析
library(enrichplot)
library(ComplexHeatmap)
library(circlize)
gene <- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/GNLY_STMN1_DEG.csv",header = T)
colnames(gene)[1] <- "gene_name"
gene$regulate <-ifelse(gene$avg_log2FC>0,"up","down")
gene_ID=bitr(gene$gene_name,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db")
colnames(gene_ID)
View(gene_ID)
table(is.na(gene_ID))
gene_df <- base::merge(gene_ID,gene,by.x="SYMBOL",by.y="gene_name")

#按照LOG_FOIDCHANGE进行排序
colnames(gene_df)
gene_df_sort <- gene_df[order(gene_df$avg_log2FC, decreasing = T),]
gene_fc <- gene_df_sort$avg_log2FC
names(gene_fc) <- gene_df_sort$ENTREZID

#GSEA分析——以GO数据库为例

Go_gseresult <- gseGO(gene_fc, 'org.Hs.eg.db', keyType = "ENTREZID", ont="all", 
                      minGSSize = 10, maxGSSize = 1000, pvalueCutoff=0.05,
                      pAdjustMethod = "BH",
                      verbose = FALSE,
                      eps = 0)

Go_ges <- setReadable(Go_gseresult,
                      OrgDb = org.Hs.eg.db,
                      keyType = "ENTREZID")

Go_ges_result <- Go_ges@result
Go_ges_result$Description

gseaplot2(Go_gseresult,18,pvalue_table = TRUE)
View(Go_ges)
my_colors <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b")
####


####p8p9正式画图----
library(GseaVis)
#指定名称：图纸
selected_gene_sets <- c(5, 8, 17)
p1<-gseaNb(object = Go_ges,
           geneSetID =selected_gene_sets)

View(Go_ges)
my_colors <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b")
# clsaasic with pvalue
p8<-gseaNb(object = Go_ges,
           geneSetID = selected_gene_sets,
           curveCol=c("#1f77b4", "#ff7f0e", "#2ca02c"),
           addPval = F,
           pvalX = 0.75,pvalY = 0.8,
           pCol = 'black',
           pHjust = 0,
           subPlot = 2)
p8
selected_gene_sets=c(25, 57, 66)
p9<-gseaNb(object = Go_ges,
           geneSetID = selected_gene_sets,
           curveCol=c("#d62728", "#9467bd", "#d3a193"),
           addPval = F,
           pvalX = 0.75,pvalY = 0.8,
           pCol = 'black',
           pHjust = 0,
           subPlot = 2)
p9
##############柱状图############################3
library(stringr)
library(ggplot2)
Go_ges_result <- Go_ges@result
Go_ges_result$Description
Go_ges_result_df <- cbind(ONTOLOGY=Go_ges_result$ONTOLOGY,ID=Go_ges_result$ID,
                          Description=Go_ges_result$Description,setSize=Go_ges_result$setSize,
                          enrichmentScore=Go_ges_result$enrichmentScore,NES=Go_ges_result$NES,
                          p.adjust=Go_ges_result$pvalue,Go_ges_result$p.adjust,
                          rank=Go_ges_result$qvalue,Go_ges_result$rank) %>% as.data.frame()
Go_ges_result_df$NES <- as.numeric(as.character(Go_ges_result_df$NES))
Go_ges_result_df <- Go_ges_result_df[order(Go_ges_result_df$NES),]
Go_ges_result_df$yax <- ifelse(Go_ges_result_df$NES >0, -0.02, 0.02)
Go_ges_result_df$col <- ifelse(Go_ges_result_df$NES > 0, "blue","red")
Go_ges_result_df$NES
View(Go_ges_result_df)
p <- 
  ggplot(Go_ges_result_df,aes(y=NES,x=reorder(Description,NES),label = Description))+
  geom_text(aes(y = yax,hjust = NES > 0),size=4)+ 
  geom_bar(stat="identity",aes(fill=col),width = 0.8)+
  geom_hline(yintercept = 0)+
  coord_flip(ylim = c(-2.5,2.5))+
  labs(y = "Normalised enrichment score", x = "",title="")+
  scale_fill_manual(values = c("#FB6542","#375E97"),name="Expression",labels=c("Up","Down"))+
  scale_x_discrete(breaks = NULL)+
  theme_classic(base_size =  12)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.y = element_blank(),
        legend.position = 'bottom')
p

#####细胞过程分析###################
library(Seurat)
library(MASS)
library(harmony)
library(clustree)
library(dplyr)
library(Matrix)
library(tidyverse)
library(openxlsx)
library(scales)
T_cells_object<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/T_cells_object.rds")
T_cells_object = UpdateSeuratObject(T_cells_object)
T_cells_object
excel_data<-read.xlsx("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/不同分组T亚群拟时序分析/T_cell_gene_list.xlsx")
View(excel_data)

#文章给的基因集的excle表格
# 2.处理数据 ------------------------------------------------------------------# 使用lapply按列生成列表

marker.list <- lapply(excel_data, function(column) column)

T_cells_object<-AddModuleScore(T_cells_object, 
                           features = marker.list,
                           ctrl = 5,
                           name = "FunctionScore")

T_cells_object$FunctionScore17
for(i in 1:length(marker.list)){  
  colnames(T_cells_object@meta.data)[colnames(T_cells_object@meta.data) == paste0("FunctionScore", i)] <- names(marker.list)[i]}
Idents(T_cells_object) <- T_cells_object$cd4_cd8_group
Differentiation <- c("Naïve", "Activation:Effector.function", "Exhaustion")
Function <- c("TCR.Signaling", "Cytotoxicity", "Cytokine:Cytokine.receptor", 
              "Chemokine:Chemokine.receptor", "Senescence", "Anergy", 
              "NFKB.Signaling", "Stress.response", "MAPK.Signaling", "Adhesion",
              "IFN.Response")
Metabolism <- c("Oxidative.phosphorylation", "Glycolysis", "Fatty.acid.metabolism")
Apoptosis <- c("Pro-apoptosis", "Anti-apoptosis")
MarkerNameVector <- c(Differentiation, Function, Metabolism, Apoptosis)
FunctionScoreMatrix <- matrix(0,  ncol = length(unique(T_cells_object$cd4_cd8_group)), 
                              nrow = length(marker.list))

View(FunctionScoreMatrix)
colnames(FunctionScoreMatrix)<-c("CD4_C1-PDCD1","CD4_C2-IL7R","CD4_C3-FOXP3","CD4_C4-FOSB",
                                      "CD8_C1-GNLY","CD8_C2-CTSW","CD8_C3-IL7R","CD8_C4-GZMK","CD8_C5-STMN1",
                                      "CD3+CD4-CD8-","NK")
rownames(FunctionScoreMatrix) <- MarkerNameVector
for(ci in 1:ncol(FunctionScoreMatrix)){  
  for(ri in 1:nrow(FunctionScoreMatrix)){    
  FunctionVec <- as_tibble(T_cells_object@meta.data) %>% pull(MarkerNameVector[ri])   
  fv <- mean(FunctionVec[T_cells_object$cd4_cd8_group == levels(as.factor(T_cells_object$cd4_cd8_group))[ci]])   
  FunctionScoreMatrix[ri, ci] <- fv }}

#### 3.热图--------------------------------------------------------------------
FunctionScoreMatrix<-t(apply(FunctionScoreMatrix,1,rescale,to=c(-1,1)))
orderC = c("CD4_C1-PDCD1","CD4_C2-IL7R","CD4_C3-FOXP3","CD4_C4-FOSB",
           "CD8_C1-GNLY","CD8_C2-CTSW","CD8_C3-IL7R","CD8_C4-GZMK","CD8_C5-STMN1",
           "CD3+CD4-CD8-","NK")
FunctionScoreMatrix <- FunctionScoreMatrix[,orderC]
my.breaks <- c(seq(-1, 0, by=0.1), seq(0.1, 1, by=0.1))
my.colors <- c(colorRampPalette(colors = c("#6DCCFD", "white"))(length(my.breaks)/2),  
               colorRampPalette(colors = c("white", "#FD9AA0"))(length(my.breaks)/2))
## cellType_col <- data.frame(cell.type = CD8_Obj_CellType)
## rownames(cellType_col) <- colnames(FunctionScoreMatrix)
signatureType_row <- data.frame(Signature.type = c(  
  rep("Differentiation", length(Differentiation)),  
  rep("Function", length(Function)),  
  rep("Metabolism", length(Metabolism)),  
  rep("Apoptosis", length(Apoptosis))))
rownames(signatureType_row) <- MarkerNameVector
library(pheatmap)
pheatmap(FunctionScoreMatrix,        
         show_colnames = T,        
         show_rownames = T,        
         ## annotation_col = cellType_col,         
         annotation_row = signatureType_row,        
         gaps_row = c(3, 14, 17),        
         cluster_rows = F,         
         cluster_cols = F,        
         breaks = my.breaks,         
         color = my.colors,        
         border_color = "NA",        
         fontsize = 8,         
         width = 5,        
         height = 3.8,         
         filename = file.path("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/不同分组T亚群拟时序分析/", paste0("CD8_FunctionScore_heatmap.pdf")))

####CD8 AUCcell评分----

BiocManager::install('AUCell')

library(AUCell)
geneSets <- list(geneSet1=c("gene1", "gene2", "gene3"))

T_exhausted<-read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/不同分组T亚群拟时序分析和特征得分/T细胞耗竭毒性所用marker.csv")


T_cell_markers<-list(T_exhausted$T.Cell.Exhuastion.signature,T_exhausted$T.Cell.Mediated.Immune.Response.to.Tumor.Cell,
                     T_exhausted$Terminally.Exhausted.CD8..T.Cell.Signature,T_exhausted$Progenitor.Exhausted.CD8..T.Cell.Signature,
                     T_exhausted$Cytotoxicity,T_exhausted$Terminally.Exhausted.CD8..T.Cell.Signature_filter,
                     T_exhausted$Progenitor.Exhausted.CD8..T.Cell.Signature_filter)

####"CD8_C1-GNLY","CD8_C2-CTSW","CD8_C3-IL7R","CD8_C4-GZMK","CD8_C5-STMN1"
CD8_object = subset(T_cells_object,idents = c("CD8_C1-GNLY","CD8_C2-CTSW","CD8_C3-IL7R","CD8_C4-GZMK","CD8_C5-STMN1"))
CD8_tex_object
pbmc<-CD8_object

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
table(pbmc$RNA_snn_res.0.8)
sel.clust <- "RNA_snn_res.0.8"

pbmc <- SetIdent(pbmc,value = sel.clust)
metadata <- pbmc@meta.data

pbmc <- RunUMAP(pbmc,dims = 1:20)#UMAP二次降维，取前20个维度
pbmc <- RunTSNE(pbmc, dims = 1:20, reduction.key = "tsne_")
#harmony分析***
pbmc <- RunHarmony(object =pbmc, reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")#harmony需要pca降维的结果，而scRNA对象此前已经做过PCA降维了，所以此处不再重复
pbmc <- RunUMAP(pbmc, reduction = "harmony", dims = 1:20,reduction.name = "umap.harmony")
pbmc <- RunTSNE(object = pbmc,reduction = "harmony", reduction.key = "tsneharmony_",dims.use = 1:20,reduction.name = "tsne.harmony")
pbmc <- FindNeighbors(pbmc, reduction = "umap.harmony", dims = 1:2)
pbmc <- FindClusters(pbmc)#标准聚类

#绘图和保存
DefaultAssay(pbmc) <- "RNA"

CD8_object<-pbmc

CD8_tex_object<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/不同分组T亚群拟时序分析和特征得分/CD8_tex_object.rds")
CD8_object<-CD8_tex_object

exprMatrix<-CD8_object@assays$RNA@counts

######cells_AUC_Exhuastion.signature----

cells_AUC_Exhuastion.signature <- AUCell_run(exprMatrix, T_exhausted$T.Cell.Exhuastion.signature)

#cells_assignments<- AUCell_exploreThresholds(cells_AUC_Exhuastion.signature,plotHist = TRUE,assign=TRUE)

AUC_Exp <- as.numeric(getAUC(cells_AUC_Exhuastion.signature))

CD8_object$AUC <- AUC_Exp

plot.df <- data.frame(CD8_object@meta.data,
                      CD8_object@reductions$umap.harmony@cell.embeddings)

names(plot.df)

p <- ggplot() + 
  geom_point(data=plot.df, aes(x=umapharmony_1,y=umapharmony_2,colour=AUC), size =0.5) +
  #viridis::scale_color_viridis(option="A") +
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text= element_text(colour= 'black',size=14),
        axis.title= element_text(size = 14),
        axis.line= element_line(colour= 'black'),
        panel.border = element_rect(size = 0.5, linetype = "solid", colour = "black"), 
        aspect.ratio = 1)
p

class_avg <- data.frame(
  umapharmony_1  = c(0.70586, -6.1612144),
  umapharmony_2  = c(1.24, -1.694),
  cd4_cd8_group = c("CD8_C1-GNLY", "CD8_C5-STMN1")
)

Exhuastion_plot<- ggplot(plot.df, aes(umapharmony_1, umapharmony_2))  +
  geom_point(aes(colour  = AUC)) + viridis::scale_color_viridis(option="A") +
  ggrepel::geom_label_repel(aes(label = cd4_cd8_group),
                            data = class_avg,
                            size = 5,
                            label.size = 0,
                            segment.color = NA,
                            nudge_x = 0.1,        # 水平调整位置
                            nudge_y = 0.2)+   
  theme(legend.position = "right") + 
  theme_bw()+
  ggtitle("T Cell Exhuastion signature (AUCell)") +
  theme(plot.title = element_text(hjust = 0.5))+
  xlab("UMAP_1") +
  ylab("UMAP_2") +
  theme(
    panel.background = element_rect(fill = "white", color = NA),  # 设置背景为白色
    panel.grid.major = element_blank(),  # 移除主网格线
    panel.grid.minor = element_blank(),  # 移除次网格线
    axis.title.x = element_text(size = 12),  # X 轴名称字体大小和样式
    axis.title.y = element_text(size = 12),
    panel.border = element_rect(color = "black", fill = NA, size = 0.8),
    axis.text.x = element_text(size = 11),  # 设置 X 轴刻度字体大小
    axis.text.y = element_text(size = 11) 
  )

Exhuastion_plot
ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/不同分组T亚群拟时序分析和特征得分/CD8_cells_exhausted_UMAP_plot.pdf",Exhuastion_plot,width=7,height=6)


######
saveRDS(CD8_tex_object,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/不同分组T亚群拟时序分析和特征得分/CD8_tex_object.rds")
CD8_tex_object = subset(T_cells_object,idents = c("CD8_C1-GNLY","CD8_C5-STMN1"))

pbmc<-CD8_tex_object

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
table(pbmc$RNA_snn_res.0.8)
sel.clust <- "RNA_snn_res.0.8"

pbmc <- SetIdent(pbmc,value = sel.clust)
metadata <- pbmc@meta.data

pbmc <- RunUMAP(pbmc,dims = 1:20)#UMAP二次降维，取前20个维度
pbmc <- RunTSNE(pbmc, dims = 1:20, reduction.key = "tsne_")
#harmony分析***
pbmc <- RunHarmony(object =pbmc, reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")#harmony需要pca降维的结果，而scRNA对象此前已经做过PCA降维了，所以此处不再重复
pbmc <- RunUMAP(pbmc, reduction = "harmony", dims = 1:20,reduction.name = "umap.harmony")
pbmc <- RunTSNE(object = pbmc,reduction = "harmony", reduction.key = "tsneharmony_",dims.use = 1:20,reduction.name = "tsne.harmony")
pbmc <- FindNeighbors(pbmc, reduction = "umap.harmony", dims = 1:2)
pbmc <- FindClusters(pbmc)#标准聚类

#绘图和保存
DefaultAssay(pbmc) <- "RNA"

CD8_tex_object<-pbmc


exprMatrix<-CD8_tex_object@assays$RNA@counts

cells_AUC_Terminally.Exhuastion.signature <- AUCell_run(exprMatrix, T_exhausted$Terminally.Exhausted.CD8..T.Cell.Signature_filter)

#cells_assignments<- AUCell_exploreThresholds(cells_AUC_Exhuastion.signature,plotHist = TRUE,assign=TRUE)

AUC_Exp <- as.numeric(getAUC(cells_AUC_Terminally.Exhuastion.signature))

CD8_tex_object$AUC.Terminally <- AUC_Exp
CD8_tex_object$cd4_cd8_group
plot.df <- data.frame(CD8_tex_object@meta.data,
                      CD8_tex_object@reductions$umap.harmony@cell.embeddings)

names(plot.df)

p <- ggplot() + 
  geom_point(data=plot.df, aes(x=umapharmony_1,y=umapharmony_2,colour=AUC.Terminally), size =0.5) +
  #viridis::scale_color_viridis(option="A") +
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text= element_text(colour= 'black',size=14),
        axis.title= element_text(size = 14),
        axis.line= element_line(colour= 'black'),
        panel.border = element_rect(size = 0.5, linetype = "solid", colour = "black"), 
        aspect.ratio = 1)
p

class_avg <- data.frame(
  umapharmony_1  = c(0.70586, -6.1612144),
  umapharmony_2  = c(1.24, -1.694),
  cd4_cd8_group = c("CD8_C1-GNLY", "CD8_C5-STMN1")
)
Terminally.Exhuastion_plot<- ggplot(plot.df, aes(umapharmony_1, umapharmony_2))  +
  geom_point(aes(colour  = AUC.Terminally)) + viridis::scale_color_viridis(option="A") +
  ggrepel::geom_label_repel(aes(label = cd4_cd8_group),
                            data = class_avg,
                            size = 5,
                            label.size = 0,
                            segment.color = NA,
                            nudge_x = 0.1,        # 水平调整位置
                            nudge_y = 0.2)+   
  theme(legend.position = "right") + 
  theme_bw()+
  ggtitle("Terminally exhausted CD8+ T cell signature (AUCell)") +
  theme(plot.title = element_text(hjust = 0.5))+
  xlab("UMAP_1") +
  ylab("UMAP_2") +
  theme(
    panel.background = element_rect(fill = "white", color = NA),  # 设置背景为白色
    panel.grid.major = element_blank(),  # 移除主网格线
    panel.grid.minor = element_blank(),  # 移除次网格线
    axis.title.x = element_text(size = 12),  # X 轴名称字体大小和样式
    axis.title.y = element_text(size = 12),
    panel.border = element_rect(color = "black", fill = NA, size = 0.8),
    axis.text.x = element_text(size = 11),  # 设置 X 轴刻度字体大小
    axis.text.y = element_text(size = 11) 
  )

Terminally.Exhuastion_plot
ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/不同分组T亚群拟时序分析和特征得分/CD8_cells_Terminally.exhausted_UMAP_plot.pdf",Terminally.Exhuastion_plot,width=7,height=6)


exprMatrix<-CD8_tex_object@assays$RNA@counts

cells_AUC_Progenitor.Exhuastion.signature <- AUCell_run(exprMatrix, T_exhausted$Progenitor.Exhausted.CD8..T.Cell.Signature_filter)

#cells_assignments<- AUCell_exploreThresholds(cells_AUC_Exhuastion.signature,plotHist = TRUE,assign=TRUE)

AUC_Exp <- as.numeric(getAUC(cells_AUC_Progenitor.Exhuastion.signature))

CD8_tex_object$AUC.Progenitor <- AUC_Exp

plot.df <- data.frame(CD8_tex_object@meta.data,
                      CD8_tex_object@reductions$umap.harmony@cell.embeddings)

names(plot.df)

p <- ggplot() + 
  geom_point(data=plot.df, aes(x=umapharmony_1,y=umapharmony_2,colour=AUC.Progenitor), size =0.5) +
  #viridis::scale_color_viridis(option="A") +
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text= element_text(colour= 'black',size=14),
        axis.title= element_text(size = 14),
        axis.line= element_line(colour= 'black'),
        panel.border = element_rect(size = 0.5, linetype = "solid", colour = "black"), 
        aspect.ratio = 1)
p

class_avg <- data.frame(
  umapharmony_1  = c(0.70586, -6.1612144),
  umapharmony_2  = c(1.24, -1.694),
  cd4_cd8_group = c("CD8_C1-GNLY", "CD8_C5-STMN1")
)

Progenitor.Exhuastion_plot<- ggplot(plot.df, aes(umapharmony_1, umapharmony_2))  +
  geom_point(aes(colour  = AUC.Progenitor)) + viridis::scale_color_viridis(option="A") +
  ggrepel::geom_label_repel(aes(label = cd4_cd8_group),
                            data = class_avg,
                            size = 5,
                            label.size = 0,
                            segment.color = NA,
                            nudge_x = 0.1,        # 水平调整位置
                            nudge_y = 0.2)+   
  theme(legend.position = "right") + 
  theme_bw()+
  ggtitle("Progenitor exhausted CD8+ T cell signature(AUCell)") +
  theme(plot.title = element_text(hjust = 0.5))+
  xlab("UMAP_1") +
  ylab("UMAP_2") +
  theme(
    panel.background = element_rect(fill = "white", color = NA),  # 设置背景为白色
    panel.grid.major = element_blank(),  # 移除主网格线
    panel.grid.minor = element_blank(),  # 移除次网格线
    axis.title.x = element_text(size = 12),  # X 轴名称字体大小和样式
    axis.title.y = element_text(size = 12),
    panel.border = element_rect(color = "black", fill = NA, size = 0.8),
    axis.text.x = element_text(size = 11),  # 设置 X 轴刻度字体大小
    axis.text.y = element_text(size = 11) 
  )

Progenitor.Exhuastion_plot
ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/不同分组T亚群拟时序分析和特征得分/CD8_cells_Progenitor.exhausted_UMAP_plot.pdf",Progenitor.Exhuastion_plot,width=7,height=6)


####CD4 AUCcell评分----

BiocManager::install('AUCell')

library(AUCell)
geneSets <- list(geneSet1=c("gene1", "gene2", "gene3"))

T_exhausted<-read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/不同分组T亚群拟时序分析和特征得分/T细胞耗竭毒性所用marker.csv")


T_cell_markers<-list(T_exhausted$T.Cell.Exhuastion.signature,T_exhausted$T.Cell.Mediated.Immune.Response.to.Tumor.Cell,
                     T_exhausted$Terminally.Exhausted.CD8..T.Cell.Signature,T_exhausted$Progenitor.Exhausted.CD8..T.Cell.Signature,
                     T_exhausted$Cytotoxicity)

####"CD8_C1-GNLY","CD8_C2-CTSW","CD8_C3-IL7R","CD8_C4-GZMK","CD8_C5-STMN1"
CD4_object = subset(T_cells_object,idents = c("CD4_C1-PDCD1","CD4_C2-IL7R","CD4_C3-FOXP3","CD4_C4-FOSB"))

pbmc<-CD4_object

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
table(pbmc$RNA_snn_res.0.8)
sel.clust <- "RNA_snn_res.0.8"

pbmc <- SetIdent(pbmc,value = sel.clust)
metadata <- pbmc@meta.data

pbmc <- RunUMAP(pbmc,dims = 1:20)#UMAP二次降维，取前20个维度
pbmc <- RunTSNE(pbmc, dims = 1:20, reduction.key = "tsne_")
#harmony分析***
pbmc <- RunHarmony(object =pbmc, reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")#harmony需要pca降维的结果，而scRNA对象此前已经做过PCA降维了，所以此处不再重复
pbmc <- RunUMAP(pbmc, reduction = "harmony", dims = 1:20,reduction.name = "umap.harmony")
pbmc <- RunTSNE(object = pbmc,reduction = "harmony", reduction.key = "tsneharmony_",dims.use = 1:20,reduction.name = "tsne.harmony")
pbmc <- FindNeighbors(pbmc, reduction = "umap.harmony", dims = 1:2)
pbmc <- FindClusters(pbmc)#标准聚类

#绘图和保存
DefaultAssay(pbmc) <- "RNA"

CD4_object<-pbmc


exprMatrix<-CD4_object@assays$RNA@counts


cells_AUC_Exhuastion.signature <- AUCell_run(exprMatrix, T_exhausted$T.Cell.Exhuastion.signature)

#cells_assignments<- AUCell_exploreThresholds(cells_AUC_Exhuastion.signature,plotHist = TRUE,assign=TRUE)

AUC_Exp <- as.numeric(getAUC(cells_AUC_Exhuastion.signature))

CD4_object$AUC <- AUC_Exp

plot.df <- data.frame(CD4_object@meta.data,
                      CD4_object@reductions$umap.harmony@cell.embeddings)

names(plot.df)

p <- ggplot() + 
  geom_point(data=plot.df, aes(x=umapharmony_1,y=umapharmony_2,colour=AUC), size =0.5) +
  #viridis::scale_color_viridis(option="A") +
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text= element_text(colour= 'black',size=14),
        axis.title= element_text(size = 14),
        axis.line= element_line(colour= 'black'),
        panel.border = element_rect(size = 0.5, linetype = "solid", colour = "black"), 
        aspect.ratio = 1)
p

class_avg <- data.frame(
  umapharmony_1  = c(0.70586, -6.1612144),
  umapharmony_2  = c(1.24, -1.694),
  cd4_cd8_group = c("CD8_C1-GNLY", "CD8_C5-STMN1")
)

Exhuastion_plot<- ggplot(plot.df, aes(umapharmony_1, umapharmony_2))  +
  geom_point(aes(colour  = AUC)) + viridis::scale_color_viridis(option="A") +
  ggrepel::geom_label_repel(aes(label = cd4_cd8_group),
                            data = class_avg,
                            size = 5,
                            label.size = 0,
                            segment.color = NA,
                            nudge_x = 0.1,        # 水平调整位置
                            nudge_y = 0.2)+   
  theme(legend.position = "right") + 
  theme_bw()+
  ggtitle("CD4+ cell Exhuastion signature (AUCell)") +
  theme(plot.title = element_text(hjust = 0.5))+
  xlab("UMAP_1") +
  ylab("UMAP_2") +
  theme(
    panel.background = element_rect(fill = "white", color = NA),  # 设置背景为白色
    panel.grid.major = element_blank(),  # 移除主网格线
    panel.grid.minor = element_blank(),  # 移除次网格线
    axis.title.x = element_text(size = 12),  # X 轴名称字体大小和样式
    axis.title.y = element_text(size = 12),
    panel.border = element_rect(color = "black", fill = NA, size = 0.8),
    axis.text.x = element_text(size = 11),  # 设置 X 轴刻度字体大小
    axis.text.y = element_text(size = 11) 
  )

Exhuastion_plot
ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/不同分组T亚群拟时序分析和特征得分/CD4_cells_exhausted_UMAP_plot.pdf",Exhuastion_plot,width=7,height=6)



####拟时序分析----

#####直接写一个功能函数----
ks_run_Monocle2 <- function(object, #seurat obj or expression matrix (建议数据格式转为matrix,如果数据量大转化为稀疏矩阵as(as.matrix(data), "sparseMatrix"))
                            layer, #used when object is a seurat obj
                            assay, #used when object is a seurat obj
                            lowerDetectionLimit = 0.1, 
                            VARgenesM=c("dispersionTable","seurat","differentialGeneTest"),
                            cellAnno=NULL, 
                            define_root=F,
                            root_state,
                            reverse=NULL
){
  
  
  if(class(object)[1] == 'Seurat') {
    
    data <- GetAssayData(object=object, layer=layer, assay=assay)#get expression matrix data
    
    pd <- new("AnnotatedDataFrame", data = object@meta.data)
    fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
    fd <- new("AnnotatedDataFrame", data = fData)
    
    if(all(data == floor(data))) {
      expressionFamily <- negbinomial.size()
    } else if(any(data < 0)){
      expressionFamily <- uninormal()
    } else {
      expressionFamily <- tobit()
    }
    
    #Creates a new CellDateSet object.
    monocle_cds <- newCellDataSet(data,
                                  phenoData = pd, 
                                  featureData = fd,
                                  lowerDetectionLimit=0.1,
                                  expressionFamily=expressionFamily)
    
  }else{
    
    print("This fucntions only apply for a seurat obj")
  }
  
  
  
  #Estimate size factors and dispersions
  #数据处理
  monocle_cds <- estimateSizeFactors(monocle_cds)#size facotr标准化细胞之间的mRNA差异
  monocle_cds <- estimateDispersions(monocle_cds)
  
  #质量控制-filter cells
  monocle_cds <- detectGenes(monocle_cds, min_expr = 0.1)
  # print(head(fData(monocle_cds)))
  # print(head(pData(monocle_cds)))
  # expressed_genes <- row.names(subset(fData(mouse_monocle), num_cells_expressed >= 10))
  monocle_cds <- monocle_cds[fData(monocle_cds)$num_cells_expressed >= 10, ]
  
  
  #select methods for VariableFeatures
  if(VARgenesM=="dispersionTable"){
    
    disp_table <- dispersionTable(monocle_cds)
    ordering_genes <- subset(disp_table,
                             mean_expression >= 0.1 &
                               dispersion_empirical >= 1.5* dispersion_fit)$gene_id
    
  }
  
  
  if(VARgenesM=="seurat"){
    
    ordering_genes <- VariableFeatures(FindVariableFeatures(object, assay = "RNA"), assay = "RNA")
    
  }
  
  
  if(VARgenesM=="differentialGeneTest"){
    
    diff_test_res <- differentialGeneTest(monocle_cds,fullModelFormulaStr = paste0("~",cellAnno))##~后面是表示对谁做差异分析的变量
    diff_test_res_sig <- diff_test_res[order(diff_test_res$qval,decreasing=F),]
    
    ordering_sce <- diff_test_res_sig[diff_test_res_sig$qval< 0.01,]
    
    if(nrow(ordering_sce)>3000){
      
      ordering_genes <- ordering_sce$gene_short_name[1:3000]
      
    }else{
      
      ordering_genes <- rdering_sce$gene_short_name
    }
    
  }
  
  #Marks genes for clustering
  monocle_cds <- setOrderingFilter(monocle_cds, ordering_genes)
  plot_ordering_genes(monocle_cds)
  
  
  #cluster
  monocle_cds <- reduceDimension(monocle_cds, max_components = 2,reduction_method = 'DDRTree')
  
  #order cells
  monocle_cds <- orderCells(monocle_cds, reverse=reverse)
  
  if(define_root){
    monocle_cds <- monocle_cds <- orderCells(monocle_cds,root_state = root_state)
  }
  
  
  return(monocle_cds)
  
}

####CD8 耗竭拟时序----
library(monocle)
library(Seurat)
library(tidyverse)
library(harmony)
library(dplyr)
library(patchwork)
library(SeuratWrappers)
# T_cells_object <- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/T_cells_object.rds")
# CD8_tex_object = subset(T_cells_object,idents = c("CD8_C1-GNLY","CD8_C5-STMN1"))
CD8_tex_object<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/不同分组T亚群拟时序分析和特征得分/CD8_tex_object.rds")
CD8_tex_object$cd4_cd8_group
sce_CDS <- ks_run_Monocle2(object =CD8_object,
                           layer = 'counts',
                           assay = "RNA",
                           VARgenesM="seurat",
                           cellAnno = "cd4_cd8_group")

sce_CDS$cd4_cd8_group
saveRDS(sce_CDS,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/不同分组T亚群拟时序分析和特征得分/CD8IL7R拟时序.rds")

plot_cell_trajectory(sce_CDS,color_by = "Pseudotime")
plot_cell_trajectory(sce_CDS,color_by = "State")

plot_cell_trajectory(sce_CDS,color_by = "cd4_cd8_group")


sce_CDS

####美化cd8_plot_cell_trajectory----
colour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",  
         "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0",
         "#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE",
         "#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")

colour=c("#83B5D3","#AEE94F")
p1 <- plot_cell_trajectory(sce_CDS, x = 1, y = 2, cell_size=2,cell_link_size = 1,
                           color_by = "cd4_cd8_group") + 
  theme(legend.position='none',panel.border = element_blank()) + #去掉第一个的legend
  scale_color_manual(values = colour) +
  theme(axis.text = element_text(color = 'black', size = 12),        
        axis.title = element_text(color = 'black', size = 14),       
        strip.text = element_text(color = 'black', size = 14),
        axis.line.x.bottom = element_line(linewidth = 0.6, color = "black"),
        axis.line.y.left = element_line(linewidth = 0.6, color = "black"),
        axis.line = element_line(linewidth = 2, color = "black")) 
p2 <- plot_complex_cell_trajectory(sce_CDS, x = 1, y = 2,cell_size=2,cell_link_size = 1,
                                   color_by = "cd4_cd8_group")+
  scale_color_manual(values = colour) +
  theme(legend.title = element_blank()) +
  theme(axis.text = element_text(color = 'black', size = 12),        
        axis.title = element_text(color = 'black', size = 14),       
        strip.text = element_text(color = 'black', size = 14),
        axis.line.x.bottom = element_line(linewidth = 0.6, color = "black"),
        axis.line.y.left = element_line(linewidth = 0.6, color = "black"),
        axis.line = element_line(linewidth = 2, color = "black"),
        legend.title = element_text(size = 14),   # 图例标题大小
        legend.text = element_text(size = 12),    # 图例文本大小
        legend.position = "bottom",               # 图例位置，例如底部
        legend.key.size = unit(1, "cm")) +
  guides(color = guide_legend(override.aes = list(size = 3)))

p3<- p1 / p2

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/不同分组T亚群拟时序分析和特征得分/monoclecd8拟时序加树.pdf",p3,width = 8,height=9)

####monocle 分支----
plot_cell_trajectory(sce_CDS, color_by = "Pseudotime", show_branch_points = TRUE)
BEAM_res <- BEAM(sce_CDS, branch_point = 1, cores = 1, progenitor_method = "duplicate")
BEAM_res <- BEAM(sce_CDS, branch_point = 1, cores = 1)
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
#只保留"gene_short_name", "pval", "qval"这3列；
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
head(BEAM_res)

disp_table <- dispersionTable(sce_CDS)
disp.genes <- subset(disp_table, mean_expression >= 0.5&dispersion_empirical >= 1*dispersion_fit)
disp.genes <- as.character(disp.genes$gene_id)
mycds_sub <- sce_CDS[disp.genes,]
plot_cell_trajectory(mycds_sub, color_by = "State")
beam_res <- BEAM(mycds_sub, branch_point = 1, cores = 8,progenitor_method == "duplicate")
beam_res <- beam_res[order(beam_res$qval),]
beam_res <- beam_res[,c("gene_short_name", "pval", "qval")]
mycds_sub_beam <- mycds_sub[row.names(subset(beam_res, qval < 1e-4)),]
plot_genes_branched_heatmap(mycds_sub_beam,  branch_point = 1, num_clusters = 3, show_rownames = T)

lung_genes <- row.names(subset(fData(sce_CDS),gene_short_name %in% c("PDCD1")))
#绘制表达趋势散点图；
plot_genes_branched_pseudotime(sce_CDS[lung_genes,],
                               branch_point = 1,
                               color_by = "Time",
                               ncol = 1)

####monocle差异基因----

library(ggplot2)
library(ggsci) 
# pData(sce_CDS)$PDCD1 = log2(exprs(sce_CDS)['PDCD1',]+1)
# plot_cell_trajectory(sce_CDS, color_by = 'PDCD1') + 
#   scale_color_gsea()

cg1 <- c("CD69",
        "PDCD1",
        "HAVCR2",
        "ENTPD1",
        "LAG3",
        "CTLA4",
        "GNLY",
        "MKI67",
        "STMN1") 


p1 <- plot_genes_in_pseudotime(sce_CDS[cg1, ], 
                              cell_size = 1, 
                              ncol = 3, 
                              color_by = "cd4_cd8_group")+
  scale_color_manual(name=NULL, values = c("#83B5D3","#AEE94F")) +
  theme(axis.text = element_text(color = 'black', size = 12),        
        axis.title = element_text(color = 'black', size = 14),       
        strip.text = element_text(color = 'black', size = 14),
        axis.line.x.bottom = element_line(linewidth = 0.6, color = "black"),
        axis.line.y.left = element_line(linewidth = 0.6, color = "black"),
        axis.line = element_line(linewidth = 2, color = "black"),
       legend.title = element_text(size = 14),   # 图例标题大小
              legend.text = element_text(size = 12),    # 图例文本大小
              legend.position = "right",               # 图例位置，例如底部
              legend.key.size = unit(1, "cm")) +
  guides(color = guide_legend(override.aes = list(size = 2)))

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/不同分组T亚群拟时序分析和特征得分/monocle基因伪时间变化曲线.pdf",p1,width=10,height = 6)

####monocle 基因变化热图----
library(monocle)
sce_CDS <- orderCells(sce_CDS)
names(pData(sce_CDS))
diff_test_res <- differentialGeneTest(sce_CDS, fullModelFormulaStr = "~ State")
sig_genes <- diff_test_res[diff_test_res$qval < 0.05, ]



T_cells_object_deg_all=FindMarkers(T_cells_object, ident.1 = c("Bacteria+"), ident.2 = "Bacteria-")


# 提取感兴趣的基因表达数据---
library(dplyr)
library(monocle)
library(ggplot2)
library(Seurat)
colnames(pData(sce_CDS))
genes <-  c("CD69",
            "PDCD1",
            "HAVCR2",
            "ENTPD1",
            "LAG3",
            "CTLA4",
            "GNLY",
            "MKI67",
            "STMN1") 

genes_exp<-list()

for(i in 1:length(genes)){  
  A <- log2(exprs(sce_CDS)[genes[i],]+1)  
  A <- as.data.frame(A)  
  genes_exp[[i]] <- A}

gene_exp <- do.call(cbind, genes_exp)
colnames(gene_exp) <- genes#将上述几个基因的拟时表达添加到
pData(sce_CDS) = cbind(pData(sce_CDS), gene_exp)
#names(pData(sce_CDS))
#pData(sce_CDS) <-pData(sce_CDS)[,-c(93:95)]

data <- pData(sce_CDS)
colnames(data)
library(reshape2)
data<-data %>% select("orig.ident","Pseudotime","MPR_Response2","CD69",
                      "PDCD1",
                      "HAVCR2",
                      "ENTPD1",
                      "LAG3",
                      "CTLA4",
                      "GNLY",
                      "MKI67",
                      "STMN1")

data_long_m<-melt(data, id.vars = c("orig.ident","Pseudotime","MPR_Response2"), #需保留的不参与聚合的变量列名                 
                  measure.vars = 4:12,#选择需要转化的列                 
                  variable.name = c('gene'),#聚合变量的新列名                 
                  value.name = 'value')#聚合值的新列名
  colnames(data_long_m)

data_long_m$MPR_Response2<- factor(data_long_m$MPR_Response2,levels = c("Pre_NMPR","Pre_MPR","Post_MPR"))          
p_group<-ggplot(data_long_m, aes(x=Pseudotime, y=value,color=MPR_Response2))+
  geom_smooth(aes(fill=MPR_Response2))+ #平滑的填充  
  xlab('Pseudotime') +   
  ylab('Relative expression') +  
  facet_wrap(~gene, scales = "free_y")+ #分面，y轴用各自数据 
  theme(axis.text = element_text(color = 'black',size = 12),        
        axis.title = element_text(color = 'black',size = 14),       
        strip.text = element_text(color = 'black',size = 14),
        panel.background = element_blank(),  # 去除面板背景
        panel.grid.major = element_blank(),   # 去除主网格线
        panel.grid.minor = element_blank())+
   scale_color_manual(name=NULL, values =  c("#4974a4","#4dae47","#f29600")) +
   scale_fill_manual(name=NULL, values =  c("#4974a4","#4dae47","#f29600"))

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/不同分组T亚群拟时序分析和特征得分/monocle基因伪时间变化曲线分组.pdf",p_group,width=10,height = 6)

####时间轴密度----

library(ggpubr) 
df <- pData(sce_CDS)  
# pData(cds)取出的是cds对象中cds@phenoData@data的内容
names(df)
density_plot <- ggplot(df, 
       aes(Pseudotime, 
           colour = cd4_cd8_group, 
           fill = cd4_cd8_group)) +   
  geom_density( # 绘制密度图
    bw = 0.5,  # 设定带宽，影响密度曲线的平滑程度
    size = 0.8, # 密度曲线的粗细
    alpha = 0.5) +  # 填充颜色的透明度
  theme_classic2()+
scale_color_manual(name=NULL, values = c("#83B5D3","#AEE94F"))+
  scale_fill_manual(name=NULL, values = c("#83B5D3","#AEE94F"))

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/不同分组T亚群拟时序分析和特征得分/monocle_密度图.pdf",density_plot,width = 7,height = 6)



####tradeSeq拟时序----
if(!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("tradeSeq")
library(Seurat)
library(tradeSeq)
library(RColorBrewer)
library(SingleCellExperiment)
library(slingshot)
library(tidyverse)
library(viridis)
library(SingleCellExperiment)
CD8_object = subset(T_cells_object,idents = c("CD8_C1-GNLY", "CD8_C2-CTSW", "CD8_C3-IL7R", "CD8_C4-GZMK","CD8_C5-STMN1"))
CD8_object = subset(T_cells_object,idents = c("CD8_C1-GNLY", "CD8_C5-STMN1"))
CD8_object = subset(T_cells_object,idents =c("CD8_C3-IL7R", "CD8_C4-GZMK"))
CD8_tex_object<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/不同分组T亚群拟时序分析和特征得分/CD8_tex_object.rds")
#CD8_tex_object<- readRDS("/home/zqwangyansu/oscc_data/CD8_tex_object.rds")

pbmc<-CD8_object

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
table(pbmc$RNA_snn_res.0.8)
sel.clust <- "RNA_snn_res.0.8"

pbmc <- SetIdent(pbmc,value = sel.clust)
metadata <- pbmc@meta.data

pbmc <- RunUMAP(pbmc,dims = 1:20)#UMAP二次降维，取前20个维度
pbmc <- RunTSNE(pbmc, dims = 1:20, reduction.key = "tsne_")
#harmony分析***
library(harmony)
pbmc <- RunHarmony(object =pbmc, reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")#harmony需要pca降维的结果，而scRNA对象此前已经做过PCA降维了，所以此处不再重复
pbmc <- RunUMAP(pbmc, reduction = "harmony", dims = 1:20,reduction.name = "umap.harmony")
pbmc <- RunTSNE(object = pbmc,reduction = "harmony", reduction.key = "tsneharmony_",dims.use = 1:20,reduction.name = "tsne.harmony")
pbmc <- FindNeighbors(pbmc, reduction = "umap.harmony", dims = 1:2)
pbmc <- FindClusters(pbmc)#标准聚类

#绘图和保存
DefaultAssay(pbmc) <- "RNA"

CD8_object<-pbmc
saveRDS(CD8_object,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/不同分组T亚群拟时序分析和特征得分/CD8_IL7R_object.rds")
unique(CD8_object$cd4_cd8_group)

#CD8_object<-CD8_tex_object
CD8_object@reductions
# Seurat 转化为 SingleCellExperiment
CD8_object_sce <-
  as.SingleCellExperiment(x = CD8_object, assay = "RNA")
reducedDims(CD8_object_sce)
# Run slingshot
CD8_object_sce <-
  slingshot(
    CD8_object_sce, # 选择的SingleCellExperiment对象
    reducedDim = "UMAP.HARMONY",  # 设置使用的降维
    clusterLabels = CD8_object_sce$cd4_cd8_group, # cluster标签
    start.clus = "CD8_C5-STMN1", # 谱系起点cluster
    # end.clus = "PS7D", # 谱系终点cluster
    approx_points = 150 # 近似点 默认设置
  )
CD8_object_sce <-
  slingshot(
    CD8_object_sce, # 选择的SingleCellExperiment对象
    reducedDim = "UMAP.HARMONY",  # 设置使用的降维
    clusterLabels = CD8_object_sce$cd4_cd8_group, # cluster标签
    # end.clus = "PS7D", # 谱系终点cluster
    approx_points = 150 # 近似点 默认设置
  )

# Add to Seurat object 将slingshot信息添加至Seurat
CD8_object_sce@colData %>%
  colnames() %>%
  grep(pattern = "^sling", x = ., value = TRUE) -> slingshot_res_colnames
slingshot_res_colnames
for (cn in slingshot_res_colnames) {
  message(cn)
  CD8_object[[cn]] <- CD8_object_sce@colData[, cn]
}

# plot
Tcell_type_color<-c("#77DCCF", "#BA4EE4")

Tcell_type_color<-c("#83B5D3","#AEE94F")

CD8_object_sce.groupByGroup.umap <- DimPlot(
  CD8_object,
  reduction = "umap.harmony",
  group.by = "cd4_cd8_group",
  label = FALSE,
  shuffle = TRUE,
  pt.size = 1.5,
  cols = c("#83B5D3","#AEE94F") # 从ggsci包中得到 
) +
  ggtitle("UMAP") +
  labs(x = "UMAP_1", y = "UMAP_2") +
  theme(title = element_text(size = 10))+
  labs(x = "UMAP 1", y = "UMAP 2") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    legend.text = element_text(size = 11),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    plot.title = element_text(size =15, hjust = 0.5),
    axis.text.x = element_text(size = 11),
    #axis.line = element_line(linewidth = 1.2, color = "black"),  # 调整轴线的粗细和颜色
    panel.border = element_rect(linewidth = 1, color = "black") # 调整四周边框的粗细和颜色
  )+
  guides(
    color = guide_legend(override.aes = list(size = 6))    # 调整图例符号的大小
  )

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/不同分组T亚群拟时序分析和特征得分/slingshot & tradeSeq拟时序分析.pdf",CD8_object_sce.groupByGroup.umap,width = 7,height = 5)

plot.df <- data.frame(CD8_object@meta.data,
                      CD8_object@reductions$umap.harmony@cell.embeddings)

names(plot.df)

plot.df_clean <- plot.df %>%
  select(-slingshot)  # 去掉 slingshot 列
plot.df_clean$cd4_cd8_group

library(dplyr)
class_avg <- plot.df_clean %>%
  group_by(cd4_cd8_group) %>%
  summarise(
    umapharmony_1 = median(umapharmony_1),
    umapharmony_2 = median(umapharmony_2)
  )



ffplot<- FeaturePlot(
  CD8_object,
  reduction = "umap.harmony",
  features = "slingPseudotime_1",
  pt.size = 1.5,
  cols = viridis_pal(begin = 0.1, end = 0.9)(2)
)+
  #lims(x = c(-0.25, 0.25), y = c(-0.25, 0.25))+
  labs(x = "UMAP_1", y = "UMAP_2") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    legend.text = element_text(size = 11),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    plot.title = element_text(size =15, hjust = 0.5),
    axis.text.x = element_text(size = 11),
    #axis.line = element_line(linewidth = 1.2, color = "black"),  # 调整轴线的粗细和颜色
    panel.border = element_rect(linewidth = 1, color = "black") # 调整四周边框的粗细和颜色
  )+ ggtitle("Pseudotime")

class_avg <- data.frame(
  umapharmony_1  = c(2.41, 3.05),
  umapharmony_2  = c(-0.345, -0.314),
  cd4_cd8_group = c("CD8_C3-IL7R", "CD8_C4-GZMK")
)

ffplot1<-ffplot+ 
  ggrepel::geom_label_repel(data = class_avg, 
                            aes(x = umapharmony_1, y = umapharmony_2, label = cd4_cd8_group),
                            size = 5, 
                            box.padding = 0.5, 
                            point.padding = 0.5, 
                            segment.color = 'grey50')+theme(
                              legend.position = "top",        # 图例放置在顶部
                              legend.direction = "horizontal" # 图例横向显示
                            )

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/不同分组T亚群拟时序分析和特征得分/slingshot & tradeSeq拟时序分析_Pseudotime.pdf",ffplot1,width =6,height = 6)

pdf("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/不同分组T亚群拟时序分析和特征得分/slingshot & tradeSeq拟时序分析_Pseudotime_重新画图.pdf", width = 6, height = 6) 
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(CD8_object_sce$slingPseudotime_1, breaks=100)]
plot(reducedDims(CD8_object_sce)$UMAP.HARMONY, col = plotcol, pch=16, asp = 1, cex = 0.8,xlab = "UMAP_1",ylab = "UMAP_2")
box(lwd = 2, col = "black")
lines(SlingshotDataSet(CD8_object_sce), lwd=2, col='black')
dev.off()

pdf("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/不同分组T亚群拟时序分析和特征得分/slingshot & tradeSeq拟时序分析_Pseudotime_重新画图legend.pdf", width = 3, height = 3) 
image(1:100, 1, as.matrix(1:100), col = colors, axes = FALSE, xlab = "", ylab = "")
axis(1, at = seq(1, 100, by = 25), labels = seq(0, 15, length.out = 4))
title(xlab = "Pseudotime", line = 2)
dev.off()

####trade 基因绘图----
pseudotime <- slingPseudotime(CD8_object_sce, na = FALSE)
cellWeights <- slingCurveWeights(CD8_object_sce)
crv <- SlingshotDataSet(CD8_object_sce)
counts<-CD8_object_sce@assays@data@listData$counts
length(rownames(counts))
library(BiocParallel)

set.seed(111)
icMat <- evaluateK(counts = counts, 
                   sds = crv, 
                   k = 3:10,    # no more than 12
                   nGenes = 200, # 每个细胞纳入分析的基因数量，默认是500，这里为了节省示例计算时间
                   verbose = T)


system.time(sce <- fitGAM(counts = counts, 
              pseudotime = pseudotime, 
              cellWeights = cellWeights,
              nknots = 6, verbose = TRUE,
              BPPARAM = MulticoreParam(20),parallel=T))

table(rowData(sce)$tradeSeq$converged)
# pdf("plotSmoothers.pdf", width = 8, height = 6)
# plotSmoothers(sce, counts, gene = "PDCD1")
# dev.off()
# 
# pdf("plotGeneCount.pdf", width = 8, height = 6)
# plotGeneCount(crv, counts, gene = sigGeneStart)
# dev.off()
# 
# plotSmoothers(gam, gene = "HAVCR2", theme_size = 12)
####探索基因表达与拟时序的相关性
assoRes <- associationTest(sce)
#write.csv(assoRes,"trade-seq_assoRes_gene.csv")
head(assoRes)
###寻找与起止点相关的基因
startRes <- startVsEndTest(sce)
#write.csv(startRes,"trade-seq_startRes_gene.csv")
head(startRes)

# 按相关性进行排序
oStart <- order(startRes$waldStat, decreasing = TRUE)
#write.csv(oStart,"trade-seq_startRes_gene排序.csv")
# 挑选相关性最强的基因，并可视化

sigGeneStart <- names(sce)[oStart[1]]

plotSmoothers(sce, counts, gene = sigGeneStart)

plotGeneCount(crv, counts, gene = sigGeneStart)

plotGeneCount(crv, counts, gene = sigGeneStart)

####美化基因表达量图

# 取celltpye和配色信息
coldata <- data.frame(celltype = CD8_object_sce@colData$cd4_cd8_group)
rownames(coldata) = colnames(CD8_object_sce)

#load("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/不同分组T亚群拟时序分析和特征得分/trade-seq基因画图/trade_seq_my_data.RData")

# 把sce中的3000个细胞对应信息取出

filter_coldata <- coldata[colnames(sce),] %>% as.data.frame()

View(filter_coldata)
# 添加拟时序信息
filter_coldata$Pseudotime = sce$crv$pseudotime.Lineage1

# top6 genes
top6 <- names(sce)[oStart[1:6]]
top6_exp = sce@assays@data$counts[top6,] 
top6_exp = log2(top6_exp + 1) %>% t()

gene_names <- c("CD69","CTLA4","ENTPD1","GNLY","HAVCR2","LAG3","MKI67","PDCD1","STMN1")
gene_names_exp = sce@assays@data$counts[gene_names,] 
gene_names_exp = log2(gene_names_exp + 1) %>% t()
# 获得最终的清洁数据
plt_data = cbind(filter_coldata, gene_names_exp)
colnames(plt_data)
# [1] "celltype"   "plotcol"    "Pseudotime" "CENPF"      "HBB"       
# [6] "HMGB2"      "CD74"       "HBA1"       "MKI67" 

names(plt_data)[1]<-"celltype"

head(plt_data)
####画图
library(ggplot2)
library(ggsci)
library(ggpubr)
library(RColorBrewer)
getPalette = colorRampPalette(brewer.pal(12, "Paired"))

mycolors = c("#83B5D3","#AEE94F")

# 为了拼图美观，我们把legend隐藏掉
plt_list = list()
for (gene in gene_names) {
  print(gene)
  p = ggscatter(data = plt_data,
                x = 'Pseudotime',
                y = gene,
                color = 'celltype',
                size = 0.6)+
    geom_smooth(se = T, color = 'orange', formula=y ~ poly(x, 3),method = NULL)+
    theme_bw()+
    scale_color_manual(values = mycolors)+
    theme(legend.position = 'none')
  plt_list[[gene]] = p
}

library(patchwork)
wrap_plots(plt_list)
ggsave(filename = 'figures/05_top6_genes.pdf', width = 9, height = 5)

# 这里单独save一张有legend的图
gene = 'HBB'
p_test = ggscatter(data = plt_data,
                   x = 'Pseudotime',
                   y = gene,
                   color = 'celltype',
                   size = 0.6)+
  geom_smooth(se = F, color = 'orange')+
  theme_bw()+
  scale_color_manual(values = mycolors)+
  theme(legend.position = 'right')+
  guides(color = guide_legend(ncol = 1,  # 将图例分为两列显示
                              override.aes = list(size = 3)) )  # 调整图例展示的点大小
p_test
ggsave(plot = p_test, filename = 'figures/05_HBB_legend.pdf', width = 5, height = 5)



gene_dynamic <- list()
genes_plot <- c("S100a8","Anxa1","Ly6g","S100a9")
for(i in 1:length(genes_plot)){
  p = plotSmoothers(sce_slinghot, assays(sce_slinghot)$counts,
                    gene =genes_plot[i], alpha = 0.6, border = T, lwd = 2)+
    ggtitle(genes_plot[i])
  gene_dynamic[[i]] <- p
}


####CytoTRACE推断细胞干性----
devtools::install_local("C:/Users/Wangy/Downloads/CytoTRACE_0.3.3.tar.gz")
devtools::install_github("digitalcytometry/cytotrace", subdir = "cytotrace_r",force=TRUE) 
library(CytoTRACE2)
CD8_tex_object<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/不同分组T亚群拟时序分析和特征得分/CD8_tex_object.rds")

table(is.na(CD8_object$cd4_cd8_group))

cytotrace2_result_sce <- cytotrace2(CD8_object, 
                                    is_seurat = TRUE, 
                                    slot_type = "counts", 
                                    species = 'human',
                                    seed = 1234)

annotation <- data.frame(phenotype = CD8_object@meta.data$cd4_cd8_group) %>% 
  set_rownames(., colnames(CD8_object))

# plotting
plots <- plotData(cytotrace2_result = cytotrace2_result_sce, 
                  annotation = annotation, 
                  is_seurat = TRUE)


# 绘制CytoTRACE2_Potency的umap图
umap_raw <- as.data.frame(CD8_object@reductions$umap.harmony@cell.embeddings)

umap_raw$umapharmony_1
plot_names <- c("CytoTRACE2_UMAP", "CytoTRACE2_Potency_UMAP", 
                "CytoTRACE2_Relative_UMAP", "Phenotype_UMAP" )

for (plot_name in plot_names) {
  if (!is.null(plots[[plot_name]][[1]]$data)) {
    plots[[plot_name]][[1]]$data$umap_1 <- umap_raw$umapharmony_1
    plots[[plot_name]][[1]]$data$umap_2 <- umap_raw$umapharmony_2
  }
}

x_limits <- range(umap_raw$umapharmony_1, na.rm = TRUE)
y_limits <- range(umap_raw$umapharmony_2, na.rm = TRUE)


plots$CytoTRACE2_UMAP[[1]] <- plots$CytoTRACE2_UMAP[[1]] + 
  scale_x_continuous(limits = x_limits) + 
  scale_y_continuous(limits = y_limits)
plots$CytoTRACE2_UMAP

plots$CytoTRACE2_Potency_UMAP[[1]] <- plots$CytoTRACE2_Potency_UMAP[[1]] + 
  coord_cartesian(xlim = x_limits, ylim = y_limits)
plots$CytoTRACE2_Potency_UMAP

plots$CytoTRACE2_Relative_UMAP[[1]] <- plots$CytoTRACE2_Relative_UMAP[[1]] + 
  coord_cartesian(xlim = x_limits, ylim = y_limits)
plots$CytoTRACE2_Relative_UMAP

plots$Phenotype_UMAP[[1]] <- plots$Phenotype_UMAP[[1]] + 
  coord_cartesian(xlim = x_limits, ylim = y_limits)
plots$Phenotype_UMAP
Idents(cytotrace2_result_sce) <- ""

p1<-FeaturePlot(cytotrace2_result_sce, "CytoTRACE2_Relative",pt.size = 1.5,
            reduction = "umap.harmony") + 
  scale_colour_gradientn(colours = 
                           (c("#9E0142", "#F46D43", "#FEE08B", "#E6F598", 
                              "#66C2A5", "#5E4FA2")), 
                         na.value = "transparent", 
                         limits = c(0, 1), 
                         breaks = seq(0, 1, by = 0.2), 
                         labels = c("0.0 (More diff.)", 
                                    "0.2", "0.4", "0.6", "0.8", "1.0 (Less diff.)"), 
                         name = "Relative\norder \n", 
                         guide = guide_colorbar(frame.colour = "black", 
                                                ticks.colour = "black")) + 
  ggtitle("CytoTRACE 2") + 
  xlab("UMAP1") + ylab("UMAP2") + 
  theme(legend.text = element_text(size = 10), 
        legend.title = element_text(size = 12), 
        axis.text = element_text(size = 13), 
        axis.title = element_text(size = 12), 
        plot.title = element_text(size = 12, 
                                  face = "bold", hjust = 0.5, 
                                  margin = margin(b = 20)),
        axis.line = element_line(linewidth = 0.7)
        ) + 
  theme(aspect.ratio = 1)

p1

class_avg <- data.frame(
  umapharmony_1  = c(0.70586, -5.1612144),
  umapharmony_2  = c(1.24, -1.694),
  cd4_cd8_group = c("CD8_C1-GNLY", "CD8_C5-STMN1")
)

p2<-p1+ 
  ggrepel::geom_label_repel(data = class_avg, 
                            aes(x = umapharmony_1, y = umapharmony_2, label = cd4_cd8_group),
                            size = 5, 
                            box.padding = 0.5, 
                            point.padding = 0.5, 
                            segment.color = 'grey50')+theme(
                              legend.position = "right",        # 图例放置在顶部
                            )

p2

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/不同分组T亚群拟时序分析和特征得分/CD8_cytoTRACE2细胞干性推断.pdf",p2,width=7,height = 5)

library(ggpubr)
c("#83B5D3","#AEE94F")
p1 <- ggboxplot(cytotrace2_result_sce@meta.data, x="cd4_cd8_group", y="CytoTRACE2_Score", width = 0.6, 
                color = "cd4_cd8_group",#轮廓颜色
                fill="cd4_cd8_group",#填充
                palette = c("#83B5D3","#AEE94F"),
                xlab = F, #不显示x轴的标签
                bxp.errorbar=T,#显示误差条
                bxp.errorbar.width=0.5, #误差条大小
                size=1, #箱型图边线的粗细
                outlier.shape=NA, #不显示outlier
                legend = "none") #图例放右边 


  p1<-p1+ theme(axis.title.x = element_blank(),        # 去掉 X 轴标题
             axis.text = element_text(size = 13), 
             axis.title = element_text(size = 12), 
             plot.title = element_text(size = 12, 
                                       face = "bold", hjust = 0.5, 
                                       margin = margin(b = 20)),
             axis.line = element_line(linewidth = 0.7)
)
###指定组比较
#my_comparisons <- list(c("Epi", "un"), c("T", "un"),c("Myeloid", "un"))
cytotrace2_result_sce@meta.data$cd4_cd8_group
my_comparisons <- list(c("CD8_C1-GNLY", "CD8_C5-STMN1"))
p2<-p1+stat_compare_means(comparisons = my_comparisons,
                      method = "wilcox.test")

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/不同分组T亚群拟时序分析和特征得分/CD8_cytoTRACE2细胞干性推断boxplot.pdf",p2,width=4,height = 5)


####tradeSeq 差异基因画图----
CD8_object_sce
head(colnames(colData(CD8_object_sce)))

plot(reducedDims(CD8_object_sce)$UMAP.HARMONY, pch=16, asp = 1)
mycolors = brewer.pal(4,"Set1")
lines(SlingshotDataSet(CD8_object_sce), lwd=2, col=mycolors)
library(tradeSeq)

# Fit negative binomial model
counts <- CD8_object_sce@assays@data$counts
crv <- SlingshotDataSet(CD8_object_sce)
set.seed(111)
icMat <- evaluateK(counts = counts, 
                   sds = crv, 
                   k = 3:10,    # no more than 12
                   nGenes = 200, # 每个细胞纳入分析的基因数量，默认是500，这里为了节省示例计算时间
                   verbose = T)

set.seed(111)
pseudotime <- slingPseudotime(crv, na = FALSE)
cellWeights <- slingCurveWeights(crv)
# fit negative binomial GAM
# 2k cells ~13 min
# system.time()这个函数可以计算运行时间
system.time({ 
  sce <- fitGAM(counts = counts, 
                pseudotime = pseudotime, 
                cellWeights = cellWeights,
                nknots = 6, 
                verbose = FALSE)
})
table(rowData(sce)$tradeSeq$converged)

oStart <- order(startRes$waldStat, decreasing = TRUE)

# 挑选相关性最强的基因，并可视化
sigGeneStart <- names(sce)[oStart[1]]
plotSmoothers(sce, counts, gene = sigGeneStart)

plotGeneCount(crv, counts, gene = sigGeneStart)

####trade-seq拟时序差异基因----
# 取celltpye和配色信息
coldata <- data.frame(celltype = CD8_object_sce@colData$cd4_cd8_group,
                      CD8_object_sce@colData$slingPseudotime_1)
rownames(coldata) = colnames(CD8_object_sce)
CD8_object_sce$slingPseudotime_1
# 把sce中的3000个细胞对应信息取出
# filter_coldata <- coldata[colnames(CD8_object_sce),]

# 添加拟时序信息
# filter_coldata$Pseudotime = sce$crv$pseudotime.Lineage1

# top6 genes

top6 <- c("HAVCR2","LAG3","PDCD1","ENTPD1","TBX21","REL","EOMES","TOX","NFKB1","FOXO1","SPRY1")
top6<-c("CD69",
        "PDCD1",
        "HAVCR2",
        "ENTPD1",
        "CTLA4",
        "LAG3")
top6_exp =CD8_object_sce@assays@data$counts[top6,] 
top6_exp = log2(top6_exp + 1) %>% t()

# 获得最终的清洁数据
plt_data = cbind(coldata, top6_exp)
colnames(plt_data)

# [1] "celltype"   "plotcol"    "Pseudotime" "CENPF"      "HBB"       
# [6] "HMGB2"      "CD74"       "HBA1"       "MKI67"  

####trade拟时序基因画图----
library(ggplot2)
library(ggsci)
library(ggpubr)
library(RColorBrewer)
getPalette = colorRampPalette(brewer.pal(12, "Paired"))
mycolors = getPalette(length(unique(plt_data$celltype)))
names(plt_data)[2]<-"Pseudotime"
# 为了拼图美观，我们把legend隐藏掉
plt_list = list()
for (gene in top6) {
  print(gene)
  p = ggscatter(data = plt_data,
                x = 'Pseudotime',
                y = gene,
                color = 'celltype',
                size = 0.6)+
    geom_smooth(se = F, color = 'orange')+
    theme_bw()+
    scale_color_manual(values = mycolors)+
    theme(legend.position = 'none')
  plt_list[[gene]] = p
}

library(patchwork)
wrap_plots(plt_list)
ggsave(filename = 'figures/05_top6_genes.pdf', width = 9, height = 5)

assoRes <- associationTest(CD8_object_sce)

library(ggplot2)
library(ggsci)
library(ggpubr)
library(RColorBrewer)
getPalette = colorRampPalette(brewer.pal(12, "Paired"))
mycolors = getPalette(length(unique(plt_data$celltype)))

# 为了拼图美观，我们把legend隐藏掉
plt_list = list()
for (gene in top6) {
  print(gene)
  p = ggscatter(data = plt_data,
                x = 'Pseudotime',
                y = gene,
                color = 'celltype',
                size = 0.6)+
    geom_smooth(se = F, color = 'orange')+
    theme_bw()+
    scale_color_manual(values = mycolors)+
    theme(legend.position = 'none')
  plt_list[[gene]] = p
}

library(patchwork)
wrap_plots(plt_list)
ggsave(filename = 'figures/05_top6_genes.pdf', width = 9, height = 5)

####
View(CD8_object)
names(pData(sce_CDS))
monocle_data <- pData(sce_CDS)[,c(1,35,91,92)] %>% as.data.frame()
names(monocle_data)
write.csv(monocle_data,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/不同分组T亚群拟时序分析和特征得分/monocle_data拟时序数据.csv")

trade_coldata <- data.frame(orig.ident = CD8_object_sce@colData$orig.ident,
                            cell_id = CD8_object_sce@colData$cell_id,
                            slingPseudotime = CD8_object_sce@colData$slingPseudotime_1)

View(CD8_tex_object)

metadata<-CD8_tex_object@meta.data

merge_data<- cbind(metadata,trade_coldata[match(metadata$cell_id,trade_coldata$cell_id),],monocle_data[match(metadata$cell_id,monocle_data$cell_id),])%>% as.data.frame()

names(merge_data)

####细胞通讯分析----


####topmarker 热图

T_cells_object$cd4_cd8_group
Idents(T_cells_object)<-"cd4_cd8_group"
pbmc.markers <- FindAllMarkers(T_cells_object,
                               only.pos = TRUE,
                               min.pct = 0.25,
                               logfc.threshold = 0.25)

# get top 10 genes
top5pbmc.markers <- pbmc.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)

col <- pal_npg()(9)
#data$rid <- factor(data$rid, levels=c("CXCL2+iCAFs","CXCR4+iCAFs","TPSAB1+iCAFs","MyCAFs","apCAFs","eCAFs"))

cell_type_color<-c("CD4_C1-PDCD1"="#DF66B0",
                   "CD4_C2-IL7R"="#9888DB",
                   "CD4_C3-FOXP3"="#DC8062",
                   "CD4_C4-FOSB"="#D8ACC9",
                   "CD8_C1-GNLY"="#83B5D3",
                   "CD8_C2-CTSW"="#D5D8C5",
                   "CD8_C3-IL7R"="#77DCCF",
                   "CD8_C4-GZMK"="#BA4EE4",
                   "CD8_C5-STMN1"="#AEE94F",
                   "CD3+CD4-CD8-"="#89DB8E",
                   "NK"="#E0D175")

# plot

tcell_marker<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/marker热图分型/T_NK_marker.csv")

DoHeatmap(T_cells_object,features =tcell_marker$gene_name,
          group.colors = cell_type_color) +
  scale_colour_npg() +
  scale_fill_gradient2(low = '#0099CC',mid = 'white',high = '#CC0033',
                       name = 'Z-score')
top5pbmc.markers$gene[50]<-"GNLY"
top5pbmc.markers$gene[51]<-"CTSW"
mean_gene_exp <- AverageExpression(T_cells_object,
                                   features = top5pbmc.markers$gene,
                                   group.by = 'cd4_cd8_group',
                                   slot = 'data') %>%
  data.frame() %>%
  as.matrix()

# add colnames

# Z-score
htdf <- t(scale(t(mean_gene_exp),scale = T,center = T))

# color
col_fun = colorRamp2(c(-2, 0, 2), c("#0099CC", "white", "#CC0033"))

anno_col <- brewer.pal(9, "Paired")
names(anno_col) <- paste('clsuter ',0:8,sep = '')
colnames(htdf)<- c("CD3+CD4-CD8-","CD4_C1-PDCD1","CD4_C2-IL7R","CD4_C3-FOXP3","CD4_C4-FOSB", 
                   "CD8_C1-GNLY","CD8_C2-CTSW","CD8_C3-IL7R","CD8_C4-GZMK","CD8_C5-STMN1",
                   "NK")
# top annotation
column_ha = HeatmapAnnotation(cluster = colnames(htdf),
                              col = list(cluster = cell_type_color))


Heatmap(htdf,
        name = "Z-score",
        cluster_columns = F,cluster_rows = F,
        row_title = "Cluster top 5 Marker genes",
        column_title = "Clusters",
        row_names_gp = gpar(fontface = 'italic',fontsize = 10),
        row_names_side = 'left',
        border = T,
        rect_gp = gpar(col = "white", lwd = 1),
        column_names_side = 'top',
        column_names_rot = 45,
        top_annotation = column_ha,
        # column_split = paste('clsuter ',0:8,sep = ''),
        col = col_fun)

Idents(T_cells_object)<-"cd4_cd8_group"
unique(Idents(T_cells_object))
CD8_object<- subset(T_cells_object,idents=c("CD8_C4-GZMK",
                                            "CD8_C2-CTSW",
                                            "CD8_C3-IL7R",
                                            "CD8_C1-GNLY",
                                            "CD8_C5-STMN1"))

gene_list<-list(c("CD8A","CD8B","PDCD1","GNLY"))
gene_list<-list(c("CD8A","CD8B","PDCD1","CTSW"))
gene_list<-c(deg_all$X %>% head(10),c("CD8A","CD8B","PDCD1"))
gene_list<-list(c("IL7R","MKI67","PDCD1","ZNF683","HOPX","ITGAE","TCF7","SLAMF6","CD8A","CD8B","PDCD1","GNLY"))
gene_list<-list(c("IL7R","GPR183","LMNA","NR4A3","TCF7","MGAT4A","CD55",
                  "AIM1","PER1","FOSL2","EGR1","TSPYL2","YPEL5","CSRNP1","REL",
                  "SKIL","PIK3R1","FOXP1","RGCC","PFKFB3","MYADM","ZFP36L2",
                  "USP36","TC2N","FAM177A1","BTG2","TSC22D2","FAM65B","STAT4","RGPD5",
                  "NEU1","IFRD1","PDE4B","NR4A1","CD8A","CD8B","PDCD1","GNLY"))

gene_list<-list(c("CD69",
                  "PDCD1",
                  "HAVCR2",
                  "ENTPD1",
                  "CTLA4",
                  "LAG3",
                  "TIGIT",
                  "TBX21"
))
imm_T2<-AddModuleScore(CD8_object,features=gene_list,name="CD8_GNLY_signature")
head(imm_T2)

library(tidyverse)
library(gghalves)
library(rstatix)
library(ggpubr)
library(ggsci)
library(Seurat)

sig_score<- cbind(imm_T2$MPR_Response2,imm_T2$CD8_GNLY_signature1) %>% as.data.frame()

names(sig_score)<- c("MPR_Response2","CD8_GNLY_signature1")

df<-sig_score
df$CD8_GNLY_signature1 <- as.numeric(df$CD8_GNLY_signature1)
names(df)
df <- na.omit(df)

df_p_val1 <- df %>%
  wilcox_test(CD8_GNLY_signature1~ MPR_Response2) %>%
  adjust_pvalue(p.col = "p", method = "bonferroni") %>%
  add_significance(p.col = "p.adj") %>%
  add_xy_position(x = "MPR_Response2", dodge = 0.8)

treatment_color <- c("#4974a4","#4dae47","#f29600","#F39B7Fcc")
treatment_color <- c("#4974a4","#f29600")
treatment_color <- c("#4974a4","#4dae47","#f29600")
unique(df$MPR_Response2)
df$MPR_Response2<- factor(df$MPR_Response2,levels = c("Pre_NMPR","Pre_MPR","Post_MPR"))
View(df)

t_exhausted<- df %>%
  ggplot(aes(MPR_Response2,CD8_GNLY_signature1)) +
  geom_half_boxplot(fill=treatment_color,color="black",side="l",errorbar.draw = T,
                    outlier.shape = NA,width=0.8,lwd= 0.7) +
  geom_half_point(color=treatment_color,side = "r",alpha=0.1,
                  transformation_params = list(height = 0,width = 0.001,seed = 2))+
  stat_pvalue_manual(df_p_val1,label = "p.adj.signif",label.size=5,hide.ns = F,bracket.size = 0.5,y.position = c(2.4,2.5,2.6))+
  #facet_wrap(.~continent,nrow=1)+
  #scale_y_continuous(limits = c(0,5),breaks = seq(0,2,5))+
  scale_fill_npg()+
  scale_color_npg()+
  labs(title = "CD8_GNLY_signature score",y='Signature score',x="")+
  theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),units=,"cm"),
        plot.title = element_text(hjust = 0.5, size = 13),
        panel.background = element_blank(),
        plot.background = element_blank(),
        strip.text = element_text(size=12),
        axis.title.y = element_text(size = 14),
        axis.line = element_line(color = "black",size = 0.4),
        axis.text.y = element_text(color="black",size=12),
        axis.text.x = element_text(margin = margin(t = 2),color="black",size=12),
        legend.position = "none",
        panel.spacing = unit(0,"lines"),
        axis.ticks.x=element_line(color="black",size=0.5,lineend = 1),
        axis.line.y=element_line(linetype=1,color="black",size=0.7),
        axis.ticks.y=element_line(color="black",size=0.5,lineend = 1),
        axis.line.x=element_line(linetype=1,color="black",size=0.7))+
  
  coord_cartesian()
t_exhausted

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/多队列验证/EGAS50000000037/特征得分/CD8_GNLY_signature.pdf",height = 5,width=5)


deg_all= FindMarkers(CD8_object,ident.1 = c("CD8_C1-GNLY"))

write.csv(deg_all,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/CD8_deg/CD8_C1-GNLY_deg.csv")

deg_all<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/CD8_deg/CD8_C1-GNLY_deg.csv")

deg_all<- deg_all[order(deg_all$pct.1,decreasing = T),]
significant_genes <- deg_all[deg_all$avg_log2FC > 2 & deg_all$p_val_adj < 0.05, ]

write.csv(significant_genes,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/CD8_deg/CD8_C1-GNLY_deg_排序.csv")

deg_all<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/CD8_deg/CD8_C1-GNLY_deg_排序.csv")


##组织驻留记忆T细胞----
Tissue_resident_memory <-c("CA10",
                           "ITGA1",
                           "ITGAE",
                           "IL2",
                           "IL10",
                           "CXCR6",
                           "CXCL13",
                           "KCNK5",
                           "RGS1",
                           "CRTAM",
                           "DUSP6",
                           "PDCD1",
                           "IL23R") 


Cytotoxicity<- c("PRF1",
                 "GZMB",
                 "GZMA",
                 "GZMH",
                 "NKG7",
                 "GNLY")

NK_cell_receptors<-c("KLRD1",
                     "FGFBP2",
                     "FCGR3A",
                     "S1PR5",
                     "KLRC1",
                     "KLRC3",
                     "KLRB1",
                     "KLRC2")


Inhibitory_Receptors<-c("PDCD1",
                        "CTLA4",
                        "HAVCR2",
                        "LAG3",
                        "TIGIT")

####Signature gene sets from (E) used to analyze bulk RNA-seq data from an independent ICB-treated patient cohort (metastatic urothelial cancer patients treated with anti-PD-L1)
Tx_E_signature<-c("KLRD1","CXCL13","GNLY","KLRC2","CD7",
                  "ITGAE","GZMB","CTLA4","ZNF683","TNFRSF18",
                  "CD226","SIRPG","KLRB1","PRF1","TIGIT","LAG3")                                                                                                              
Tx_NE_signature<-c("GZMK","CS77","EOMES","GZMM","SAMD3","LIME1","CXCR3","SH2D1A","TRAT1","CD5")  


m1m2_pws<- list(Tx_E_signature,Tx_NE_signature)

imm_anno<-AddModuleScore(object=T_cells_object,features=m1m2_pws,name=c("Tx_E_signature","Tx_NE_signature"),nbin=12)

head(imm_anno)

use_colors<-c(
  "CD4_C1-PDCD1"="#DF66B0", 
  "CD4_C2-IL7R"="#9888DB", 
  "CD4_C3-FOXP3"="#DC8062", 
  "CD4_C4-FOSB"="#D8ACC9",
  "CD8_C1-GNLY"="#83B5D3",
  "CD8_C2-CTSW"="#D5D8C5", 
  "CD8_C3-IL7R"="#77DCCF", 
  "CD8_C4-GZMK"="#BA4EE4", 
  "CD8_C5-STMN1"="#AEE94F", 
  "CD3+CD4-CD8-"="#89DB8E",
  "NK"="#E0D175")

library(randomcoloR)
Tcell_type_color = distinctColorPalette(11)
Tcell_type_color<-c("#DF66B0","#9888DB","#DC8062","#D8ACC9","#83B5D3","#D5D8C5", "#77DCCF", "#BA4EE4",
                    "#AEE94F","#89DB8E","#E0D175")

imm_anno$cd4_cd8_group<- factor(imm_anno$cd4_cd8_group,levels=c("CD4_C1-PDCD1", "CD4_C2-IL7R","CD4_C3-FOXP3","CD4_C4-FOSB","CD8_C1-GNLY",
                                                                "CD8_C2-CTSW","CD8_C3-IL7R","CD8_C4-GZMK","CD8_C5-STMN1","CD3+CD4-CD8-","NK"))
imm_anno_df <- imm_anno@meta.data
names(imm_anno_df)
imm_anno_filter<-imm_anno_df  %>% select(cd4_cd8_group,Tx_E_signature1,Tx_NE_signature2,MPR_Response2)

library(tidyverse)
library(gghalves)
library(rstatix)
library(ggpubr)
library(ggsci)
df<-imm_anno_filter

df_p_val1 <- df %>%
  wilcox_test(Tx_E_signature1~ cd4_cd8_group) %>%
  adjust_pvalue(p.col = "p", method = "bonferroni") %>%
  add_significance(p.col = "p.adj") %>% 
  add_xy_position(x = "Tx_E_signature1", dodge = 0.8) 

View(df_p_val1)

treatment_color <- c("#4974a4","#4dae47","#f29600")
cell_type_color<- c("#DF66B0","#9888DB","#DC8062","#D8ACC9","#83B5D3","#D5D8C5", "#77DCCF", "#BA4EE4",
                    "#AEE94F","#89DB8E","#E0D175")

df$cd4_cd8_group<- factor(df$cd4_cd8_group,levels=c("CD4_C1-PDCD1", "CD4_C2-IL7R","CD4_C3-FOXP3","CD4_C4-FOSB","CD8_C1-GNLY",
                                                    "CD8_C2-CTSW","CD8_C3-IL7R","CD8_C4-GZMK","CD8_C5-STMN1","CD3+CD4-CD8-","NK"))
df$MPR_Response2<- factor(df$MPR_Response2,levels = c("Pre_NMPR","Pre_MPR","Post_MPR"))
Tissue_resident<- df %>%
  ggplot(aes(cd4_cd8_group,Tx_E_signature1)) +
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
  labs(title = "Tissue-resident memory",y='Signature score',x="")+
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
Tissue_resident


##CD8 组织驻留拟时序----
library(monocle)
library(Seurat)
library(tidyverse)
library(harmony)
library(dplyr)
library(patchwork)
library(SeuratWrappers)
# T_cells_object <- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/T_cells_object.rds")
# CD8_tex_object = subset(T_cells_object,idents = c("CD8_C1-GNLY","CD8_C5-STMN1"))

CD8_res_object = subset(T_cells_object,idents =c("CD8_C3-IL7R", "CD8_C4-GZMK"))


pbmc<- CD8_res_object
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
table(pbmc$RNA_snn_res.1)
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
pbmc <- FindNeighbors(pbmc, reduction = "umap.harmony", dims = 1:10)
pbmc <- FindClusters(pbmc)#标准聚类

#绘图和保存
p1 <- DimPlot(pbmc, group.by = "orig.ident",reduction = "umap.harmony", pt.size=0.1) +
  ggtitle("Integrated by harmony")#去批次后的可视化
p1

p1 <- DimPlot(pbmc, group.by = "cd4_cd8_group",reduction = "umap.harmony", pt.size=0.1) +
  ggtitle("Integrated by harmony")#去批次后的可视化
p1


p2 <- DimPlot(pbmc, group.by = "orig.ident",reduction = "umap", pt.size=0.1,label = T) +
  ggtitle("Integrated by harmony")#去批次后的可视化
p2

p3 <- DimPlot(pbmc, group.by = "cd4_cd8_group",reduction = "tsne.harmony", pt.size=0.1,label = T) +
  ggtitle("Integrated by harmony")#去批次后的可视化
p3

DefaultAssay(pbmc) <- "RNA"
markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
all.markers =markers %>% dplyr::select(gene, everything()) %>% subset(p_val<0.05)
top20 =all.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
Idents(pbmc) <- pbmc$RNA_snn_res.0.8
p_heatmap<-DoHeatmap(pbmc,features = top20$gene)+ #热图复现
  theme(text = element_text(size = 5))

CD8_res_object<-pbmc


saveRDS(CD8_res_object,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/不同分组T亚群拟时序分析和特征得分/CD8_res_object.rds")
CD8_res_object<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/不同分组T亚群拟时序分析和特征得分/CD8_res_object.rds")
CD8_res_object$cd4_cd8_group

####monocle2--------

ks_run_Monocle2 <- function(object, #seurat obj or expression matrix (建议数据格式转为matrix,如果数据量大转化为稀疏矩阵as(as.matrix(data), "sparseMatrix"))
                            layer, #used when object is a seurat obj
                            assay, #used when object is a seurat obj
                            lowerDetectionLimit = 0.1, 
                            VARgenesM=c("dispersionTable","seurat","differentialGeneTest"),
                            cellAnno=NULL, 
                            define_root=F,
                            root_state,
                            reverse=NULL
){
  
  
  if(class(object)[1] == 'Seurat') {
    
    data <- GetAssayData(object=object, layer=layer, assay=assay)#get expression matrix data
    
    pd <- new("AnnotatedDataFrame", data = object@meta.data)
    fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
    fd <- new("AnnotatedDataFrame", data = fData)
    
    if(all(data == floor(data))) {
      expressionFamily <- negbinomial.size()
    } else if(any(data < 0)){
      expressionFamily <- uninormal()
    } else {
      expressionFamily <- tobit()
    }
    
    #Creates a new CellDateSet object.
    monocle_cds <- newCellDataSet(data,
                                  phenoData = pd, 
                                  featureData = fd,
                                  lowerDetectionLimit=0.1,
                                  expressionFamily=expressionFamily)
    
  }else{
    
    print("This fucntions only apply for a seurat obj")
  }
  
  
  
  #Estimate size factors and dispersions
  #数据处理
  monocle_cds <- estimateSizeFactors(monocle_cds)#size facotr标准化细胞之间的mRNA差异
  monocle_cds <- estimateDispersions(monocle_cds)
  
  #质量控制-filter cells
  monocle_cds <- detectGenes(monocle_cds, min_expr = 0.1)
  # print(head(fData(monocle_cds)))
  # print(head(pData(monocle_cds)))
  # expressed_genes <- row.names(subset(fData(mouse_monocle), num_cells_expressed >= 10))
  monocle_cds <- monocle_cds[fData(monocle_cds)$num_cells_expressed >= 10, ]
  
  
  #select methods for VariableFeatures
  if(VARgenesM=="dispersionTable"){
    
    disp_table <- dispersionTable(monocle_cds)
    ordering_genes <- subset(disp_table,
                             mean_expression >= 0.1 &
                               dispersion_empirical >= 1.5* dispersion_fit)$gene_id
    
  }
  
  
  if(VARgenesM=="seurat"){
    
    ordering_genes <- VariableFeatures(FindVariableFeatures(object, assay = "RNA"), assay = "RNA")
    
  }
  
  
  if(VARgenesM=="differentialGeneTest"){
    
    diff_test_res <- differentialGeneTest(monocle_cds,fullModelFormulaStr = paste0("~",cellAnno))##~后面是表示对谁做差异分析的变量
    diff_test_res_sig <- diff_test_res[order(diff_test_res$qval,decreasing=F),]
    
    ordering_sce <- diff_test_res_sig[diff_test_res_sig$qval< 0.01,]
    
    if(nrow(ordering_sce)>3000){
      
      ordering_genes <- ordering_sce$gene_short_name[1:3000]
      
    }else{
      
      ordering_genes <- rdering_sce$gene_short_name
    }
    
  }
  
  #Marks genes for clustering
  monocle_cds <- setOrderingFilter(monocle_cds, ordering_genes)
  plot_ordering_genes(monocle_cds)
  
  
  #cluster
  monocle_cds <- reduceDimension(monocle_cds, max_components = 2,reduction_method = 'DDRTree')
  
  #order cells
  monocle_cds <- orderCells(monocle_cds, reverse=reverse)
  
  if(define_root){
    monocle_cds <- monocle_cds <- orderCells(monocle_cds,root_state = root_state)
  }
  
  
  return(monocle_cds)
  
}

####CD8 序----
library(monocle)
library(Seurat)
library(tidyverse)
library(harmony)
library(dplyr)
library(patchwork)
library(SeuratWrappers)
# T_cells_object <- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/T_cells_object.rds")
# CD8_tex_object = subset(T_cells_object,idents = c("CD8_C1-GNLY","CD8_C5-STMN1"))
sce_CDS <- ks_run_Monocle2(object =CD8_object,
                           layer = 'counts',
                           assay = "RNA",
                           VARgenesM="seurat",
                           cellAnno = "cd4_cd8_group")

sce_CDS$cd4_cd8_group
saveRDS(sce_CDS,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/不同分组T亚群拟时序分析和特征得分/CD8IL7R拟时序.rds")

plot_cell_trajectory(sce_CDS,color_by = "Pseudotime")
plot_cell_trajectory(sce_CDS,color_by = "State")

plot_cell_trajectory(sce_CDS,color_by = "cd4_cd8_group")


sce_CDS

####美化cd8_plot_cell_trajectory----

colour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",  
         "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0",
         "#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE",
         "#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")

colour=c("#83B5D3","#AEE94F")

p1 <- plot_cell_trajectory(sce_CDS, x = 1, y = 2, cell_size=2,cell_link_size = 1,
                           color_by = "cd4_cd8_group") + 
  theme(legend.position='none',panel.border = element_blank()) + #去掉第一个的legend
  scale_color_manual(values = colour) +
  theme(axis.text = element_text(color = 'black', size = 12),        
        axis.title = element_text(color = 'black', size = 14),       
        strip.text = element_text(color = 'black', size = 14),
        axis.line.x.bottom = element_line(linewidth = 0.6, color = "black"),
        axis.line.y.left = element_line(linewidth = 0.6, color = "black"),
        axis.line = element_line(linewidth = 2, color = "black")) 
p2 <- plot_complex_cell_trajectory(sce_CDS, x = 1, y = 2,cell_size=2,cell_link_size = 1,
                                   color_by = "cd4_cd8_group")+
  scale_color_manual(values = colour) +
  theme(legend.title = element_blank()) +
  theme(axis.text = element_text(color = 'black', size = 12),        
        axis.title = element_text(color = 'black', size = 14),       
        strip.text = element_text(color = 'black', size = 14),
        axis.line.x.bottom = element_line(linewidth = 0.6, color = "black"),
        axis.line.y.left = element_line(linewidth = 0.6, color = "black"),
        axis.line = element_line(linewidth = 2, color = "black"),
        legend.title = element_text(size = 14),   # 图例标题大小
        legend.text = element_text(size = 12),    # 图例文本大小
        legend.position = "bottom",               # 图例位置，例如底部
        legend.key.size = unit(1, "cm")) +
  guides(color = guide_legend(override.aes = list(size = 3)))

p3<- p1 / p2

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/不同分组T亚群拟时序分析和特征得分/monoclecd8拟时序加树.pdf",p3,width = 8,height=9)





####monocle3------------

library(monocle3)
library(Seurat)
library(ggplot2)
library(dplyr)
install.packages("devtools")
devtools::install_github('cole-trapnell-lab/monocle3')

CD8_object<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/T_cells_CD8_object.rds")

CD8_res_object<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/CD8_res_object.rds")

mouse_data<-CD8_object

expression_matrix <- as(as.matrix(mouse_data@assays$RNA@counts), 'sparseMatrix')
cell_metadata <- mouse_data@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(expression_matrix))
rownames(gene_annotation) <- rownames(expression_matrix)

#构建cds对象
cds <- new_cell_data_set(expression_matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)

cds <- preprocess_cds(cds, num_dim = 100)#method默认为PCA
plot_pc_variance_explained(cds)#展示PC数，和seurat降维一摸一样
cds <- reduce_dimension(cds, reduction_method = "UMAP", preprocess_method = "PCA")
View(mouse_data)
cds.embed <- cds@int_colData$reducedDims$UMAP#monocle3中UMAP信息
int.embed <- Embeddings(mouse_data, reduction = "umap.harmony")#seurat中UMAP降维信息
int.embed <- int.embed[rownames(cds.embed),]#替换
cds@int_colData$reducedDims$UMAP <- int.embed #替换

plot_cells(cds, color_cells_by="cd4_cd8_group",
           cell_size=0.5,group_label_size=4) 

mycds <- cds

plot_cells(mycds, color_cells_by="cd4_cd8_group",
           cell_size=0.5,group_label_size=4) 


?plot_cells

mycds <- cluster_cells(mycds, reduction_method = "UMAP")
?cluster_cells
mycds <- learn_graph(mycds,
                     verbose=T,
                     learn_graph_control=list(minimal_branch_len=10,#在图修剪过程中要保留的分支直径路径的最小长度。默认值是10。
                                              euclidean_distance_ratio=1#生成树中两个末端节点的欧氏距离与生成树上允许连接的任何连接点之间的最大距离之比。默认值为1。
                     ))

?learn_graph

View(mycds)
plot_cells(mycds, 
           color_cells_by = 'cd4_cd8_group',
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           cell_size=0.5,group_label_size=4)


plot_cells(mycds, 
           color_cells_by = "cd4_cd8_group", 
           label_groups_by_cluster=FALSE,
           label_leaves=TRUE, 
           label_branch_points=TRUE,
           graph_label_size=4)

mycds1 <- mycds
install.packages("htmltools")
remotes::install_github("rstudio/htmltools")
packageVersion("htmltools")
.rs.restartR()
remove.packages("htmltools")
install.packages("htmltools")


mycds1 <- order_cells(mycds1)#在交互界面中选择

#可视化拟时图
plot_cells(mycds1, label_cell_groups = F, 
           color_cells_by = "pseudotime", 
           label_leaves = F, 
           label_branch_points = F, 
           graph_label_size = 2, 
           cell_size=0.2, 
           trajectory_graph_segment_size = 2)

#或者函数推断
mycds2 <- mycds
get_earliest_principal_node <- function(cds){
  cell_ids <- dim(cds)[2]
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))]
  root_pr_nodes
}

mycds2 <- order_cells(mycds2, root_pr_nodes=get_earliest_principal_node(mycds2))

monocle_plot1<- plot_cells(mycds2, label_cell_groups = F, 
           color_cells_by = "pseudotime", 
           label_leaves = F, 
           label_branch_points = F, 
           graph_label_size = 3, 
           cell_size=0.3, 
           trajectory_graph_segment_size = 1)

monocle_plot2<-plot_cells(mycds2, label_cell_groups = F, 
           color_cells_by = "cd4_cd8_group", 
           label_leaves = F, 
           label_branch_points = F, 
           graph_label_size = 3, 
           cell_size=0.3, 
           trajectory_graph_segment_size = 1)

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/组织驻留T细胞/monocle/CD8整体拟时序monocle_plot1.pdf",monocle_plot1,width=6,height=5)

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/组织驻留T细胞/monocle/CD8整体拟时序monocle_plot2.pdf",monocle_plot2,width=6,height=5)


save(mycds, mycds1, mycds2, file = 'D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/组织驻留T细胞/monocle/Monocle3_cds_peu.RData')

load('D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/组织驻留T细胞/monocle/Monocle3_cds_peu.RData')

####组织驻留T细胞拟时序----

CD8_res_object<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/CD8_res_object.rds")

mouse_data<-CD8_res_object

expression_matrix <- as(as.matrix(mouse_data@assays$RNA@counts), 'sparseMatrix')
cell_metadata <- mouse_data@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(expression_matrix))
rownames(gene_annotation) <- rownames(expression_matrix)

#构建cds对象
cds <- new_cell_data_set(expression_matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)

cds <- preprocess_cds(cds, num_dim = 100)#method默认为PCA
plot_pc_variance_explained(cds)#展示PC数，和seurat降维一摸一样
cds <- reduce_dimension(cds, reduction_method = "UMAP", preprocess_method = "PCA")
View(mouse_data)
cds.embed <- cds@int_colData$reducedDims$UMAP#monocle3中UMAP信息
int.embed <- Embeddings(mouse_data, reduction = "umap.harmony")#seurat中UMAP降维信息
int.embed <- int.embed[rownames(cds.embed),]#替换
cds@int_colData$reducedDims$UMAP <- int.embed #替换

plot_cells(cds, color_cells_by="cd4_cd8_group",
           cell_size=0.5,group_label_size=4) 

mycds <- cds

plot_cells(mycds, color_cells_by="cd4_cd8_group",
           cell_size=0.5,group_label_size=4) 


?plot_cells

mycds <- cluster_cells(mycds, reduction_method = "UMAP")
?cluster_cells
mycds <- learn_graph(mycds,
                     verbose=T,
                     learn_graph_control=list(minimal_branch_len=10,#在图修剪过程中要保留的分支直径路径的最小长度。默认值是10。
                                              euclidean_distance_ratio=1#生成树中两个末端节点的欧氏距离与生成树上允许连接的任何连接点之间的最大距离之比。默认值为1。
                     ))

?learn_graph

View(mycds)
plot_cells(mycds, 
           color_cells_by = 'cd4_cd8_group',
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           cell_size=0.5,group_label_size=4)


plot_cells(mycds, 
           color_cells_by = "cd4_cd8_group", 
           label_groups_by_cluster=FALSE,
           label_leaves=TRUE, 
           label_branch_points=TRUE,
           graph_label_size=4)

mycds1 <- mycds
install.packages("htmltools")
remotes::install_github("rstudio/htmltools")
packageVersion("htmltools")
.rs.restartR()
remove.packages("htmltools")
install.packages("htmltools")


mycds1 <- order_cells(mycds1)#在交互界面中选择

#可视化拟时图
plot_cells(mycds1, label_cell_groups = F, 
           color_cells_by = "pseudotime", 
           label_leaves = F, 
           label_branch_points = F, 
           graph_label_size = 2, 
           cell_size=0.2, 
           trajectory_graph_segment_size = 2)

monocle_plot1<- plot_cells(mycds1, label_cell_groups = F, 
                           color_cells_by = "pseudotime", 
                           label_leaves = F, 
                           label_branch_points = F, 
                           graph_label_size = 3, 
                           cell_size=0.3, 
                           trajectory_graph_segment_size = 1)


colors <- c('CD8_C3-IL7R' = '#FF7F24', 'CD8_C4-GZMK' = '#9400D3')

# 确保分组名称与颜色向量中的名称一致
# 更新 cds 对象的元数据以包含颜色信息
mycds1$cd4_cd8_group_color <- colors[as.character(mycds1$cd4_cd8_group)]

# 使用 plot_cells 函数并传递 color_cells_by 参数
monocle_plot2 <- plot_cells(mycds1, 
                            label_cell_groups = FALSE,
                            color_cells_by = "cd4_cd8_group", # 注意这里使用了新的颜色变量
                            label_leaves = FALSE,
                            label_branch_points = FALSE,
                            graph_label_size = 3,
                            cell_size = 0.3,
                            trajectory_graph_segment_size = 1)


monocle_plot2<-plot_cells(mycds1, label_cell_groups = F, 
                          color_cells_by = "cd4_cd8_group", 
                          label_leaves = F, 
                          label_branch_points = F, 
                          graph_label_size = 3, 
                          cell_size=0.3, 
                          trajectory_graph_segment_size = 1,
                          group_colors = c('#FF7F24',"#9400D3"))

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/组织驻留T细胞/monocle/组织驻留T_monocle_plot1.pdf",monocle_plot1,width=6,height=5)

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/组织驻留T细胞/monocle/组织驻留T_monocle_plot2.pdf",monocle_plot2,width=6,height=5)


#或者函数推断
mycds2 <- mycds
get_earliest_principal_node <- function(cds){
  cell_ids <- dim(cds)[2]
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))]
  root_pr_nodes
}

mycds2 <- order_cells(mycds2, root_pr_nodes=get_earliest_principal_node(mycds2))

monocle_plot1<- plot_cells(mycds2, label_cell_groups = F, 
                           color_cells_by = "pseudotime", 
                           label_leaves = F, 
                           label_branch_points = F, 
                           graph_label_size = 3, 
                           cell_size=0.3, 
                           trajectory_graph_segment_size = 1)

monocle_plot2<-plot_cells(mycds2, label_cell_groups = F, 
                          color_cells_by = "cd4_cd8_group", 
                          label_leaves = F, 
                          label_branch_points = F, 
                          graph_label_size = 3, 
                          cell_size=0.3, 
                          trajectory_graph_segment_size = 1)

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/组织驻留T细胞/monocle/组织驻留T_monocle_plot1.pdf",monocle_plot1,width=6,height=5)

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/组织驻留T细胞/monocle/组织驻留T_monocle_plot2.pdf",monocle_plot2,width=6,height=5)


####耗竭T细胞monocle 拟时序----

CD8_tex_object<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/不同分组T亚群拟时序分析和特征得分/CD8_tex_object.rds")

mouse_data<-CD8_tex_object

expression_matrix <- as(as.matrix(mouse_data@assays$RNA@counts), 'sparseMatrix')
cell_metadata <- mouse_data@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(expression_matrix))
rownames(gene_annotation) <- rownames(expression_matrix)

#构建cds对象
cds <- new_cell_data_set(expression_matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)

cds <- preprocess_cds(cds, num_dim = 100)#method默认为PCA
plot_pc_variance_explained(cds)#展示PC数，和seurat降维一摸一样
cds <- reduce_dimension(cds, reduction_method = "UMAP", preprocess_method = "PCA")
View(mouse_data)
cds.embed <- cds@int_colData$reducedDims$UMAP#monocle3中UMAP信息
int.embed <- Embeddings(mouse_data, reduction = "umap.harmony")#seurat中UMAP降维信息
int.embed <- int.embed[rownames(cds.embed),]#替换
cds@int_colData$reducedDims$UMAP <- int.embed #替换

plot_cells(cds, color_cells_by="cd4_cd8_group",
           cell_size=0.5,group_label_size=4) 

mycds <- cds

plot_cells(mycds, color_cells_by="cd4_cd8_group",
           cell_size=0.5,group_label_size=4) 


?plot_cells

mycds <- cluster_cells(mycds, reduction_method = "UMAP")
?cluster_cells
mycds <- learn_graph(mycds,
                     verbose=T,
                     learn_graph_control=list(minimal_branch_len=10,#在图修剪过程中要保留的分支直径路径的最小长度。默认值是10。
                                              euclidean_distance_ratio=1#生成树中两个末端节点的欧氏距离与生成树上允许连接的任何连接点之间的最大距离之比。默认值为1。
                     ))

?learn_graph

View(mycds)
plot_cells(mycds, 
           color_cells_by = 'cd4_cd8_group',
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           cell_size=0.5,group_label_size=4)


plot_cells(mycds, 
           color_cells_by = "cd4_cd8_group", 
           label_groups_by_cluster=FALSE,
           label_leaves=TRUE, 
           label_branch_points=TRUE,
           graph_label_size=4)

mycds1 <- mycds
install.packages("htmltools")
remotes::install_github("rstudio/htmltools")
packageVersion("htmltools")
.rs.restartR()
remove.packages("htmltools")
install.packages("htmltools")


mycds1 <- order_cells(mycds1)#在交互界面中选择

#可视化拟时图
plot_cells(mycds1, label_cell_groups = F, 
           color_cells_by = "pseudotime", 
           label_leaves = F, 
           label_branch_points = F, 
           graph_label_size = 2, 
           cell_size=0.2, 
           trajectory_graph_segment_size = 2)

#或者函数推断
mycds2 <- mycds
get_earliest_principal_node <- function(cds){
  cell_ids <- dim(cds)[2]
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))]
  root_pr_nodes
}

mycds2 <- order_cells(mycds2, root_pr_nodes=get_earliest_principal_node(mycds2))

monocle_plot1<- plot_cells(mycds2, label_cell_groups = F, 
                           color_cells_by = "pseudotime", 
                           label_leaves = F, 
                           label_branch_points = F, 
                           graph_label_size = 3, 
                           cell_size=0.5, 
                           trajectory_graph_segment_size = 1)

monocle_plot2<-plot_cells(mycds2, label_cell_groups = F, 
                          color_cells_by = "cd4_cd8_group", 
                          label_leaves = F, 
                          label_branch_points = F, 
                          graph_label_size = 3, 
                          cell_size=0.3, 
                          trajectory_graph_segment_size = 1)

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/组织驻留T细胞/monocle/耗竭T_monocle_plot1.pdf",monocle_plot1,width=6,height=5)

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/组织驻留T细胞/monocle/耗竭T_monocle_plot2.pdf",monocle_plot2,width=6,height=5)


#===============================setp6：Monocle3拟时结果可视化====================
#将拟时结果添加到seurat对象

mouse_data<- CD8_res_object
pd <- pseudotime(mycds1, reduction_method = 'UMAP')
mouse_data <- AddMetaData(mouse_data,metadata = pd,col.name = 'pseudotime')#添加拟时


mycds1_res <- graph_test(mycds1, 
                         neighbor_graph="principal_graph", cores=4)

write.csv(mycds1_res, file = "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/组织驻留T细胞/monocle/mycds1_res.csv")

mycds1_res<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/组织驻留T细胞/monocle/mycds1_res.csv",header = T,row.names = 1)
head(mycds1_res)
res_ids <- row.names(subset(mycds1_res, q_value < 0.01))#显著性基因
gene_module_df <- find_gene_modules(mycds1[res_ids,], 
                                    resolution=c(10^seq(-6,-1)))#find_gene_modules函数是将这些显著基因聚类成细胞中共表达的moducle
write.csv(gene_module_df, file = "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/组织驻留T细胞/monocle/gene_module_df.csv")#保存moducle文件，可以查看具体的moducle中的基因

cell_group_df <- tibble::tibble(cell=row.names(colData(mycds1)), 
                                cell_group=colData(mycds1)$cd4_cd8_group)#细胞信息


agg_mat <- aggregate_gene_expression(mycds1, 
                                     gene_module_df, 
                                     cell_group_df)#计算基因平均表达量
row.names(agg_mat) <- stringr::str_c("Module", row.names(agg_mat))
pheatmap::pheatmap(agg_mat,
                   scale="column", clustering_method="ward.D2")


#选择top10的基因可视化，也可以自己指定基因
genes_sig <- mycds1_res %>% top_n(n=10, morans_I) %>%
  pull(gene_short_name) %>% as.character()
#显示celltype
plot_genes_in_pseudotime(mycds1[genes_sig,], color_cells_by="cd4_cd8_group", 
                         min_expr=0.5, ncol = 2)

#显示拟时
plot_genes_in_pseudotime(mycds1[genes_sig,], color_cells_by="pseudotime", 
                         min_expr=0.5, ncol = 2)


#直接将基因表达用UMAP显示
plot_cells(mycds1, genes=genes_sig, show_trajectory_graph=FALSE,
           label_cell_groups=FALSE,  label_leaves=FALSE)

#######可视化显示轨迹
library(viridis)
monocle_plot3_1<- plot_cells(mycds1,
           genes=c("IL7R","CXCR4","FOS","CD69","CCL4","CCL4L2","GZMK","ANXA1","FOSB"),           
           label_cell_groups=F,
           show_trajectory_graph=T, 
           cell_size=1, trajectory_graph_color = "black", 
           label_branch_points = F, 
           label_roots = F, label_leaves = F)+
  scale_color_viridis(option="inferno")

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/组织驻留T细胞/monocle/gene_拟时序monocle_plot3_1.pdf",monocle_plot3_1,width=7,height=5)


sel_genes <- c("IL7R","CXCR4","FOS","CD69","CCL4","CCL4L2","GZMK","ANXA1","FOSB")
monocle_plot3<- plot_genes_in_pseudotime(mycds1[rowData(mycds1)$gene_short_name %in% sel_genes,],
                         color_cells_by="cd4_cd8_group",
                         min_expr=0.5, ncol = 3)

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/组织驻留T细胞/monocle/gene_拟时序monocle_plot3.pdf",monocle_plot3,width=7,height=5)

#---------------------monocle3可视化3  不同分组的拟时展示-----------------------------
#我们这里展示不同的分组展示拟时

#提取graph坐标信息
ica_space_df <- t(mycds1@principal_graph_aux$UMAP$dp_mst) %>% as.data.frame() %>% dplyr::mutate(sample_name = rownames(.),sample_state = rownames(.))
pgraph <- principal_graph(mycds1)$UMAP
edge_df <- pgraph %>% igraph::as_data_frame() %>% dplyr::select_(source = "from", target = "to") %>%
  dplyr::left_join(ica_space_df %>%
                     dplyr::select_(
                       source="sample_name",
                       source_prin_graph_dim_1="umapharmony_1",
                       source_prin_graph_dim_2="umapharmony_2"),
                   by = "source") %>%
  dplyr::left_join(ica_space_df %>%
                     dplyr::select_(
                       target="sample_name",
                       target_prin_graph_dim_1="umapharmony_1",
                       target_prin_graph_dim_2="umapharmony_2"),
                   by = "target")


#设置一些参数，也可以在函数中设置
trajectory_graph_segment_size = 1
trajectory_graph_color = "red"
#提取数据
custom.plot.data <- data.frame(reducedDims(mycds1)[["UMAP"]], 
                               condition = factor(mycds1@colData$MPR_Response2, levels = c("Pre_NMPR", "Pre_MPR","Post_MPR")),
                               celltype = mycds1 @colData$cd4_cd8_group,
                               pseudotime = pseudotime(mycds1, reduction_method = "UMAP")
)
colnames(custom.plot.data)[1:2] <- c("umapharmony_1", "umapharmony_2")

#作图
library(ggplot2)
monocle_plot4 <- ggplot(custom.plot.data, aes(x = umapharmony_1, y = umapharmony_2, colour = celltype)) +
  geom_point(size = 0.2) +
  theme_classic() +
  geom_segment(aes_string(x = "source_prin_graph_dim_1", 
                          y = "source_prin_graph_dim_2", xend = "target_prin_graph_dim_1", 
                          yend = "target_prin_graph_dim_2"), size = trajectory_graph_segment_size, 
               color = I(trajectory_graph_color), linetype = "solid", na.rm = TRUE, data = edge_df) +
  guides(color = guide_legend(title = 'group', 
                              override.aes = list(size = 3)))+
  scale_colour_manual(values = c('#FF7F24',"#9400D3"))


ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/组织驻留T细胞/monocle/monocle_plot4.pdf",monocle_plot4,width=7,height=5)


monocle_plot5 <- ggplot(custom.plot.data, aes(x = umapharmony_1, y = umapharmony_2, colour = condition)) +
  geom_point(size = 0.2) +
  theme_classic() +
  geom_segment(aes_string(x = "source_prin_graph_dim_1", 
                          y = "source_prin_graph_dim_2", xend = "target_prin_graph_dim_1", 
                          yend = "target_prin_graph_dim_2"), size = trajectory_graph_segment_size, 
               color = I(trajectory_graph_color), linetype = "solid", na.rm = TRUE, data = edge_df) +
  guides(color = guide_legend(title = 'group', 
                              override.aes = list(size = 3)))+
  scale_colour_manual(values = c("#4974a4","#4dae47","#f29600"))


####样本基因比例随拟时变化相关性
library(Hmisc)
library(dplyr)
mycds_order <- mycds1

pData(mycds_order)$pseudotime <- pseudotime(mycds_order)#拟时结果
mycds_order <- mycds_order[,is.finite(pseudotime(mycds_order))]#去除inf值,因为有可能有些细胞没有纳入拟时也就是之前分析时一些灰色细胞，所以它的拟时时inf
a <- as.numeric(mycds_order$pseudotime)#拟时时间轴


#为了做后期相关分析，我们将拟时时间轴划分成不同的区间，这里我们使用0.15作为区间。
#也可以使用0.1或者0.2等等
mycds_order$pseudotime_bin <-  cut2(as.numeric(mycds_order$pseudotime), 
                                    seq(min(a),max(a), by = 0.15)) #区间划分
table(mycds_order@colData@listData[["pseudotime_bin"]])#查看区间，这里很重要，注意数据
table(mycds_order@colData@listData[["MPR_Response2"]],
      mycds_order@colData@listData[["pseudotime_bin"]])#这里查看不同分组，每个拟时区间中的细胞数量

table(mycds_order@colData@listData[["cd4_cd8_group"]],
      mycds_order@colData@listData[["pseudotime_bin"]])#这里查看不同分组，每个拟时区间中的细胞数量

#接下来使用循环计算每个区间细胞比例
bins <- levels(mycds_order$pseudotime_bin)#区间
proportion_df <- data.frame()#先建立一个空的数据框
for(i in 1:length(bins)){
  bin = bins[i]
  meta <- subset(as.data.frame(mycds_order@colData@listData), pseudotime_bin == bin)
  if(nrow(meta) == 0) { #因为我的数据有些拟时不是区间，而是一些数字，所以需要这一步筛选，自己构建
    df <- data.frame(MPR_Response2 = c("Pre_NMPR", "Pre_MPR"),
                     proportion = c(0,0),
                     bin_num = c(i, i),
                     bin= c(bins[i], bins[i]))
    proportion_df <- rbind(proportion_df, df)
  }
  
  else{
    df <- data.frame(table(meta$MPR_Response2) / nrow(meta))
    df <- df %>% dplyr::rename(c(MPR_Response2 = Var1, proportion = Freq))
    df$bin_num <- i
    df$bin <- bin
    proportion_df <- rbind(proportion_df, df)
  }
  
}

#不同组细胞比例与伪时间线的pearson相关性
proportion_cor <- cor.test(
  subset(proportion_df[,], MPR_Response2 == 'Pre_NMPR') %>% .$bin_num,
  subset(proportion_df[,], MPR_Response2 == 'Pre_NMPR') %>% .$proportion,
  method='pearson')

#作图并进行回归拟合
library(ggpubr)
p_pre_nmpr <- ggscatter(
  subset(proportion_df, MPR_Response2 == 'Pre_NMPR'),
  color = "#4974a4",
  x = 'bin_num', y = 'proportion',
  add ='reg.line',
  add.params = list(color = '#5c4b4d', fill = 'lightgray'),
  conf.int=TRUE,
  cor.coef=TRUE,
  cor.coeff.args = list(method = 'pearson', label.sep = '\n'),
  alpha = 0.8,
  size = 2.5, 
  title = "Propotion of  cells in Pre_NMPR group \n along the pseudotime") + 
  xlab('pseudotime') +
  ylab('proportion')
p_pre_nmpr


proportion_cor <- cor.test(
  subset(proportion_df[,], MPR_Response2 == 'Pre_MPR') %>% .$bin_num,
  subset(proportion_df[,], MPR_Response2 == 'Pre_MPR') %>% .$proportion,
  method='pearson')

#作图并进行回归拟合
library(ggpubr)
p_pre_mpr <- ggscatter(
  subset(proportion_df, MPR_Response2 == 'Pre_MPR'),
  color = "#4dae47",
  x = 'bin_num', y = 'proportion',
  add ='reg.line',
  add.params = list(color = '#5c4b4d', fill = 'lightgray'),
  conf.int=TRUE,
  cor.coef=TRUE,
  cor.coeff.args = list(method = 'pearson', label.sep = '\n'),
  alpha = 0.8,
  size = 2.5, 
  title = "Propotion of  cells in Pre_MPR group \n along the pseudotime") + 
  xlab('pseudotime') +
  ylab('proportion')
p_pre_mpr
names(proportion_df)
proportion_cor <- cor.test(
  subset(proportion_df[,], MPR_Response2 == 'Post_MPR') %>% .$bin_num,
  subset(proportion_df[,], MPR_Response2 == 'Post_MPR') %>% .$proportion,
  method='pearson')

#c("#4974a4","#4dae47","#f29600")
#作图并进行回归拟合
library(ggpubr)
p_post_mpr <- ggscatter(
  subset(proportion_df, MPR_Response2 == 'Post_MPR'),
  color = "#f29600",
  x = 'bin_num', y = 'proportion',
  add ='reg.line',
  add.params = list(color = '#5c4b4d', fill = 'lightgray'),
  conf.int=TRUE,
  cor.coef=TRUE,
  cor.coeff.args = list(method = 'pearson', label.sep = '\n'),
  alpha = 0.8,
  size = 2.5, 
  title = "Propotion of  cells in Post_MPR group \n along the pseudotime") + 
  xlab('pseudotime') +
  ylab('proportion')
p_post_mpr

pp_group<-p_pre_nmpr+p_pre_mpr+p_post_mpr

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/组织驻留T细胞/monocle/group细胞比例.pdf",pp_group,width=10,height=4)



bins <- levels(mycds_order$pseudotime_bin)#区间
proportion_df <- data.frame()#先建立一个空的数据框
for(i in 1:length(bins)){
  bin = bins[i]
  meta <- subset(as.data.frame(mycds_order@colData@listData), pseudotime_bin == bin)
  if(nrow(meta) == 0) { #因为我的数据有些拟时不是区间，而是一些数字，所以需要这一步筛选，自己构建
    df <- data.frame(cd4_cd8_group = c("CD8_C3-IL7R",
                                       "CD8_C4-GZMK"),
                     proportion = c(0,0),
                     bin_num = c(i, i),
                     bin= c(bins[i], bins[i]))
    proportion_df <- rbind(proportion_df, df)
  }
  
  else{
    df <- data.frame(table(meta$cd4_cd8_group) / nrow(meta))
    df <- df %>% dplyr::rename(c(cd4_cd8_group = Var1, proportion = Freq))
    df$bin_num <- i
    df$bin <- bin
    proportion_df <- rbind(proportion_df, df)
  }
  
}


proportion_cor <- cor.test(
  subset(proportion_df[,], cd4_cd8_group == 'CD8_C3-IL7R') %>% .$bin_num,
  subset(proportion_df[,], cd4_cd8_group == 'CD8_C3-IL7R') %>% .$proportion,
  method='pearson')

#作图并进行回归拟合
library(ggpubr)
p_IL7R <- ggscatter(
  subset(proportion_df, cd4_cd8_group == 'CD8_C3-IL7R'),
  color ='#FF7F24',
  x = 'bin_num', y = 'proportion',
  add ='reg.line',
  add.params = list(color = '#5c4b4d', fill = 'lightgray'),
  conf.int=TRUE,
  cor.coef=TRUE,
  cor.coeff.args = list(method = 'pearson', label.sep = '\n'),
  alpha = 0.8,
  size = 2.5, 
  title = "Propotion of CD8_C3-IL7R cells \n along the pseudotime") + 
  xlab('pseudotime') +
  ylab('proportion')
p_IL7R


proportion_cor <- cor.test(
  subset(proportion_df[,], cd4_cd8_group == 'CD8_C4-GZMK') %>% .$bin_num,
  subset(proportion_df[,], cd4_cd8_group == 'CD8_C4-GZMK') %>% .$proportion,
  method='pearson')

#作图并进行回归拟合
library(ggpubr)
p_GZMK <- ggscatter(
  subset(proportion_df, cd4_cd8_group == 'CD8_C4-GZMK'),
  color = "#9400D3",
  x = 'bin_num', y = 'proportion',
  add ='reg.line',
  add.params = list(color = '#5c4b4d', fill = 'lightgray'),
  conf.int=TRUE,
  cor.coef=TRUE,
  cor.coeff.args = list(method = 'pearson', label.sep = '\n'),
  alpha = 0.8,
  size = 2.5, 
  title = "Propotion of CD8_C4-GZMK cells \n along the pseudotime") + 
  xlab('pseudotime') +
  ylab('proportion')
p_GZMK
pp<- p_IL7R+p_GZMK
pp
ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/组织驻留T细胞/monocle/IL7R_GZMK细胞数随着拟时间变化.pdf",pp, width = 8,height=4)


####gene score或者module随拟时的变化

#将拟时添加到seurat对象
mouse_sc <- mouse_data[, rownames(mycds_order@colData)]#mycds_order就是上面分析好拟时的cds对象
mouse_sc$pseudotime <- mycds_order@colData$pseudotime
mouse_sc$pseudotime_bin <- mycds_order@colData$pseudotime_bin
DefaultAssay(mouse_sc) <- "RNA"
mouse_data<-CD8_res_object

#这里我选择了一些炎症基因
PMN_state_gene <- c("CA10",
                    "ITGA1",
                    "ITGAE",
                    "IL2",
                    "IL10",
                    "CXCR6",
                    "CXCL13",
                    "KCNK5",
                    "RGS1",
                    "CRTAM",
                    "DUSP6",
                    "PDCD1",
                    "IL23R"
)


PMN_state_gene <- c("PRF1",
"GZMB",
"GZMA",
"GZMH",
"NKG7",
"GNLY")

T_exhausted <-c("PDCD1","LAYN","HAVCR2","LAG3","CTLA4","TIGIT","TOX","RBPJ","VCAM1","GZMB","MYO7A","CD244","VSIR","BTLA","ENTPD1","CD160","LAIR1")
cytotoxicity <-c("GZMA","GZMB","GZMH","GZMK","GZMH","GNLY","PRF1","IFNG","TNF","SERPINB1",
"SERPINB6","SERPINB9","CTSA","CTSB","CTSC","CTSD","CTSW","CST3","CST7",	
"CSTB","LAMP1","LAMP3","CAPN2","KLRK1","KLRK1","KLRB1","NKG7")
Mediated_Immune_Response<-c("HLA-A","HLA-DRB1","HLA-DRB3","MR1","HMGB1","HSPD1","FBXO38","SLC22A13")
Inhibitory_Receptors <-c("PDCD1",
"CTLA4",
"HAVCR2",
"LAG3",
"TIGIT")



PMN_state_gene<-cytotoxicity
PMN_state_gene<-Mediated_Immune_Response
PMN_state_gene<-c("CA10","ITGA1","ITGAE","IL2","IL10","CXCR6","CXCL13","KCNK5","RGS1","CRTAM","DUSP6","PDCD1",
                  "IL23R",
                  "PDCD1","LAYN","HAVCR2","LAG3","CTLA4","TIGIT","TOX","RBPJ","VCAM1","GZMB","MYO7A","CD244","VSIR","BTLA","ENTPD1","CD160","LAIR1",
                  "GZMA","GZMB","GZMH","GZMK","GZMH","GNLY","PRF1","IFNG","TNF","SERPINB1",
                  "SERPINB6","SERPINB9","CTSA","CTSB","CTSC","CTSD","CTSW","CST3","CST7",	
                  "CSTB","LAMP1","LAMP3","CAPN2","KLRK1","KLRK1","KLRB1","NKG7",
                  "HLA-A","HLA-DRB1","MR1","HMGB1","HSPD1","FBXO38","SLC22A13",
                  "PDCD1","CTLA4","HAVCR2","LAG3","TIGIT")



#计算炎症分数，并添加到seurat对象。
mouse_sc <- AddModuleScore(mouse_sc,
                           features=list('ISC' = PMN_state_gene),
                           pool = rownames(mouse_sc), k=F, nbin=24,
                           name = c('ISC'))

#计算随着拟时变化，各个拟时区间gene score的平均表达
library(dplyr)
modules <- select(mouse_sc@meta.data, c(pseudotime_bin, ISC1))#只选择需要的两列即可
# 确保元数据是 data.frame
meta_data <- as.data.frame(mouse_sc@meta.data)

# 使用 dplyr::select 提取所需列
modules <- dplyr::select(meta_data, pseudotime_bin, ISC1)


modules$pseudotime_bin_num <- as.numeric(modules$pseudotime_bin)
features <- c('ISC1')

tmp <- lapply(1:max(modules$pseudotime_bin_num), function(i){
  df <- modules %>% subset(pseudotime_bin_num == i)
  data.frame(
    value = as.numeric(sum(df[,'ISC1']) / nrow(df)),
    bin_num = i,
    feature = features)
}) #计算出来是list

plot_df <- Reduce(rbind, tmp)#list转为dataframe
plot_df <- na.omit(plot_df)#去除NA

#作图拟合曲线
library(ggplot2)
p_cytotoxicity<-ggplot(plot_df, aes(bin_num, value)) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'gray', alpha = 0.75) +
  geom_point(size = 2) +
  geom_smooth() +
  xlab('pseudotime') + ylab('module score') +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) + 
  ggtitle('ISC along the pseudotime')


p_cytotoxicity

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/组织驻留T细胞/monocle/毒力细胞得分随着拟时间变化.pdf",p_cytotoxicity,width = 5,height = 4)

p_Immune_Response<-ggplot(plot_df, aes(bin_num, value)) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'gray', alpha = 0.75) +
  geom_point(size = 2) +
  geom_smooth() +
  xlab('pseudotime') + ylab('module score') +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) + 
  ggtitle('ISC along the pseudotime')


p_Immune_Response

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/组织驻留T细胞/monocle/免疫细胞得分随着拟时间变化.pdf",p_Immune_Response,width = 5,height = 4)


#从这个例子可以看出，随着拟时的发展，ISC score是明显增加的，而结合之前的相关分析
#male的拟时相关是显著的。所以可说，是不是male的PMN随着拟时。。。。。。有什么生物学意义
#随后，我们可以看看两组之间的相关基因表达差异。

#至于具体的生物学意义就需要各位的想象力了！！！


VlnPlot(mouse_sc, features = 'CA10', group.by = 'MPR_Response2', pt.size = 1)


####Monocle3个性化可视化3---基因表达时序热图

#monocle2中，就有基因拟时热图，monocle3我们依然可以这样呈现。
#前面通过graph_test函数已经得到拟时各阶段基因变化
#我们筛选拟时显著的基因进行可视化，至于选择什么，选择权在你手上
mycds1_res <- graph_test(mycds1, 
                         neighbor_graph="principal_graph", cores=4)#拟时差异基因
pseudotime_table <- subset(mycds1_res, q_value < 0.001)#显著的基因
rownames(pseudotime_table) <- pseudotime_table$gene_short_name
pseudotime_genes <- as.character(rownames(pseudotime_table))

# 拟时区间基因平均表达
#将拟时加入到seurat对象
mouse_sec <- mouse_data[, rownames(mycds_order@colData)]
mouse_sec$pseudotime <- mycds_order@colData$pseudotime
mouse_sec$pseudotime_bin <- mycds_order@colData$pseudotime_bin
#计算拟时区间基因平均表达量
Idents(mouse_sec) <- mouse_sec$pseudotime_bin
average_exp <- AverageExpression(mouse_sec, slot='data', assay='RNA', features=pseudotime_genes)
average_exp <- average_exp$RNA
# 数据缩放
library(scales)
library(dplyr)
scaled_exp <- sapply(1:nrow(average_exp), function(i) rescale(as.numeric(average_exp[i,]), to=c(-1,1))) %>% t %>% as.data.frame
rownames(scaled_exp) <- rownames(average_exp)
average_exp <- scaled_exp
write.csv(average_exp, file = 'D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/组织驻留T细胞/monocle/pseudotime_bin_average_exp.csv')

average_exp<- read.csv('D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/组织驻留T细胞/monocle/pseudotime_bin_average_exp.csv',header = T,row.names = 1)

#如果直接用平均值缩放的数据进行热土绘制，没有拟时的感觉，很差，所以对数据进行处理一下
library(future.apply)

range01 <- function(x){
  cur <- average_exp[x,]
  (cur-min(cur))/(max(cur)-min(cur))
}#这个函数的意思就是用每个点的表达减去最小值后除以极差，这样就获得了高低表达

scaled <- lapply(rownames(average_exp), range01)
scaled <- do.call(rbind, scaled)
rownames(scaled) <- rownames(average_exp)

# 对基因表达排序，按照拟时，从小到大
ordering <- future_lapply(1:nrow(scaled), function(i){
  match(names(scaled[i,])[scaled[i,] >= 0.99][1], colnames(scaled))
})
ordered <- scaled[rownames(scaled)[order(unlist(ordering))],]
ordered <- ordered[rownames(ordered),]

#现在做热图就有那个效果了
library(viridis)
library(ComplexHeatmap)
Heatmap(as.matrix(ordered),
        show_column_names = F, 
        show_row_names=F,
        col = viridis(256),
        cluster_rows=FALSE,
        cluster_columns=FALSE,
        use_raster = TRUE)

#----------------------------------
#我们可以可视化一些主要的基因，至于选择什么基因，按照自己实际情况
grouped <- scaled
grouped$group <- unlist(ordering)#pseudotime_bin group
table(grouped$group)#我这个数据有80个group
grouped <- grouped[rownames(grouped)[order(unlist(ordering))],]

#group分的太细了，我们归纳一下，归纳按照自己的要求处理即可
#我这里粗糙的归为4类
grouped$Group <- 1
grouped["Group"][grouped["group"] > 20] <- 2
grouped["Group"][grouped["group"] > 40] <- 3
grouped["Group"][grouped["group"] > 60] <- 4
module = grouped$Group
module_col = c("1" = "#009E73", "2" = "#F0E442", "3" = "#E69F00", "4" = "#D55E00")
write.csv(grouped, 'grouped.csv')

#选择需要展示的基因
fit_m_genes <- c("CA10",
                    "ITGA1",
                    "ITGAE",
                    "IL2",
                    "IL10",
                    "CXCR6",
                    "CXCL13",
                    "KCNK5",
                    "RGS1",
                    "CRTAM",
                    "DUSP6",
                    "PDCD1",
                    "IL23R")

T_exhausted <-c("PDCD1","LAYN","HAVCR2","LAG3","CTLA4","TIGIT","TOX","RBPJ","VCAM1","GZMB","MYO7A","CD244","VSIR","BTLA","ENTPD1","CD160","LAIR1")
cytotoxicity <-c("GZMA","GZMB","GZMH","GZMK","GZMH","GNLY","PRF1","IFNG","TNF","SERPINB1",
                 "SERPINB6","SERPINB9","CTSA","CTSB","CTSC","CTSD","CTSW","CST3","CST7",	
                 "CSTB","LAMP1","LAMP3","CAPN2","KLRK1","KLRK1","KLRB1","NKG7")
Mediated_Immune_Response<-c("HLA-A","HLA-DRB1","HLA-DRB3","MR1","HMGB1","HSPD1","FBXO38","SLC22A13")
Inhibitory_Receptors <-c("PDCD1",
                         "CTLA4",
                         "HAVCR2",
                         "LAG3",
                         "TIGIT")
fit_m_genes<-c("CA10","ITGA1","ITGAE","IL2","IL10","CXCR6","CXCL13","KCNK5","RGS1","CRTAM","DUSP6",
                  "IL23R",
                  "GZMA","GZMB","GZMH","GZMK","GZMH","GNLY","PRF1","IFNG","TNF","SERPINB1",
                  "SERPINB6","SERPINB9","CTSA","CTSB","CTSC","CTSD","CTSW","CST3","CST7",	
                  "CSTB","LAMP1","LAMP3","CAPN2","KLRK1","KLRK1","KLRB1","NKG7",
                  "PDCD1","CTLA4","HAVCR2","LAG3","TIGIT")


fit_m_genes <- as.data.frame(fit_m_genes)
pdf('D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/组织驻留T细胞/monocle/pseudotime_hetamap2.pdf', width=4, height=8)
ht_list <-Heatmap(module, col = module_col,
                  cluster_rows = FALSE,
                  show_row_names = FALSE,
                  show_column_names = FALSE,
                  row_title = NULL,
                  column_title = NULL,
                  heatmap_legend_param = list(title = "Module"))+
  Heatmap(
    as.matrix(ordered),
    show_column_names = FALSE, 
    show_row_names=FALSE,
    col = viridis(256),
    cluster_rows=FALSE,
    cluster_columns=FALSE,
    use_raster = TRUE,
    heatmap_legend_param = list(title = "%Max")
  )+
  rowAnnotation(link = anno_mark(at = which(rownames(ordered) %in% fit_m_genes$fit_m_genes), 
                                 labels = fit_m_genes$fit_m_genes, labels_gp = gpar(fontsize = 8)))


draw(ht_list, row_km = 4, row_split = module,
     column_title = "Heatmap of pseudotime DEGs", 
     column_title_gp = gpar(fontsize = 12, fontface = "bold"), 
     merge_legends = TRUE, heatmap_legend_side = "right") 
dev.off()

##########################################################################
####拟时连续热图----
load('D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/组织驻留T细胞/monocle/Monocle3_cds_peu.RData')

library(pheatmap)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(monocle3)
#mycds1_res <- graph_test(mycds1, neighbor_graph="principal_graph", cores=4)

mycds1_res<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/组织驻留T细胞/monocle/mycds1_res.csv",header = T,row.names = 1)
genes <- row.names(subset(mycds1_res, q_value< 0.01 & morans_I > 0.2))

plot_matrix <- exprs(mycds1)[match(genes,#挑选我们前面确定需要展示的基因
                                   rownames(rowData(mycds1))),
                             order(pseudotime(mycds1))]#按照拟时对基因排序

#数据拟合平滑处理smooth.spline。计算Z-score
plot_matrix <- t(apply(plot_matrix,1,function(x){smooth.spline(x,df=3)$y}))
plot_matrix <- t(apply(plot_matrix,1,function(x){(x-mean(x))/sd(x)}))
rownames(plot_matrix) <- genes;
dim(plot_matrix)

plot_matrix_combin <- list()
for (i in 1:length(seq(1, 6025-50, 50))){
  num <- seq(1, 6025-50, 50)
  A <- plot_matrix[,num[i]:(50+num[i]-1)]
  a <- rowMeans(A)
  a <- as.data.frame(a)
  a <- a$a
  plot_matrix_combin[[i]] <- a
}

length(plot_matrix_combin)
tailed <- as.data.frame(rowMeans(plot_matrix[,6001:6025]))
plot_matrix_combin[[121]] <- tailed[,1]

plot_matrix_combin <- do.call(cbind, plot_matrix_combin)
rownames(plot_matrix_combin) <- rownames(plot_matrix)

cutColumn_Means <- function(data_exp,#需要分割数据
                            cut#需要分割的列数
){
  plot_matrix_combin <- list()
  nums <- ncol(data_exp)/cut
  if (nums-round(nums, 0)==0){
    
    for (i in 1:length(seq(1, ncol(data_exp), cut))){
      num <- seq(1, ncol(data_exp), cut)
      A <- as.data.frame(rowMeans(data_exp[,num[i]:(cut+num[i]-1)]))[,1]
      plot_matrix_combin[[i]] <- A
      
    }
    plot_matrix_combin <- do.call(cbind, plot_matrix_combin)
    rownames(plot_matrix_combin) <- rownames(data_exp)
    colnames(plot_matrix_combin) <- seq(1,ncol(plot_matrix_combin),1)
    return(plot_matrix_combin)
    
  }else{
    
    for (i in 1:length(seq(1, ncol(data_exp)-cut, cut))){
      num <- seq(1, ncol(data_exp)-cut, cut)
      A <- as.data.frame(rowMeans(data_exp[,num[i]:(cut+num[i]-1)]))[,1]
      plot_matrix_combin[[i]] <- A
    }
    
    plot_matrix_combin[[length(seq(1, ncol(data_exp)-cut, cut))+1]] <- as.data.frame(rowMeans(data_exp[,(max(num)+cut):ncol(data_exp)]))                       
    plot_matrix_combin <- do.call(cbind, plot_matrix_combin)
    rownames(plot_matrix_combin) <- rownames(data_exp)
    colnames(plot_matrix_combin) <- seq(1,ncol(plot_matrix_combin),1)
    return(plot_matrix_combin)
  }
  
  
}


plot_test <- cutColumn_Means(plot_matrix,cut = 25)


#得到矩阵就可以做热图了，很简单
##排序设置
callback = function(hc, mat){
  sv = svd(t(mat))$v[,1]
  dend = reorder(as.dendrogram(hc), wts = sv)
  as.hclust(dend)
}


p1 <- pheatmap::pheatmap(plot_matrix_combin, 
                         useRaster = T,
                         cluster_cols=FALSE, 
                         cluster_rows=T, 
                         show_rownames=F, 
                         show_colnames=F, 
                         clustering_method = "ward.D2",
                         cutree_rows=4,
                         filename=NA,
                         border_color = NA,
                         fontsize_row = 8,
                         color=colorRampPalette(c('#1A5592','white',"#B83D3D"))(100),
                         clustering_callback = callback)


p2 <- pheatmap::pheatmap(plot_test, 
                         useRaster = T,
                         cluster_cols=FALSE, 
                        cluster_rows=T, 
                         show_rownames=F, 
                         show_colnames=F, 
                         clustering_method = "ward.D2",
                         cutree_rows=4,
                         filename=NA,
                         border_color = NA,
                         fontsize_row = 8,
                         color=colorRampPalette(c('#1A5592','white',"#B83D3D"))(100),
                         clustering_callback = callback)


#行注释
annotation_row <- data.frame(Cluster=factor(cutree(p1$tree_row, 4)))
row.names(annotation_row) <- rownames(plot_matrix_combin)

rowcolor <- c("#85B22E","#E29827","#922927",'#57C3F3') 
names(rowcolor) <- c("1","2","3","4") #类型颜色
#注释颜色设置
ann_colors <- list(Cluster=rowcolor) #颜色设置

p3 <- pheatmap::pheatmap(plot_matrix_combin, 
                         cluster_cols=FALSE, 
                         cluster_rows=T, 
                         show_rownames=F, 
                         show_colnames=F, 
                         clustering_method = "ward.D2",
                         filename=NA,
                         cutree_rows=4,
                         border_color = NA,
                         fontsize_row = 8,
                         color=colorRampPalette(c('#1A5592','white',"#B83D3D"))(100),
                         annotation_colors=ann_colors,
                         annotation_row = annotation_row,
                         clustering_callback = callback,
                         annotation_names_col = F,
                         annotation_names_row = F,
                         main="Pseudotime")

#展示需要的基因
gene <- c("CA10",
          "ITGA1",
          "ITGAE",
          "IL2",
          "IL10",
          "CXCR6",
          "CXCL13",
          "KCNK5",
          "RGS1",
          "CRTAM",
          "DUSP6",
          "PDCD1",
          "IL23R")
#这里我们使用add.flag.R函数添加，需要注意的是添加flag要成功，不要添加cutree_rows参数
p4 <- pheatmap::pheatmap(plot_matrix_combin, 
                         cluster_cols=FALSE, 
                         cluster_rows=T, 
                         show_rownames=T, 
                         show_colnames=F, 
                         clustering_method = "ward.D2",
                         filename=NA,
                         border_color = NA,
                         fontsize_row = 8,
                         color=colorRampPalette(c('#1A5592','white',"#B83D3D"))(100),
                         annotation_colors=ann_colors,
                         #annotation_row = annotation_row,
                         clustering_callback = callback,
                         annotation_names_col = F,
                         #annotation_names_row = F,
                         main="Pseudotime")
ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/组织驻留T细胞/monocle/连续基因分析.pdf",p4,width=6,height=7)

source('F:/生信分析代码/KS-10.5合集/KS-account-codes/1-单细胞转录组（scRNA）/2-scRNA 拟时分析/Monocle3/3-monocle3_GO富集分析/add.flag.R')

add.flag(p4,kept.labels = gene,repel.degree = 0.2)

####拟时热图不同module基因GO分析

library(clusterProfiler)
library(org.Mm.eg.db)

###首先提取热图中各个module的基因
module_gene <- as.data.frame(cutree(p3$tree_row, k=4))
colnames(module_gene) <- "Module"
module_gene$gene <- rownames(module_gene)
module_gene$gene <- paste0(toupper(substr(module_gene$gene, 1, 1)), tolower(substr(module_gene$gene, 2, nchar(module_gene$gene))))
#我们这里展示GO结果
Module_GO=data.frame()

for (i in unique(module_gene$Module)) {
  
  data=filter(module_gene,module_gene$Module==i)
  df=bitr(data$gene, 
          fromType="SYMBOL",
          toType=c("ENTREZID"), 
          OrgDb="org.Mm.eg.db")#Symbol转化为ID
  
  go <- enrichGO(gene= unique(df$ENTREZID),
                 OrgDb= org.Mm.eg.db,
                 keyType= 'ENTREZID',
                 ont= "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff= 0.05,
                 qvalueCutoff= 0.05,
                 readable= TRUE)
  go_res=go@result
  
  if (dim(go_res)[1] != 0) {
    go_res$cluster=i
    Module_GO=rbind(Module_GO,go_res)
  }
}

#筛选显著的Terms
Module_GO <- Module_GO[which(Module_GO$qvalue <= 0.05),]
Module_GO <- Module_GO[,c("ID","Description","qvalue","cluster")]
write.csv(Module_GO, file = 'D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/组织驻留T细胞/monocle/Module_GO.csv')
#挑选需要的Terms，使用AI添加到热图上即可
#当然，也可以用ggplot做前面的热图，然后与GO作图合并，但我觉得很麻烦，没有必要。还是AI修饰效果显著
View(Module_GO)



#=======================================================================================
#                             3、monocle2基因随拟时变化趋势
#=======================================================================================
#这个就更简单了，和我们之前讲过的monocle2的趋势一样
#提取拟时基因表达，这里我们用标准化的数据
#前面数据的处理和做热图一样，只不过不需要缩放计算Z-score了
#这里我们提取log标准化的表达数据
peuGene_exp <- normalized_counts(mycds2, norm_method = "log")
peuGene_exp <- as.data.frame(peuGene_exp)
dim(peuGene_exp)
# [1] 11199  6025

pseudotime <- as.data.frame(pseudotime(mycds2))#提取pseudotime
colnames(pseudotime) <- "pseudotime"
#去除Inf
pseudotime <- pseudotime[rowSums(abs(pseudotime)==Inf)==0,,drop=F]
peuGene_exp <- peuGene_exp[,colnames(peuGene_exp) %in% rownames(pseudotime)]
dim(peuGene_exp)
# [1] 11199  5951

#因为我的数据是包含两个分组的，细胞名前面在构建seurat的时候就添加了group，这里分开
peuGene_exp_Female <- peuGene_exp[,grep(pattern="^Female",colnames(peuGene_exp))]
pseudotime_Female <- pseudotime[colnames(peuGene_exp_Female),]

peuGene_exp_Male <- peuGene_exp[,grep(pattern="^Male",colnames(peuGene_exp))]
pseudotime_Male <- pseudotime[colnames(peuGene_exp_Male),]

#这里以一个基因为例，我们做一下拟时的趋势图
Female_anxal1 <- peuGene_exp_Female['CXCL13',]
Female_anxal1 <- as.data.frame(t(Female_anxal1))
Female_anxal1$pseudotime <- pseudotime_Female
Female_anxal1$group <- 'Female'

Male_anxal1 <- peuGene_exp_Male['CXCL13',]
Male_anxal1 <- as.data.frame(t(Male_anxal1))
Male_anxal1$pseudotime <- pseudotime_Male
Male_anxal1$group <- 'Male'

peu_trend <- rbind(Female_anxal1, Male_anxal1)
#ggplot作图
ggplot(peu_trend, aes(x=pseudotime, y=Anxa1, color=group))+
  geom_smooth(aes(fill=group))+ #平滑的填充
  xlab('pseudotime') + 
  ylab('Relative expression') +
  ggtitle('Anxa1')+
  theme_classic(base_size = 12)+
  theme(axis.text = element_text(color = 'black',size = 12),
        axis.title = element_text(color = 'black',size = 14))+
  scale_color_manual(name=NULL, values = c("#089E86","#3D5387"))+#修改颜色
  scale_fill_manual(name=NULL, values = c("#089E86","#3D5387"))#修改颜色







####################################################

####CytoTRACE推断细胞干性----
devtools::install_local("C:/Users/Wangy/Downloads/CytoTRACE_0.3.3.tar.gz")

devtools::install_local("/home/zqwangyansu/CytoTRACE_0.3.3.tar.gz")

devtools::install_github("digitalcytometry/cytotrace", subdir = "cytotrace_r",force=TRUE) 
library(CytoTRACE2)
CD8_res_object<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/CD8_res_object.rds")

table(is.na(CD8_object$cd4_cd8_group))

devtools::install_github("digitalcytometry/cytotrace2", subdir = "cytotrace2_r") 
library(CytoTRACE2)


CD8_object<-readRDS("/home/zqwangyansu/oscc_data/T_cells_CD8_object.rds")
options(future.globals.maxSize = 1024^3 * 50)
cytotrace2_result_sce <- cytotrace2(CD8_res_object, 
                                    is_seurat = TRUE, 
                                    slot_type = "counts", 
                                    species = 'human',
                                    seed = 1234)

annotation <- data.frame(phenotype = CD8_object@meta.data$cd4_cd8_group) %>% 
  set_rownames(., colnames(CD8_object))

# plotting
plots <- plotData(cytotrace2_result = cytotrace2_result_sce, 
                  annotation = annotation, 
                  is_seurat = TRUE)


# 绘制CytoTRACE2_Potency的umap图
umap_raw <- as.data.frame(CD8_res_object@reductions$umap.harmony@cell.embeddings)
umap_raw <- as.data.frame(CD8_res_object@reductions$umap.harmony@cell.embeddings)
umap_raw$umapharmony_1
plot_names <- c("CytoTRACE2_UMAP", "CytoTRACE2_Potency_UMAP", 
                "CytoTRACE2_Relative_UMAP", "Phenotype_UMAP" )

for (plot_name in plot_names) {
  if (!is.null(plots[[plot_name]][[1]]$data)) {
    plots[[plot_name]][[1]]$data$umap_1 <- umap_raw$umapharmony_1
    plots[[plot_name]][[1]]$data$umap_2 <- umap_raw$umapharmony_2
  }
}

x_limits <- range(umap_raw$umapharmony_1, na.rm = TRUE)
y_limits <- range(umap_raw$umapharmony_2, na.rm = TRUE)


plots$CytoTRACE2_UMAP[[1]] <- plots$CytoTRACE2_UMAP[[1]] + 
  scale_x_continuous(limits = x_limits) + 
  scale_y_continuous(limits = y_limits)
plots$CytoTRACE2_UMAP

plots$CytoTRACE2_Potency_UMAP[[1]] <- plots$CytoTRACE2_Potency_UMAP[[1]] + 
  coord_cartesian(xlim = x_limits, ylim = y_limits)
plots$CytoTRACE2_Potency_UMAP

plots$CytoTRACE2_Relative_UMAP[[1]] <- plots$CytoTRACE2_Relative_UMAP[[1]] + 
  coord_cartesian(xlim = x_limits, ylim = y_limits)
plots$CytoTRACE2_Relative_UMAP

plots$Phenotype_UMAP[[1]] <- plots$Phenotype_UMAP[[1]] + 
  coord_cartesian(xlim = x_limits, ylim = y_limits)
plots$Phenotype_UMAP
Idents(cytotrace2_result_sce) <- ""

p1<-FeaturePlot(cytotrace2_result_sce, "CytoTRACE2_Relative",pt.size = 1.5,
                reduction = "umap.harmony") + 
  scale_colour_gradientn(colours = 
                           (c("#9E0142", "#F46D43", "#FEE08B", "#E6F598", 
                              "#66C2A5", "#5E4FA2")), 
                         na.value = "transparent", 
                         limits = c(0, 1), 
                         breaks = seq(0, 1, by = 0.2), 
                         labels = c("0.0 (More diff.)", 
                                    "0.2", "0.4", "0.6", "0.8", "1.0 (Less diff.)"), 
                         name = "Relative\norder \n", 
                         guide = guide_colorbar(frame.colour = "black", 
                                                ticks.colour = "black")) + 
  ggtitle("CytoTRACE 2") + 
  xlab("UMAP1") + ylab("UMAP2") + 
  theme(legend.text = element_text(size = 10), 
        legend.title = element_text(size = 12), 
        axis.text = element_text(size = 13), 
        axis.title = element_text(size = 12), 
        plot.title = element_text(size = 12, 
                                  face = "bold", hjust = 0.5, 
                                  margin = margin(b = 20)),
        axis.line = element_line(linewidth = 0.7)
  ) + 
  theme(aspect.ratio = 1)

p1

class_avg <- data.frame(
  umapharmony_1  = c(0.70586, -5.1612144),
  umapharmony_2  = c(1.24, -1.694),
  cd4_cd8_group = c("CD8_C3-IL7R", "CD8_C4-GZMK")
)

p2<-p1+ 
  ggrepel::geom_label_repel(data = class_avg, 
                            aes(x = umapharmony_1, y = umapharmony_2, label = cd4_cd8_group),
                            size = 5, 
                            box.padding = 0.5, 
                            point.padding = 0.5, 
                            segment.color = 'grey50')+theme(
                              legend.position = "right",        # 图例放置在顶部
                            )

p2

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/组织驻留T细胞/细胞干性推断.pdf",p2,width=7,height = 5)

library(ggpubr)
c("#83B5D3","#AEE94F")
p1 <- ggboxplot(cytotrace2_result_sce@meta.data, x="cd4_cd8_group", y="CytoTRACE2_Score", width = 0.6, 
                color = "cd4_cd8_group",#轮廓颜色
                fill="cd4_cd8_group",#填充
                palette = c("#83B5D3","#AEE94F"),
                xlab = F, #不显示x轴的标签
                bxp.errorbar=T,#显示误差条
                bxp.errorbar.width=0.5, #误差条大小
                size=1, #箱型图边线的粗细
                outlier.shape=NA, #不显示outlier
                legend = "none") #图例放右边 


p1<-p1+ theme(axis.title.x = element_blank(),        # 去掉 X 轴标题
              axis.text = element_text(size = 13), 
              axis.title = element_text(size = 12), 
              plot.title = element_text(size = 12, 
                                        face = "bold", hjust = 0.5, 
                                        margin = margin(b = 20)),
              axis.line = element_line(linewidth = 0.7)
)
###指定组比较
#my_comparisons <- list(c("Epi", "un"), c("T", "un"),c("Myeloid", "un"))
cytotrace2_result_sce@meta.data$cd4_cd8_group
my_comparisons <- list(c("CD8_C1-GNLY", "CD8_C5-STMN1"))
p2<-p1+stat_compare_means(comparisons = my_comparisons,
                          method = "wilcox.test")

ggsave("细胞干性推断boxplot.pdf",p2,width=4,height = 5)





####组织驻留T细胞通讯分析----

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

unique(merged_obj@meta.data$merged_group)

unique(merged_obj$MPR_Response2)

### CD8_C3-IL7R
### CD8_C4-GZMK

merged_obj_filter<- subset(merged_obj,subset=merged_group %in% c("CD8_C4-GZMK","CD8_C2-CTSW",                     
                                                                "CD4_C3-FOXP3","NK","CD3+CD4-CD8-",                    
                                                                 "CD4_C2-IL7R","CD8_C3-IL7R","CD8_C1-GNLY",                     
                                                                 "CD8_C5-STMN1","CD4_C4-FOSB","CD4_C1-PDCD1"))


pre_nmpr_myeloid_t_chat<- subset(merged_obj_filter, subset=MPR_Response2=="Pre_NMPR")
unique(pre_nmpr_myeloid_t_chat$MPR_Response2)
pre_mpr_myeloid_t_chat<- subset(merged_obj_filter, subset=MPR_Response2=="Pre_MPR")
unique(pre_mpr_myeloid_t_chat$MPR_Response2)
post_mpr_myeloid_t_chat<- subset(merged_obj_filter, subset=MPR_Response2=="Post_MPR")
unique(post_mpr_myeloid_t_chat$MPR_Response2)
unique(post_mpr_myeloid_t_chat$CRNCR_Response2)

pre_nmpr.cellchat <- KS_cellchat(pre_nmpr_myeloid_t_chat, assay = 'RNA', group.by = "merged_group",
                                 workers=2, species='human')

saveRDS(pre_nmpr.cellchat, file = "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/细胞通讯分析/pre_nmpr.cellchatT细胞.rds")



options(future.globals.maxSize = 1e+09) 
pre_mpr.cellchat <- KS_cellchat(pre_mpr_myeloid_t_chat, assay = 'RNA', group.by = "merged_group",
                                workers=2, species='human')

saveRDS(pre_mpr.cellchat , file = "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/细胞通讯分析/pre_mpr.cellchatT细胞.rds")


post_mpr.cellchat <- KS_cellchat(post_mpr_myeloid_t_chat, assay = 'RNA', group.by = "merged_group",
                                 workers=2, species='human')

saveRDS(post_mpr.cellchat, file = "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/细胞通讯分析/post_mpr.cellchatt细胞.rds")



####cellchat多组比较分析----

pre_nmpr.cellchat<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/细胞通讯分析/pre_nmpr.cellchatT细胞.rds")

pre_mpr.cellchat<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/细胞通讯分析/pre_mpr.cellchatt细胞.rds")

post_mpr.cellchat<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/细胞通讯分析/post_mpr.cellchatT细胞.rds")


object.list <- list(Pre_NMPR= pre_nmpr.cellchat, Pre_MPR= pre_mpr.cellchat,Post_MPR=post_mpr.cellchat)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))


mat<-cellchat@net$Pre_NMPR$weight

Pre_NMPR<-pheatmap(as.matrix(mat),
                   cluster_rows=F, 
                   cluster_cols=F, 
                   show_colnames=TRUE, 
                   main="Heatmap of Interaction Strength")

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/中央记忆细胞和耗竭性细胞T细胞/细胞通讯分析/和B细胞通讯分析/Pre_NMPR_heatmap.pdf",Pre_NMPR,width = 7,height = 7)

mat<-cellchat@net$Pre_MPR$weight

Pre_MPR<-pheatmap(as.matrix(mat),
                  cluster_rows=F, 
                  cluster_cols=F, 
                  show_colnames=TRUE, 
                  main="Heatmap of Interaction Strength")

mat<-cellchat@net$Post_MPR$weight

Post_MPR<-pheatmap(as.matrix(mat),
                   cluster_rows=F, 
                   cluster_cols=F, 
                   show_colnames=TRUE, 
                   main="Heatmap of Interaction Strength")



# object.list2 <- list(Post_MPR_CR=post_mpr_cr.cellchat,Post_MPR_NCR=post_mpr_ncr.cellchat_ncr)
# cellchat2 <- mergeCellChat(object.list2, add.names = names(object.list2),cell.prefix = TRUE)


gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2,3))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2,3), measure = "weight")
p1<-gg1 + gg2

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/细胞通讯分析/分组比较柱状图.pdf",p1,width=6,height=4)


gg1 <- netVisual_heatmap(cellchat)
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
pdf(file = "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/细胞通讯分析/分组比较热图.pdf", width = 10, height = 6)
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
ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/细胞通讯分析/combined_plots.pdf", p_combined, width = 10, height = 4) 


?rankNet
treatment_color <- c("#4974a4","#4dae47","#f29600")
gg1 <- rankNet(cellchat, mode = "comparison",comparison = c(1, 2), stacked = T,font.size = 10, sources.use = "CD4_C1-PDCD1",targets.use = c("CD8_C4-GZMK","CD8_C2-CTSW",                     
                                                                                                                                           "CD4_C3-FOXP3","NK","CD3+CD4-CD8-",                    
                                                                                                                                           "CD4_C2-IL7R","CD8_C3-IL7R","CD8_C1-GNLY",                     
                                                                                                                                           "CD8_C5-STMN1","CD4_C4-FOSB"),do.stat = TRUE, color.use = c("#4974a4","#4dae47"))

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
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all", signaling = pathway.union, title = names(object.list)[i], width = 6, height = 12)
object.list[[i+1]] <- netAnalysis_computeCentrality(object.list[[i+1]])
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "all", signaling = pathway.union, title = names(object.list)[i+1], width = 6, height = 12)
object.list[[i+2]] <- netAnalysis_computeCentrality(object.list[[i+2]])
ht3 = netAnalysis_signalingRole_heatmap(object.list[[i+2]], pattern = "all", signaling = pathway.union, title = names(object.list)[i+2], width = 6, height = 12)
draw(ht1 + ht2+ht3, ht_gap = unit(0.5, "cm"))



pbubble<-netVisual_bubble(cellchat, sources.use = "CD4_C1-PDCD1", targets.use =c("CD8_C4-GZMK","CD8_C2-CTSW",                     
                                                                                "CD4_C3-FOXP3","NK","CD3+CD4-CD8-",                    
                                                                                "CD4_C2-IL7R","CD8_C3-IL7R","CD8_C1-GNLY",                     
                                                                                "CD8_C5-STMN1","CD4_C4-FOSB"),comparison = c(1,2,3), angle.x = 45,font.size = 12,color.text =c("#4974a4","#4dae47","#f29600"))
?netVisual_bubble

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/组织驻留T细胞/IL7R_受配体_combined_pbubble.pdf", pbubble, width = 9, height = 9) 


pathways.show <- c("PDCD1")
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show)
par(mfrow = c(1,3), xpd=TRUE)

pathways.show <- c("PD-L1")
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 1, signaling.name = paste(pathways.show, names(object.list)[i]))
}

netVisual_aggregate(object.list[[2]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[2], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[2]))
netVisual_aggregate(object.list[[3]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[3], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[3]))

pathways.show <- c("MHC-I")
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show)
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







####CD4 细胞拟时序分析----

