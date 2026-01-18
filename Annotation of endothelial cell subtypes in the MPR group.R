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

setwd("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/内皮细胞亚群注释")
edo_cells_object<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/edo_cells_object.rds")
setwd("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/内皮细胞亚群注释")

edo_cells_object= readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/edo_cells_object.rds")


####提取myeloid cells subject----
Idents(MPRNMPR_object_miMPRobe)<-"cell_type_new"
levels(Idents(MPRNMPR_object_miMPRobe))
Idents(myeloid_cells_object) <-"RNA_snn_res.0.4"
levels(Idents(myeloid_cells_object))

edo_cells_object=subset(MPRNMPR_object_miMPRobe,idents = c("Endothelial cells"))


table(Idents(edo_cells_object))

colnames(edo_cells_object@meta.data)

pbmc<- edo_cells_object

pbmc <- NormalizeData(pbmc,normalization.method = "LogNormalize",
                      scale.factor = 1e4)
GetAssay(pbmc,assay = "RNA")
pbmc <- FindVariableFeatures(pbmc,
                             selection.method = "vst",nfeatures = 5000)
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

p6 <- DimPlot(pbmc, group.by = "orig.ident",reduction = "tsne", pt.size=0.1)
p6

p7 <- DimPlot(pbmc, group.by = "RNA_snn_res.0.4",reduction = "umap", pt.size=0.1,label = T)

p7

p8 <- DimPlot(pbmc, group.by = "RNA_snn_res.0.4",reduction = "umap", pt.size=0.1)

p8

DefaultAssay(pbmc) <- "RNA"

edo_cells_object <- pbmc
Idents(edo_cells_object) <-"RNA_snn_res.0.4"
table(Idents(edo_cells_object))

markers <- FindAllMarkers(edo_cells_object, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
all.markers =markers %>% dplyr::select(gene, everything()) %>% subset(p_val<0.05)
top20 =all.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)

write.csv(all.markers,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/内皮细胞亚群注释/edo_cells_marker_all.markers.csv")
write.csv(top20,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/内皮细胞亚群注释/edo_cells_marker_top20.csv")

#保存RDS
saveRDS(edo_cells_object,"D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/edo_cells_object.rds")

####
####marker基因小提琴图----

setwd("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/内皮细胞亚群注释")
edo_cells_object<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/edo_cells_object.rds")

cell_marker<- read.csv("marker2.csv")

######删除找不到的基因

nofinder<- c("R2F1", "LY6C1", "LY6A", "NR", "GKN3", "MPR2", "CAR4", "ENFB2", "ENDRB", "LRG","CAR14")

cell_marker_filtered <- cell_marker[!cell_marker$gene_name %in% nofinder, ]

Idents(edo_cells_object)<-"RNA_snn_res.0.4"
levels(Idents(edo_cells_object))
unique(cell_marker_filtered$gene_name)
library(randomcoloR)

my36colors = distinctColorPalette(84)

vln.dat=FetchData(edo_cells_object,c(unique(cell_marker_filtered$gene_name),"RNA_snn_res.0.4"))
colnames(vln.dat)
#宽转长
vln.dat.melt <- reshape2::melt(vln.dat, id.vars = c("RNA_snn_res.0.4"), 
                               measure.vars = unique(cell_marker_filtered$gene_name),
                               variable.name = "gene", 
                               value.name = "Expr") 

vln.dat.melt$gene<- factor(vln.dat.melt$gene,levels = unique(vln.dat.melt$gene))

View(vln.dat.melt)

library(cowplot)
vln.dat.melt$RNA_snn_res.0.4<- factor(vln.dat.melt$RNA_snn_res.0.4,levels = c("0","1","2","3","4","5","6","7","8","9","10","11","12",
                                                                              "13"))

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

ggsave("内皮细胞marker小提琴图RNA_snn_res.0.4.jpg",p1,width = 7,height = 7)

DoHeatmap(edo_cells_object, features = unique(cell_marker_filtered$gene_name), size = 3)

FeaturePlot(object=edo_cells_object,reduction = "tsne.harmony",raster=FALSE,features =c("HSPG2",
                                                                                        "IGFBP3",
                                                                                        "SPRY1"
)) 


VlnPlot(edo_cells_object, features = c("PECAM1",
                                       "CDH5",
                                       "FLT1"
                                       
                                       
                                       
)
        ,raster=FALSE)


####cluster 注释----

old_to_new_groups <- c("0"="ACKR1+Blood ECs",
                       "1"="ACKR1+Blood ECs",
                       "3"="ACKR1+Blood ECs",
                       "4"="ACKR1+Blood ECs", 
                       "7"="ACKR1+Blood ECs", 
                       "11"="ACKR1+Blood ECs",
                       "13"="ACKR1+Blood ECs",
                       "5"="Capillary-arterial ECs",
                       "2" ="Lymphatic ECs",
                       "6"="PTPRC+Blood ECs",
                       "9"="MKI67+Blood ECs",
                       "8" = "SPP1+ECs",
                       "12" ="SPP1+ECs",
                       "10"="SPP1+ECs")


edo_cells_object@meta.data$edo_cells_type <- plyr::mapvalues(edo_cells_object@meta.data$RNA_snn_res.0.4,
                                                        from = names(old_to_new_groups), 
                                                        to = old_to_new_groups)

table(edo_cells_object@meta.data$edo_cells_type)

DimPlot(edo_cells_object, group.by = "edo_cells_type",reduction = "umap.harmony", pt.size=0.1)

saveRDS(edo_cells_object,"D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/edo_cells_object.rds")

Idents(edo_cells_object)<-"edo_cells_type"

markers <- FindAllMarkers(edo_cells_object, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
all.markers =markers %>% dplyr::select(gene, everything()) %>% subset(p_val<0.05)
top20 =all.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)

write.csv(all.markers,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/内皮细胞亚群注释/edo_cells_marker_all.markers_cluster.csv")
write.csv(top20,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/内皮细胞亚群注释/edo_cells_marker_top20_cluster.csv")

####marker热图----
####dotplot数据画热图----
cell_marker<- read.csv("marker画图.csv")
p1 <- DotPlot(edo_cells_object, features =unique(cell_marker$gene_name),
              assay='RNA' )
p1
sub_markers_count<- p1$data

names(sub_markers_count)
sub_markers_count$id
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
        axis.text.x = element_text(size = 16,color = 'black',angle =40,hjust =1,vjust=1),
        legend.title = element_text(size=15, color = "black"),
        legend.text = element_text(size=15,color = "black",angle =0),
        legend.position = "right") +
  scale_y_discrete(position = "left")
p1_heat


ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/内皮细胞亚群注释/细胞亚型分类热图.jpg",p1_heat,height =10,width = 7)

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/内皮细胞亚群注释/细胞亚型分类热图.pdf",p1_heat,height =10,width = 7)

#############################
##############################

cell_marker<- read.csv("marker画图.csv")
markers <-cell_marker$gene_name
markers <- as.data.frame(markers)

vln.dat=FetchData(edo_cells_object,c(unique(cell_marker$gene_name),"edo_cells_type"))
colnames(vln.dat)
#宽转长
vln.dat.melt <- reshape2::melt(vln.dat, id.vars = c("edo_cells_type"), 
                               measure.vars = unique(cell_marker$gene_name),
                               variable.name = "gene", 
                               value.name = "Expr") 
vln.dat.melt_mean <- aggregate(vln.dat.melt$Expr, by = list(vln.dat.melt$edo_cells_type, vln.dat.melt$gene), FUN = mean)

colnames(vln.dat.melt_mean) <- c("edo_cells_type","gene","value")
vln.dat.melt<-vln.dat.melt_mean
aver_dt_df_long<- vln.dat.melt
colnames(aver_dt_df_long)
View(aver_dt_df_long)
write.csv(aver_dt_df_long,"edo_cells_object_filter_marker热图数据.csv")
####热图顶部色块
library(tidyverse)
library(ggnewscale)
library(MetBrewer)
library(patchwork)
library(ggtext)

marker_group_list<-c("Lymphatic ECs","Capillary-arterial ECs","ACKR1+Blood ECs","PTPRC+Blood ECs","MKI67+Blood ECs","SPP1+ECs")
aver_dt_df_long$gene <- factor(aver_dt_df_long$gene,levels =rev(unique(aver_dt_df_long$gene)))
aver_dt_df_long$edo_cells_type <- factor(aver_dt_df_long$edo_cells_type,levels =marker_group_list)


p1_heat<-ggplot(aver_dt_df_long) +
  geom_tile(aes(edo_cells_type,gene,fill=value), colour = "white", linewidth = 0.5)+
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
        axis.text.y = element_text(size = 17,color = 'black',
        ),
        panel.border=element_rect(fill=NA,color="black",size=0.5,linetype="solid"),
        axis.text.x = element_text(size = 17,color = 'black',angle =40,hjust =1,vjust=1),
        legend.title = element_text(size=17, color = "black"),
        legend.text = element_text(size=17,color = "black",angle =0),
        legend.position = "right") +
  scale_y_discrete(position = "left")
p1_heat


ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/内皮细胞亚群注释/内皮细胞亚型分类热图.jpg",p1_heat,height =10,width = 7)

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/内皮细胞亚群注释/内皮细胞亚型分类热图.pdf",p1_heat,height =10,width = 7)

####umap 注释之后画图----

library(tidyverse)
edo_cells.harmonyredu<-edo_cells_object@reductions$umap.harmony@cell.embeddings%>%
  as.data.frame() %>% 
  cbind(cell_type =edo_cells_object@meta.data$edo_cells_type)%>% 
  cbind(Patients = edo_cells_object@meta.data$Patients) %>% 
  cbind(Treatments = edo_cells_object@meta.data$Treatments) %>% 
  cbind(MPR_Response = edo_cells_object@meta.data$MPR_Response) %>%
  cbind(MPR_Response2 =edo_cells_object@meta.data$MPR_Response2) %>%
  cbind(cell_id = edo_cells_object@meta.data$cell_id) %>%
  cbind(orig.ident = edo_cells_object@meta.data$orig.ident)

names(edo_cells.harmonyredu)

colnames(edo_cells.harmonyredu) <- c("umapharmony_1","umapharmony_2","cell_type","Patients",
                                      "Treatments","MPR_Response", "MPR_Response2","cell_id",
                                      "orig.ident")

write.csv(edo_cells.harmonyredu,"CAFs_cells_umap降维数据.csv")

###颜色设置###画小坐标轴####https://cloud.tencent.com/developer/article/1924260
library(randomcoloR)
cell_type_color = distinctColorPalette(6)
cell_type_color<- c("#00bfff","#ff1493","#66CDAA","#FFa07a","#0000FF","#FF00ff")

unique(edo_cells.harmonyredu$cell_type)

edo_cells.harmonyredu$cell_type <- factor(edo_cells.harmonyredu$cell_type,levels= c("Lymphatic ECs","Capillary-arterial ECs","ACKR1+Blood ECs","PTPRC+Blood ECs","MKI67+Blood ECs","SPP1+ECs"))
library(ggplot2)

edo_cell_type <- ggplot(edo_cells.harmonyredu,aes(x= umapharmony_1, y =umapharmony_2,color = cell_type)) +  
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
  geom_segment(aes(x = min(edo_cells.harmonyredu$umapharmony_1) , y = min(edo_cells.harmonyredu$umapharmony_2) ,
                   xend = min(edo_cells.harmonyredu$umapharmony_1) +3, yend = min(edo_cells.harmonyredu$umapharmony_2) ),
               colour = "black", size=1,arrow = arrow(length = unit(0.5,"cm")))+ 
  geom_segment(aes(x = min(edo_cells.harmonyredu$umapharmony_1)  , y = min(edo_cells.harmonyredu$umapharmony_2)  ,
                   xend = min(edo_cells.harmonyredu$umapharmony_1) , yend = min(edo_cells.harmonyredu$umapharmony_2) +3),
               colour = "black", size=1,arrow = arrow(length = unit(0.5,"cm"))) +
  annotate("text", x = min(edo_cells.harmonyredu$umapharmony_1) +1.5, y = min(edo_cells.harmonyredu$umapharmony_2)-1, label = "UMPA_1",
           color="black",size = 4, fontface="bold" ) + 
  annotate("text", x = min(edo_cells.harmonyredu$umapharmony_1)-1, y = min(edo_cells.harmonyredu$umapharmony_2) + 1.5, label = "UMAP_2",
           color="black",size = 4, fontface="bold" ,angle=90)

edo_cell_type

edo_cell_type_med <- edo_cells.harmonyredu %>%
  group_by(cell_type) %>%
  summarise(
    umapharmony_1 = median(umapharmony_1),
    umapharmony_2 = median(umapharmony_2)
  )

library(ggrepel)
edo_cell_type_med

edo_cell_type2<-edo_cell_type +geom_label_repel(aes(label=cell_type),size=4,color="black",fontface="bold",data = edo_cell_type_med,
                                                  point.padding=unit(0.5, "lines"),fill = alpha(c("white"),0.5),
                                                  segment.size=0.5, nudge_x=1, nudge_y = 0.5,direction="y",max.overlaps=50)+
  theme(legend.position = "none")


ggsave('./edo_cell_type.jpg',edo_cell_type2,width =6,height = 6)
ggsave('./edo_cell_type.pdf',edo_cell_type2,width =6,height = 6)







