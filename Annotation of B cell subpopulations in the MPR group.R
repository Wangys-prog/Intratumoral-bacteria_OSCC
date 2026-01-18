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

B_cells_object<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/B_cells_object_combined使用.rds")
B_cells_object = UpdateSeuratObject(B_cells_object)
B_cells_object
####提取Bcell subject----
Idents(MPRNMPR_object_miMPRobe)<-"cell_type_new"

B_cells_object=subset(MPRNMPR_object_miMPRobe,idents = c('B cells',"Plasma cells"))

colnames(B_cells_object@meta.data)
pbmc<- B_cells_object
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

# cell_cluster <- data.frame(cell_ID=rownames(metadata), cluster_ID=metadata$RNA_snn_res.0.8) #将聚类结果另存为data.frame
#
# write.csv(cell_cluster,'D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/tcell_cluster.csv',row.names = F, quote = F)

pbmc <- RunUMAP(pbmc,dims = 1:20)#UMAP二次降维，取前20个维度
pbmc <- RunTSNE(pbmc,dims = 1:20)
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
pbmc <- FindNeighbors(pbmc, reduction = "tsne.harmony", dims = 1:2)
pbmc <- FindClusters(pbmc)#标准聚类

#绘图和保存
p1 <- DimPlot(pbmc, group.by = "orig.ident",reduction = "umap.harmony", pt.size=0.1) +
  ggtitle("Integrated by harmony")#去批次后的可视化
p1
p2 <- DimPlot(pbmc, group.by = "RNA_snn_res.0.4",reduction = "umap.harmony", pt.size=0.1,label = T) +
  ggtitle("Integrated by harmony")#去批次后的可视化
p2

p3 <- DimPlot(pbmc, group.by = "RNA_snn_res.0.4",reduction = "tsne.harmony", pt.size=0.1,label = T) +
  ggtitle("Integrated by harmony")#去批次后的可视化
p3

p4 <- DimPlot(pbmc, group.by = "RNA_snn_res.0.4",reduction = "tsne.harmony", pt.size=0.1,label = T) +
  ggtitle("Integrated by harmony")#去批次后的可视化
p4

p5 <- DimPlot(pbmc, group.by = "orig.ident",reduction = "tsne.harmony", pt.size=0.1)
p5

DefaultAssay(pbmc) <- "RNA"
markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
all.markers =markers %>% dplyr::select(gene, everything()) %>% subset(p_val<0.05)
top20 =all.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
Idents(pbmc) <- pbmc$RNA_snn_res.0.4
p_heatmap<-DoHeatmap(pbmc,features = top20$gene)+ #热图复现
  theme(text = element_text(size = 5))

write.csv(top20,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/B_cell注释/Bcell_marker_top20.csv")
ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/B_cell注释/Bcell_topmarkergene_heatmap.jpg",p_heatmap, width = 10,height = 15)
B_cells_object<-pbmc
#保存RDS
saveRDS(B_cells_object,"D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/B_cells_object.rds")

####marker基因小提琴图----
setwd("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/B_cell注释")
B_cells_object<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/B_cells_object.rds")

Bcell_marker<- read.csv("B_marker.csv")
markers <-Bcell_marker$gene_name
markers <- as.data.frame(markers)
levels(Idents(B_cells_object))
View(B_cells_object)
Idents(B_cells_object)<-"RNA_snn_res.0.4"

b_dimplot<-DimPlot(B_cells_object, group.by = "RNA_snn_res.0.4",reduction = "tsne.harmony", pt.size=0.1,label = T)
ggsave("B_celldimplot.jpg",b_dimplot,width = 8,height = 9)

b_dimplot<-DimPlot(B_cells_object, group.by = "RNA_snn_res.0.4",reduction = "umap.harmony", pt.size=0.1,label = T)
ggsave("B_celldimplot_umpp.jpg",b_dimplot,width = 8,height = 9)


FeaturePlot(object=B_cells_object,reduction = "tsne.harmony",raster=FALSE,features =c("FLT3",
                                                                                      "CD1C"
)) 
library(randomcoloR)

my36colors = distinctColorPalette(47)


my36colors <-c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
               '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
               '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
               '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
               '#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
               '#968175'
)

VlnPlot(object = B_cells_object, features = "IGHA1")

VlnPlot(B_cells_object, features = unique(Bcell_marker$gene_name), 
        stack = TRUE, 
        sort = TRUE, 
        cols = my36colors,
        split.by ="RNA_snn_res.0.4" , #每种cluster 一个颜色
        flip = TRUE) +
  theme(legend.position = "none") + 
  ggtitle("Identity on x-axis")

Idents(B_cells_object)<- "B_cells_type"
markers <- FindAllMarkers(B_cells_object, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
all.markers =markers %>% dplyr::select(gene, everything()) %>% subset(p_val<0.05)
all.markers$cluster
top20 =all.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)

write.csv(top20,"Bcells_marker_top20.csv")

####ggplot2画小提琴图----

vln.dat=FetchData(B_cells_object,c(unique(Bcell_marker$gene_name),"RNA_snn_res.0.4"))
colnames(vln.dat)
#宽转长
vln.dat.melt <- reshape2::melt(vln.dat, id.vars = c("RNA_snn_res.0.4"), 
                               measure.vars = unique(Bcell_marker$gene_name),
                               variable.name = "gene", 
                               value.name = "Expr") 

#vln.dat.melt_mean <- aggregate(vln.dat.melt$Expr, by = list(vln.dat.melt$RNA_snn_res.0.4, vln.dat.melt$gene), FUN = mean)


# ggplot(vln.dat.melt, aes(factor(RNA_snn_res.0.4), Expr, fill = gene)) +
#   geom_violin(scale = "width", adjust = 1, trim = TRUE) +
#   facet_grid(rows = vars(gene), scales = "free", switch = "y")

vln.dat.melt$gene<- factor(vln.dat.melt$gene,levels = unique(vln.dat.melt$gene))
library(cowplot)

p1 <- ggplot(vln.dat.melt, aes(gene, Expr, fill = gene)) +
  geom_violin(scale = "width", adjust = 1, trim = TRUE) +
  scale_y_continuous(expand = c(0, 0), position="right", labels = function(x)
    c(rep(x = "", times = length(x)-2), x[length(x) - 1], "")) +
  facet_grid(rows = vars(RNA_snn_res.0.4), scales = "free", switch = "y") +
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

ggsave("B细胞marker小提琴图.jpg",p1,width = 15,height = 14)


featureplot

FeaturePlot(object =B_cells_object,reduction = "tsne.harmony",raster=FALSE,features =c("CD34")) 


####命名cluster----

old_to_new_groups <- c("12" = "Naive B cells",
                        "3" = "Primarily mature naive B cells",
                       "7" = "Pre-B cells",
                       "13" = "Immature B cells",
                         "1" = "Plasma cells",
                         "4" ="Plasma cells",
                        "6" ="Plasma cells",
                        "10"="Plasma cells",
                        "14"="Plasma cells",
                        "5"= "Naive plasma cell",
                        "0"="Active B cells",
                        "8"="Active B cells",
                        "2"="Active B cells",
                        "9"="Active B cells",
                        "11" ="Active B cells")


B_cells_object@meta.data$B_cells_type <- plyr::mapvalues(B_cells_object@meta.data$RNA_snn_res.0.4,
                                                          from = names(old_to_new_groups), 
                                                          to = old_to_new_groups)

table(B_cells_object@meta.data$B_cells_type)

DimPlot(B_cells_object, group.by = "B_cells_type",reduction = "tsne.harmony", pt.size=0.1)

saveRDS(B_cells_object,"D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/B_cells_object.rds")

####细胞注释正式画图----

setwd("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/B_cell注释")
B_cells_object<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/B_cells_object.rds")

library(tidyverse)
B_cells_object.harmonyredu<- B_cells_object@reductions$tsne.harmony@cell.embeddings%>%
  as.data.frame() %>% 
  cbind(cell_type = B_cells_object@meta.data$B_cells_type)%>% 
  cbind(Patients = B_cells_object@meta.data$Patients) %>% 
  cbind(Treatments = B_cells_object@meta.data$Treatments) %>% 
  cbind(MPR_Response = B_cells_object@meta.data$MPR_Response) %>%
  cbind(MPR_Response2 = B_cells_object@meta.data$MPR_Response2) %>%
  cbind(cell_id = B_cells_object@meta.data$cell_id) %>%
  cbind(orig.ident = B_cells_object@meta.data$orig.ident)

names(B_cells_object.harmonyredu)

colnames(B_cells_object.harmonyredu) <- c("tsneharmony_1","tsneharmony_2","cell_type","Patients",
                                          "Treatments","MPR_Response", "MPR_Response2","cell_id",
                                          "orig.ident")

write.csv(B_cells_object.harmonyredu,"B_cells_tsne降维数据.csv")


###颜色设置###画小坐标轴####https://cloud.tencent.com/developer/article/1924260
library(randomcoloR)
cell_type_color = distinctColorPalette(7)
cell_type_color<- c("#C453DC","#00ffff","#0000FF","#FF6347","#8EDCC1","#7fff00","#ffd700")

B_cells_object.harmonyredu$cell_type <- factor(B_cells_object.harmonyredu$cell_type,levels= c("Naive B cells","Pre-B cells","Primarily mature naive B cells","Immature B cells","Active B cells",
                                                                                              "Naive plasma cell","Plasma cells"))
library(ggplot2)

Bcells_cell_type <- ggplot(B_cells_object.harmonyredu,aes(x= tsneharmony_1, y = tsneharmony_2,color = cell_type)) +  
  geom_point(size = 0.5, alpha =1)+  
  scale_color_manual(values = cell_type_color)+
  xlab("t-SNE1")+
  ylab("t-SNE2")+
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
  geom_segment(aes(x = min(B_cells_object.harmonyredu$tsneharmony_1) , y = min(B_cells_object.harmonyredu$tsneharmony_2) ,
                   xend = min(B_cells_object.harmonyredu$tsneharmony_1) +10, yend = min(B_cells_object.harmonyredu$tsneharmony_2) ),
               colour = "black", size=1,arrow = arrow(length = unit(0.5,"cm")))+ 
  geom_segment(aes(x = min(B_cells_object.harmonyredu$tsneharmony_1)  , y = min(B_cells_object.harmonyredu$tsneharmony_2)  ,
                   xend = min(B_cells_object.harmonyredu$tsneharmony_1) , yend = min(B_cells_object.harmonyredu$tsneharmony_2) +10),
               colour = "black", size=1,arrow = arrow(length = unit(0.5,"cm"))) +
  annotate("text", x = min(B_cells_object.harmonyredu$tsneharmony_1) +3.5, y = min(B_cells_object.harmonyredu$tsneharmony_2)-3, label = "t-SNE_1",
           color="black",size = 4, fontface="bold" ) + 
  annotate("text", x = min(B_cells_object.harmonyredu$tsneharmony_1)-3, y = min(B_cells_object.harmonyredu$tsneharmony_2) + 3.5, label = "t-SNE_2",
           color="black",size = 4, fontface="bold" ,angle=90)

Bcells_cell_type

B_cell_type_med <- B_cells_object.harmonyredu %>%
  group_by(cell_type) %>%
  summarise(
    tsneharmony_1 = median(tsneharmony_1),
    tsneharmony_2 = median(tsneharmony_2)
  )

library(ggrepel)
B_cell_type_med 

Bcells_cell_type2<-Bcells_cell_type +geom_label_repel(aes(label=cell_type),size=4,color="black",fontface="bold",data = B_cell_type_med,
                                                        point.padding=unit(0.5, "lines"),fill = alpha(c("white"),0.5),
                                                        segment.size=0.5, nudge_x=1, nudge_y = 0.5,direction="y",max.overlaps=50)+
  theme(legend.position = "none")


ggsave('./Bcells_cell_type2.jpg',Bcells_cell_type2,width =6,height = 6)
ggsave('./Bcells_cell_type2.pdf',Bcells_cell_type2,width =6,height = 6)

####重新选择marker 画小提琴图----

setwd("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/B_cell注释")
B_cells_object<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/B_cells_object.rds")

Bcell_marker<- read.csv("B_marker.csv")
markers <-Bcell_marker$gene_name
markers <- as.data.frame(markers)
levels(Idents(B_cells_object))
View(B_cells_object)
Idents(B_cells_object)<-"RNA_snn_res.0.4"
unique(Bcell_marker$gene_name)
library(randomcoloR)

my36colors = distinctColorPalette(40)


my36colors <-c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
               '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
               '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
               '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
               '#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
               '#968175'
)


####ggplot2画小提琴图----

vln.dat=FetchData(B_cells_object,c(unique(Bcell_marker$gene_name),"B_cells_type"))
colnames(vln.dat)
#宽转长
vln.dat.melt <- reshape2::melt(vln.dat, id.vars = c("B_cells_type"), 
                               measure.vars = unique(Bcell_marker$gene_name),
                               variable.name = "gene", 
                               value.name = "Expr") 

vln.dat.melt$gene<- factor(vln.dat.melt$gene,levels = unique(vln.dat.melt$gene))
vln.dat.melt$B_cells_type<- factor(vln.dat.melt$B_cells_type,levels = c("Pre-B cells","Immature B cells","Primarily mature naive B cells","Naive B cells","Active B cells",
                                                                        "Naive plasma cell","Plasma cells"))
View(vln.dat.melt)

library(cowplot)

p1 <- ggplot(vln.dat.melt, aes(gene,Expr,fill = gene)) +
  geom_violin(scale = "width", adjust = 2, trim = TRUE) +
  scale_y_continuous(expand = c(0, 0), position="right", labels = function(x)
    c(rep(x = "", times = length(x)-2), x[length(x) - 1], "")) +
  # facet_grid(rows = vars(B_cells_type), scales = "free", switch = "y") +
  facet_grid(rows = vars(B_cells_type),switch = "y")+
  # scale_y_continuous(limits = c(0, 7),breaks = seq(0,6,4))+
  scale_fill_manual(values = my36colors) + 
  theme_cowplot(font_size = 12) +
  theme(legend.position = "none", panel.spacing = unit(0, "lines"),
        plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = NA, color = "black"),
        panel.border = element_rect(linewidth=1),
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

ggsave("B细胞marker小提琴图.jpg",p1,width = 15,height = 14)


#### 细胞注释marker 热图----
Bcell_marker<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/B_cell注释/B_marker.csv")
markers <-Bcell_marker$gene_name
markers <- as.data.frame(markers)

vln.dat=FetchData(B_cells_object,c(unique(Bcell_marker$gene_name),"B_cells_type"))
colnames(vln.dat)
#宽转长
vln.dat.melt <- reshape2::melt(vln.dat, id.vars = c("B_cells_type"), 
                               measure.vars = unique(Bcell_marker$gene_name),
                               variable.name = "gene", 
                               value.name = "Expr") 

vln.dat.melt_mean <- aggregate(vln.dat.melt$Expr, by = list(vln.dat.melt$B_cells_type, vln.dat.melt$gene), FUN = mean)

colnames(vln.dat.melt_mean) <- c("B_cells_type","gene","value")
vln.dat.melt<-vln.dat.melt_mean

aver_dt_df_long<- vln.dat.melt
View(aver_dt_df_long)
write.csv(aver_dt_df_long,"B细胞小类marker画图数据.csv")
####热图顶部色块
library(tidyverse)
library(ggnewscale)
library(MetBrewer)
library(patchwork)
library(ggtext)
marker_group_list<-c("Pre-B cells","Immature B cells","Primarily mature naive B cells","Naive B cells","Active B cells",
                     "Naive plasma cell","Plasma cells")

aver_dt_df_long$gene <- factor(aver_dt_df_long$gene,levels =rev(unique(aver_dt_df_long$gene)))
aver_dt_df_long$Bcell_group <- factor(aver_dt_df_long$B_cells_type,levels =marker_group_list)

p1_heat<-aver_dt_df_long %>% ggplot(aes(Bcell_group,gene)) +
  geom_tile(data=aver_dt_df_long,aes(fill=scale(value)), colour = "white", linewidth = 0.5)+
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
        axis.text.x = element_text(size = 17,color = 'black',angle =40,hjust =1,vjust=1),
        legend.title = element_text(size=15, color = "black"),
        legend.text = element_text(size=15,color = "black",angle =0),
        legend.position = "right") +
  scale_y_discrete(position = "left")
p1_heat


ggsave("Bcell_细胞亚型分类热图.jpg",p1_heat,height = 11,width = 6.5)
ggsave("Bcell_细胞亚型分类热图.pdf",p1_heat,height = 11,width = 6.5)


####dotplot 转heatmap----
cell_marker <- Bcell_marker
p1 <- DotPlot(B_cells_object, features =unique(cell_marker$gene_name),
              assay='RNA' )
p1
sub_markers_count<- p1$data

names(sub_markers_count)

sub_markers_count$id<- factor(sub_markers_count$id, levels = c("Pre-B cells","Immature B cells","Primarily mature naive B cells","Naive B cells","Active B cells",
                                                               "Naive plasma cell","Plasma cells"))

colnames(sub_markers_count)
sub_markers_count<- sub_markers_count %>% data.frame()
sub_markers_count$features.plot <- factor(sub_markers_count$features.plot,levels =rev(unique(sub_markers_count$features.plot)))
sub_markers_count$id <- factor(sub_markers_count$id,levels =c("Pre-B cells","Immature B cells","Primarily mature naive B cells","Naive B cells","Active B cells",
                                                              "Naive plasma cell","Plasma cells"))
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
        axis.text.y = element_text(size = 10,color = 'black',
        ),
        panel.border=element_rect(fill=NA,color="black",size=0.5,linetype="solid"),
        axis.text.x = element_text(size = 10,color = 'black',angle =40,hjust =1,vjust=1),
        legend.title = element_text(size=10, color = "black"),
        legend.text = element_text(size=10,color = "black",angle =0),
        legend.position = "right") +
  scale_y_discrete(position = "left")
p1_heat

# ggsave("myeloid_cells_细胞亚型分类热图.jpg",p1_heat,height = 15,width = 8)
# ggsave("myeloid_cells_细胞亚型分类热图.pdf",p1_heat,height = 11,width = 6.5)


####B 细胞细胞比例----
library(ggplot2)
library(dplyr)
library(ggpubr)
library(tidyverse)
library(ggh4x)
library(rstatix)
library(ggpubr)
cell_types <- B_cells_object@meta.data$B_cells_type %>% as.data.frame()

microbe_group<-cbind(B_cells_object@meta.data$cell_id) %>%
  cbind(B_cells_object@meta.data$orig.ident)%>%
  cbind(cell_type = as.character(B_cells_object@meta.data$B_cells_type))%>% 
  cbind(Patients = B_cells_object@meta.data$Patients) %>% 
  cbind(Treatments = B_cells_object@meta.data$Treatments) %>% 
  cbind(MPR_Response2 = B_cells_object@meta.data$MPR_Response2) %>%as.data.frame()
  
colnames(microbe_group) <- c("cell_id","samples","cell_type","Patients","Treatments","MPR_Response2")


Cellratio <- prop.table(table(Idents(B_cells_object),B_cells_object@meta.data$orig.ident), margin = 2)
Cellratio <- as.data.frame(Cellratio)
colnames(Cellratio) <- c("cell_type","sample_id","Freq")
View(Cellratio)
Cellratio$samples<- microbe_group[match(Cellratio$sample_id,microbe_group$samples),2]
Cellratio$Patients<- microbe_group[match(Cellratio$sample_id,microbe_group$samples),4]
Cellratio$MPR_Response2<- microbe_group[match(Cellratio$sample_id,microbe_group$samples),6]
View(Cellratio)

Cellratio$MPR_Response2<- factor(Cellratio$MPR_Response2,levels=c("Pre_NMPR","Pre_MPR","Post_MPR"))
Cellratio$cell_type<- factor(Cellratio$cell_type,levels=c("Pre-B cells","Immature B cells","Primarily mature naive B cells","Naive B cells","Active B cells",
                                                            "Naive plasma cell","Plasma cells"))


View(Cellratio)

Cellratio_2<- Cellratio%>% filter(grepl("Plasma cells",cell_type)) %>%as.data.frame()
View(Cellratio_2)
##Immature B cells Pre_MPR XL_SC_pre	XYH_SC_pre

####
df_p_val2 <-Cellratio_2 %>% group_by(cell_type)%>%
  wilcox_test(Freq~MPR_Response2) %>%
  adjust_pvalue(p.col="p",method="bonferroni") %>%
  add_significance(p.col="p.adj")%>% 
  add_xy_position(x="MPR_Response2")

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

####合并髓系细胞cluster20----
pbmc.combined <- merge(B_cells_object, y = cluster20, add.cell.ids = c("B cells", "cluster20"), project = "PBMC12K")

pbmc<- pbmc.combined
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
pbmc <- RunTSNE(pbmc,dims = 1:20)
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
pbmc <- FindNeighbors(pbmc, reduction = "tsne.harmony", dims = 1:2)
pbmc <- FindClusters(pbmc)#标准聚类

#绘图和保存
p1 <- DimPlot(pbmc, group.by = "orig.ident",reduction = "umap.harmony", pt.size=0.1) +
  ggtitle("Integrated by harmony")#去批次后的可视化
p1
p2 <- DimPlot(pbmc, group.by = "RNA_snn_res.0.4",reduction = "umap.harmony", pt.size=0.1,label = T) +
  ggtitle("Integrated by harmony")#去批次后的可视化
p2

p3 <- DimPlot(pbmc, group.by = "RNA_snn_res.0.4",reduction = "tsne.harmony", pt.size=0.1,label = T) +
  ggtitle("Integrated by harmony")#去批次后的可视化
p3

p4 <- DimPlot(pbmc, group.by = "RNA_snn_res.0.4",reduction = "tsne.harmony", pt.size=0.1,label = T) +
  ggtitle("Integrated by harmony")#去批次后的可视化
p4

p5 <- DimPlot(pbmc, group.by = "orig.ident",reduction = "tsne.harmony", pt.size=0.1)
p5

DefaultAssay(pbmc) <- "RNA"
markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
all.markers =markers %>% dplyr::select(gene, everything()) %>% subset(p_val<0.05)
top20 =all.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
Idents(pbmc) <- pbmc$RNA_snn_res.0.4
p_heatmap<-DoHeatmap(pbmc,features = top20$gene)+ #热图复现
  theme(text = element_text(size = 5))

write.csv(top20,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/B_cell注释/Bcell_marker_top20.csv")
ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/B_cell注释/Bcell_topmarkergene_heatmap.jpg",p_heatmap, width = 10,height = 15)
B_cells_object<-pbmc
#保存RDS
saveRDS(pbmc.combined,"D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/B_cells_object_combined.rds")


####marker基因小提琴图----
pbmc.combined <-pbmc

setwd("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/B_cell注释")

cell_marker<- read.csv("B_marker2.csv")
levels(Idents(pbmc.combined))
View(pbmc.combined)
Idents(pbmc.combined)<-"RNA_snn_res.0.4"

library(randomcoloR)

my36colors = distinctColorPalette(49)

nofinder<- c( "CD20", "CD138", "IGHA")

cell_marker_filtered <- cell_marker[!cell_marker$gene_name %in% nofinder, ]

vln.dat=FetchData(pbmc.combined,c(unique(cell_marker_filtered$gene_name),"RNA_snn_res.0.4"))
colnames(vln.dat)
#宽转长
vln.dat.melt <- reshape2::melt(vln.dat, id.vars = c("RNA_snn_res.0.4"), 
                               measure.vars = unique(cell_marker_filtered$gene_name),
                               variable.name = "gene", 
                               value.name = "Expr") 

vln.dat.melt$gene<- factor(vln.dat.melt$gene,levels = unique(vln.dat.melt$gene))
library(cowplot)

p1 <- ggplot(vln.dat.melt, aes(gene, Expr, fill = gene)) +
  geom_violin(scale = "width", adjust = 1, trim = TRUE) +
  scale_y_continuous(expand = c(0, 0), position="right", labels = function(x)
    c(rep(x = "", times = length(x)-2), x[length(x) - 1], "")) +
  facet_grid(rows = vars(RNA_snn_res.0.4), scales = "free", switch = "y") +
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

ggsave("B细胞marker小提琴图.jpg",p1,width = 15,height = 14)


featureplot

FeaturePlot(object =B_cells_object,reduction = "tsne.harmony",raster=FALSE,features =c("CD34")) 

####合并之后重新命名----

old_to_new_groups <- c("13" = "Naive B cells",
                       "3" = "Primarily mature naive B cells",
                       "8" = "Pre-B cells",
                       "2" = "Plasma cells",
                       "5" ="Plasma cells",
                       "7" ="Plasma cells",
                       "11"="Plasma cells",
                       "14"="Plasma cells",
                       "15"="Plasma cells",
                       "6"= "Naive plasma cell",
                       "0"="Active B cells",
                       "1"="Active B cells",
                       "4"="Active B cells",
                       "9"="Active B cells",
                       "12" ="Active B cells",
                       "10" ="Active B cells")


pbmc.combined@meta.data$B_cells_type <- plyr::mapvalues(pbmc.combined@meta.data$RNA_snn_res.0.4,
                                                         from = names(old_to_new_groups), 
                                                         to = old_to_new_groups)

table(pbmc.combined@meta.data$B_cells_type)

DimPlot(pbmc.combined, group.by = "B_cells_type",reduction = "tsne.harmony", pt.size=0.1)

saveRDS(pbmc.combined,"D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/B_cells_object_combined.rds")

####合并轴细胞注释正式画图----

setwd("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/B_cell注释")
pbmc.combined<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/B_cells_object_combined.rds")

library(tidyverse)
pbmc.combined.harmonyredu<- pbmc.combined@reductions$tsne.harmony@cell.embeddings%>%
  as.data.frame() %>% 
  cbind(cell_type = pbmc.combined@meta.data$B_cells_type)%>% 
  cbind(Patients = pbmc.combined@meta.data$Patients) %>% 
  cbind(Treatments = pbmc.combined@meta.data$Treatments) %>% 
  cbind(MPR_Response = pbmc.combined@meta.data$MPR_Response) %>%
  cbind(MPR_Response2 = pbmc.combined@meta.data$MPR_Response2) %>%
  cbind(cell_id = pbmc.combined@meta.data$cell_id) %>%
  cbind(orig.ident = B_cells_object@meta.data$orig.ident)

names(pbmc.combined.harmonyredu)

colnames(pbmc.combined.harmonyredu) <- c("tsneharmony_1","tsneharmony_2","cell_type","Patients",
                                          "Treatments","MPR_Response", "MPR_Response2","cell_id",
                                          "orig.ident")

write.csv(pbmc.combined.harmonyredu,"pbmc.combined_B_cells_tsne降维数据.csv")


###颜色设置###画小坐标轴####https://cloud.tencent.com/developer/article/1924260
library(randomcoloR)
cell_type_color = distinctColorPalette(7)
cell_type_color<- c("#C453DC","#00ffff","#0000FF","#8EDCC1","#7fff00","#ffd700")

unique(pbmc.combined.harmonyredu$cell_type)
pbmc.combined.harmonyredu$cell_type <- factor(pbmc.combined.harmonyredu$cell_type,levels= c("Naive B cells","Pre-B cells","Primarily mature naive B cells","Active B cells",
                                                                                              "Naive plasma cell","Plasma cells"))
library(ggplot2)

Bcells_cell_type <- ggplot(pbmc.combined.harmonyredu,aes(x= tsneharmony_1, y = tsneharmony_2,color = cell_type)) +  
  geom_point(size = 0.5, alpha =1)+  
  scale_color_manual(values = cell_type_color)+
  xlab("t-SNE1")+
  ylab("t-SNE2")+
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
        legend.text = element_text(size=15), #设置legend标签的大小
        legend.key.size=unit(1,'cm'),
        legend.position = "None") +  # 设置legend标签之间的大小
  guides(color = guide_legend(override.aes = list(size=5)))+ #设置legend中 点的大小 
  geom_segment(aes(x = min(B_cells_object.harmonyredu$tsneharmony_1) , y = min(B_cells_object.harmonyredu$tsneharmony_2) ,
                   xend = min(B_cells_object.harmonyredu$tsneharmony_1) +10, yend = min(B_cells_object.harmonyredu$tsneharmony_2) ),
               colour = "black", size=1,arrow = arrow(length = unit(0.5,"cm")))+ 
  geom_segment(aes(x = min(B_cells_object.harmonyredu$tsneharmony_1)  , y = min(B_cells_object.harmonyredu$tsneharmony_2)  ,
                   xend = min(B_cells_object.harmonyredu$tsneharmony_1) , yend = min(B_cells_object.harmonyredu$tsneharmony_2) +10),
               colour = "black", size=1,arrow = arrow(length = unit(0.5,"cm"))) +
  annotate("text", x = min(B_cells_object.harmonyredu$tsneharmony_1) +3.5, y = min(B_cells_object.harmonyredu$tsneharmony_2)-3, label = "t-SNE_1",
           color="black",size = 4, fontface="bold" ) + 
  annotate("text", x = min(B_cells_object.harmonyredu$tsneharmony_1)-3, y = min(B_cells_object.harmonyredu$tsneharmony_2) + 3.5, label = "t-SNE_2",
           color="black",size = 4, fontface="bold" ,angle=90)

Bcells_cell_type

B_cell_type_med <- pbmc.combined.harmonyredu %>%
  group_by(cell_type) %>%
  summarise(
    tsneharmony_1 = median(tsneharmony_1),
    tsneharmony_2 = median(tsneharmony_2)
  )

library(ggrepel)
B_cell_type_med 

Bcells_cell_type2<-Bcells_cell_type +geom_label_repel(aes(label=cell_type),size=4,color="black",fontface="bold",data = B_cell_type_med,
                                                      point.padding=unit(0.5, "lines"),fill = alpha(c("white"),0.5),
                                                      segment.size=0.5, nudge_x=1, nudge_y = 0.5,direction="y",max.overlaps=50)+
  theme(legend.position = "none")


ggsave('./Bcells_cell_type合并.jpg',Bcells_cell_type2,width =6,height = 6)
ggsave('./Bcells_cell_type合并.pdf',Bcells_cell_type2,width =6,height = 6)

####合并之后热图----

Bcell_marker<- read.csv("B_marker.csv")
markers <-Bcell_marker$gene_name
markers <- as.data.frame(markers)

vln.dat=FetchData(pbmc.combined,c(unique(Bcell_marker$gene_name),"B_cells_type"))
colnames(vln.dat)
#宽转长
vln.dat.melt <- reshape2::melt(vln.dat, id.vars = c("B_cells_type"), 
                               measure.vars = unique(Bcell_marker$gene_name),
                               variable.name = "gene", 
                               value.name = "Expr") 




vln.dat.melt_mean <- aggregate(vln.dat.melt$Expr, by = list(vln.dat.melt$B_cells_type, vln.dat.melt$gene), FUN = mean)

colnames(vln.dat.melt_mean) <- c("B_cells_type","gene","value")
vln.dat.melt<-vln.dat.melt_mean

aver_dt_df_long<- vln.dat.melt
View(aver_dt_df_long)
write.csv(aver_dt_df_long,"B细胞小类marker画图数据合并.csv")
####热图顶部色块
library(tidyverse)
library(ggnewscale)
library(MetBrewer)
library(patchwork)
library(ggtext)
marker_group_list<-c("Pre-B cells","Primarily mature naive B cells","Naive B cells","Active B cells",
                     "Naive plasma cell","Plasma cells")

aver_dt_df_long$gene <- factor(aver_dt_df_long$gene,levels =rev(unique(aver_dt_df_long$gene)))
aver_dt_df_long$Bcell_group <- factor(aver_dt_df_long$B_cells_type,levels =marker_group_list)

p1_heat<-aver_dt_df_long %>% ggplot(aes(Bcell_group,gene)) +
  geom_tile(data=aver_dt_df_long,aes(fill=scale(value)), colour = "white", linewidth = 0.5)+
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
        axis.text.x = element_text(size = 17,color = 'black',angle =40,hjust =1,vjust=1),
        legend.title = element_text(size=15, color = "black"),
        legend.text = element_text(size=15,color = "black",angle =0),
        legend.position = "right") +
  scale_y_discrete(position = "left")
p1_heat


ggsave("Bcell_细胞亚型分类热图合并.jpg",p1_heat,height = 11,width = 6.5)
ggsave("Bcell_细胞亚型分类热图合并.pdf",p1_heat,height = 11,width = 6.5)


####有无微生物细胞差异分析----
MPRNMPR_object_miMPRobe_remove<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/MPRNMPR_object_miMPRobe_remove.rds")
MPRNMPR_object_miMPRobe_remove = UpdateSeuratObject(MPRNMPR_object_miMPRobe_remove)
MPRNMPR_object_miMPRobe_remove

Idents(MPRNMPR_object_miMPRobe_remove) <- "cell_type_new"
MPRNMPR_object_miMPRobe_remove@meta.data$MPR_Response2
B_cells_object=subset(MPRNMPR_object_miMPRobe_remove,idents =  c("B cells",
                                                                       "Plasma cells"))

B_cells_object@meta.data$group_microbe2
Idents(B_cells_object) <- "group_microbe2"
B_cells_object_deg_all=FindMarkers(B_cells_object, ident.1 = c("Bacteria+"), ident.2 = "Bacteria-")
dim(B_cells_object_deg_all)
names(B_cells_object_deg_all)
write.csv(B_cells_object_deg_all,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/免疫细胞/B_cells_object_deg_all.csv")
B_cells_object_deg_all<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/免疫细胞/B_cells_object_deg_all.csv",header = T)

####根据过滤掉pt 小于0.1 的行

names(B_cells_object_deg_all)

B_cells_object_deg_all_filter<- B_cells_object_deg_all %>% filter(pct.1>0.1 &pct.2>0.1)

####读取免疫细胞相关基因，从差异基因集中筛选

gene_name<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/免疫细胞/免疫细胞gene_name.csv")

intersection <- intersect(gene_name$Gene.symbol, B_cells_object_deg_all_filter$X)

B_cells_object_deg_all_filter_inter<- B_cells_object_deg_all_filter %>% filter(B_cells_object_deg_all_filter$X %in%intersection )

df <- B_cells_object_deg_all_filter_inter
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
  ggtitle("B cells DEG (Bacteria+ vs Bacteria-)") + #标题
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
ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/免疫细胞/B细胞火山图.jpg",degp,width = 10,height=7)
ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/免疫细胞/B细胞火山图.pdf",degp,width = 10,height=7)

####OR指数----

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
out.prefix <- "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/B_cell注释/"
B_cells_object<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/B_cells_object_combined使用.rds")

library(dplyr)

# 转换元数据为数据框
metadata <- as.data.frame(B_cells_object@meta.data)

# 替换 B_cells_type 列中的值
metadata <- metadata %>%
  mutate(B_cells_type = recode(B_cells_type, 
                               `Naive plasma cell` = "Naive plasma cells"))

# 将修改后的元数据更新回 Seurat 对象
B_cells_object@meta.data <- metadata
B_cells_object$B_cells_type

saveRDS(B_cells_object,"D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/B_cells_object_combined使用.rds")

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

meta.tb <- B_cells_object@meta.data

names(meta.tb)
View(meta.tb)

A <- do.tissueDist(cellInfo.tb = meta.tb,
                   meta.cluster = meta.tb$B_cells_type,
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
data$rid <- factor(data$rid, levels=c("Pre-B cells", "Primarily mature naive B cells", 
                                      "Naive B cells", "Active B cells", "Naive plasma cell", 
                                      "Plasma cells"))
head(data)
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


ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/B_cell注释/OR指数/OR_plot_MPR_treatments.pdf",OR_plot,width=9,height=5)


####细胞比例计算----

Cellratio <- prop.table(table(B_cells_object@meta.data$B_cells_type,B_cells_object@meta.data$orig.ident), margin = 2)

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


cell_type_groups = c("Pre-B cells", "Primarily mature naive B cells", 
                     "Naive B cells", "Active B cells", "Naive plasma cells", 
                     "Plasma cells")
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
cells_<-plot_grid(pplist[['Pre-B cells']],
                  pplist[['Primarily mature naive B cells']],
                  pplist[['Naive B cells']],
                  pplist[['Active B cells']],
                  pplist[['Naive plasma cells']],
                  align = "h",  #axis = 'l',
                  nrow =2)

ggsave('D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/B_cell注释/细胞百分比/B_cells_MPR_percentage.pdf',cells_,width =10,height =8,limitsize = FALSE)

####B 细胞分组火山图----
B_cells_object$UMIs_counts

library(Seurat)
library(ggplot2)
library(cowplot)
library(dplyr)
library(tidydr)
library(stringr)
library(viridis)
library(CellChat)

pal <- viridis(n = 15, option = "D", direction = -1)
FeaturePlot(object =B_cells_object,features = "UMIs_counts",reduction = "umap.harmony",raster=FALSE, cols = pal, order = T)


####T 细胞B细胞细胞通讯分析----

T_cells_object<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/T_cells_object.rds")
T_cells_object = UpdateSeuratObject(T_cells_object)
T_cells_object


B_cells_object<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/B_cells_object_combined使用.rds")
B_cells_object = UpdateSeuratObject(B_cells_object)
B_cells_object

library(CellChat)
library(tidyverse)
library(Seurat)
options(stringsAsFactors = FALSE)
####合并T细胞和B细胞的seurat
merged_obj <- merge(T_cells_object, B_cells_object, 
                    add.cell.ids = c("T_cells_object", "B_cells_object"), 
                    project = "Combined_Samples")
head(merged_obj@meta.data)

merged_obj@meta.data$merged_group <- paste(merged_obj@meta.data$cd4_cd8_group, 
                                            merged_obj@meta.data$B_cells_type, sep = "_")

# 查看更新后的 meta.data
head(merged_obj@meta.data)

merged_obj@meta.data$merged_group <- gsub("^NA_", "", merged_obj@meta.data$merged_group)
merged_obj@meta.data$merged_group <- gsub("^NA_", "", merged_obj@meta.data$merged_group)
merged_obj@meta.data$merged_group <- gsub("_NA$", "", merged_obj@meta.data$merged_group)

unique(merged_obj@meta.data$merged_group)

####创建CellChat对象
####可以使用矩阵数据、Seurat或SingleCellExperiment 对象中创建CellChat对象。如果是Seurat或SingleCellExperiment 对象，meta信息会默认使用meta.data（当然也可以指定），通过 group.by 定义分组。

cellchat <- createCellChat(object = merged_obj,
                           meta = merged_obj@meta.data,
                           group.by = "merged_group")
cellchat

####设置参考数据库
####根据分析数据的物种，可选CellChatDB.human, 或者 CellChatDB.mouse 。通过showDatabaseCategory函数可以查看该数据库的情况
CellChatDB <- CellChatDB.human  
showDatabaseCategory(CellChatDB)
# use a subset of CellChatDB for cell-cell communication analysis
#CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") 
# set the used database in the object
#cellchat@DB <- CellChatDB.use

cellchat@DB <- CellChatDB

#人的数据包括61.8%的旁分泌/自分泌信号互作、21.7%的细胞外基质(ECM)受体互作和16.5%的细胞-细胞通讯互作。

#注1：如果想用全部的用于cellchat分析，不进行subsetDB，直接指定cellchat@DB <- CellChatDB 即可。

#注2：如果你有关心的配受体对 且 不在该数据库中，也可以自行添加上。大概步骤就是下载对应的csv（数据库），在对应的列上添加上你的配受体对信息，保存后重新读取新的csv即可，详细见https://htmlpreview.github.io/?https://github.com/sqjin/CellChat/blob/master/tutorial/Update-CellChatDB.html。

####表达数据的预处理

# This step is necessary even if using the whole database
cellchat <- subsetData(cellchat) 
# do parallel ，根据配置设置
#future::plan("multiprocess", workers = 1) 

#识别过表达基因
cellchat <- identifyOverExpressedGenes(cellchat)
#识别过表达配体受体对
cellchat <- identifyOverExpressedInteractions(cellchat)

#project gene expression data onto PPI (Optional: when running it, USER should set `raw.use = FALSE` in the function `computeCommunProb()` in order to use the projected data)
cellchat <- projectData(cellchat, PPI.human)
cellchat@data.project[1:4,1:4]

####耗时
####推断 cell-cell communication network
cellchat <- computeCommunProb(cellchat, raw.use = TRUE, population.size = TRUE) 
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)

####注1：raw.use = TRUE 表示使用raw数据，而不使用上一步projectData后的结果。

####注2：在假设细胞数较多的群 往往比 细胞数较少的群发送更强的信号的前提下，当population.size = TRUE时候，CellChat可以在概率计算中考虑每个细胞群中细胞比例的影响。
####保存结果
####根据需求保存结果，默认保存全部结果（推荐），设置slot.name = "netP" 保存显著的结果，指定sources.use和targets.use则能获取指定配体受体方向的结果，signaling获取指定signaling的结果。


#all the inferred cell-cell communications at the level of ligands/receptors
df.net <- subsetCommunication(cellchat)
write.csv(df.net, "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/B_cell注释/细胞通讯BT/cell-cell_communications_ligands_receptors.csv")

#access the the inferred communications at the level of signaling pathways
df.net1 <- subsetCommunication(cellchat,slot.name = "netP")
write.csv(df.net1, "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/B_cell注释/细胞通讯BT/cell-cell_communications_signaling_pathways.csv")

#gives the inferred cell-cell communications sending from cell groups 1 and 2 to cell groups 4 and 5.
levels(cellchat@idents)

df.net2 <- subsetCommunication(cellchat, sources.use = c("CD8_C1-GNLY"), targets.use = c("Naive B cells" ,             
                                                                                           "Naive plasma cells","Active B cells" ,"Plasma cells",                  
                                                                                            "Pre-B cells","Primarily mature naive B cells")) 

write.csv(df.net2, "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/B_cell注释/细胞通讯BT/cell-cell_communications_CD8_C1-GNLY.csv")

df.net_ctsw <- subsetCommunication(cellchat, sources.use = c("CD8_C2-CTSW"), targets.use = c("Naive B cells" ,             
                                                                                         "Naive plasma cells","Active B cells" ,"Plasma cells",                  
                                                                                         "Pre-B cells","Primarily mature naive B cells")) 

write.csv(df.net_ctsw, "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/B_cell注释/细胞通讯BT/cell-cell_communications_CD8_ctsw.csv")



#gives the inferred cell-cell communications mediated by signaling WNT and TGFb.
df.net3 <- subsetCommunication(cellchat, signaling = c("CCL", "TGFb"))

write.csv(df.net3, "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/B_cell注释/细胞通讯BT/cell-cell_communications_WNT_TGFb.csv")

#3，计算cell-cell communication
#使用computeCommunProbPathway计算每个信号通路的所有配体-受体相互作用的通信结果，结存存放在net 和 netP中 。

#然后使用aggregateNet计算细胞类型间整合的细胞通讯结果

#计算每个信号通路相关的所有配体-受体相互作用的通信结果
cellchat <- computeCommunProbPathway(cellchat)
#计算整合的细胞类型之间通信结果
cellchat <- aggregateNet(cellchat)

####CellChat 可视化
####celltype之间通讯结果
####1）根据使用netVisual_circle显示任意两个celltype之间的通讯次数（左）或总通讯强度(右)

groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, 
                 weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, 
                 weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

####
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, 
                 weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, 
                 weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

netVisual_circle(df.net2@net$weight, vertex.weight = groupSize, 
                 weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
df.net2

mat_count <- cellchat@net$count
write.csv(mat_count,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/B_cell注释/细胞通讯BT/cell-cell_communications_mat_count.csv")
mat_weight<- cellchat@net$weight
write.csv(mat_weight,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/B_cell注释/细胞通讯BT/cell-cell_communications_mat_weight.csv")
mat_prob<-cellchat@netP$prob
write.csv(mat_prob,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/B_cell注释/细胞通讯BT/cell-cell_communications_mat_prob.csv")
print(mat)
par(mfrow = c(2,3), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
colnames(mat2)
mat2[8,c(1,12:17)] <- mat[8,c(1,12:17)]

rownames(mat2)

netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat2), title.name = rownames(mat2)[2])

saveRDS(cellchat,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/B_cell注释/细胞通讯BT/cellchat.rds")


if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("CCPlotR")

## or for development version:

devtools::install_github("Sarah145/CCPlotR")
library(CCPlotR)
data(toy_data, toy_exp, package = 'CCPlotR')
head(toy_data)
head(toy_exp)

#'A', 'B', 'CellPhoneDB' or 'Liana'
cc_heatmap(toy_data,option = "A", n_top_ints = 30)
?cc_heatmap
head(df.net2)

######绘制环图----

df.GNLY<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/B_cell注释/细胞通讯BT/cell-cell_communications_CD8_C1-GNLY.csv")
head(df.GNLY)
df<- df.GNLY[,-1]
head(df)
colnames(df)[5]<- "score"
pdf(file = "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/B_cell注释/细胞通讯BT/cc_heatmap_GNLY.pdf", width = 10, height = 5)
cc_heatmap(df)
dev.off()

pdf(file = "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/B_cell注释/细胞通讯BT/cc_heatmap_GNLY受配体.pdf", width = 10, height = 5)
cc_heatmap(df, option = "B", n_top_ints = 10)
dev.off()

df.ctsw<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/B_cell注释/细胞通讯BT/cell-cell_communications_CD8_ctsw.csv")
head(df.ctsw)
df.net_ctsw<- df.ctsw[,-1]
head(df.net_ctsw)
colnames(df.net_ctsw)[5]<- "score"
cc_heatmap(df.net_ctsw)
cc_heatmap(df.net_ctsw, option = "B", n_top_ints = 10)
cc_heatmap(df.net_ctsw, option = "CellPhoneDB")
cc_dotplot(toy_data)
cc_dotplot(toy_data, option = "B", n_top_ints = 10)
cc_dotplot(toy_data, option = "Liana", n_top_ints = 15)
cc_network(df.net_ctsw)
cc_network(df.net_ctsw, colours = c("orange", "cornflowerblue", "hotpink"), option = "B")
cc_circos(df.net_ctsw)
cc_circos(toy_data, option = "B", n_top_ints = 10, cex = 0.75)

merged_obj$merged_group
avg_expression<-AverageExpression(merged_obj, group.by = "merged_group")

library(tidyr)

# 将矩阵转换为数据框
avg_expression_df <- as.data.frame(avg_expression$RNA)

# 将行名（基因名）转换为一列
avg_expression_df$gene <- rownames(avg_expression$RNA)

# 使用pivot_longer将列转换为长格式
long_data <- pivot_longer(avg_expression_df, 
                          cols = -gene,  # 除了Gene列以外的所有列
                          names_to = "cell_type", 
                          values_to = "mean_exp")

head(long_data)
long_data <- long_data[, c(2, 1, 3:ncol(long_data))]
head(long_data)
library(randomcoloR)
colors<- distinctColorPalette(6)
# "#D18CB5" "#8BE07A" "#AE66DD" "#D5D4BB" "#82CDD4" "#E2C55B"
df$target
unique(df$target)
scales::show_col(paletteMartin())
dev.new()
pdf(file = "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/B_cell注释/细胞通讯BT/cc_circos_p_GNLY.pdf", width = 8, height = 8)
 cc_circos(df, option = "C", n_top_ints = 15, exp_df = long_data,
              cell_cols =c("CD8_C1-GNLY"= "hotpink", "Naive B cells" = "orange", "Naive plasma cells" = "cornflowerblue",
                           "Active B cells"="#8BE07A", "Plasma cells"="#82CDD4","Pre-B cells"="#D5D4BB","Primarily mature naive B cells" ="#AE66DD"), palette = "PuRd", cex = 0.75)
dev.off()
graphics.off()

p_ctsw<- cc_circos(df.net_ctsw, option = "C", n_top_ints = 15, exp_df = long_data,
          cell_cols =c("CD8_C2-CTSW"= "hotpink", `Naive B cells` = "orange", `Naive plasma cells` = "cornflowerblue",
                       "Active B cells"="#8BE07A", "Plasma cells"="#82CDD4","Pre-B cells"="#D5D4BB","Primarily mature naive B cells" ="#AE66DD"), palette = "PuRd", cex = 0.75)

?cc_circos
dev.off()  # 关闭当前的绘图设备
plot.new() # 

####分组细胞通讯分析----

######pre-nmpr----
library(CCPlotR)
library(CellChat)
library(tidyverse)
library(Seurat)

T_cells_object<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/T_cells_object.rds")
T_cells_object = UpdateSeuratObject(T_cells_object)
T_cells_object

B_cells_object<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/B_cells_object_combined使用.rds")
B_cells_object = UpdateSeuratObject(B_cells_object)
B_cells_object

merged_obj <- merge(T_cells_object, B_cells_object, 
                    add.cell.ids = c("T_cells_object", "B_cells_object"), 
                    project = "Combined_Samples")

head(merged_obj@meta.data)

merged_obj@meta.data$merged_group <- paste(merged_obj@meta.data$cd4_cd8_group, 
                                           merged_obj@meta.data$B_cells_type, sep = "_")

# 查看更新后的 meta.data
head(merged_obj@meta.data)
merged_obj@meta.data$merged_group <- gsub("^NA_", "", merged_obj@meta.data$merged_group)
merged_obj@meta.data$merged_group <- gsub("^NA_", "", merged_obj@meta.data$merged_group)
merged_obj@meta.data$merged_group <- gsub("_NA$", "", merged_obj@meta.data$merged_group)

unique(merged_obj@meta.data$merged_group)

unique(merged_obj$MPR_Response2)


####相关性热图----
merged_obj_filter<- subset(merged_obj,subset=merged_group %in% c("CD8_C2-CTSW","CD8_C1-GNLY","Active B cells",                
                                                                 "Naive plasma cells","Plasma cells","Primarily mature naive B cells",
                                                                 "Naive B cells","Pre-B cells"))

merged_obj_filter<- subset(merged_obj,subset=merged_group %in% c("CD8_C3-IL7R","CD8_C4-GZMK","CD8_C1-GNLY","Active B cells",                
                                                                 "Naive plasma cells","Plasma cells","Primarily mature naive B cells",
                                                                 "Naive B cells","Pre-B cells"))

merged_obj_filter<- subset(merged_obj,subset=merged_group %in% c("CD8_C3-IL7R","CD8_C4-GZMK","Active B cells",                
                                                                 "Naive plasma cells","Plasma cells","Primarily mature naive B cells",
                                                                 "Naive B cells","Pre-B cells"))

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

# 数据清洗
varechem <- read_tsv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/B_cell注释/细胞通讯BT/相关性热图/varechem.xls") %>% 
  column_to_rownames(var="id")
table1 <- varechem %>% rownames_to_column(var="id") %>% 
  pivot_longer(-id) %>% 
  pivot_wider(values_from = value) %>% 
  column_to_rownames(var="id")
varespec <- read_tsv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/B_cell注释/细胞通讯BT/相关性热图/varespec.xls") %>% column_to_rownames(var="id") 

mantel <- mantel_test(varespec, table1,
                      spec_select = list(Spec01 = 1:7,
                                         Spec02 = 8:18,
                                         Spec03 = 19:20)) %>% 
  mutate(rd = cut(r, breaks = c(-Inf, 0.2, 0.4, Inf),
                  labels = c("< 0.2", "0.2 - 0.4", ">= 0.4")),
         pd = cut(p, breaks = c(-Inf, 0.01, 0.05, Inf),
                  labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))

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
  pivot_longer(cols = c(`CD8-C1-GNLY`, `CD8-C2-CTSW`), 
               names_to = "Group", 
               values_to = "Value") %>% as.data.frame()


long_df <- cor_df %>%
  pivot_longer(cols = c(`CD8-C1-GNLY`,"CD8-C3-IL7R","CD8-C4-GZMK"), 
               names_to = "Group", 
               values_to = "Value") %>% as.data.frame()

long_df <- cor_df %>%
  pivot_longer(cols = c(`CD8-C1-GNLY`, `CD8-C2-CTSW`,"CD8-C3-IL7R","CD8-C4-GZMK"), 
               names_to = "Group", 
               values_to = "Value") %>% as.data.frame()


long_df <- cor_df %>%
  pivot_longer(cols = c("CD8-C3-IL7R","CD8-C4-GZMK"), 
               names_to = "Group", 
               values_to = "Value") %>% as.data.frame()

long_df <- long_df[, c("Group", "Other_Group", "Value")]


names(long_df)<-c("spec","env","r")


p_df<- as.data.frame(res1$p)
p_df$Other_Group <- rownames(p_df)


p_long_df <- p_df %>%
  pivot_longer(cols = c(`CD8-C1-GNLY`, `CD8-C2-CTSW`), 
               names_to = "Group", 
               values_to = "Value") %>% as.data.frame()


p_long_df <- p_df %>%
  pivot_longer(cols = c(`CD8-C1-GNLY`, "CD8-C3-IL7R","CD8-C4-GZMK"), 
               names_to = "Group", 
               values_to = "Value") %>% as.data.frame()

p_long_df <- p_df %>%
  pivot_longer(cols = c("CD8-C3-IL7R","CD8-C4-GZMK"), 
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

# 使用 Pearson 方法计算表格 table1 的相关系数矩阵
cor_matrix_filtered[lower.tri(cor_matrix_filtered1)] <- NA
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
ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/B_cell注释/细胞通讯BT/相关性热图/BT相关性热图_GNLY_IL7R_GZMK.pdf",bt_plot,width=6,height =6)

####细胞通讯----

pre_nmpr <- subset(merged_obj, subset =MPR_Response2 == 'Pre_NMPR')

pre_nmpr

cellchat <- createCellChat(object = pre_nmpr,
                           meta = pre_nmpr@meta.data,
                           group.by = "merged_group")
cellchat
####设置参考数据库
####根据分析数据的物种，可选CellChatDB.human, 或者 CellChatDB.mouse 。通过showDatabaseCategory函数可以查看该数据库的情况
CellChatDB <- CellChatDB.human  
showDatabaseCategory(CellChatDB)
cellchat@DB <- CellChatDB
cellchat <- subsetData(cellchat) 
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
cellchat@data.project[1:4,1:4]
cellchat <- computeCommunProb(cellchat, raw.use = TRUE, population.size = TRUE) 
cellchat <- filterCommunication(cellchat, min.cells = 10)
df.net <- subsetCommunication(cellchat)
write.csv(df.net, "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/B_cell注释/细胞通讯BT/pre_nmpr_cell-cell_communications_ligands_receptors.csv")
df.net1 <- subsetCommunication(cellchat,slot.name = "netP")
write.csv(df.net1, "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/B_cell注释/细胞通讯BT/pre_nmpr_cell-cell_communications_signaling_pathways.csv")
levels(cellchat@idents)
df.net2 <- subsetCommunication(cellchat, sources.use = c("CD8_C1-GNLY"), targets.use = c("Naive B cells" ,             
                                                                                         "Naive plasma cells","Active B cells" ,"Plasma cells",                  
                                                                                         "Pre-B cells","Primarily mature naive B cells")) 

write.csv(df.net2, "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/B_cell注释/细胞通讯BT/pre_nmpr_cell-cell_communications_CD8_C1-GNLY.csv")

df.net_ctsw <- subsetCommunication(cellchat, sources.use = c("CD8_C2-CTSW"), targets.use = c("Naive B cells" ,             
                                                                                             "Naive plasma cells","Active B cells" ,"Plasma cells",                  
                                                                                             "Pre-B cells","Primarily mature naive B cells")) 

write.csv(df.net_ctsw, "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/B_cell注释/细胞通讯BT/pre_nmpr_cell-cell_communications_CD8_ctsw.csv")

df.net3 <- subsetCommunication(cellchat, signaling = c("CCL", "TGFb"))
write.csv(df.net3, "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/B_cell注释/细胞通讯BT/pre_nmpr_cell-cell_communications_WNT_TGFb.csv")

cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
saveRDS(cellchat,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/B_cell注释/细胞通讯BT/pre_nmpr_cellchat.rds")

df_pre_nmpr_gnly<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/B_cell注释/细胞通讯BT/pre_nmpr_cell-cell_communications_CD8_C1-GNLY.csv")
df.net_gnly<- df_pre_nmpr_gnly[,-1]
colnames(df.net_gnly)[5]<- "score"
head(df.net_gnly)

df_pre_nmpr_ctsw<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/B_cell注释/细胞通讯BT/pre_nmpr_cell-cell_communications_CD8_ctsw.csv")
df.net_ctsw<- df_pre_nmpr_ctsw[,-1]
colnames(df.net_ctsw)[5]<- "score"
head(df.net_ctsw)

avg_expression<-AverageExpression(merged_obj, group.by = "merged_group")

library(tidyr)
# 将矩阵转换为数据框
avg_expression_df <- as.data.frame(avg_expression$RNA)

# 将行名（基因名）转换为一列
avg_expression_df$gene <- rownames(avg_expression$RNA)

# 使用pivot_longer将列转换为长格式
long_data <- pivot_longer(avg_expression_df, 
                          cols = -gene,  # 除了Gene列以外的所有列
                          names_to = "cell_type", 
                          values_to = "mean_exp")

head(long_data)
long_data <- long_data[, c(2, 1, 3:ncol(long_data))]
head(long_data)
cc_heatmap(df.net_gnly)
cc_heatmap(df.net_gnly, option = "B", n_top_ints = 10)

cc_heatmap(df.net_ctsw)
cc_heatmap(df.net_ctsw, option = "B", n_top_ints = 10)

dev.new()
cc_circos(df.net_gnly, option = "C", n_top_ints = 15, exp_df = long_data,
          cell_cols =c("CD8_C1-GNLY"= "hotpink", "Naive B cells" = "orange", "Naive plasma cells" = "cornflowerblue",
                       "Active B cells"="#8BE07A", "Plasma cells"="#82CDD4","Pre-B cells"="#D5D4BB","Primarily mature naive B cells" ="#AE66DD"), palette = "PuRd", cex = 0.75)
dev.off()


cc_circos(df.net_ctsw, option = "C", n_top_ints = 15, exp_df = long_data,
                   cell_cols =c("CD8_C2-CTSW"= "hotpink", `Naive B cells` = "orange", `Naive plasma cells` = "cornflowerblue",
                                "Active B cells"="#8BE07A", "Plasma cells"="#82CDD4","Pre-B cells"="#D5D4BB","Primarily mature naive B cells" ="#AE66DD"), palette = "PuRd", cex = 0.75)


######pre-mpr----
library(CCPlotR)
library(CellChat)
library(tidyverse)
library(Seurat)

T_cells_object<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/T_cells_object.rds")
T_cells_object = UpdateSeuratObject(T_cells_object)
T_cells_object

B_cells_object<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/B_cells_object_combined使用.rds")
B_cells_object = UpdateSeuratObject(B_cells_object)
B_cells_object

merged_obj <- merge(T_cells_object, B_cells_object, 
                    add.cell.ids = c("T_cells_object", "B_cells_object"), 
                    project = "Combined_Samples")

head(merged_obj@meta.data)

merged_obj@meta.data$merged_group <- paste(merged_obj@meta.data$cd4_cd8_group, 
                                           merged_obj@meta.data$B_cells_type, sep = "_")

# 查看更新后的 meta.data
head(merged_obj@meta.data)
merged_obj@meta.data$merged_group <- gsub("^NA_", "", merged_obj@meta.data$merged_group)
merged_obj@meta.data$merged_group <- gsub("^NA_", "", merged_obj@meta.data$merged_group)
merged_obj@meta.data$merged_group <- gsub("_NA$", "", merged_obj@meta.data$merged_group)

unique(merged_obj@meta.data$merged_group)

unique(merged_obj$MPR_Response2)

pre_mpr <- subset(merged_obj, subset =MPR_Response2 == 'Pre_MPR')

pre_mpr

cellchat <- createCellChat(object = pre_mpr,
                           meta = pre_mpr@meta.data,
                           group.by = "merged_group")
cellchat
####设置参考数据库
####根据分析数据的物种，可选CellChatDB.human, 或者 CellChatDB.mouse 。通过showDatabaseCategory函数可以查看该数据库的情况
CellChatDB <- CellChatDB.human  
showDatabaseCategory(CellChatDB)
cellchat@DB <- CellChatDB
cellchat <- subsetData(cellchat) 
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
cellchat@data.project[1:4,1:4]
cellchat <- computeCommunProb(cellchat, raw.use = TRUE, population.size = TRUE) 
cellchat <- filterCommunication(cellchat, min.cells = 10)
df.net <- subsetCommunication(cellchat)
write.csv(df.net, "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/B_cell注释/细胞通讯BT/pre_mpr_cell-cell_communications_ligands_receptors.csv")
df.net1 <- subsetCommunication(cellchat,slot.name = "netP")
write.csv(df.net1, "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/B_cell注释/细胞通讯BT/pre_mpr_cell-cell_communications_signaling_pathways.csv")
levels(cellchat@idents)
df.net2 <- subsetCommunication(cellchat, sources.use = c("CD8_C1-GNLY"), targets.use = c("Naive B cells" ,             
                                                                                         "Naive plasma cells","Active B cells" ,"Plasma cells",                  
                                                                                         "Pre-B cells","Primarily mature naive B cells")) 

write.csv(df.net2, "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/B_cell注释/细胞通讯BT/pre_mpr_cell-cell_communications_CD8_C1-GNLY.csv")

df.net_ctsw <- subsetCommunication(cellchat, sources.use = c("CD8_C2-CTSW"), targets.use = c("Naive B cells" ,             
                                                                                             "Naive plasma cells","Active B cells" ,"Plasma cells",                  
                                                                                             "Pre-B cells","Primarily mature naive B cells")) 

write.csv(df.net_ctsw, "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/B_cell注释/细胞通讯BT/pre_mpr_cell-cell_communications_CD8_ctsw.csv")

df.net3 <- subsetCommunication(cellchat, signaling = c("CCL", "TGFb"))
write.csv(df.net3, "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/B_cell注释/细胞通讯BT/pre_mpr_cell-cell_communications_WNT_TGFb.csv")

cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
saveRDS(cellchat,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/B_cell注释/细胞通讯BT/pre_mpr_cellchat.rds")

df_pre_mpr_gnly<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/B_cell注释/细胞通讯BT/pre_mpr_cell-cell_communications_CD8_C1-GNLY.csv")
df.net_gnly<- df_pre_mpr_gnly[,-1]
colnames(df.net_gnly)[5]<- "score"
head(df.net_gnly)

df_pre_mpr_ctsw<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/B_cell注释/细胞通讯BT/pre_mpr_cell-cell_communications_CD8_ctsw.csv")
df.net_ctsw<- df_pre_mpr_ctsw[,-1]
colnames(df.net_ctsw)[5]<- "score"
head(df.net_ctsw)

avg_expression<-AverageExpression(merged_obj, group.by = "merged_group")

library(tidyr)
# 将矩阵转换为数据框
avg_expression_df <- as.data.frame(avg_expression$RNA)

# 将行名（基因名）转换为一列
avg_expression_df$gene <- rownames(avg_expression$RNA)

# 使用pivot_longer将列转换为长格式
long_data <- pivot_longer(avg_expression_df, 
                          cols = -gene,  # 除了Gene列以外的所有列
                          names_to = "cell_type", 
                          values_to = "mean_exp")

head(long_data)
long_data <- long_data[, c(2, 1, 3:ncol(long_data))]
head(long_data)
cc_heatmap(df.net_gnly)
cc_heatmap(df.net_gnly, option = "B", n_top_ints = 10)

cc_heatmap(df.net_ctsw)
cc_heatmap(df.net_ctsw, option = "B", n_top_ints = 10)

dev.new()
cc_circos(df.net_gnly, option = "C", n_top_ints = 15, exp_df = long_data,
          cell_cols =c("CD8_C1-GNLY"= "hotpink", "Naive B cells" = "orange", "Naive plasma cells" = "cornflowerblue",
                       "Active B cells"="#8BE07A", "Plasma cells"="#82CDD4","Pre-B cells"="#D5D4BB","Primarily mature naive B cells" ="#AE66DD"), palette = "PuRd", cex = 0.75)
dev.off()


cc_circos(df.net_ctsw, option = "C", n_top_ints = 15, exp_df = long_data,
                   cell_cols =c("CD8_C2-CTSW"= "hotpink", `Naive B cells` = "orange", `Naive plasma cells` = "cornflowerblue",
                                "Active B cells"="#8BE07A", "Plasma cells"="#82CDD4","Pre-B cells"="#D5D4BB","Primarily mature naive B cells" ="#AE66DD"), palette = "PuRd", cex = 0.75)



######post-mpr----
library(CCPlotR)
library(CellChat)
library(tidyverse)
library(Seurat)

T_cells_object<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/T_cells_object.rds")
T_cells_object = UpdateSeuratObject(T_cells_object)
T_cells_object

B_cells_object<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/B_cells_object_combined使用.rds")
B_cells_object = UpdateSeuratObject(B_cells_object)
B_cells_object

merged_obj <- merge(T_cells_object, B_cells_object, 
                    add.cell.ids = c("T_cells_object", "B_cells_object"), 
                    project = "Combined_Samples")

head(merged_obj@meta.data)

merged_obj@meta.data$merged_group <- paste(merged_obj@meta.data$cd4_cd8_group, 
                                           merged_obj@meta.data$B_cells_type, sep = "_")

# 查看更新后的 meta.data
head(merged_obj@meta.data)
merged_obj@meta.data$merged_group <- gsub("^NA_", "", merged_obj@meta.data$merged_group)
merged_obj@meta.data$merged_group <- gsub("^NA_", "", merged_obj@meta.data$merged_group)
merged_obj@meta.data$merged_group <- gsub("_NA$", "", merged_obj@meta.data$merged_group)

unique(merged_obj@meta.data$merged_group)

unique(merged_obj$MPR_Response2)

post_mpr <- subset(merged_obj, subset =MPR_Response2 == 'Post_MPR')

post_mpr

cellchat <- createCellChat(object = post_mpr,
                           meta = post_mpr@meta.data,
                           group.by = "merged_group")
cellchat
####设置参考数据库
####根据分析数据的物种，可选CellChatDB.human, 或者 CellChatDB.mouse 。通过showDatabaseCategory函数可以查看该数据库的情况
CellChatDB <- CellChatDB.human  
showDatabaseCategory(CellChatDB)
cellchat@DB <- CellChatDB
cellchat <- subsetData(cellchat) 
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
cellchat@data.project[1:4,1:4]
cellchat <- computeCommunProb(cellchat, raw.use = TRUE, population.size = TRUE) 
cellchat <- filterCommunication(cellchat, min.cells = 10)
df.net <- subsetCommunication(cellchat)
write.csv(df.net, "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/B_cell注释/细胞通讯BT/post_mpr_cell-cell_communications_ligands_receptors.csv")
df.net1 <- subsetCommunication(cellchat,slot.name = "netP")
write.csv(df.net1, "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/B_cell注释/细胞通讯BT/post_mpr_cell-cell_communications_signaling_pathways.csv")
levels(cellchat@idents)
df.net2 <- subsetCommunication(cellchat, sources.use = c("CD8_C1-GNLY"), targets.use = c("Naive B cells" ,             
                                                                                         "Naive plasma cells","Active B cells" ,"Plasma cells",                  
                                                                                         "Pre-B cells","Primarily mature naive B cells")) 

write.csv(df.net2, "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/B_cell注释/细胞通讯BT/post_mpr_cell-cell_communications_CD8_C1-GNLY.csv")

df.net_ctsw <- subsetCommunication(cellchat, sources.use = c("CD8_C2-CTSW"), targets.use = c("Naive B cells" ,             
                                                                                             "Naive plasma cells","Active B cells" ,"Plasma cells",                  
                                                                                             "Pre-B cells","Primarily mature naive B cells")) 

write.csv(df.net_ctsw, "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/B_cell注释/细胞通讯BT/post_mpr_cell-cell_communications_CD8_ctsw.csv")

df.net3 <- subsetCommunication(cellchat, signaling = c("CCL", "TGFb"))
write.csv(df.net3, "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/B_cell注释/细胞通讯BT/post_mpr_cell-cell_communications_WNT_TGFb.csv")

cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
saveRDS(cellchat,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/B_cell注释/细胞通讯BT/post_mpr_cellchat.rds")

df_post_mpr_gnly<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/B_cell注释/细胞通讯BT/post_mpr_cell-cell_communications_CD8_C1-GNLY.csv")
df.net_gnly<- df_post_mpr_gnly[,-1]
colnames(df.net_gnly)[5]<- "score"
head(df.net_gnly)

df_post_mpr_ctsw<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/B_cell注释/细胞通讯BT/post_mpr_cell-cell_communications_CD8_ctsw.csv")
df.net_ctsw<- df_post_mpr_ctsw[,-1]
colnames(df.net_ctsw)[5]<- "score"
head(df.net_ctsw)

avg_expression<-AverageExpression(merged_obj, group.by = "merged_group")

library(tidyr)
# 将矩阵转换为数据框
avg_expression_df <- as.data.frame(avg_expression$RNA)

# 将行名（基因名）转换为一列
avg_expression_df$gene <- rownames(avg_expression$RNA)

# 使用pivot_longer将列转换为长格式
long_data <- pivot_longer(avg_expression_df, 
                          cols = -gene,  # 除了Gene列以外的所有列
                          names_to = "cell_type", 
                          values_to = "mean_exp")

head(long_data)
long_data <- long_data[, c(2, 1, 3:ncol(long_data))]
head(long_data)
cc_heatmap(df.net_gnly)
cc_heatmap(df.net_gnly, option = "B", n_top_ints = 10)

cc_heatmap(df.net_ctsw)
cc_heatmap(df.net_ctsw, option = "B", n_top_ints = 10)

dev.new()
cc_circos(df.net_gnly, option = "C", n_top_ints = 15, exp_df = long_data,
          cell_cols =c("CD8_C1-GNLY"= "hotpink", "Naive B cells" = "orange", "Naive plasma cells" = "cornflowerblue",
                       "Active B cells"="#8BE07A", "Plasma cells"="#82CDD4","Pre-B cells"="#D5D4BB","Primarily mature naive B cells" ="#AE66DD"), palette = "PuRd", cex = 0.75)
dev.off()

cc_circos(df.net_ctsw, option = "C", n_top_ints = 15, exp_df = long_data,
          cell_cols =c("CD8_C2-CTSW"= "hotpink", `Naive B cells` = "orange", `Naive plasma cells` = "cornflowerblue",
                       "Active B cells"="#8BE07A", "Plasma cells"="#82CDD4","Pre-B cells"="#D5D4BB","Primarily mature naive B cells" ="#AE66DD"), palette = "PuRd", cex = 0.75)


######KS 细胞通讯方法----
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

library(CCPlotR)


T_cells_object<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/T_cells_object.rds")
T_cells_object = UpdateSeuratObject(T_cells_object)
T_cells_object

B_cells_object<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/B_cells_object_combined使用.rds")
B_cells_object = UpdateSeuratObject(B_cells_object)
B_cells_object

merged_obj <- merge(T_cells_object, B_cells_object, 
                    add.cell.ids = c("T_cells_object", "B_cells_object"), 
                    project = "Combined_Samples")

head(merged_obj@meta.data)

merged_obj@meta.data$merged_group <- paste(merged_obj@meta.data$cd4_cd8_group, 
                                           merged_obj@meta.data$B_cells_type, sep = "_")

# 查看更新后的 meta.data
head(merged_obj@meta.data)
merged_obj@meta.data$merged_group <- gsub("^NA_", "", merged_obj@meta.data$merged_group)
merged_obj@meta.data$merged_group <- gsub("^NA_", "", merged_obj@meta.data$merged_group)
merged_obj@meta.data$merged_group <- gsub("_NA$", "", merged_obj@meta.data$merged_group)

unique(merged_obj@meta.data$merged_group)

unique(merged_obj$MPR_Response2)

merged_obj_filter<- subset(merged_obj,subset=merged_group %in% c("Pre-B cells", "Primarily mature naive B cells", 
                                                                 "Naive B cells", "Active B cells", "Naive plasma cells", 
                                                                 "Plasma cells","CD8_C2-CTSW","CD8_C1-GNLY","CD8_C3-IL7R","CD8_C4-GZMK"))

pre_nmpr_bt <- subset(merged_obj_filter, subset =MPR_Response2 == 'Pre_NMPR')
unique(pre_nmpr_bt$MPR_Response2)
pre_mpr_bt <- subset(merged_obj_filter, subset =MPR_Response2 == 'Pre_MPR')
unique(pre_mpr_bt$MPR_Response2)
post_mpr_bt <- subset(merged_obj_filter, subset =MPR_Response2 == 'Post_MPR')
unique(post_mpr_bt$MPR_Response2)

pre_nmpr.bt.cellchat <- KS_cellchat(pre_nmpr_bt, assay = 'RNA', group.by = "merged_group",
                                 workers=2, species='human')

saveRDS(pre_nmpr.bt.cellchat, file = "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/B_cell注释/细胞通讯BT/KS的代码重新分析画图/pre_nmpr.bt.cellchat_all.rds")



pre_mpr.bt.cellchat <- KS_cellchat(pre_mpr_bt, assay = 'RNA', group.by = "merged_group",
                                    workers=2, species='human')

saveRDS(pre_mpr.bt.cellchat, file = "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/B_cell注释/细胞通讯BT/KS的代码重新分析画图/pre_mpr.bt.cellchat_all.rds")


post_mpr.bt.cellchat <- KS_cellchat(post_mpr_bt, assay = 'RNA', group.by = "merged_group",
                                   workers=2, species='human')

saveRDS(post_mpr.bt.cellchat, file = "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/B_cell注释/细胞通讯BT/KS的代码重新分析画图/post_mpr.bt.cellchat_all.rds")

object.list <- list(Pre_NMPR= pre_nmpr.bt.cellchat, Pre_MPR= pre_mpr.bt.cellchat,Post_MPR=post_mpr.bt.cellchat)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))


gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2,3))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2,3), measure = "weight")
p1<-gg1 + gg2

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/B_cell注释/细胞通讯BT/KS的代码重新分析画图/分组比较柱状图_记忆和耗竭.pdf",p1,width=6,height=4)



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
ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/B_cell注释/细胞通讯BT/KS的代码重新分析画图/combined_plots记忆和耗竭.pdf", p_combined, width = 10, height = 4) 

treatment_color <- c("#4974a4","#4dae47","#f29600")


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
pathway.union

pbubble<-netVisual_bubble(cellchat, sources.use = "CD8_C1-GNLY", targets.use =c("Pre-B cells", "Primarily mature naive B cells", 
                                                                                "Naive B cells", "Active B cells", "Naive plasma cells", 
                                                                                "Plasma cells"),comparison = c(1,2,3), angle.x = 45,font.size = 12,color.text =c("#4974a4","#4dae47","#f29600"))
?netVisual_bubble

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/B_cell注释/细胞通讯BT/KS的代码重新分析画图/GNLY_受配体_combined_pbubble.pdf", pbubble, width = 9, height = 9) 




merged_obj_filter.bt.cellchat <- KS_cellchat(merged_obj_filter, assay = 'RNA', group.by = "merged_group",
                                    workers=2, species='human')

saveRDS(merged_obj_filter.bt.cellchat, file = "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/B_cell注释/细胞通讯BT/KS的代码重新分析画图/merged_obj_filter.bt.cellchat记忆和耗竭.rds")

df.net <- subsetCommunication(merged_obj_filter.bt.cellchat)
write.csv(df.net, "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/B_cell注释/细胞通讯BT/KS的代码重新分析画图/cell-cell_communications_ligands_receptors记忆和耗竭.csv")
df.net1 <- subsetCommunication(merged_obj_filter.bt.cellchat,slot.name = "netP")
write.csv(df.net1, "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/B_cell注释/细胞通讯BT/KS的代码重新分析画图/cell-cell_communications_signaling_pathways记忆和耗竭.csv")
levels(merged_obj_filter.bt.cellchat@idents)
df.net2 <- subsetCommunication(merged_obj_filter.bt.cellchat, sources.use = c("CD8_C1-GNLY"), targets.use = c("Naive B cells" ,             
                                                                                         "Naive plasma cells","Active B cells" ,"Plasma cells",                  
                                                                                         "Pre-B cells","Primarily mature naive B cells")) 
write.csv(df.net2, "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/B_cell注释/细胞通讯BT/KS的代码重新分析画图/cell-cell_communications_CD8_C1-GNLY记忆和耗竭.csv")

df.net_ctsw <- subsetCommunication(merged_obj_filter.bt.cellchat, sources.use = c("CD8_C2-CTSW"), targets.use = c("Naive B cells" ,             
                                                                                             "Naive plasma cells","Active B cells" ,"Plasma cells",                  
                                                                                             "Pre-B cells","Primarily mature naive B cells")) 
write.csv(df.net_ctsw, "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/B_cell注释/细胞通讯BT/KS的代码重新分析画图/cell-cell_communications_CD8_ctsw记忆和耗竭.csv")

merged_obj_filter.bt.cellchat<- computeCommunProbPathway(merged_obj_filter.bt.cellchat)
merged_obj_filter.bt.cellchat <- aggregateNet(merged_obj_filter.bt.cellchat)

df_pre_nmpr_gnly<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/B_cell注释/细胞通讯BT/KS的代码重新分析画图/cell-cell_communications_CD8_C1-GNLY记忆和耗竭.csv")
df.net_gnly<- df_pre_nmpr_gnly[,-1]
colnames(df.net_gnly)[5]<- "score"
head(df.net_gnly)

df_pre_nmpr_ctsw<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/B_cell注释/细胞通讯BT/KS的代码重新分析画图/cell-cell_communications_CD8_ctsw记忆和耗竭.csv")
df.net_ctsw<- df_pre_nmpr_ctsw[,-1]
colnames(df.net_ctsw)[5]<- "score"
head(df.net_ctsw)

avg_expression<-AverageExpression(merged_obj_filter, group.by = "merged_group")

library(tidyr)
avg_expression_df <- as.data.frame(avg_expression$RNA)

avg_expression_df$gene <- rownames(avg_expression$RNA)

long_data <- pivot_longer(avg_expression_df, 
                          cols = -gene,  # 除了Gene列以外的所有列
                          names_to = "cell_type", 
                          values_to = "mean_exp")

head(long_data)
long_data <- long_data[, c(2, 1, 3:ncol(long_data))]
head(long_data)
cc_heatmap(df.net_gnly)
cc_heatmap(df.net_gnly, option = "B", n_top_ints = 10)

cc_heatmap(df.net_ctsw)
cc_heatmap(df.net_ctsw, option = "B", n_top_ints = 10)

dev.new()
cc_circos(df.net_gnly, option = "C", n_top_ints = 15, exp_df = long_data,
          cell_cols =c("CD8_C1-GNLY"= "hotpink", "Naive B cells" = "orange", "Naive plasma cells" = "cornflowerblue",
                       "Active B cells"="#8BE07A", "Plasma cells"="#82CDD4","Pre-B cells"="#D5D4BB","Primarily mature naive B cells" ="#AE66DD"), palette = "PuRd", cex = 0.75)
dev.off()


cc_circos(df.net_ctsw, option = "C", n_top_ints = 15, exp_df = long_data,
          cell_cols =c("CD8_C2-CTSW"= "hotpink", `Naive B cells` = "orange", `Naive plasma cells` = "cornflowerblue",
                       "Active B cells"="#8BE07A", "Plasma cells"="#82CDD4","Pre-B cells"="#D5D4BB","Primarily mature naive B cells" ="#AE66DD"), palette = "PuRd", cex = 0.75)


##################################
df.net <- subsetCommunication(pre_mpr.bt.cellchat)
write.csv(df.net, "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/B_cell注释/细胞通讯BT/KS的代码重新分析画图/pre_mpr.cell-cell_communications_ligands_receptors记忆和耗竭.csv")
df.net1 <- subsetCommunication(pre_mpr.bt.cellchat,slot.name = "netP")
write.csv(df.net1, "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/B_cell注释/细胞通讯BT/KS的代码重新分析画图/pre_mpr.cell-cell_communications_signaling_pathways记忆和耗竭.csv")
levels(pre_mpr.bt.cellchat)
df.net2 <- subsetCommunication(pre_mpr.bt.cellchat, sources.use = c("CD8_C1-GNLY"), targets.use = c("Naive B cells" ,             
                                                                                                              "Naive plasma cells","Active B cells" ,"Plasma cells",                  
                                                                                                              "Pre-B cells","Primarily mature naive B cells")) 
write.csv(df.net2, "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/B_cell注释/细胞通讯BT/KS的代码重新分析画图/pre_mpr.cell-cell_communications_CD8_C1-GNLY.csv")

df.net_ctsw <- subsetCommunication(pre_mpr.bt.cellchat, sources.use = c("CD8_C2-CTSW"), targets.use = c("Naive B cells" ,             
                                                                                                                  "Naive plasma cells","Active B cells" ,"Plasma cells",                  
                                                                                                                  "Pre-B cells","Primarily mature naive B cells")) 
write.csv(df.net_ctsw, "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/B_cell注释/细胞通讯BT/KS的代码重新分析画图/pre_mpr.cell-cell_communications_CD8_ctsw.csv")

pre_mpr.bt.cellchat<- computeCommunProbPathway(pre_mpr.bt.cellchat)
pre_mpr.bt.cellchat <- aggregateNet(pre_mpr.bt.cellchat)

df_pre_mpr_gnly<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/B_cell注释/细胞通讯BT/KS的代码重新分析画图/pre_mpr.cell-cell_communications_CD8_C1-GNLY.csv")
df.net_gnly<- df_pre_mpr_gnly[,-1]
colnames(df.net_gnly)[5]<- "score"
head(df.net_gnly)

df_pre_mpr_ctsw<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/B_cell注释/细胞通讯BT/KS的代码重新分析画图/pre_mpr.cell-cell_communications_CD8_ctsw.csv")
df.net_ctsw<- df_pre_mpr_ctsw[,-1]
colnames(df.net_ctsw)[5]<- "score"
head(df.net_ctsw)

avg_expression<-AverageExpression(pre_nmpr_bt, group.by = "merged_group")

library(tidyr)
avg_expression_df <- as.data.frame(avg_expression$RNA)

avg_expression_df$gene <- rownames(avg_expression$RNA)

long_data <- pivot_longer(avg_expression_df, 
                          cols = -gene,  # 除了Gene列以外的所有列
                          names_to = "cell_type", 
                          values_to = "mean_exp")

head(long_data)
long_data <- long_data[, c(2, 1, 3:ncol(long_data))]
head(long_data)
cc_heatmap(df.net_gnly)
cc_heatmap(df.net_gnly, option = "B", n_top_ints = 10)

cc_heatmap(df.net_ctsw)
cc_heatmap(df.net_ctsw, option = "B", n_top_ints = 10)

dev.new()
cc_circos(df.net_gnly, option = "C", n_top_ints = 15, exp_df = long_data,
          cell_cols =c("CD8_C1-GNLY"= "hotpink", "Naive B cells" = "orange", "Naive plasma cells" = "cornflowerblue",
                       "Active B cells"="#8BE07A", "Plasma cells"="#82CDD4","Pre-B cells"="#D5D4BB","Primarily mature naive B cells" ="#AE66DD"), palette = "PuRd", cex = 0.75)
dev.off()


cc_circos(df.net_ctsw, option = "C", n_top_ints = 15, exp_df = long_data,
          cell_cols =c("CD8_C2-CTSW"= "hotpink", `Naive B cells` = "orange", `Naive plasma cells` = "cornflowerblue",
                       "Active B cells"="#8BE07A", "Plasma cells"="#82CDD4","Pre-B cells"="#D5D4BB","Primarily mature naive B cells" ="#AE66DD"), palette = "PuRd", cex = 0.75)



######KS 细胞通讯方法中央记忆细胞和耗竭性细胞----
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

library(CCPlotR)


T_cells_object<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/T_cells_object.rds")
T_cells_object = UpdateSeuratObject(T_cells_object)
T_cells_object

B_cells_object<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/B_cells_object_combined使用.rds")
B_cells_object = UpdateSeuratObject(B_cells_object)
B_cells_object

merged_obj <- merge(T_cells_object, B_cells_object, 
                    add.cell.ids = c("T_cells_object", "B_cells_object"), 
                    project = "Combined_Samples")

head(merged_obj@meta.data)

merged_obj@meta.data$merged_group <- paste(merged_obj@meta.data$cd4_cd8_group, 
                                           merged_obj@meta.data$B_cells_type, sep = "_")

# 查看更新后的 meta.data
head(merged_obj@meta.data)
merged_obj@meta.data$merged_group <- gsub("^NA_", "", merged_obj@meta.data$merged_group)
merged_obj@meta.data$merged_group <- gsub("^NA_", "", merged_obj@meta.data$merged_group)
merged_obj@meta.data$merged_group <- gsub("_NA$", "", merged_obj@meta.data$merged_group)

unique(merged_obj@meta.data$merged_group)

unique(merged_obj$MPR_Response2)

merged_obj_filter<- subset(merged_obj,subset=merged_group %in% c("Pre-B cells", "Primarily mature naive B cells", 
                                                                 "Naive B cells", "Active B cells", "Naive plasma cells", 
                                                                 "Plasma cells","CD8_C1-GNLY","CD8_C3-IL7R","CD8_C4-GZMK"))

pre_nmpr_bt <- subset(merged_obj_filter, subset =MPR_Response2 == 'Pre_NMPR')
unique(pre_nmpr_bt$MPR_Response2)
pre_mpr_bt <- subset(merged_obj_filter, subset =MPR_Response2 == 'Pre_MPR')
unique(pre_mpr_bt$MPR_Response2)
post_mpr_bt <- subset(merged_obj_filter, subset =MPR_Response2 == 'Post_MPR')
unique(post_mpr_bt$MPR_Response2)

pre_nmpr.bt.cellchat <- KS_cellchat(pre_nmpr_bt, assay = 'RNA', group.by = "merged_group",
                                    workers=2, species='human')

saveRDS(pre_nmpr.bt.cellchat, file = "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/中央记忆细胞和耗竭性细胞T细胞/细胞通讯分析/和B细胞通讯分析/pre_nmpr.bt.cellchat记忆和耗竭.rds")



pre_mpr.bt.cellchat <- KS_cellchat(pre_mpr_bt, assay = 'RNA', group.by = "merged_group",
                                   workers=2, species='human')

saveRDS(pre_mpr.bt.cellchat, file = "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/中央记忆细胞和耗竭性细胞T细胞/细胞通讯分析/和B细胞通讯分析/pre_mpr.bt.cellchat记忆和耗竭.rds")


post_mpr.bt.cellchat <- KS_cellchat(post_mpr_bt, assay = 'RNA', group.by = "merged_group",
                                    workers=2, species='human')

saveRDS(post_mpr.bt.cellchat, file = "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/中央记忆细胞和耗竭性细胞T细胞/细胞通讯分析/和B细胞通讯分析/post_mpr.bt.cellchat记忆和耗竭.rds")

####有无细菌----
B.T.bacteria_yes_chat_df<- subset(merged_obj_filter, subset=group_microbe2=="Bacteria+")
B.T.bacteria_yes_pre_nmpr<- subset(B.T.bacteria_yes_chat, subset=MPR_Response2=="Pre_NMPR")
unique(B.T.bacteria_yes_pre_nmpr$MPR_Response2)
B.T.bacteria_yes_pre_mpr<- subset(B.T.bacteria_yes_chat, subset=MPR_Response2=="Pre_MPR")
unique(B.T.bacteria_yes_pre_mpr$MPR_Response2)
B.T.bacteria_yes_post_mpr<- subset(B.T.bacteria_yes_chat, subset=MPR_Response2=="Post_MPR")
unique(B.T.bacteria_yes_post_mpr$MPR_Response2)

B.T.bacteria_yes_chat<- KS_cellchat(B.T.bacteria_yes_chat, assay = 'RNA', group.by = "merged_group",
                                       workers=2, species='human')

saveRDS(B.T.bacteria_yes_chat,
        file = "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/中央记忆细胞和耗竭性细胞T细胞/细胞通讯分析/和B细胞通讯分析/B.T.bacteria_yes_chat.rds")

B.T.bacteria_yes_pre_nmpr_chat<- KS_cellchat(B.T.bacteria_yes_pre_nmpr, assay = 'RNA', group.by = "merged_group",
                                                workers=2, species='human')

B.T.bacteria_yes_pre_mpr_chat<- KS_cellchat(B.T.bacteria_yes_pre_mpr, assay = 'RNA', group.by = "merged_group",
                                               workers=2, species='human')

B.T.bacteria_yes_post_mpr_chat<- KS_cellchat(B.T.bacteria_yes_post_mpr, assay = 'RNA', group.by = "merged_group",
                                                workers=2, species='human')


saveRDS(B.T.bacteria_yes_pre_nmpr_chat,
        file = "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/中央记忆细胞和耗竭性细胞T细胞/细胞通讯分析/和B细胞通讯分析/B.T.bacteria_yes_pre_nmpr_chat.rds")

saveRDS(B.T.bacteria_yes_pre_mpr_chat,
        file = "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/中央记忆细胞和耗竭性细胞T细胞/细胞通讯分析/和B细胞通讯分析/B.T.bacteria_yes_pre_mpr_chat.rds")

saveRDS(B.T.bacteria_yes_post_mpr_chat,
        file = "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/中央记忆细胞和耗竭性细胞T细胞/细胞通讯分析/和B细胞通讯分析/B.T.bacteria_yes_post_mpr_chat.rds")



B.T.bacteria_no_chat_df<- subset(merged_obj_filter, subset=group_microbe2=="Bacteria-")
B.T.bacteria_no_pre_nmpr<- subset(B.T.bacteria_no_chat, subset=MPR_Response2=="Pre_NMPR")
unique(B.T.bacteria_no_pre_nmpr$MPR_Response2)
B.T.bacteria_no_pre_mpr<- subset(B.T.bacteria_no_chat, subset=MPR_Response2=="Pre_MPR")
unique(B.T.bacteria_no_pre_mpr$MPR_Response2)
B.T.bacteria_no_post_mpr<- subset(B.T.bacteria_no_chat, subset=MPR_Response2=="Post_MPR")
unique(B.T.bacteria_no_post_mpr$MPR_Response2)

options(future.globals.maxSize = 3e+09) 
B.T.bacteria_no_chat<- KS_cellchat(B.T.bacteria_no_chat, assay = 'RNA', group.by = "merged_group",
                                      workers=2, species='human')

saveRDS(B.T.bacteria_no_chat,
        file = "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/中央记忆细胞和耗竭性细胞T细胞/细胞通讯分析/和B细胞通讯分析/B.T.bacteria_no_chat.rds")


B.T.bacteria_no_pre_nmpr_chat<- KS_cellchat(B.T.bacteria_no_pre_nmpr, assay = 'RNA', group.by = "merged_group",
                                               workers=2, species='human')

B.T.bacteria_no_pre_mpr_chat<- KS_cellchat(B.T.bacteria_no_pre_mpr, assay = 'RNA', group.by = "merged_group",
                                              workers=2, species='human')

B.T.bacteria_no_post_mpr_chat<- KS_cellchat(B.T.bacteria_no_post_mpr, assay = 'RNA', group.by = "merged_group",
                                               workers=2, species='human')


saveRDS(B.T.bacteria_no_pre_nmpr_chat,
        file = "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/中央记忆细胞和耗竭性细胞T细胞/细胞通讯分析/和B细胞通讯分析/B.T.bacteria_no_pre_nmpr_chat.rds")

saveRDS(B.T.bacteria_no_pre_mpr_chat,
        file = "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/中央记忆细胞和耗竭性细胞T细胞/细胞通讯分析/和B细胞通讯分析/B.T.bacteria_no_pre_mpr_chat.rds")

saveRDS(B.T.bacteria_no_post_mpr_chat,
        file = "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/中央记忆细胞和耗竭性细胞T细胞/细胞通讯分析/和B细胞通讯分析/B.T.bacteria_no_post_mpr_chat.rds")



####cellchat多组比较分析----

pre_nmpr.bt.cellchat<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/中央记忆细胞和耗竭性细胞T细胞/细胞通讯分析/和B细胞通讯分析/pre_nmpr.bt.cellchat记忆和耗竭.rds")
pre_mpr.bt.cellchat<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/中央记忆细胞和耗竭性细胞T细胞/细胞通讯分析/和B细胞通讯分析/pre_mpr.bt.cellchat记忆和耗竭.rds")
post_mpr.bt.cellchat<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/中央记忆细胞和耗竭性细胞T细胞/细胞通讯分析/和B细胞通讯分析/post_mpr.bt.cellchat记忆和耗竭.rds")





object.list <- list(bacteria_yes= B.T.bacteria_yes_chat, bacteria_no= B.T.bacteria_no_chat)

object.list <- list(Pre_NMPR= pre_nmpr.bt.cellchat, Pre_MPR= pre_mpr.bt.cellchat,Post_MPR=post_mpr.bt.cellchat)

object.list <- list(Pre_NMPR= B.T.bacteria_yes_pre_nmpr_chat, 
                    Pre_MPR= B.T.bacteria_yes_pre_mpr_chat,
                    Post_MPR=B.T.bacteria_yes_post_mpr_chat)

object.list <- list(Pre_NMPR= B.T.bacteria_no_pre_nmpr_chat, 
                    Pre_MPR= B.T.bacteria_no_pre_mpr_chat,
                    Post_MPR=B.T.bacteria_no_post_mpr_chat)



object.list <- list(Pre_NMPR= pre_nmpr.bt.cellchat, Pre_MPR= pre_mpr.bt.cellchat,Post_MPR=post_mpr.bt.cellchat)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2,3))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2,3), measure = "weight")
p1<-gg1 + gg2

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/中央记忆细胞和耗竭性细胞T细胞/细胞通讯分析/和B细胞通讯分析/分组比较柱状图记忆和耗竭有无细菌.pdf",p1,width=6,height=4)


gg1 <- netVisual_heatmap(cellchat)
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
pdf(file = "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/中央记忆细胞和耗竭性细胞T细胞/细胞通讯分析/和B细胞通讯分析/分组比较热图记忆和耗竭有无细菌.pdf", width = 10, height = 6)
gg1 + gg2
dev.off()

gg1 <- netVisual_heatmap(B.T.bacteria_yes_chat)
gg2 <- netVisual_heatmap(B.T.bacteria_yes_chat, measure = "weight")
gg1 + gg2

gg1 <- netVisual_heatmap(B.T.bacteria_no_chat)
gg2 <- netVisual_heatmap(B.T.bacteria_no_chat, measure = "weight")
gg1 + gg2



groupSize<-as.numeric(table(pre_nmpr.bt.cellchat@idents))
p1<-netVisual_circle(pre_nmpr.bt.cellchat@net$weight, vertex.weight=groupSize, sources.use = c("CD8_C3-IL7R","CD8_C4-GZMK","CD8_C1-GNLY"), targets.use =c("Pre-B cells", "Primarily mature naive B cells", 
                                                                                                                          "Naive B cells", "Active B cells", "Naive plasma cells", 
                                                                                                                          "Plasma cells"),weight.scale=T, label.edge=F, title.name="Interaction weights/strength")
p1
groupSize<-as.numeric(table(pre_mpr.bt.cellchat@idents))
p2<-netVisual_circle(pre_mpr.bt.cellchat@net$weight, vertex.weight=groupSize, sources.use = c("CD8_C3-IL7R","CD8_C4-GZMK","CD8_C1-GNLY"), targets.use =c("Pre-B cells", "Primarily mature naive B cells", 
                                                                                                                          "Naive B cells", "Active B cells", "Naive plasma cells", 
                                                                                                                          "Plasma cells"),weight.scale=T, label.edge=F, title.name="Interaction weights/strength")
p2


groupSize<-as.numeric(table(post_mpr.bt.cellchat@idents))
p3<-netVisual_circle(post_mpr.bt.cellchat@net$weight, vertex.weight=groupSize, sources.use =  c("CD8_C3-IL7R","CD8_C4-GZMK","CD8_C1-GNLY"), targets.use =c("Pre-B cells", "Primarily mature naive B cells", 
                                                                                                                         "Naive B cells", "Active B cells", "Naive plasma cells", 
                                                                                                                         "Plasma cells"),weight.scale=T, label.edge=F, title.name="Interaction weights/strength")
p3



groupSize<-as.numeric(table(pre_nmpr.bt.cellchat@idents))
p4<-netVisual_circle(pre_nmpr.bt.cellchat@net$weight, vertex.weight=groupSize, sources.use =c("Pre-B cells", "Primarily mature naive B cells", 
                                                                                                "Naive B cells", "Active B cells", "Naive plasma cells", 
                                                                                                "Plasma cells"), targets.use =c("CD8_C3-IL7R","CD8_C4-GZMK","CD8_C1-GNLY"),weight.scale=T, label.edge=F, title.name="Interaction weights/strength")
p4

groupSize<-as.numeric(table(pre_mpr.bt.cellchat@idents))
p5<-netVisual_circle(pre_mpr.bt.cellchat@net$weight, vertex.weight=groupSize, sources.use = c("Pre-B cells", "Primarily mature naive B cells", 
                                                                                              "Naive B cells", "Active B cells", "Naive plasma cells", 
                                                                                              "Plasma cells"), targets.use =c("CD8_C3-IL7R","CD8_C4-GZMK","CD8_C1-GNLY"),weight.scale=T, label.edge=F, title.name="Interaction weights/strength")
p5


groupSize<-as.numeric(table(post_mpr.bt.cellchat@idents))
p6<-netVisual_circle(post_mpr.bt.cellchat@net$weight, vertex.weight=groupSize, sources.use = c("Pre-B cells", "Primarily mature naive B cells", 
                                                                                               "Naive B cells", "Active B cells", "Naive plasma cells", 
                                                                                               "Plasma cells"), targets.use =c("CD8_C3-IL7R","CD8_C4-GZMK","CD8_C1-GNLY"),weight.scale=T, label.edge=F, title.name="Interaction weights/strength")
p6




df.net<-subsetCommunication(cellchat)
mat<-cellchat@net$Pre_NMPR$weight

Pre_NMPR<-pheatmap(as.matrix(mat),
         cluster_rows=F, 
         cluster_cols=F, 
         show_colnames=TRUE, 
         main="Heatmap of Interaction Strength")

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/中央记忆细胞和耗竭性细胞T细胞/细胞通讯分析/和B细胞通讯分析/Pre_NMPR_heatmap.pdf",Pre_NMPR,width = 7,height = 7)

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
ppp<-bacteria_yes+bacteria_no


ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/中央记忆细胞和耗竭性细胞T细胞/细胞通讯分析/和B细胞通讯分析/bacteria_yes_bacteria_no_heatmap.pdf",ppp,width = 8,height = 8)

num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
num.link <- unlist(num.link)
num.link <- as.numeric(unlist(num.link))
weight.MinMax <- c(min(num.link), max(num.link))
#weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  # 计算网络中心性得分
  object.list[[i]] <- netAnalysis_computeCentrality(object.list[[i]], slot.name = "netP") 
  # 然后生成散点图
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}
patchwork::wrap_plots(plots = gg)

p_combined <- patchwork::wrap_plots(plots = gg)
ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/中央记忆细胞和耗竭性细胞T细胞/细胞通讯分析/和B细胞通讯分析/combined_plots记忆和耗竭无细菌.pdf", p_combined, width = 10, height = 4) 

treatment_color <- c("#4974a4","#4dae47","#f29600")


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
pathway.union



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
pathway.union


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
pathway.union




pbubble<-netVisual_bubble(cellchat, sources.use = "CD8_C3-IL7R", targets.use =c("Pre-B cells", "Primarily mature naive B cells", 
                                                                                "Naive B cells", "Active B cells", "Naive plasma cells", 
                                                                                "Plasma cells"),comparison = c(1,2,3), angle.x = 45,font.size = 12,color.text =c("#4974a4","#4dae47","#f29600"))
?netVisual_bubble

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/中央记忆细胞和耗竭性细胞T细胞/细胞通讯分析/和B细胞通讯分析/IL7R_受配体_combined_pbubble.pdf", pbubble, width = 9, height = 9) 

pbubble<-netVisual_bubble(cellchat, sources.use = "CD8_C4-GZMK", targets.use =c("Pre-B cells", "Primarily mature naive B cells", 
                                                                                "Naive B cells", "Active B cells", "Naive plasma cells", 
                                                                                "Plasma cells"),comparison = c(1,2,3), angle.x = 45,font.size = 12,color.text =c("#4974a4","#4dae47","#f29600"))
?netVisual_bubble

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/中央记忆细胞和耗竭性细胞T细胞/细胞通讯分析/和B细胞通讯分析/GZMK_受配体_combined_pbubble.pdf", pbubble, width = 9, height = 9) 

#"CD8_C1-GNLY","CD8_C2-CTSW"

pbubble<-netVisual_bubble(cellchat, sources.use = "CD8_C1-GNLY", targets.use =c("Pre-B cells", "Primarily mature naive B cells", 
                                                                                "Naive B cells", "Active B cells", "Naive plasma cells", 
                                                                                "Plasma cells"),comparison = c(1,2,3), angle.x = 45,font.size = 12,color.text =c("#4974a4","#4dae47","#f29600"))
?netVisual_bubble

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/中央记忆细胞和耗竭性细胞T细胞/细胞通讯分析/和B细胞通讯分析/GNLY_受配体_combined_pbubble.pdf", pbubble, width = 9, height = 9) 

pbubble<-netVisual_bubble(cellchat, sources.use = "CD8_C2-CTSW", targets.use =c("Pre-B cells", "Primarily mature naive B cells", 
                                                                                "Naive B cells", "Active B cells", "Naive plasma cells", 
                                                                                "Plasma cells"),comparison = c(1,2,3), angle.x = 45,font.size = 12,color.text =c("#4974a4","#4dae47","#f29600"))
?netVisual_bubble

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/中央记忆细胞和耗竭性细胞T细胞/细胞通讯分析/和B细胞通讯分析/CTSW_受配体_combined_pbubble.pdf", pbubble, width = 9, height = 9) 




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


options(future.globals.maxSize = 5e9)
merged_obj_filter.bt.cellchat <- KS_cellchat(merged_obj_filter, assay = 'RNA', group.by = "merged_group",
                                             workers=2, species='human')

saveRDS(merged_obj_filter.bt.cellchat, file = "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/中央记忆细胞和耗竭性细胞T细胞/细胞通讯分析/和B细胞通讯分析/merged_obj_filter.bt.cellchat记忆和耗竭.rds")

df.net <- subsetCommunication(merged_obj_filter.bt.cellchat)
write.csv(df.net, "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/中央记忆细胞和耗竭性细胞T细胞/细胞通讯分析/和B细胞通讯分析/cell-cell_communications_ligands_receptors记忆和耗竭.csv")
df.net1 <- subsetCommunication(merged_obj_filter.bt.cellchat,slot.name = "netP")
write.csv(df.net1, "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/中央记忆细胞和耗竭性细胞T细胞/细胞通讯分析/和B细胞通讯分析/cell-cell_communications_signaling_pathways记忆和耗竭.csv")
levels(merged_obj_filter.bt.cellchat@idents)
df.net2 <- subsetCommunication(merged_obj_filter.bt.cellchat, sources.use = c("CD8_C3-IL7R"), targets.use = c("Naive B cells" ,             
                                                                                                              "Naive plasma cells","Active B cells" ,"Plasma cells",                  
                                                                                                              "Pre-B cells","Primarily mature naive B cells")) 
write.csv(df.net2, "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/中央记忆细胞和耗竭性细胞T细胞/细胞通讯分析/和B细胞通讯分析/cell-cell_communications_CD8_C3-IL7R记忆和耗竭.csv")



####B.T.bacteria_yes_chat, bacteria_no= B.T.bacteria_no_chat
df.net_GZMK <- subsetCommunication(B.T.bacteria_yes_chat, sources.use = c("CD8_C4-GZMK"), targets.use = c("Naive B cells" ,             
                                                                                                                  "Naive plasma cells","Active B cells" ,"Plasma cells", "Pre-B cells","Primarily mature naive B cells")) 
write.csv(df.net_GZMK, "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/中央记忆细胞和耗竭性细胞T细胞/细胞通讯分析/和B细胞通讯分析/cell-cell_communications_CD8_GZMK记忆和耗竭有细菌存在.csv")

df.net_IL7R <- subsetCommunication(B.T.bacteria_yes_chat, sources.use = c("CD8_C3-IL7R"), targets.use = c("Naive B cells" ,             
                                                                                                                  "Naive plasma cells","Active B cells" ,"Plasma cells",  "Pre-B cells","Primarily mature naive B cells")) 

write.csv(df.net_IL7R, "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/中央记忆细胞和耗竭性细胞T细胞/细胞通讯分析/和B细胞通讯分析/cell-cell_communications_CD8_IL7R记忆和耗竭有细菌存在.csv")

df.net_GNLY <- subsetCommunication(B.T.bacteria_yes_chat, sources.use = c("CD8_C1-GNLY"), targets.use = c("Naive B cells" ,             
                                                                                                                    "Naive plasma cells","Active B cells" ,"Plasma cells",  "Pre-B cells","Primarily mature naive B cells")) 

write.csv(df.net_GNLY, "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/中央记忆细胞和耗竭性细胞T细胞/细胞通讯分析/和B细胞通讯分析/cell-cell_communications_CD8_GNLY记忆和耗竭有细菌存在.csv")

merged_obj_filter.bt.cellchat<- computeCommunProbPathway(merged_obj_filter.bt.cellchat)
merged_obj_filter.bt.cellchat <- aggregateNet(merged_obj_filter.bt.cellchat)

df_pre_nmpr_IL7R<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/中央记忆细胞和耗竭性细胞T细胞/细胞通讯分析/和B细胞通讯分析/cell-cell_communications_CD8_C3-IL7R记忆和耗竭.csv")
df.net_IL7R<- df_pre_nmpr_IL7R[,-1]
colnames(df.net_IL7R)[5]<- "score"
head(df.net_IL7R)

df_pre_nmpr_GZMK<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/中央记忆细胞和耗竭性细胞T细胞/细胞通讯分析/和B细胞通讯分析/cell-cell_communications_CD8_GZMK记忆和耗竭.csv")
df.net_GZMK<- df_pre_nmpr_GZMK[,-1]
colnames(df.net_GZMK)[5]<- "score"
head(df.net_GZMK)

df_pre_nmpr_GNLY<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/中央记忆细胞和耗竭性细胞T细胞/细胞通讯分析/和B细胞通讯分析/cell-cell_communications_CD8_GNLY记忆和耗竭.csv")
df.net_GNLY<- df_pre_nmpr_GNLY[,-1]
colnames(df.net_GNLY)[5]<- "score"
head(df.net_GNLY)

avg_expression<-AverageExpression(B.T.bacteria_yes_chat_df, group.by = "merged_group")

library(tidyr)
avg_expression_df <- as.data.frame(avg_expression$RNA)

avg_expression_df$gene <- rownames(avg_expression$RNA)

long_data <- pivot_longer(avg_expression_df, 
                          cols = -gene,  # 除了Gene列以外的所有列
                          names_to = "cell_type", 
                          values_to = "mean_exp")

head(long_data)
long_data <- long_data[, c(2, 1, 3:ncol(long_data))]
head(long_data)
cc_heatmap(df.net_IL7R)
cc_heatmap(df.net_IL7R, option = "B", n_top_ints = 10)

cc_heatmap(df.net_GZMK)
cc_heatmap(df.net_GZMK, option = "B", n_top_ints = 10)

dev.new()

pdf("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/中央记忆细胞和耗竭性细胞T细胞/细胞通讯分析/和B细胞通讯分析/CD8_C3-IL7R_Bcell记忆和耗竭.pdf",width=7,height=7)

cc_circos(df.net_IL7R, option = "C", n_top_ints = 15, exp_df = long_data,
          cell_cols =c("CD8_C3-IL7R"= "hotpink", "Naive B cells" = "orange", "Naive plasma cells" = "cornflowerblue",
                       "Active B cells"="#8BE07A", "Plasma cells"="#82CDD4","Pre-B cells"="#D5D4BB","Primarily mature naive B cells" ="#AE66DD"), palette = "PuRd", cex = 0.75)
dev.off()

pdf("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/中央记忆细胞和耗竭性细胞T细胞/细胞通讯分析/和B细胞通讯分析/CD8_C4-GZMK_Bcell记忆和耗竭.pdf",width=7,height=7)

cc_circos(df.net_GZMK, option = "C", n_top_ints = 15, exp_df = long_data,
          cell_cols =c("CD8_C4-GZMK"= "hotpink", `Naive B cells` = "orange", `Naive plasma cells` = "cornflowerblue",
                       "Active B cells"="#8BE07A", "Plasma cells"="#82CDD4","Pre-B cells"="#D5D4BB","Primarily mature naive B cells" ="#AE66DD"), palette = "PuRd", cex = 0.75)

dev.off()


pdf("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/中央记忆细胞和耗竭性细胞T细胞/细胞通讯分析/和B细胞通讯分析/CD8_C5-STMN1_Bcell记忆和耗竭.pdf",width=7,height=7)

cc_circos(df.net_STMN1, option = "C", n_top_ints = 15, exp_df = long_data,
          cell_cols =c("CD8_C5-STMN1"= "hotpink", "Naive B cells" = "orange", "Naive plasma cells" = "cornflowerblue",
                       "Active B cells"="#8BE07A", "Plasma cells"="#82CDD4","Pre-B cells"="#D5D4BB","Primarily mature naive B cells" ="#AE66DD"), palette = "PuRd", cex = 0.75)
dev.off()

pdf("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/中央记忆细胞和耗竭性细胞T细胞/细胞通讯分析/和B细胞通讯分析/CD8_C1-GNLY_Bcell记忆和耗竭.pdf",width=7,height=7)

cc_circos(df.net_GNLY, option = "C", n_top_ints = 15, exp_df = long_data,
          cell_cols =c("CD8_C1-GNLY"= "hotpink", `Naive B cells` = "orange", `Naive plasma cells` = "cornflowerblue",
                       "Active B cells"="#8BE07A", "Plasma cells"="#82CDD4","Pre-B cells"="#D5D4BB","Primarily mature naive B cells" ="#AE66DD"), palette = "PuRd", cex = 0.75)

dev.off()


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

T_cells_object<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/T_cells_object.rds")
T_cells_object = UpdateSeuratObject(T_cells_object)
T_cells_object

B_cells_object<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/B_cells_object_combined使用.rds")
B_cells_object = UpdateSeuratObject(B_cells_object)
B_cells_object

merged_obj <- merge(T_cells_object, B_cells_object, 
                    add.cell.ids = c("T_cells_object", "B_cells_object"), 
                    project = "Combined_Samples")

head(merged_obj@meta.data)

merged_obj@meta.data$merged_group <- paste(merged_obj@meta.data$cd4_cd8_group, 
                                           merged_obj@meta.data$B_cells_type, sep = "_")

# 查看更新后的 meta.data
head(merged_obj@meta.data)
merged_obj@meta.data$merged_group <- gsub("^NA_", "", merged_obj@meta.data$merged_group)
merged_obj@meta.data$merged_group <- gsub("^NA_", "", merged_obj@meta.data$merged_group)
merged_obj@meta.data$merged_group <- gsub("_NA$", "", merged_obj@meta.data$merged_group)

unique(merged_obj@meta.data$merged_group)

unique(merged_obj$MPR_Response2)

merged_obj_filter<- subset(merged_obj,subset=merged_group %in% c("Pre-B cells", "Primarily mature naive B cells", 
                                                                 "Naive B cells", "Active B cells", "Naive plasma cells", 
                                                                 "Plasma cells","CD8_C1-GNLY","CD8_C3-IL7R","CD8_C4-GZMK"))

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
sender_celltypes <- c("Pre-B cells", "Primarily mature naive B cells", 
                      "Naive B cells", "Active B cells", "Naive plasma cells", 
                      "Plasma cells")
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
####"CD8_C1-GNLY","CD8_C3-IL7R","CD8_C4-GZMK"


receiver="CD8_C1-GNLY"  
sender=c("Pre-B cells", "Primarily mature naive B cells", 
         "Naive B cells", "Active B cells", "Naive plasma cells", 
         "Plasma cells")

sender=c("ZBTB16+TAMs")

receiver="CD8_C1-GNLY"

sender=c("Pre-B cells")

unique(B.T.bacteria_no_chat_df$merged_group)
##B.T.bacteria_no_chat_df
B.T.bacteria_no_chat_df$cell_cluster = ifelse(B.T.bacteria_no_chat_df$merged_group==receiver, receiver,"other_cells")
DEG_Fib<-FindMarkers(B.T.bacteria_no_chat_df, ident.1=receiver, ident.2="other_cells", only.pos=TRUE, logfc.threshold=0.2, group.by = "cell_cluster")
expressed_genes_receiver = get_expressed_genes(receiver, B.T.bacteria_no_chat_df, pct = 0.1) %>% .[. %in% rownames(DEG_Fib)] 

#sender expressed_genes
list_expressed_genes_sender <- sender %>% unique() %>% lapply(get_expressed_genes,B.T.bacteria_no_chat_df, 0.1)
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
                                                        celltype_col = "merged_group",
                                                        sender = sender) 

receiver="CD8_C1-GNLY"  
sender=c("Pre-B cells", "Primarily mature naive B cells", 
         "Naive B cells", "Active B cells", "Naive plasma cells", 
         "Plasma cells")

#Define the ligand-target links of interest
active_ligand_target_links_df$target_type <- "CD8_C1-GNLY" # needed for joining tables
circos_links <- get_ligand_target_links_oi(ligand_type_indication_df,
                                           active_ligand_target_links_df,
                                           cutoff = 0) 


unique(circos_links$ligand_type)
unique(circos_links$target_type)

ligand_colors <- c("General" = "#377EB8", "Naive B cells" = "#C453DC", "Pre-B cells" = "#00ffff", "Primarily mature naive B cells" = "#0000FF",
                   "Active B cells" = "#8EDCC1", "Naive plasma cells" = "#7fff00","Plasma cells"="#ffd700") 

target_colors <- c("CD8_C1-GNLY" = "#999999") 

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
  #rename(ligand=from, target=to) %>% 
  mutate(target_type = "CD8_C1-GNLY") %>% 
  inner_join(ligand_type_indication_df)

receptor_colors <- c("CD8_C1-GNLY" = "#E41A1C")

vis_circos_receptor_obj <- prepare_circos_visualization(lr_network_top_df,
                                                        ligand_colors = ligand_colors,
                                                        target_colors = receptor_colors) 


make_circos_plot(vis_circos_receptor_obj, transparency = T,link.visible = TRUE,  args.circos.text = list(cex = 0.8)) 







##################################
df.net <- subsetCommunication(pre_mpr.bt.cellchat)
write.csv(df.net, "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/组织驻留T细胞/细胞通讯分析/和B细胞通讯分析/pre_mpr.cell-cell_communications_ligands_receptors.csv")

df.net1 <- subsetCommunication(pre_mpr.bt.cellchat,slot.name = "netP")
write.csv(df.net1, "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/组织驻留T细胞/细胞通讯分析/和B细胞通讯分析/pre_mpr.cell-cell_communications_signaling_pathways.csv")
levels(pre_mpr.bt.cellchat)
df.net2 <- subsetCommunication(pre_mpr.bt.cellchat, sources.use = c("CD8_C3-IL7R"), targets.use = c("Naive B cells" ,             
                                                                                                    "Naive plasma cells","Active B cells" ,"Plasma cells",                  
                                                                                                    "Pre-B cells","Primarily mature naive B cells")) 
write.csv(df.net2, "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/组织驻留T细胞/细胞通讯分析/和B细胞通讯分析/pre_mpr.cell-cell_communications_CD8_C3-IL7R.csv")

df.net_GZMK <- subsetCommunication(pre_mpr.bt.cellchat, sources.use = c("CD8_C4-GZMK"), targets.use = c("Naive B cells" ,             
                                                                                                        "Naive plasma cells","Active B cells" ,"Plasma cells",                  
                                                                                                        "Pre-B cells","Primarily mature naive B cells")) 
write.csv(df.net_GZMK, "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/组织驻留T细胞/细胞通讯分析/和B细胞通讯分析/pre_mpr.cell-cell_communications_CD8_C4-GZMK.csv")

pre_mpr.bt.cellchat<- computeCommunProbPathway(pre_mpr.bt.cellchat)
pre_mpr.bt.cellchat <- aggregateNet(pre_mpr.bt.cellchat)

df_pre_mpr_IL7R<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/组织驻留T细胞/细胞通讯分析/和B细胞通讯分析/pre_mpr.cell-cell_communications_CD8_C3-IL7R.csv")
df.net_IL7R<- df_pre_mpr_IL7R[,-1]
colnames(df.net_IL7R)[5]<- "score"
head(df.net_IL7R)

df_pre_mpr_GZMK<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/组织驻留T细胞/细胞通讯分析/和B细胞通讯分析/pre_mpr.cell-cell_communications_CD8_C4-GZMK.csv")
df.net_GZMK<- df_pre_mpr_GZMK[,-1]
colnames(df.net_GZMK)[5]<- "score"
head(df.net_GZMK)

avg_expression<-AverageExpression(pre_nmpr_bt, group.by = "merged_group")

library(tidyr)
avg_expression_df <- as.data.frame(avg_expression$RNA)

avg_expression_df$gene <- rownames(avg_expression$RNA)

long_data <- pivot_longer(avg_expression_df, 
                          cols = -gene,  # 除了Gene列以外的所有列
                          names_to = "cell_type", 
                          values_to = "mean_exp")

head(long_data)
long_data <- long_data[, c(2, 1, 3:ncol(long_data))]
head(long_data)
cc_heatmap(df.net_IL7R)

pdf("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/组织驻留T细胞/细胞通讯分析/和B细胞通讯分析/IL7R_heatmap.pdf",width=7,height=7)
cc_heatmap(df.net_IL7R, option = "B", n_top_ints = 10)
dev.off()

pdf("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/组织驻留T细胞/细胞通讯分析/和B细胞通讯分析/GZMK_heatmap.pdf",width=7,height=7)
cc_heatmap(df.net_GZMK, option = "B", n_top_ints = 10)
dev.off()

# 
# 
# pdf("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/组织驻留T细胞/细胞通讯分析/和B细胞通讯分析/GZMK_heatmap.pdf",width=7,height=7)
# cc_circos(df.net_IL7R, option = "C", n_top_ints = 15, exp_df = long_data,
#           cell_cols =c("CD8_C1-IL7R"= "hotpink", "Naive B cells" = "orange", "Naive plasma cells" = "cornflowerblue",
#                        "Active B cells"="#8BE07A", "Plasma cells"="#82CDD4","Pre-B cells"="#D5D4BB","Primarily mature naive B cells" ="#AE66DD"), palette = "PuRd", cex = 0.75)
# dev.off()
# 
# 
# cc_circos(df.net_GZMK, option = "C", n_top_ints = 15, exp_df = long_data,
#           cell_cols =c("CD8_C2-GZMK"= "hotpink", `Naive B cells` = "orange", `Naive plasma cells` = "cornflowerblue",
#                        "Active B cells"="#8BE07A", "Plasma cells"="#82CDD4","Pre-B cells"="#D5D4BB","Primarily mature naive B cells" ="#AE66DD"), palette = "PuRd", cex = 0.75)
# 




