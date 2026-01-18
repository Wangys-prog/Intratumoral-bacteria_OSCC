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

setwd("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/成纤维细胞/")

MPRNMPR_object_miMPRobe_remove<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/MPRNMPR_object_miMPRobe_remove.rds")
MPRNMPR_object_miMPRobe_remove = UpdateSeuratObject(MPRNMPR_object_miMPRobe_remove)
MPRNMPR_object_miMPRobe_remove


CAFs_cells_object<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/CAFs_cells_object.rds")
CAFs_cells_object = UpdateSeuratObject(CAFs_cells_object)
CAFs_cells_object


setwd("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/成纤维细胞")

####提取myeloid cells subject----
Idents(MPRNMPR_object_miMPRobe)<-"cell_type_new"

Idents(myeloid_cells_object) <-"RNA_snn_res.0.4"
levels(Idents(myeloid_cells_object))

cluster9_17<- subset(myeloid_cells_object,idents = c("9","17"))

Fibroblasts_cells_object=subset(MPRNMPR_object_miMPRobe,idents = c("Fibroblasts"))

CAFs_cells_object = merge(Fibroblasts_cells_object, y = cluster9_17, project = "CombinedProject")

#CAFs_cells_object <-Fibroblasts_cells_object
table(Idents(CAFs_cells_object))

colnames(CAFs_cells_object@meta.data)

pbmc<- CAFs_cells_object

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

CAFs_cells_object <- pbmc
Idents(CAFs_cells_object) <-"RNA_snn_res.0.4"
table(Idents(CAFs_cells_object))

markers <- FindAllMarkers(CAFs_cells_object, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
all.markers =markers %>% dplyr::select(gene, everything()) %>% subset(p_val<0.05)
top20 =all.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
p_heatmap<-DoHeatmap(CAFs_cells_object,features = top20$gene)+ #热图复现
  theme(text = element_text(size = 5))

write.csv(all.markers,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/成纤维细胞/CAFs_cells_marker_all.markers去除.csv")
write.csv(top20,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/成纤维细胞/CAFs_cells_marker_top20去除.csv")

#保存RDS
saveRDS(CAFs_cells_object,"D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/CAFs_cells_object.rds")

####marker基因小提琴图----
setwd("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/成纤维细胞/")
CAFs_cells_object<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/CAFs_cells_object.rds")

cell_marker<- read.csv("marker1.csv")
######删除找不到的基因
nofinder<- c("IL8", "SRFP4", "IGK1", "CTGF")

cell_marker_filtered <- cell_marker[!cell_marker$gene_name %in% nofinder, ]

Idents(CAFs_cells_object)<-"RNA_snn_res.0.4"
levels(Idents(CAFs_cells_object))
unique(cell_marker_filtered$gene_name)
library(randomcoloR)

my36colors = distinctColorPalette(180)

vln.dat=FetchData(CAFs_cells_object,c(unique(cell_marker_filtered$gene_name),"RNA_snn_res.0.4"))
colnames(vln.dat)
#宽转长
vln.dat.melt <- reshape2::melt(vln.dat, id.vars = c("RNA_snn_res.0.4"), 
                               measure.vars = unique(cell_marker_filtered$gene_name),
                               variable.name = "gene", 
                               value.name = "Expr") 

vln.dat.melt$gene<- factor(vln.dat.melt$gene,levels = unique(vln.dat.melt$gene))

View(vln.dat.melt)

library(cowplot)
vln.dat.melt$RNA_snn_res.0.4<- factor(vln.dat.melt$RNA_snn_res.0.4,levels = c("0","1","2","3","4","5",
                                                                              "6","7","8","9","10","11",
                                                                              "12","13","14","15","16"))

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

ggsave("成纤维细胞_marker小提琴图RNA_snn_res.0.4.jpg",p1,width = 10,height = 7)



markers <- FindAllMarkers(CAFs_cells_object, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
all.markers =markers %>% dplyr::select(gene, everything()) %>% subset(p_val<0.05)
top20 =all.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
Idents(CAFs_cells_object) <- CAFs_cells_object$RNA_snn_res.0.4

write.csv(all.markers,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/成纤维细胞/CAF_cells_marker_all.markers.csv")

write.csv(top20,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/成纤维细胞/CAF_cells_marker_top20.csv")

DoHeatmap(CAFs_cells_object, features = unique(cell_marker_filtered$gene_name), size = 3)

FeaturePlot(object=CAFs_cells_object,reduction = "umap.harmony",raster=FALSE,features =c("IFN1")) 


VlnPlot(CAFs_cells_object, features = "PTPRC",raster=FALSE)

####注释----
old_to_new_groups <- c("1"="CXCL2+iCAFs",
                       "3"="CXCL2+iCAFs",
                       "4"="CXCL2+iCAFs",
                       "7"="CXCL2+iCAFs",
                       "11"="CXCL2+iCAFs", 
                       "9"="CXCR4+iCAFs",
                       "15"="TPSAB1+iCAFs",
                       "0"="MyCAFs",
                       "13"="MyCAFs",
                       "6"="MyCAFs",
                       "14"="MyCAFs",
                       "16"="apCAFs",
                       "8"="apCAFs",
                       "2"="apCAFs",
                       "5"="eCAFs",
                       "12"="eCAFs",
                       "10"="eCAFs")

CAFs_cells_object@meta.data$CAFs_cells_type <- plyr::mapvalues(CAFs_cells_object@meta.data$RNA_snn_res.0.4,
                                                                            from = names(old_to_new_groups), 
                                                                            to = old_to_new_groups)

table(CAFs_cells_object@meta.data$CAFs_cells_type)
unique(CAFs_cells_object@meta.data$CAFs_cells_type)

DimPlot(CAFs_cells_object, group.by = "CAFs_cells_type",reduction = "umap.harmony", pt.size=0.1)
Idents(CAFs_cells_object)<-"CAFs_cells_type"

DoHeatmap(CAFs_cells_object, features = unique(cell_marker_filtered$gene_name), size = 3)

DotPlot(CAFs_cells_object, features = unique(cell_marker_filtered$gene_name)) + RotatedAxis()

FeaturePlot(object=CAFs_cells_object,reduction = "umap.harmony",raster=FALSE,features =c("VEGFA"
)) 

saveRDS(CAFs_cells_object,"D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/CAFs_cells_object.rds")

####cluster 之后重新画图----

cell_marker<- read.csv("marker画图.csv")

nofinder<- c("IL8", "SRFP4", "IGK1", "CTGF")

cell_marker_filtered <- cell_marker[!cell_marker$gene_name %in% nofinder, ]

Idents(CAFs_cells_object)<-"CAFs_cells_type"
levels(Idents(CAFs_cells_object))
unique(cell_marker_filtered$gene_name)

library(tidyverse)
CAFs_cells.harmonyredu<- CAFs_cells_object@reductions$umap.harmony@cell.embeddings%>%
  as.data.frame() %>% 
  cbind(cell_type =CAFs_cells_object@meta.data$CAFs_cells_type)%>% 
  cbind(Patients = CAFs_cells_object@meta.data$Patients) %>% 
  cbind(Treatments = CAFs_cells_object@meta.data$Treatments) %>% 
  cbind(MPR_Response = CAFs_cells_object@meta.data$MPR_Response) %>%
  cbind(MPR_Response2 = CAFs_cells_object@meta.data$MPR_Response2) %>%
  cbind(cell_id = CAFs_cells_object@meta.data$cell_id) %>%
  cbind(orig.ident = CAFs_cells_object@meta.data$orig.ident)

names(CAFs_cells.harmonyredu)

colnames(CAFs_cells.harmonyredu) <- c("umapharmony_1","umapharmony_2","cell_type","Patients",
                                         "Treatments","MPR_Response", "MPR_Response2","cell_id",
                                         "orig.ident")

write.csv(CAFs_cells.harmonyredu,"CAFs_cells_umap降维数据.csv")

###颜色设置###画小坐标轴####https://cloud.tencent.com/developer/article/1924260
library(randomcoloR)
cell_type_color = distinctColorPalette(6)
cell_type_color<- c("#00bfff","#ff1493","#66CDAA","#FFa07a","#0000FF","#FF00ff")

unique(CAFs_cells.harmonyredu$cell_type)

CAFs_cells.harmonyredu$cell_type <- factor(CAFs_cells.harmonyredu$cell_type,levels= c("CXCL2+iCAFs","CXCR4+iCAFs","TPSAB1+iCAFs","MyCAFs","apCAFs","eCAFs"))
library(ggplot2)

CAFs_cell_type <- ggplot(CAFs_cells.harmonyredu,aes(x= umapharmony_1, y =umapharmony_2,color = cell_type)) +  
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
  geom_segment(aes(x = min(CAFs_cells.harmonyredu$umapharmony_1) , y = min(CAFs_cells.harmonyredu$umapharmony_2) ,
                   xend = min(CAFs_cells.harmonyredu$umapharmony_1) +3, yend = min(CAFs_cells.harmonyredu$umapharmony_2) ),
               colour = "black", size=1,arrow = arrow(length = unit(0.5,"cm")))+ 
  geom_segment(aes(x = min(CAFs_cells.harmonyredu$umapharmony_1)  , y = min(CAFs_cells.harmonyredu$umapharmony_2)  ,
                   xend = min(CAFs_cells.harmonyredu$umapharmony_1) , yend = min(CAFs_cells.harmonyredu$umapharmony_2) +3),
               colour = "black", size=1,arrow = arrow(length = unit(0.5,"cm"))) +
  annotate("text", x = min(CAFs_cells.harmonyredu$umapharmony_1) +1.5, y = min(CAFs_cells.harmonyredu$umapharmony_2)-1, label = "UMPA_1",
           color="black",size = 4, fontface="bold" ) + 
  annotate("text", x = min(CAFs_cells.harmonyredu$umapharmony_1)-1, y = min(CAFs_cells.harmonyredu$umapharmony_2) + 1.5, label = "UMAP_2",
           color="black",size = 4, fontface="bold" ,angle=90)

CAFs_cell_type

CAFs_cell_type_med <- CAFs_cells.harmonyredu %>%
  group_by(cell_type) %>%
  summarise(
    umapharmony_1 = median(umapharmony_1),
    umapharmony_2 = median(umapharmony_2)
  )

library(ggrepel)
CAFs_cell_type_med

CAFs_cell_type2<-CAFs_cell_type +geom_label_repel(aes(label=cell_type),size=4,color="black",fontface="bold",data = CAFs_cell_type_med,
                                                        point.padding=unit(0.5, "lines"),fill = alpha(c("white"),0.5),
                                                        segment.size=0.5, nudge_x=1, nudge_y = 0.5,direction="y",max.overlaps=50)+
  theme(legend.position = "none")


ggsave('./CAFs_cell_type.jpg',CAFs_cell_type2,width =6,height = 6)
ggsave('./CAFs_cell_type.pdf',CAFs_cell_type2,width =6,height = 6)


####dotplot数据画热图----
cell_marker<- read.csv("marker画图.csv")
p1 <- DotPlot(CAFs_cells_object, features =unique(cell_marker$gene_name),
              assay='RNA' )
p1
sub_markers_count<- p1$data
marker_group_list<-c("CXCL2+iCAFs","CXCR4+iCAFs","TPSAB1+iCAFs","MyCAFs","apCAFs","eCAFs")
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
        axis.text.x = element_text(size = 16,color = 'black',angle =40,hjust =1,vjust=1),
        legend.title = element_text(size=15, color = "black"),
        legend.text = element_text(size=15,color = "black",angle =0),
        legend.position = "right") +
  scale_y_discrete(position = "left")
p1_heat


ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/成纤维细胞/成纤维细胞细胞亚型分类热图.jpg",p1_heat,height =10,width = 7)

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/成纤维细胞/成纤维细胞细胞亚型分类热图.pdf",p1_heat,height =10,width = 7)


####
####成纤维细胞OR指数----

library("sscVis")
library("data.table")
library("grid")
library("cowplot")
library("ggrepel")
library("readr")
library("plyr")
library("ggpubr")
library("ggplot2")
library(dplyr)

#设置图片输出目录
CAFs_cells_object
out.prefix <- "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/成纤维细胞/OR指数/"
CAFs_cells_object<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/CAFs_cells_object.rds")
CAFs_cells_object$CAFs_cells_type
Idents(CAFs_cells_object)<-"CAFs_cells_type"

# 转换元数据为数据框
metadata <- as.data.frame(CAFs_cells_object@meta.data)
colnames(metadata)
unique(metadata$CAFs_cells_type)


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

CAFs_cells_object$CAFs_cells_type

meta.tb <-CAFs_cells_object@meta.data


A <- do.tissueDist(cellInfo.tb = meta.tb,
                   meta.cluster = meta.tb$CAFs_cells_type,
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
data$rid <- factor(data$rid, levels=c("CXCL2+iCAFs","CXCR4+iCAFs","TPSAB1+iCAFs","MyCAFs","apCAFs","eCAFs"))

cell_type_color<- c("#00bfff","#ff1493","#66CDAA","#FFa07a","#0000FF","#FF00ff")


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

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/成纤维细胞/OR指数和细胞比例/成纤维_OR_plot_MPR.pdf",OR_plot,width=9,height=5)



####髓系细胞比例计算----

CAFs_cells_object$CAFs_cells_type

Cellratio <- prop.table(table(CAFs_cells_object@meta.data$CAFs_cells_type,CAFs_cells_object@meta.data$orig.ident), margin = 2)

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

cell_type_groups = c("CXCL2+iCAFs","CXCR4+iCAFs","TPSAB1+iCAFs","MyCAFs","apCAFs","eCAFs")

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


#c("CXCL2+iCAFs","CXCR4+iCAFs","TPSAB1+iCAFs","MyCAFs","apCAFs","eCAFs")

cells_<-plot_grid(pplist[['CXCL2+iCAFs']],
                  pplist[['CXCR4+iCAFs']],
                  pplist[['TPSAB1+iCAFs']],
                  pplist[['MyCAFs']],
                  pplist[['apCAFs']],
                  pplist[['eCAFs']],
                  align = "h",  #axis = 'l',
                  nrow =2)

ggsave('D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/成纤维细胞/OR指数和细胞比例/CAF_cells_MPR_percentage.pdf',cells_,width =9,height =8,limitsize = FALSE)


####成纤维细胞与免疫细胞相关性热图--
library(tidyverse)
library(corrplot)

CAFs_cells_object$CAFs_cells_type

myeloid_cells_object=readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/myeloid_cells_object_filter使用.rds")
myeloid_cells_object = UpdateSeuratObject(myeloid_cells_object)
myeloid_cells_object

B_cells_object<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/B_cells_object_combined使用.rds")
B_cells_object = UpdateSeuratObject(B_cells_object)
B_cells_object

T_cells_object<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/T_cells_object.rds")
T_cells_object = UpdateSeuratObject(T_cells_object)
T_cells_object

merged_obj <- merge(CAFs_cells_object, y = c(T_cells_object,B_cells_object,myeloid_cells_object), add.cell.ids = c("CAFs_cells_object","T_cells_object","B_cells_object","myeloid_cells_object"))

View(merged_obj)
merged_obj$CAFs_cells_type
merged_obj$cd4_cd8_group
merged_obj$B_cells_type
merged_obj$myeloid_cells_type

merged_obj@meta.data$merged_group <- paste(merged_obj$CAFs_cells_type,
                                           merged_obj$cd4_cd8_group,
                                           merged_obj$B_cells_type,
                                           merged_obj$myeloid_cells_type,
                                           sep = "_")

# 查看更新后的 meta.data
head(merged_obj@meta.data)
merged_obj@meta.data$merged_group <- gsub("^NA_", "", merged_obj@meta.data$merged_group)
merged_obj@meta.data$merged_group <- gsub("_NA$", "", merged_obj@meta.data$merged_group)
merged_obj@meta.data$merged_group <- gsub("_NA_NA_[0-9]+$", "",merged_obj@meta.data$merged_group)
merged_obj@meta.data$merged_group <- gsub("^NA_", "", merged_obj@meta.data$merged_group)
merged_obj@meta.data$merged_group <- gsub("_NA$", "", merged_obj@meta.data$merged_group)
merged_obj@meta.data$merged_group <- gsub("_NA_NA_[0-9]+$", "",merged_obj@meta.data$merged_group)
merged_obj@meta.data$merged_group <- gsub("^NA_", "", merged_obj@meta.data$merged_group)
merged_obj@meta.data$merged_group <- gsub("_NA$", "", merged_obj@meta.data$merged_group)
merged_obj@meta.data$merged_group <- gsub("_NA_NA_[0-9]+$", "",merged_obj@meta.data$merged_group)

unique(merged_obj@meta.data$merged_group)

unique(merged_obj$MPR_Response2)
merged_obj$Patients
aver_dt<-AverageExpression(merged_obj,
                           group.by = c("merged_group"),
                           assays = "RNA")

#aver_dt <- as.data.frame(aver_dt$RNA)
aver_dt=aver_dt[[1]]
head(aver_dt)
cg=names(tail(sort(apply(aver_dt, 1, sd)),1000))
aver_dt_matrix <- as.matrix(aver_dt)

# corr_imm <- cor(aver_dt_matrix[cg,],method = 'spearman')
# pmat_imm <- cor.mtest(aver_dt_matrix[cg,],method = 'spearman')


library(psych)
cor_result <- corr.test(aver_dt_matrix[cg,], method = "pearson")
# 提取相关系数矩阵
correlation_matrix <- cor_result$r
# 提取p值矩阵
p_value_matrix <- cor_result$p

my_col <- colorRampPalette(rev(brewer.pal(7, "RdBu")))
#corrlpot绘制相关性热图
corrplot(correlation_matrix, p.mat = p_value_matrix, method = 'color', diag = T, type = 'full',
         sig.level = c(0.001, 0.01, 0.05), pch.cex = 0.7,
         insig = 'label_sig',  order = 'hclust', tl.col="black", tl.cex = 0.7, col=my_col(200))

pheatmap::pheatmap(cor(aver_dt_matrix[cg,],method = 'pearson'))

corrplot(correlation_matrix, p.mat = p_value_matrix, sig.level = 0.05, insig = "blank", type = "full", tl.col = "black")

library(ggcorrplot)
library(corrplot)
test_cor <- cor(aver_dt_matrix[cg,], method = "spearman") #pearson
test_p<- round(cor_pmat(aver_dt_matrix[cg,],method = "spearman"),3)
corrplot(corr =test_cor,type="upper",tl.pos="lt",tl.col="black",p.mat = test_p, insig = "label_sig", 
         sig.level = c(.01, .05),pch.cex=2,pch.col = "black",order = "original")#order = "AOE"
corrplot(corr = test_cor,type="lower",add=TRUE,method="number",
         tl.pos="n",tl.col="black",
         col="black",tl.cex=1.2,diag=FALSE, cl.pos="n",pch.col = "black",
         number.cex = 0.8,order = "original")

if (F) {
  library(ComplexHeatmap)
  library(circlize)
  library(ggcorrplot)
  
  test_cor <- cor(aver_dt_matrix[cg,], method = "spearman") #pearson
  test_p<- round(cor_pmat(aver_dt_matrix[cg,],method = "spearman"),3)
  diag(test_p) <- 1
  Heatmap(
    test_cor, 
    name = " ", 
    border = "white",
    col = colorRamp2(c(0,0.5,1), c("#5296cc", "yellow","#f064af")),
    #colorRampPalette(c("#417e46", "#f7f8f9", "#6f2f7e"))(100),
    rect_gp = gpar(col = "white", lwd = 2),
    row_order = rownames(test_cor),
    row_names_side = "left", 
    row_names_gp = gpar(fontsize = 10),
    column_order = colnames(test_cor),
    column_names_side = "top",  
    column_names_gp = gpar(
      fontsize = 10,fontface="italic",
      col=c(rep("#ff339f",14),rep("#5facee",20))),
    
    cell_fun = function(j, i, x, y, width, height, fill) {
      if (test_p[i,j] < 0.001) {
        grid.text("***",  x, y, gp = gpar(col = 'white',fontsize = 10))
      } else if (test_p[i,j]<0.01) {
        grid.text("**", x, y, gp = gpar(fontsize = 10))
      } else if (test_p[i,j]<0.05) {
        grid.text("*",  x, y, gp = gpar(fontsize = 10))
      }else {
        grid.text("", x, y)
      }
      
    }
  )
  
}
####不同成纤维细胞差异基因表达---
library(dplyr)
library(Seurat)
library(ggplot2)
library(patchwork)
library(ggsci)
library(circlize)
library(RColorBrewer)
library(ComplexHeatmap)

CAFs_cells_object$CAFs_cells_type
Idents(CAFs_cells_object)<-"CAFs_cells_type"
pbmc.markers <- FindAllMarkers(CAFs_cells_object,
                               only.pos = TRUE,
                               min.pct = 0.25,
                               logfc.threshold = 0.25)

# get top 10 genes
top5pbmc.markers <- pbmc.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)

col <- pal_npg()(9)
data$rid <- factor(data$rid, levels=c("CXCL2+iCAFs","CXCR4+iCAFs","TPSAB1+iCAFs","MyCAFs","apCAFs","eCAFs"))

cell_type_color<- c(
  "CXCL2+iCAFs" ="#00bfff",
  "CXCR4+iCAFs" ="#ff1493",
  "TPSAB1+iCAFs"="#66CDAA",
  "MyCAFs"="#FFa07a",
  "apCAFs"="#0000FF",
  "eCAFs"="#FF00ff")

# plot
top5pbmc.markers$gene[7]<-"CXCL2"

markers <- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/成纤维细胞/marker画图.csv",header = T)

DoHeatmap(CAFs_cells_object,features = markers$gene_name,
          group.colors = cell_type_color) +
  scale_colour_npg() +
  scale_fill_gradient2(low = '#0099CC',mid = 'white',high = '#CC0033',
                       name = 'Z-score')




####不同成纤维细胞在不同分组中的差异基因通路分析----
library(gggsea)
library(ggplot2)
library(GSVA)
library(enrichplot)
library(GSEABase)
library(fgsea)
library(clusterProfiler)
unique(CAFs_cells_object$MPR_Response2)
Idents(CAFs_cells_object)<-"MPR_Response2"
deg_all=FindMarkers(CAFs_cells_object, ident.1 = "Pre_MPR", ident.2 = "Pre_NMPR")
write.csv(deg_all, "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/成纤维细胞/gsea通路分析/Pre_MPR_Pre_NMPR_deg_all.csv", row.names = T)
deg_all2=FindMarkers(CAFs_cells_object, ident.1 = "Post_MPR", ident.2 = "Pre_MPR")
write.csv(deg_all2, "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/成纤维细胞/gsea通路分析/Post_MPR_Pre_MPR_deg_all.csv", row.names = T)
deg_all3=FindMarkers(CAFs_cells_object, ident.1 = "Post_MPR", ident.2 = "Pre_NMPR")
write.csv(deg_all3, "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/成纤维细胞/gsea通路分析/Post_MPR_Pre_NMPR_deg_all.csv", row.names = T)


####pre_nmpr_pre_mpr_deg_gsea----

deg_all<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/成纤维细胞/gsea通路分析/Pre_MPR_Pre_NMPR_deg_all.csv")
colnames(deg_all)[1]<- "gene_name"
deg<- deg_all
head(deg)

gmtfile <- "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/h.all.v2024.1.Hs.symbols.gmt"
#txtfile<- "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/髓系细胞亚型/M1M2特征得分和GESA分析/特征得分marker.txt"

hallmark <- read.gmt(gmtfile)
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
            nrow = 2,
            ncol = 5
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

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/成纤维细胞/gsea通路分析/Pre_MPR_Pre_NMPR_gsea.pdf",p_nes,width=13,height=6)

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

####成纤维细胞和免疫细胞的通讯分析all----

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


library(CellChat)
library(tidyverse)
library(Seurat)

CAFs_cells_object$CAFs_cells_type

CAFs_cells_object<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/CAFs_cells_object.rds")
CAFs_cells_object = UpdateSeuratObject(CAFs_cells_object)
CAFs_cells_object


myeloid_cells_object=readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/myeloid_cells_object_filter使用.rds")
myeloid_cells_object = UpdateSeuratObject(myeloid_cells_object)
myeloid_cells_object

B_cells_object<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/B_cells_object_combined使用.rds")
B_cells_object = UpdateSeuratObject(B_cells_object)
B_cells_object

T_cells_object<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/T_cells_object.rds")
T_cells_object = UpdateSeuratObject(T_cells_object)
T_cells_object

merged_obj <- merge(CAFs_cells_object, y = c(T_cells_object,B_cells_object,myeloid_cells_object), add.cell.ids = c("CAFs_cells_object","T_cells_object","B_cells_object","myeloid_cells_object"))

merged_obj <- merge(CAFs_cells_object, y = c(T_cells_object), add.cell.ids = c("CAFs_cells_object","T_cells_object"))


View(merged_obj)
merged_obj$CAFs_cells_type
merged_obj$cd4_cd8_group
merged_obj$B_cells_type
merged_obj$myeloid_cells_type

merged_obj@meta.data$merged_group <- paste(merged_obj$CAFs_cells_type,
                                           merged_obj$cd4_cd8_group,
                                           merged_obj$B_cells_type,
                                           merged_obj$myeloid_cells_type,
                                           sep = "_")


merged_obj@meta.data$merged_group <- paste(merged_obj$CAFs_cells_type,
                                           merged_obj$cd4_cd8_group,
                                           sep = "_")

# 查看更新后的 meta.data
head(merged_obj@meta.data)
merged_obj@meta.data$merged_group <- gsub("^NA_", "", merged_obj@meta.data$merged_group)
merged_obj@meta.data$merged_group <- gsub("_NA$", "", merged_obj@meta.data$merged_group)
merged_obj@meta.data$merged_group <- gsub("_NA_NA_[0-9]+$", "",merged_obj@meta.data$merged_group)
merged_obj@meta.data$merged_group <- gsub("^NA_", "", merged_obj@meta.data$merged_group)
merged_obj@meta.data$merged_group <- gsub("_NA$", "", merged_obj@meta.data$merged_group)
merged_obj@meta.data$merged_group <- gsub("_NA_NA_[0-9]+$", "",merged_obj@meta.data$merged_group)
merged_obj@meta.data$merged_group <- gsub("^NA_", "", merged_obj@meta.data$merged_group)
merged_obj@meta.data$merged_group <- gsub("_NA$", "", merged_obj@meta.data$merged_group)
merged_obj@meta.data$merged_group <- gsub("_NA_NA_[0-9]+$", "",merged_obj@meta.data$merged_group)



merged_obj@meta.data$merged_group <- gsub("^NA_", "", merged_obj@meta.data$merged_group)
merged_obj@meta.data$merged_group <- gsub("_NA$", "", merged_obj@meta.data$merged_group)
unique(merged_obj@meta.data$merged_group)


merged_obj_filter<- subset(merged_obj,subset=merged_group %in% c("MyCAFs","apCAFs","CXCL2+iCAFs","CXCR4+iCAFs","eCAFs","TPSAB1+iCAFs",
                                                                 "CD8_C1-GNLY","CD8_C4-GZMK","CD8_C3-IL7R"))


pre_nmpr_CAFs_chat<- subset(merged_obj_filter, subset=MPR_Response2=="Pre_NMPR")
unique(pre_nmpr_CAFs_chat$MPR_Response2)
pre_mpr_CAFs_chat<- subset(merged_obj_filter, subset=MPR_Response2=="Pre_MPR")
unique(pre_mpr_CAFs_chat$MPR_Response2)
post_mpr_CAFs_chat<- subset(merged_obj_filter, subset=MPR_Response2=="Post_MPR")
unique(post_mpr_CAFs_chat$MPR_Response2)
options(future.globals.maxSize = 3e+09)
pre_nmpr.cellchat <- KS_cellchat(pre_nmpr_CAFs_chat, assay = 'RNA', group.by = "merged_group",
                                 workers=2, species='human')

saveRDS(pre_nmpr.cellchat, file = "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/成纤维细胞/细胞通讯分析/免疫细胞/pre_nmpr.CAFs.T.cellchat.rds")

options(future.globals.maxSize = 3e+09) 
pre_mpr.cellchat <- KS_cellchat(pre_mpr_CAFs_chat, assay = 'RNA', group.by = "merged_group",
                                workers=2, species='human')

saveRDS(pre_mpr.cellchat , file = "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/成纤维细胞/细胞通讯分析/免疫细胞/pre_mpr.CAFs.T.cellchat.rds")

options(future.globals.maxSize = 3e+09) 
post_mpr.cellchat <- KS_cellchat(post_mpr_CAFs_chat, assay = 'RNA', group.by = "merged_group",
                                 workers=2, species='human')

saveRDS(post_mpr.cellchat, file = "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/成纤维细胞/细胞通讯分析/免疫细胞/post_mpr.CAFs.T.cellchat.rds")


pre_nmpr.cellchat<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/成纤维细胞/细胞通讯分析/免疫细胞/pre_nmprCAFs.T.cellchat.rds")

pre_mpr.cellchat<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/成纤维细胞/细胞通讯分析/免疫细胞/pre_mpr.CAFs.T.cellchat.rds")

post_mpr.cellchat<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/成纤维细胞/细胞通讯分析/免疫细胞/post_mpr.CAFs.T.cellchat.rds")


####有无细菌----

CAFs.T.bacteria_yes_chat<- subset(merged_obj_filter, subset=group_microbe2=="Bacteria+")
CAFs.T.bacteria_yes_pre_nmpr<- subset(CAFs.T.bacteria_yes_chat, subset=MPR_Response2=="Pre_NMPR")
unique(CAFs.T.bacteria_yes_pre_nmpr$MPR_Response2)
CAFs.T.bacteria_yes_pre_mpr<- subset(CAFs.T.bacteria_yes_chat, subset=MPR_Response2=="Pre_MPR")
unique(CAFs.T.bacteria_yes_pre_mpr$MPR_Response2)
CAFs.T.bacteria_yes_post_mpr<- subset(CAFs.T.bacteria_yes_chat, subset=MPR_Response2=="Post_MPR")
unique(CAFs.T.bacteria_yes_post_mpr$MPR_Response2)

CAFs.T.bacteria_yes_chat<- KS_cellchat(CAFs.T.bacteria_yes_chat, assay = 'RNA', group.by = "merged_group",
                                                workers=2, species='human')

saveRDS(CAFs.T.bacteria_yes_chat,
        file = "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/髓系细胞亚型/细胞通讯分析/CAFs.T.bacteria_yes_chat.rds")

CAFs.T.bacteria_yes_pre_nmpr_chat<- KS_cellchat(CAFs.T.bacteria_yes_pre_nmpr, assay = 'RNA', group.by = "merged_group",
                                                             workers=2, species='human')

CAFs.T.bacteria_yes_pre_mpr_chat<- KS_cellchat(CAFs.T.bacteria_yes_pre_mpr, assay = 'RNA', group.by = "merged_group",
                                                            workers=2, species='human')

CAFs.T.bacteria_yes_post_mpr_chat<- KS_cellchat(CAFs.T.bacteria_yes_post_mpr, assay = 'RNA', group.by = "merged_group",
                                                             workers=2, species='human')


saveRDS(CAFs.T.bacteria_yes_pre_nmpr_chat,
        file = "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/髓系细胞亚型/细胞通讯分析/CAFs.T.bacteria_yes_pre_nmpr_chat.rds")

saveRDS(CAFs.T.bacteria_yes_pre_mpr_chat,
        file = "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/髓系细胞亚型/细胞通讯分析/CAFs.T.bacteria_yes_pre_mpr_chat.rds")

saveRDS(CAFs.T.bacteria_yes_post_mpr_chat,
        file = "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/髓系细胞亚型/细胞通讯分析/CAFs.T.bacteria_yes_post_mpr_chat.rds")



CAFs.T.bacteria_no_chat<- subset(merged_obj_filter, subset=group_microbe2=="Bacteria-")
CAFs.T.bacteria_no_pre_nmpr<- subset(CAFs.T.bacteria_no_chat, subset=MPR_Response2=="Pre_NMPR")
unique(CAFs.T.bacteria_no_pre_nmpr$MPR_Response2)
CAFs.T.bacteria_no_pre_mpr<- subset(CAFs.T.bacteria_no_chat, subset=MPR_Response2=="Pre_MPR")
unique(CAFs.T.bacteria_no_pre_mpr$MPR_Response2)
CAFs.T.bacteria_no_post_mpr<- subset(CAFs.T.bacteria_no_chat, subset=MPR_Response2=="Post_MPR")
unique(CAFs.T.bacteria_no_post_mpr$MPR_Response2)

options(future.globals.maxSize = 3e+09) 
CAFs.T.bacteria_no_chat<- KS_cellchat(CAFs.T.bacteria_no_chat, assay = 'RNA', group.by = "merged_group",
                                       workers=2, species='human')

saveRDS(CAFs.T.bacteria_no_chat,
        file = "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/髓系细胞亚型/细胞通讯分析/CAFs.T.bacteria_no_chat.rds")


CAFs.T.bacteria_no_pre_nmpr_chat<- KS_cellchat(CAFs.T.bacteria_no_pre_nmpr, assay = 'RNA', group.by = "merged_group",
                                                workers=2, species='human')

CAFs.T.bacteria_no_pre_mpr_chat<- KS_cellchat(CAFs.T.bacteria_no_pre_mpr, assay = 'RNA', group.by = "merged_group",
                                               workers=2, species='human')

CAFs.T.bacteria_no_post_mpr_chat<- KS_cellchat(CAFs.T.bacteria_no_post_mpr, assay = 'RNA', group.by = "merged_group",
                                                workers=2, species='human')


saveRDS(CAFs.T.bacteria_no_pre_nmpr_chat,
        file = "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/髓系细胞亚型/细胞通讯分析/CAFs.T.bacteria_no_pre_nmpr_chat.rds")

saveRDS(CAFs.T.bacteria_no_pre_mpr_chat,
        file = "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/髓系细胞亚型/细胞通讯分析/CAFs.T.bacteria_no_pre_mpr_chat.rds")

saveRDS(CAFs.T.bacteria_no_post_mpr_chat,
        file = "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/髓系细胞亚型/细胞通讯分析/CAFs.T.bacteria_no_post_mpr_chat.rds")



#############################################################

pheatmap::pheatmap(pre_nmpr.cellchat@net$count, border_color = "black", 
                   cluster_cols = F, fontsize = 10, cluster_rows = F,
                   display_numbers = T,number_color="black",number_format = "%.0f")


object.list <- list(bacteria_yes= CAFs.T.bacteria_yes_chat, bacteria_no=CAFs.T.bacteria_no_chat)

object.list <- list(Pre_NMPR= pre_nmpr.cellchat, Pre_MPR= pre_mpr.cellchat,Post_MPR=post_mpr.cellchat)

object.list <- list(Pre_NMPR= CAFs.T.bacteria_yes_pre_nmpr_chat, 
                    Pre_MPR= CAFs.T.bacteria_yes_pre_mpr_chat,
                    Post_MPR=CAFs.T.bacteria_yes_post_mpr_chat)

object.list <- list(Pre_NMPR= CAFs.T.bacteria_no_pre_nmpr_chat, 
                    Pre_MPR= CAFs.T.bacteria_no_pre_mpr_chat,
                    Post_MPR=CAFs.T.bacteria_no_post_mpr_chat)

cellchat <- mergeCellChat(object.list, add.names = names(object.list))



gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2,3))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2,3), measure = "weight")
p1<-gg1 + gg2

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/成纤维细胞/细胞通讯分析/免疫细胞/分组比较柱状图_IL7R_GNLY_GZMK有细菌存在分组.pdf",p1,width=6,height=4)

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
ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/成纤维细胞/细胞通讯分析/成纤维和T细胞通讯/combined_plots有细菌存在三个分组.pdf", p_combined, width = 10, height = 4) 


gg4 <- rankNet(cellchat, mode = "comparison",comparison = c(1, 2), stacked = T,font.size = 10, sources.use = "CD8_C1-GNLY",targets.use = c("CXCL2+iCAFs","CXCR4+iCAFs","TPSAB1+iCAFs",
                                                                                                                                           "MyCAFs","apCAFs","eCAFs","CD8_C2-CTSW","CD8_C1-GNLY"),do.stat = TRUE, color.use = c("#4974a4","#4dae47"))

gg5 <- rankNet(cellchat, mode = "comparison",comparison = c(1, 3), stacked = T, font.size = 10,sources.use = "CD8_C1-GNLY",targets.use = c("CXCL2+iCAFs","CXCR4+iCAFs","TPSAB1+iCAFs",
                                                                                                                                           "MyCAFs","apCAFs","eCAFs","CD8_C2-CTSW","CD8_C1-GNLY"),do.stat = TRUE,color.use = c("#4974a4","#f29600"))

gg6 <- rankNet(cellchat, mode = "comparison",comparison = c(2, 3), stacked = T, font.size = 10,sources.use = "CD8_C1-GNLY",targets.use = c("CXCL2+iCAFs","CXCR4+iCAFs","TPSAB1+iCAFs",
                                                                                                                                           "MyCAFs","apCAFs","eCAFs","CD8_C2-CTSW","CD8_C1-GNLY"),do.stat = TRUE,color.use = c("#4dae47","#f29600"))

gg4+gg5+gg6



pbubble<-netVisual_bubble(cellchat, sources.use = "CD8_C1-GNLY", targets.use =c("CXCL2+iCAFs","CXCR4+iCAFs","TPSAB1+iCAFs",
                                                                                "MyCAFs","apCAFs","eCAFs"),comparison = c(1,2,3), angle.x = 45,font.size = 12,color.text =c("#4974a4","#4dae47","#f29600"))



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




pathways.show <- c("MHC-I")
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show)
par(mfrow = c(1,3), xpd=TRUE)

pathways.show <- c("MHC-I")
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}


pheatmap::pheatmap(pre_nmpr.cellchat@net$count, border_color = "black", 
                   cluster_cols = F, fontsize = 10, cluster_rows = F,
                   display_numbers = T,number_color="black",number_format = "%.0f")

pathways.show <- c("MHC-I") 
netVisual_aggregate(pre_nmpr.cellchat, signaling = pathways.show)
pairLR <- extractEnrichedLR(pre_nmpr.cellchat, signaling = pathways.show, geneLR.return = FALSE)
netVisual_individual(pre_nmpr.cellchat, signaling = pathways.show, pairLR.use = "HLA-DPA1_CD4", layout = "circle")

gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2


gg1 <- netVisual_heatmap(cellchat)
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
pdf(file = "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/成纤维细胞/细胞通讯分析/成纤维和T细胞通讯/netVisual_heatmap.pdf", width = 10, height = 6)
gg1 + gg2
dev.off()

#ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/成纤维细胞/细胞通讯分析/成纤维和T细胞通讯/netVisual_heatmap.pdf",p,width=7,height=5)









#####不同细胞群之间相互作用数量或相互作用强度的差异.红色表示与第一个数据集相比，互作增强，蓝色对应的是减弱

object.list <- list(Pre_NMPR= pre_nmpr.cellchat, Pre_MPR= pre_mpr.cellchat)
object.list <- list(Pre_NMPR= pre_nmpr.cellchat, Pre_MPR= pre_mpr.cellchat,Post_MPR=post_mpr.cellchat)

cellchat <- mergeCellChat(object.list, add.names = names(object.list))

pdf(file = "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/成纤维细胞/细胞通讯分析/成纤维和T细胞通讯/pre_nmpr_pre_mpr_netVisual_diffInteraction.pdf", width = 10, height = 6)
par(mfrow = c(1,2), xpd=TRUE)
p11<-netVisual_diffInteraction(cellchat, weight.scale = T)
p12<-netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")
dev.off()

object.list <- list(Pre_MPR= pre_mpr.cellchat,Post_MPR=post_mpr.cellchat)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

pdf(file = "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/成纤维细胞/细胞通讯分析/成纤维和T细胞通讯/pre_mpr_post_mpr_netVisual_diffInteraction.pdf", width = 10, height = 6)
par(mfrow = c(1,2), xpd=TRUE)
p11<-netVisual_diffInteraction(cellchat, weight.scale = T)
p12<-netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")
dev.off()


cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")
cellchat <- netClustering(cellchat, type = "functional")

library(ComplexHeatmap)
i = 1
# combining all the identified signaling pathways from different datasets 
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)

ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 6, height = 12)

ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 6, height = 12)

ht3 = netAnalysis_signalingRole_heatmap(object.list[[i+2]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 6, height = 12)
draw(ht1 + ht2+ht3, ht_gap = unit(0.5, "cm"))

netVisual_bubble(cellchat, sources.use = "CD8_C1-GNLY",comparison = c(1,2,3), angle.x = 45,font.size = 12,color.text =c("#4974a4","#4dae47","#f29600"))



####和髓系通讯----

library(CellChat)
library(tidyverse)
library(Seurat)

CAFs_cells_object$CAFs_cells_type

CAFs_cells_object<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/CAFs_cells_object.rds")
CAFs_cells_object = UpdateSeuratObject(CAFs_cells_object)
CAFs_cells_object

myeloid_cells_object=readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/myeloid_cells_object_filter使用.rds")
myeloid_cells_object = UpdateSeuratObject(myeloid_cells_object)
myeloid_cells_object

merged_obj <- merge(CAFs_cells_object, y = myeloid_cells_object, add.cell.ids = c("CAFs_cells_object","myeloid_cells_object"))

View(merged_obj)
merged_obj$CAFs_cells_type
merged_obj$myeloid_cells_type


merged_obj@meta.data$merged_group <- paste(merged_obj$CAFs_cells_type,
                                           merged_obj$myeloid_cells_type,
                                           sep = "_")

# 查看更新后的 meta.data
head(merged_obj@meta.data)
unique(merged_obj@meta.data$merged_group)

merged_obj@meta.data$merged_group <- gsub("^NA_", "", merged_obj@meta.data$merged_group)
merged_obj@meta.data$merged_group <- gsub("_NA$", "", merged_obj@meta.data$merged_group)
merged_obj@meta.data$merged_group <- gsub("_NA_NA_[0-9]+$", "",merged_obj@meta.data$merged_group)
merged_obj@meta.data$merged_group <- gsub("_[0-9]+$", "",merged_obj@meta.data$merged_group)

unique(merged_obj@meta.data$merged_group)

merged_obj_filter<- subset(merged_obj,subset=merged_group %in% c("CXCL2+iCAFs","CXCR4+iCAFs","TPSAB1+iCAFs","MyCAFs","apCAFs","eCAFs",
                                                                 "Mast cells", "SPP1+TAMs", "FOLR2+TAMs", 
                                                                 "ZBTB16+TAMs", "CXCL10+TAMs", "Neutrophil", 
                                                                 "Monocyte-derived dendritic cells", "cDC1", 
                                                                 "cDC2", "Migratory cDCs"))

unique(merged_obj_filter@meta.data$merged_group)


pre_nmpr_CAFs_chat<- subset(merged_obj_filter, subset=MPR_Response2=="Pre_NMPR")
unique(pre_nmpr_CAFs_chat$MPR_Response2)
pre_mpr_CAFs_chat<- subset(merged_obj_filter, subset=MPR_Response2=="Pre_MPR")
unique(pre_mpr_CAFs_chat$MPR_Response2)
post_mpr_CAFs_chat<- subset(merged_obj_filter, subset=MPR_Response2=="Post_MPR")
unique(post_mpr_CAFs_chat$MPR_Response2)
options(future.globals.maxSize = 3e+09)
pre_nmpr.cellchat <- KS_cellchat(pre_nmpr_CAFs_chat, assay = 'RNA', group.by = "merged_group",
                                 workers=2, species='human')


saveRDS(pre_nmpr.cellchat, file = "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/成纤维细胞/细胞通讯分析/成纤维细胞和髓系细胞通讯/pre_nmpr.CAFs.cellchat.all.rds")

options(future.globals.maxSize = 3e+09) 
pre_mpr.cellchat <- KS_cellchat(pre_mpr_CAFs_chat, assay = 'RNA', group.by = "merged_group",
                                workers=2, species='human')

saveRDS(pre_mpr.cellchat , file = "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/成纤维细胞/细胞通讯分析/成纤维细胞和髓系细胞通讯/pre_mpr.CAFs.cellchat.all.rds")

options(future.globals.maxSize = 3e+09) 
post_mpr.cellchat <- KS_cellchat(post_mpr_CAFs_chat, assay = 'RNA', group.by = "merged_group",
                                 workers=2, species='human')

saveRDS(post_mpr.cellchat, file = "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/成纤维细胞/细胞通讯分析/成纤维细胞和髓系细胞通讯/post_mpr.CAFs.cellchat.all.rds")


pre_nmpr.cellchat<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/成纤维细胞/细胞通讯分析/成纤维细胞和髓系细胞通讯/pre_nmpr.CAFs.cellchat.rds")

pre_mpr.cellchat<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/成纤维细胞/细胞通讯分析/成纤维细胞和髓系细胞通讯/pre_mpr.CAFs.cellchat.rds")

post_mpr.cellchat<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/成纤维细胞/细胞通讯分析/成纤维细胞和髓系细胞通讯/post_mpr.CAFs.cellchat.rds")


object.list <- list(Pre_NMPR= pre_nmpr.cellchat, Pre_MPR= pre_mpr.cellchat,Post_MPR=post_mpr.cellchat)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))



