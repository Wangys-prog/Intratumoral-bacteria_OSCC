library(Seurat)
library(tidyverse)
library(harmony)
library(dplyr)
library(patchwork)
library(SeuratWrappers)

MPRNMPR_object_miMPRobe<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/MPRNMPR_object_miMPRobe_remove.rds")
MPRNMPR_object_miMPRobe = UpdateSeuratObject(MPRNMPR_object_miMPRobe)
MPRNMPR_object_miMPRobe
View(MPRNMPR_object_miMPRobe)

MPRNMPR_object_miMPRobe_remove <- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/MPRNMPR_object_miMPRobe_remove.rds")
MPRNMPR_object_miMPRobe_remove = UpdateSeuratObject(MPRNMPR_object_miMPRobe_remove)
MPRNMPR_object_miMPRobe_remove

setwd("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/上皮细胞亚群注释")
setwd("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/上皮细胞亚群注释")

epi_cells_object= readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/epi_cells_object_nomerged.rds")
epi_cells_object = UpdateSeuratObject(epi_cells_object)

##上皮细胞使用不merged----

####提取myeloid cells subject----
Idents(MPRNMPR_object_miMPRobe)<-"cell_type_new"
levels(Idents(MPRNMPR_object_miMPRobe))
Idents(myeloid_cells_object) <-"RNA_snn_res.0.4"
levels(Idents(myeloid_cells_object))

#cluster8_12_13<- subset(myeloid_cells_object,idents = c("8","12","13"))

epi_cells_object=subset(MPRNMPR_object_miMPRobe,idents = c("Epithelial cells"))

#epi_cells_object = merge(epi_cells_object, y = cluster8_12_13, project = "CombinedProject")

table(Idents(epi_cells_object))

colnames(epi_cells_object@meta.data)

pbmc<- epi_cells_object

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
p1 <- DimPlot(pbmc, group.by = "orig.ident",reduction = "umap.harmony", pt.size=0.1,label = T) +
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

p6 <- DimPlot(pbmc, group.by = "orig.ident",reduction = "tsne", pt.size=0.1,label = T)
p6

p7 <- DimPlot(pbmc, group.by = "orig.ident",reduction = "umap", pt.size=0.1,label = T)

p7

p8 <- DimPlot(pbmc, group.by = "RNA_snn_res.0.4",reduction = "umap", pt.size=0.1,label = T)

p8

DefaultAssay(pbmc) <- "RNA"

epi_cells_object <- pbmc
Idents(epi_cells_object) <-"RNA_snn_res.0.4"
table(Idents(epi_cells_object))
DimPlot(epi_cells_object, group.by = "RNA_snn_res.0.4",reduction = "umap.harmony", pt.size=0.1,label = T) +
  ggtitle("RNA_snn_res.0.4")

DimPlot(epi_cells_object_merged, group.by = "RNA_snn_res.0.4",reduction = "umap.harmony", pt.size=0.1,label = T) +
  ggtitle("RNA_snn_res.0.4")

markers <- FindAllMarkers(epi_cells_object, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
all.markers =markers %>% dplyr::select(gene, everything()) %>% subset(p_val<0.05)
top20 =all.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
p_heatmap<-DoHeatmap(epi_cells_object,features = top20$gene)+ #热图复现
  theme(text = element_text(size = 5))

write.csv(all.markers,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/上皮细胞亚群注释/epi_cells_marker_all.markers.csv")
write.csv(top20,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/上皮细胞亚群注释/epi_cells_marker_top20.csv")

#保存RDS
saveRDS(epi_cells_object,"D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/epi_cells_object_nomerged.rds")

####使用不merged的细胞继续分析----


p1 <- DimPlot(epi_cells_object, group.by = "orig.ident",reduction = "umap.harmony", pt.size=0.1) +
  ggtitle("Integrated by harmony")#去批次后的可视化

p1
p2 <- DimPlot(epi_cells_object, group.by = "RNA_snn_res.0.4",reduction = "umap.harmony", pt.size=0.1,label = T) +
  ggtitle("RNA_snn_res.0.4")#去批次后的可视化

p2

p3 <- DimPlot(epi_cells_object, group.by = "RNA_snn_res.0.4",reduction = "tsne.harmony", pt.size=0.1,label = T) +
  ggtitle("RNA_snn_res.0.4")#去批次后的可视化
p3

p4 <- DimPlot(epi_cells_object, group.by = "RNA_snn_res.0.4",reduction = "tsne.harmony", pt.size=0.1,label = T) +
  ggtitle("Integrated by harmony")#去批次后的可视化
p4

####
####marker基因小提琴图----

setwd("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/上皮细胞亚群注释")

cell_marker<- read.csv("marker1.csv")

######删除找不到的基因

nofinder<- c("DNAIC2", "OXJ1", "CYP2F2", "TRP63", "POPN", "RTNK2", "AGE")

cell_marker_filtered <- cell_marker[!cell_marker$gene_name %in% nofinder, ]
 
Idents(epi_cells_object)<-"RNA_snn_res.0.4"
levels(Idents(epi_cells_object))
unique(cell_marker_filtered$gene_name)
library(randomcoloR)

my36colors = distinctColorPalette(84)

vln.dat=FetchData(epi_cells_object,c(unique(cell_marker_filtered$gene_name),"RNA_snn_res.0.4"))
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
                                                                              "13","14","15","16","17","18","19","20","21","22",
                                                                              "23","24","25","26","27","28"))

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

ggsave("上皮细胞marker小提琴图RNA_snn_res.0.4.jpg",p1,width = 7,height = 7)

####合并公司infercnv结果数据----
#sys.setenv(JAGS_HOME= "C:/Users/Wangy/AppData/Local/Programs/JAGS/JAGS-4.3.1")

library("rjags")
library(infercnv)
library(Seurat)
library(tidyverse)
library(harmony)
library(dplyr)
library(patchwork)
dir.create("~/SeuratV4")
install.packages('Seurat', repos = c('https://satijalab.r-universe.dev'), lib = "~/SeuratV4")
.libPaths(c("~/SeuratV4", .libPaths()))
detach("package:infercnv", unload = TRUE)
detach("package:Seurat", unload = TRUE)
unloadNamespace("Seurat")
.rs.restartR()
.libPaths(c("~/SeuratV4", .libPaths()))
library(Seurat)
packageVersion("Seurat")

####infercnv feature plot 绘制----
####读取全部样本的rds----
MPRNMPR_object_miMPRobe<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/MPRNMPR_object_miMPRobe_remove.rds")
MPRNMPR_object_miMPRobe = UpdateSeuratObject(MPRNMPR_object_miMPRobe)
MPRNMPR_object_miMPRobe

cnv_level_all <-read.table("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/口腔鳞癌CNV分析数据/公司infercnv结果/DZOE2024041345_cnv中间文件/vis_all/cnv_result_all.xls",sep ="\t",header = T)
head(cnv_level_all)
unique(cnv_level_all$clusters)
MPRNMPR_object_miMPRobe@meta.data$cnv_level <-cnv_level_all[match(MPRNMPR_object_miMPRobe@meta.data$cell_id,cnv_level_all$barcode_inuse),6]
MPRNMPR_object_miMPRobe@meta.data$cnv_group <-cnv_level_all[match(MPRNMPR_object_miMPRobe@meta.data$cell_id,cnv_level_all$barcode_inuse),7]
table(MPRNMPR_object_miMPRobe@meta.data$cnv_group)
saveRDS(MPRNMPR_object_miMPRobe, file = "D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/MPRNMPR_object_miMPRobe_remove.rds")

p1 = FeaturePlot(MPRNMPR_object_miMPRobe, reduction = 'umap.harmony',
                 features = "cnv_level",raster=FALSE) +
  xlab("UMAP_1") + ylab("UMAP_2")
p1
p2 = VlnPlot(seurat_ob , features = "cnv_level" , group.by = "orig.ident",pt.size = 0)
p2



MPRNMPR.umap.harmonyredu<- MPRNMPR_object_miMPRobe@reductions$umap.harmony@cell.embeddings%>%
  as.data.frame() %>% 
  cbind(cell_type = MPRNMPR_object_miMPRobe@meta.data$cell_type_new)%>% 
  cbind(cnv_level = MPRNMPR_object_miMPRobe@meta.data$cnv_level) %>%
  cbind(cnv_group = MPRNMPR_object_miMPRobe@meta.data$cnv_group)

names(MPRNMPR.umap.harmonyredu)

MPRNMPR.umap.harmonyredu_clean <- na.omit(MPRNMPR.umap.harmonyredu)
head(MPRNMPR.umap.harmonyredu_clean )


MPRNMPR_cell_type <- ggplot(MPRNMPR.umap.harmonyredu_clean ,aes(x= umapharmony_1, y = umapharmony_2,color = cnv_level)) +  
  geom_point(size = 0.5, alpha =0.5)  +  
  scale_color_gradientn(colors = c("green4", "yellow", "red"),
                        values = scales::rescale(c(1, 30, 100), to = c(0, 1))) +
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
  theme(legend.title = element_text(size=13),
        legend.key=element_rect(fill='white'), #
        legend.text = element_text(size=12), #设置legend标签的大小
        legend.key.size=unit(0.8,'cm') ) +  # 设置legend标签之间的大小
  geom_segment(aes(x = min(MPRNMPR.umap.harmonyredu_clean$umapharmony_1) , y = min(MPRNMPR.umap.harmonyredu_clean$umapharmony_2) ,
                   xend = min(MPRNMPR.umap.harmonyredu_clean$umapharmony_1) +3, yend = min(MPRNMPR.umap.harmonyredu_clean$umapharmony_2) ),
               colour = "black", size=1,arrow = arrow(length = unit(0.5,"cm")))+ 
  geom_segment(aes(x = min(MPRNMPR.umap.harmonyredu_clean$umapharmony_1)  , y = min(MPRNMPR.umap.harmonyredu_clean$umapharmony_2)  ,
                   xend = min(MPRNMPR.umap.harmonyredu_clean$umapharmony_1) , yend = min(MPRNMPR.umap.harmonyredu_clean$umapharmony_2) + 3),
               colour = "black", size=1,arrow = arrow(length = unit(0.5,"cm"))) +
  annotate("text", x = min(MPRNMPR.umap.harmonyredu_clean$umapharmony_1) +1.5, y = min(MPRNMPR.umap.harmonyredu_clean$umapharmony_2) -1, label = "UMAP_1",
           color="black",size = 5, fontface="bold" ) + 
  annotate("text", x = min(MPRNMPR.umap.harmonyredu_clean$umapharmony_1) -1, y = min(MPRNMPR.umap.harmonyredu_clean$umapharmony_2) + 1.5, label = "UMAP_2",
           color="black",size = 5, fontface="bold" ,angle=90)

MPRNMPR_cell_type

ggsave('D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/上皮细胞亚群注释/infercnv_plot/infercnv_feature.pdf',MPRNMPR_cell_type,width =7,height = 6)

####小提琴图----

sum(is.na(MPRNMPR_object_miMPRobe$cnv_level))

# 过滤掉包含缺失值的细胞
MPRNMPR_object_miMPRobe_clean <- subset(MPRNMPR_object_miMPRobe, cells = Cells(MPRNMPR_object_miMPRobe)[!is.na(MPRNMPR_object_miMPRobe$cnv_level)])
MPRNMPR_object_miMPRobe_clean$cell_type_new <- factor(MPRNMPR_object_miMPRobe_clean$cell_type_new, 
                                               levels = c('T cells', 'B cells','Plasma cells', 'Myeloid cells',  
                                                          'Mast cells',"Endothelial cells","Fibroblasts","Smooth muscle cells",
                                                          "Pericytes","Epithelial cells"))

cell_type_color <- c("#FF6347","#ffc556","#1E90FF","#0000FF","#9370DB","#00FFFF","#32cd32","#006400","#FF00ff","#20B2AA")
p2 = VlnPlot(MPRNMPR_object_miMPRobe_clean, features = "cnv_level", group.by = "cell_type_new",pt.size = 0,raster = FALSE)+
  scale_fill_manual(values = cell_type_color)+xlab(NULL) +  # 去除 X 轴标题
  geom_hline(yintercept = 25, linetype = "dashed", color = "black",size=0.8) + 
  ylab("CNV level") +  # 设置 Y 轴标题为 "CNV level"
  ggtitle(NULL) +  # 去除整体标题
  theme(
    axis.title.x = element_blank(),  # 去除 X 轴标题
    plot.title = element_blank()  # 去除整体标题
  )


df <- MPRNMPR_object_miMPRobe_clean@meta.data %>%
  dplyr::select(cell_type_new, cnv_level) %>%
  filter(!is.na(cnv_level))

df$cell_type_new <- factor(df$cell_type_new, levels = c('T cells', 'B cells','Plasma cells', 'Myeloid cells',  
                                                                 'Mast cells',"Endothelial cells","Fibroblasts","Smooth muscle cells",
                                                                 "Pericytes","Epithelial cells"))
p2 <- ggplot(df, aes(x = cell_type_new, y = cnv_level, fill = cell_type_new)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.15, outlier.shape = NA) +
  stat_summary(fun = mean, geom = "point", size = 1, color = "red") +
  geom_hline(yintercept = 25, linetype = "dashed", color = "black", size = 0.8) +
  scale_fill_manual(values = cell_type_color) +
  ylab("CNV level") +
  theme_classic() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"   # ← 去掉图例
  )

p2

cnv_celltype_data <- FetchData(
  MPRNMPR_object_miMPRobe_remove, 
  vars = c("cnv_level", "cell_type_new")
) %>%
  filter(!is.na(cell_type_new), !is.na(cnv_level))  # 去除NA值

# 手动指定细胞类型顺序（按需修改，避免默认字母序混乱）
# 先查看所有细胞类型名称：unique(cnv_celltype_data$cell_type_new)
cell_type_order <- unique(cnv_celltype_data$cell_type_new)  # 保留原始顺序，可手动调整
cnv_celltype_data$cell_type_new <- factor(cnv_celltype_data$cell_type_new, levels = cell_type_order)

# 批量计算所有细胞类型的两两Wilcoxon检验（非参数检验，适配单细胞数据）
df_p_val_celltype <- cnv_celltype_data %>%
  pairwise_wilcox_test(
    cnv_level ~ cell_type_new,  # 因变量：cnv_level；自变量：细胞类型
    p.adjust.method = "fdr",    # 多重检验校正（控制假发现率，适合多组对比）
    paired = FALSE
  ) %>%
  # 自动计算标注位置（step.increase控制多组对比横杠的垂直间距，避免重叠）
  add_xy_position(
    x = "cell_type_new", 
    dodge = 0.6, 
    step.increase = 0.15  # 细胞类型越多，可适当增大（如0.2）
  )


ggsave('D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/上皮细胞亚群注释/infercnv_plot/infercnv_vlnplot.pdf',p2,width =8,height = 6)

ggsave("G:/OSCC_manuscripts/oscc_process/MPR_NMPR_group/细胞亚群注释/上皮细胞亚群注释/infercnv_plot/infercnv_vlnplot.pdf",p2,width =8,height = 6)


####根据公司结果提取恶性和非恶性细胞分组----

epi_infercnv<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/口腔鳞癌CNV分析数据/公司infercnv结果/DZOE2024041345_cnv中间文件/vis_epithelial/epi_cnv_seurat.rds")
epi_infercnv = UpdateSeuratObject(epi_infercnv)

View(epi_infercnv)
epi_infercnv@meta.data$clusters
unique(epi_infercnv@meta.data$cell_type_new)

table(epi_infercnv@meta.data$epi_infercnv_group)

write.csv(epi_infercnv@meta.data,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/上皮细胞亚群注释/infercnv_plot/epi_infercnv_group.csv")

epi.umap.harmonyredu<- epi_infercnv@reductions$umap.harmony@cell.embeddings%>%
  as.data.frame() %>% 
  cbind(cell_type = epi_infercnv@meta.data$cell_type_new)%>% 
  cbind(cnv_level = epi_infercnv@meta.data$cnv_level) %>%
  cbind(cnv_group = epi_infercnv@meta.data$cnv_group) %>%
  cbind(epi_infercnv_group= epi_infercnv@meta.data$epi_infercnv_group)

names(epi.umap.harmonyredu)
head(epi.umap.harmonyredu)

p1 = FeaturePlot(epi_infercnv, features = "cnv_level" )
p1
p2 = VlnPlot(epi_infercnv, features = "cnv_level" , group.by = "clusters",pt.size = 0)
p2
####使用自己的epi数据RNA_snn_res.0.4 作图vlnplot----

epi_cells_object= readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/epi_cells_object_nomerged.rds")
epi_cnv_level <-read.table("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/口腔鳞癌CNV分析数据/公司infercnv结果/DZOE2024041345_cnv中间文件/vis_epithelial/epi_cnv_result.xls",sep="\t",header = T)
head(epi_cnv_level)
epi_cells_object@meta.data$cnv_level <-cnv_level_all[match(epi_cells_object@meta.data$cell_id,epi_cnv_level$barcode_inuse),6]
epi_cells_object@meta.data$cnv_group <-cnv_level_all[match(epi_cells_object@meta.data$cell_id,epi_cnv_level$barcode_inuse),7]

unique(epi_cells_object@meta.data$cnv_group)
unique(epi_cells_object@meta.data$RNA_snn_res.0.3)

epi_cells_object@meta.data$epi_infercnv_group <- epi_infercnv@meta.data[match(epi_cells_object@meta.data$cell_id,epi_infercnv@meta.data$cell_id),"epi_infercnv_group"]


library(randomcoloR)

colors_epi<- distinctColorPalette(29)

p3 = VlnPlot(epi_cells_object, features = "cnv_level" , group.by = "RNA_snn_res.0.4",pt.size = 0)+
  scale_fill_manual(values = colors_epi)+xlab(NULL) +  # 去除 X 轴标题
  ylab("CNV level") +  # 设置 Y 轴标题为 "CNV level"
  geom_hline(yintercept = 25, linetype = "dashed", color = "black",size=0.8) + 
  ggtitle(NULL) +  # 去除整体标题
  theme(
    axis.title.x = element_blank(),  # 去除 X 轴标题
    plot.title = element_blank()  # 去除整体标题
  )
p3

df <- epi_cells_object@meta.data %>%
  dplyr::select(RNA_snn_res.0.4, cnv_level) %>%
  filter(!is.na(cnv_level))

unique(df$RNA_snn_res.0.4)

p3 <- ggplot(df, aes(x = RNA_snn_res.0.4, y = cnv_level, fill = RNA_snn_res.0.4)) +
  geom_violin(trim = FALSE) +
  #geom_boxplot(width = 0.15, outlier.shape = NA) +
  stat_summary(fun = mean, geom = "point", size = 1, color = "red") +
  geom_hline(yintercept = 25, linetype = "dashed", color = "black", size = 0.8) +
  scale_fill_manual(values = colors_epi) +
  ylab("CNV level") +
  theme_classic() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 16),  # 调大 X 轴字体
    axis.text.y = element_text(size = 14),                         # 调大 Y 轴字体
    axis.title.y = element_text(size = 14),                        # 调大 Y 轴标题
    legend.position = "none"
  )

p3

ggsave('D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/上皮细胞亚群注释/infercnv_plot/epi_infercnv_vlnplot.pdf',p3,width =9,height = 6)

ggsave("G:/OSCC_manuscripts/oscc_process/MPR_NMPR_group/细胞亚群注释/上皮细胞亚群注释/infercnv_plot/epi_infercnv_vlnplot_cluster.pdf",p3,width =9,height = 6)


colorlist <- c("#a4cde1","#67a4cc","#277fb8","#549da3","#96cb8f","#8bc96d","#4dae47","#5c9e43","#b79973","#f38989","#ec5051","#e32427","#ef6a45","#f9b769","#f9a341","#f48521","#ee8e46","#d4a6a8","#af93c4","#8660a8","#815e99","#c6b598","#f6f28f","#d4a55b","#b05a28")
colorlist <- c("#4974a4","#4dae47","#f29600")
epi_cells_object@meta.data$MPR_Response2 <- factor(epi_cells_object@meta.data$MPR_Response2,levels = c("Pre_NMPR","Pre_MPR","Post_MPR"))

p4 = VlnPlot(epi_cells_object, features = "cnv_level" , group.by = "MPR_Response2",pt.size = 0)+
  scale_fill_manual(values = colorlist)+xlab(NULL) +  # 去除 X 轴标题
  ylab("CNV level") +  # 设置 Y 轴标题为 "CNV level"
  geom_hline(yintercept = 25, linetype = "dashed", color = "black",size=0.8) + 
  ggtitle(NULL) +  # 去除整体标题
  theme(
    axis.title.x = element_blank(),  # 去除 X 轴标题
    plot.title = element_blank()  # 去除整体标题
  )
p4


df <- epi_cells_object@meta.data %>%
  dplyr::select(MPR_Response2, cnv_level) %>%
  filter(!is.na(cnv_level))

unique(df$MPR_Response2)

df$MPR_Response2 <- factor(df$MPR_Response2, levels = c('Pre_NMPR', 'Pre_MPR','Post_MPR'))
p4 <- ggplot(df, aes(x = MPR_Response2, y = cnv_level, fill = MPR_Response2)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.15, outlier.shape = NA) +
  stat_summary(fun = mean, geom = "point", size = 1, color = "red") +
  geom_hline(yintercept = 25, linetype = "dashed", color = "black", size = 0.8) +
  scale_fill_manual(values = colorlist) +
  ylab("CNV level") +
  theme_classic() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 16),  # 调大 X 轴字体
    axis.text.y = element_text(size = 14),                         # 调大 Y 轴字体
    axis.title.y = element_text(size = 14),                        # 调大 Y 轴标题
    legend.position = "none"
  )

p4

ggsave('G:/OSCC_manuscripts/oscc_process/MPR_NMPR_group/细胞亚群注释/上皮细胞亚群注释/infercnv_plot/epi_infercnv_vlnplot_分组平均值.pdf',p4,width =7,height = 6)


library(ggplot2)
library(dplyr)
library(rstatix)  # 用于pairwise_wilcox_test批量计算两两p值
library(ggpubr)   # 用于stat_pvalue_manual绘制显著性标注

# ---------------------- 1. 数据准备+批量计算两两对比p值 ----------------------
cnv_data <- FetchData(epi_cells_object, vars = c("cnv_level", "MPR_Response2")) %>%
  filter(!is.na(MPR_Response2), !is.na(cnv_level)) %>%
  # 重命名列名（与你示例一致，避免函数报错，可选）
  rename(MPR_response2 = MPR_Response2, cnv_percent = cnv_level)

# 手动指定分组顺序（按需修改，确保x轴顺序符合预期）
group_order <- c("Pre_NMPR", "Pre_MPR", "Post_MPR")  # 替换为你的实际分组
cnv_data$MPR_response2 <- factor(cnv_data$MPR_response2, levels = group_order)

# 批量计算所有两两分组的Wilcoxon检验p值（Bonferroni校正）
df_p_val <- cnv_data %>%
  pairwise_wilcox_test(
    cnv_percent ~ MPR_response2,  # 因变量~自变量
    p.adjust.method = "bonferroni",  # 多重检验校正（与你示例一致）
    paired = FALSE
  ) %>%
  # 自动计算标注的x轴位置（dodge=0.6控制横杠宽度，与小提琴图适配）
  add_xy_position(x = "MPR_response2", dodge = 0.6, step.increase = 0.1)

# 查看计算的p值和标注位置（可选，确认对比组对正确）
print(df_p_val[, c("group1", "group2", "p.adj", "p.adj.signif", "y.position")])

# ---------------------- 2. 最终美化版：两两分组显著性标注（复刻你的示例风格） ----------------------
colorlist <- c("#4974a4","#4dae47","#f29600")  # 与你示例一致的颜色（可替换为你的颜色）

p4 <- ggplot(cnv_data, aes(x = MPR_response2, y = cnv_percent, fill = MPR_response2)) +
  # 1. 小提琴图：加粗边框+优化填充（复刻boxplot的加粗风格）
  geom_violin(
    scale = "width", 
    alpha = 0.9,
    linewidth = 1.2,  # 加粗边框（与你示例的lwd=1一致，稍作优化）
    color = "black"
  ) +
  # 2. 颜色映射（与你示例的colorlist一致）
  scale_fill_manual(values = colorlist) +
  # 3. 轴标签与标题（放大+加粗，适配学术图）
  xlab(NULL) +
  ylab("CNV level") +
  ggtitle(NULL) +
  # 4. 阈值线：加粗虚线（保留之前的需求）
  geom_hline(
    yintercept = 25, 
    linetype = "dashed", 
    color = "black", 
    linewidth = 1.0
  ) +
  # 5. 核心：两两分组显著性标注（完全复刻你的示例逻辑）
  stat_pvalue_manual(
    df_p_val,  # 之前计算的两两对比p值数据框
    label = "p.adj.signif",  # 显示校正后的星号
    hide.ns = FALSE,  # 显示不显著（ns），若想隐藏设为TRUE
    tip.length = 0,  # 横杠两端无箭头（与你示例一致）
    label.size = 7,  # 星号大小（放大，适配整体字体）
    color = "black",  # 星号颜色（与你示例一致）
    fontface = "bold",  # 星号加粗
    linewidth = 1.2  # 对比横杠加粗（与小提琴边框呼应）
  ) +
  # 6. 主题：无背景+加粗边框+放大字体（复刻theme_cowplot的简洁风格）
  theme_cowplot() +
  theme(
    # 坐标轴文字（放大+加粗，比你示例更大，更醒目）
    axis.text.x = element_text(
      size = 18, 
      face = "bold", 
      color = "black",
      margin = margin(t = 5)
    ),
    axis.text.y = element_text(
      size = 16, 
      face = "bold", 
      color = "black",
      margin = margin(r = 5)
    ),
    # 轴标题（放大+加粗）
    axis.title.y = element_text(
      size = 18, 
      face = "bold", 
      color = "black",
      margin = margin(r = 10)
    ),
    # 坐标轴轴线+刻度（加粗，与你示例一致）
    axis.line.x = element_line(color = "black", linewidth = 1.2),
    axis.line.y = element_line(color = "black", linewidth = 1.2),
    axis.ticks.x = element_line(color = "black", linewidth = 2),
    axis.ticks.y = element_line(color = "black", linewidth = 2),
    # 隐藏图例（x轴已标注分组，与你示例一致）
    legend.position = "none",
    # 图边距（避免标注超出画布）
    plot.margin = margin(10, 15, 10, 15)
  )

# 显示最终图
p4
ggsave('D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/上皮细胞亚群注释/infercnv_plot/epi_infercnv_vlnplot_分组显著性.pdf',p4,width =9,height = 6)


############################################


epi_cells_object@meta.data$epi_infercnv_group <- ifelse(epi_cells_object@meta.data$cnv_level<25,
                                                    "Non-malignant epithelial cells", 
                                                    "Malignant epithelial cells")
table(epi_cells_object@meta.data$epi_infercnv_group)

epi.umap.harmonyredu<- epi_cells_object@reductions$umap.harmony@cell.embeddings%>%
  as.data.frame() %>% 
  cbind(cell_type = epi_cells_object@meta.data$cell_type_new)%>% 
  cbind(cnv_level = epi_cells_object@meta.data$cnv_level) %>%
  cbind(cnv_group = epi_cells_object@meta.data$cnv_group) %>%
  cbind(epi_infercnv_group= epi_cells_object@meta.data$epi_infercnv_group)

head(epi.umap.harmonyredu)

epi.umap.harmonyredu_clean <- na.omit(epi.umap.harmonyredu)
head(epi.umap.harmonyredu_clean)

epi_cell_type <- ggplot(epi.umap.harmonyredu_clean ,aes(x= umapharmony_1, y = umapharmony_2,color = cnv_level)) +  
  geom_point(size = 0.5, alpha =0.5)  +  
  # scale_color_gradientn(colors = c("green4", "yellow", "red"),
  #                       values = scales::rescale(c(1,30, 50,100), to = c(0, 1))) +
  scale_color_gradientn(
    colors = c("darkgreen", "yellow", "orange", "red"),  # 更鲜艳的颜色过渡
    values = scales::rescale(c(1, 25, 50, 100), to = c(0, 1))) +
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
  theme(legend.title = element_text(size=13),
        legend.key=element_rect(fill='white'), #
        legend.text = element_text(size=12), #设置legend标签的大小
        legend.key.size=unit(0.8,'cm') ) +  # 设置legend标签之间的大小
  geom_segment(aes(x = min(epi.umap.harmonyredu_clean$umapharmony_1) , y = min(epi.umap.harmonyredu_clean$umapharmony_2) ,
                   xend = min(epi.umap.harmonyredu_clean$umapharmony_1) +3, yend = min(epi.umap.harmonyredu_clean$umapharmony_2) ),
               colour = "black", size=1,arrow = arrow(length = unit(0.5,"cm")))+ 
  geom_segment(aes(x = min(epi.umap.harmonyredu_clean$umapharmony_1)  , y = min(epi.umap.harmonyredu_clean$umapharmony_2)  ,
                   xend = min(epi.umap.harmonyredu_clean$umapharmony_1) , yend = min(epi.umap.harmonyredu_clean$umapharmony_2) + 3),
               colour = "black", size=1,arrow = arrow(length = unit(0.5,"cm"))) +
  annotate("text", x = min(epi.umap.harmonyredu_clean$umapharmony_1) +1.5, y = min(epi.umap.harmonyredu_clean$umapharmony_2) -1, label = "UMAP_1",
           color="black",size = 5, fontface="bold" ) + 
  annotate("text", x = min(epi.umap.harmonyredu_clean$umapharmony_1) -1, y = min(epi.umap.harmonyredu_clean$umapharmony_2) + 1.5, label = "UMAP_2",
           color="black",size = 5, fontface="bold" ,angle=90)

epi_cell_type

ggsave('D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/上皮细胞亚群注释/infercnv_plot/epi_infercnv_feature.pdf',epi_cell_type,width =7,height = 6)


colorlist <- c("#a4cde1","#67a4cc","#277fb8","#549da3","#96cb8f","#8bc96d","#4dae47","#5c9e43","#b79973","#f38989","#ec5051","#e32427","#ef6a45","#f9b769","#f9a341","#f48521","#ee8e46","#d4a6a8","#af93c4","#8660a8","#815e99","#c6b598","#f6f28f","#d4a55b","#b05a28")
colorlist <- c("#f29600","#4974a4","#b81316")

epi_infercnv_group_color <- c("#f29600","#a4cde1")
epi_cell_type2 <- ggplot(epi.umap.harmonyredu,aes(x= umapharmony_1, y = umapharmony_2,color = epi_infercnv_group)) +  
  geom_point(size = 0.2, alpha =0.5)  +  
  scale_color_manual(values =epi_infercnv_group_color)+
  #scale_color_brewer(palette = "Set1")+
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
  theme(legend.title =element_blank(),
        legend.key=element_rect(fill='white'), #
        legend.text = element_text(size=12), #设置legend标签的大小
        legend.key.size=unit(0.8,'cm') ) +  # 设置legend标签之间的大小
  guides(color = guide_legend(override.aes = list(size=5)))+ #设置legend中 点的大小
  geom_segment(aes(x = min(epi.umap.harmonyredu_clean$umapharmony_1) , y = min(epi.umap.harmonyredu_clean$umapharmony_2) ,
                   xend = min(epi.umap.harmonyredu_clean$umapharmony_1) +3, yend = min(epi.umap.harmonyredu_clean$umapharmony_2) ),
               colour = "black", size=0.7,arrow = arrow(length = unit(0.5,"cm")))+ 
  geom_segment(aes(x = min(epi.umap.harmonyredu_clean$umapharmony_1)  , y = min(epi.umap.harmonyredu_clean$umapharmony_2)  ,
                   xend = min(epi.umap.harmonyredu_clean$umapharmony_1) , yend = min(epi.umap.harmonyredu_clean$umapharmony_2) + 3),
               colour = "black", size=0.7,arrow = arrow(length = unit(0.5,"cm"))) +
  annotate("text", x = min(epi.umap.harmonyredu_clean$umapharmony_1) +1.5, y = min(epi.umap.harmonyredu_clean$umapharmony_2) -1, label = "UMAP_1",
           color="black",size = 4, fontface="bold" ) + 
  annotate("text", x = min(epi.umap.harmonyredu_clean$umapharmony_1) -1, y = min(epi.umap.harmonyredu_clean$umapharmony_2) + 1.5, label = "UMAP_2",
           color="black",size = 4, fontface="bold" ,angle=90)
  
epi_cell_type2

ggsave('D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/上皮细胞亚群注释/infercnv_plot/epi_infercnv_malignant.pdf',epi_cell_type2,width =8.5,height = 6)


####恶性非恶性细胞在不同患者中分布情况----
####恶性非恶性箱装图----
library(tidyverse)
library(gghalves)
library(rstatix)
library(ggpubr)
library(ggsci)

df <- read_tsv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/上皮细胞亚群注释/infercnv_plot/ggplot2优雅绘制半边箱线图20240427/data.xls") %>% 
  mutate(year=as.character(year))
head(df)
######差异分析----
df_p_val1 <- df %>% group_by(continent) %>%
  wilcox_test(lifeExp  ~ year) %>%
  adjust_pvalue(p.col = "p", method = "bonferroni") %>%
  add_significance(p.col = "p.adj") %>% 
  add_xy_position(x = "year", dodge = 0.8) 

df_p_val_right <-p2_filter %>% group_by(group)%>%
  pairwise_wilcox_test(
    value ~ type, 
    p.adjust.method = "bonferroni"  # 调整方法可以选择 "bonferroni", "holm", "fdr" 等
  ) %>% 
  add_xy_position(x="group",dodge=0.6)


###数据准

colorlist <- c("#4974a4","#4dae47","#f29600")

epi_malignant_data<- epi_cells_object@meta.data$orig.ident%>%
  as.data.frame() %>% 
  cbind(MPR_response2 = epi_cells_object@meta.data$MPR_Response2) %>%
  cbind(epi_infercnv_group= epi_cells_object@meta.data$epi_infercnv_group)

names(epi_malignant_data) <-c("patient_id","MPR_response2","epi_infercnv_group")

cell_count <- epi_malignant_data %>%
  group_by(patient_id,MPR_response2, epi_infercnv_group) %>%
  summarise(cell_n = n()) %>%
  ungroup()

write.csv(cell_count,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/上皮细胞亚群注释/infercnv_plot/恶性细胞治疗后下降数据.csv")

library(rstatix)
library(ggpubr)
names(cell_count)
cell_count_filter <- cell_count %>% filter(epi_infercnv_group =="Malignant epithelial cells")
cell_count_filter$MPR_response2<- factor(cell_count_filter$MPR_response2,levels=c("Pre_NMPR","Pre_MPR","Post_MPR"))
df_p_val1 <-cell_count_filter%>%
  pairwise_wilcox_test(
    cell_n ~ MPR_response2, 
    p.adjust.method = "bonferroni"  # 调整方法可以选择 "bonferroni", "holm", "fdr" 等
  ) %>% 
  add_xy_position(x="MPR_response2",dodge=0.6)

colorlist <- c("#f29600","#4dae47","#4974a4")

colorlist <- c("#4974a4","#4dae47","#f29600")


library(cowplot)  

pp1 = ggplot(cell_count_filter,aes(x=MPR_response2,y=cell_n)) + #ggplot作图
  geom_boxplot(aes(fill=MPR_response2),outlier.shape = NA,lwd= 1)+
  scale_fill_manual(values = colorlist)+
  stat_pvalue_manual(df_p_val1,label="p.adj.signif",hide.ns = F,
                     tip.length = 0,label.size = 6,color="black")+
  theme_cowplot() +
  theme(axis.text = element_text(size = 12),
        axis.line.x=element_line(linetype=1,color="black",size=1),
        axis.ticks.x=element_line(color="black",size=2,lineend = 1),
        axis.line.y=element_line(linetype=1,color="black",size=1),
        axis.ticks.y=element_line(color="black",size=2,lineend = 1),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        axis.text.x = element_text(vjust=0,hjust=0.5,angle=0,size=16),
        axis.text.y = element_text(vjust=1,hjust=1,angle=0,size=15),
        legend.title = element_text(size = 12),
        plot.title = element_text(size = 18),
        legend.position = 'none') + 
  labs(title = "Malignant epithelial cells",y='',x="") +
  guides(fill = guide_legend(title = NULL))

pp1


pp2=ggplot(cell_count_filter, aes(x = MPR_response2, y = cell_n)) + 
  geom_boxplot(aes(fill = MPR_response2), outlier.shape = NA, lwd = 1) +
  geom_jitter(width = 0.2, size = 3, alpha = 0.8) +  # ← 添加散点
  scale_fill_manual(values = colorlist) +
  stat_pvalue_manual(df_p_val1, label = "p.adj.signif", hide.ns = FALSE,
                     tip.length = 0, label.size = 6, color = "black") +
  cowplot::theme_cowplot() +  # 推荐加命名空间，避免未加载时报错
  theme(
    axis.text = element_text(size = 12),
    axis.line.x = element_line(linetype = 1, color = "black", size = 1),
    axis.ticks.x = element_line(color = "black", size = 2, lineend = 1),
    axis.line.y = element_line(linetype = 1, color = "black", size = 1),
    axis.ticks.y = element_line(color = "black", size = 2, lineend = 1),
    axis.title = element_text(size = 12),
    legend.text = element_text(size = 12),
    axis.text.x = element_text(vjust = 0, hjust = 0.5, angle = 0, size = 16),
    axis.text.y = element_text(vjust = 1, hjust = 1, angle = 0, size = 15),
    legend.title = element_text(size = 12),
    plot.title = element_text(size = 18),
    legend.position = 'none'
  ) + 
  labs(title = "Malignant epithelial cells", y = '', x = "") +
  guides(fill = guide_legend(title = NULL))

ggsave('D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/上皮细胞亚群注释/infercnv_plot/epi_malignant分组箱状图加散点.pdf',pp2,width =6,height = 6)

saveRDS(epi_cells_object, file = "D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/epi_cells_object_nomerged.rds")



library(dplyr)
library(ggplot2)
library(rstatix)
library(ggpubr)
library(cowplot)

#-----------------------------------------------------------
# ① 准备数据
#-----------------------------------------------------------
epi_malignant_data <- epi_cells_object@meta.data$orig.ident %>%
  as.data.frame() %>% 
  cbind(MPR_response2 = epi_cells_object@meta.data$MPR_Response2) %>%
  cbind(epi_infercnv_group = epi_cells_object@meta.data$epi_infercnv_group)

names(epi_malignant_data) <- c("patient_id", "MPR_response2", "epi_infercnv_group")
View(epi_malignant_data)

write.csv(epi_malignant_data,"G:/OSCC_manuscripts/manuscripts/Cell reports/cell_reports_R2/投稿资料/3-绘图数据/Figure S2/epi_malignant_data.csv")
#-----------------------------------------------------------
# ② 每个 patient_id × MPR_response2 的总上皮细胞数
#-----------------------------------------------------------
total_cells <- epi_malignant_data %>%
  group_by(patient_id, MPR_response2) %>%
  summarise(total_n = n(), .groups = "drop")

#-----------------------------------------------------------
# ③ 每个 patient_id × MPR_response2 的恶性上皮细胞数
#-----------------------------------------------------------
malignant_cells <- epi_malignant_data %>%
  filter(epi_infercnv_group == "Malignant epithelial cells") %>%
  group_by(patient_id, MPR_response2) %>%
  summarise(malignant_n = n(), .groups = "drop")

#-----------------------------------------------------------
# ④ 合并并计算恶性细胞比例（百分比）
#-----------------------------------------------------------
tumor_pct <- total_cells %>%
  left_join(malignant_cells, by = c("patient_id", "MPR_response2")) %>%
  mutate(
    malignant_n = ifelse(is.na(malignant_n), 0, malignant_n),
    tumor_percent = malignant_n / total_n * 100,
    MPR_response2 = factor(MPR_response2, levels = c("Pre_NMPR", "Pre_MPR", "Post_MPR"))
  )

#-----------------------------------------------------------
# ⑤ 统计学检验（Wilcoxon）
#-----------------------------------------------------------
df_p_val1 <- tumor_pct %>%
  pairwise_wilcox_test(
    tumor_percent ~ MPR_response2,
    p.adjust.method = "bonferroni"
  ) %>%
  add_xy_position(x="MPR_response2", dodge = 0.6)

#-----------------------------------------------------------
# ⑥ 绘图（boxplot + scatter）
#-----------------------------------------------------------
colorlist <- c("#4974a4","#4dae47","#f29600")

pp <- ggplot(tumor_pct, aes(x=MPR_response2, y=tumor_percent)) +
  geom_boxplot(aes(fill=MPR_response2), outlier.shape = NA, lwd=1) +
  geom_jitter(width=0.2, size=3, alpha=0.8) +
  scale_fill_manual(values = colorlist) +
  stat_pvalue_manual(df_p_val1, label="p.adj.signif",
                     hide.ns = FALSE, tip.length = 0,
                     label.size = 6, color="black") +
  theme_cowplot() +
  theme(
    axis.text = element_text(size = 12),
    axis.line.x = element_line(color="black", size=1),
    axis.ticks.x = element_line(color="black", size=2),
    axis.line.y = element_line(color="black", size=1),
    axis.ticks.y = element_line(color="black", size=2),
    axis.title = element_text(size = 12),
    legend.position = 'none',
    axis.text.x = element_text(size=16),
    axis.text.y = element_text(size=15)
  ) +
  labs(title = "Malignant epithelial cells (%)", x="", y="")

pp

#-----------------------------------------------------------
# ⑦ 保存图形
#-----------------------------------------------------------
ggsave(
  "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/上皮细胞亚群注释/infercnv_plot/epi_malignant_percent_box_scatter.pdf",
  pp, width = 6, height = 6
)





####cnv结果写入到seurat项目中----
## 4. CNV LEVEL 写入seurat中
# infercnv_obj = readRDS(file.path(output_dir,"run.final.infercnv_obj"))
seurat_ob[['CNV']] = CreateAssayObject(data = infercnv_obj@expr.data)
# 计算每行的和，并去除缺失值
infercnv_level <- apply(as.data.frame(t(infercnv_obj@expr.data)), 1, function(x) {
  x[is.na(x)] <- 0
  return(sum(x))
})
# 求均值，并放大100倍，取整
infercnv_level <- round(scales::rescale(infercnv_level / nrow(infercnv_obj@expr.data), c(1, 100)), 0)
# 写入seurat
infercnv_level <- infercnv_level[Cells(seurat_ob)]
seurat_ob@meta.data$cnv_level = infercnv_level
## 5. 绘图
p1 = FeaturePlot(seurat_ob , features = "cnv_level")
p1
p2 = VlnPlot(seurat_ob , features = "cnv_level" , group.by = "orig.ident",pt.size = 0)
p2
ggsave(p1,file = file.path(output_dir,"featureplot.png"),height = 5,width = 5)
ggsave(p2,file = file.path(output_dir,"vlnplot.png"),height = 4,width = 7)


####infercnv 热图绘制----
library(infercnv)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library("RColorBrewer")

MPRNMPR_object_miMPRobe_remove <- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/MPRNMPR_object_miMPRobe_remove.rds")
MPRNMPR_object_miMPRobe_remove@meta.data$epi_infercnv_group

MPRNMPR_object_miMPRobe_remove@meta.data$epi_infercnv_group <-epi_cells_object@meta.data[match(MPRNMPR_object_miMPRobe_remove@meta.data$cell_id,epi_cells_object@meta.data$cell_id),"epi_infercnv_group"] 
saveRDS(MPRNMPR_object_miMPRobe_remove, file = "D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/MPRNMPR_object_miMPRobe_remove.rds")
infercnv_obj = readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/口腔鳞癌CNV分析数据/公司infercnv结果/DZOE2024041345_cnv中间文件/公司_run.final.infercnv_obj")
expr <- infercnv_obj@expr.data
head(expr)
dim(expr)
expr[1:4,1:4]
expr<-t(expr)
expr[1:4,1:4]

metadata<-MPRNMPR_object_miMPRobe_remove@meta.data$orig.ident %>% as.data.frame() %>%
  cbind(MPRNMPR_object_miMPRobe_remove@meta.data$MPR_Response) %>% cbind(MPRNMPR_object_miMPRobe_remove@meta.data$MPR_Response2) %>% 
  cbind(MPRNMPR_object_miMPRobe_remove@meta.data$cell_type_new) %>% cbind(MPRNMPR_object_miMPRobe_remove@meta.data$cnv_level) %>% 
  cbind(MPRNMPR_object_miMPRobe_remove@meta.data$cnv_group)%>% cbind(MPRNMPR_object_miMPRobe_remove@meta.data$cell_id) %>% 
  cbind(MPRNMPR_object_miMPRobe_remove@meta.data$epi_infercnv_group)

names(metadata)<- c("orig.ident","MPR_Response","MPR_Response2","cell_type_new","cnv_level","cnv_group","cell_id","epi_infercnv_group")

unique(metadata$cell_type_new)

ref_cell_id <- metadata %>% filter(cell_type_new %in% c("T cells", "B cells","Fibroblasts","Smooth muscle cells","Epithelial cells"))
#ref_cell_id <- metadata %>% filter(cell_type_new %in% c("Epithelial cells"))
names(ref_cell_id)
ref_cells_expr<- expr[match(ref_cell_id$cell_id,rownames(expr)),] %>% as.data.frame()
ref_cells_expr <- na.omit(ref_cells_expr)
ref_cells_expr[1:4,1:4]


ncol(ref_cells_expr)
nrow(ref_cells_expr)
ref_cells_expr <-t(ref_cells_expr)

gn <- rownames(ref_cells_expr)
geneFile <- read.table('D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/上皮细胞亚群注释/geneFile.txt')
length(geneFile$V1);length(gn)
sub_geneFile <-  geneFile[geneFile$V1%in%gn,]
dim(sub_geneFile)
gn2<-intersect(rownames(ref_cells_expr),sub_geneFile$V1)
ref_cells_expr <- ref_cells_expr[gn2,]
gn3 <- rownames(ref_cells_expr)
length(geneFile$V1);length(gn3)
sub_geneFile <-  geneFile[geneFile$V1%in%gn3,]
dim(sub_geneFile)

ref_cells_expr<- ref_cells_expr[match(sub_geneFile$V1,rownames(ref_cells_expr)), ]
identical(rownames(ref_cells_expr),sub_geneFile$V1)
ref_cells_expr[1:4,1:4]

# selected_columns <- sample(ncol(ref_cells_expr), 1000)
# # 提取随机选择的列
# ref_cells_expr_selected <- ref_cells_expr[, selected_columns]

matched_indices <- match(rownames(ref_cells_expr), sub_geneFile$V1)
# Extract corresponding values from V2
chr_values <- sub_geneFile$V2[matched_indices]
# Convert matrix to data frame and add the new column
df_matrix <- as.data.frame(ref_cells_expr)
df_matrix$chr <- chr_values
df_matrix<- df_matrix %>% filter(!chr %in% c("chrM","chrX","chrY"))
nrow(df_matrix)
ncol(df_matrix)

nrow(ref_cells_expr_selected)
ncol(ref_cells_expr_selected)

######

n <- t(df_matrix[,-ncol(df_matrix)])
ncol(n)
nrow(n)

df1<- metadata[,c("cell_id","cell_type_new")]
# df2<- metadata[,c("cell_id","epi_infercnv_group")]
# names(df2)
# names(df1)
# df_combined <- df1 %>%
#   left_join(df2, by = "cell_id") %>%
#   mutate(cell_type_new2 = ifelse(cell_type_new == "Epithelial cells", epi_infercnv_group, cell_type_new)) %>%
#   select(cell_id, cell_type_new2) # 移除不需要的列
# unique(df_combined$cell_type_new2)

matched_indices <- match(rownames(n), df1$cell_id)

# Extract corresponding values from V2
cell_type <- df1$cell_type_new[matched_indices]

# Convert matrix to data frame and add the new column
m_matrix <- as.data.frame(n)
m_matrix$group_cell_type <- cell_type

table(m_matrix$group_cell_type)

sample_sizes <- c("T cells" = 100, "B cells" = 100, "Fibroblasts"=100,"Smooth muscle cells"=100,"Epithelial cells" = 600)

m_matrix <- m_matrix %>%
  group_by(group_cell_type) %>%
  group_modify(~ slice_sample(.x, n = sample_sizes[.y$group_cell_type])) %>%
  ungroup()

unique(m_matrix$group_cell_type)
ncol(m_matrix)
nrow(m_matrix)

library(ComplexHeatmap)
# 调整配色
library("RColorBrewer")
# 来自函数：infercnv::plot_cnv
get_group_color_palette <- function () {
  return(colorRampPalette(RColorBrewer::brewer.pal(12, "Set3")))
}

color1 <- get_group_color_palette()(length(unique(df_matrix$chr)))
library(randomcoloR)
color1<- distinctColorPalette(22)
# 列标题的颜色框
unique(df_matrix$chr)
df_matrix$chr <- factor(df_matrix$chr, levels = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13",
                                                  "chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22" ))
top_color <- HeatmapAnnotation(
  cluster = anno_block(
    gp = gpar(fill = color1),  # Set fill colors for blocks
    labels = levels(df_matrix$chr),  # Set the labels
    labels_gp = gpar(cex = 0.9, col = "black")  # Set label size and color
  ))

unique_groups1 <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13",
                   "chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22" )


color2 <- c("#FF6347","#ffc556","#32cd32","#006400","#20B2AA")

unique(m_matrix$group_cell_type)

m_matrix$group_cell_type <- factor(m_matrix$group_cell_type, levels = c("T cells","B cells","Fibroblasts","Smooth muscle cells","Epithelial cells"))

unique_groups <- c("T cells","B cells","Fibroblasts","Smooth muscle cells","Epithelial cells")
names(color2) <- unique_groups
ha_row <- rowAnnotation(
  Group = m_matrix$group_cell_type,  # Add the group information
  col = list(Group = color2),  # Set the colors for the groups
  annotation_legend_param = list(
    Group = list(
      title = "Cell type",  # Title for the legend
      at = unique_groups,  # Unique group labels
      labels = unique_groups,  # Labels to show in the legend
      legend_gp = gpar(fill = color2[unique_groups])  # Colors corresponding to the groups
    )
  )
)

m_matrix$group_cell_type
m <- m_matrix %>%
  select(-group_cell_type)
ncol(m)
nrow(m)
max(m)
min(m)
m1<-m
m1[m1>1.2]=1.6
m1[m1<0.9]=0.8
m1[0.9<m1 & m1<1.2]=1

View(m)

write.csv(m,"G:/OSCC_manuscripts/manuscripts/Cell reports/cell_reports_R2/投稿资料/3-绘图数据/Figure S2/CNV_heatmap.csv")




library("circlize")
col_fun = colorRamp2(c(0.9, 1, 1.2), c("blue", "white", "red"))

pdf("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/上皮细胞亚群注释/infercnv_plot/infercnv_heatmap3.pdf",width = 11,height = 10)

ht_ref = Heatmap(as.matrix(m),
                   cluster_rows = F,
                   cluster_columns = F,
                   show_column_names = F,
                   show_row_names = F,
                   #col = col_fun,
                   column_split = df_matrix$chr,
                   row_split = m_matrix$group_cell_type,
                   heatmap_legend_param = list(
                     title = "Modified Expression",
                     title_position = "leftcenter-rot", # 图例标题位置
                     at=c(0.4,1.6), #图例范围
                     legend_height = unit(3, "cm") #图例长度
                   ),
                   top_annotation = top_color,
                   left_annotation = ha_row,
                   row_title = "Observations (Cells)",
                   row_title_side = c("right"),
                   column_title = "Genomic Region",
                   column_title_side = c("bottom"),
                  column_gap = unit(0.5, "mm"),
                 row_gap = unit(0.3, "mm")
                 )

draw(ht_ref, heatmap_legend_side = "right") # 图例位置

dev.off()

table(df_matrix$chr)
chr_info <- as.numeric(table(df_matrix$chr))
x=cumsum(chr_info)/chr_info[1]

pdf("heatmap_with_vline.pdf", width = 8, height = 6)

draw(ht_ref)
draw_x_line=function(x,ht1,color){
  decorate_heatmap_body(ht1, {
    grid.lines(c(x, x), c(0, 1), gp = gpar(lty = 2, lwd = 2, col=color))
  })
}

sapply(x, function(x1)draw_x_line(x1,ht1=ht_ref,color="black"))

dev.off()

####恶性和非恶性细胞feature dotplot----

epi_cells_object= readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/epi_cells_object_nomerged.rds")
epi_cells_object = UpdateSeuratObject(epi_cells_object)
epi_cells_object@meta.data$epi_infercnv_group

genes_of_interest <- c("EPCAM", "CDH1","SFN","KRT5","KRT8","KRT14","KRT16","KRT17", 
                       "KRT18","KRT6B")

dp<-DotPlot(epi_cells_object, features = genes_of_interest,group.by = "epi_infercnv_group") + 
  scale_color_gradient(low = "lightblue", high = "darkblue") +
  xlab(NULL)+
  ylab(NULL)+
  ggtitle("Typical epithelial marker genes") +
  theme(
    plot.title = element_text(size = 13, face = "bold", color = "black", hjust = 0.5)  # 自定义标题字体
  ) +
  RotatedAxis()

ggsave(dp,file ="D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/上皮细胞亚群注释/infercnv_plot/恶性和非恶性dotplot.pdf",height = 3.5,width = 9)



####CD45+冷热肿瘤对比----

####携带内源性细菌的 T 细胞或巨噬细胞与未携带细菌的细胞进行对比----

####绘制每位患者的CNV 小提琴图----

unique(epi_cells_object$MPR_Response2)

epi_cells_object$Patients

epi_cells_object$Response_Patient <- paste(
  epi_cells_object$MPR_Response2,
  epi_cells_object$Patients,
  sep = "_"
)

epi_cells_object$cnv_level

epi_cells_object$Response_Patient <- factor(
  epi_cells_object$Response_Patient,
  levels = c("Pre_NMPR_P3", "Pre_NMPR_P5", "Pre_NMPR_P6", "Pre_NMPR_P10",
             "Pre_MPR_P1","Pre_MPR_P2","Pre_MPR_P4","Pre_MPR_P7","Pre_MPR_P8","Pre_MPR_P9",
             "Post_MPR_P1","Post_MPR_P4","Post_MPR_P7","Post_MPR_P8")
)

VlnPlot(
  object = epi_cells_object,
  features = "cnv_level",
  group.by = "Response_Patient"
)

df <- epi_cells_object@meta.data %>%
  dplyr::select(Response_Patient, cnv_level) %>%
  filter(!is.na(cnv_level))

unique(df$Response_Patient)

df$Response_Patient <- factor(
  df$Response_Patient,
  levels = c("Pre_NMPR_P3", "Pre_NMPR_P5", "Pre_NMPR_P6", "Pre_NMPR_P10",
             "Pre_MPR_P1","Pre_MPR_P2","Pre_MPR_P4","Pre_MPR_P7","Pre_MPR_P8","Pre_MPR_P9",
             "Post_MPR_P1","Post_MPR_P4","Post_MPR_P7","Post_MPR_P8")
)

colorlist <- c(
  "#332288",  # 深紫
  "#88CCEE",  # 浅蓝
  "#44AA99",  # 蓝绿
  "#117733",  # 深绿
  "#999933",  # 暖黄
  "#DDCC77",  # 浅黄
  "#CC6677",  # 粉红
  "#882255",  # 紫红
  "#AA4499",  # 紫粉
  "#DDDDDD",  # 浅灰
  "#4477AA",  # 蓝
  "#66CCEE",  # 青
  "#228833",  # 绿
  "#EE6677"   # 珊瑚红
)

p5 <- ggplot(df, aes(x = Response_Patient, y = cnv_level, fill = Response_Patient)) +
  geom_violin(trim = FALSE) +
  #geom_boxplot(width = 0.15, outlier.shape = NA) +
  stat_summary(fun = mean, geom = "point", size = 1, color = "red") +
  geom_hline(yintercept = 25, linetype = "dashed", color = "black", size = 0.8) +
  scale_fill_manual(values = colorlist) +
  ylab("CNV level") +
  theme_classic() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 16),  # 调大 X 轴字体
    axis.text.y = element_text(size = 14),                         # 调大 Y 轴字体
    axis.title.y = element_text(size = 14),                        # 调大 Y 轴标题
    legend.position = "none"
  )

p5

ggsave("G:/OSCC_manuscripts/oscc_process/MPR_NMPR_group/细胞亚群注释/上皮细胞亚群注释/infercnv_plot/每位患者CNV平均值.pdf",p5,width=9,height=6)


