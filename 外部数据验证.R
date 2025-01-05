
library(Seurat)

# 创建 Seurat 对象
###P18 P23 P24 P27 P29 P32 P13 P14 P15 P16 P18
###P19 P20 P21 P22 P23 P24 P25 P26 P27 P28 P29 P30 P31 P32

sample_paths <- list(
  "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/多队列验证/GSE200996/GSE200996_RAW/GSM6048142_raw_feature_bc_matrix_P18_pre-Tx_GEX_sc_tumor.h5",
  "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/多队列验证/GSE200996/GSE200996_RAW/GSM6048143_raw_feature_bc_matrix_P23_pre-Tx_GEX_sc_tumor.h5",
  "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/多队列验证/GSE200996/GSE200996_RAW/GSM6048144_raw_feature_bc_matrix_P24_pre-Tx_GEX_sc_tumor.h5",
  "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/多队列验证/GSE200996/GSE200996_RAW/GSM6048145_raw_feature_bc_matrix_P27_pre-Tx_GEX_sc_tumor.h5",
  "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/多队列验证/GSE200996/GSE200996_RAW/GSM6048146_raw_feature_bc_matrix_P29_pre-Tx_GEX_sc_tumor.h5",
  "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/多队列验证/GSE200996/GSE200996_RAW/GSM6048147_raw_feature_bc_matrix_P32_pre-Tx_GEX_sc_tumor.h5",
  "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/多队列验证/GSE200996/GSE200996_RAW/GSM6048148_raw_feature_bc_matrix_P13_post-Tx_GEX_sc_tumor.h5",
  "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/多队列验证/GSE200996/GSE200996_RAW/GSM6048149_raw_feature_bc_matrix_P14_post-Tx_GEX_sc_tumor.h5",
  "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/多队列验证/GSE200996/GSE200996_RAW/GSM6048150_raw_feature_bc_matrix_P15_post-Tx_GEX_sc_tumor.h5",
  "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/多队列验证/GSE200996/GSE200996_RAW/GSM6048151_raw_feature_bc_matrix_P16_post-Tx_GEX_sc_tumor.h5",
  "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/多队列验证/GSE200996/GSE200996_RAW/GSM6048152_raw_feature_bc_matrix_P18_post-Tx_GEX_sc_tumor.h5",
  "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/多队列验证/GSE200996/GSE200996_RAW/GSM6048153_raw_feature_bc_matrix_P19_post-Tx_GEX_sc_tumor.h5",
  "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/多队列验证/GSE200996/GSE200996_RAW/GSM6048154_raw_feature_bc_matrix_P20_post-Tx_GEX_sc_tumor.h5",
  "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/多队列验证/GSE200996/GSE200996_RAW/GSM6048155_raw_feature_bc_matrix_P21_post-Tx_GEX_sc_tumor.h5",
  "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/多队列验证/GSE200996/GSE200996_RAW/GSM6048156_raw_feature_bc_matrix_P22_post-Tx_GEX_sc_tumor.h5",
  "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/多队列验证/GSE200996/GSE200996_RAW/GSM6048157_raw_feature_bc_matrix_P23_post-Tx_GEX_sc_tumor.h5",
  "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/多队列验证/GSE200996/GSE200996_RAW/GSM6048158_raw_feature_bc_matrix_P24_post-Tx_GEX_sc_tumor.h5",
  "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/多队列验证/GSE200996/GSE200996_RAW/GSM6048159_raw_feature_bc_matrix_P25_post-Tx_GEX_sc_tumor.h5",
  "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/多队列验证/GSE200996/GSE200996_RAW/GSM6048160_raw_feature_bc_matrix_P26_post-Tx_GEX_sc_tumor.h5",
  "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/多队列验证/GSE200996/GSE200996_RAW/GSM6048161_raw_feature_bc_matrix_P27_post-Tx_GEX_sc_tumor.h5",
  "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/多队列验证/GSE200996/GSE200996_RAW/GSM6048162_raw_feature_bc_matrix_P28_post-Tx_GEX_sc_tumor.h5",
  "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/多队列验证/GSE200996/GSE200996_RAW/GSM6048163_raw_feature_bc_matrix_P29_post-Tx_GEX_sc_tumor.h5",
  "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/多队列验证/GSE200996/GSE200996_RAW/GSM6048164_raw_feature_bc_matrix_P30_post-Tx_GEX_sc_tumor.h5",
  "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/多队列验证/GSE200996/GSE200996_RAW/GSM6048165_raw_feature_bc_matrix_P31_post-Tx_GEX_sc_tumor.h5",
  "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/多队列验证/GSE200996/GSE200996_RAW/GSM6048166_raw_feature_bc_matrix_P32_post-Tx_GEX_sc_tumor.h5"
)


# 定义样本名称，帮助后续标注数据
sample_names <- c("P18-pre","P23-pre","P24-pre","P27-pre","P29-pre","P32-pre","P13-post","P14-post","P15-post","P16-post","P18-post",
                  "P19-post","P20-post","P21-post","P22-post", "P23-post","P24-post","P25-post","P26-post",
                  "P27-post","P28-post","P29-post","P30-post","P31-post","P32-post")


sample_paths <- c(
  "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/多队列验证/GSE200996/GSE200996_RAW/GSM6048142_raw_feature_bc_matrix_P18_pre-Tx_GEX_sc_tumor.h5",
  "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/多队列验证/GSE200996/GSE200996_RAW/GSM6048143_raw_feature_bc_matrix_P23_pre-Tx_GEX_sc_tumor.h5",
  "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/多队列验证/GSE200996/GSE200996_RAW/GSM6048144_raw_feature_bc_matrix_P24_pre-Tx_GEX_sc_tumor.h5",
  "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/多队列验证/GSE200996/GSE200996_RAW/GSM6048145_raw_feature_bc_matrix_P27_pre-Tx_GEX_sc_tumor.h5",
  "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/多队列验证/GSE200996/GSE200996_RAW/GSM6048147_raw_feature_bc_matrix_P32_pre-Tx_GEX_sc_tumor.h5",
  "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/多队列验证/GSE200996/GSE200996_RAW/GSM6048148_raw_feature_bc_matrix_P13_post-Tx_GEX_sc_tumor.h5",
  "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/多队列验证/GSE200996/GSE200996_RAW/GSM6048152_raw_feature_bc_matrix_P18_post-Tx_GEX_sc_tumor.h5",
  "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/多队列验证/GSE200996/GSE200996_RAW/GSM6048155_raw_feature_bc_matrix_P21_post-Tx_GEX_sc_tumor.h5",
  "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/多队列验证/GSE200996/GSE200996_RAW/GSM6048157_raw_feature_bc_matrix_P23_post-Tx_GEX_sc_tumor.h5",
  "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/多队列验证/GSE200996/GSE200996_RAW/GSM6048158_raw_feature_bc_matrix_P24_post-Tx_GEX_sc_tumor.h5",
  "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/多队列验证/GSE200996/GSE200996_RAW/GSM6048161_raw_feature_bc_matrix_P27_post-Tx_GEX_sc_tumor.h5",
  "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/多队列验证/GSE200996/GSE200996_RAW/GSM6048166_raw_feature_bc_matrix_P32_post-Tx_GEX_sc_tumor.h5"
)

sample_names <- c("P18-pre","P23-pre","P24-pre","P27-pre","P32-pre","P13-post","P18-post",
                  "P21-post","P23-post","P24-post","P27-post","P32-post")




# # # 读取每个样本并存储到列表
# seurat_objects <- list()
# for(i in 1:length(sample_paths)){
#   print(i)
#   counts <- Read10X_h5(sample_paths[i])
#   seurat_objects[[i]] <- CreateSeuratObject(counts, min.cells=1) #和方法一一样，只是一个样本一个样本构建seurat对象
# } #把每个样本的数据创建一个seurat对象，并存放到列表scRNAlist里
# # 
# #seurat对象合并
# #scRNA2 <- merge(seurat_objects[[1]], y=seurat_objects[[2]]) #使用merge函数将2个seurat对象合并成一个seurat对象 #2个以上的样本合并seurat对象，需要写一个循环进行merge，因为merge合并只针对2个对象合并
# merged_obj <- seurat_objects[[1]]
# for (i in 2:length(seurat_objects)) {
#   merged_obj <- merge(merged_obj, seurat_objects[[i]],
#                       add.cell.ids = sample_names[i],
#                       project = "Combined_Samples")
# }

# merged_obj2 <- merge(
#   x = seurat_objects[[1]],
#   y = seurat_objects[-1],
#   add.cell.ids = sample_names, # 传递完整的 cell ids 向量
#   project = "Combined_Samples"
# )

# 
# sce.big <- merge(seurat_objects[[1]],
#                  y = c(seurat_objects[[2]],seurat_objects[[3]],seurat_objects[[4]],
#                        seurat_objects[[5]],seurat_objects[[6]],
#                        seurat_objects[[7]],seurat_objects[[8]],seurat_objects[[9]],seurat_objects[[10]],
#                        seurat_objects[[11]],seurat_objects[[12]]),
#                  add.cell.ids = sample_names,
#                  project = "Combined_Samples")


for (i in 1:length(sample_paths)) {
  print(i)
  # 读取 h5 数据
  data <- Read10X_h5(sample_paths[[i]])

  # 创建 Seurat 对象，并为每个样本添加样本信息
  seurat_obj <- CreateSeuratObject(counts = data, project = sample_names[i])

  # 添加样本标注信息，以便区分不同样本
  seurat_obj$sample <- sample_names[i]

  # 存入列表
  seurat_objects[[i]] <- seurat_obj
}
# 
# # 将多个 Seurat 对象合并为一个
combined_seurat_object <- merge(
  seurat_objects[[1]],
  y = seurat_objects[2:length(seurat_objects)],
  add.cell.ids = sample_names,  # 添加样本ID到细胞名称中
  project = "Combined_Samples"
)

combined_seurat_object<-sce.big
table(combined_seurat_object@meta.data$orig.ident)
combined_seurat_object@meta.data$nCount_RNA
rownames(combined_seurat_object@meta.data)
combined_seurat_object@meta.data$cell_id <- rownames(combined_seurat_object@meta.data)

saveRDS(combined_seurat_object,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/多队列验证/GSE200996/GSE200996_RAW/combined_seurat_object减少样本.rds")
# 标准化数据

####加入metadata
metadata1<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/多队列验证/GSE200996/metadata.csv")
head(metadata1)
metadata1$Sample_id
metadata<- FetchData(combined_seurat_object,"sample")
metadata$sample <- gsub("-", "_", metadata$sample)  # 将 - 替换为 _
head(metadata)
metadata_merged<- left_join(x=metadata,y=metadata1,by = join_by(sample==Sample_id))
head(metadata_merged)
rownames(metadata_merged)<-metadata_merged$cell_id

combined_seurat_object<- AddMetaData(combined_seurat_object,metadata = metadata_merged)

View(combined_seurat_object)

saveRDS(combined_seurat_object,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/多队列验证/GSE200996/GSE200996_RAW/combined_seurat_object合并.rds")

combined_seurat_object<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/多队列验证/GSE200996/GSE200996_RAW/combined_seurat_object合并.rds")
combined_seurat_object<- readRDS("./combined_seurat_object减少样本.rds")
pbmc<- combined_seurat_object
scRNA1 = subset(pbmc, downsample=100000,seed=42)
pbmc<-combined_seurat_object

pbmc[["percent.mt"]] = PercentageFeatureSet(pbmc,pattern = "^MT-") #因为人的线粒体基因是MT-开头的，所以可以通过这种方式匹配，如果是小鼠的话，是以“mt-”开头的。其他特殊物种可以直接提供线粒体基因
#percent.HB计算***
HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ") #直接提供血红细胞基因
HB_m <- match(HB.genes, rownames(pbmc@assays$RNA)) #血红细胞基因和seurat对象的基因匹配，返回对应的基因定位信息

HB.genes <- rownames(pbmc@assays$RNA)[HB_m] #将位置信息定位回去
HB.genes <- HB.genes[!is.na(HB.genes)] #去除NA信息
pbmc[["percent.HB"]]<-PercentageFeatureSet(pbmc, features=HB.genes) #根据提供的红细胞基因列表计算红细胞基因表达的比例

#数据过滤前的可视化***
beforeQC_vlnplot = VlnPlot(pbmc, 
                           features = c("nFeature_RNA", 
                                        "nCount_RNA", 
                                        "percent.mt",
                                        "percent.HB"), 
                           raster=FALSE,
                           ncol = 4, 
                           pt.size = 0) #设置对应的图片排版和点的大小
beforeQC_vlnplot


minGene=500 #这里设置了比较严格的条件，一般可以设置200，保留较多的细胞用于后续分析，保留细胞类型多样性，由于今天演示，为保证细胞的基因表达数量足够，设置得较高。
maxGene=3000 
pctMT=10
pbmc = subset(pbmc,subset = nFeature_RNA > minGene & nFeature_RNA < maxGene & percent.mt < pctMT) #取符合条件的子集

afterQC_vlnplot = VlnPlot(pbmc, 
                          features = c("nFeature_RNA", 
                                       "nCount_RNA", 
                                       "percent.mt",
                                       "percent.HB"),
                          raster=FALSE,
                          ncol = 4, 
                          pt.size = 0) #设置对应的图片排版和点的大小

afterQC_vlnplot

pbmc <- NormalizeData(pbmc,normalization.method = "LogNormalize",
                      scale.factor = 1e4)
pbmc <- FindVariableFeatures(pbmc,
                             selection.method = "vst",nfeatures = 2000)

pbmc<- ScaleData(pbmc,features = VariableFeatures(pbmc))
pbmc <- RunPCA(object = pbmc,pc.genes = VariableFeatures(pbmc))

options(future.globals.maxSize = 8 * 1024^3)  # 将最大大小设置为8 GiB
pbmc <- FindNeighbors(pbmc,dims = 1:20)
for (res in c(0.01, 0.05, 0.1, 0.2 ,0.3, 0.5, 0.6,0.8, 1)) {
  pbmc = FindClusters(pbmc,
                      resolution = res, algorithm = 1)
}
pbmc = FindClusters(pbmc,
                    resolution = 0.4, algorithm = 1)
sel.clust <- "RNA_snn_res.0.4"
pbmc <- SetIdent(pbmc,value = sel.clust)
library(harmony)
pbmc <- RunUMAP(pbmc,dims = 1:20)#UMAP二次降维，取前20个维度
pbmc <- RunHarmony(object =pbmc, reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")#harmony需要pca降维的结果，而scRNA对象此前已经做过PCA降维了，所以此处不再重复
pbmc <- RunUMAP(pbmc, reduction = "harmony", dims = 1:20,reduction.name = "umap.harmony")
pbmc <- RunTSNE(object = pbmc,reduction = "harmony", reduction.key = "tsneharmony_",dims.use = 1:20,reduction.name = "tsne.harmony")
pbmc <- FindNeighbors(pbmc, reduction = "harmony", dims = 1:20)
pbmc <- FindClusters(pbmc)#标准聚类


combined_seurat_object<-pbmc
#saveRDS(combined_seurat_object,"./combined_seurat_object减少样本_sampledown10万.rds")
saveRDS(combined_seurat_object,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/多队列验证/GSE200996/GSE200996_RAW/combined_seurat_object合并12个样本.rds")
plot3 = DimPlot(pbmc, reduction = "umap", label=T) #查看clusters在UMAP降维图中的分布
plot3
plot4 = DimPlot(pbmc, reduction = "umap", group.by='orig.ident') #查看每个样本在UMAP降维图中的分布
plot5 = DimPlot(pbmc, reduction = "umap", split.by='orig.ident') #查看每个样本在UMAP降维图中的分面图
plotc <- plot4+plot5
plotc
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

combined_seurat_object<-pbmc

#CD3D CD3E CD3G CD4 CD8A CD8B
combined_seurat_object<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/多队列验证/GSE200996/GSE200996_RAW/combined_seurat_object合并12个样本.rds")

DimPlot(combined_seurat_object, reduction = "umap.harmony", group.by='orig.ident')

p1<-DimPlot(combined_seurat_object, reduction = "umap.harmony", group.by='RNA_snn_res.0.4',label = TRUE, repel = TRUE)+ 
  scale_color_manual(values = palette)+
  xlab("UMAP_1")+
  ylab("UMAP_2")+NoLegend()

p2<- FeaturePlot(combined_seurat_object,features ="CD3D",reduction = "umap.harmony",raster=FALSE)+
  xlab("UMAP_1")+
  ylab("UMAP_2")

p3<-FeaturePlot(combined_seurat_object,features ="CD3E",reduction = "umap.harmony",raster=FALSE)+
  xlab("UMAP_1")+
  ylab("UMAP_2")

p4<-FeaturePlot(combined_seurat_object,features ="CD3G",reduction = "umap.harmony",raster=FALSE)+
  xlab("UMAP_1")+
  ylab("UMAP_2")

p5<-FeaturePlot(combined_seurat_object,features ="CD4",reduction = "umap.harmony",raster=FALSE)+
  xlab("UMAP_1")+
  ylab("UMAP_2")

p6<-FeaturePlot(combined_seurat_object,features ="CD8A",reduction = "umap.harmony",raster=FALSE)+
  xlab("UMAP_1")+
  ylab("UMAP_2")
p7<-FeaturePlot(combined_seurat_object,features ="CD8B",reduction = "umap.harmony",raster=FALSE)+
  xlab("UMAP_1")+
  ylab("UMAP_2")

p<- p1+p2+p3+p4+p5+p6+p7

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/多队列验证/GSE200996/细胞注释/UMAP.pdf",p,width=12,height=12)

library(randomcoloR)
library(ggplot2)
palette <- distinctColorPalette(15)
# > palette
# [1] "#6C8A8C" "#B63FE2" "#E53F91" "#7DAFDE" "#D4D59E" "#86E3D9" "#D8D5D9" "#76E29A" "#D86DCB" "#746DD7" "#DADD5B" "#D78E9A" "#86EA58" "#CDA4DB"
# [15] "#D77F4D"
Idents(combined_seurat_object)<-"RNA_snn_res.0.4"
T_object <- subset(combined_seurat_object,idents=c("1","2","5","0","4","13","6"))
DimPlot(T_object, reduction = "umap.harmony", group.by='RNA_snn_res.0.4',label = TRUE, repel = TRUE)+ scale_color_manual(values = palette)+
  xlab("UMAP_1")+
  ylab("UMAP_2")+NoLegend()

saveRDS(T_object,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/多队列验证/GSE200996/GSE200996_RAW/T_object合并12个样本.rds")


####T_object质控----
pbmc<-T_object

pbmc <- NormalizeData(pbmc,normalization.method = "LogNormalize",
                      scale.factor = 1e4)
pbmc <- FindVariableFeatures(pbmc,
                             selection.method = "vst",nfeatures = 2000)

pbmc<- ScaleData(pbmc,features = VariableFeatures(pbmc))
pbmc <- RunPCA(object = pbmc,pc.genes = VariableFeatures(pbmc))

options(future.globals.maxSize = 8 * 1024^3)  # 将最大大小设置为8 GiB
pbmc <- FindNeighbors(pbmc,dims = 1:20)
for (res in c(0.01, 0.05, 0.1, 0.2 ,0.3, 0.5, 0.6,0.8, 1)) {
  pbmc = FindClusters(pbmc,
                      resolution = res, algorithm = 1)
}
pbmc = FindClusters(pbmc,
                    resolution = 0.4, algorithm = 1)
sel.clust <- "RNA_snn_res.0.8"
pbmc <- SetIdent(pbmc,value = sel.clust)
library(harmony)
pbmc <- RunUMAP(pbmc,dims = 1:20)#UMAP二次降维，取前20个维度
pbmc <- RunHarmony(object =pbmc, reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")#harmony需要pca降维的结果，而scRNA对象此前已经做过PCA降维了，所以此处不再重复
pbmc <- RunUMAP(pbmc, reduction = "harmony", dims = 1:20,reduction.name = "umap.harmony")
pbmc <- RunTSNE(object = pbmc,reduction = "harmony", reduction.key = "tsneharmony_",dims.use = 1:20,reduction.name = "tsne.harmony")
pbmc <- FindNeighbors(pbmc, reduction = "harmony", dims = 1:20)
pbmc <- FindClusters(pbmc)#标准聚类

T_object<-pbmc

saveRDS(T_object,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/多队列验证/GSE200996/GSE200996_RAW/T_object合并12个样本重新聚类0.8.rds")

p1<-DimPlot(T_object, reduction = "umap.harmony", group.by='RNA_snn_res.0.8',label = TRUE, repel = TRUE)+ 
  scale_color_manual(values = palette)+
  xlab("UMAP_1")+
  ylab("UMAP_2")+NoLegend()

p2<- FeaturePlot(T_object,features ="CD3D",reduction = "umap.harmony",raster=FALSE)+
  xlab("UMAP_1")+
  ylab("UMAP_2")

p3<-FeaturePlot(T_object,features ="CD3E",reduction = "umap.harmony",raster=FALSE)+
  xlab("UMAP_1")+
  ylab("UMAP_2")

p4<-FeaturePlot(T_object,features ="CD3G",reduction = "umap.harmony",raster=FALSE)+
  xlab("UMAP_1")+
  ylab("UMAP_2")

p5<-FeaturePlot(T_object,features ="CD4",reduction = "umap.harmony",raster=FALSE)+
  xlab("UMAP_1")+
  ylab("UMAP_2")

p6<-FeaturePlot(T_object,features ="CD8A",reduction = "umap.harmony",raster=FALSE)+
  xlab("UMAP_1")+
  ylab("UMAP_2")
p7<-FeaturePlot(T_object,features ="CD8B",reduction = "umap.harmony",raster=FALSE)+
  xlab("UMAP_1")+
  ylab("UMAP_2")

p<- p1+p2+p3+p4+p5+p6+p7

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/多队列验证/GSE200996/细胞注释/T_object_UMAP.pdf",p,width=12,height=12)

unique(T_object$sample)

CD8_object<-subset(T_object,idents=c("10","7","11","0","14","4","13")) 

saveRDS(CD8_object,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/多队列验证/GSE200996/GSE200996_RAW/CD8_object合并12个样本.rds")

layer_names <- c("counts.P18-pre", "counts.P23-pre", "counts.P24-pre", "counts.P27-pre",
                 "counts.P32-pre", "counts.P13-post", "counts.P18-post", "counts.P21-post",
                 "counts.P23-post", "counts.P24-post", "counts.P27-post", "counts.P32-post")

# 提取各个层的矩阵，并存入列表
count_matrices <- lapply(layer_names, function(layer) {
  GetAssayData(CD8_object[["RNA"]], slot = layer)
})

# 将所有计数矩阵按列合并
merged_counts <- do.call(cbind, count_matrices)

CD8_object[["RNA"]] <- CreateAssayObject(counts = merged_counts)
View(CD8_object)
CD8_object <- NormalizeData(CD8_object)


####添加metadata信息----
library(dplyr)
group<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/多队列验证/GSE200996/metadata.csv")
names(group)
View(CD8_object)
metadata<- FetchData(CD8_object,"sample")
metadata$sample <- gsub("-", "_", metadata$sample)
metadata<- left_join(x=metadata,y=group,by = join_by(sample==Sample_id))
CD8_object<- AddMetaData(CD8_object,metadata = metadata)
View(CD8_object)
table(CD8_object@meta.data$sample %in% group$Sample_id)
head(CD8_object@meta.data)

CD8_object$Path_response
table(CD8_object$Path_response)
Assays(CD8_object)

T_exhausted<-read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/不同分组T亚群拟时序分析和特征得分/T细胞耗竭毒性所用marker.csv")



T_cell_markers<-list(T_exhausted$T.Cell.Exhuastion.signature,T_exhausted$T.Cell.Mediated.Immune.Response.to.Tumor.Cell,
                     T_exhausted$Terminally.Exhausted.CD8..T.Cell.Signature,T_exhausted$Terminally.Exhausted.CD8..T.Cell.Signature_filter,T_exhausted$Progenitor.Exhausted.CD8..T.Cell.Signature,
                     T_exhausted$Cytotoxicity,T_exhausted$Progenitor.Exhausted.CD8..T.Cell.Signature_filter)

# 查看 `CD8_object` 中 RNA assay 的特征名称
available_genes <- rownames(GetAssayData(CD8_object, assay = "RNA", slot = "scale.data"))

# 查看哪些基因不在 `CD8_object` 中
missing_genes <- lapply(T_cell_markers, function(x) setdiff(x, available_genes))
missing_genes
T_cell_markers <- lapply(T_cell_markers, function(x) x[x %in% available_genes])

imm_T<-AddModuleScore(CD8_object,features=T_cell_markers,name=c("T.Cell.Exhuastion.signature",
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

VlnPlot(imm_T,features=c("Terminally.Exhausted.CD8..T.Cell.Signature_filter4","Progenitor.Exhausted.CD8..T.Cell.Signature_filter7"),pt.size=0,group.by="Path_response",ncol=1)


###CD8A,CD8B,PDCD1, and GNLY

CD8_object <- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/多队列验证/GSE200996/GSE200996_RAW/CD8_object合并12个样本.rds")
CD8_object@meta.data$Stage
CD8_object@meta.data$Cohort
unique(CD8_object@meta.data$Cohort)
CD8_object_Mono <- subset(CD8_object, subset = Cohort == "Mono")
unique(CD8_object_Mono@meta.data$Cohort)

CD8_object_Combo <- subset(CD8_object, subset = Cohort == "Combo")

unique(CD8_object_Combo@meta.data$Cohort)


gene_list<-list(c("CD8A","CD8B","IL7R","CD69","FOS","CXCR4"))
gene_list<-list(c("IL7R","MKI67","PDCD1","ZNF683","HOPX","ITGAE","TCF7","SLAMF6","CD8A","CD8B","PDCD1","GNLY"))

gene_list<-list(c("CD8A","CD8B","IL7R"))
gene_list<-list(c("CD8A","CD8B","GZMK"))
gene_list<-list(c("CD8A","CD8B","GNLY"))
gene_list<-list(c("CD8A","CD8B","CTSW"))
gene_list<-list(c("CD8A","CD8B","STMN1"))

gene_list<-list(c("IL7R","GPR183","LMNA","NR4A3","TCF7","MGAT4A","CD55","AIM1","PER1",
"FOSL2","EGR1","TSPYL2","YPEL5","CSRNP1","REL","SKIL","PIK3R1","FOXP1",
"RGCC","PFKFB3","MYADM","ZFP36L2","USP36","TC2N","FAM177A1","BTG2","TSC22D2","FAM65B",
"STAT4","RGPD5","NEU1","IFRD1","PDE4B","NR4A1")) ####前体耗竭

gene_list<-list(c("PDCD1",
                  "LAYN",
                  "HAVCR2",
                  "LAG3",
                  "CTLA4",
                  "TIGIT",
                  "TOX",
                  "RBPJ",
                  "VCAM1",
                  "GZMB",
                  "MYO7A",
                  "CD244",
                  "VSIR",
                  "BTLA",
                  "ENTPD1",
                  "CD160",
                  "LAIR1"
))



imm_T2<-AddModuleScore(CD8_object_Mono,features=gene_list,name="CD8_IL7R_signature")

imm_T2<-AddModuleScore(CD8_object,features=gene_list,name="CD8_IL7R_signature")

imm_T2<-AddModuleScore(CD8_object_Combo,features=gene_list,name="CD8_IL7R_signature")

head(imm_T2)

library(tidyverse)
library(gghalves)
library(rstatix)
library(ggpubr)
library(ggsci)

sig_score<- cbind(imm_T2$sample,imm_T2$Path_response,imm_T2$Stage,imm_T2$CD8_IL7R_signature1) %>% as.data.frame()

names(sig_score)<- c("sample","Path_response","Stage","CD8_IL7R_signature1")
unique(sig_score$Path_response)
df<-sig_score
df$CD8_IL7R_signature1 <- as.numeric(df$CD8_IL7R_signature1)
names(df)

df_p_val1 <- df %>% 
  wilcox_test(CD8_IL7R_signature1~ Path_response) %>%
  adjust_pvalue(p.col = "p", method = "bonferroni") %>%
  add_significance(p.col = "p.adj") %>%
  add_xy_position(x = "Path_response", dodge = 0.8)


response_color <- c("#FF9a9b","#9998ff","#c99bff")

#treatment_color <- c("#4974a4","#f29600")
unique(df$Stage)
unique(df$Path_response)
df$Path_response<- factor(df$Path_response,levels = c("Low","Medium","High"))
df$Stage<- factor(df$Stage,levels = c("Pre-Tx","Post-Tx"))
df_p_val1 <- df %>% group_by(Stage) %>%
  wilcox_test(CD8_IL7R_signature1~ Path_response) %>%
  adjust_pvalue(p.col = "p", method = "bonferroni") %>%
  add_significance(p.col = "p.adj") %>%
  add_xy_position(x = "Path_response", dodge = 0.8)
View(df)
t_exhausted<- df %>%
  ggplot(aes(Path_response,CD8_IL7R_signature1)) +
  geom_half_boxplot(aes(fill=Path_response),color="black",side="l",errorbar.draw = T,
                    outlier.shape = NA,width=0.8,lwd= 0.7) +
  geom_half_point(aes(),side = "r")+
  stat_pvalue_manual(df_p_val1,label = "p.adj.signif",label.size=5,hide.ns = F,bracket.size = 0.5,y.position = c(2.6,2.8,3.2))+
  facet_wrap(.~Stage,nrow=1)+
  #scale_y_continuous(limits = c(0,5),breaks = seq(0,2,5))+
  scale_fill_manual(values= c("High" = "#FF9a9b",
                              "Medium" = "#9998ff",
                              "Low"="#c99bff"))+
  
  # scale_fill_manual(values= c("Medium" = "#9998ff",
  #                             "Low"="#c99bff"))+
  scale_color_npg()+
  labs(title = "CD8_STMN1_signature score",y='Signature score',x="")+
  #labs(title ="Progenitor exhausted CD8 T cell signature",x="")+
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
ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/多队列验证/GSE200996/特征得分/GSE200996_CD8_STMN1_combo_signature治疗前后分面.pdf",height = 5,width=7)


df_p_val1 <- df %>% 
  wilcox_test(CD8_IL7R_signature1~ Stage) %>%
  adjust_pvalue(p.col = "p", method = "bonferroni") %>%
  add_significance(p.col = "p.adj") %>%
  add_xy_position(x = "Stage", dodge = 0.8)
treatment_color<-c("#FF9a9b","#c99bff")

df$Stage<- factor(df$Stage,levels=c("Pre-Tx","Post-Tx"))

t_exhausted<- df %>%
  ggplot(aes(Stage,CD8_IL7R_signature1)) +
  geom_half_boxplot(aes(fill=Stage),color="black",side="l",errorbar.draw = T,
                    outlier.shape = NA,width=0.8,lwd= 0.7) +
  geom_half_point(color=treatment_color,side = "r",alpha=0.1,
                  transformation_params = list(height = 0,width = 0.001,seed = 2))+
  stat_pvalue_manual(df_p_val1,label = "p.adj.signif",label.size=5,hide.ns = F,bracket.size = 0.5,y.position = c(2.5))+
  #scale_y_continuous(limits = c(0,5),breaks = seq(0,2,5))+
  scale_fill_manual(values= c("Pre-Tx" = "#FF9a9b",
                              "Post-Tx"="#c99bff"))+
  
  # scale_fill_manual(values= c("Medium" = "#9998ff",
  #                             "Low"="#c99bff"))+
  scale_color_npg()+
  labs(title = "CD8_STMN1_signature score",y='Signature score',x="")+
  #labs(title ="Progenitor exhausted CD8 T cell signature",x="")+
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

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/多队列验证/GSE200996/特征得分/GSE200996_CD8_STMN1_combo_signature.pdf",height = 5,width=5)


#####CTSW----

CD8_object <- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/多队列验证/GSE200996/GSE200996_RAW/CD8_object合并12个样本.rds")

gene_list<-list(c("CD8A","CD8B","IL7R","CD69","FOS","CXCR4"))
gene_list<-list(c("IL7R","MKI67","PDCD1","ZNF683","HOPX","ITGAE","TCF7","SLAMF6","CD8A","CD8B","PDCD1","GNLY"))

gene_list<-list(c("CD8A","CD8B","PDCD1","GNLY"))
imm_T2<-AddModuleScore(CD8_object,features=gene_list,name="CD8_GNLY_signature")
head(imm_T2)

library(tidyverse)
library(gghalves)
library(rstatix)
library(ggpubr)
library(ggsci)

sig_score<- cbind(imm_T2$sample,imm_T2$Path_response,imm_T2$Stage,imm_T2$CD8_GNLY_signature1) %>% as.data.frame()

names(sig_score)<- c("sample","Path_response","Stage","CD8_GNLY_signature1")

df<-sig_score
df$CD8_GNLY_signature1 <- as.numeric(df$CD8_GNLY_signature1)
names(df)
df_p_val1 <- df %>% group_by(Stage) %>%
  wilcox_test(CD8_GNLY_signature1~ Path_response) %>%
  adjust_pvalue(p.col = "p", method = "bonferroni") %>%
  add_significance(p.col = "p.adj") %>%
  add_xy_position(x = "Path_response", dodge = 0.8)
df_p_val1 <- df %>% 
  wilcox_test(CD8_GNLY_signature1~ Path_response) %>%
  adjust_pvalue(p.col = "p", method = "bonferroni") %>%
  add_significance(p.col = "p.adj") %>%
  add_xy_position(x = "Path_response", dodge = 0.8)


response_color <- c("#FF9a9b","#9998ff","#c99bff")

treatment_color <- c("#4974a4","#f29600")
unique(df$Stage)
df$Path_response<- factor(df$Path_response,levels = c("Low","Medium","High"))
df$Stage<- factor(df$Stage,levels = c("Pre-Tx","Post-Tx"))

View(df)
t_exhausted<- df %>%
  ggplot(aes(Path_response,CD8_GNLY_signature1)) +
  geom_half_boxplot(aes(fill=Path_response),color="black",side="l",errorbar.draw = T,
                    outlier.shape = NA,width=0.8,lwd= 0.7) +
  geom_half_point(aes(),side = "r")+
  stat_pvalue_manual(df_p_val1,label = "p.adj.signif",label.size=5,hide.ns = F,bracket.size = 0.5,y.position = c(2.4,2.6,2.8))+
  facet_wrap(.~Stage,nrow=1)+
  #scale_y_continuous(limits = c(0,5),breaks = seq(0,2,5))+
  scale_fill_manual(values= c("High" = "#FF9a9b",
                              "Medium" = "#9998ff",
                              "Low"="#c99bff"))+
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
ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/多队列验证/GSE200996/特征得分/CD8_GNLY_signature.pdf",height = 5,width=7)



######TCGA 生存分析数据下载----
install.packages("BiocManager")
BiocManager::install("TCGAbiolinks")
library(TCGAbiolinks)

# 查询HNSCC RNA-Seq数据
query_bulk <- GDCquery(project = "TCGA-HNSC",
                       data.category = "Transcriptome Profiling",
                       data.type = "Gene Expression Quantification",
                       workflow.type = "STAR - Counts")

# 下载数据
GDCdownload(query_bulk)
bulk_data <- GDCprepare(query_bulk)

# 查看数据
head(bulk_data)
View(bulk_data)
bulk_data<- saveRDS(bulk_data,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/多队列验证/生存分析/bulk_data.rds")



# 查询HNSCC临床数据
query_clinical <- GDCquery(project = "TCGA-HNSC",
                           data.category = "Clinical",
                           data.type = "Clinical Supplement")

# 下载数据
GDCdownload(query_clinical)
clinical_data <- GDCprepare(query_clinical)

# 查看临床数据
head(clinical_data)

clinical_data<- saveRDS(clinical_data,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/多队列验证/生存分析/clinical_data.rds")
clinical_data<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/多队列验证/生存分析/clinical_data.rds")
View(clinical_data)

# 使用GEO或Single Cell Portal
# 访问GEO： 使用GEO数据库搜索HNSCC相关的scRNA-seq数据，网址：GEO
# 
# 访问Single Cell Portal： Single Cell Portal是一个很好的资源，可以搜索与HNSCC相关的scRNA-seq数据。

#获取到的bulk RNA和临床数据可以用来进行生存分析，使用生存分析包（如survival）进行分析和可视化。


##########EGAS50000000037----
###读取矩阵数据
library(Seurat)
data_dir <- "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/多队列验证/EGAS50000000037/5577-HNSCC_Tissue/"  # 数据目录
sc_data <- Read10X(data.dir = data_dir)

# 创建 Seurat 对象
seurat_obj <- CreateSeuratObject(counts = sc_data)

View(seurat_obj)

pbmc<- seurat_obj 

pbmc[["percent.mt"]] = PercentageFeatureSet(pbmc,pattern = "^MT-") #因为人的线粒体基因是MT-开头的，所以可以通过这种方式匹配，如果是小鼠的话，是以“mt-”开头的。其他特殊物种可以直接提供线粒体基因
#percent.HB计算***
HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ") #直接提供血红细胞基因
HB_m <- match(HB.genes, rownames(pbmc@assays$RNA)) #血红细胞基因和seurat对象的基因匹配，返回对应的基因定位信息

HB.genes <- rownames(pbmc@assays$RNA)[HB_m] #将位置信息定位回去
HB.genes <- HB.genes[!is.na(HB.genes)] #去除NA信息
pbmc[["percent.HB"]]<-PercentageFeatureSet(pbmc, features=HB.genes) #根据提供的红细胞基因列表计算红细胞基因表达的比例

#数据过滤前的可视化***
beforeQC_vlnplot = VlnPlot(pbmc, 
                           features = c("nFeature_RNA", 
                                        "nCount_RNA", 
                                        "percent.mt",
                                        "percent.HB"), 
                           raster=FALSE,
                           ncol = 4, 
                           pt.size = 0) #设置对应的图片排版和点的大小
beforeQC_vlnplot


minGene=500 #这里设置了比较严格的条件，一般可以设置200，保留较多的细胞用于后续分析，保留细胞类型多样性，由于今天演示，为保证细胞的基因表达数量足够，设置得较高。
maxGene=3000 
pctMT=10
pbmc = subset(pbmc,subset = nFeature_RNA > minGene & nFeature_RNA < maxGene & percent.mt < pctMT) #取符合条件的子集

afterQC_vlnplot = VlnPlot(pbmc, 
                          features = c("nFeature_RNA", 
                                       "nCount_RNA", 
                                       "percent.mt",
                                       "percent.HB"),
                          raster=FALSE,
                          ncol = 4, 
                          pt.size = 0) #设置对应的图片排版和点的大小

afterQC_vlnplot

pbmc <- NormalizeData(pbmc,normalization.method = "LogNormalize",
                      scale.factor = 1e4)
pbmc <- FindVariableFeatures(pbmc,
                             selection.method = "vst",nfeatures = 2000)

pbmc<- ScaleData(pbmc,features = VariableFeatures(pbmc))
pbmc <- RunPCA(object = pbmc,pc.genes = VariableFeatures(pbmc))

options(future.globals.maxSize = 8 * 1024^3)  # 将最大大小设置为8 GiB
pbmc <- FindNeighbors(pbmc,dims = 1:20)
for (res in c(0.01, 0.05, 0.1, 0.2 ,0.3, 0.5, 0.6,0.8, 1)) {
  pbmc = FindClusters(pbmc,
                      resolution = res, algorithm = 1)
}
pbmc = FindClusters(pbmc,
                    resolution = 0.4, algorithm = 1)
sel.clust <- "RNA_snn_res.0.4"
pbmc <- SetIdent(pbmc,value = sel.clust)
library(harmony)
pbmc <- RunUMAP(pbmc,dims = 1:20)#UMAP二次降维，取前20个维度
library(dplyr)
library(tidyr)
pbmc$sampleid <- sub("(^[^_]*_[^_]*).*", "\\1", rownames(pbmc@meta.data))

pbmc <- RunHarmony(object =pbmc, reduction = "pca",group.by.vars = "sampleid",reduction.save = "harmony")#harmony需要pca降维的结果，而scRNA对象此前已经做过PCA降维了，所以此处不再重复
pbmc <- RunUMAP(pbmc, reduction = "harmony", dims = 1:20,reduction.name = "umap.harmony")
pbmc <- RunTSNE(object = pbmc,reduction = "harmony", reduction.key = "tsneharmony_",dims.use = 1:20,reduction.name = "tsne.harmony")
pbmc <- FindNeighbors(pbmc, reduction = "harmony", dims = 1:20)
pbmc <- FindClusters(pbmc)#标准聚类

seurat_obj<-pbmc

saveRDS(seurat_obj,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/多队列验证/EGAS50000000037/seurat_obj.rds")

DimPlot(seurat_obj, reduction = "umap.harmony", group.by='orig.ident',raster=FALSE)

library(ggplot2)
library(randomcoloR)
palette<- distinctColorPalette(22)

p1<-DimPlot(seurat_obj, reduction = "umap.harmony", group.by='RNA_snn_res.0.4',label = TRUE, raster=FALSE,repel = TRUE)+ 
  scale_color_manual(values = palette)+
  xlab("UMAP_1")+
  ylab("UMAP_2")+NoLegend()

p2<- FeaturePlot(seurat_obj,features ="CD3D",reduction = "umap.harmony",raster=FALSE)+
  xlab("UMAP_1")+
  ylab("UMAP_2")

p3<-FeaturePlot(seurat_obj,features ="CD3E",reduction = "umap.harmony",raster=FALSE)+
  xlab("UMAP_1")+
  ylab("UMAP_2")

p4<-FeaturePlot(seurat_obj,features ="CD3G",reduction = "umap.harmony",raster=FALSE)+
  xlab("UMAP_1")+
  ylab("UMAP_2")

p5<-FeaturePlot(seurat_obj,features ="CD4",reduction = "umap.harmony",raster=FALSE)+
  xlab("UMAP_1")+
  ylab("UMAP_2")

p6<-FeaturePlot(seurat_obj,features ="CD8A",reduction = "umap.harmony",raster=FALSE)+
  xlab("UMAP_1")+
  ylab("UMAP_2")
p7<-FeaturePlot(seurat_obj,features ="CD8B",reduction = "umap.harmony",raster=FALSE)+
  xlab("UMAP_1")+
  ylab("UMAP_2")

p<- p1+p2+p3+p4+p5+p6+p7

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/多队列验证/EGAS50000000037/细胞注释/UMAP.pdf",p,width=12,height=12)


metadata1<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/多队列验证/EGAS50000000037/6010-Immunity_metadata2.csv")
head(metadata1)
metadata1$Sample
metadata<- FetchData(seurat_obj,"sampleid")
rownames(metadata)
metadata_merged<- left_join(x=metadata,y=metadata1,by = join_by(sampleid==Sample))
head(metadata_merged)
rownames(metadata_merged)<-rownames(metadata)
seurat_obj<- AddMetaData(seurat_obj,metadata = metadata_merged)
View(seurat_obj)
saveRDS(seurat_obj,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/多队列验证/EGAS50000000037/seurat_obj.rds")

head(seurat_obj@meta.data)

seurat_obj<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/多队列验证/EGAS50000000037/seurat_obj.rds")

names(seurat_obj@meta.data)

Idents(seurat_obj)<-"RNA_snn_res.0.4"

sub_object_9 = subset(seurat_obj,idents=c("9"))

sub_object_9$RNA_snn_res.0.4 = ifelse(sub_object_9@assays$RNA$counts['CD3D',]>0& sub_object_9@assays$RNA$counts['CD3E',]>0,'9_T','9_neg')

seurat_obj$RNA_snn_res.0.4 <- as.character(seurat_obj$RNA_snn_res.0.4)
sub_object_9$RNA_snn_res.0.4 <- as.character(sub_object_9$RNA_snn_res.0.4)

seurat_obj$RNA_snn_res.0.4[match(colnames(sub_object_9),colnames(seurat_obj))]=sub_object_9$RNA_snn_res.0.4

table(seurat_obj$RNA_snn_res.0.4)

Idents(seurat_obj)<-"RNA_snn_res.0.4"


sub_object_1 = subset(seurat_obj,idents=c("1"))

sub_object_1$RNA_snn_res.0.4 = ifelse(sub_object_1@assays$RNA$counts['CD3D',]>0& sub_object_1@assays$RNA$counts['CD3E',]>0,'1_T','1_neg')

seurat_obj$RNA_snn_res.0.4 <- as.character(seurat_obj$RNA_snn_res.0.4)
sub_object_1$RNA_snn_res.0.4 <- as.character(sub_object_1$RNA_snn_res.0.4)

seurat_obj$RNA_snn_res.0.4[match(colnames(sub_object_1),colnames(seurat_obj))]=sub_object_1$RNA_snn_res.0.4

table(seurat_obj$RNA_snn_res.0.4)

Idents(seurat_obj)<-"RNA_snn_res.0.4"

Idents(seurat_obj)<-"RNA_snn_res.0.4"
T_object <- subset(seurat_obj,idents=c("0","3","1_T","9_T"))
DimPlot(T_object, reduction = "umap.harmony", group.by='RNA_snn_res.0.4',raster=FALSE,label = TRUE, repel = TRUE)+ scale_color_manual(values = palette)+
  xlab("UMAP_1")+
  ylab("UMAP_2")+NoLegend()

saveRDS(T_object,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/多队列验证/EGAS50000000037/T_object.rds")


pbmc<-T_object

pbmc <- NormalizeData(pbmc,normalization.method = "LogNormalize",
                      scale.factor = 1e4)
pbmc <- FindVariableFeatures(pbmc,
                             selection.method = "vst",nfeatures = 2000)

pbmc<- ScaleData(pbmc,features = VariableFeatures(pbmc))
pbmc <- RunPCA(object = pbmc,pc.genes = VariableFeatures(pbmc))

options(future.globals.maxSize = 8 * 1024^3)  # 将最大大小设置为8 GiB
pbmc <- FindNeighbors(pbmc,dims = 1:20)
for (res in c(0.01, 0.05, 0.1, 0.2 ,0.3, 0.5, 0.6,0.8, 1)) {
  pbmc = FindClusters(pbmc,
                      resolution = res, algorithm = 1)
}
pbmc = FindClusters(pbmc,
                    resolution = 0.4, algorithm = 1)
sel.clust <- "RNA_snn_res.0.8"
pbmc <- SetIdent(pbmc,value = sel.clust)
library(harmony)
pbmc <- RunUMAP(pbmc,dims = 1:20)#UMAP二次降维，取前20个维度
pbmc <- RunHarmony(object =pbmc, reduction = "pca",group.by.vars = "sampleid",reduction.save = "harmony")#harmony需要pca降维的结果，而scRNA对象此前已经做过PCA降维了，所以此处不再重复
pbmc <- RunUMAP(pbmc, reduction = "harmony", dims = 1:20,reduction.name = "umap.harmony")
pbmc <- RunTSNE(object = pbmc,reduction = "harmony", reduction.key = "tsneharmony_",dims.use = 1:20,reduction.name = "tsne.harmony")
pbmc <- FindNeighbors(pbmc, reduction = "harmony", dims = 1:20)
pbmc <- FindClusters(pbmc)#标准聚类

T_object<-pbmc

saveRDS(T_object,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/多队列验证/EGAS50000000037/T_object.rds")


p1<-DimPlot(T_object, reduction = "umap.harmony", group.by='RNA_snn_res.0.8',label = TRUE, repel = TRUE)+ 
  scale_color_manual(values = palette)+
  xlab("UMAP_1")+
  ylab("UMAP_2")+NoLegend()

p2<- FeaturePlot(T_object,features ="CD3D",reduction = "umap.harmony",raster=FALSE)+
  xlab("UMAP_1")+
  ylab("UMAP_2")

p3<-FeaturePlot(T_object,features ="CD3E",reduction = "umap.harmony",raster=FALSE)+
  xlab("UMAP_1")+
  ylab("UMAP_2")

p4<-FeaturePlot(T_object,features ="CD3G",reduction = "umap.harmony",raster=FALSE)+
  xlab("UMAP_1")+
  ylab("UMAP_2")

p5<-FeaturePlot(T_object,features ="CD4",reduction = "umap.harmony",raster=FALSE)+
  xlab("UMAP_1")+
  ylab("UMAP_2")

p6<-FeaturePlot(T_object,features ="CD8A",reduction = "umap.harmony",raster=FALSE)+
  xlab("UMAP_1")+
  ylab("UMAP_2")
p7<-FeaturePlot(T_object,features ="CD8B",reduction = "umap.harmony",raster=FALSE)+
  xlab("UMAP_1")+
  ylab("UMAP_2")

p<- p1+p2+p3+p4+p5+p6+p7

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/多队列验证/EGAS50000000037/细胞注释/T_object_UMAP.jpg",p,width=12,height=12)
ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/多队列验证/EGAS50000000037/细胞注释/T_object_UMAP.pdf",p,width=12,height=12)

unique(T_object$sample)


sub_object_12 = subset(T_object,idents=c("12"))

sub_object_12$RNA_snn_res.0.8 = ifelse(sub_object_12@assays$RNA$counts['CD8A',]>0& sub_object_12@assays$RNA$counts['CD8B',]>0,'12_cd8','12_neg')

T_object$RNA_snn_res.0.8 <- as.character(T_object$RNA_snn_res.0.8)
sub_object_12$RNA_snn_res.0.8 <- as.character(sub_object_12$RNA_snn_res.0.8)

T_object$RNA_snn_res.0.8[match(colnames(sub_object_12),colnames(T_object))]=sub_object_12$RNA_snn_res.0.8

table(T_object$RNA_snn_res.0.8)

Idents(T_object)<-"RNA_snn_res.0.8"


sub_object_20 = subset(T_object,idents=c("20"))

sub_object_20$RNA_snn_res.0.8 = ifelse(sub_object_20@assays$RNA$counts['CD8A',]>0& sub_object_20@assays$RNA$counts['CD8B',]>0,'20_cd8','20_neg')

T_object$RNA_snn_res.0.8 <- as.character(T_object$RNA_snn_res.0.8)
sub_object_20$RNA_snn_res.0.8 <- as.character(sub_object_20$RNA_snn_res.0.8)

T_object$RNA_snn_res.0.8[match(colnames(sub_object_20),colnames(T_object))]=sub_object_20$RNA_snn_res.0.8

table(T_object$RNA_snn_res.0.8)

Idents(T_object)<-"RNA_snn_res.0.8"

CD8_object<- subset(T_object,idents=c("11","20_cd8","19","6","12_cd8","10","2","17"))
DimPlot(CD8_object, reduction = "umap.harmony", group.by='RNA_snn_res.0.8',raster=FALSE,label = TRUE, repel = TRUE)+ scale_color_manual(values = palette)+
  xlab("UMAP_1")+
  ylab("UMAP_2")+NoLegend()



pbmc<-CD8_object

pbmc <- NormalizeData(pbmc,normalization.method = "LogNormalize",
                      scale.factor = 1e4)
pbmc <- FindVariableFeatures(pbmc,
                             selection.method = "vst",nfeatures = 2000)

pbmc<- ScaleData(pbmc,features = VariableFeatures(pbmc))
pbmc <- RunPCA(object = pbmc,pc.genes = VariableFeatures(pbmc))

options(future.globals.maxSize = 8 * 1024^3)  # 将最大大小设置为8 GiB
pbmc <- FindNeighbors(pbmc,dims = 1:20)
for (res in c(0.01, 0.05, 0.1, 0.2 ,0.3, 0.5, 0.6,0.8, 1)) {
  pbmc = FindClusters(pbmc,
                      resolution = res, algorithm = 1)
}
pbmc = FindClusters(pbmc,
                    resolution = 0.4, algorithm = 1)
sel.clust <- "RNA_snn_res.0.8"
pbmc <- SetIdent(pbmc,value = sel.clust)
library(harmony)
pbmc <- RunUMAP(pbmc,dims = 1:20)#UMAP二次降维，取前20个维度
pbmc <- RunHarmony(object =pbmc, reduction = "pca",group.by.vars = "sampleid",reduction.save = "harmony")#harmony需要pca降维的结果，而scRNA对象此前已经做过PCA降维了，所以此处不再重复
pbmc <- RunUMAP(pbmc, reduction = "harmony", dims = 1:20,reduction.name = "umap.harmony")
pbmc <- RunTSNE(object = pbmc,reduction = "harmony", reduction.key = "tsneharmony_",dims.use = 1:20,reduction.name = "tsne.harmony")
pbmc <- FindNeighbors(pbmc, reduction = "harmony", dims = 1:20)
pbmc <- FindClusters(pbmc)#标准聚类

CD8_object<-pbmc

saveRDS(CD8_object,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/多队列验证/EGAS50000000037/CD8_object.rds")

CD8_object<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/多队列验证/EGAS50000000037/CD8_object.rds")


T_object <- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/多队列验证/EGAS50000000037/T_object.rds")

###CD8A,CD8B,PDCD1, and GNLY

gene_list<-list(c("CD8A","CD8B","PDCD1","GNLY"))

gene_list<-list(c("IL7R","MKI67","PDCD1","ZNF683","HOPX","ITGAE","TCF7","SLAMF6","CD8A","CD8B","PDCD1","GNLY"))

gene_list<-list(c("CD69",
                  "PDCD1",
                  "HAVCR2",
                  "ENTPD1",
                  "CTLA4",
                  "LAG3",
                  "TIGIT",
                  "TBX21",
                  "CD8A","CD8B","PDCD1","GNLY"))

gene_list<-list(c("CD8A","CD8B","IL7R","CD69","FOS","CXCR4"))

gene_list<-list(c("CD8A","CD8B","IL7R"))

gene_list<-list(c("CD8A","CD8B","GZMK"))

gene_list<-list(c("CD8A","CD8B","IL7R"))
gene_list<-list(c("CD8A","CD8B","GZMK"))

gene_list<-list(c("CD8A","CD8B","GNLY"))

gene_list<-list(c("CD8A","CD8B","CTSW"))

gene_list<-list(c("CD8A","CD8B","STMN1"))


imm_T2<-AddModuleScore(CD8_object,features=gene_list,name="CD8_IL7R_signature")
head(imm_T2)

library(tidyverse)
library(gghalves)
library(rstatix)
library(ggpubr)
library(ggsci)

sig_score<- cbind(imm_T2$sampleid,imm_T2$Tissue,imm_T2$Immunotherapy,imm_T2$Patient,imm_T2$CD8_IL7R_signature1) %>% as.data.frame()

names(sig_score)<- c("sampleid","Tissue","Immunotherapy","Patient","CD8_IL7R_signature1")

df<-sig_score
df$CD8_IL7R_signature1 <- as.numeric(df$CD8_IL7R_signature1)
names(df)
unique(df$Immunotherapy)
df <- na.omit(df)
df1<-df %>% filter(Tissue %in% c("Tumor pre-treatment","Tumor on-treatment")& Immunotherapy %in% c("Durvalumab + Tremelimumab") )

df_p_val1 <- df1 %>%
  wilcox_test(CD8_IL7R_signature1~ Tissue) %>%
  adjust_pvalue(p.col = "p", method = "bonferroni") %>%
  add_significance(p.col = "p.adj") %>%
  add_xy_position(x = "Tissue", dodge = 0.8)
treatment_color <- c("#4974a4","#4dae47","#f29600","#F39B7Fcc")
treatment_color <- c("#4974a4","#f29600")
treatment_color <- c("#4974a4","#4dae47","#f29600")
unique(df1$Tissue)
df1$Tissue<- factor(df1$Tissue,levels = c("Tumor pre-treatment","Tumor on-treatment","Normal adjacent","tdLN"))
View(df1)

t_exhausted<- df1 %>%
  ggplot(aes(Tissue,CD8_IL7R_signature1)) +
  geom_half_boxplot(fill=treatment_color,color="black",side="l",errorbar.draw = T,
                    outlier.shape = NA,width=0.8,lwd= 0.7) +
  geom_half_point(color=treatment_color,side = "r",alpha=0.1,
                  transformation_params = list(height = 0,width = 0.001,seed = 2))+
  stat_pvalue_manual(df_p_val1,label = "p.adj.signif",label.size=5,hide.ns = F,bracket.size = 0.5,y.position = c(2.7))+
  #facet_wrap(.~continent,nrow=1)+
  #scale_y_continuous(limits = c(0,5),breaks = seq(0,2,5))+
  scale_fill_npg()+
  scale_color_npg()+
  labs(title = "CD8_STMN1_signature score(Durvalumab + Tremelimumab)",y='Signature score',x="")+
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

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/多队列验证/EGAS50000000037/特征得分/CD8_STMN1_signature_Durvalumab + Tremelimumab.pdf",t_exhausted,height = 5,width=5)


df1<-df %>% filter(Tissue %in% c("Tumor pre-treatment","Tumor on-treatment")& Immunotherapy %in% c("Durvalumab + Tremelimumab") )

df_p_val1 <- df1 %>%
  wilcox_test(CD8_IL7R_signature1~ Tissue) %>%
  adjust_pvalue(p.col = "p", method = "bonferroni") %>%
  add_significance(p.col = "p.adj") %>%
  add_xy_position(x = "Tissue", dodge = 0.8)
treatment_color <- c("#4974a4","#4dae47","#f29600","#F39B7Fcc")
treatment_color <- c("#4974a4","#f29600")
treatment_color <- c("#4974a4","#4dae47","#f29600")
unique(df1$Tissue)
df1$Tissue<- factor(df1$Tissue,levels = c("Tumor pre-treatment","Tumor on-treatment","Normal adjacent","tdLN"))
View(df1)


df1$Tissue<- factor(df1$Tissue,levels = c("Tumor pre-treatment","Tumor on-treatment","Normal adjacent","tdLN"))
View(df1)

t_p2<- df1 %>%
  ggplot(aes(Tissue,CD8_IL7R_signature1)) +
  geom_half_boxplot(fill=treatment_color,color="black",side="l",errorbar.draw = T,
                    outlier.shape = NA,width=0.8,lwd= 0.7) +
  geom_half_point(color=treatment_color,side = "r",alpha=0.1,
                  transformation_params = list(height = 0,width = 0.001,seed = 2))+
  stat_pvalue_manual(df_p_val1,label = "p.adj.signif",label.size=5,hide.ns = F,bracket.size = 0.5,y.position = c(3))+
  #facet_wrap(.~continent,nrow=1)+
  #scale_y_continuous(limits = c(0,5),breaks = seq(0,2,5))+
  scale_fill_npg()+
  scale_color_npg()+
  labs(title = "CD8_GZMK_signature score(Durvalumab + Tremelimumab)",y='Signature score',x="")+
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
t_p2

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/多队列验证/EGAS50000000037/特征得分/CD8_GZMK_signature_Durvalumab_Tremelimumab.pdf",t_p2,height = 5,width=5)


#####CTSW----
gene_list<-list(c("CD8A","CD8B","PDCD1","GNLY"))
gene_list<-list(c("IL7R","MKI67","PDCD1","ZNF683","HOPX","ITGAE","TCF7","SLAMF6","CD8A","CD8B","PDCD1","GNLY"))

imm_T2<-AddModuleScore(CD8_object,features=gene_list,name="CD8_GNLY_signature")
head(imm_T2)

library(tidyverse)
library(gghalves)
library(rstatix)
library(ggpubr)
library(ggsci)

sig_score<- cbind(imm_T2$sampleid,imm_T2$Tissue,imm_T2$Immunotherapy,imm_T2$Patient,imm_T2$CD8_GNLY_signature1) %>% as.data.frame()

names(sig_score)<- c("sampleid","Tissue","Immunotherapy","Patient","CD8_GNLY_signature1")

df<-sig_score
df$CD8_GNLY_signature1 <- as.numeric(df$CD8_GNLY_signature1)
names(df)
df <- na.omit(df)
df1<-df %>% filter(Tissue %in% c("Tumor pre-treatment","Tumor on-treatment")& Immunotherapy %in% c("Durvalumab") )

df_p_val1 <- df1 %>%
  wilcox_test(CD8_GNLY_signature1~ Tissue) %>%
  adjust_pvalue(p.col = "p", method = "bonferroni") %>%
  add_significance(p.col = "p.adj") %>%
  add_xy_position(x = "Tissue", dodge = 0.8)
treatment_color <- c("#4974a4","#4dae47","#f29600","#F39B7Fcc")
treatment_color <- c("#4974a4","#f29600")
treatment_color <- c("#4974a4","#4dae47","#f29600")
unique(df1$Tissue)
df1$Tissue<- factor(df1$Tissue,levels = c("Tumor pre-treatment","Tumor on-treatment","Normal adjacent","tdLN"))
View(df1)

t_exhausted<- df1 %>%
  ggplot(aes(Tissue,CD8_GNLY_signature1)) +
  geom_half_boxplot(fill=treatment_color,color="black",side="l",errorbar.draw = T,
                    outlier.shape = NA,width=0.8,lwd= 0.7) +
  geom_half_point(color=treatment_color,side = "r",alpha=0.1,
                  transformation_params = list(height = 0,width = 0.001,seed = 2))+
  stat_pvalue_manual(df_p_val1,label = "p.adj.signif",label.size=5,hide.ns = F,bracket.size = 0.5,y.position = c(2.1))+
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

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/多队列验证/GSE200996/特征得分/CD8_GNLY_signature.pdf",height = 5,width=5)



