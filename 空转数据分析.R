
dir.create("~/SeuratV4")
# 然后安装的时候，指定安装目录
install.packages('Seurat', repos = c('https://satijalab.r-universe.dev'), lib = "~/SeuratV4")

.libPaths(c("~/SeuratV4", .libPaths()))
.libPaths(c("~/SeuratV5", .libPaths()))
library(Seurat)
packageVersion("Seurat")

####加载R 包----
library(Seurat)
remove.packages("Seurat")
packageVersion("Seurat")
library(ggplot2)
library(patchwork)
library(data.table)
library(clusterProfiler)
#BiocManager::install("org.Mm.eg.db")
library(org.Mm.eg.db)
library(fgsea)
library(enrichplot)
library(pheatmap)
library(igraph)
library(ggraph)
#BiocManager::install("SPOTlight")
library(SPOTlight)
library(RColorBrewer)
#devtools::install_github("dmcable/spacexr", build_vignettes = FALSE)
library(spacexr)
library(doParallel)
## 如果 glmGamPoi没有安装上也没关系，不加载此安装包即可
#BiocManager::install("glmGamPoi")
library(glmGamPoi)
library(tidyverse)
library(SeuratData)
library(dplyr)
devtools::install_github("satijalab/seurat-data")

st_object<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/seurat_microbe_genus.rds")
st_object = UpdateSeuratObject(st_object)
unique(st_object$sampleid)

st_object<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/st_microbe_genus修改.rds")
st_object = UpdateSeuratObject(st_object)
unique(st_object$sampleid)

View(st_object)

metadata_group<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/group.csv")

metadata<- st_object@meta.data
rownames(st_object@meta.data)
metadata$MPR_response2<- metadata_group[match(metadata$sampleid,metadata_group$sampleid),2]

#left_join(x=metadata,y=metadata_group,by = join_by(sampleid==sampleid))

st_object<- AddMetaData(st_object,metadata = metadata)

names(st_object@meta.data)

View(st_object)

saveRDS(st_object,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/st_microbe_genus.rds")

####批量展示Vlnplot
#DefaultAssay(st_object) <- "Spatial"
#st_object[["Spatial"]]
library(patchwork)
plot1 <- VlnPlot(st_object,features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(st_object, features = "nCount_Spatial") 
wrap_plots(plot1, plot2)
Idents(st_object)<-"orig.ident"

SpatialFeaturePlot(st_object_LXD_pre, features = "nCount_Spatial") 
SpatialPlot(st_object_LXD_pre,
            ncol = 4,
            features="nCount_Spatial",
            pt.size.factor = 1.4,
            alpha = 1)

plot_list = list()
for (i in c("PDCD1", "STMN1")){
  a= SpatialFeaturePlot(st_object_LXD,
                        features = i,
                        alpha = 1,
                        ncol = 2,
                        pt.size.factor = 1.6,
                        min.cutoff = "q0",
                        max.cutoff = "q100")
  plot_list[[i]] = a
}
p1 = wrap_plots(plot_list,ncol =6)
p1

####拆分每个样本rds----
###合并之后的拆分样本保存在st_data_microbe---
sampleid<-c("LXD_ST_pre","LXD_ST_post","XL_ST_pre","XL_ST_post","XYH_ST_post",
            "XYH_ST_pre","YMS_ST_post","YMS_ST_pre","YXJ_ST_pre","LM_ST_pre")

current_microbe <- subset(st_object, subset = sampleid == "LXD_ST_pre")
saveRDS(current_microbe,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/st_data_microbe/LXD_ST_pre_microbe.rds")
LXD_ST_pre_microbe<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/st_data_microbe/LXD_ST_pre_microbe.rds")

current_microbe <- subset(st_object, subset = sampleid == "LXD_ST_post")
saveRDS(current_microbe,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/st_data_microbe/LXD_ST_post_microbe.rds")
LXD_ST_post_microbe<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/st_data_microbe/LXD_ST_post_microbe.rds")

current_microbe <- subset(st_object, subset = sampleid == "XL_ST_pre")
saveRDS(current_microbe,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/st_data_microbe/XL_ST_pre_microbe.rds")
XL_ST_pre_microbe<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/st_data_microbe/XL_ST_pre_microbe.rds")

current_microbe <- subset(st_object, subset = sampleid == "XL_ST_post")
saveRDS(current_microbe,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/st_data_microbe/XL_ST_post_microbe.rds")
XL_ST_post_microbe<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/st_data_microbe/XL_ST_post_microbe.rds")

current_microbe <- subset(st_object, subset = sampleid == "XYH_ST_pre")
saveRDS(current_microbe,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/st_data_microbe/XYH_ST_pre_microbe.rds")
XYH_ST_pre_microbe<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/st_data_microbe/XYH_ST_pre_microbe.rds")

current_microbe <- subset(st_object, subset = sampleid == "XYH_ST_post")
saveRDS(current_microbe,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/st_data_microbe/XYH_ST_post_microbe.rds")
XYH_ST_post_microbe<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/st_data_microbe/XYH_ST_post_microbe.rds")

current_microbe <- subset(st_object, subset = sampleid == "YMS_ST_post")
saveRDS(current_microbe,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/st_data_microbe/YMS_ST_post_microbe.rds")
YMS_ST_post_microbe<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/st_data_microbe/YMS_ST_post_microbe.rds")

current_microbe <- subset(st_object, subset = sampleid == "YMS_ST_pre")
saveRDS(current_microbe,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/st_data_microbe/YMS_ST_pre_microbe.rds")
YMS_ST_pre_microbe<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/st_data_microbe/YMS_ST_pre_microbe.rds")

current_microbe <- subset(st_object, subset = sampleid == "YXJ_ST_pre")
saveRDS(current_microbe,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/st_data_microbe/YXJ_ST_pre_microbe.rds")
YXJ_ST_pre_microbe<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/st_data_microbe/YXJ_ST_pre_microbe.rds")

current_microbe <- subset(st_object, subset = sampleid == "LM_ST_pre")
saveRDS(current_microbe,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/st_data_microbe/LM_ST_pre_microbe.rds")
LM_ST_pre_microbe<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/st_data_microbe/LM_ST_pre_microbe.rds")


# for (i in seq_along(sampleid)) {
#   # 动态创建子集
#   current_microbe <- subset(st_object, subset = sampleid == sampleid[i])
#   
#   # 动态保存 RDS 文件，使用 paste0() 函数拼接文件名
#   saveRDS(current_microbe, paste0("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/st_data_microbe/", sampleid[i], "_microbe.rds"))
#   
#   # 打印进度或检查信息
#   print(paste0("Processed and saved: ", sampleid[i]))
# }

####读取空转微生物rds----
####读取rds
current_microbe <- subset(st_object, subset = sampleid == "LXD_ST_post")
saveRDS(current_microbe,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/st_data_microbe/LXD_ST_post_microbe.rds")
LXD_ST_post_microbe<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/st_data_microbe/LXD_ST_post_microbe.rds")

for (i in seq_along(sampleid)) {
  # 动态创建子集
  assign(paste0(sampleid[i],"_microbe"),readRDS(paste0("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/st_data_microbe/",sampleid[i],"_microbe.rds")))
  # 读取 RDS 文件，使用 paste0() 函数拼接文件名
  print(paste0("ReadRDS_microbe: ", sampleid[i]))
}
LXD_ST_pre_microbe@meta.data$sampleid
LXD_ST_post_microbe@meta.data$orig.ident
LXD_ST_post_microbe@meta.data$sampleid
unique(XL_ST_pre_microbe$sampleid)
####分别获取微生物的metadata
for (i in seq_along(sampleid)) {
  # 获取动态创建的对象的元数据
  # 使用 get() 函数从动态变量名中获取对象
  boundary_object <- get(paste0(sampleid[i], "_microbe"))
  
  # 提取该对象的元数据并存储为新变量
  assign(paste0(sampleid[i], "_microbe_metadata"), boundary_object@meta.data)
  
  metadata <- get(paste0(sampleid[i], "_microbe_metadata"))
  
  metadata$splot_id <- paste(rownames(metadata))
  
  metadata$splot_id <- gsub("-", "_", metadata$splot_id)
  
  boundary_object@meta.data <- metadata
  assign(paste0(sampleid[i], "_microbe"),boundary_object)
  assign(paste0(sampleid[i], "_microbe_metadata"),metadata)
  # 打印进度信息
  print(paste0("metadata_microbe: ", sampleid[i]))
}

names(LXD_ST_pre_microbe@meta.data)
LXD_ST_post_microbe@meta.data$splot_id
XL_ST_pre_microbe$sampleid

####获取肿瘤边界rds----
####批量读取rds
for (i in seq_along(sampleid)) {
  # 动态创建子集
  assign(paste0(sampleid[i],"_BoundaryDefine"),readRDS(paste0("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/Cottrazm肿瘤边界识别/Cottrazm结果/",sampleid[i],"_BoundaryDefine.rds")))
  # 读取 RDS 文件，使用 paste0() 函数拼接文件名
  print(paste0("ReadRDS_BoundaryDefine: ", sampleid[i]))
}

unique(LXD_ST_pre_BoundaryDefine$orig.ident)
unique(LXD_ST_post_BoundaryDefine$orig.ident)
unique(LM_ST_pre_BoundaryDefine$orig.ident)


####分别获取metadata

for (i in seq_along(sampleid)) {
  # 获取动态创建的对象的元数据
  # 使用 get() 函数从动态变量名中获取对象
  boundary_object <- get(paste0(sampleid[i], "_BoundaryDefine"))
  
  # 提取该对象的元数据并存储为新变量
  assign(paste0(sampleid[i], "_BoundaryDefine_metadata"), boundary_object@meta.data)
  
  metadata <- get(paste0(sampleid[i], "_BoundaryDefine_metadata"))
  
  metadata$splot_id <- paste(metadata$orig.ident, rownames(metadata), sep = "_")
  metadata$splot_id <- sub("-1$", "", metadata$splot_id)
  boundary_object@meta.data <- metadata
  assign(paste0(sampleid[i], "_BoundaryDefine"),boundary_object)
  assign(paste0(sampleid[i], "_BoundaryDefine_metadata"),metadata)
  # 打印进度信息
  print(paste0("metadata_BoundaryDefine: ", sampleid[i]))
}
#names(LXD_ST_pre_BoundaryDefine@meta.data)
# names(LXD_ST_pre_BoundaryDefine_metadata)
# LXD_ST_post_BoundaryDefine_metadata$splot_id <- paste(LXD_ST_pre_BoundaryDefine_metadata$orig.ident, rownames(LXD_ST_pre_BoundaryDefine_metadata), sep = "_")


####合并微生物和肿瘤边界metadata数据---
microbe_object<-LXD_ST_pre_microbe
BoundaryDefine_object<- LXD_ST_pre_BoundaryDefine
metadata1 <- microbe_object@meta.data
metadata2 <- BoundaryDefine_object@meta.data
#interaction(metadata1$splot_id,metadata2$splot_id)
merged_metadata <- left_join(metadata1, metadata2, by = "splot_id")
merged_metadata <- merged_metadata %>%
  rename_with(~ gsub("\\.x$", "", .), ends_with(".x"))
merged_metadata <- merged_metadata %>%
  select(-ends_with(".y"))
microbe_object<- AddMetaData(microbe_object,metadata = merged_metadata)
microbe_object@meta.data

for (i in seq_along(sampleid)){
  microbe_object<- get(paste0(sampleid[i], "_microbe")) 
  BoundaryDefine_object<- get(paste0(sampleid[i], "_BoundaryDefine")) 
  metadata1 <- microbe_object@meta.data
  metadata2 <- BoundaryDefine_object@meta.data
  merged_metadata <- left_join(metadata1, metadata2, by = "splot_id")
  merged_metadata <- merged_metadata %>%
    rename_with(~ gsub("\\.x$", "", .), ends_with(".x"))
  merged_metadata <- merged_metadata %>%
    select(-ends_with(".y"))
  microbe_object<- AddMetaData(microbe_object,metadata = merged_metadata)
  saveRDS(microbe_object,paste0("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/st_data_microbe/",sampleid[i],"_microbe.rds"))
  print(paste0("saveRDS: ", sampleid[i]))
}
for (i in seq_along(sampleid)) {
  # 动态创建子集
  assign(paste0(sampleid[i],"_microbe"),readRDS(paste0("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/st_data_microbe/",sampleid[i],"_microbe.rds")))
  # 读取 RDS 文件，使用 paste0() 函数拼接文件名
  print(paste0("ReadRDS_microbe: ", sampleid[i]))
}

LXD_ST_pre_microbe
LXD_ST_pre_BoundaryDefine
names(LXD_ST_pre_microbe@meta.data)

LXD_ST_post_microbe@meta.data
LXD_ST_pre_microbe@meta.data

######纵向合并几个样本的metadata----
###拆分的样本保存
folder_path <- "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/st_data_microbe/"

# 获取文件夹中的所有 .rds 文件
file_list <- list.files(path = folder_path, pattern = "*.rds", full.names = TRUE)
metadata_list <- list()
# 加载 Seurat 对象并存储到列表中
for (i in seq_along(file_list)) {
  # 读取 .rds 文件并将其分配给一个临时变量 object
  object <- readRDS(file_list[i])
  # 提取 metadata
  metadata <- object@meta.data
  # 将 metadata 添加到列表中
  metadata_list[[i]] <- metadata
}

merged_metadata1 <- bind_rows(metadata_list)
nrow(merged_metadata1)
print(head(merged_metadata1))
merged_metadata1$CNVScores

st_object<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/seurat_microbe_genus.rds")

other_metadata <- st_object@meta.data
other_metadata$splot_id<- row.names(other_metadata)
other_metadata$splot_id <- gsub("-", "_", other_metadata$splot_id)
nrow(other_metadata)
names(other_metadata)
rownames(other_metadata)
# 根据 splot_id 进行合并，保留所有行
final_metadata <- left_join(merged_metadata1, other_metadata, by = "splot_id")
# 查看合并后的元数据
final_metadata <- final_metadata %>%
  rename_with(~ gsub("\\.x$", "", .), ends_with(".x"))
# 去除多余的 .y 列
final_metadata <- final_metadata %>%
  select(-ends_with(".y"))
# 查看处理后的元数据
nrow(final_metadata)
nrow(st_object@meta.data)
names(final_metadata)
st_object@meta.data$splot_id<- row.names(st_object@meta.data)
st_object@meta.data$splot_id <- gsub("-", "_", st_object@meta.data$splot_id)
st_object@meta.data$Mito.percent<- final_metadata[match(st_object@meta.data$splot_id,final_metadata$splot_id),79]
st_object@meta.data$Morph_snn_res.0.4<- final_metadata[match(st_object@meta.data$splot_id,final_metadata$splot_id),78]
st_object@meta.data$NormalScore<- final_metadata[match(st_object@meta.data$splot_id,final_metadata$splot_id),81]
st_object@meta.data$CNVLabel<- final_metadata[match(st_object@meta.data$splot_id,final_metadata$splot_id),82]

st_object@meta.data$cnv_score<- final_metadata[match(st_object@meta.data$splot_id,final_metadata$splot_id),83]

st_object@meta.data$CNVScores<- final_metadata[match(st_object@meta.data$splot_id,final_metadata$splot_id),84]
st_object@meta.data$Location<- final_metadata[match(st_object@meta.data$splot_id,final_metadata$splot_id),85]


print(head(final_metadata))
saveRDS(st_object,"D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/st_microbe_genus修改.rds")


####获取解卷积细胞类型rds----
######读取cell_type_new

sampleid<-c("LXD_ST_pre","LXD_ST_post","XL_ST_pre","XL_ST_post","XYH_ST_post",
            "XYH_ST_pre","YMS_ST_post","YMS_ST_pre","YXJ_ST_pre","LM_ST_pre")

for (i in seq_along(sampleid)) {
  bacteria_object <- readRDS(paste0("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/RCTD解卷积/cell_type_new/", sampleid[i], "_RCTD_object.rds"))
  # 读取 RDS 文件，使用 paste0() 函数拼接文件名
  colnames(bacteria_object@meta.data)[c(95, 96)] <- c("top1_cell_type_new", "top2_cell_typw_new")
  assign(paste0(sampleid[i], "_cell_type_new"), bacteria_object)
  print(paste0("ReadRDS_microbe: ", sampleid[i]))
}

names(LXD_ST_pre_cell_type_new@meta.data)
LXD_ST_post_cell_type_new
XL_ST_pre_cell_type_new
XL_ST_post_cell_type_new

for (i in seq_along(sampleid)) {
  # 动态创建子集
  bacteria_object <- readRDS(paste0("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/RCTD解卷积/cell_subtype/", sampleid[i], "_RCTD_object_subtype.rds"))
  # 读取 RDS 文件，使用 paste0() 函数拼接文件名
  colnames(bacteria_object@meta.data)[c(128, 129)] <- c("top1_celltype_subtype", "top2_celltype_subtype")
  assign(paste0(sampleid[i], "_subtype"), bacteria_object)
  print(paste0("ReadRDS_microbe: ", sampleid[i]))
}

names(LXD_ST_pre_subtype@meta.data)
LXD_ST_post_subtype
XL_ST_pre_subtype
XL_ST_post_subtype

for (i in seq_along(sampleid)) {
  # 动态创建子集
  bacteria_object <- readRDS(paste0("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/RCTD解卷积/cell_type_bacteria_group/", sampleid[i], "_RCTD_object_bacteria.rds"))
  # 读取 RDS 文件，使用 paste0() 函数拼接文件名
  colnames(bacteria_object@meta.data)[c(105, 106)] <- c("top1_celltype_bacteria", "top2_celltype_bacteria")
  assign(paste0(sampleid[i], "_cell_type_bacteria"), bacteria_object)
  print(paste0("ReadRDS_microbe: ", sampleid[i]))
}

names(LXD_ST_pre_cell_type_bacteria@meta.data)
names(LXD_ST_post_cell_type_bacteria@meta.data)
XL_ST_pre_cell_type_bacteria


for (i in seq_along(sampleid)){
  cell_type_new_object<- get(paste0(sampleid[i], "_cell_type_new")) 
  subtype_object<- get(paste0(sampleid[i], "_subtype")) 
  cell_type_bacteria_object<- get(paste0(sampleid[i], "_cell_type_bacteria"))
  
  metadata1 <- cell_type_new_object@meta.data
  metadata2 <- subtype_object@meta.data
  metadata3 <- cell_type_bacteria_object@meta.data
  
  merged_metadata1 <- left_join(metadata1, metadata2, by = "splot_id")
  merged_metadata1 <- merged_metadata1 %>%
    rename_with(~ gsub("\\.x$", "", .), ends_with(".x"))
  merged_metadata1 <- merged_metadata1 %>%
    select(-ends_with(".y"))
  
  merged_metadata2 <- left_join(merged_metadata1, metadata3, by = "splot_id")
  merged_metadata2 <- merged_metadata2 %>%
    rename_with(~ gsub("\\.x$", "", .), ends_with(".x"))
  merged_metadata2 <- merged_metadata2 %>%
    select(-ends_with(".y"))
  
  cell_type_new_object<- AddMetaData(cell_type_new_object,metadata = merged_metadata2)
  saveRDS(cell_type_new_object,paste0("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/RCTD解卷积/解卷积合并rds/",sampleid[i],"_RCTD合并.rds"))
  print(paste0("saveRDS: ", sampleid[i]))
  
}

for (i in seq_along(sampleid)) {
  # 动态创建子集
  assign(paste0(sampleid[i],"_RCTD合并"),readRDS(paste0("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/RCTD解卷积/解卷积合并rds/",sampleid[i],"_RCTD合并.rds")))
  # 读取 RDS 文件，使用 paste0() 函数拼接文件名
  print(paste0("ReadRDS_RCTD合并: ", sampleid[i],"_RCTD合并"))
}

LXD_ST_pre_RCTD合并@meta.data

####合并解卷积Fusobacterium_group Streptococcus_group----

sampleid<-c("LXD_ST_pre","LXD_ST_post","XL_ST_pre","XL_ST_post","XYH_ST_post",
            "XYH_ST_pre","YMS_ST_post","YMS_ST_pre","YXJ_ST_pre","LM_ST_pre")

bacteria_object <- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/RCTD解卷积/Fusobacterium_group/LXD_ST_pre_RCTD_object.rds")

names(bacteria_object@meta.data)

for (i in seq_along(sampleid)) {
  bacteria_object <- readRDS(paste0("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/RCTD解卷积/解卷积合并rds/", sampleid[i], "_RCTD合并.rds"))
  # 读取 RDS 文件，使用 paste0() 函数拼接文件名
  assign(paste0(sampleid[i], "_合并"), bacteria_object)
  print(paste0("ReadRDS_合并: ", sampleid[i]))
}
LXD_ST_pre_合并
LXD_ST_post_合并
names(LXD_ST_post_合并@meta.data)

for (i in seq_along(sampleid)) {
  bacteria_object <- readRDS(paste0("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/RCTD解卷积/Fusobacterium_group/", sampleid[i], "_RCTD_object.rds"))
  # 读取 RDS 文件，使用 paste0() 函数拼接文件名
  colnames(bacteria_object@meta.data)[c(87, 88)] <- c("top1_Fusobacterium_group", "top2_Fusobacterium_group")
  assign(paste0(sampleid[i], "_Fusobacterium_group"), bacteria_object)
  print(paste0("ReadRDS_microbe: ", sampleid[i]))
}

names(LXD_ST_pre_Fusobacterium_group@meta.data)
LXD_ST_post_Fusobacterium_group
XL_ST_pre_Fusobacterium_group
XL_ST_post_Fusobacterium_group

bacteria_object <- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/RCTD解卷积/Streptococcus_group/LXD_ST_pre_RCTD_object.rds")
names(bacteria_object@meta.data)
for (i in seq_along(sampleid)) {
  # 动态创建子集
  bacteria_object <- readRDS(paste0("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/RCTD解卷积/Streptococcus_group/", sampleid[i], "_RCTD_object.rds"))
  # 读取 RDS 文件，使用 paste0() 函数拼接文件名
  colnames(bacteria_object@meta.data)[c(87,88)] <- c("top1_Streptococcus_group", "top2_Streptococcus_group")
  assign(paste0(sampleid[i], "_Streptococcus_group"), bacteria_object)
  print(paste0("ReadRDS_microbe: ", sampleid[i]))
}

names(LXD_ST_pre_Streptococcus_group@meta.data)
LXD_ST_post_Streptococcus_group
XL_ST_pre_Streptococcus_group
XL_ST_post_Streptococcus_group
XL_ST_post_Streptococcus_group@meta.data

#######合并

for (i in seq_along(sampleid)){
  hebing_object<- get(paste0(sampleid[i], "_合并")) 
  Fusobacterium_group_object<- get(paste0(sampleid[i], "_Fusobacterium_group")) 
  Streptococcus_group_object<- get(paste0(sampleid[i], "_Streptococcus_group")) 
  
  metadata1 <- hebing_object@meta.data
  metadata2 <- Fusobacterium_group_object@meta.data
  metadata3 <- Streptococcus_group_object@meta.data
  
  merged_metadata1 <- left_join(metadata1, metadata2, by = "splot_id")
  merged_metadata1 <- merged_metadata1 %>%
    rename_with(~ gsub("\\.x$", "", .), ends_with(".x"))
  merged_metadata1 <- merged_metadata1 %>%
    select(-ends_with(".y"))
  
  merged_metadata2 <- left_join(merged_metadata1, metadata3, by = "splot_id")
  merged_metadata2 <- merged_metadata2 %>%
    rename_with(~ gsub("\\.x$", "", .), ends_with(".x"))
  merged_metadata2 <- merged_metadata2 %>%
    select(-ends_with(".y"))
  
  hebing_object<- AddMetaData( hebing_object,metadata = merged_metadata2)
  saveRDS(hebing_object,paste0("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/RCTD解卷积/解卷积合并rds2/",sampleid[i],"_RCTD合并.rds"))
  print(paste0("saveRDS: ", sampleid[i]))
  
}

for (i in seq_along(sampleid)) {
  # 动态创建子集
  assign(paste0(sampleid[i],"_RCTD合并"),readRDS(paste0("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/RCTD解卷积/解卷积合并rds2/",sampleid[i],"_RCTD合并.rds")))
  # 读取 RDS 文件，使用 paste0() 函数拼接文件名
  print(paste0("ReadRDS_RCTD合并: ", sampleid[i],"_RCTD合并"))
}

LXD_ST_pre_RCTD合并@meta.data
LXD_ST_post_RCTD合并@meta.data
XL_ST_pre_RCTD合并@meta.data


############################################################################
####rds 合并解卷积细胞类型结果----
folder_path <- "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/RCTD解卷积/解卷积合并rds2/"

# 获取文件夹中的所有 .rds 文件
file_list <- list.files(path = folder_path, pattern = "*.rds", full.names = TRUE)
metadata_list <- list()
# 加载 Seurat 对象并存储到列表中
for (i in seq_along(file_list)) {
  # 读取 .rds 文件并将其分配给一个临时变量 object
  object <- readRDS(file_list[i])
  # 提取 metadata
  metadata <- object@meta.data
  # 将 metadata 添加到列表中
  metadata_list[[i]] <- metadata
}

merged_metadata1 <- bind_rows(metadata_list)
nrow(merged_metadata1)
print(head(merged_metadata1))
merged_metadata1$top2_celltype_bacteria

st_object<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/st_microbe_genus修改.rds")

other_metadata <- st_object@meta.data
other_metadata$splot_id<- row.names(other_metadata)
other_metadata$splot_id <- gsub("-", "_", other_metadata$splot_id)
nrow(other_metadata)
names(other_metadata)
rownames(other_metadata)
# 根据 splot_id 进行合并，保留所有行
final_metadata <- left_join(merged_metadata1, other_metadata, by = "splot_id")
# 查看合并后的元数据
final_metadata <- final_metadata %>%
  rename_with(~ gsub("\\.x$", "", .), ends_with(".x"))
# 去除多余的 .y 列
final_metadata <- final_metadata %>%
  select(-ends_with(".y"))
# 查看处理后的元数据
nrow(final_metadata)
nrow(st_object@meta.data)
names(final_metadata)

st_object<- AddMetaData(st_object,metadata = final_metadata)
names(st_object@meta.data)
st_object@meta.data$top2_celltype_bacteria
saveRDS(st_object,"D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/st_microbe_genus修改.rds")


for (i in seq_along(sampleid)) {
  # 动态创建子集
  assign(paste0(sampleid[i],"_RCTD合并"),readRDS(paste0("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/RCTD解卷积/解卷积合并rds2/",sampleid[i],"_RCTD合并.rds")))
  # 读取 RDS 文件，使用 paste0() 函数拼接文件名
  print(paste0("ReadRDS_RCTD合并: ", sampleid[i],"_RCTD合并"))
}


####重新合并肿瘤边界注释----

for (i in seq_along(sampleid)) {
  # 动态创建子集
  assign(paste0(sampleid[i],"_BoundaryDefine"),readRDS(paste0("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/Cottrazm肿瘤边界识别重新注释/",sampleid[i],"_BoundaryDefine.rds")))
  # 读取 RDS 文件，使用 paste0() 函数拼接文件名
  print(paste0("ReadRDS_BoundaryDefine: ", sampleid[i]))
}

unique(LXD_ST_pre_BoundaryDefine$orig.ident)
unique(LXD_ST_post_BoundaryDefine$orig.ident)
unique(LM_ST_pre_BoundaryDefine$orig.ident)


for (i in seq_along(sampleid)) {
  # 动态创建子集
  assign(paste0(sampleid[i],"_RCTD合并"),readRDS(paste0("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/RCTD解卷积/解卷积合并rds2/",sampleid[i],"_RCTD合并.rds")))
  # 读取 RDS 文件，使用 paste0() 函数拼接文件名
  print(paste0("ReadRDS_RCTD合并: ", sampleid[i],"_RCTD合并"))
}

LXD_ST_pre_RCTD合并<-readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/RCTD解卷积/解卷积合并rds2/LXD_ST_pre_RCTD合并.rds")
LXD_ST_pre_RCTD合并@meta.data <- LXD_ST_pre_RCTD合并@meta.data[, !(colnames(LXD_ST_pre_RCTD合并@meta.data) %in% c("NormalScore", "CNVLable", "cnv_score", "CNVScores", "Location"))]
names(LXD_ST_pre_RCTD合并@meta.data )

for (i in seq_along(sampleid)){
  RCTD合并_object<- get(paste0(sampleid[i], "_RCTD合并")) 
  BoundaryDefine_object<- get(paste0(sampleid[i], "_BoundaryDefine"))
  RCTD合并_object@meta.data <- RCTD合并_object@meta.data[, !(colnames(RCTD合并_object@meta.data) %in% c("NormalScore", "CNVLabel", "cnv_score", "CNVScores", "Location"))]
  metadata1 <- RCTD合并_object@meta.data
  metadata2 <- BoundaryDefine_object@meta.data
  merged_metadata <- left_join(metadata1, metadata2, by = "splot_id")
  merged_metadata <- merged_metadata %>%
    select(-matches("\\.x$"), -matches("\\.y$"))
  RCTD合并_object<- AddMetaData(RCTD合并_object,metadata = merged_metadata)
  saveRDS(RCTD合并_object,paste0("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/RCTD解卷积/解卷积合并rds2/",sampleid[i],"_RCTD合并.rds"))
  print(paste0("saveRDS: ", sampleid[i]))
}

View(st_object)


LXD_ST_pre<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/RCTD解卷积/解卷积合并rds2/LXD_ST_pre_microbe.rds")

View(LXD_ST_pre)
####写出微生物数据----

####根据单细胞数据反卷积----
####Cell2location能更加准确的预测细胞类型的空间分布----。
library(reticulate)
py_config()
reticulate::py_install("cell2location")
#pip install scvi-tools
reticulate::py_install("scvi")
scvi <- import("scvi")
cell2location <- import("cell2location")
scanpy <- import ("scanpy")
anndata <- import ("anndata")
pandas <- import ("pandas")
numpy <- import ("numpy")
matplotlib.pyplot <- import ("matplotlib.pyplot")
matplotlib <- import ("matplotlib")
scvi <- import("scvi")
rcParams <- matplotlib$rcParams
rcParams['pdf.fonttype'] = 42 # enables correct plotting of text for PDFs
os<- import("os")

install.packages("SeuratDisk")
devtools::install_github("mojaveazure/seurat-disk")
# 加载 SeuratDisk
library(SeuratDisk)
# 将 Seurat 对象转换为 h5ad 格式（AnnData）
SaveH5Seurat(seurat_object, filename = "seurat_object.h5Seurat")
Convert("seurat_object.h5Seurat", dest = "h5ad")

####SPOTlight----
library(SPOTlight)
library(Seurat)
library(ggplot2)
library(SingleCellExperiment)
#BiocManager::install('SpatialExperiment')
#library(SpatialExperiment)
# BiocManager::install('scater')
# BiocManager::install('scran')
library(scater)
library(scran)
library(tidyverse)
library(harmony)
library(dplyr)
library(patchwork)
library(SeuratWrappers)
library(data.table)
library(clusterProfiler)
library(org.Hs.eg.db)
library(fgsea)
library(enrichplot)
library(pheatmap)
library(igraph)
library(ggraph)
library(RColorBrewer)
library(spacexr)
library(doParallel)
## 如果 glmGamPoi没有安装上也没关系，不加载此安装包即可
library(glmGamPoi)
####读取单细胞数据----
MPRNMPR_object_miMPRobe_remove<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/MPRNMPR_object_miMPRobe_remove.rds")
#MPRNMPR_object_miMPRobe_remove<- readRDS("/home/zqwangyansu/oscc_data/MPRNMPR_object_miMPRobe_remove.rds")

MPRNMPR_object_miMPRobe_remove = UpdateSeuratObject(MPRNMPR_object_miMPRobe_remove)
MPRNMPR_object_miMPRobe_remove
MPRNMPR_object_miMPRobe_remove$group_microbe2
MPRNMPR_object_miMPRobe_remove$cell_type_bacteria_group <- paste(MPRNMPR_object_miMPRobe_remove$cell_type_new, MPRNMPR_object_miMPRobe_remove$group_microbe2, sep = "_")
unique(MPRNMPR_object_miMPRobe_remove$cell_type_bacteria_group)
MPRNMPR_object_miMPRobe_remove$cell_subtype

saveRDS(MPRNMPR_object_miMPRobe_remove,"D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/MPRNMPR_object_miMPRobe_remove.rds")

scrna <- MPRNMPR_object_miMPRobe_remove
scrna$cell_type_new
View(scrna)
p1 = DimPlot(scrna, 
             group.by = "cell_type_new",
             reduction = "umap.harmony", 
             label = TRUE,
             raster=FALSE)
p1

Idents(object = scrna) <- scrna@meta.data$cell_type_new
# cluster_markers_all <-  FindAllMarkers(scrna,
#                                        assay = "RNA",
#                                        only.pos = TRUE,
#                                        test.use = "wilcox",
#                                        verbose = TRUE)

# write.csv(cluster_markers_all,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/cell_type_FindAllMarkers.csv")

write.csv(cluster_markers_all,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/CD8_FindAllMarkers.csv")

# cluster_markers_all<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/cell_type_FindAllMarkers.csv",header=T)


head(cluster_markers_all)
View(scrna)

common_genes <- intersect(rownames(scrna), rownames(st_object_LXD_pre))
scrna <- scrna[common_genes, ]
st_object_LXD_pre <- st_object_LXD_pre[common_genes, ]

sce <- as.SingleCellExperiment(scrna)

sce <- logNormCounts(sce)

### Variance modelling
# 去掉核糖体和线粒体基因
genes <- !grepl(pattern = "^RP[L|S]|MT", x = rownames(sce))
dec <- modelGeneVar(sce , subset.row = genes)

# 计算高变基因
hvg <- getTopHVGs(dec, n = 2000)

# 加上细胞注释信息
colLabels(sce) <- colData(sce)$cell_type_new
table(sce$cell_type_new)
# Compute marker genes
mgs <- scoreMarkers(sce, subset.row = genes)
# 计算标记基因，表明该基因对该细胞类型的重要性
# 保留最相关的marker基因
mgs_fil <- lapply(names(mgs), function(i) {
  x <- mgs[[i]]
  # 筛选并保留相关标记基因，即AUC>0.8的基因
  x <- x[x$mean.AUC > 0, ]
  # 将基因从高到低的权重排序
  x <- x[order(x$mean.AUC, decreasing = TRUE), ]
  # 将基因和聚类ID添加到数据框中
  x$gene <- rownames(x)
  x$cluster <- i
  data.frame(x)
})
mgs_df <- do.call(rbind, mgs_fil)
####使用lapply函数批量得到每种celltype的marker 基因，这里是根据AUC的阈值（0.8）进行筛选，其中0.8 可以根据需要自行更改。

####加载空转数据
st_object<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/st_microbe_genus.rds")
st_object_LXD_pre <- subset(st_object, subset =sampleid =="LXD_ST_pre")
#st_object<- readRDS("/home/zqwangyansu/oscc_data/st_microbe_genus.rds")
View(st_object_LXD_pre)
res <- SPOTlight(
  x = scrna,
  y = st_object_LXD_pre,
  groups = as.character(scrna$cell_type_new),
  mgs = mgs_df,
  hvg = hvg,
  weight_id = "mean.AUC",
  group_id = "cluster",
  gene_id = "gene")

#Extract data from `SPOTlight`:
str(res) #查看结果类型 
decon_mtrx <- res$mat
#这里重命名，非必要
colnames(decon_mtrx) <- paste(colnames(decon_mtrx),"Spotlight",sep = "_")

###结果得到的是每个spot的 各个celltype的占比，这里重命名是为了区分（比如后面使用其他方法进行注释）。

# 计算排名第一及第二的细胞类型
norm_weights = normalize_weights(decon_mtrx)%>% as.data.frame
maxn = function(n) function(x) order(x, decreasing = TRUE)[!is.na(x)][n]
decon_mtrx_1 = norm_weights %>% 
  dplyr::mutate(
    top1_celltype = apply(., 1, function(x) names(x)[maxn(1)(x)]),
    top2_celltype = apply(., 1, function(x) names(x)[maxn(2)(x)])
  )
head(decon_mtrx_1)


#和Seurat一致，也可以2种保存方式

#（1）添加至metadata.

#（2）虽然res 不是SeuratObject ，但是也可以构建为新的slot。

#meta 方式添加
st_object@meta.data <- cbind(st_object@meta.data, decon_mtrx)
head(st_object)
saveRDS(st_object,file="/home/zqwangyansu/oscc_data/st_object_SPOTlight.rds")
#构建新的 slot
# st_object[["SPOTlight"]] <- CreateAssayObject(t(res$mat))
# DefaultAssay(st_object) <- "SPOTlight"

# 绘制spots统计饼图
## 只保留细胞类型名
cell_types_all = colnames(decon_mtrx)[which(colnames(decon_mtrx) != "res_ss")]
head(cell_types_all)

## img_path 样本 tissue_lowres_image.png 所在
## pie_scale 表示饼图的大小
st_object$image@
st_object$cell <- names(st_object$orig.ident)
st_object <- subset(st_object, cell %in% rownames(decon_mtrx))
p = plotSpatialScatterpie(st_object, decon_mtrx, scatterpie_alpha = 1, pie_scale = 0.3)
nrow(st_object)
nrow(decon_mtrx)

p = plotSpatialScatterpie(x,
                          y,
                          cell_types = colnames(y),
                          img = FALSE,
                          slice = NULL,
                          scatterpie_alpha = 1,
                          pie_scale = 0.4,
                          degrees = NULL,
                          axis = NULL,)
p
# 保存数据
ggsave("./result/5-SPOTlight/Celltype_SPOTlight.pdf", 
       plot = p, 
       device = "pdf", 
       width = 14, 
       height = 12, 
       dpi = 300)

# 绘制 Astro 细胞类型的分布情况
p = SpatialFeaturePlot(object = Anterior1_ob,
                       features = "Astro",
                       alpha = c(0.1, 1))
p
# 保存数据
ggsave("./result/5-SPOTlight/Astro_SPOTlight.pdf",
       plot = p, 
       device = "pdf", 
       width = 14,
       height = 12,
       dpi = 300)

# 谷氨酸能神经元亚类
p = SpatialFeaturePlot(Anterior1_ob, 
                       features = c("L6b"),
                       pt.size.factor = 1.6,
                       ncol = 1)
p
# 保存数据
ggsave("./result/5-SPOTlight/L6b_SPOTlight.pdf",
       plot = p, 
       device = "pdf", 
       width = 14,
       height = 12,
       dpi = 300)

####SpatialDWLS----

####RCTD算法----
rm(list = ls())
library(Seurat)
library(ggplot2)
library(patchwork)
detach("package:data.table", unload = TRUE)
library(data.table)
library(clusterProfiler)
library(org.Mm.eg.db)
library(fgsea)
library(enrichplot)
library(pheatmap)
library(igraph)
library(ggraph)
library(SPOTlight)
library(igraph)
library(RColorBrewer)
library(spacexr)
library(doParallel)
## 如果 glmGamPoi没有安装上也没关系，不加载此安装包即可
library(glmGamPoi)
library(tidyverse)

dir.create("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/RCTD解卷积")
dir.create("./RCTD解卷积")

# 为了避免后续出现内存不足，我们这里重新读取rds
rm(list = ls())
MPRNMPR_object_miMPRobe_remove<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/MPRNMPR_object_miMPRobe_remove.rds")
#MPRNMPR_object_miMPRobe_remove<- readRDS("./MPRNMPR_object_miMPRobe_remove.rds")

MPRNMPR_object_miMPRobe_remove = UpdateSeuratObject(MPRNMPR_object_miMPRobe_remove)
MPRNMPR_object_miMPRobe_remove
scrna <- MPRNMPR_object_miMPRobe_remove
st_object<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/st_microbe_genus修改.rds")

#st_object<- readRDS("./st_microbe_genus修改.rds")

dir.create("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/RCTD解卷积/cell_type_new")
library(stringr)

sampleid<-c("LXD_ST_pre","LXD_ST_post","XL_ST_pre","XL_ST_post","XYH_ST_post",
            "XYH_ST_pre","YMS_ST_post","YMS_ST_pre","YXJ_ST_pre","LM_ST_pre")

st_subobject <- subset(st_object, subset =sampleid =="LXD_ST_post")

####RCTD cell_type_new ----
####
for (i in seq_along(sampleid)){
  ###抽取
  st_subobject <- subset(st_object, subset =sampleid ==sampleid[i])
  # 单细胞基因表达谱数据 
  counts_ref = GetAssayData(scrna, 
                            assay = "RNA",
                            layer = "counts")
  head(counts_ref)
  # 细胞类型信息
  ## 为了避免出问题，利用dplyr::mutate将细胞类型名中的"-"替换成"_"
  cell_types_ref <-scrna@meta.data %>%  dplyr::select("cell_type_new") %>%
    dplyr::mutate(cell_type2 =str_replace_all(!!!rlang::syms("cell_type_new"), "-", "_") %>%   
                    as.factor()) %>% .$cell_type2
  # 向量或列表中的元素命名
  cell_types_ref <- setNames(cell_types_ref,
                             scrna@meta.data  %>% rownames)
  head(cell_types_ref)
  
  # nUMI信息
  nUMI_ref <- scrna@meta.data %>% dplyr::select("nCount_RNA") %>% .$nCount_RNA
  nUMI_ref <- setNames(nUMI_ref,
                       scrna@meta.data %>% rownames)
  head(nUMI_ref)
  
  # 存储空间基因表达数据的参考信息
  reference = spacexr::Reference(counts_ref,
                                 cell_types_ref,
                                 nUMI_ref)
  
  str(reference)
  
  ####空转数据转换成 RCTD 对象格式
  seurat_to_SpatialRNA = function(st_ob, st_assay) {
    ## 空间坐标信息
    images = Seurat::Images(st_ob)  
    coords = do.call(rbind, lapply(images, function(slice) {
      spatial_coord = data.frame(st_ob[[slice]]@coordinates)
      colnames(spatial_coord) = c("in_tissue", "y", "x", "pxl_col_in_fullres", "pxl_row_in_fullres")
      spatial_coord$y = -1 * (spatial_coord$y) + 78
      spatial_coord = spatial_coord %>% dplyr::select(c(x, y))
      return(spatial_coord)
    }))
    ## 空间spots基因表达谱
    counts = Seurat::GetAssayData(st_ob, assay = st_assay, slot = "counts")
    ## 存储空间转录组数据的信息
    puck = spacexr::SpatialRNA(coords, counts) 
    return(puck)}
  
  spatialRNA = seurat_to_SpatialRNA(st_subobject, "Spatial")
  str(spatialRNA)
  
  #### 创建空间转录组数据的参考模板
  RCTD_ob = spacexr::create.RCTD(
    spatialRNA=spatialRNA,      ## 空间转录组
    reference=reference,        ## 参考单细胞转录组
    test_mode = FALSE,          ## 是否测试模式 
    gene_cutoff = 1.25e-4,      ## 平台效应归一化过程中：基因表达阈值
    fc_cutoff = 0.5,            ## 平台效应归一化过程中：差异倍数阈值（细胞类型间）
    gene_cutoff_reg = 2e-04,    ## 基因表达阈值
    fc_cutoff_reg = 0.75,       ## 差异倍数阈值（细胞类型间）                          
    UMI_min = 100,              ## 每个spot点的总UMI数最小阈值
    UMI_max = 2e+05,            ## 每个spot点的总UMI数最大阈值
    counts_MIN = 10,            ## 每个spot点的基因数目最小阈值   
    UMI_min_sigma = 300,        ## choose_sigma_c函数中对每个spot最小UMI数的阈值设定
    class_df = NULL,            ## 一个dataframe列表，用于对细胞类型进行分类
    CELL_MIN_INSTANCE = 1 ,     ## 每个细胞类型所需要的细胞数目
    cell_type_names = NULL,     ## 细胞类型名称列表，可以选取部分，默认使用所有细胞类型
    cell_type_profiles = NULL,  ## gene-celltype矩阵，如果提供，reference将不再使用
    MAX_MULTI_TYPES = 4,        ## multi模式时生效，用于控制每个spot中细胞数目
    keep_reference = F,         ## 是否保留所有单细胞转录组参考细胞
    CONFIDENCE_THRESHOLD = 10,  ## 细胞类型推断时置信度阈值
    DOUBLET_THRESHOLD = 25,     ## spot点为doublet非singlet时给予的罚分权重
    max_cores = 20 )            ## 并行处理，使用的线程数
  
  ######## 3.4 运行RCTD软件(full模式,运行时间较长）
  
  ## doublet_mode：RCTD方法可选 "full"，"doublet"，"multi"
  myRCTD = spacexr::run.RCTD(RCTD_ob,
                             doublet_mode = "full")
  head(myRCTD@results$weights)
  
  save(myRCTD, file = paste0("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/RCTD解卷积/cell_type_new/",sampleid[i],"_myRCTD.RData"))
  
  ####计算排名第一及第二的细胞类型
  # normalize_weights函数对权重进行归一化，使得每个spots中的各种细胞类型比例加起来等于1，并获取占比排名第一及第二的细胞类型，并将其存入seurat对象中，便于后续使用
  norm_weights = normalize_weights(myRCTD@results$weights)%>% as.data.frame
  st_subobject@misc$RCTD_result = norm_weights
  maxn = function(n) function(x) order(x, decreasing = TRUE)[!is.na(x)][n]
  metadata = norm_weights %>% 
    dplyr::mutate(
      top1_celltype = apply(., 1, function(x) names(x)[maxn(1)(x)]),
      top2_celltype = apply(., 1, function(x) names(x)[maxn(2)(x)])
    )
  st_subobject = Seurat::AddMetaData(st_subobject, metadata = metadata)
  head(st_subobject@meta.data)
  # 保存
  saveRDS(st_subobject,file=paste0("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/RCTD解卷积/cell_type_new/",sampleid[i],"_RCTD_object.rds"))
  print(paste0("saveRDS: ", sampleid[i]))
}


for (i in seq_along(sampleid)) {
  # 动态创建子集
  assign(paste0(sampleid[i],"_RCTD_object"),readRDS(paste0("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/Cottrazm肿瘤边界识别/Cottrazm结果/",sampleid[i],"_RCTD_object.rds")))
  # 读取 RDS 文件，使用 paste0() 函数拼接文件名
  paste0(sampleid[i],"_RCTD_object")
  print(paste0("ReadRDS_RCTD_object: ", sampleid[i]))
}

####验证跑出的rds 正不正确

LXD_pre<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/RCTD解卷积/cell_type_new/LXD_ST_pre_RCTD_object.rds")

LXD_pre

LXD_post<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/RCTD解卷积/cell_type_new/LXD_ST_post_RCTD_object.rds")
LXD_post

XL_ST_pre<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/RCTD解卷积/cell_type_new/XL_ST_pre_RCTD_object.rds")
XL_ST_pre

XL_ST_post<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/RCTD解卷积/cell_type_new/XL_ST_post_RCTD_object.rds")
XL_ST_post

XYH_ST_post<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/RCTD解卷积/cell_type_new/XYH_ST_post_RCTD_object.rds")
XYH_ST_post

XYH_ST_pre<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/RCTD解卷积/cell_type_new/XYH_ST_pre_RCTD_object.rds")
XYH_ST_pre

YMS_ST_post<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/RCTD解卷积/cell_type_new/YMS_ST_post_RCTD_object.rds")
YMS_ST_post

YXJ_ST_pre<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/RCTD解卷积/cell_type_new/YXJ_ST_pre_RCTD_object.rds")
YXJ_ST_pre

LM_ST_pre<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/RCTD解卷积/cell_type_new/LM_ST_pre_RCTD_object.rds")
LM_ST_pre

sampleid <- c("LXD_ST_pre","LXD_ST_post","XL_ST_pre","XL_ST_post","XYH_ST_post",
              "XYH_ST_pre","YMS_ST_post","YMS_ST_pre","YXJ_ST_pre","LM_ST_pre")

LM_ST_pre$top1_celltype

####RCTD cell_subtype ----

LXD_pre<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/RCTD解卷积/cell_subtype/LXD_ST_pre_RCTD_object_subtype.rds")

LXD_pre

LXD_post<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/RCTD解卷积/cell_subtype/LXD_ST_post_RCTD_object_subtype.rds")
LXD_post

XL_ST_pre<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/RCTD解卷积/cell_subtype/XL_ST_pre_RCTD_object_subtype.rds")
XL_ST_pre

XL_ST_post<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/RCTD解卷积/cell_subtype/XL_ST_post_RCTD_object_subtype.rds")
XL_ST_post

XYH_ST_post<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/RCTD解卷积/cell_subtype/XYH_ST_post_RCTD_object_subtype.rds")
XYH_ST_post

XYH_ST_pre<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/RCTD解卷积/cell_subtype/XYH_ST_pre_RCTD_object_subtype.rds")
XYH_ST_pre

YMS_ST_post<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/RCTD解卷积/cell_subtype/YMS_ST_post_RCTD_object_subtype.rds")
YMS_ST_post

YXJ_ST_pre<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/RCTD解卷积/cell_subtype/YXJ_ST_pre_RCTD_object_subtype.rds")
YXJ_ST_pre

LM_ST_pre<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/RCTD解卷积/cell_subtype/LM_ST_pre_RCTD_object_subtype.rds")
LM_ST_pre

####数据可视化
# 绘制饼图
mat=normalize_weights(myRCTD@results$weights)%>%as.data.frame()
?plotSpatialScatterpie
?spatial_scatterpie
p = SPOTlight::spatial_scatterpie(se_obj = st_object_LXD_pre,
                                  cell_types_all = colnames(mat),
                                  img_path = "./RCTD_st_object_LXD_pre/RCTD_st_object_LXD_pre_tissue_lowres_image.png",
                                  pie_scale = 0.4)
p
# 保存数据
ggsave("Celltype_RCTD.pdf", 
       plot = p,
       device = "pdf", 
       width = 14, 
       height = 12,
       dpi = 300)

# 绘制 Astro 细胞类型的分布情况

p = SpatialFeaturePlot(object = st_object_LXD_pre,
                       features = "cell_type_new",
                       alpha = c(0.1, 1))
p
#保存数据
ggsave("Astro_RCTD.pdf", 
       plot = p,
       device = "pdf", 
       width = 14, 
       height = 12,
       dpi = 300)

# 谷氨酸能神经元亚类
p = SpatialFeaturePlot(st_object_LXD_pre, 
                       features = c("L6b"),
                       pt.size.factor = 1.6,
                       ncol = 1)
p
#保存数据
ggsave("L6b_RCTD.pdf",
       plot = p, 
       device = "pdf", 
       width = 14,
       height = 12,
       dpi = 300)

# top1_celltype
p = SpatialDimPlot(Anterior1_ob,group.by = "top1_celltype")+
  theme( legend.key = element_rect(fill = NA, color = NA),
         legend.text = element_text(size = 10),
         plot.title = element_text(face = "bold", hjust = 0.5))

p
#保存数据
ggsave("top1_celltype_RCTD.pdf",
       plot = p, 
       device = "pdf", 
       width = 14,
       height = 12,
       dpi = 300)

####CRAD算法----
devtools::install_github('xuranw/MuSiC') #安装依赖包
devtools::install_github('YingMa0107/CARD') #安装CARD反卷积R包
#BiocManager::install("TOAST")

library(CARD)
library(MuSiC)
library(Seurat)
library(patchwork)
library(tidyverse)

####载入空转数据
st_object<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/st_microbe_genus.rds")
st_object_LXD_pre <- subset(st_object, subset =sampleid =="LXD_ST_pre")

###获取空转的counts表达矩阵
spatial_count <-  st_object_LXD_pre@assays$Spatial@counts
nrow(spatial_count)
ncol(spatial_count)
spatial_count[1:4,1:4]
######获取LXD_ST_pre空转的空间位置矩阵----
spatial_loca <- st_object@images$LXD_ST_pre@coordinates
nrow(spatial_loca)
ncol(spatial_loca)
spatial_location <- spatial_loca[,2:3]
#名字必须是x y ，否则后面CARD_deconvolution会报错
colnames(spatial_location) <- c("x","y")
spatial_location[1:3,]

#                   x   y
#AAACAAGTATCTCCCA-1 50 102
#AAACACCAATAACTGC-1 59  19
#AAACAGAGCGACTCCT-1 14  94

#载入‘单细胞多样本整合分析’教程 保存的 GBM-scRNA 单细胞数据
MPRNMPR_object_miMPRobe_remove<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/MPRNMPR_object_miMPRobe_remove.rds")
MPRNMPR_object_miMPRobe_remove = UpdateSeuratObject(MPRNMPR_object_miMPRobe_remove)
MPRNMPR_object_miMPRobe_remove
scRNA <- MPRNMPR_object_miMPRobe_remove

#获取单细胞counts矩阵
sc_count <- scRNA@assays$RNA@counts

#获取单细胞细胞注释矩阵
sc_meta <- scRNA@meta.data %>% 
  dplyr::select(cell_id,orig.ident,cell_type_new) %>% as.data.frame()

head(sc_meta)
#                                                 cellID orig.ident   celltype
#GSM4119531_AAACCTGAGTCAAGGC GSM4119531_AAACCTGAGTCAAGGC GSM4119531   MES like
#GSM4119531_AAACCTGTCAGGCAAG GSM4119531_AAACCTGTCAGGCAAG GSM4119531 Macrophage
#GSM4119531_AAACCTGTCCTGCCAT GSM4119531_AAACCTGTCCTGCCAT GSM4119531   OPC like


##构建CARD对象，并进行空间细胞成分反卷积
CARD_obj = createCARDObject( 
  sc_count = sc_count, 
  sc_meta = sc_meta, 
  spatial_count = spatial_count, 
  spatial_location = spatial_location, 
  ct.varname = "cell_type_new", 
  ct.select = unique(sc_meta$cell_type_new), #细胞类型列名
  sample.varname = "orig.ident")
## QC on scRNASeq dataset! ...
## QC on spatially-resolved dataset! ...

#CARD 解卷积
CARD_obj = CARD_deconvolution(CARD_object = CARD_obj)

## create reference matrix from scRNASeq...
## Select Informative Genes! ...
## Deconvolution Starts! ...
## Deconvolution Finish! ...
##以上步骤就完成了空间细胞成分反卷积

#下面对反卷积结果进行可视化
#CARD-spot 可视化spot的细胞类型分布饼图

######修改颜色顺序与图1对应----
cell_type_color1 <- c("#ffc556",#B cells
                      "#00FFFF",#"Endothelial cells"
                      "#20B2AA",#"Epithelial cells"
                      "#32cd32",#"Fibroblasts"
                      "#9370DB",#'Mast cells
                      "#0000FF",#'Myeloid cells'
                      "#FF00ff",#"Pericytes"
                      "#1E90FF",#'Plasma cells'
                      "#006400",#"Smooth muscle cells"
                      "#FF6347")#T cells)

cell_type_color1 <- c("grey90",#B cells
                      "grey90",#"Endothelial cells"
                      "grey90",#"Epithelial cells"
                      "grey90",#"Fibroblasts"
                      "grey90",#'Mast cells
                      "#FF6347",#'Myeloid cells'
                      "grey90",#"Pericytes"
                      "grey90",#'Plasma cells'
                      "grey90",#"Smooth muscle cells"
                      "grey90")#T cells)

head(CARD_obj@Proportion_CARD)
head(CARD_obj@spatial_location)
p1<-CARD.visualize.pie(proportion = CARD_obj@Proportion_CARD,
                       spatial_location = CARD_obj@spatial_location, 
                       colors = cell_type_color1)

p2<-SpatialPlot(st_object_LXD_pre,group.by = 'microbe_detected_status',cols = c('Microbe detected'='#007799','No microbe detected'='#AA0000'))
p1+p2
unique(st_object_LXD_pre$microbe_detected_status)

#左图反卷积结果显示，Tumor区域主要是OPC like细胞

#选择一些感兴趣的细胞类型进行可视化，查看细胞类型比例的空间分布
#选择一些感兴趣的细胞类型分别进行可视化
ct.visualize = c("OPC like","AC like","MES like")

p3 <- CARD.visualize.prop(proportion = CARD_obj@Proportion_CARD,        
                          spatial_location = CARD_obj@spatial_location, 
                          ct.visualize = ct.visualize,                
                          colors = c("lightblue","lightyellow","red"), 
                          NumCols = 3,pointSize = 1)#图中spot大小
p3

#图片向右旋转了90°，可以自行对CARD.visualize.prop函数进行更改。
#原始的基因表达可视化
p4 <- CARD.visualize.gene(
  spatial_expression = CARD_obj@spatial_countMat,
  spatial_location = CARD_obj@spatial_location,
  gene.visualize = c("PDCD1","SYT1","JUNB"),
  colors = NULL,
  NumCols =3)
p4

#同时可视化两种细胞类型
p5 = CARD.visualize.prop.2CT(
  proportion = CARD_obj@Proportion_CARD,                             ### Cell type proportion estimated by CARD
  spatial_location = CARD_obj@spatial_location,                  ### two cell types you want to visualize
  ct2.visualize = c("Myeloid cells","B cells"),
  colors = list(c("lightblue","lightyellow","red"),c("lightblue","lightyellow","black")))       ### two color scales    

p6 = CARD.visualize.prop.2CT(
  proportion = CARD_obj@Proportion_CARD,                             ### Cell type proportion estimated by CARD
  spatial_location = CARD_obj@spatial_location,                  ### two cell types you want to visualize
  ct2.visualize = c("MES like","Oligo"),
  colors = list(c("lightblue","lightyellow","red"),c("lightblue","lightyellow","black")))       ### two color scales    

p5+p6

#细胞类型比例相关图热图
p6 <- CARD.visualize.Cor(CARD_obj@Proportion_CARD,colors = NULL) 
p6

#保存反卷积结果
save(CARD_obj,file = 'CARD_obj.rdata')



####Cottrazm 识别肿瘤边界----

library(RColorBrewer)
library(patchwork)
library(ggtree)
library(BiocGenerics)
library(readr)
#BiocManager::install("rtracklayer")
library(rtracklayer)
library(infercnv)
library(phylogram)
library(utils)
library(dendextend)
library(assertthat)
library(reticulate)
library(openxlsx)
library(scatterpie)
library(cowplot)
library(stats)
library(quadprog)
library(data.table)
library(Rfast)
library(ggrepel)
library(tibble)
library(clusterProfiler)
library(utils)
library(org.Hs.eg.db)
devtools::install_local("C:/Users/Wangy/Downloads/Cottrazm-main.zip")
library(CellChat)
remotes::install_github("sqjin/CellChat")
source("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/肿瘤边界cottrazm_rewrite_functions.R")
library(Cottrazm)

####合并
LM_ST_pre<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/Cottrazm肿瘤边界识别/Cottrazm结果/LM_ST_pre/LM_ST_pre_BoundaryDefine.rds/LM_ST_pre_BoundaryDefine.rds")
View(LM_ST_pre)
LM_ST_pre$Location
library(Seurat)
library(ggplot2)
SpatialDimPlot(object = LM_ST_pre, group.by = "CNVLabel", pt.size = 2 ) +
  ggtitle("Spatial Distribution by Location")

SpatialFeaturePlot(object = LM_ST_pre, features = c("PDCD1"), ncol = 1) + 
  ggtitle("Spatial Expression of Multiple Genes")

SpatialFeaturePlot(object = LM_ST_pre, features = "PDCD1", alpha = c(0.1, 1),pt.size.factor = 4) +
  ggtitle("Spatial Expression with Adjusted Transparency")

VlnPlot(object = LM_ST_pre, features = "GAPDH")

SpatialPlot(LM_ST_pre, features = "MS4A1",pt.size = 2)
SpatialFeaturePlot(LM_ST_pre, features = "MS4A1",pt.size = 2)


####空转微生物组成分析----
.libPaths(c("~/SeuratV4", .libPaths()))
library(Seurat)
packageVersion("Seurat")
st_object<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/st_microbe_genus修改.rds")
st_metadata$species
group<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/group.csv")
View(group)
library(dplyr)
metadata<- FetchData(st_object,"sampleid")
metadata$cell_id <- rownames(metadata)
metadata<- left_join(x=metadata,y=group,by = join_by(sampleid==sampleid))
rownames(metadata)<-metadata$cell_id
View(metadata)
st_object<- AddMetaData(st_object,metadata = metadata)
st_metadata<- st_object@meta.data

write.csv(st_metadata,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/空转微生物SAHMI/st_metadata.csv")

library(microeco)
library(ggplot2)
library(dplyr) 

######MPR----
#genus特征表
genus<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/空转微生物SAHMI/st_metadata.csv",header = T)
names(genus)
genus$group4<- paste0(genus$splot_id,"-",genus$MPR_Response2)
genus$group4 <- genus$sampleid
genus0<- genus[,c(88,28:72)] %>% as.data.frame()
names(genus0)
genus1 <- genus0 %>%
  group_by(!!sym("group4")) %>%  # 将 V53 替换为你的实际列名
  summarise(across(1:45, sum, na.rm = TRUE))  %>% as.data.frame()
names(genus1)
genus1$group4

genus2<- t(genus1)
rownames(genus2)
View(genus2)
colnames(genus2)<- genus2[1,]
genus2<- genus2[-1,] %>% as.data.frame
View(genus2)
genus2[] <- lapply(genus2, function(x) as.numeric(as.character(x)))

###分组信息
sample <- genus[,c(19,85,87,88)] %>% as.data.frame()
View(sample)
otu_cols <- colnames(genus2)
# 去除重复的 sampleid，并保持顺序
filtered_df <- sample %>%
  distinct(group4, .keep_all = TRUE) %>%
  arrange(match(group4, otu_cols))
# 确保过滤后的数据框与 otu_table 的列名顺序一致
final_df <- filtered_df[filtered_df$group4 %in% otu_cols, ]
# 如果需要调整列名顺序
final_df <- final_df[match(otu_cols, final_df$group4), ]
sample1<- final_df
View(sample1)
rownames(sample1)<- sample1$group4
#物种分类表
tax<- rownames(genus2) %>% as.data.frame()
View(tax)
colnames(tax)<- "genus"
rownames(tax) <- tax$genus

df<-microtable$new(sample_table=sample1,
                     otu_table=genus2,
                     tax_table=tax,
                     auto_tidy=F)
df
colnames(genus2)
rownames(sample1)
sample1$Location
sample$MPR_Response2
genus_df<-trans_abund$new(dataset=df,taxrank="genus",ntaxa=30,groupmean = "MPR_Response2")
genus_df$data_abund$Sample<- factor(genus_df$data_abund$Sample,levels=c("Pre_NMPR","Pre_MPR","Post_MPR"))
p1<-genus_df$plot_bar(others_color="grey70",#剩余分类的填充色
                #facet="MPR_Response2",#根据组进行分面
                xtext_keep=T,#是否显示样本名称
                legend_text_italic=F)+
  theme(
    axis.text.x = element_text(size = 12,face ="bold"),
    axis.text.y = element_text(size = 12), 
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 13, face = "bold") 
  )
 
p1
ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/空转微生物SAHMI/空转微生物_MPR_Response2.pdf",p1,width = 9,height=8)

######箱状图----
?trans_abund
#input_taxaname
#color_values
#order_x
#x_axis_name
genus_df$data_abund$MPR_Response2
genus_df<-trans_abund$new(dataset=df,taxrank="genus",input_taxaname=c("Streptococcus","Fusobacterium"))
genus_df$data_abund$MPR_Response2<- factor(genus_df$data_abund$MPR_Response2,levels=c("Pre_NMPR","Pre_MPR","Post_MPR"))
p3<- genus_df$plot_box(color_values=c("#4974a4","#4dae47","#f29600"),group="MPR_Response2", 
                  position = position_dodge(0.9),
                  xtext_angle=30,
                  show_point=FALSE,
                  point_size = 6)+
  scale_y_continuous(limits = c(0, 3))+
  theme_minimal()+
  theme(
    panel.background = element_blank(),  # 去除面板背景
    plot.background = element_blank(),    # 去除整个图形的背景
    panel.grid.major = element_blank(),   # 去除主要网格线
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 14,angle=0,vjust = 0.5,hjust=0.5),
    axis.text.y = element_text(size = 12), 
    axis.title.y = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 13, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA,linewidth = 1)  # 添加边框
  )
p3

library(ggplot2)
library(ggpubr)
library(rstatix)

df_p_val1 <- genus_df1 %>%group_by(Taxonomy) %>%
  wilcox_test(Abundance~MPR_response2) %>%
  adjust_pvalue(p.col="p",method="bonferroni") %>%
  add_significance(p.col="p.adj") %>%
  add_xy_position(x="MPR_response2",dodge=0.8)

stat_pvalue_manual(df_p_val1,label="p.adj.signif",hide.ns=T,
                   tip.length = 0,label.size = 5,color="black")

names(genus_df$data_abund)
genus_df1<- genus_df$data_abund %>% filter(Taxonomy %in% c("Streptococcus", "Fusobacterium"))
head(genus_df1)
p3<- ggplot(genus_df1, aes(x =Taxonomy, y = Abundance,fill=MPR_Response2)) +
  geom_boxplot(outlier.shape = NA) +  
  labs(title = "", x = "", y = "Relative abundance(%)") +  # 添加标题和轴标签
  theme_minimal() +  # 使用简单主题
  scale_fill_manual(values = c("#c62d17", "#04d486","#f6c619"))+
  theme_minimal()+
  theme(
    panel.background = element_blank(),  # 去除面板背景
    plot.background = element_blank(),    # 去除整个图形的背景
    panel.grid.major = element_blank(),   # 去除主要网格线
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 14,angle=0,vjust = 0.5,hjust=0.5),
    axis.text.y = element_text(size = 12), 
    axis.title.y = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 13, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA,linewidth = 1),  # 添加边框
  )

p3

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/空转微生物SAHMI/空转微生物_MPR_boxplot.pdf",p3,width = 6,height=5)

library(ggradar)
devtools::install_github("ricardo-bion/ggradar")

genus_df<-trans_abund$new(dataset=df,taxrank="genus",input_taxaname=c("Streptococcus","Fusobacterium"),groupmean = "MPR_Response2")
genus_df<-trans_abund$new(dataset=df,taxrank="genus",ntaxa=10,groupmean = "MPR_Response2")
genus_df<-trans_abund$new(dataset=df,taxrank="genus",ntaxa=10)
genus_df$data_abund$Sample<- factor(genus_df$data_abund$Sample,levels=c("Pre_NMPR","Pre_MPR","Post_MPR"))
genus_df$data_abund$Sample<- factor(genus_df$data_abund$Sample,levels=c("YXJ_ST_pre","LM_ST_pre",
                                                                        "LXD_ST_pre","XL_ST_pre","XYH_ST_pre","YMS_ST_pre",
                                                                        "LXD_ST_post","XL_ST_post","XYH_ST_post","YMS_ST_post"))
genus_df$data_abund$MPR_Response2<- factor(genus_df$data_abund$MPR_Response2,levels=c("Pre_NMPR","Pre_MPR","Post_MPR"))
p4<- genus_df$plot_heatmap(
  color_values = rev(RColorBrewer::brewer.pal(n = 11, name = "RdYlBu")),
  #facet ="MPR_Response2",
  x_axis_name = NULL,
  order_x = NULL,
  xtext_size = 14,
  ytext_size = 14,
  xtext_angle = 30)
p4
ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/空转微生物SAHMI/空转微生物_MPR_heatmap_top30.pdf",p4,width =8,height=7)

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/空转微生物SAHMI/空转微生物_MPR_heatmap_top10每个患者.pdf",p4,width =8,height=8)


######location----
genus<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/空转微生物SAHMI/st_metadata.csv",header = T)
names(genus)
head(genus)
genus$group3<- paste0(genus$sampleid,"-",genus$Location)
genus$group3<- paste0(genus$MPR_Response2,"-",genus$Location)
genus0<- genus[,c(88,28:72)] %>% as.data.frame()
names(genus0)

genus1 <- genus0 %>%
  group_by(!!sym("group3")) %>%  # 将 V53 替换为你的实际列名
  summarise(across(1:45, sum, na.rm = TRUE))  %>% as.data.frame()
names(genus1)
genus1$group3
genus2<- t(genus1)
rownames(genus2)
View(genus2)
colnames(genus2)<- genus2[1,]
genus2<- genus2[-1,] %>% as.data.frame
View(genus2)
genus2[] <- lapply(genus2, function(x) as.numeric(as.character(x)))
colnames(genus2)
genus2<- genus2[,-27]
###分组信息
sample <- genus[,c(19,85,87,88)] %>% as.data.frame()
View(sample)
otu_cols <- colnames(genus2)
# 去除重复的 sampleid，并保持顺序
filtered_df <- sample %>%
  distinct(group3, .keep_all = TRUE) %>%
  arrange(match(group3, otu_cols))
# 确保过滤后的数据框与 otu_table 的列名顺序一致
final_df <- filtered_df[filtered_df$group3 %in% otu_cols, ]
# 如果需要调整列名顺序
final_df <- final_df[match(otu_cols, final_df$group3), ]
sample1<- final_df
View(sample1)
rownames(sample1)<- sample1$group3
#物种分类表
tax<- rownames(genus2) %>% as.data.frame()
View(tax)
colnames(tax)<- "genus"
rownames(tax) <- tax$genus

df<-microtable$new(sample_table=sample1,
                   otu_table=genus2,
                   tax_table=tax,
                   auto_tidy=F)
df
colnames(genus2)
rownames(sample1)
sample1$Location
sample$MPR_Response2

location_df<-trans_abund$new(dataset=df,taxrank="genus",ntaxa=10,groupmean = "Location")
location_df$data_abund$Sample<- factor(location_df$data_abund$Sample,levels=c("Mal","Bdy","nMal"))
p2<-location_df$plot_bar(others_color="grey70",#剩余分类的填充色
                         #facet="MPR_Response2",#根据组进行分面
                         xtext_keep=T,#是否显示样本名称
                         legend_text_italic=F)+
  theme(
    axis.text.x = element_text(size = 14,face ="bold"),
    axis.text.y = element_text(size = 12), 
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 13, face = "bold") 
  )

p2

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/空转微生物SAHMI/空转微生物_恶行细胞非恶性细胞.pdf",p2,width = 9,height=8)

location_df$data_abund$Location
location_df<-trans_abund$new(dataset=df,taxrank="genus",input_taxaname=c("Streptococcus","Fusobacterium"))
location_df$data_abund$Location<- factor(location_df$data_abund$Location,levels=c("Mal","Bdy","nMal"))
write.csv(location_df$data_abund,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/空转微生物SAHMI/恶性非恶性boxplot数据.csv")
data<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/空转微生物SAHMI/恶性非恶性boxplot数据.csv",header = T,row.names = 1)
p5<- location_df$plot_box(color_values=c("#c62d17", "#04d486","#f6c619"),group="Location", 
                       position = position_dodge(0.9),
                       xtext_angle=30,
                       point_size = 6)+
  scale_y_continuous(limits = c(0, 3))+
  theme_minimal()+
  theme(
    panel.background = element_blank(),  # 去除面板背景
    plot.background = element_blank(),    # 去除整个图形的背景
    panel.grid.major = element_blank(),   # 去除主要网格线
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 14,angle=0,vjust = 0.5,hjust=0.5),
    axis.text.y = element_text(size = 12), 
    axis.title.y = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 13, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA,linewidth = 1),  # 添加边框
  )
p5
names(data)
data1 <- location_df$data_abund %>% filter(Taxonomy %in% c("Streptococcus", "Fusobacterium"))
data1$Taxonomy<- factor(data1$Taxonomy,levels=c("Streptococcus","Fusobacterium"))
p5<- ggplot(data1, aes(x =Taxonomy, y = Abundance,fill=Location)) +
  geom_boxplot(outlier.shape = NA) + 
  labs(title = "", x = "", y = "Relative abundance(%)") +  # 添加标题和轴标签
  theme_minimal() +  # 使用简单主题
  scale_fill_manual(values = c("#c62d17", "#04d486","#f6c619"))+
  scale_y_continuous(limits = c(0, 6))+
  theme_minimal()+
  theme(
    panel.background = element_blank(),  # 去除面板背景
    plot.background = element_blank(),    # 去除整个图形的背景
    panel.grid.major = element_blank(),   # 去除主要网格线
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 14,angle=0,vjust = 0.5,hjust=0.5),
    axis.text.y = element_text(size = 12), 
    axis.title.y = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 13, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA,linewidth = 1),  # 添加边框
  )

p5
ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/空转微生物SAHMI/空转微生物_Mal_boxplot.pdf",p5,width = 6,height=5)

library(ggradar)
devtools::install_github("ricardo-bion/ggradar")

location_df<-trans_abund$new(dataset=df,taxrank="genus",input_taxaname=c("Streptococcus","Fusobacterium"),groupmean = "Location")
location_df<-trans_abund$new(dataset=df,taxrank="genus",ntaxa=11,groupmean = "Location")
location_df$data_abund$Sample<- factor(location_df$data_abund$Sample,levels=c("Mal","Bdy","nMal"))
location_df<-trans_abund$new(dataset=df,taxrank="genus",ntaxa=11)
location_df$data_abund$MPR_Response2<- factor(location_df$data_abund$MPR_Response2,levels=c("Pre_NMPR","Pre_MPR","Post_MPR"))
View(location_df$data_abund)
p6<- location_df$plot_heatmap(
  color_values = rev(RColorBrewer::brewer.pal(n = 11, name = "RdYlBu")),
  facet ="MPR_Response2",
  x_axis_name = NULL,
  order_x = NULL,
  xtext_size = 10,
  ytext_size = 14,
  xtext_angle = 30)
p6
ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/空转微生物SAHMI/空转微生物_Mal_heatmap_top10.pdf",p6,width =8,height=8)


####自己做微生物图----
abundance<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/空转微生物SAHMI/st_metadata.csv",header = T,row.names = 1)
View(abundance)
colnames(abundance)

abundance_mat <- abundance[,c(27:71)]
View(abundance_mat)
abundance_mat$group <- paste0(sub("(-.*)", "", rownames(abundance_mat)),"-",genus$Location)
# 对每个分组进行按列求和
result <- abundance_mat %>%
  group_by(group) %>%
  summarise(across(everything(), sum, .names = "sum_{col}"))
abundance_mat2<- result

tax_name <- colnames(abundance_mat) %>% as.data.frame()
View(tax_name)
names(tax_name)<- "genus"
head(abundance_mat)
head(tax_name)

group <- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/group.csv", header=T)
names(group)
######genusheatmap----
library (pheatmap)
abundance_mat3 <- t(abundance_mat2) %>% as.data.frame()
View(abundance_mat3)
colnames(abundance_mat3) <-abundance_mat3[1,]
abundance_mat3<- abundance_mat3[-1,]
rownames(abundance_mat3)
rownames(abundance_mat3) <- gsub("^sum_", "", rownames(abundance_mat3))

View(abundance_mat3)
df3<- abundance_mat3
df3_numeric <- apply(abundance_mat3, 2, as.numeric)
df3_numeric <- apply(df3_numeric,2,function(x) x/sum(x)) %>% as.data.frame()
colnames(df3_numeric)
df3_numeric<-df3_numeric[,-27]
View(df3_numeric)
# 计算行和
df3_numeric$rowsum <- rowSums(df3_numeric)
rownames(df3_numeric)<- rownames(abundance_mat3)

View(df3_numeric)
df4 <- df3_numeric[order (df3_numeric$rowsum,decreasing=TRUE),]
colnames(df4)
df5 = df4[,-31]#删除求和列
#求物种相对丰度
View(df5)

#取前10行
df7 <-  df5[1:11,]
df8 <- 1-apply(df7, 2, sum) #计算剩下物种的总丰度
#合并数据
df9 <- rbind(df7,df8)
row.names(df9)[12]="Others"
View(df9)
#导出数据
write.table (df9, file ="genus_x.csv",sep =",", quote =FALSE)

bdy_avg <- df9 %>%
  select(contains("Bdy")) %>%
  rowMeans() 

Mal_avg <- df9 %>%
  select(contains("Mal")) %>%
  rowMeans() 

nMal_avg <- df9 %>%
  select(contains("nMal")) %>%
  rowMeans()

df10<- scale(cbind(bdy_avg,Mal_avg,nMal_avg) )%>% as.data.frame()

# 行列注释信息
group1<- gsub(".*-", "", colnames(abundance_mat3))[-27]
group1<- colnames(df10)
annotation_col = data.frame(group1)#行注释矩阵
rownames(annotation_col) <- annotation_col$group1
rownames(annotation_col) = colnames(df10)
#annotation_row = data.frame( g2 = factor(rep(c("c", "d", "e"), c(3, 4, 4))))#列注释矩阵
#rownames(annotation_row) = rownames(df9)
#分组颜色
colors = list(group1 = c(Bdy = "#1B9E77", nMal = "#D95F02",Mal = "red"))
              #g2 = c( c= "blue", d = "#E7298A", e = "#66A61E"))
#绘图
pheatmap(df9, 
         display_numbers = matrix(ifelse(df9 > 0.5, "*", ""), nrow(df9)), #设置条件并进行展示
         annotation_col = annotation_col,#加入行列注释信息
         annotation_colors = colors, #行列注释颜色设置
         angle_col = "45", #行标签倾斜角度
         cellwidth=25, cellheight=15, #格子大小
         cluster_rows=F, treeheight_col = 15, #树的高度
         legend_breaks=c(0,0.2,0.4,0.6,0.8) #图注标签展示
)

View(df9)

####空转微生物图片展示----
####读取空转数据
sampleid <- c("LXD_ST_pre","LXD_ST_post","XL_ST_pre","XL_ST_post","XYH_ST_post",
              "XYH_ST_pre","YMS_ST_post","YMS_ST_pre","YXJ_ST_pre","LM_ST_pre")

st_object<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/st_microbe_genus修改.rds")

for (i in seq_along(sampleid)) {
  # 动态创建子集
  assign(paste0(sampleid[i],"_microbe"),readRDS(paste0("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/st_data_microbe/",sampleid[i],"_microbe.rds")))
  # 读取 RDS 文件，使用 paste0() 函数拼接文件名
  print(paste0("ReadRDS_microbe: ", sampleid[i],"_microbe"))
}


library(ggplot2)
library(randomcoloR)
colors<- distinctColorPalette(24)
#c("#c62d17", "#04d486", "#f6c619")

YXJ_ST_pre_microbe$Streptococcus_group <- ifelse(YXJ_ST_pre_microbe$Streptococcus > 0, 
                                                 "Streptococcus+", 
                                                 "Streptococcus-")
YXJ_ST_pre_microbe$Fusobacterium_group <- ifelse(YXJ_ST_pre_microbe$Fusobacterium > 0, 
                                                 "Fusobacterium+", 
                                                 "Fusobacterium-")

st_object
current_microbe <- subset(st_object, subset = sampleid == "LXD_ST_pre")
View(current_microbe)
SpatialPlot(current_microbe, 
            ncol = 1,
            group.by = "Location",
            crop=T,
            image.alpha =0.3,
            pt.size.factor = 1.5)+
  scale_fill_manual(values = c("#c62d17", "#04d486", "#f6c619")) +
  labs(title=paste0(sampleid[i]))+
  theme_void()+
  theme(
    plot.title = element_text(size = 15, hjust = 0.5),
    panel.background = element_blank(),  # 去除面板背景
    plot.background = element_blank(),    # 去除整个图形的背景
    panel.grid.major = element_blank(),   # 去除主要网格线
    panel.grid.minor = element_blank(),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 13, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA,linewidth = 1)  # 添加边框
  )+
  guides(
    fill = guide_legend(
      override.aes = list(size = 3)  # Adjust the size of the legend points/icons
    )
  )  

#######绘制肿瘤边界----

for (i in seq_along(sampleid)) {
  # 动态创建子集
  assign(paste0(sampleid[i],"_RCTD合并"),readRDS(paste0("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/RCTD解卷积/解卷积合并rds2/",sampleid[i],"_RCTD合并.rds")))
  # 读取 RDS 文件，使用 paste0() 函数拼接文件名
  print(paste0("ReadRDS_RCTD合并: ", sampleid[i],"_RCTD合并"))
}

for (i in seq_along(sampleid)){
  seurat_object <- get(paste0(sampleid[i], "_RCTD合并")) 
  p<- SpatialPlot(seurat_object, 
                  ncol = 1,
                  group.by = "Location",
                  crop=T,
                  image.alpha =0.3,
                  pt.size.factor = 1.5)+
    scale_fill_manual(values = c("#c62d17", "#04d486", "#f6c619")) +
    labs(title=paste0(sampleid[i]))+
    theme_void()+
    theme(
      plot.title = element_text(size = 15, hjust = 0.5),
      panel.background = element_blank(),  # 去除面板背景
      plot.background = element_blank(),    # 去除整个图形的背景
      panel.grid.major = element_blank(),   # 去除主要网格线
      panel.grid.minor = element_blank(),
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 13, face = "bold"),
      panel.border = element_rect(color = "black", fill = NA,linewidth = 1)  # 添加边框
    )+
    guides(
      fill = guide_legend(
        override.aes = list(size = 3)  # Adjust the size of the legend points/icons
      )
    )  
  ggsave(paste0("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/Cottrazm肿瘤边界识别重新注释/绘图/",sampleid[i],"_Cottrazm.pdf"),p, width=8,height = 8)
  print(paste0("Cottrazm绘图:",sampleid[i]))
  
}

######Streptococcus----
for (i in seq_along(sampleid)){
seurat_object <- get(paste0(sampleid[i], "_RCTD合并")) 
 p2 <- SpatialFeaturePlot(
  object = seurat_object,
  features = c("Streptococcus"),  # 替换为实际的特征名称
  pt.size.factor = 1.5,
  image.alpha =0.3
) +
  scale_fill_gradientn(
    colours = c("#B6d7ce","blue","yellow", "red"),   # 定义从灰色到红色的颜色梯度                                # 未检测到的点设置为灰色
    values = scales::rescale(c(0,1,7.5,15)),                               # 设置数值范围，根据您的数据调整
    name = "Bacteria UMIs"                            # 调整图例的标题
  ) +
  theme_void() + 
  labs(title=paste0(sampleid[i],"(Streptococcus)"))+# 去除背景网格线
  theme(legend.position = "right",     
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 12),
        panel.background = element_blank(),  # 去除面板背景
        plot.background = element_blank(),    # 去除整个图形的背景
        panel.grid.major = element_blank(),   # 去除主要网格线
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 14),
        panel.border = element_rect(color = "black", fill = NA,linewidth = 1) 
        )
 ggsave(paste0("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/空转微生物SAHMI/SAHMI_Streptococcus/",sampleid[i],"_Streptococcus.pdf"),p2, width=8,height = 8)
 print(paste0("SAHMI绘图:",sampleid[i]))
}


p2 <- SpatialFeaturePlot(
  object =XYH_ST_pre_RCTD合并,
  features = c("Streptococcus"),  # 替换为实际的特征名称
  pt.size.factor = 2.6,
  image.alpha =0.3
) +
  scale_fill_gradientn(
    colours = c("#B6d7ce","blue","yellow", "red"),   # 定义从灰色到红色的颜色梯度                                # 未检测到的点设置为灰色
    values = scales::rescale(c(0,1,7.5,15)),                               # 设置数值范围，根据您的数据调整
    name = "Bacteria UMIs"                            # 调整图例的标题
  ) +
  theme_void() + 
  labs(title="XYH_ST_pre(Streptococcus)")+# 去除背景网格线
  theme(legend.position = "right",     
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 12),
        panel.background = element_blank(),  # 去除面板背景
        plot.background = element_blank(),    # 去除整个图形的背景
        panel.grid.major = element_blank(),   # 去除主要网格线
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 14),
        panel.border = element_rect(color = "black", fill = NA,linewidth = 1) 
  )
p2
ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/空转微生物SAHMI/SAHMI_Streptococcus/XYH_ST_pre_Streptococcus.pdf",p2, width=8,height = 8)


######Fusobacterium----
for (i in seq_along(sampleid)){
  seurat_object <- get(paste0(sampleid[i], "_RCTD合并")) 
  p2 <- SpatialFeaturePlot(
    object = seurat_object,
    features = c("Fusobacterium"),  # 替换为实际的特征名称
    pt.size.factor = 1.5,
    image.alpha =0.3
  ) +
    scale_fill_gradientn(
      colours = c("#B6d7ce","blue","yellow", "red"),   # 定义从灰色到红色的颜色梯度                                # 未检测到的点设置为灰色
      values = scales::rescale(c(0,1,7.5,15)),                               # 设置数值范围，根据您的数据调整
      name = "Bacteria UMIs"                            # 调整图例的标题
    ) +
    theme_void() + 
    labs(title=paste0(sampleid[i]," (Fusobacterium)"))+# 去除背景网格线
    theme(legend.position = "right",     
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 12),
          panel.background = element_blank(),  # 去除面板背景
          plot.background = element_blank(),    # 去除整个图形的背景
          panel.grid.major = element_blank(),   # 去除主要网格线
          panel.grid.minor = element_blank(),
          legend.text = element_text(size = 14),
          panel.border = element_rect(color = "black", fill = NA,linewidth = 1) 
    )
  ggsave(paste0("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/空转微生物SAHMI/SAHMI_Fusobacterium/",sampleid[i],"_Fusobacterium.pdf"),p2, width=8,height = 8)
  print(paste0("SAHMI绘图:",sampleid[i]))
}

p2 <- SpatialFeaturePlot(
  object =XYH_ST_pre_RCTD合并,
  features = c("Fusobacterium"),  # 替换为实际的特征名称
  pt.size.factor = 2.6,
  image.alpha =0.3
) +
  scale_fill_gradientn(
    colours = c("#B6d7ce","blue","yellow", "red"),   # 定义从灰色到红色的颜色梯度                                # 未检测到的点设置为灰色
    values = scales::rescale(c(0,1,7.5,15)),                               # 设置数值范围，根据您的数据调整
    name = "Bacteria UMIs"                            # 调整图例的标题
  ) +
  theme_void() + 
  labs(title="XYH_ST_pre(Fusobacterium)")+# 去除背景网格线
  theme(legend.position = "right",     
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 12),
        panel.background = element_blank(),  # 去除面板背景
        plot.background = element_blank(),    # 去除整个图形的背景
        panel.grid.major = element_blank(),   # 去除主要网格线
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 14),
        panel.border = element_rect(color = "black", fill = NA,linewidth = 1) 
  )
p2
ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/空转微生物SAHMI/SAHMI_Fusobacterium/XYH_ST_pre_Fusobacterium.pdf",p2, width=8,height = 8)


p2 <- SpatialFeaturePlot(
  object =YMS_ST_pre_RCTD合并,
  features = c("Streptococcus"),  # 替换为实际的特征名称
  pt.size.factor = 2.4,
  image.alpha =1
) +
  scale_fill_gradientn(
    colours = c("#B6d7ce","blue","yellow", "red"),   # 定义从灰色到红色的颜色梯度                                # 未检测到的点设置为灰色
    values = scales::rescale(c(0,1,7.5,15)),                               # 设置数值范围，根据您的数据调整
    name = "Bacteria UMIs"                            # 调整图例的标题
  ) +
  theme_void() + 
  labs(title="YMS_ST_pre(Streptococcus)")+# 去除背景网格线
  theme(legend.position = "right",     
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 12),
        panel.background = element_blank(),  # 去除面板背景
        plot.background = element_blank(),    # 去除整个图形的背景
        panel.grid.major = element_blank(),   # 去除主要网格线
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 14),
        panel.border = element_rect(color = "black", fill = NA,linewidth = 1) 
  )
p2
ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/空转微生物SAHMI/SAHMI_Streptococcus/YMS_ST_pre_Streptococcus.pdf",p2, width=8,height = 8)


p2 <- SpatialFeaturePlot(
  object =LXD_ST_post_RCTD合并,
  features = c("Streptococcus"),  # 替换为实际的特征名称
  pt.size.factor = 2.3,
  image.alpha =0.3
) +
  scale_fill_gradientn(
    colours = c("#B6d7ce","blue","yellow", "red"),   # 定义从灰色到红色的颜色梯度                                # 未检测到的点设置为灰色
    values = scales::rescale(c(0,1,7.5,15)),                               # 设置数值范围，根据您的数据调整
    name = "Bacteria UMIs"                            # 调整图例的标题
  ) +
  theme_void() + 
  labs(title="LXD_ST_post(Streptococcus)")+# 去除背景网格线
  theme(legend.position = "right",     
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 12),
        panel.background = element_blank(),  # 去除面板背景
        plot.background = element_blank(),    # 去除整个图形的背景
        panel.grid.major = element_blank(),   # 去除主要网格线
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 14),
        panel.border = element_rect(color = "black", fill = NA,linewidth = 1) 
  )
p2
ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/空转微生物SAHMI/SAHMI_Streptococcus/LXD_ST_post_Streptococcus.pdf",p2, width=8,height = 8)


p2 <- SpatialFeaturePlot(
  object =XL_ST_post_RCTD合并,
  features = c("Streptococcus"),  # 替换为实际的特征名称
  pt.size.factor = 2.1,
  image.alpha =0.3
) +
  scale_fill_gradientn(
    colours = c("#B6d7ce","blue","yellow", "red"),   # 定义从灰色到红色的颜色梯度                                # 未检测到的点设置为灰色
    values = scales::rescale(c(0,1,7.5,15)),                               # 设置数值范围，根据您的数据调整
    name = "Bacteria UMIs"                            # 调整图例的标题
  ) +
  theme_void() + 
  labs(title="XL_ST_post(Streptococcus)")+# 去除背景网格线
  theme(legend.position = "right",     
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 12),
        panel.background = element_blank(),  # 去除面板背景
        plot.background = element_blank(),    # 去除整个图形的背景
        panel.grid.major = element_blank(),   # 去除主要网格线
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 14),
        panel.border = element_rect(color = "black", fill = NA,linewidth = 1) 
  )
p2
ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/空转微生物SAHMI/SAHMI_Streptococcus/XL_ST_post_Streptococcus.pdf",p2, width=8,height = 8)


p2 <- SpatialFeaturePlot(
  object =XYH_ST_post_RCTD合并,
  features = c("Streptococcus"),  # 替换为实际的特征名称
  pt.size.factor = 2,
  image.alpha =0.3
) +
  scale_fill_gradientn(
    colours = c("#B6d7ce","blue","yellow", "red"),   # 定义从灰色到红色的颜色梯度                                # 未检测到的点设置为灰色
    values = scales::rescale(c(0,1,7.5,15)),                               # 设置数值范围，根据您的数据调整
    name = "Bacteria UMIs"                            # 调整图例的标题
  ) +
  theme_void() + 
  labs(title="XYH_ST_post(Streptococcus)")+# 去除背景网格线
  theme(legend.position = "right",     
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 12),
        panel.background = element_blank(),  # 去除面板背景
        plot.background = element_blank(),    # 去除整个图形的背景
        panel.grid.major = element_blank(),   # 去除主要网格线
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 14),
        panel.border = element_rect(color = "black", fill = NA,linewidth = 1) 
  )
p2
ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/空转微生物SAHMI/SAHMI_Streptococcus/XYH_ST_post_Streptococcus.pdf",p2, width=8,height = 8)




######cell_type_new----
YXJ_ST_pre_RCTD合并$top1_cell_type_new

cell_type_color <- c("#FF6347","#ffc556","#1E90FF","#0000FF","#9370DB","#00FFFF","#32cd32","#006400","#FF00ff","#20B2AA")

cell_type_color1 <- c("B cells"="#ffc556",#B cells
                      "Endothelial cells"= "#00FFFF",#"Endothelial cells"
                      "Epithelial cells"="#20B2AA",#"Epithelial cells"
                      "Fibroblasts"="#32cd32",#"Fibroblasts"
                      "Mast cells"="#9370DB",#"Mast cells"
                      "Myeloid cells"="#0000FF",#"Myeloid cells"
                      "Pericytes"="#FF00ff",#"Pericytes"
                      "Plasma cells"="#1E90FF",#"Plasma cells"
                      "Smooth muscle cells"="#006400",#"Smooth muscle cells"
                      "T cells"="#FF6347")#"T cells")

library(randomcoloR)
colors<- distinctColorPalette(17)

for (i in seq_along(sampleid)){
  seurat_object <- get(paste0(sampleid[i], "_RCTD合并")) 
  p<- SpatialPlot(seurat_object, 
                  ncol = 1,
                  group.by = "top1_cell_type_new",
                  crop=T,
                  image.alpha =0.3,
                  pt.size.factor = 1.5)+
    scale_fill_manual(values = cell_type_color1) +
    labs(title=paste0(sampleid[i]))+
    theme_void()+
    theme(
      plot.title = element_text(size = 15, hjust = 0.5),
      panel.background = element_blank(),  # 去除面板背景
      plot.background = element_blank(),    # 去除整个图形的背景
      panel.grid.major = element_blank(),   # 去除主要网格线
      panel.grid.minor = element_blank(),
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 13, face = "bold"),
      panel.border = element_rect(color = "black", fill = NA,linewidth = 1)  # 添加边框
    )+
    guides(
      fill = guide_legend(
        override.aes = list(size = 3)  # Adjust the size of the legend points/icons
      )
    )  
  ggsave(paste0("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/空转微生物SAHMI/top1_cell_type_new/",sampleid[i],"_top1_cell_type_new.pdf"),p, width=8,height = 8)
  print(paste0("绘图:",sampleid[i]))
  
}

#XYH_ST_pre_RCTD合并
XYH_ST_pre_RCTD合并$top1_cell_type_new <- replace(XYH_ST_pre_RCTD合并$top1_cell_type_new, is.na(XYH_ST_pre_RCTD合并$top1_cell_type_new), "Smooth muscle cells")

p<- SpatialPlot(XYH_ST_pre_RCTD合并, 
                ncol = 1,
                group.by = "top1_cell_type_new",
                crop=T,
                image.alpha =0.3,
                pt.size.factor = 2.9)+
  scale_fill_manual(values = cell_type_color1) +
  labs(title="XYH_ST_pre")+
  theme_void()+
  theme(
    plot.title = element_text(size = 15, hjust = 0.5),
    panel.background = element_blank(),  # 去除面板背景
    plot.background = element_blank(),    # 去除整个图形的背景
    panel.grid.major = element_blank(),   # 去除主要网格线
    panel.grid.minor = element_blank(),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 13, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA,linewidth = 1)  # 添加边框
  )+
  guides(
    fill = guide_legend(
      override.aes = list(size = 3)  # Adjust the size of the legend points/icons
    )
  )  

p
unique(XYH_ST_pre_RCTD合并$top1_cell_type_new)
ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/空转微生物SAHMI/top1_cell_type_new/XYH_ST_pre_top1_cell_type_new.pdf",p, width=8,height = 8)

YMS_ST_pre_RCTD合并$top1_cell_type_new <- replace(YMS_ST_pre_RCTD合并$top1_cell_type_new, is.na(YMS_ST_pre_RCTD合并$top1_cell_type_new), "Smooth muscle cells")

p<- SpatialPlot(YMS_ST_pre_RCTD合并, 
                ncol = 1,
                group.by = "top1_cell_type_new",
                crop=T,
                image.alpha =0.3,
                pt.size.factor = 2.4)+
  scale_fill_manual(values = cell_type_color1) +
  labs(title="YMS_ST_pre")+
  theme_void()+
  theme(
    plot.title = element_text(size = 15, hjust = 0.5),
    panel.background = element_blank(),  # 去除面板背景
    plot.background = element_blank(),    # 去除整个图形的背景
    panel.grid.major = element_blank(),   # 去除主要网格线
    panel.grid.minor = element_blank(),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 13, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA,linewidth = 1)  # 添加边框
  )+
  guides(
    fill = guide_legend(
      override.aes = list(size = 3)  # Adjust the size of the legend points/icons
    )
  )  

p
unique(YMS_ST_pre_RCTD合并$top1_cell_type_new)
ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/空转微生物SAHMI/top1_cell_type_new/YMS_ST_pre_top1_cell_type_new.pdf",p, width=8,height = 8)


p<- SpatialPlot(LXD_ST_post_RCTD合并, 
                ncol = 1,
                group.by = "top1_cell_type_new",
                crop=T,
                image.alpha =0.3,
                pt.size.factor = 2.3)+
  scale_fill_manual(values = cell_type_color1) +
  labs(title="LXD_ST_post")+
  theme_void()+
  theme(
    plot.title = element_text(size = 15, hjust = 0.5),
    panel.background = element_blank(),  # 去除面板背景
    plot.background = element_blank(),    # 去除整个图形的背景
    panel.grid.major = element_blank(),   # 去除主要网格线
    panel.grid.minor = element_blank(),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 13, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA,linewidth = 1)  # 添加边框
  )+
  guides(
    fill = guide_legend(
      override.aes = list(size = 3)  # Adjust the size of the legend points/icons
    )
  )  

p
unique(YMS_ST_pre_RCTD合并$top1_cell_type_new)
ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/空转微生物SAHMI/top1_cell_type_new/LXD_ST_post_top1_cell_type_new.pdf",p, width=8,height = 8)

p<- SpatialPlot(XL_ST_post_RCTD合并, 
                ncol = 1,
                group.by = "top1_cell_type_new",
                crop=T,
                image.alpha =0.3,
                pt.size.factor = 2.2)+
  scale_fill_manual(values = cell_type_color1) +
  labs(title="XL_ST_post")+
  theme_void()+
  theme(
    plot.title = element_text(size = 15, hjust = 0.5),
    panel.background = element_blank(),  # 去除面板背景
    plot.background = element_blank(),    # 去除整个图形的背景
    panel.grid.major = element_blank(),   # 去除主要网格线
    panel.grid.minor = element_blank(),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 13, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA,linewidth = 1)  # 添加边框
  )+
  guides(
    fill = guide_legend(
      override.aes = list(size = 3)  # Adjust the size of the legend points/icons
    )
  )  

p
ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/空转微生物SAHMI/top1_cell_type_new/XL_ST_post_top1_cell_type_new.pdf",p, width=8,height = 8)

XYH_ST_post_RCTD合并$top1_cell_type_new <- replace(XYH_ST_post_RCTD合并$top1_cell_type_new, is.na(XYH_ST_post_RCTD合并$top1_cell_type_new), "Fibroblasts")


p<- SpatialPlot(XYH_ST_post_RCTD合并, 
                ncol = 1,
                group.by = "top1_cell_type_new",
                crop=T,
                image.alpha =0.3,
                pt.size.factor = 2)+
  scale_fill_manual(values = cell_type_color1) +
  labs(title=paste0("XYH_ST_post"))+
  theme_void()+
  theme(
    plot.title = element_text(size = 15, hjust = 0.5),
    panel.background = element_blank(),  # 去除面板背景
    plot.background = element_blank(),    # 去除整个图形的背景
    panel.grid.major = element_blank(),   # 去除主要网格线
    panel.grid.minor = element_blank(),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 13, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA,linewidth = 1)  # 添加边框
  )+
  guides(
    fill = guide_legend(
      override.aes = list(size = 3)  # Adjust the size of the legend points/icons
    )
  )  

p
ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/空转微生物SAHMI/top1_cell_type_new/XYH_ST_post_top1_cell_type_new.pdf",p, width=8,height = 8)


XL_ST_pre_RCTD合并$top1_cell_type_new <- replace(XL_ST_pre_RCTD合并$top1_cell_type_new, is.na(XL_ST_pre_RCTD合并$top1_cell_type_new), "Smooth muscle cells")


p<- SpatialPlot(XL_ST_pre_RCTD合并, 
                ncol = 1,
                group.by = "top1_cell_type_new",
                crop=T,
                image.alpha =0.3,
                pt.size.factor = 1.6)+
  scale_fill_manual(values = cell_type_color1) +
  labs(title=paste0("XL_ST_pre"))+
  theme_void()+
  theme(
    plot.title = element_text(size = 15, hjust = 0.5),
    panel.background = element_blank(),  # 去除面板背景
    plot.background = element_blank(),    # 去除整个图形的背景
    panel.grid.major = element_blank(),   # 去除主要网格线
    panel.grid.minor = element_blank(),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 13, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA,linewidth = 1)  # 添加边框
  )+
  guides(
    fill = guide_legend(
      override.aes = list(size = 3)  # Adjust the size of the legend points/icons
    )
  )  

p
ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/空转微生物SAHMI/top1_cell_type_new/XL_ST_pre_top1_cell_type_new.pdf",p, width=8,height = 8)


#########################################


#XYH_ST_pre_RCTD合并
XYH_ST_pre_RCTD合并$top1_cell_type_new <- replace(XYH_ST_pre_RCTD合并$top1_cell_type_new, is.na(XYH_ST_pre_RCTD合并$top1_cell_type_new), "Smooth muscle cells")

p<- SpatialPlot(XYH_ST_pre_RCTD合并, 
                ncol = 1,
                group.by = "Location",
                crop=T,
                image.alpha =0.3,
                pt.size.factor = 2.9)+
  scale_fill_manual(values =  c("#c62d17", "#04d486", "#f6c619")) +
  labs(title="XYH_ST_pre")+
  theme_void()+
  theme(
    plot.title = element_text(size = 15, hjust = 0.5),
    panel.background = element_blank(),  # 去除面板背景
    plot.background = element_blank(),    # 去除整个图形的背景
    panel.grid.major = element_blank(),   # 去除主要网格线
    panel.grid.minor = element_blank(),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 13, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA,linewidth = 1)  # 添加边框
  )+
  guides(
    fill = guide_legend(
      override.aes = list(size = 3)  # Adjust the size of the legend points/icons
    )
  )  

p
ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/Cottrazm肿瘤边界识别重新注释/绘图/XYH_ST_pre_Cottrazm.pdf",p, width=8,height = 8)



p<- SpatialPlot(YMS_ST_pre_RCTD合并, 
                ncol = 1,
                group.by = "Location",
                crop=T,
                image.alpha =0.3,
                pt.size.factor = 2.4)+
  scale_fill_manual(values = c("#c62d17", "#04d486", "#f6c619")) +
  labs(title="YMS_ST_pre")+
  theme_void()+
  theme(
    plot.title = element_text(size = 15, hjust = 0.5),
    panel.background = element_blank(),  # 去除面板背景
    plot.background = element_blank(),    # 去除整个图形的背景
    panel.grid.major = element_blank(),   # 去除主要网格线
    panel.grid.minor = element_blank(),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 13, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA,linewidth = 1)  # 添加边框
  )+
  guides(
    fill = guide_legend(
      override.aes = list(size = 3)  # Adjust the size of the legend points/icons
    )
  )  

p
unique(YMS_ST_pre_RCTD合并$top1_cell_type_new)
ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/Cottrazm肿瘤边界识别重新注释/绘图/YMS_ST_pre_Cottrazm.pdf",p, width=8,height = 8)


p<- SpatialPlot(LXD_ST_post_RCTD合并, 
                ncol = 1,
                group.by = "Location",
                crop=T,
                image.alpha =0.3,
                pt.size.factor = 2.3)+
  scale_fill_manual(values = c("#c62d17", "#04d486", "#f6c619")) +
  labs(title="LXD_ST_post")+
  theme_void()+
  theme(
    plot.title = element_text(size = 15, hjust = 0.5),
    panel.background = element_blank(),  # 去除面板背景
    plot.background = element_blank(),    # 去除整个图形的背景
    panel.grid.major = element_blank(),   # 去除主要网格线
    panel.grid.minor = element_blank(),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 13, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA,linewidth = 1)  # 添加边框
  )+
  guides(
    fill = guide_legend(
      override.aes = list(size = 3)  # Adjust the size of the legend points/icons
    )
  )  

p
unique(YMS_ST_pre_RCTD合并$top1_cell_type_new)
ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/Cottrazm肿瘤边界识别重新注释/绘图/LXD_ST_post_Cottrazm.pdf",p, width=8,height = 8)

p<- SpatialPlot(XL_ST_post_RCTD合并, 
                ncol = 1,
                group.by = "Location",
                crop=T,
                image.alpha =0.3,
                pt.size.factor = 2.2)+
  scale_fill_manual(values = c("#c62d17", "#04d486", "#f6c619")) +
  labs(title="XL_ST_post")+
  theme_void()+
  theme(
    plot.title = element_text(size = 15, hjust = 0.5),
    panel.background = element_blank(),  # 去除面板背景
    plot.background = element_blank(),    # 去除整个图形的背景
    panel.grid.major = element_blank(),   # 去除主要网格线
    panel.grid.minor = element_blank(),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 13, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA,linewidth = 1)  # 添加边框
  )+
  guides(
    fill = guide_legend(
      override.aes = list(size = 3)  # Adjust the size of the legend points/icons
    )
  )  

p
ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/Cottrazm肿瘤边界识别重新注释/绘图/XL_ST_post_Cottrazm.pdf",p, width=8,height = 8)

XYH_ST_post_RCTD合并$top1_cell_type_new <- replace(XYH_ST_post_RCTD合并$top1_cell_type_new, is.na(XYH_ST_post_RCTD合并$top1_cell_type_new), "Fibroblasts")


p<- SpatialPlot(XYH_ST_post_RCTD合并, 
                ncol = 1,
                group.by = "Location",
                crop=T,
                image.alpha =0.3,
                pt.size.factor = 2)+
  scale_fill_manual(values = c("#c62d17", "#04d486", "#f6c619")) +
  labs(title=paste0("XYH_ST_post"))+
  theme_void()+
  theme(
    plot.title = element_text(size = 15, hjust = 0.5),
    panel.background = element_blank(),  # 去除面板背景
    plot.background = element_blank(),    # 去除整个图形的背景
    panel.grid.major = element_blank(),   # 去除主要网格线
    panel.grid.minor = element_blank(),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 13, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA,linewidth = 1)  # 添加边框
  )+
  guides(
    fill = guide_legend(
      override.aes = list(size = 3)  # Adjust the size of the legend points/icons
    )
  )  

p
ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/Cottrazm肿瘤边界识别重新注释/绘图/XYH_ST_post_Cottrazm.pdf",p, width=8,height = 8)


XL_ST_pre_RCTD合并$top1_cell_type_new <- replace(XL_ST_pre_RCTD合并$top1_cell_type_new, is.na(XL_ST_pre_RCTD合并$top1_cell_type_new), "Smooth muscle cells")


p<- SpatialPlot(XL_ST_pre_RCTD合并, 
                ncol = 1,
                group.by = "Location",
                crop=T,
                image.alpha =0.3,
                pt.size.factor = 1.6)+
  scale_fill_manual(values = c("#c62d17", "#04d486", "#f6c619")) +
  labs(title=paste0("XL_ST_pre"))+
  theme_void()+
  theme(
    plot.title = element_text(size = 15, hjust = 0.5),
    panel.background = element_blank(),  # 去除面板背景
    plot.background = element_blank(),    # 去除整个图形的背景
    panel.grid.major = element_blank(),   # 去除主要网格线
    panel.grid.minor = element_blank(),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 13, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA,linewidth = 1)  # 添加边框
  )+
  guides(
    fill = guide_legend(
      override.aes = list(size = 3)  # Adjust the size of the legend points/icons
    )
  )  

p
ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/Cottrazm肿瘤边界识别重新注释/绘图/XL_ST_pre_Cottrazm.pdf",p, width=8,height = 8)



######bacteria----
YXJ_ST_pre_RCTD合并$top1_celltype_bacteria
unique(YXJ_ST_pre_RCTD合并$top1_celltype_bacteria)
colors<- c("#f4a460","#6ca6cd")

colors <-c("Endothelial cells_Bacteria+"="#00FFFF",#"Endothelial cells"
           "Epithelial cells_Bacteria+" ="#20B2AA",#"Epithelial cells"
           "Fibroblasts_Bacteria+"="#32cd32",#"Fibroblasts"
           "Myeloid cells_Bacteria+"="#0000FF",#'Myeloid cells'
           "Pericytes_Bacteria+"="#FF00ff",#"Pericytes"
           "Plasma cells_Bacteria+" = "#1E90FF",#'Plasma cells'
           "Smooth muscle cells_Bacteria+"="#006400",#"Smooth muscle cells"
           "T cells_Bacteria+"= "#FF6347",#T cells
           "Bacteria-"="gray")


colors<- distinctColorPalette(9)

for (i in seq_along(sampleid)){
  seurat_object <- get(paste0(sampleid[i], "_RCTD合并")) 
  seurat_object$microbe_group_scrna <- ifelse(grepl("Bacteria[+]", seurat_object$top1_celltype_bacteria),
                                              seurat_object$top1_celltype_bacteria, 
                                              "Bacteria-")
  p<- SpatialPlot(seurat_object, 
                  ncol = 1,
                  group.by = "microbe_group_scrna",
                  crop=T,
                  image.alpha =0.3,
                  pt.size.factor = 1.5)+
    scale_fill_manual(values = colors) +
    labs(title=paste0(sampleid[i]))+
    theme_void()+
    theme(
      plot.title = element_text(size = 15, hjust = 0.5),
      panel.background = element_blank(),  # 去除面板背景
      plot.background = element_blank(),    # 去除整个图形的背景
      panel.grid.major = element_blank(),   # 去除主要网格线
      panel.grid.minor = element_blank(),
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 13, face = "bold"),
      panel.border = element_rect(color = "black", fill = NA,linewidth = 1)  # 添加边框
    )+
    guides(
      fill = guide_legend(
        override.aes = list(size = 3)  # Adjust the size of the legend points/icons
      )
    )  
  ggsave(paste0("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/空转微生物SAHMI/top1_microbe_group_scrna/",sampleid[i],"_top1_microbe_group_scrna.pdf"),p, width=8,height = 8)
  print(paste0("绘图:",sampleid[i]))
  
}

#####
miMPRobe_data<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/PDCD1每个细胞以及微生物分组表达情况/16s_filter.combined.genus.umi.matrix.csv")

####
######Streptococcus_group----
YXJ_ST_pre_RCTD合并$top1_Streptococcus_group
unique(YXJ_ST_pre_RCTD合并$top1_Streptococcus_group)

colors<- c("#ffdc91b2","#91d1c2b2")

colors<- distinctColorPalette(9)

for (i in seq_along(sampleid)){
  seurat_object <- get(paste0(sampleid[i], "_RCTD合并")) 
  p<- SpatialPlot(seurat_object, 
                  ncol = 1,
                  group.by = "top1_Streptococcus_group",
                  crop=T,
                  image.alpha =0.3,
                  pt.size.factor = 1.5)+
    scale_fill_manual(values = colors) +
    labs(title=paste0(sampleid[i]))+
    theme_void()+
    theme(
      plot.title = element_text(size = 15, hjust = 0.5),
      panel.background = element_blank(),  # 去除面板背景
      plot.background = element_blank(),    # 去除整个图形的背景
      panel.grid.major = element_blank(),   # 去除主要网格线
      panel.grid.minor = element_blank(),
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 13, face = "bold"),
      panel.border = element_rect(color = "black", fill = NA,linewidth = 1)  # 添加边框
    )+
    guides(
      fill = guide_legend(
        override.aes = list(size = 3)  # Adjust the size of the legend points/icons
      )
    )  
  ggsave(paste0("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/空转微生物SAHMI/top1_Streptococcus_group/",sampleid[i],"_top1_Streptococcus_group.pdf"),p, width=8,height = 8)
  print(paste0("绘图:",sampleid[i]))
  
}


######Fusobacterium_group----
YXJ_ST_pre_RCTD合并$top1_Fusobacterium_group
unique(YXJ_ST_pre_RCTD合并$top1_Fusobacterium_group)

colors<- c("#ffdc91b2","#AB82FF")

colors<- distinctColorPalette(9)

for (i in seq_along(sampleid)){
  seurat_object <- get(paste0(sampleid[i], "_RCTD合并")) 
  p<- SpatialPlot(seurat_object, 
                  ncol = 1,
                  group.by = "top1_Fusobacterium_group",
                  crop=T,
                  image.alpha =0.3,
                  pt.size.factor = 1.5)+
    scale_fill_manual(values = colors) +
    labs(title=paste0(sampleid[i]))+
    theme_void()+
    theme(
      plot.title = element_text(size = 15, hjust = 0.5),
      panel.background = element_blank(),  # 去除面板背景
      plot.background = element_blank(),    # 去除整个图形的背景
      panel.grid.major = element_blank(),   # 去除主要网格线
      panel.grid.minor = element_blank(),
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 13, face = "bold"),
      panel.border = element_rect(color = "black", fill = NA,linewidth = 1)  # 添加边框
    )+
    guides(
      fill = guide_legend(
        override.aes = list(size = 3)  # Adjust the size of the legend points/icons
      )
    )  
  ggsave(paste0("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/空转微生物SAHMI/top1_Fusobacterium_group/",sampleid[i],"_top1_Fusobacterium_group.pdf"),p, width=8,height = 8)
  print(paste0("绘图:",sampleid[i]))
  
}


#####Microbe_UMIs_featureplot----
YXJ_ST_pre_microbe$Microbe_UMIs
YXJ_ST_pre_microbe$CNVLabel


for (i in seq_along(sampleid)){
  seurat_object <- get(paste0(sampleid[i], "_microbe")) 
  p2 <- SpatialFeaturePlot(
    object = seurat_object,
    features = c("Microbe_UMIs"),  # 替换为实际的特征名称
    pt.size.factor = 1.5,
    image.alpha =0.3
  ) +
    scale_fill_gradientn(
      colours = c("#B6d7ce","blue","yellow", "red"),   # 定义从灰色到红色的颜色梯度                                # 未检测到的点设置为灰色
      values = scales::rescale(c(0,1,7.5,15)),                               # 设置数值范围，根据您的数据调整
      name = "Bacteria UMIs"                            # 调整图例的标题
    ) +
    theme_void() + 
    labs(title=paste0(sampleid[i],"(microbe_UMIs)"))+# 去除背景网格线
    theme(legend.position = "right",     
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 12),
          panel.background = element_blank(),  # 去除面板背景
          plot.background = element_blank(),    # 去除整个图形的背景
          panel.grid.major = element_blank(),   # 去除主要网格线
          panel.grid.minor = element_blank(),
          legend.text = element_text(size = 14),
          panel.border = element_rect(color = "black", fill = NA,linewidth = 1) 
    )
  ggsave(paste0("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/空转微生物SAHMI/SAHMI_Microbe_UMIs/",sampleid[i],"_Microbe_UMIs.pdf"),p2, width=8,height = 8)
  print(paste0("SAHMI绘图:",sampleid[i]))
}


XYH_ST_pre_microbe$Microbe_UMIs<- as.numeric(XYH_ST_pre_microbe$Microbe_UMIs)

p<- SpatialPlot(XYH_ST_pre_microbe, 
                ncol = 1,
                features = "Microbe_UMIs",
                image.alpha =0.3,
                pt.size.factor = 2.8)+
  scale_fill_gradientn(
    colours = c("#B6d7ce","blue","yellow", "red"),   # 定义从灰色到红色的颜色梯度                                # 未检测到的点设置为灰色
    values = scales::rescale(c(0,1,7.5,15)),                               # 设置数值范围，根据您的数据调整
    name = "Bacteria UMIs"                            # 调整图例的标题
  ) +
  labs(title="XYH_ST_pre")+
  theme_void()+
  theme(
    plot.title = element_text(size = 15, hjust = 0.5),
    panel.background = element_blank(),  # 去除面板背景
    plot.background = element_blank(),    # 去除整个图形的背景
    panel.grid.major = element_blank(),   # 去除主要网格线
    panel.grid.minor = element_blank(),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 13, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA,linewidth = 1)  # 添加边框
  )

p
ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/空转微生物SAHMI/SAHMI_Microbe_UMIs/XYH_ST_pre_Microbe_UMIs.pdf",p, width=8,height = 8)



p<- SpatialPlot(YMS_ST_pre_microbe, 
                ncol = 1,
                features = "Microbe_UMIs",
                crop=T,
                image.alpha =0.3,
                pt.size.factor = 2.4)+
  scale_fill_gradientn(
    colours = c("#B6d7ce","blue","yellow", "red"),   # 定义从灰色到红色的颜色梯度                                # 未检测到的点设置为灰色
    values = scales::rescale(c(0,1,7.5,15)),                               # 设置数值范围，根据您的数据调整
    name = "Bacteria UMIs"                            # 调整图例的标题
  ) +
  labs(title="YMS_ST_pre")+
  theme_void()+
  theme(
    plot.title = element_text(size = 15, hjust = 0.5),
    panel.background = element_blank(),  # 去除面板背景
    plot.background = element_blank(),    # 去除整个图形的背景
    panel.grid.major = element_blank(),   # 去除主要网格线
    panel.grid.minor = element_blank(),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 13, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA,linewidth = 1)  # 添加边框
  )

p

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/空转微生物SAHMI/SAHMI_Microbe_UMIs/YMS_ST_pre_Microbe_UMIs.pdf",p, width=8,height = 8)


p<- SpatialPlot(LXD_ST_post_microbe, 
                ncol = 1,
                features = "Microbe_UMIs",
                crop=T,
                image.alpha =0.3,
                pt.size.factor = 2.3)+
  scale_fill_gradientn(
    colours = c("#B6d7ce","blue","yellow", "red"),   # 定义从灰色到红色的颜色梯度                                # 未检测到的点设置为灰色
    values = scales::rescale(c(0,1,7.5,15)),                               # 设置数值范围，根据您的数据调整
    name = "Bacteria UMIs"                            # 调整图例的标题
  ) +
  labs(title="LXD_ST_post")+
  theme_void()+
  theme(
    plot.title = element_text(size = 15, hjust = 0.5),
    panel.background = element_blank(),  # 去除面板背景
    plot.background = element_blank(),    # 去除整个图形的背景
    panel.grid.major = element_blank(),   # 去除主要网格线
    panel.grid.minor = element_blank(),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 13, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA,linewidth = 1)  # 添加边框
  )

p
ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/空转微生物SAHMI/SAHMI_Microbe_UMIs/LXD_ST_post_Microbe_UMIs.pdf",p, width=8,height = 8)

p<- SpatialPlot(XL_ST_post_microbe, 
                ncol = 1,
                features = "Microbe_UMIs",
                crop=T,
                image.alpha =0.3,
                pt.size.factor = 2.2)+
  scale_fill_gradientn(
    colours = c("#B6d7ce","blue","yellow", "red"),   # 定义从灰色到红色的颜色梯度                                # 未检测到的点设置为灰色
    values = scales::rescale(c(0,1,7.5,15)),                               # 设置数值范围，根据您的数据调整
    name = "Bacteria UMIs"                            # 调整图例的标题
  ) +
  labs(title="XL_ST_post")+
  theme_void()+
  theme(
    plot.title = element_text(size = 15, hjust = 0.5),
    panel.background = element_blank(),  # 去除面板背景
    plot.background = element_blank(),    # 去除整个图形的背景
    panel.grid.major = element_blank(),   # 去除主要网格线
    panel.grid.minor = element_blank(),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 13, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA,linewidth = 1)  # 添加边框
  )
p
ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/空转微生物SAHMI/SAHMI_Microbe_UMIs/XL_ST_post_Microbe_UMIs.pdf",p, width=8,height = 8)


p<- SpatialPlot(XYH_ST_post_microbe, 
                ncol = 1,
                features = "Microbe_UMIs",
                crop=T,
                image.alpha =0.3,
                pt.size.factor = 2)+
  scale_fill_gradientn(
    colours = c("#B6d7ce","blue","yellow", "red"),   # 定义从灰色到红色的颜色梯度                                # 未检测到的点设置为灰色
    values = scales::rescale(c(0,1,7.5,15)),                               # 设置数值范围，根据您的数据调整
    name = "Bacteria UMIs"                            # 调整图例的标题
  ) +
  labs(title=paste0("XYH_ST_post"))+
  theme_void()+
  theme(
    plot.title = element_text(size = 15, hjust = 0.5),
    panel.background = element_blank(),  # 去除面板背景
    plot.background = element_blank(),    # 去除整个图形的背景
    panel.grid.major = element_blank(),   # 去除主要网格线
    panel.grid.minor = element_blank(),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 13, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA,linewidth = 1)  # 添加边框
  )  

p
ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/空转微生物SAHMI/SAHMI_Microbe_UMIs/XYH_ST_post_Microbe_UMIs.pdf",p, width=8,height = 8)




p<- SpatialPlot(XL_ST_pre_microbe, 
                ncol = 1,
                features = "Microbe_UMIs",
                crop=T,
                image.alpha =0.3,
                pt.size.factor = 1.6)+
  scale_fill_gradientn(
    colours = c("#B6d7ce","blue","yellow", "red"),   # 定义从灰色到红色的颜色梯度                                # 未检测到的点设置为灰色
    values = scales::rescale(c(0,1,7.5,15)),                               # 设置数值范围，根据您的数据调整
    name = "Bacteria UMIs"                            # 调整图例的标题
  ) +
  labs(title=paste0("XL_ST_pre"))+
  theme_void()+
  theme(
    plot.title = element_text(size = 15, hjust = 0.5),
    panel.background = element_blank(),  # 去除面板背景
    plot.background = element_blank(),    # 去除整个图形的背景
    panel.grid.major = element_blank(),   # 去除主要网格线
    panel.grid.minor = element_blank(),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 13, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA,linewidth = 1)  # 添加边框
  )

p
ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/空转微生物SAHMI/SAHMI_Microbe_UMIs//XL_ST_pre_Microbe_UMIs.pdf",p, width=8,height = 8)

####组织绘图----

for (i in seq_along(sampleid)){
  seurat_object <- get(paste0(sampleid[i], "_microbe")) 
  p2 <- SpatialFeaturePlot(
    object = seurat_object,
    features = c("Microbe_UMIs"),  # 替换为实际的特征名称
    pt.size.factor = 0,
    image.alpha =1
  ) +
    theme_void() + 
    labs(title=paste0(sampleid[i],"(microbe_UMIs)"))+# 去除背景网格线
    theme(    
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 12),
          panel.background = element_blank(),  # 去除面板背景
          plot.background = element_blank(),    # 去除整个图形的背景
          panel.grid.major = element_blank(),   # 去除主要网格线
          panel.grid.minor = element_blank(),
          legend.text = element_text(size = 14),
          panel.border = element_blank(),
          legend.position ="none"
    )
  ggsave(paste0("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/空转微生物SAHMI/SAHMI_组织图/",sampleid[i],"_组织图.pdf"),p2, width=8,height = 8)
  print(paste0("SAHMI绘图:",sampleid[i]))
}

XYH_ST_post_microbe <- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/st_data_microbe/XYH_ST_post_microbe.rds")

p2 <- SpatialFeaturePlot(
  object = XYH_ST_post_microbe,
  features = c("Microbe_UMIs"),  # 替换为实际的特征名称
  pt.size.factor = 0,
  image.alpha =1
) +
  theme_void() + 
  labs(title="XYH_ST_post")+# 去除背景网格线
  theme(    
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    panel.background = element_blank(),  # 去除面板背景
    plot.background = element_blank(),    # 去除整个图形的背景
    panel.grid.major = element_blank(),   # 去除主要网格线
    panel.grid.minor = element_blank(),
    legend.text = element_text(size = 14),
    panel.border = element_blank(),
    legend.position ="none"
  )

p2

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/空转微生物SAHMI/SAHMI_组织图/XYH_ST_post_组织图.pdf",p2, width=8,height = 8)


p<- SpatialPlot(XL_ST_pre_microbe, 
                ncol = 1,
                features = "Microbe_UMIs",
                crop=T,
                image.alpha =0.3,
                pt.size.factor = 0) +
  labs(title=paste0("XL_ST_pre"))+
  theme_void()+
  theme(
    plot.title = element_text(size = 15, hjust = 0.5),
    panel.background = element_blank(),  # 去除面板背景
    plot.background = element_blank(),    # 去除整个图形的背景
    panel.grid.major = element_blank(),   # 去除主要网格线
    panel.grid.minor = element_blank(),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 13, face = "bold"),
    panel.border = element_blank(),
    legend.position ="none",# 添加边框
  )

p


####统计微生物占比----
View(st_object)

metadata<- st_object@meta.data

names(metadata)

metadata$Total_Microbe_UMIs_Per_Sample

metadata$dominant_genus

ncol(st_object)
metadata$status
library(dplyr)
metadata %>%
  summarise(status2 = mean(status == "Microbe detected") * 100)

names(metadata)
metadata_microbe<- metadata[,c(27:71)]

metadata_microbe$sum <- rowSums(metadata_microbe, na.rm = TRUE)  # 对每一行求和，去除最后一列（如果最后一列是group列）
metadata_microbe$group <- ifelse(metadata_microbe$sum > 3, "detected", "non-detected")
count_detected <- sum(metadata_microbe$group == "detected") / nrow(metadata_microbe) * 100
max(metadata_microbe$sum)
min(metadata_microbe$sum)
max(metadata$Microbe_UMIs)
median(metadata$Microbe_UMIs)

sampleid<-c("LXD_ST_pre","LXD_ST_post","XL_ST_pre","XL_ST_post","XYH_ST_post",
            "XYH_ST_pre","YMS_ST_post","YMS_ST_pre","YXJ_ST_pre","LM_ST_pre")


count_detected <- sum(LXD_ST_pre_microbe$microbe_detected_status == "Microbe detected")/ncol(LXD_ST_pre_microbe) *100


sampleid <- c("LXD_ST_pre","LXD_ST_post","XL_ST_pre","XL_ST_post","XYH_ST_post",
              "XYH_ST_pre","YMS_ST_post","YMS_ST_pre","YXJ_ST_pre","LM_ST_pre")

for (i in seq_along(sampleid)) {
  # 动态创建子集
  assign(paste0(sampleid[i],"_microbe"),readRDS(paste0("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/空转数据/st_data_microbe/",sampleid[i],"_microbe.rds")))
  # 读取 RDS 文件，使用 paste0() 函数拼接文件名
  print(paste0("ReadRDS_microbe: ", sampleid[i]))
}

sample_data_list <- list(LXD_ST_pre_microbe, LXD_ST_post_microbe, XL_ST_pre_microbe, 
                         XL_ST_post_microbe, XYH_ST_post_microbe, XYH_ST_pre_microbe, 
                         YMS_ST_post_microbe, YMS_ST_pre_microbe, YXJ_ST_pre_microbe, 
                         LM_ST_pre_microbe)

# 循环遍历每个样本数据框
for (i in 1:length(sampleid)) {
  # 获取当前样本的数据框
  current_sample_data <- sample_data_list[[i]]
  
  # 计算 Microbe detected 的比例
  count_detected <- sum(current_sample_data$microbe_detected_status == "Microbe detected") / ncol(current_sample_data) * 100
  
  # 打印结果
  cat(paste("Sample:", sampleid[i], " - Microbe detected percentage:", count_detected, "%\n"))
}
