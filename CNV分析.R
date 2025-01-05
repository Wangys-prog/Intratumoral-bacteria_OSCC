####
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("infercnv",force = TRUE)
# 
# options(repos='http://cran.rstudio.com/')
library(infercnv)
library(Seurat)
library(tidyverse)
library(harmony)
library(dplyr)
library(patchwork)
library(SeuratWrappers)
library(SeuratData)

MPRNMPR.SC.subset<- readRDS("./MPRNMPR.SC.subset_v1_harmony.rds")
MPRNMPR.SC.subset = UpdateSeuratObject(MPRNMPR.SC.subset)


##单细胞表达量的原始矩阵
counts_matrix<- Seurat::GetAssayData(MPRNMPR.SC.subset, slot="counts")
##细胞类型注释文件

cellanno<-Seurat:: FetchData(MPRNMPR.SC.subset,vars="cell_type_new")%>%tibble::rownames_to_column(var="cellbarcode")
write.table(cellanno,"cnv_celltype_group.xls",sep="\t",col.names=F,row.names=F,quote=F)

####基因排序文件自己构建
dat <-counts_matrix
library(AnnoProbe)
geneInfor=annoGene(rownames(dat),"SYMBOL","human")
colnames(geneInfor)
geneInfor=geneInfor[with(geneInfor, order(chr, start)),c(1,4:6)]
geneInfor=geneInfor[!duplicated(geneInfor[,1]),]
length(unique(geneInfor[,1]))
# head(geneInfor)
# 这里可以去除性染色体
#也可以把染色体排序方式改变
dat=dat[rownames(dat) %in% geneInfor[,1],]
dat=dat[match( geneInfor[,1], rownames(dat) ),] 
# dim(dat)
# head(geneInfor)
geneFile='./geneFile.txt'
write.table(geneInfor,file = geneFile,sep = '\t',quote = F,col.names = F,row.names = F)


# ##基因排序文件
# gtf="./refdata-gex-GRCh38-2020-A/genes.gtf"
# gtf=plyranges::read_gff(gtf)
# gene.chr=gtf%>%
#   plyranges::filter(type=="gene"&gene_name%in%rownames(MPRNMPR.SC.subset))%>%
#   as.data.frame()%>%
#   dplyr::select(gene_name,seqnames,start,end)%>%
#   dplyr::distinct(gene_name,.keep_all=T)
# write.table(gene.chr,"gene_order_file.xls",col.names=F,row.names=F,sep="\t",quote=F)

###创建inferCNV对象

infercnv_obj=CreateInfercnvObject(raw_counts_matrix=counts_matrix,#可以直接提供矩阵对象
                                  annotations_file="cnv_celltype_group.xls",
                                  delim="\t",
                                  gene_order_file="geneFile.txt",
                                  ref_group_names=NULL)# 用于设置参考组，正常的细胞类型组别
								  

sampled_cellmeta = seurat_ob@meta.data[which(seurat_ob@meta.data[,'celltype'] %in% 'T_cells'),] %>% tibble::rownames_to_column() %>%
                            dplyr::group_by( .dots= "seurat_clusters" ) %>%
                            dplyr::sample_frac( size = 0.75,replace = F) %>%
                            tibble::column_to_rownames()
subsetcell=setdiff(colnames(seurat_ob),rownames(sampled_cellmeta))
seurat_ob = seurat_ob[, subsetcell]
								  
infercnv_obj=CreateInfercnvObject(raw_counts_matrix=counts_matrix,#可以直接提供矩阵对象
                                  annotations_file="cnv_celltype_group.xls",
                                  delim="\t",
                                  gene_order_file="geneFile.txt",
                                  ref_group_names=NULL)# 用于设置参考组，正常的细胞类型组别								  
								  

###运行流程
infercnv_obj=infercnv::run(infercnv_obj,
                           cutoff=0.1,#use1 for smart-seq, 0.1 for 10x-genomics 
                           out_dir="infercnv_output_dir2",# dir is auto-created for storing outputs 
                           analysis_mode="samples", 
                           cluster_by_groups=T,# 根据细胞类型对肿瘤细胞执行单独的聚类
                           denoise=T,# 去噪处理
                           num_threads=10,#设置线程数,多线程运行，加快计算速度
                           HMM=F)#是否基于HMM预测CNV,选择F会加快运算速度
