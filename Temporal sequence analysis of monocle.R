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

####CD4拟时序----
T_cells_object <- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/T_cells_object.rds")
CD4_object = subset(T_cells_object,idents = c("CD4_C1-PDCD1","CD4_C2-IL7R","CD4_C3-FOXP3","CD4_C4-FOSB"))
CD4_object$RNA_snn_res.0.8
saveRDS(CD4_object,"D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/T_cells_CD4_object.rds")
saveRDS(CD8_object,"D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/T_cells_CD8_object.rds")
CD4_object$cd4_cd8_group
sce_CDS <- ks_run_Monocle2(object = CD4_object,
                           layer = 'counts',
                           assay = "RNA",
                           VARgenesM="seurat",
                           cellAnno = "cd4_cd8_group")

saveRDS(sce_CDS,"D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/T_cells_CD4拟时序分析cds数据.rds")

sce_CDS<-readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/T_cells_CD4拟时序分析cds数据.rds")

plot_cell_trajectory(sce_CDS,color_by = "cd4_cd8_group")
dev.off()
plot_cell_trajectory(sce_CDS, color_by = "Massilia")+
  scale_color_gradientn(colours = colorRampPalette(c('blue','green','yellow','red'))(100))

pdf(file="D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/拟时序分析/T细胞拟时序分析/CD8_拟时序/af_trajectory.Pseudotime.pdf",width=6.5,height=6)
plot_cell_trajectory(sce_CDS,color_by = "Pseudotime")
dev.off()
#State轨迹分布图
pdf(file="D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/拟时序分析/T细胞拟时序分析/CD8_拟时序/af_trajectory.state.pdf",width=6.5,height=6)
p1 = plot_cell_trajectory(cds, color_by = "State")
#分面展示
p2 = plot_cell_trajectory(cds, color_by = "State") + facet_wrap(~State, nrow = 1)
p1|p2
dev.off()
sce_CDS$Massilia
sce_CDS$Herbaspirillum
sce_CDS$Fusobacterium
sce_CDS$Herbaspirillum_group
plot_cell_trajectory(sce_CDS,color_by = "Herbaspirillum_group",grouping = "State")

s.genes <- c("Fusobacterium","Massilia","Herbaspirillum")
# 点图（抖动）
p1 <- plot_genes_jitter(sce_CDS$Fusobacterium, grouping = "State", color_by = "State")
# 小提琴图
p2 <- plot_genes_violin(sce_CDS[s.genes,], grouping = "State", color_by = "State")
# 伪时间图
p3 <- plot_genes_in_pseudotime(cds[s.genes,], color_by = "State")
plotc <- p1|p2|p3

####CD8拟时序----
BiocManager::install("monocle",force = TRUE)
library(monocle)
library(Seurat)
library(tidyverse)
library(harmony)
library(dplyr)
library(patchwork)
library(SeuratWrappers)
# 读入文件
T_cells_object <- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/T_cells_object.rds")
# CD8_object = subset(T_cells_object,idents = c("CD8_C1-GNLY","CD8_C2-CTSW","CD8_C3-IL7R","CD8_C4-GZMK","CD8_C5-STMN1"))
CD8_object$cd4_cd8_group
sce_CDS_cd8 <- ks_run_Monocle2(object = CD8_object,
                           layer = 'counts',
                           assay = "RNA",
                           VARgenesM="seurat",
                           cellAnno = "cd4_cd8_group",
                           root_state="CD8_C5-STMN1")

saveRDS(sce_CDS_cd8,"D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/T_cells_CD8拟时序分析cds数据.rds")

plot_cell_trajectory(sce_CDS_cd8,color_by = "Pseudotime")
plot_cell_trajectory(sce_CDS_cd8,color_by = "State")

plot_cell_trajectory(sce_CDS_cd8,color_by = "cd4_cd8_group")


sce_CDS_cd8

####美化cd8_plot_cell_trajectory----
colour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",  
         "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0",
         "#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE",
         "#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")
p1 <- plot_cell_trajectory(sce_CDS_cd8, x = 1, y = 2, color_by = "cd4_cd8_group") + 
  theme(legend.position='none',panel.border = element_blank()) + #去掉第一个的legend
  scale_color_manual(values = colour) 
p2 <- plot_complex_cell_trajectory(sce_CDS_cd8, x = 1, y = 2,
                                   color_by = "cd4_cd8_group")+
  scale_color_manual(values = colour) +
  theme(legend.title = element_blank()) 

p3<- p1 / p2

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/拟时序分析/T细胞拟时序分析/CD8_拟时序/cd8拟时序加树.pdf",p3,width = 8,height=9)

####CytoTRACE2分析----
devtools::install_github("guldenemrah/CytoTRACE")
library(CytoTRACE2)

cytotrace2_result_sce <- cytotrace2(CD8_object, 
                                    is_seurat = TRUE, 
                                    slot_type = "counts", 
                                    species = 'human',
                                    seed = 1234)
cytotrace2_result_sce

saveRDS(cytotrace2_result_sce,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/拟时序分析/T细胞拟时序分析/CD8_拟时序/cytotrace2/cytotrace2_result_sce.rds")

# making an annotation dataframe that matches input requirements for plotData function
annotation <- data.frame(phenotype =CD8_object@meta.data$cd4_cd8_group) %>% 
  set_rownames(., colnames(CD8_object))

# plotting
plots <- plotData(cytotrace2_result = cytotrace2_result_sce, 
                  annotation = annotation, 
                  is_seurat = TRUE)
# 绘制CytoTRACE2_Potency的umap图
p1 <- plots$CytoTRACE2_UMAP
# 绘制CytoTRACE2_Potency的umap图
p2 <- plots$CytoTRACE2_Potency_UMAP
# 绘制CytoTRACE2_Relative的umap图 ，v1 
p3 <- plots$CytoTRACE2_Relative_UMAP 
# 绘制各细胞类型CytoTRACE2_Score的箱线图
p4 <- plots$CytoTRACE2_Boxplot_byPheno
ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/拟时序分析/T细胞拟时序分析/CD8_拟时序/cytotrace2/cytotrace2_CytoTRACE2_Boxplot_byPheno.pdf",p4,width =7,height=8)
(p1+p2+p3+p4) + plot_layout(ncol = 2)

#添加cytotrace结果到cds对象中
CytoTRACE2_Score <- cytotrace2_result_sce$CytoTRACE2_Relative %>% as.data.frame()
rownames(CytoTRACE2_Score)
names(CytoTRACE2_Score)<-"Score"
sce_CDS_cd8$cytotrace=CytoTRACE2_Score[match(sce_CDS_cd8$cell_id,rownames(CytoTRACE2_Score)),"Score"]

sce_CDS_cd8$Pseudotime
plot_cell_trajectory(sce_CDS_cd8,color_by = "cytotrace")
cols <- rev(RColorBrewer::brewer.pal(11, "Spectral"))
pdf("cytotrace.pdf",6,6)
plot_cell_trajectory(sce_CDS_cd8, color_by = "cytotrace") +
  scale_colour_gradientn(colours = cols)
dev.off()



##################################################


plot_cell_trajectory(sce_CDS_cd8,color_by = "MPR_Response2")
sce_CDS_cd8$MPR_Response2
p1 = plot_cell_trajectory(sce_CDS_cd8, color_by = "State")
p2 = plot_cell_trajectory(sce_CDS_cd8, color_by = "State") + facet_wrap(~State, nrow = 1)
p1|p2

plot_cell_trajectory(sce_CDS_cd8,color_by = "Herbaspirillum",grouping = "State")
plot_cell_trajectory(sce_CDS_cd8,color_by = "Massilia",grouping = "State")
sce_CDS_cd8
s.genes <- c("Fusobacterium","Massilia","Herbaspirillum")

s.genes <- c("ITGB1","CCR7","KLRB1","GNLY")
# 点图（抖动）
p1 <- plot_genes_jitter(sce_CDS_cd8[s.genes,], grouping = "State", color_by = "State")
# 小提琴图
p2 <- plot_genes_violin(sce_CDS_cd8[s.genes,], grouping = "State", color_by = "State")
# 伪时间图
p3 <- plot_genes_in_pseudotime(sce_CDS_cd8[s.genes,], color_by = "State")
plotc <- p1|p2|p3




####RNAvelocity速率分析----
##https://lishensuo.github.io/posts/bioinfo/024%E5%8D%95%E7%BB%86%E8%83%9E%E5%88%86%E6%9E%90%E5%B7%A5%E5%85%B7--rnavelocity%E9%80%9F%E7%8E%87%E5%88%86%E6%9E%90/

#tree -L 1 ./201008_D0_scRNA_fastq/outs/
  # ./201008_D0_scRNA_fastq/outs/
  # ├── cellsorted_possorted_genome_bam.bam
  # ├── filtered_feature_bc_matrix
  # ├── filtered_feature_bc_matrix.h5
  # ├── filtered_feature_bc_matrix.tar.gz
  # ├── metrics_summary.csv
  # ├── molecule_info.h5
  # ├── possorted_genome_bam.bam
  # ├── possorted_genome_bam.bam.bai
  # ├── raw_feature_bc_matrix
  # ├── raw_feature_bc_matrix.h5
# └── web_summary.html

# (2) 准备三份数据
## 10X比对结果文件夹 
# cellranger_outDir=/path/to/201008_D8Ven_scRNA_fastq
# ## 10X参考数据集中的gtf文件
# cellranger_gtf=/path/to/refdata-gex-GRCh38-2020-A/genes/genes.gtf
# ## 重复序列注释文件(optiona;)
# rmsk_gtf=/path/to/hg38_repeat_rmsk.gtf

#rmsk_gtf下载链接：https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=611454127_NtvlaW6xBSIRYJEBI0iRDEWisITa&clade=mammal&org=Human&db=0&hgta_group=allTracks&hgta_track=rmsk&hgta_table=rmsk&hgta_regionType=genome&position=&hgta_outputType=gff&hgta_outFileName=mm10_rmsk.gtf


#####################################
pdf(file="af_cluster.trajectory.pdf",width=6.5,height=6)
plot_cell_trajectory(cds, color_by = "cd4_cd8_group")+
  scale_color_gradientn(colours = colorRampPalette(c('blue','green','yellow','red'))(100))
dev.off()
pdf(file="D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/拟时序分析/T细胞拟时序分析/CD8_拟时序/af_cellType.trajectory.pdf",width=6.5,height=6)
plot_cell_trajectory(cds,color_by = "cd4_cd8_group")
dev.off()
pdf(file="D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/拟时序分析/T细胞拟时序分析/CD8_拟时序/af_trajectory.Pseudotime.pdf",width=6.5,height=6)
plot_cell_trajectory(cds,color_by = "Pseudotime")
dev.off()
#State轨迹分布图
pdf(file="D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/拟时序分析/T细胞拟时序分析/CD8_拟时序/af_trajectory.state.pdf",width=6.5,height=6)
p1 = plot_cell_trajectory(cds, color_by = "State")
#分面展示
p2 = plot_cell_trajectory(cds, color_by = "State") + facet_wrap(~State, nrow = 1)
p1|p2
dev.off()
cds

plot_cell_trajectory(sce_CDS,show_branch_points = F,color_by = "cluster",cell_size = 3)+
  theme_dr(arrow = grid::arrow(length = unit(0, "inches")))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  scale_color_manual(values = c("#E69253", "#EDB931", "#E4502E", "#4378A0"))+
  geom_curve(aes(x=9,y=2,yend=0,xend=5),
             arrow=arrow(length=unit(0.03,"npc")),
             size=1,curvature=0,
             color="black")+
  annotate("text", x=6.5, y=1.5, label = "bold(YSMP)", parse = TRUE, angle=45, size=4)+
  geom_curve(aes(x=-4,y=-1,yend=3,xend=-6),
             arrow=arrow(length=unit(0.03,"npc")),
             size=1,curvature=0,
             color="black")+
  annotate("text", x=-4, y=1, label = "bold(Monocyte~Fate)", parse = TRUE, angle=-75, size=4)+
  geom_curve(aes(x=-5,y=-3,yend=-5,xend=-6),
             arrow=arrow(length=unit(0.03,"npc")),
             size=1,curvature=0,
             color="black")+
  annotate("text", x=-6.5, y=-4, label = "bold(Neu~Fate)", parse = TRUE, angle=75, size=4)

#结合可视化cytroTRACE2结果
plot_cell_trajectory(sce_CDS,show_branch_points = F,color_by = "CytoTRACE2_Relative",cell_size = 3)+
  scale_colour_gradientn(colours = (c( "#000004FF", "#3B0F70FF", "#8C2981FF", "#DE4968FF", "#FE9F6DFF", "#FCFDBFFF")),
                         na.value = "transparent",
                         limits=c(0,1),
                         breaks = c(0,1),
                         labels=c("(More diff.)", "(Less diff.)"),
                         name = "Relative\norder \n" ,
                         guide = guide_colorbar(frame.colour = "black", ticks.colour = "black"))+
  theme_dr(arrow = grid::arrow(length = unit(0, "inches")))+#坐标轴主题修改
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

####输出拟时序差异基因
groups=subset(pData(cds),select='State')
pbmc=AddMetaData(object=af, metadata=groups, col.name="group")
geneList=list()
for(i in levels(factor(groups$State))){
  pbmc.markers=FindMarkers(pbmc, ident.1 = i, group.by = 'group')
  sig.markers=pbmc.markers[(abs(as.numeric(as.vector(pbmc.markers$avg_log2FC)))>logFCfilter & as.numeric(as.vector(pbmc.markers$p_val_adj))<adjPvalFilter),]
  sig.markers=cbind(Gene=row.names(sig.markers), sig.markers)
  write.table(sig.markers,file=paste0("05.monocleDiff.", i, ".txt"),sep="\t",row.names=F,quote=F)
  geneList[[i]]=row.names(sig.markers)
}
#保存交集基因
unionGenes=Reduce(union,geneList)
write.table(file="monocleDiff.union.txt",unionGenes,sep="\t",quote=F,col.names=F,row.names=F)

####得到这些基因后我们可以随便查看一个基因在时序中的位置
pdf(file = "HSPA1B_af_cds_geneselect.pdf",width =8,height = 7)
pData(cds)[,'HSPA1B'] = scRNA_tpm@assays$RNA@scale.data['HSPA1B',]
plot_cell_trajectory(cds, color_by = 'HSPA1B') + scale_color_gsea()
dev.off()


####拟时相关基因聚类热图
####monocle包的differentialGeneTest()函数可以按条件进行差异分析，将相关参数设为fullModelFormulaStr = "~sm.ns(Pseudotime)"时，可以找到与拟时先关的差异基因。
####为了提高效率，可直接使用cluster差异基因或者高变基因，进行筛选，再绘制热图。以高变基因为例：
disp_table <- dispersionTable(cds)
disp.genes <- subset(disp_table, mean_expression >= 0.5&dispersion_empirical >= 1*dispersion_fit)
disp.genes <- as.character(disp.genes$gene_id)
diff_test <- differentialGeneTest(mycds[disp.genes,], cores = 4, 
                                  fullModelFormulaStr = "~sm.ns(Pseudotime)")
#找到与拟时先关的差异基因
sig_gene_names <- row.names(subset(diff_test, qval < 1e-04))
#以qval为指标，挑选显著性
plot_pseudotime_heatmap(cds[sig_gene_names,], num_clusters=5,
                        show_rownames=T, return_heatmap=T)

####单细胞轨迹通常包含分支（如本例），主要是因为细胞的表达模式受各种因素影响而发生变化；
####BEAM（Branched expression analysis modeling）统计方法可用于寻找以依赖分支方式来实现调控的基因
disp_table <- dispersionTable(mycds)
disp.genes <- subset(disp_table, mean_expression >= 0.5&dispersion_empirical >= 1*dispersion_fit)
disp.genes <- as.character(disp.genes$gene_id)
cds_sub <- cds[disp.genes,]
plot_cell_trajectory(mycds_sub, color_by = "State")
beam_res <- BEAM(mycds_sub, branch_point = 1, cores = 8)
beam_res <- beam_res[order(beam_res$qval),]
beam_res <- beam_res[,c("gene_short_name", "pval", "qval")]
mycds_sub_beam <- mycds_sub[row.names(subset(beam_res, qval < 1e-4)),]
plot_genes_branched_heatmap(mycds_sub_beam,  branch_point = 1, num_clusters = 3, show_rownames = T)

saveRDS(cds, file="D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/T_cells_object拟时序分析cds数据.rds")

#拟时差异基因
expressed_genes=row.names(subset(fData(sce_CDS),use_for_ordering=="TRUE"))
peu_gene <- differentialGeneTest(sce_CDS[expressed_genes,],fullModelFormulaStr = "~sm.ns(Pseudotime)",cores = 1)
# write.csv(peu_gene,file='peu_gene.csv')#保存好文件，这个分析过程挺费时间
peu_gene <- peu_gene[which(peu_gene$qval<0.01 & peu_gene$num_cells_expressed>100),]#筛选显著是的
# peu_gene %>% arrange(qval)  -> peu_gene#按照qval排个序
# peu_gene <- peu_gene[1:100,] #这里我们取前100个基因演示
#改造热图函数，展示需要的基因
source('./new_heatmap.R')
library(pheatmap)
library(grid)
p1 <-plot_pseudotime_heatmap(sce_CDS[peu_gene$gene_short_name,],
                             num_clusters = 4,
                             cores = 2,
                             show_rownames = T,return_heatmap =T,
                             hmcols = viridis(256),
                             use_gene_short_name = T,
                             clustercolor  = c("#d2981a", "#a53e1f", "#457277", "#8f657d"))


p1

#只展示感兴趣的基因
genes <- c("C1QB","C1QC","C1QA","MRC1","LGMN","MS4A7","MAF","FOLR2",
           "HLA-DPA1","CLEC10A","IL10RA","CD163","KCTD12","CLEC7A","MS4A6A","CD14",
           "ITM2A","CYTL1","MDK","SELP","CD24",
           "S100A8","S100A9","S100A12")

source('./add.flag.R')
add.flag(p1, kept.labels = genes, repel.degree = 0.2)


