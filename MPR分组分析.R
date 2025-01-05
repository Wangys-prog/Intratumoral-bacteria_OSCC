#########################################################
#读取数据
library(Seurat)
library(tidyverse)
library(harmony)
library(dplyr)
library(patchwork)
library(SeuratWrappers)
####癌旁数据
object_24_all<- readRDS("D:/wangys_uestc_data2/OSCC_data/24个样本整合分析报告/24oscc_scRNA_data.rds")
object_24_all = UpdateSeuratObject(object_24_all)
object_24_all

View(object_24_all)
Idents(object_24_all)<-"orig.ident"
seurat_sub_sp <- subset(object_24_all, idents = c("LXD_SP_pre",
                                                  "LXD_SP_post",
                                                  "QHH_SP_pre",
                                                  "HTX_SP_pre",
                                                  "XL_SP_pre",
                                                  "XL_SP_post",
                                                  "YXJ_SP_pre",
                                                  "LM_SP_pre",
                                                  "XYH_SP_pre",
                                                  "YXC_SP_pre"))

View(seurat_sub_sp)

library(dplyr)

names(seurat_sub_sp@meta.data)

seurat_sub_sp@meta.data <- seurat_sub_sp@meta.data %>% select(orig.ident,
                                                              nCount_RNA,
                                                              nFeature_RNA,
                                                              rawbc,
                                                              sampleid,
                                                              percent.mito,
                                                              percent.HB,
                                                              log10GenesPerUMI,
                                                              is_metric_outlier,
                                                              is_valid,
                                                              is_doublets,
                                                              RNA_snn_res.0.4)

saveRDS(seurat_sub_sp,"D:/wangys_uestc_data2/OSCC_data/不同治疗癌旁数据-公司分析/rds数据/seurat_sub_sp.rds")

####################################################

setwd("./OSCC_data/MPR_NMPR_group")
OSCC_object_anno<- readRDS("../OSCC_object_anno.rds")
OSCC_object_anno = UpdateSeuratObject(OSCC_object_anno)
OSCC_object_anno
View(OSCC_object_anno)


MPRNMPR.SC.subset<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/MPRNMPR_object_miMPRobe_remove.rds")
MPRNMPR.SC.subset = UpdateSeuratObject(MPRNMPR.SC.subset)
MPRNMPR.SC.subset

MPRNMPR_object_miMPRobe<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/MPRNMPR_object_miMPRobe_v1_harmony.rds")
MPRNMPR_object_miMPRobe = UpdateSeuratObject(MPRNMPR_object_miMPRobe)
MPRNMPR_object_miMPRobe
View(MPRNMPR_object_miMPRobe)

MPRNMPR_object_miMPRobe_remove<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/MPRNMPR_object_miMPRobe_remove.rds")
MPRNMPR_object_miMPRobe_remove = UpdateSeuratObject(MPRNMPR_object_miMPRobe_remove)
MPRNMPR_object_miMPRobe_remove
View(MPRNMPR_object_miMPRobe_remove)
saveRDS(MPRNMPR_object_miMPRobe_remove, file = "D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/MPRNMPR_object_miMPRobe_remove.rds")

##########注释画图#########################

cell_type_new_all<-DimPlot(OSCC_object_anno, reduction = 'umap.harmony', group.by = 'cell_type_new',
                           label = TRUE,raster=FALSE) + NoLegend()

ggsave("细胞注释所有.pdf",cell_type_new_all,width =9, height =7)

#########分组注释##########

library(tidyverse)
library(dplyr)

###读取分组信息
metadata_res<-data.table::fread("MPR_NMPR.csv",header = TRUE)  
metadata<- FetchData(OSCC_object_anno,"orig.ident")
metadata$cell_id <- rownames(metadata)
metadata<- left_join(x=metadata,y=metadata_res,by="orig.ident")
rownames(metadata)<-metadata$cell_id
OSCC_object_anno <- AddMetaData(OSCC_object_anno,metadata = metadata)
table(OSCC_object_anno@meta.data$MPR_MPR_Response2)

metadata_tp<-data.table::fread("SC_SP.csv",header = TRUE)
metadata<- FetchData(OSCC_object_anno,"orig.ident")
metadata$cell_id <- rownames(metadata)
metadata<- left_join(x=metadata,y=metadata_tp,by="orig.ident")
rownames(metadata)<-metadata$cell_id
OSCC_object_anno <- AddMetaData(OSCC_object_anno,metadata = metadata)
table(OSCC_object_anno@meta.data$TPgroup4)


###umap效果较好###
p3 <- DimPlot(OSCC_object_anno, reduction = "umap.harmony",group.by = "MPR_Response2",label = F,raster=FALSE)
p4 <- DimPlot(OSCC_object_anno, reduction = "tsne.harmony",group.by = "MPR_Response2", label =T,raster=FALSE)

saveRDS(OSCC_object_anno, file = "OSCC_object_anno.rds")

##########################################################################################

MPRNMPR.SC.subset<-subset(x =OSCC_object_anno, subset = (orig.ident == "XL_SC_pre"|
                                                         orig.ident == "YMS_SC_pre"|
                                                         orig.ident == "WRY_SC_pre"|
                                                         orig.ident == "HTX_SC_pre"|
                                                         orig.ident == "LM_SC_pre"|
                                                         orig.ident == "QHH_SC_pre"|
                                                         orig.ident == "LXD_SC_pre"|
                                                         orig.ident == "YXJ_SC_pre"|
                                                         orig.ident == "YXC_SC_pre"|
                                                         orig.ident == "XYH_SC_pre"|
                                                         orig.ident == "XL_SC_post"|
                                                         orig.ident == "YMS_SC_post"|
                                                         orig.ident == "LXD_SC_post"|
                                                         orig.ident == "XYH_SC_post"))

table(MPRNMPR.SC.subset$Patients)
table(MPRNMPR.SC.subset$MPR_Response)
saveRDS(MPRNMPR.SC.subset, file = "MPRNMPR.SC.subset_v1_harmony.rds")
MPRNMPR.SC.subset<- readRDS("./MPRNMPR.SC.subset_v1_harmony.rds")
MPRNMPR.SC.subset = UpdateSeuratObject(MPRNMPR.SC.subset)
MPRNMPR.SC.subset

View(MPRNMPR.SC.subset)
MPRNMPR_cell_type<-DimPlot(MPRNMPR.SC.subset, reduction = 'umap.harmony', group.by = 'cell_type_new',
                          label = T,raster=FALSE)+ NoLegend()

MPRNMPR_Patients<-DimPlot(MPRNMPR.SC.subset, reduction = 'umap.harmony', group.by = 'Patients',
                         label = F,raster=FALSE)
MPRNMPR_Treatments<-DimPlot(MPRNMPR.SC.subset, reduction = 'umap.harmony', group.by = 'Treatments',
                           label = F,raster=FALSE)

MPRNMPR_MPR_Response<-DimPlot(MPRNMPR.SC.subset, reduction = 'umap.harmony', group.by = 'MPR_Response',
                         label = F,raster=FALSE)

save(MPRNMPR.SC.subset,MPRNMPR_cell_type,MPRNMPR_Patients,
     MPRNMPR_MPR_Response,
     file = "MPRNMPR.SC.subset_注释.RData")
load(file = "MPRNMPR.SC.subset_注释.RData")
plot_grid(plot.mpg, plot.diamonds, labels = c("A", "B"))

###提取降维信息用ggplot 画图####

library(tidyverse)
MPRNMPR.umap.harmonyredu<- MPRNMPR.SC.subset@reductions$umap.harmony@cell.embeddings%>%
  as.data.frame() %>% 
  cbind(cell_type = MPRNMPR.SC.subset@meta.data$cell_type_new)%>% 
  cbind(Patients = MPRNMPR.SC.subset@meta.data$Patients) %>% 
  cbind(Treatments = MPRNMPR.SC.subset@meta.data$Treatments) %>% 
  cbind(MPR_Response = MPRNMPR.SC.subset@meta.data$MPR_Response)

names(MPRNMPR.umap.harmonyredu)

###颜色设置###画小坐标轴####https://cloud.tencent.com/developer/article/1924260

cell_type_color<- c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080",
                    "#1E90FF","#7CFC00","#FFFF00", "#808000","#FF00FF","#FA8072","#7B68EE",
                    "#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0",
                    "#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4",
                    "#32CD32","#F0E68C","#FFFFE0","#EE82EE","#FF6347","#6A5ACD","#9932CC","#8B008B"
                    ,"#8B4513","#DEB887")

cell_type_color <- c("#FF6347","#ffc556","#1E90FF","#0000FF","#9370DB","#00FFFF","#32cd32","#006400","#FF00ff","#20B2AA")

MPRNMPR.umap.harmonyredu$cell_type<- factor(MPRNMPR.umap.harmonyredu$cell_type, levels = c('T cells', 'B cells','Plasma cells', 'Myeloid cells',  
                                                       'Mast cells',"Endothelial cells","Fibroblasts","Smooth muscle cells",
                                                       "Pericytes","Epithelial cells"))


# library(RColorBrewer)
# library(ggplot2)
# par(mar=c(3,4,2,2))
# display.brewer.all()
# br_pal <- brewer.pal(10,"Paired")



MPRNMPR_cell_type <- ggplot(MPRNMPR.umap.harmonyredu,aes(x= umapharmony_1, y = umapharmony_2,color = cell_type)) +  
  geom_point(size = 0.5, alpha =0.5)  +  
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
        legend.key.size=unit(1,'cm') ) +  # 设置legend标签之间的大小
  guides(color = guide_legend(override.aes = list(size=5)))+ #设置legend中 点的大小 
  geom_segment(aes(x = min(MPRNMPR.umap.harmonyredu$umapharmony_1) , y = min(MPRNMPR.umap.harmonyredu$umapharmony_2) ,
                   xend = min(MPRNMPR.umap.harmonyredu$umapharmony_1) +3, yend = min(MPRNMPR.umap.harmonyredu$umapharmony_2) ),
               colour = "black", size=1,arrow = arrow(length = unit(0.5,"cm")))+ 
  geom_segment(aes(x = min(MPRNMPR.umap.harmonyredu$umapharmony_1)  , y = min(MPRNMPR.umap.harmonyredu$umapharmony_2)  ,
                   xend = min(MPRNMPR.umap.harmonyredu$umapharmony_1) , yend = min(MPRNMPR.umap.harmonyredu$umapharmony_2) + 3),
               colour = "black", size=1,arrow = arrow(length = unit(0.5,"cm"))) +
  annotate("text", x = min(MPRNMPR.umap.harmonyredu$umapharmony_1) +1.5, y = min(MPRNMPR.umap.harmonyredu$umapharmony_2) -1, label = "UMAP_1",
           color="black",size = 7, fontface="bold" ) + 
  annotate("text", x = min(MPRNMPR.umap.harmonyredu$umapharmony_1) -1, y = min(MPRNMPR.umap.harmonyredu$umapharmony_2) + 1.5, label = "UMAP_2",
           color="black",size = 7, fontface="bold" ,angle=90)

MPRNMPR_cell_type

MPRNMPR_cell_type_med <- MPRNMPR.umap.harmonyredu %>%
  group_by(cell_type) %>%
  summarise(
    umapharmony_1 = median(umapharmony_1),
    umapharmony_2 = median(umapharmony_2)
  )

library(ggrepel)
MPRNMPR_cell_type_med

MPRNMPR_cell_type2<-MPRNMPR_cell_type +geom_label_repel(aes(label=cell_type),size=7,color="black",fontface="bold",data = MPRNMPR_cell_type_med,
                                                      point.padding=unit(0.5, "lines"),fill = alpha(c("white"),0.5),
                                                      segment.size=0.5, nudge_x=0.5, nudge_y = 0,direction="y",max.overlaps=50)+
  theme(legend.position = "none")


ggsave('D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞注释_dotplot_featureplot/第二次umap图/MPRNMPR_cell_type2大类umap.jpg',MPRNMPR_cell_type2,width =10,height = 10)



#####Patients

Patients_color <- c("#40e0d0","#ee82ee","#7ccd7c","#551a8b","#ffc0cb","#ff6347",
                    "#f4a460","#b2dfee","#00ff7f","#63b8ff")
##调整legend 顺序

MPRNMPR.umap.harmonyredu$Patients <- factor(MPRNMPR.umap.harmonyredu$Patients, levels=c('P1', 'P2', 'P3', 'P4', 
                                                                        'P5',"P6","P7","P8","P9","P10"))

MPRNMPR_Patients <- ggplot(MPRNMPR.umap.harmonyredu,aes(x= umapharmony_1, y = umapharmony_2,color = Patients)) +  
  geom_point(size = 0.2, alpha =1)  +  
  scale_color_manual(values = Patients_color)+
  xlab("UMAP1")+
  ylab("UMAP2")+
  theme(legend.position = "bottom")+
  guides(color = guide_legend(nrow = 4))+
  theme(panel.grid.major = element_blank(), #主网格线
        panel.grid.minor = element_blank(), #次网格线
        #panel.border = element_blank(), #边框
        axis.title = element_blank(),  #轴标题
        axis.text = element_blank(), # 文本
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = 'white', color="black",linewidth=1), #背景色
        plot.background=element_rect(fill="white"))+
  theme(legend.title = element_blank(), #去掉legend.title 
        legend.key=element_rect(fill='white'), #
        legend.text = element_text(size=20), #设置legend标签的大小
        legend.key.size=unit(1,'cm') ) +  # 设置legend标签之间的大小
  guides(color = guide_legend(override.aes = list(size=5)))+ #设置legend中 点的大小 
  geom_segment(aes(x = min(MPRNMPR.umap.harmonyredu$umapharmony_1) , y = min(MPRNMPR.umap.harmonyredu$umapharmony_2) ,
                   xend = min(MPRNMPR.umap.harmonyredu$umapharmony_1) +3, yend = min(MPRNMPR.umap.harmonyredu$umapharmony_2) ),
               colour = "black", size=0.5,arrow = arrow(length = unit(0.3,"cm")))+ 
  geom_segment(aes(x = min(MPRNMPR.umap.harmonyredu$umapharmony_1)  , y = min(MPRNMPR.umap.harmonyredu$umapharmony_2)  ,
                   xend = min(MPRNMPR.umap.harmonyredu$umapharmony_1) , yend = min(MPRNMPR.umap.harmonyredu$umapharmony_2) + 3),
               colour = "black", size=0.5,arrow = arrow(length = unit(0.3,"cm"))) +
  annotate("text", x = min(MPRNMPR.umap.harmonyredu$umapharmony_1) +1.5, y = min(MPRNMPR.umap.harmonyredu$umapharmony_2) -1, label = "UMAP_1",
           color="black",size = 3, fontface="bold" ) + 
  annotate("text", x = min(MPRNMPR.umap.harmonyredu$umapharmony_1) -1, y = min(MPRNMPR.umap.harmonyredu$umapharmony_2) + 1.5, label = "UMAP_2",
           color="black",size = 3, fontface="bold" ,angle=90)

MPRNMPR_Patients

ggsave('D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞注释_dotplot_featureplot/第二次umap图/MPRNMPR_cell_type2大类病人分组.jpg',MPRNMPR_Patients,width =5,height = 5)

Treatments_color <- c("#32cd32","#b9d3ee")

MPRNMPR_Treatments <- ggplot(MPRNMPR.umap.harmonyredu,aes(x= umapharmony_1, y = umapharmony_2,color = Treatments)) +  
  geom_point(size = 0.2, alpha =1)  +  
  scale_color_manual(values = Treatments_color)+
  #xlab("UMAP1")+
  #ylab("UMAP2")+
  theme(legend.position = "bottom")+
  guides(color = guide_legend(nrow = 1))+
  theme(panel.grid.major = element_blank(), #主网格线
        panel.grid.minor = element_blank(), #次网格线
        #panel.border = element_blank(), #边框
        axis.title = element_blank(),  #轴标题
        axis.text = element_blank(), # 文本
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = 'white', color="black",linewidth=1), #背景色
        plot.background=element_rect(fill="white"))+
  theme(legend.title = element_blank(), #去掉legend.title 
        legend.key=element_rect(fill='white'), #
        legend.text = element_text(size=20), #设置legend标签的大小
        legend.key.size=unit(1,'cm') ) +  # 设置legend标签之间的大小
  guides(color = guide_legend(override.aes = list(size=5)))+ #设置legend中 点的大小 
  geom_segment(aes(x = min(MPRNMPR.umap.harmonyredu$umapharmony_1) , y = min(MPRNMPR.umap.harmonyredu$umapharmony_2) ,
                   xend = min(MPRNMPR.umap.harmonyredu$umapharmony_1) +3, yend = min(MPRNMPR.umap.harmonyredu$umapharmony_2) ),
               colour = "black", size=0.5,arrow = arrow(length = unit(0.3,"cm")))+ 
  geom_segment(aes(x = min(MPRNMPR.umap.harmonyredu$umapharmony_1)  , y = min(MPRNMPR.umap.harmonyredu$umapharmony_2)  ,
                   xend = min(MPRNMPR.umap.harmonyredu$umapharmony_1) , yend = min(MPRNMPR.umap.harmonyredu$umapharmony_2) + 3),
               colour = "black", size=0.5,arrow = arrow(length = unit(0.3,"cm"))) +
  annotate("text", x = min(MPRNMPR.umap.harmonyredu$umapharmony_1) +1.5, y = min(MPRNMPR.umap.harmonyredu$umapharmony_2) -1, label = "UMAP_1",
           color="black",size = 3, fontface="bold" ) + 
  annotate("text", x = min(MPRNMPR.umap.harmonyredu$umapharmony_1) -1, y = min(MPRNMPR.umap.harmonyredu$umapharmony_2) + 1.5, label = "UMAP_2",
           color="black",size = 3, fontface="bold" ,angle=90)

MPRNMPR_Treatments

ggsave('D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞注释_dotplot_featureplot/第二次umap图/MPRNMPR_cell_type2大类治疗分组.jpg',MPRNMPR_Treatments,width =5,height = 5)

MPR_Response_color <- c("#ff8247","#a020f0")


MPRNMPR_MPR_Response <- ggplot(MPRNMPR.umap.harmonyredu,aes(x= umapharmony_1, y = umapharmony_2,color = MPR_Response)) +  
  geom_point(size = 0.2, alpha =1)  +  
  scale_color_manual(values = MPR_Response_color)+
  #xlab("UMAP1")+
  #ylab("UMAP2")+
  theme(legend.position = "bottom")+
  guides(color = guide_legend(nrow = 1))+
  theme(panel.grid.major = element_blank(), #主网格线
        panel.grid.minor = element_blank(), #次网格线
        #panel.border = element_blank(), #边框
        axis.title = element_blank(),  #轴标题
        axis.text = element_blank(), # 文本
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = 'white', color="black",linewidth=1), #背景色
        plot.background=element_rect(fill="white"))+
  theme(legend.title = element_blank(), #去掉legend.title 
        legend.key=element_rect(fill='white'), #
        legend.text = element_text(size=20), #设置legend标签的大小
        legend.key.size=unit(1,'cm') ) +  # 设置legend标签之间的大小
  guides(color = guide_legend(override.aes = list(size=5)))+ #设置legend中 点的大小 
  geom_segment(aes(x = min(MPRNMPR.umap.harmonyredu$umapharmony_1) , y = min(MPRNMPR.umap.harmonyredu$umapharmony_2) ,
                   xend = min(MPRNMPR.umap.harmonyredu$umapharmony_1) +3, yend = min(MPRNMPR.umap.harmonyredu$umapharmony_2) ),
               colour = "black", size=0.5,arrow = arrow(length = unit(0.3,"cm")))+ 
  geom_segment(aes(x = min(MPRNMPR.umap.harmonyredu$umapharmony_1)  , y = min(MPRNMPR.umap.harmonyredu$umapharmony_2)  ,
                   xend = min(MPRNMPR.umap.harmonyredu$umapharmony_1) , yend = min(MPRNMPR.umap.harmonyredu$umapharmony_2) + 3),
               colour = "black", size=0.5,arrow = arrow(length = unit(0.3,"cm"))) +
  annotate("text", x = min(MPRNMPR.umap.harmonyredu$umapharmony_1) +1.5, y = min(MPRNMPR.umap.harmonyredu$umapharmony_2) -1, label = "UMAP_1",
           color="black",size = 3, fontface="bold" ) + 
  annotate("text", x = min(MPRNMPR.umap.harmonyredu$umapharmony_1) -1, y = min(MPRNMPR.umap.harmonyredu$umapharmony_2) + 1.5, label = "UMAP_2",
           color="black",size = 3, fontface="bold" ,angle=90)

MPRNMPR_MPR_Response

ggsave('D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞注释_dotplot_featureplot/第二次umap图/MPRNMPR_cell_type2大类反应分组.jpg',MPRNMPR_MPR_Response,width =5,height = 5)

library(cowplot)

bottom_row <- plot_grid(MPRNMPR_Patients,MPRNMPR_Treatments,MPRNMPR_MPR_Response, nrow =1, align = 'h', rel_widths = c(1,1,1),rel_heights = c(1,1,1))

MPRNMPR_cell_type_merged <-plot_grid(MPRNMPR_cell_type2,bottom_row,
                             ncol = 1,
                             #labels = c('A','空图','B'),
                             rel_widths = c(1,1.3),
                             rel_heights = c(1,0.7))

ggsave('./MPRNMPR_cell_type_merged.jpg',MPRNMPR_cell_type_merged,width =15,height = 15)
ggsave('./MPRNMPR_cell_type_merged.pdf',MPRNMPR_cell_type_merged,width =15,height = 15)

save(MPRNMPR.SC.subset,MPRNMPR_cell_type_med, MPRNMPR_cell_type2,MPRNMPR_Patients,
     MPRNMPR_MPR_Response,MPRNMPR_Treatments,MPRNMPR.umap.harmonyredu,
     file = "MPRNMPR.SC.subset_注释.RData")

load(file="MPRNMPR.SC.subset_注释.RData")


saveRDS(MPRNMPR.SC.subset, file = "MPRNMPR.SC.subset_v1_harmony.rds")


####不同细胞在不同分组中的比例###https://zhuanlan.zhihu.com/p/653120609
library(reshape2)
library(ggplot2)
library(dplyr)
library(Seurat)
library(RColorBrewer)
cell_type_color <- c("#FF6347","#ffc556","#1E90FF","#0000FF","#9370DB","#00FFFF","#32cd32","#006400","#FF00ff","#20B2AA")

Cell_number <- table(MPRNMPR.SC.subset@meta.data$cell_type_new,MPRNMPR.SC.subset@meta.data$MPR_Response2) %>% melt()
colnames(Cell_number) <- c("Cluster","Group","Number")
Cell_number$Cluster <- factor(Cell_number$Cluster)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
Cell_number_color <- col_vector[1:10]

Cell_number$Cluster <- factor(Cell_number$Cluster, levels=c('T cells', 'B cells','Plasma cells', 'Myeloid cells',  
                                                        'Mast cells',"Endothelial cells","Fibroblasts","Smooth muscle cells",
                                                        "Pericytes","Epithelial cells"))
View(Cell_number)

ggplot(data = Cell_number, aes(x =Number, y = Group, fill = Cluster)) +
  geom_bar(stat = "identity", width=0.8,position="fill")+
  scale_fill_manual(values = cell_type_color)+
  theme_bw()+
  theme(panel.grid =element_blank()) +
  labs(x="Ratio",y="")+
  ####用来将y轴移动位置
  theme(axis.text.y = element_text(size=12, colour = "black"))+
  theme(axis.text.x = element_text(size=12, colour = "black"))+
  theme(
    axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 45)
  )


Cellratio <- prop.table(table(MPRNMPR.SC.subset@meta.data$cell_type_new,MPRNMPR.SC.subset@meta.data$MPR_Response2), margin = 2)

Cellratio <- as.data.frame(Cellratio)
colnames(Cellratio) <- c("Cluster","Group","Freq")

colourCount = length(unique(Cellratio$Cluster))
Cellratio$Cluster <- factor(Cellratio$Cluster, levels=c('T cells', 'B cells','Plasma cells', 'Myeloid cells',  
                                                            'Mast cells',"Endothelial cells","Fibroblasts","Smooth muscle cells",
                                                            "Pericytes","Epithelial cells"))
ggplot(data = Cellratio, aes(x =Group, y = Freq, fill = Cluster)) +
  geom_bar(stat = "identity", width=0.8,position="fill")+
  scale_fill_manual(values = cell_type_color )+
  theme_bw()+
  theme(panel.grid =element_blank()) +
  labs(x="",y="Ratio")+
  ####用来将y轴移动位置
  theme(axis.text.y = element_text(size=12, colour = "black"))+
  theme(axis.text.x = element_text(size=12, colour = "black"))+
  theme(
    axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 45)
  ) 

###绘制饼图
Cellratio$Group

Pre_NMPR_Cellratio<- Cellratio[Cellratio$Group %in% "Pre_NMPR", ]


Pre_NMPR_bp<- ggplot(Pre_NMPR_Cellratio, aes(x="", y=Freq, fill=Cluster))+
  geom_bar(width = 1, stat = "identity")+
  coord_polar("y", start=0)+
  scale_fill_manual(values = cell_type_color)+
  theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
  )
Pre_NMPR_bp

Pre_MPR_Cellratio<- Cellratio[Cellratio$Group %in% "Pre_MPR", ]


Pre_MPR_bp<- ggplot(Pre_MPR_Cellratio, aes(x="", y=Freq, fill=Cluster))+
  geom_bar(width = 1, stat = "identity")+
  coord_polar("y", start=0)+
  scale_fill_manual(values = cell_type_color)+
  theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
  )
Pre_MPR_bp


Post_NMPR_Cellratio<- Cellratio[Cellratio$Group %in% "Post_NMPR", ]


Pre_MPR_bp<- ggplot(Pre_MPR_Cellratio, aes(x="", y=Freq, fill=Cluster))+
  geom_bar(width = 1, stat = "identity")+
  coord_polar("y", start=0)+
  scale_fill_manual(values = cell_type_color)+
  theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
  )
Pre_MPR_bp


#######细胞比例差异统计图######################
Cellratio <- prop.table(table(MPRNMPR.SC.subset@meta.data$cell_type_new,MPRNMPR.SC.subset@meta.data$orig.ident), margin = 2)

Cellratio <- table(MPRNMPR.SC.subset@meta.data$cell_type_new,MPRNMPR.SC.subset@meta.data$orig.ident)

Cellratio <- as.data.frame(Cellratio)
colnames(Cellratio) <- c("Cluster","Group","Freq")
library(reshape2)
cellper <- dcast(Cellratio,Group~Cluster, value.var = "Freq")#长数据转为宽数据
write.csv(cellper,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞比例scCODA/细胞大类比例.csv")
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


cell_type_groups = c('T cells', 'B cells','Plasma cells', 'Myeloid cells',  
                     'Mast cells',"Endothelial cells","Fibroblasts","Smooth muscle cells",
                     "Pericytes","Epithelial cells")
library(ggplot2)
library(dplyr)
library(ggpubr)
library(cowplot)
treatment_color <- c("#f4a460","#6ca6cd")

cell_type <- "T cells"
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
  
  cellper_cell_type$MPR_Response2 <- factor(cellper_cell_type$MPR_Response2,levels=c("Pre_MPR","Post_MPR","Pre_NMPR","Post_NMPR"))
  
  pp1 = ggplot(cellper_cell_type,aes(x=MPR_Response2,y=percent)) + #ggplot作图
    geom_boxplot(aes(fill=Treatments),outlier.shape = NA,lwd= 1)+
    #geom_boxplot(outlier.colour="red", outlier.shape=7,outlier.size=1) +
    scale_fill_manual(values = treatment_color)+
    #geom_jitter(shape = 21,aes(fill=Treatments$Treatments),width = 0.25) + 
    #stat_summary(fun=mean, geom="point", color="grey60") +
    theme_cowplot() +
    theme(axis.text = element_text(size = 25),
          axis.line.x=element_line(linetype=1,color="black",size=1),
          axis.ticks.x=element_line(color="black",size=2,lineend = 1),
          axis.line.y=element_line(linetype=1,color="black",size=1),
          axis.ticks.y=element_line(color="black",size=2,lineend = 1),
          axis.title = element_text(size = 25),
          legend.text = element_text(size = 25),
          axis.text.x = element_text(vjust=1,hjust=1,angle=30,size=25),
          axis.text.y = element_text(vjust=1,hjust=1,angle=0,size=20),
          legend.title = element_text(size = 25),
          plot.title = element_text(size = 25),
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
                               method="wilcox.test",size=10,bracket.size = 1)
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

# library(cowplot)
cells_<-plot_grid(pplist[['T cells']],
                  pplist[['B cells']],
                  pplist[['Plasma cells']],
                  pplist[['Myeloid cells']],
                  pplist[['Mast cells']],
                  pplist[['Endothelial cells']],
                  pplist[['Fibroblasts']],
                  pplist[['Pericytes']],
                  pplist[['Smooth muscle cells']],
                  pplist[['Epithelial cells']],
                  align = "h",  #axis = 'l',
                  nrow =2)


plt1<- ggplot(cellper_cell_type,aes(x=MPR_Response,y=percent)) + #ggplot作图
  geom_boxplot(aes(fill=Treatments),outlier.shape = NA,lwd= 1)+
  #geom_boxplot(outlier.colour="red", outlier.shape=7,outlier.size=1) +
  scale_fill_manual(values = treatment_color)+
  #geom_jitter(shape = 21,aes(fill=Treatments$Treatments),width = 0.25) + 
  #stat_summary(fun=mean, geom="point", color="grey60") +
  theme_cowplot() +
  theme(axis.text = element_text(size = 25),
        legend.key.size = unit(50, "pt"),
        axis.title = element_text(size = 25,face="bold"),
        legend.text = element_text(size = 20,face="bold"),
        axis.text.x = element_text(vjust=1,hjust=1,angle=30,size=25,face="bold"),
        axis.text.y = element_text(vjust=1,hjust=1,angle=0,size=20,face="bold"),
        legend.title = element_text(size = 20),
        plot.title = element_text(size = 25,face="bold"),
        legend.position = "right") + 
  labs(title = cell_type,y='Percentage',x="") +
  guides(fill = guide_legend(title = NULL))

legend <- get_legend(plt1+
                       #guides(color = guide_legend(nrow = 1)) +
                       theme(legend.position = "right"))

cells_pre2<- plot_grid(cells_,legend,nrow=2,rel_widths=c(1,0.5))

ggsave('MPR_cells_percentage.jpg',cells_,width =30,height =25,limitsize = FALSE)

ggsave('MPR_cells_percentage.pdf',cells_,width =30,height =25,limitsize = FALSE)

save(MPRNMPR.SC.subset,MPRNMPR_cell_type_med, MPRNMPR_cell_type2,MPRNMPR_Patients,
     MPRNMPR_MPR_Response,MPRNMPR_Treatments,MPRNMPR.umap.harmonyredu,cell_type_color,
     Patients_color,MPR_Response_color,Treatments_color,Cellratio,Pre_NMPR,Pre_NMPR_Cellratio,
     Pre_NMPR_bp,
     file = "MPRNMPR.SC.subset_注释.RData")


##第一列样本或者分组信息，第一行物种名称，细胞名称
###整理数据为固定格式数据


####细胞比例九宫格图----
#画图使用数据为长数据，第一列group，第二列物种或者细胞类型，第三列为value
Cellratio <- prop.table(table(MPRNMPR.SC.subset@meta.data$cell_type_new,MPRNMPR.SC.subset@meta.data$MPR_Response2), margin = 2)
Cellratio <- as.data.frame(Cellratio)
View(Cellratio)
colnames(Cellratio) <- c("Cluster","Response","Freq")

write.csv(Cellratio,"Cellratio_细胞比例饼图数据.csv")
###加载
library(tidyverse)
library(ggtext)
library(glue)
# devtools::install_github("AllanCameron/VoronoiPlus")
library(VoronoiPlus) 

au_vor1 <- voronoi_treemap(Freq ~ Cluster,
                           data = Cellratio %>% filter(Response=="Pre_NMPR")) %>% 
  mutate(type="Pre_NMPR")

View(au_vor1)

au_vor2 <- voronoi_treemap(Freq ~ Cluster,
                           data = Cellratio %>% filter(Response=="Pre_MPR")) %>% 
  mutate(type="Pre_MPR")

au_vor3 <- voronoi_treemap(Freq ~ Cluster,
                           data = Cellratio %>% filter(Response=="Post_MPR")) %>% 
  mutate(type="Post_MPR")

groups <- au_vor1 %>% bind_rows(au_vor2,au_vor3)

groups$type<- factor(groups$type, levels = c('Pre_NMPR', 'Pre_MPR','Post_MPR'))
groups$group<- factor(groups$group, levels = c('T cells', 'B cells','Plasma cells', 'Myeloid cells',  
                                                                 'Mast cells',"Endothelial cells","Fibroblasts","Smooth muscle cells",
                                                                 "Pericytes","Epithelial cells"))
cell_type_color <- c("#FF6347","#ffc556","#1E90FF","#0000FF","#9370DB","#00FFFF","#32cd32","#006400","#FF00ff","#20B2AA")

label_group<-groups%>% group_by(group,type,value) %>% summarise(x_median=median(x),y_median=median(y))
nrow(value_group)
nrow(label_group)
View(label_group)

label_group$value <- paste(round(label_group$value *100,2),sep="","%")
pie_plot<- ggplot() +
  geom_polygon(data = groups,mapping = aes(x = x, y = y, group = group, fill = group),
               colour = "white",linewidth =0.8) +
  facet_wrap(.~type,scale="free")+
  scale_fill_manual(values =cell_type_color)+
  geom_text(aes(x_median,y_median,label=value),size=5,color="black",fontface="bold",data = label_group)+
  labs(x=NULL,y=NULL)+
  theme_test()+
  theme(axis.text=element_blank(),
        axis.ticks = element_blank(),
        plot.background = element_rect(fill = "white", colour = "white"),
        panel.background = element_rect(fill = "white", colour = "white"),
        panel.spacing.x=unit(-0.1,"cm"),
        panel.spacing.y=unit(-0.1,"cm"),
        plot.margin = margin(2,2,0,0),
        legend.key.height = unit(0.2,"in"), 
        legend.key.width = unit(0.2,"in"),
        legend.title = element_blank(),
        legend.text = element_text(margin = margin(l=0,unit="cm"),size=12),
        strip.background = element_rect(fill="grey96"),
        strip.text = element_text(face="bold",size=12))

pie_plot


ggsave("MPR_pie_细胞比例2.jpg",pie_plot,width=16.5,height=5)

ggsave("MPR_pie_细胞比例2.pdf",pie_plot,width=16.5,height=5)

write.csv(groups,"细胞比例饼图长数据groups.csv")

####细胞比例scCODA差异统计----
setwd("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group")
Idents(MPRNMPR_object_miMPRobe)<-"orig.ident"
metadata_res<-data.table::fread("MPR_NMPR.csv",header = TRUE)[,c(1,4)]
View(metadata_res)
metadata<- FetchData(MPRNMPR_object_miMPRobe,"orig.ident")
metadata$cell_id <- rownames(metadata)
metadata<- left_join(x=metadata,y=metadata_res,by="orig.ident")
rownames(metadata)<-metadata$cell_id
MPRNMPR_object_miMPRobe <- AddMetaData(MPRNMPR_object_miMPRobe,metadata = metadata)
table(MPRNMPR_object_miMPRobe@meta.data$MPR_Response3)

saveRDS(MPRNMPR_object_miMPRobe, file = "D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/MPRNMPR.SC.subset_v1_harmony.rds")

MPRNMPR_object_miMPRobe<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/MPRNMPR_object_miMPRobe_v1_harmony.rds")
MPRNMPR_object_miMPRobe = UpdateSeuratObject(MPRNMPR_object_miMPRobe)
MPRNMPR_object_miMPRobe

setwd("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞比例scCODA")

Cellratio <- table(MPRNMPR.SC.subset@meta.data$cell_type_new,MPRNMPR.SC.subset@meta.data$orig.ident)

Cellratio <- as.data.frame(Cellratio)
colnames(Cellratio) <- c("Cluster","Group","Freq")
library(reshape2)
cellper <- dcast(Cellratio,Group~Cluster, value.var = "Freq")#长数据转为宽数据
write.csv(cellper,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞比例scCODA/细胞大类比例.csv")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ALDEx2")
library(ALDEx2)

##Seurat 初始点图marker dotplot_featureplot----


devtools::install_github("junjunlab/jjAnno",force = TRUE)
library(jjAnno)

markers <- c("TRAC","CD3E","CD2","CD7", # Tcells
             "MS4A1","CD19","CD79A","BANK1", #B cells
             "IGLC2","MZB1","DERL3","PDK1" ,#Plasma cell
             "C1QB","C1QA","GCA","CD14", # Myeloid cells
             "TPSB2","TPSAB1","IL1RL1","MS4A2", # Mast cells
             "FLT4","PROX1","RASIP1","NOTCH4",#Endothelial cells
             "COL12A1","MFAP5","PDGFRA","FBN2", # Fibroblasts
             "ACTA1","MYL1","MYBPC1","MYH2", #Smooth muscle cells
             "HIGD1B","HEYL","STEAP4","ABCC9",#Pericytes
             "CDH1","CD24","KRT15","KRT13") # Endothelial cells

cluster<- c('T cells', 'B cells','Plasma cells', 'Myeloid cells',  
            'Mast cells',"Endothelial cells","Fibroblasts","Smooth muscle cells",
            "Pericytes","Epithelial cells")

p1 <- DotPlot(MPRNMPR.SC.subset, features = markers ,
              assay='RNA' )

sub_markers_count<- p1$data

names(sub_markers_count)

sub_markers_count$id<- factor(sub_markers_count$id, levels = c('T cells', 'B cells','Plasma cells', 'Myeloid cells',  
                                                               'Mast cells',"Endothelial cells","Fibroblasts","Smooth muscle cells",
                                                               "Pericytes","Epithelial cells"))


cell_type_color <- c("#FF6347","#ffc556","#1E90FF","#0000FF","#9370DB","#00FFFF","#32cd32","#006400","#FF00ff","#20B2AA")
pdot <-
  ggplot(sub_markers_count,aes(x = features.plot,y = id)) +
  geom_point(aes(fill =avg.exp.scaled,size = pct.exp),
             color = 'black',
             shape = 21) +
  theme_bw(base_size = 14) +
  xlab('') + ylab('') +
  scale_fill_gradient2(low = 'white',mid = '#EB1D36',high = '#990000',
                       midpoint = 5,
                       limits = c(0,10),
                       #breaks = seq(0,2,10),
                       #labels = seq(0,2,10),
                       name = 'Average Expression') +
  scale_size(range = c(1,5),
             limits = c(0, 100),
             breaks = seq(20,100,20),
             labels = seq(20,100,20)
  ) +
  theme(panel.grid = element_blank(),
        legend.title = element_text(size=12,face="bold"),
        legend.text=element_text(size=9,face="bold"),
        axis.text = element_text(color = 'black'),
        aspect.ratio = 0.5,
        plot.margin = margin(t = 0.5,r = 0.5,b = 0.5,l = 0.5,unit = 'cm'),
        axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5, face = 'bold.italic',size=15),
        axis.text.y = element_text(color = "black",size=15,face="bold")
  ) +
  coord_cartesian(clip = 'off') +
  guides(
    fill = guide_colorbar(
      direction = "horizontal",
      title.position = "top",
      barwidth = unit(5, "cm")
    ),
    size = guide_legend(
      title = "Percent Expressed",
      direction = "horizontal",
      title.position = "top",
      label.position = "bottom",
      override.aes = list(
        color = "black",
        fill = "grey"
      )
    )
  )

pdot
library(RColorBrewer)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
#处理后有73种差异还比较明显的颜色，基本够用
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))) 


# p1 <- annoSegment(object = pdot,
#                   annoPos = 'top', # 将注释条，添加到顶部
#                   xPosition = c(2.5,6.5,10.5,14.5,18.5,22.5,26.5,30.5,34.5,38.5), # 注释条的x轴的位置
#                   yPosition = c(10.8,12), # 注释条的y轴位置，需要多次个性化调整
#                   segWidth = 3, # 注释条的宽度，需要多次个性化调整
#                   pCol = c(cell_type_color) # 注释条的颜色
# )
# 
# p1

p2 <- annoRect(object = pdot,
               annoPos = 'left', # 添加方框的位置
               annoManual = T, # 可以手动设置
               # 设置方框的位置，由于刻度是整数，从下往上依次是1，2，3，4
               # 所以设置方框的位置是0.5，1.5，等，上下包含要框的内容
               # list 表示添加控件的个数，比如下面表示添加 [0.5,4.5]和[4.5,8.5]
               yPosition = list(c(0.5,1.5),c(1.5,2.5)),
               # 个性化多次调整，多尝试，找一个合适的数据值
               xPosition = c(-12,0.3),
               pCol = rep('white',10), # 边框颜色
               pFill = c("#FF6347","#ffc556"), # 填充颜色
               alpha = 0.2 # 填充颜色透明度
)

p2

p2_1 <- annoRect(object = p2,
                 annoPos = 'left', # 添加方框的位置
                 annoManual = T, # 可以手动设置
                 # 设置方框的位置，由于刻度是整数，从下往上依次是1，2，3，4
                 # 所以设置方框的位置是0.5，1.5，等，上下包含要框的内容
                 # list 表示添加控件的个数，比如下面表示添加 [0.5,4.5]和[4.5,8.5]
                 yPosition = list(c(2.5,3.5),c(3.5,4.5)),
                 # 个性化多次调整，多尝试，找一个合适的数据值
                 xPosition = c(-12,0.3),
                 pCol = rep('white',10), # 边框颜色
                 pFill = c("#1E90FF","#0000FF"), # 填充颜色
                 alpha = 0.2 # 填充颜色透明度
)

p2_1

p2_2 <- annoRect(object = p2_1,
                 annoPos = 'left', # 添加方框的位置
                 annoManual = T, # 可以手动设置
                 # 设置方框的位置，由于刻度是整数，从下往上依次是1，2，3，4
                 # 所以设置方框的位置是0.5，1.5，等，上下包含要框的内容
                 # list 表示添加控件的个数，比如下面表示添加 [0.5,4.5]和[4.5,8.5]
                 yPosition = list(c(4.5,5.5),c(5.5,6.5)),
                 # 个性化多次调整，多尝试，找一个合适的数据值
                 xPosition = c(-12,0.3),
                 pCol = rep('white',10), # 边框颜色
                 pFill = c("#9370DB","#00FFFF"), # 填充颜色
                 alpha = 0.2 # 填充颜色透明度
)


p2_2

p2_3 <- annoRect(object = p2_2,
                 annoPos = 'left', # 添加方框的位置
                 annoManual = T, # 可以手动设置
                 # 设置方框的位置，由于刻度是整数，从下往上依次是1，2，3，4
                 # 所以设置方框的位置是0.5，1.5，等，上下包含要框的内容
                 # list 表示添加控件的个数，比如下面表示添加 [0.5,4.5]和[4.5,8.5]
                 yPosition = list(c(6.5,7.5),c(7.5,8.5)),
                 # 个性化多次调整，多尝试，找一个合适的数据值
                 xPosition = c(-12,0.3),
                 pCol = rep('white',10), # 边框颜色
                 pFill = c("#32cd32","#006400"), # 填充颜色
                 alpha = 0.2 # 填充颜色透明度
)

p2_3

p2_4 <- annoRect(object = p2_3,
                 annoPos = 'left', # 添加方框的位置
                 annoManual = T, # 可以手动设置
                 # 设置方框的位置，由于刻度是整数，从下往上依次是1，2，3，4
                 # 所以设置方框的位置是0.5，1.5，等，上下包含要框的内容
                 # list 表示添加控件的个数，比如下面表示添加 [0.5,4.5]和[4.5,8.5]
                 yPosition = list(c(8.5,9.5),c(9.5,10.5)),
                 # 个性化多次调整，多尝试，找一个合适的数据值
                 xPosition = c(-12,0.3),
                 pCol = rep('white',10), # 边框颜色
                 pFill = c("#FF00ff","#20B2AA"), # 填充颜色
                 alpha = 0.2 # 填充颜色透明度
)

p2_4

p3<- annoSegment(object = p2_4,
                 annoPos = 'top',
                 xPosition = c(2.5,6.5,10.5,14.5,18.5,22.5,26.5,30.5,34.5,38.5),
                 yPosition = c(10.8,12),
                 segWidth = 3,
                 #addBranch = T,
                 #lwd = 4,
                 #branDirection = -1,
                 addText = T,
                 textLabel = c('T cells', 'B cells','Plasma cells', 'Myeloid cells',  
                               'Mast cells',"Endothelial cells","Fibroblasts","Smooth muscle cells",
                               "Pericytes","Epithelial cells"),
                 textRot = 45,
                 hjust = 0,
                 vjust = 0,
                 textCol = c(rep("black", 10)),
                 fontface = "bold",
                 pCol = c(cell_type_color),
                 textSize = 12)
p3

p4_1 <- annoRect(object = p3,
                 annoPos = 'left',
                 annoManual = T,
                 yPosition = c(0.5,10.5),
                 xPosition = c(0.5,4.5),
                 pCol = 'black',
                 pFill = 'transparent',
                 lty = 'dashed',
                 lwd = 1.5)

p4

p4_2 <- annoRect(object = p4_1,
                 annoPos = 'left',
                 annoManual = T,
                 yPosition = c(0.5,10.5),
                 xPosition = c(4.5,8.5),
                 pCol = 'black',
                 pFill = 'transparent',
                 lty = 'dashed',
                 lwd = 1.5)

p4_3 <- annoRect(object = p4_2,
                 annoPos = 'left',
                 annoManual = T,
                 yPosition = c(0.5,10.5),
                 xPosition = c(8.5,12.5),
                 pCol = 'black',
                 pFill = 'transparent',
                 lty = 'dashed',
                 lwd = 1.5)

p4_4 <- annoRect(object = p4_3,
                 annoPos = 'left',
                 annoManual = T,
                 yPosition = c(0.5,10.5),
                 xPosition = c(12.5,16.5),
                 pCol = 'black',
                 pFill = 'transparent',
                 lty = 'dashed',
                 lwd = 1.5)

p4_5 <- annoRect(object = p4_4,
                 annoPos = 'left',
                 annoManual = T,
                 yPosition = c(0.5,10.5),
                 xPosition = c(16.5,20.5),
                 pCol = 'black',
                 pFill = 'transparent',
                 lty = 'dashed',
                 lwd = 1.5)

p4_6 <- annoRect(object = p4_5,
                 annoPos = 'left',
                 annoManual = T,
                 yPosition = c(0.5,10.5),
                 xPosition = c(20.5,24.5),
                 pCol = 'black',
                 pFill = 'transparent',
                 lty = 'dashed',
                 lwd = 1.5)

p4_7 <- annoRect(object = p4_6,
                 annoPos = 'left',
                 annoManual = T,
                 yPosition = c(0.5,10.5),
                 xPosition = c(24.5,28.5),
                 pCol = 'black',
                 pFill = 'transparent',
                 lty = 'dashed',
                 lwd = 1.5)

p4_8 <- annoRect(object = p4_7,
                 annoPos = 'left',
                 annoManual = T,
                 yPosition = c(0.5,10.5),
                 xPosition = c(28.5,32.5),
                 pCol = 'black',
                 pFill = 'transparent',
                 lty = 'dashed',
                 lwd = 1.5)

p4_9 <- annoRect(object = p4_8,
                 annoPos = 'left',
                 annoManual = T,
                 yPosition = c(0.5,10.5),
                 xPosition = c(32.5,36.5),
                 pCol = 'black',
                 pFill = 'transparent',
                 lty = 'dashed',
                 lwd = 1.5)

p4_10 <- annoRect(object = p4_9,
                  annoPos = 'left',
                  annoManual = T,
                  yPosition = c(0.5,10.5),
                  xPosition = c(36.5,40.5),
                  pCol = 'black',
                  pFill = 'transparent',
                  lty = 'dashed',
                  lwd = 1.5)

p4_10

ggsave("MPR_marker_dotplot.jpg",p4_10,width =14,height =9)

ggsave("MPR_marker_dotplot.pdf",p4_10,width =14,height =9)

####marker feature plot####

library(viridis)
library(scCustomize)
marker <- c("TRAC","CD3E","CD2","CD7", # Tcells
            "MS4A1","CD19","CD79A","BANK1", #B cells
            "IGLC2","MZB1","DERL3","PDK1" ,#Plasma cell
            "C1QB","C1QA","GCA","CD14", # Myeloid cells
            "TPSB2","TPSAB1","IL1RL1","MS4A2", # Mast cells
            "FLT4","PROX1","RASIP1","NOTCH4",#Endothelial cells
            "COL12A1","MFAP5","PDGFRA","FBN2", # Fibroblasts
            "ACTA1","MYL1","MYBPC1","MYH2", #Smooth muscle cells
            "HIGD1B","HEYL","STEAP4","ABCC9",#Pericytes
            "CDH1","CD24","KRT15","KRT13") # Endothelial cells

marker_gene_list <- c("TRAC", # Tcells
                      "MS4A1", #B cells
                      "IGLC2",#Plasma cell
                      "C1QB", # Myeloid cells
                      "TPSB2", # Mast cells
                      "PROX1",#Endothelial cells
                      "VEGFC",#Endothelial cells
                      "COL12A1", # Fibroblasts
                      "ACTA1", #Smooth muscle cells
                      "HIGD1B",#Pericytes
                      "KRT15") # Endothelial cells

pal <- viridis(n = 10, option = "C")
# "#0D0887FF" "#47039FFF" "#7301A8FF" "#9C179EFF" "#BD3786FF" "#D8576BFF" "#ED7953FF" "#FA9E3BFF" "#FDC926FF" "#F0F921FF"
pal <- viridis(n = 15, option = "D", direction = -1)
FeaturePlot_scCustom(seurat_object = MPRNMPR.SC.subset, reduction = 'umap.harmony',
                         colors_use = colorRampPalette(c("#3288BD", "white", "#D53E4F" ))(50), 
                         features = "KRT16")

plots=list()
for (i in marker_gene_list){
  plots[[i]]=FeaturePlot_scCustom(seurat_object = MPRNMPR.SC.subset, reduction = 'umap.harmony',
                                  colors_use = colorRampPalette(c("#A9A9A9", "white", "#D53E4F" ))(50), 
                                  features = i)+NoAxes()+NoLegend()+
    theme(text=element_text(size=12))+ 
    theme(legend.text=element_text(size=8))#+
  #theme(panel.border = element_rect(fill = NA,color = "black",
  #size=1.5,linetype = "solid"))
}

library(cowplot)
library(patchwork)

plots_<-plot_grid(plots[['TRAC']],
                  plots[['MS4A1']],
                  plots[['IGLC2']],
                  plots[['C1QB']],
                  plots[['TPSB2']],
                  plots[['PROX1']],
                  plots[['VEGFC']],
                  plots[['COL12A1']],
                  plots[['ACTA1']],
                  plots[['HIGD1B']],
                  plots[['KRT15']],
                  align = "h",  #axis = 'l',
                  nrow =3)

# featureplot<-wrap_plots(plots, ncol = 5)
# featureplot


plt_featureplot<- FeaturePlot_scCustom(seurat_object = MPRNMPR.SC.subset, reduction = 'umap.harmony',
                                       colors_use = colorRampPalette(c("#A9A9A9", "white", "#D53E4F" ))(50), 
                                       features ='KRT15' )+NoAxes()+ 
  theme(text=element_text(size=12))+ 
  theme(text=element_text(face = "bold"))+
  theme(legend.text=element_text(size=8))

legend <- get_legend(plt_featureplot+
                       #guides(color = guide_legend(nrow = 1)) +
                       theme(legend.position = "left"))
plt_featureplot2<- plot_grid(plots_,legend, ncol=2,rel_widths = c(1,0.1))


ggsave(plt_featureplot2,file="MPR_featureplot.jpg",width =12,height = 7)

ggsave(plt_featureplot2,file="MPR_featureplot.pdf",width =12,height = 7)

######重新featureplot----

setwd("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞注释_dotplot_featureplot")

library(Seurat)
library(ggplot2)
library(cowplot)
library(dplyr)
library(tidydr)
library(stringr)
library(viridis)

markers <- c("TRAC", # Tcells
             "MS4A1", #B cells
             "IGLC2",#Plasma cell
             "C1QB", # Myeloid cells
             "TPSB2", # Mast cells
             "PROX1","VEGFC",#Endothelial cells
             "COL12A1", # Fibroblasts
             "ACTA1", #Smooth muscle cells
             "HIGD1B",#Pericytes
             "KRT15",
             "PDCD1") # Endothelial cells

pal <- viridis(n = 10, option = "C")
# "#0D0887FF" "#47039FFF" "#7301A8FF" "#9C179EFF" "#BD3786FF" "#D8576BFF" "#ED7953FF" "#FA9E3BFF" "#FDC926FF" "#F0F921FF"
pal <- viridis(n = 15, option = "D", direction = -1)

p<-FeaturePlot(object =MPRNMPR_object_miMPRobe,features = markers,reduction = "umap.harmony",raster=FALSE, cols = pal, order = T)

ggsave(file="D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞注释_dotplot_featureplot/MPR细胞大类markerfeatureplot.jpg",p,width=15,height=9)
ggsave(file="D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞注释_dotplot_featureplot/MPR细胞大类markerfeatureplot.pdf",p,width=15,height=9)


i=1
plots=list()
for (i in 1:length(markers)){
  plots[[i]]=FeaturePlot(object=MPRNMPR_object_miMPRobe,features = markers[i],reduction = "umap.harmony",raster=FALSE, cols = pal, order = T)+
    xlab("UMAP_1")+
    ylab("UMAP_1")+
    theme(legend.position = "none")
}
library(patchwork)
p<-wrap_plots(plots, ncol = 4)+plot_annotation(tag_levels = "A");p

p_legend<-FeaturePlot(object=MPRNMPR_object_miMPRobe,features = "PDCD1",reduction = "umap.harmony",raster=FALSE, cols = pal, order = T)+
  xlab("UMAP_1")+
  ylab("UMAP_1")+
  theme(legend.position = "right")

legend <- get_legend(p_legend+
                       #guides(color = guide_legend(nrow = 1)) +
                       theme(legend.position = "right"))
plt_featureplot2<- plot_grid(p,legend, ncol=2,rel_widths = c(1,0.1))

ggsave(file="D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞注释_dotplot_featureplot/MPR细胞大类markerfeatureplot_legend.jpg",plt_featureplot2,width=15,height=9)
ggsave(file="D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞注释_dotplot_featureplot/MPR细胞大类markerfeatureplot_legend.pdf",plt_featureplot2,width=15,height=9)

#ggsave("T_cells_CD3E.jpg",device = "pdf",width = 10,height = 10.5,units = "cm")

####组合dotplot 和featureplot#######

library(imager)
imgjpg <- load.image("MPR_marker_dotplot.jpg")

imgjpg2<-ggplot2::ggplot() + ggplot2::annotation_custom(grid::rasterGrob(imgjpg,
                                                                         width=ggplot2::unit(1,"npc"),
                                                                         height=ggplot2::unit(1,"npc")),
                                                        -Inf, Inf, -Inf, Inf)


marker_dotplot_featureplot<- plot_grid(imgjpg2,
                                       plt_featureplot2,
                                       nrow=2,
                                       rel_heights = c(1,1),
                                       rel_widths =c(0.5,1))

ggsave(marker_dotplot_featureplot,file="MPR_marker_dotplot_featureplot.jpg",width =10,height = 10)

ggsave(marker_dotplot_featureplot,file="MPR_marker_dotplot_featureplot.pdf",width =10,height = 10)


####OR组织偏好性##----

#OR指数比较单细胞亚群的组织偏好----

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
out.prefix <- "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/MPR_OR"

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

meta.tb <- MPRNMPR_object_miMPRobe_remove@meta.data

names(meta.tb)
View(meta.tb)

A <- do.tissueDist(cellInfo.tb = meta.tb,
                   meta.cluster = meta.tb$cell_type_new,
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
data$rid <- factor(data$rid, levels=c('T cells', 'B cells','Plasma cells', 'Myeloid cells',  
                                      'Mast cells',"Endothelial cells","Fibroblasts","Smooth muscle cells",
                                      "Pericytes","Epithelial cells"))
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
data$cid <- factor(data$cid,levels=c("Pre_NMPR","Pre_MPR","Post_MPR"))
OR_plot2<-ggplot(data, aes(cid,rid)) + 
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
OR_plot2

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞比例_OR值计算/OR_plot_MPR_treatments2.pdf",OR_plot2,width=5,height=7)

ggsave("OR_plot_MPR.jpg",OR_plot,width=5,height=7)
ggsave("OR_plot_MPR.pdf",OR_plot,width=5,height=7)


ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞比例_OR值计算/OR_plot_MPR_treatments.pdf",OR_plot,width=9,height=5)

####RO/e
setwd("./cell_OR")
install.packages("devtools")
devtools::install_github("Japrin/sscVis",force = TRUE)
devtools::install_github("Japrin/sscClust",force = TRUE)
source("Roe.R")

distribution_Roe(
  meta_data = metainfo,
  celltype_column = "majorCluster",
  condition_column = "tissue",
  add_label = "number",
  tile_color = "grey",
  tile_fill = "D"
)

####重新画PDCD1在不同分组、不同细胞类型表达#####----
library(viridis)
library(RColorBrewer)

markers <- c("PDCD1")
markers <- as.data.frame(markers)
markerdata <- ScaleData(MPRNMPR_object_miMPRobe_remove, features = as.character(unique(markers$markers)), assay = "RNA")

aver_dt<-AverageExpression(markerdata,
                             features = "PDCD1",
                             group.by =c('cell_type_new','MPR_Response2'))


aver_dt <- as.data.frame(aver_dt$RNA)

aver_dt<-as.data.frame(t(aver_dt))
colnames(aver_dt)<- "Average_Expression"
aver_dt$cell_type <- rownames(aver_dt)

aver_dt <- aver_dt %>% separate(cell_type, c('cell_type', 'response'),sep="_P")

aver_dt$response <- paste("P",aver_dt$response,sep = "")
# aver_dt$response <- paste("PDCD1",aver_dt$response,sep = "")

View(aver_dt)
library(RColorBrewer)
cell_type_order<-c('T cells', 'B cells','Plasma cells', 'Myeloid cells',  
              'Mast cells',"Endothelial cells","Fibroblasts","Smooth muscle cells",
              "Pericytes","Epithelial cells")

aver_dt$cell_type<-factor(aver_dt$cell_type,levels=rev(c('T cells', 'B cells','Plasma cells', 'Myeloid cells',  
                                                     'Mast cells',"Endothelial cells","Fibroblasts","Smooth muscle cells",
                                                     "Pericytes","Epithelial cells")))

write.csv(aver_dt,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/PDCD1每个细胞以及微生物分组表达情况/PDCD1不同细胞类型热图数据.csv")
df1<-aver_dt
PDCD1_plot<-ggplot(df1, aes(cell_type,response)) + 
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
        axis.text.y = element_text(size = 12,color = 'black',
        ),
        axis.text.x = element_text(size = 15,color = 'black',angle =40,hjust =1,vjust=1),
        legend.title = element_text(size=12, color = "black"), 
        legend.text = element_text(size=10,color = "black",angle =45),
        legend.position = "top") + 
  scale_y_discrete(position = "right")
  

PDCD1_plot

ggsave("PDCD1_response_MPR.jpg",PDCD1_plot,width=8,height=5)
ggsave("PDCD1_response_MPR.pdf",PDCD1_plot,width=8,height=5)

setwd("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞亚群注释/T_cell注释/PDCD1有无微生物表达对比")

#####PDCD1 微生物分裂提琴图####----
{markers <- c("PDCD1")
markers <- as.data.frame(markers)

markerdata <- ScaleData(MPRNMPR_object_miMPRobe_remove, features = as.character(unique(markers$markers)), assay = "RNA")

Idents(MPRNMPR_object_miMPRobe) <- "cell_id"

Idents(MPRNMPR_object_miMPRobe)

aver_dt<-AverageExpression(T_cell_markerdata,
                           features = "PDCD1",
                           group.by = "cell_id")


aver_dt <- as.data.frame(aver_dt$RNA)
aver_dt<-as.data.frame(t(aver_dt))
colnames(aver_dt)<- "Average_Expression"

group<-cbind(MPRNMPR_object_miMPRobe_remove@meta.data$cell_id)%>%
  cbind(MPRNMPR_object_miMPRobe_remove@meta.data$orig.ident) %>%
  cbind(MPRNMPR_object_miMPRobe_remove@meta.data$group_microbe2) %>%
  cbind(MPRNMPR_object_miMPRobe_remove@meta.data$MPR_Response2) %>%
  as.data.frame()
colnames(group) <- c("cell_id","orig.ident","group_microbe2","MPR_Response2")

rownames(aver_dt)<- gsub('-', '_', rownames(aver_dt), fixed = TRUE)

aver_dt$microbe_group <- group[match(rownames(aver_dt),group$cell_id),3]
aver_dt$MPR_Response2 <- group[match(rownames(aver_dt),group$cell_id),4]
colnames(aver_dt)

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
ggsave("PDCD1_微生物分裂小提琴图.jpg",a,width = 7,height = 5 )
ggsave("PDCD1_微生物分裂小提琴图.pdf",a,width = 7,height = 5 )
}

######PDCD1 微生物热图----
library(viridis)
library(RColorBrewer)
Idents(markerdata)
T_cell_markerdata<- subset(markerdata,idents = 'T cells')
aver_dt<-AverageExpression(T_cell_markerdata,
                           features = "PDCD1",
                           group.by =c('group_microbe2','MPR_Response2'))


aver_dt <- as.data.frame(aver_dt$RNA)

aver_dt<-as.data.frame(t(aver_dt))
colnames(aver_dt)<- "Average_Expression"
aver_dt$cell_type <- rownames(aver_dt)

aver_dt <- aver_dt %>% separate(cell_type, c('cell_type', 'response'),sep="_P")

aver_dt$response <- paste("P",aver_dt$response,sep = "")
# aver_dt$response <- paste("PDCD1",aver_dt$response,sep = "")

df2<-aver_dt
write.csv(df2,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/PDCD1每个细胞以及微生物分组表达情况/PDCD1微生物热图.csv")


PDCD1_plot2<-ggplot(df2, aes(cell_type,response)) + 
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
        axis.text.y = element_text(size = 12,color = 'black',
        ),
        axis.text.x = element_text(size = 15,color = 'black',angle =40,hjust =1,vjust=1),
        legend.title = element_text(size=12, color = "black"), 
        legend.text = element_text(size=10,color = "black",angle =45),
        legend.position = "top") + 
  scale_y_discrete(position = "right")

PDCD1_plot2

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/PDCD1每个细胞以及微生物分组表达情况/PDCD1_response_microbe.jpg",PDCD1_plot2,width=3,height=5)
ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/PDCD1每个细胞以及微生物分组表达情况/PDCD1_response_microbe_T_cells.pdf",PDCD1_plot2,width=5,height=3)

######PDCD1表达和微生物相关性----
setwd("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/PDCD1每个细胞以及微生物分组表达情况")
PDCD1_expression_matrix <- FetchData(MPRNMPR_object_miMPRobe_remove, vars = "PDCD1")
PDCD1_expression_matrix$cell_id <- rownames(PDCD1_expression_matrix)
names(PDCD1_expression_matrix)
microbe<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/PDCD1每个细胞以及微生物分组表达情况/16s_filter.combined.genus.umi.matrix.csv")
colnames(microbe)
PDCD1_microbe<- left_join(x=PDCD1_expression_matrix,y=microbe,by = join_by(cell_id==barcode)) 
rownames(PDCD1_microbe)<-PDCD1_microbe$cell_id
colnames(PDCD1_microbe)[2]

PDCD1_microbe<-PDCD1_microbe[,-2] %>% filter(complete.cases(.))
View(PDCD1_microbe)
write.csv(PDCD1_microbe,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/PDCD1每个细胞以及微生物分组表达情况/PDCD1_微生物数据.csv")

PDCD1_microbe<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/PDCD1每个细胞以及微生物分组表达情况/PDCD1_微生物数据.csv",header = T,row.names = 1)
names(PDCD1_microbe)[1]

microbe_group <- cbind(MPRNMPR_object_miMPRobe_remove@meta.data$cell_id,
                       MPRNMPR_object_miMPRobe_remove@meta.data$orig.ident,
                       MPRNMPR_object_miMPRobe_remove@meta.data$cell_type_new,
                       MPRNMPR_object_miMPRobe_remove@meta.data$group_microbe2,
                       MPRNMPR_object_miMPRobe_remove@meta.data$Patients,
                       MPRNMPR_object_miMPRobe_remove@meta.data$MPR_Response2,
                       MPRNMPR_object_miMPRobe_remove@meta.data$Reads_counts,
                       MPRNMPR_object_miMPRobe_remove@meta.data$UMIs_counts) %>%
  as.data.frame()

colnames(microbe_group)<-c("cell_id","orig.ident","cell_type_new","group_microbe2","Patients",
                           "MPR_Response2","Reads_counts","UMIs_counts")


names(microbe_group)

PDCD1_microbe$Patients<- microbe_group[match(rownames(PDCD1_microbe),microbe_group$cell_id),5]
PDCD1_microbe$MPR_Response2<- microbe_group[match(rownames(PDCD1_microbe),microbe_group$cell_id),6]
PDCD1_microbe$orig.ident <- microbe_group[match(rownames(PDCD1_microbe),microbe_group$cell_id),2]
ncol(PDCD1_microbe)
colnames(PDCD1_microbe)
head(PDCD1_microbe)

library(dplyr)
PDCD1_microbe_group_df <- PDCD1_microbe %>%
  group_by(orig.ident) %>%
  summarise(across(1:514, mean, .names = "mean_{col}"))

PDCD1_microbe_group_df <- PDCD1_microbe %>%
  group_by(Patients, MPR_Response2) %>%
  summarise(across(1:514, mean, .names = "mean_{col}"))

PDCD1_microbe_group_df$mean_PDCD1
PDCD1_microbe_group_df$group<-paste(PDCD1_microbe_group_df$MPR_Response2,PDCD1_microbe_group_df$Patients,sep="-")
colnames(PDCD1_microbe_group_df)[516]
 
PDCD1_microbe_group_df<- PDCD1_microbe_group_df[,-c(1,2)] %>% as.data.frame()
rownames(PDCD1_microbe_group_df) <- PDCD1_microbe_group_df$group

colnames(PDCD1_microbe_group_df)
PDCD1_microbe_group_df<- PDCD1_microbe_group_df[,-c(515)] %>% as.data.frame()

####相关性棒棒图----
library(psych)
PDCD1 <- PDCD1_microbe_group_df$mean_PDCD1
spe <- PDCD1_microbe_group_df[,-c(1,2)]
View(spe)
colnames(spe)[1]
View(spearman)
spearman <- corr.test(PDCD1, spe, method = 'spearman', adjust = 'bonferroni')
r <- data.frame(spearman$r)  #pearson 相关系数矩阵
p <- data.frame(spearman$p.adj)  #p 值矩阵
p <- data.frame(spearman$p)
#结果整理以便于作图
r$PDCD1 <- rownames(r)
p$PDCD1 <- rownames(p)
r <- melt(r, id = 'PDCD1')
p <- melt(p, id = 'PDCD1')
spearman <- cbind(r, p$value)
colnames(spearman) <- c('PDCD1', 'spe', 'spearman_correlation', 'p.adjusted.value')
spearman_filtered <- subset(spearman, p.adjusted.value< 0.05)

spearman_filtered$spe <- factor(spearman_filtered$spe, levels = colnames(spe))
head(spearman_filtered$spe)  #整理好的环境变量和物种丰度的 pearson 相关性统计表


col_set =c("mean_Duganella" = "#ffb25f",
           "mean_Rothia" = "#bab4d0",
           "mean_Streptococcus" = "#b57aaa",
           "mean_Veillonella" = "#7bbebd",
           "mean_Blattabacterium" = "#f7c9d9",
           "mean_Gemella" = "#aac567",
           "mean_Subdoligranulum" = "#f67e6f",
           "mean_Actinomyces" = "#fdeeaf",
           "mean_Eubacterium" = "#c4dbbc",
           "mean_Leptotrichia" = "#ffdf71")

p4 <- ggplot(spearman_filtered,
            aes(x = spe, y = spearman_correlation)) +
  # 棒棒图的连线
  geom_segment(aes(x = spe, xend = spe, y = 0, yend = spearman_correlation),  
               linetype = "solid", # 实线
               size = 1, # 连线的粗细 
               color = "gray40" # 连线的颜色
  ) + 
  # y轴0刻度的水平线
  geom_hline( 
    yintercept = 0,  # 水平线位置
    linetype = "dashed",  # 虚线
    size = 1, # 连线的粗细 
    colour="gray40" # 连线的颜色
  ) +
  # 绘制点
  geom_point(aes(color = spe),  
             color = col_set,
             size = 19) +   
  # 添加相关性R值标签
  geom_text(aes(label = round(spearman_correlation, 2)), 
            #color = ifelse(round(spearman_correlation, 2) != 0.96, "black", 'red'), 
            size = 5) +
  # 添加相关性P值标签,这一步最精华的一点是，根据r值正负调整水平移动的位置。
  geom_text(aes(label = paste("p=",round(p.adjusted.value, 2),sep="")),
            hjust =-1,
              #ifelse(spearman_filtered$spearman_correlation >= 0, 1.5, -0.5),
            vjust = -0.2,
            angle = 90,
            fontface = 'italic',
            #color = ifelse(data$p != 'p<0.001', "black", 'red'), 
            size = 5) +
  # y轴刻度设置
  scale_y_continuous( #设置y轴
    limits = c(-0.8, 0),
    breaks = c(-0.8, -0.4, 0),
    labels = c(-0.8, -0.4, 0)) +
  # 坐标title和图形title设置
  labs(y = "Spearman correlation coefficient",
    title = paste0("Correlations between PDCD1 expression \n with cell-associated bacteria")) +
  theme_classic() +
  theme(axis.line = element_line(size = 0.7),
    plot.title = element_text(size = 20, hjust = 0.5),
    axis.text.x = element_text(size = 16, angle = 45, hjust=1,colour = "black"), 
    axis.text.y = element_text(size = 16),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 20),
    plot.margin = unit(c(1, 1, 1, 1), "cm")
  )

p4

ggsave("PDCD1和微生物相关性棒棒图.jpg",p4,width = 10,height = 6)
ggsave("PDCD1和微生物相关性棒棒图.pdf",p4,width = 10,height = 6)

####随机森林模型可解释率

#随机森林的模型显著性
#构建模型查看模型分类准确性
library(randomForest)
set.seed(123)
treat_rf <- randomForest(mean_PDCD1~., data =PDCD1_microbe_group_df,importance=TRUE,proximity=TRUE)
treat_perm <- rf.significance(treat_rf, PDCD1_microbe_group_df, nperm=99, ntree=500)
treat_rf
#通过置换检验进一步查看模型的显著性水平
set.seed(123)

treat_perm

# > treat_perm
# Number of permutations:  99 
# p-value:  0.01 
# Model signifiant at p = 0.01 
# Model R-square:  0.879014 
# Random R-square:  -0.1447458 
# Random R-square variance:  0.01448512

#保存至ppt
library(eoffice)
topptx(filename = "组合图.pptx",height=6,width=3)

library(vegan)
library(dplyr)
library(ggcor)
library(ggplot2)

PDCD1_microbe<- read.csv("PDCD1_微生物数据.csv",header = T,row.names = 1)
row_sums <- rowSums(PDCD1_microbe)

# 创建逻辑向量，表示行总和不为 0 的行
non_zero_sum_rows <- row_sums != 0

# 使用逻辑向量过滤数据
data_filtered <-PDCD1_microbe[non_zero_sum_rows, ]
nrow(data_filtered)
PDCD1 <- data_filtered$PDCD1
spe <- data_filtered[,-1]

df_mantel <- mantel_test(PDCD1, spe, mantel.fun = 'mantel',
                         spec.dist.method = 'euclidean', 
                         env.dist.method = 'euclidean',
                         )#将群落数据按组进行分开

df_mantel <- df_mantel %>%
  mutate(df_r = cut(r, breaks = c(-Inf, 0.1, 0.2, 0.4, Inf),
                    labels = c("< 0.1", "0.1 - 0.2", "0.2 - 0.4", ">= 0.4")),#定义Mantel的R值范围标签
         df_p = cut(p.value, breaks = c(-Inf, 0.01, 0.05, Inf),
                    labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))#定义Mantel的P值范围标签

quickcor(env,method = "spearman", type = "upper", cor.test = T, cluster.type = "all") +#环境因子之间的相关性热图
  geom_square() +#相关性显示形式
  geom_mark(r = NA,sig.thres = 0.05, size = 3.5, colour = "black")+#显著性标签
  scale_fill_gradient2( high = 'red', mid = 'white',low = 'blue') + #颜色设置
  anno_link(df_mantel, aes(color = df_p,
                           size = df_r))+
  scale_size_manual(values = c(0.5, 1, 1.5, 2))+#连线粗细设置
  scale_color_manual(values = c("green","blue","orange"))+#线条颜色设置
  guides(fill = guide_colorbar(title = "correlation", order = 1),#图例相关设置
         size = guide_legend(title = "Mantel's r",order = 2),
         color = guide_legend(title = "Mantel's p", order = 3),
         linetype = "none")


###组合图片#####
# library(cowplot)
# PDCD1_MPR<-plot_grid(PDCD1_featureplot,PDCD1_plot,ncol=2,rel_heights = c(1,1.5),rel_widths = c(0.5,1))
# 
# ggsave("PDCD1_MPR.jpg",PDCD1_MPR,width=5,height=3)


#####T cell

#########加入tMPR bMPR 到metadata##############

{tMPR <- read.csv("./TMPR_BMPR_16s_table/tMPR_merged_vdj_contig_annotation.csv")
tMPR_clono <- read.csv("./TMPR_BMPR_16s_table/tMPR_clonotype_count_by_umi_for_each_sample.csv")
tMPR_vdj<-read.csv("./TMPR_BMPR_16s_table/tMPR_merged.csv")
View(tMPR_clono)
tMPR<- left_join(x=tMPR,y=tMPR_vdj,by ="barcode")
tMPR<- left_join(x=tMPR,y=tMPR_clono ,by ="clonotype_id")
View(tMPR)
rownames(tMPR) <- tMPR[,1]
rownames(tMPR)
metadata<- FetchData(OSCC_object,"cell_id")
metadata$cell_id <- rownames(metadata)
metadata<- left_join(x=metadata,y=tMPR,by = join_by(cell_id==barcode))
metadata<- left_join(x=metadata,y=tMPR_vdj,by = join_by(cell_id==barcode))
rownames(metadata)<-metadata$cell_id
OSCC_object_tMPR<- AddMetaData(OSCC_object,metadata = metadata)
View(OSCC_object_tMPR)
table(rownames(OSCC_object_tMPR@meta.data) %in% rownames(tMPR))
OSCC_object_tMPR<- subset(OSCC_object_tMPR,cells =rownames(tMPR))
p2 <-DimPlot(OSCC_object_tMPR, reduction = 'umap.harmony',group.by = "cdr3s_aa",raster=FALSE)
p2
DimPlot(OSCC_object_tMPR, reduction = 'umap.harmony',raster=FALSE)
}

####bMPR加到metadata#######

{bMPR <- read.csv("./TMPR_BMPR_16s_table/bMPR_merged_vdj_contig_annotation.csv")
bMPR_clono <- read.csv("./TMPR_BMPR_16s_table/bMPR_clonotype_count_by_umi_for_each_sample.csv")
head(tMPR_clono)
bMPR<- left_join(x=bMPR,y=bMPR_clono ,by ="clonotype_id")
rownames(bMPR) <- bMPR[,1]
rownames(bMPR)
metadata<- FetchData(OSCC_object,"cell_id")
metadata$cell_id <- rownames(metadata)
metadata<- left_join(x=metadata,y=bMPR,by = join_by(cell_id==barcode))
rownames(metadata)<-metadata$cell_id
OSCC_object_bMPR<- AddMetaData(OSCC_object,metadata = metadata)
View(OSCC_object_bMPR)
table(rownames(OSCC_object_bMPR@meta.data) %in% rownames(bMPR))
OSCC_object_bMPR<- subset(OSCC_object_bMPR,cells =rownames(bMPR))
p2 <-DimPlot(OSCC_object_bMPR, reduction = 'umap.harmony',group.by = "cdr3s_aa",raster=FALSE)
p2

DimPlot(OSCC_object_tMPR, reduction = 'umap.harmony',raster=FALSE)
}
####有无微生物存在对不同细胞类型功能影响####----
setwd("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对细胞大类功能影响")
GetAssay(MPRNMPR_object_miMPRobe_remove,assay = "RNA")
Idents(MPRNMPR_object_miMPRobe_remove)
Idents(MPRNMPR_object_miMPRobe_remove)<-"cell_type_new"
MPRNMPR_object_miMPRobe_remove$group_microbe2
Idents(MPRNMPR_object_miMPRobe_remove)<-"group_microbe2"

deg_all=FindMarkers(MPRNMPR_object_miMPRobe_remove, ident.1 = c("Bacteria+"), ident.2 = "Bacteria-")
dim(deg_all)

cell_type_list <- c('T cells', 'B cells','Plasma cells', 'Myeloid cells',  
  'Mast cells',"Endothelial cells","Fibroblasts","Smooth muscle cells",
  "Pericytes","Epithelial cells")

for (cell in cell_type_list){
  print(cell)
  Idents(MPRNMPR_object_miMPRobe)<-"cell_type_new"
  sub_object=subset(MPRNMPR_object_miMPRobe,idents = cell)
  Idents(sub_object)<-"group_microbe2"
  deg_all=FindMarkers(sub_object, ident.1 = c("Bacteria+"), ident.2 = "Bacteria-")
  new_filename <- paste(cell, "_deg_all",".csv", sep="")
  print(new_filename)
  write.csv(deg_all,new_filename)
  }

# 设置阈值
#logFC_threshold <- 0.25
pct.1_threshold<- 0.5

# 获取所有csv文件列表

files <- list.files(path = getwd(), pattern = "\\.csv$", recursive = TRUE)

# 创建一个空数据框用于存储合并后的数据
merged_data <- data.frame()

# 循环读取每个csv文件并按照某一列降序排序，进行阈值筛选，并合并到merged_data中
for (file in files) {
  data <- read.csv(file)
  sorted_data <- data[order(-data$pct.1), ]
  file<-basename(file)
  print(file)
  file_name <- rep(substr(file, 1, nchar(file)-12), nrow(filtered_data))
  list1<- c("B cells","Plasma cells")
  if(unique(file_name) %in% list1){
    filtered_data <- sorted_data[sorted_data$pct.1 >= 0.25, ]}
  else{
    filtered_data <- sorted_data[sorted_data$pct.1 >= 0.25, ] 
  }
  file_name <- rep(substr(file, 1, nchar(file)-12), nrow(filtered_data))
  filtered_data$group <- file_name
  merged_data <- rbind(merged_data, filtered_data)
}

# 导出合并后的数据为新的csv文件
write.csv(merged_data, "未分组_filtered_merged_data.csv", row.names = FALSE)

######根据差异基因画多组火山图###----
library(ggplot2) # Create Elegant Data Visualisations Using the Grammar of Graphics
library(dplyr) # A Grammar of Data Manipulation
library(RColorBrewer) # ColorBrewer Palettes
library(grid) # The Grid Graphics Package
library(scales)

df <- read.table("未分组_filtered_merged_data.csv", header = 1, check.names = F, sep = ",")
df$group <- factor(df$group, levels = c('T cells', 'B cells','Plasma cells', 'Myeloid cells',  
                                        'Mast cells',"Endothelial cells","Fibroblasts","Smooth muscle cells",
                                        "Pericytes","Epithelial cells"))


df$group2<-as.factor(ifelse(df$p_val_adj < 0.05 & abs(df$avg_log2FC) >= 0.25, 
                            ifelse(df$avg_log2FC>= 0.25 ,'up','down'),'ns'))

#三分组
ifelse(x <= 10, 1,
       ifelse(x <= 20, 2, 3))
df$group2 <- factor(df$group2, levels = c("up", "down", "ns"))
##确定添加标签的数据
df$label<-ifelse(df$p_val_adj<0.05&abs(df$avg_log2FC)>=0.25,"Y","N")
df$label<-ifelse(df$label == 'Y', as.character(df$X), '')
write.csv(df,"未分组_filtered_merged_data_filter.csv")

df_bg <- df %>%
  group_by(group) %>%
  summarize(max_log2FC = max(avg_log2FC),min_log2FC = min(avg_log2FC))

df_bg$group<- factor(df_bg$group, levels = c('T cells', 'B cells','Plasma cells', 'Myeloid cells',  
                                          'Mast cells',"Endothelial cells","Fibroblasts","Smooth muscle cells",
                                          "Pericytes","Epithelial cells"))


library(ggrepel)

cell_type_color <- c("#FF6347","#ffc556","#1E90FF","#0000FF","#9370DB","#00FFFF","#32cd32","#006400","#FF00ff","#20B2AA")

vol_plot<-ggplot()+
  ##y轴正半轴的灰色背景
  geom_col(data = df_bg, 
           mapping = aes(group,max_log2FC),
           fill = "grey85", width = 0.8, alpha = 0.5) +
  ##y轴负半轴的灰色背景
  geom_col(data = df_bg, 
           mapping = aes(group, min_log2FC),
           fill = "grey85", width = 0.8, alpha = 0.5)+
  #添加各组的数据点
  geom_jitter(data = df,
              mapping = aes(x = group, y = avg_log2FC, color = group2),
              size= 4,width = 0.4, alpha = 0.7)+
  # 通过在X=0的位置添加方块进行展示分组信息，采用geom_col方法添加
  geom_col(data = df_bg,
           mapping = aes(x= group, y = 0.1, fill = group),
           width = 0.8)+
  geom_col(data = df_bg,
           mapping = aes(x= group, y = -0.1, fill = group),
           width = 0.8)+
  # 在方块中添加分组的文字信息
  geom_text(data=df_bg,
            mapping = aes(x=group, y=0, label=group),
            size = 3.2, color ="white",fontface = "bold")+
  #根据需要决定是添加辅助线
  # geom_hline(yintercept = 4, lty=2, color = '#ae63e4', lwd=0.8)+
  # geom_hline(yintercept = -4, lty=2, color = '#ae63e4', lwd=0.8)+
  #颜色设置
  scale_color_manual(values = c("#e42313", "#0061d5", "#8b8c8d"))+
  scale_fill_manual(values = cell_type_color)+
  #添加显著性标签
  geom_text_repel(data = df,
                  mapping = aes(x = group, y = avg_log2FC, label = label),
                  max.overlaps = 10000,
                  size=3,
                  box.padding=unit(0.3,'lines'),
                  point.padding=unit(0.3, 'lines'),
                  segment.color='black',
                  show.legend=FALSE)+
  # 主题设置
  theme_classic()+
  theme(axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.y = element_line(linewidth = 0.8),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 14, color = "black"),
        axis.ticks.y = element_line(linewidth = 0.8))+
  labs(x = "Bacteria+ vs Bacteria-", y = "avg_log2FC", fill= NULL, color = NULL)+
  #调整图例
  guides(color=guide_legend(override.aes = list(size=7,alpha=1)))

ggsave("未分组火山图.jpg",vol_plot,width = 15,height = 8)

#背景色
color <- colorRampPalette(brewer.pal(11,"BrBG"))(30)
#添加背景
grid.raster(alpha(color, 0.2), 
            width = unit(1, "npc"), 
            height = unit(1,"npc"),
            interpolate = T)


####Pre_NMPR 分组有无微生物差异基因----
Idents(MPRNMPR_object_miMPRobe) <- "MPR_Response2"
MPRNMPR_object_miMPRobe@meta.data$MPR_Response2
Pre_NMPR_object=subset(MPRNMPR_object_miMPRobe,idents = "Pre_NMPR")

Idents(Pre_NMPR_object)
unique(Pre_NMPR_object@meta.data$MPR_Response2)


{cell_type_list <- c('T cells', 'B cells','Plasma cells', 'Myeloid cells',  
                    'Mast cells',"Endothelial cells","Fibroblasts","Smooth muscle cells",
                    "Pericytes","Epithelial cells")
dir.create("./Pre_NMPR_deg_all")
for (cell in cell_type_list){
  print(cell)
  Idents(Pre_NMPR_object)<-"cell_type_new"
  sub_object=subset(Pre_NMPR_object,idents = cell)
  Idents(sub_object)<-"group_microbe2"
  deg_all=FindMarkers(sub_object, ident.1 = c("Bacteria+"), ident.2 = "Bacteria-")
  new_filename <- paste(cell, "_deg_all",".csv", sep="")
  print(new_filename)
  write.csv(deg_all,paste("./Pre_NMPR_deg_all/",new_filename,sep=""))
}

# 获取所有csv文件列表

files <- list.files(path = getwd(), pattern = "\\.csv$", recursive = TRUE)

files <- list.files(path = "./Pre_NMPR_deg_all/", pattern = "\\.csv$", recursive = TRUE)

# 创建一个空数据框用于存储合并后的数据
merged_data <- data.frame()
# 循环读取每个csv文件并按照某一列降序排序，进行阈值筛选，并合并到merged_data中
for (file in files) {
  data <- read.csv(paste("./Pre_NMPR_deg_all/",file,sep=""))
  sorted_data <- data[order(-data$pct.1), ]
  file<-basename(file)
  print(file)
  file_name <- rep(substr(file, 1, nchar(file)-12), nrow(filtered_data))
  list1<- c("B cells","Plasma cells")
  if(unique(file_name) %in% list1){
    filtered_data <- sorted_data[sorted_data$pct.1 >= 0.25, ]}
  else{
    filtered_data <- sorted_data[sorted_data$pct.1 >= 0.25, ] 
  }
  file_name <- rep(substr(file, 1, nchar(file)-12), nrow(filtered_data))
  filtered_data$group <- file_name
  merged_data <- rbind(merged_data, filtered_data)
}

# 导出合并后的数据为新的csv文件
write.csv(merged_data, "Pre_NMPR_filtered_merged_data.csv", row.names = FALSE)


######根据差异基因画多组火山图###----
library(ggplot2) # Create Elegant Data Visualisations Using the Grammar of Graphics
library(dplyr) # A Grammar of Data Manipulation
library(RColorBrewer) # ColorBrewer Palettes
library(grid) # The Grid Graphics Package
library(scales)

df <- read.table("Pre_NMPR_filtered_merged_data.csv", header = 1, check.names = F, sep = ",")
df$group <- factor(df$group, levels = c('T cells', 'B cells','Plasma cells', 'Myeloid cells',  
                                        'Mast cells',"Endothelial cells","Fibroblasts","Smooth muscle cells",
                                        "Pericytes","Epithelial cells"))


df$group2<-as.factor(ifelse(df$p_val_adj < 0.05 & abs(df$avg_log2FC) >= 0.25, 
                            ifelse(df$avg_log2FC>= 0.25 ,'up','down'),'ns'))
df$group2 <- factor(df$group2, levels = c("up", "down", "ns"))
##确定添加标签的数据
df$label<-ifelse(df$p_val_adj<0.05&abs(df$avg_log2FC)>=3,"Y","N")
df$label<-ifelse(df$label == 'Y', as.character(df$X), '')
write.csv(df,"Pre_NMPR_filtered_merged_data_filter_label.csv")

df_bg <- df %>%
  group_by(group) %>%
  summarize(max_log2FC = max(avg_log2FC),min_log2FC = min(avg_log2FC))

df_bg$group<- factor(df_bg$group, levels = c('T cells', 'B cells','Plasma cells', 'Myeloid cells',  
                                             'Mast cells',"Endothelial cells","Fibroblasts","Smooth muscle cells",
                                             "Pericytes","Epithelial cells"))


library(ggrepel)

cell_type_color <- c("#FF6347","#ffc556","#1E90FF","#0000FF","#9370DB","#00FFFF","#32cd32","#006400","#FF00ff","#20B2AA")

vol_plot_pre_nmpr<-ggplot()+
  ##y轴正半轴的灰色背景
  geom_col(data = df_bg, 
           mapping = aes(group,max_log2FC),
           fill = "grey85", width = 0.8, alpha = 0.5) +
  ##y轴负半轴的灰色背景
  geom_col(data = df_bg, 
           mapping = aes(group, min_log2FC),
           fill = "grey85", width = 0.8, alpha = 0.5)+
  #添加各组的数据点
  geom_jitter(data = df,
              mapping = aes(x = group, y = avg_log2FC, color = group2),
              size= 2,width = 0.4, alpha = 0.7)+
  # 通过在X=0的位置添加方块进行展示分组信息，采用geom_col方法添加
  geom_col(data = df_bg,
           mapping = aes(x= group, y = 0.2, fill = group),
           width = 0.8)+
  geom_col(data = df_bg,
           mapping = aes(x= group, y = -0.2, fill = group),
           width = 0.8)+
  # 在方块中添加分组的文字信息
  geom_text(data=df_bg,
            mapping = aes(x=group, y=0, label=group),
            size = 3.2, color ="white",fontface = "bold")+
  #根据需要决定是添加辅助线
  # geom_hline(yintercept = 4, lty=2, color = '#ae63e4', lwd=0.8)+
  # geom_hline(yintercept = -4, lty=2, color = '#ae63e4', lwd=0.8)+
  #颜色设置
  scale_color_manual(values = c("#e42313", "#0061d5", "#8b8c8d"))+
  scale_fill_manual(values = cell_type_color)+
  #添加显著性标签
  geom_text_repel(data = df,
                  mapping = aes(x = group, y = avg_log2FC, label = label),
                  max.overlaps = 10000,
                  size=4,
                  box.padding=unit(0.3,'lines'),
                  point.padding=unit(0.3, 'lines'),
                  segment.color='black',
                  show.legend=FALSE)+
  # 主题设置
  theme_classic()+
  theme(axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.y = element_line(linewidth = 0.8),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 14, color = "black"),
        axis.ticks.y = element_line(linewidth = 0.8))+
  labs(x = "Pre_NMPR Bacteria+ vs Bacteria-", y = "avg_log2FC", fill= NULL, color = NULL)+
  #调整图例
  guides(color=guide_legend(override.aes = list(size=7,alpha=1)))
}
ggsave("Pre_NMPR_火山图.jpg",vol_plot_pre_nmpr,width = 15,height = 8)

####Pre_MPR分组有无微生物基因差异分析----

Idents(MPRNMPR_object_miMPRobe) <- "MPR_Response2"
MPRNMPR_object_miMPRobe@meta.data$MPR_Response2
Pre_MPR_object=subset(MPRNMPR_object_miMPRobe,idents = "Pre_MPR")

Idents(Pre_MPR_object)
unique(Pre_MPR_object@meta.data$MPR_Response2)


{cell_type_list <- c('T cells', 'B cells','Plasma cells', 'Myeloid cells',  
                     'Mast cells',"Endothelial cells","Fibroblasts","Smooth muscle cells",
                     "Pericytes","Epithelial cells")
  dir.create("./Pre_MPR_deg_all")
  for (cell in cell_type_list){
    print(cell)
    Idents(Pre_MPR_object)<-"cell_type_new"
    sub_object=subset(Pre_MPR_object,idents = cell)
    Idents(sub_object)<-"group_microbe2"
    deg_all=FindMarkers(sub_object, ident.1 = c("Bacteria+"), ident.2 = "Bacteria-")
    new_filename <- paste(cell, "_deg_all",".csv", sep="")
    print(new_filename)
    write.csv(deg_all,paste("./Pre_MPR_deg_all/",new_filename,sep=""))
  }
  
  # 获取所有csv文件列表
  
  files <- list.files(path = getwd(), pattern = "\\.csv$", recursive = TRUE)
  
  files <- list.files(path = "./Pre_MPR_deg_all/", pattern = "\\.csv$", recursive = TRUE)
  
  # 创建一个空数据框用于存储合并后的数据
  merged_data <- data.frame()
  # 循环读取每个csv文件并按照某一列降序排序，进行阈值筛选，并合并到merged_data中
  for (file in files) {
    data <- read.csv(paste("./Pre_MPR_deg_all/",file,sep=""))
    sorted_data <- data[order(-data$pct.1), ]
    file<-basename(file)
    print(file)
    file_name <- rep(substr(file, 1, nchar(file)-12), nrow(filtered_data))
    list1<- c("B cells","Plasma cells")
    if(unique(file_name) %in% list1){
      filtered_data <- sorted_data[sorted_data$pct.1 >= 0.25, ]}
    else{
      filtered_data <- sorted_data[sorted_data$pct.1 >= 0.25, ] 
    }
    file_name <- rep(substr(file, 1, nchar(file)-12), nrow(filtered_data))
    filtered_data$group <- file_name
    merged_data <- rbind(merged_data, filtered_data)
  }
  
  # 导出合并后的数据为新的csv文件
  write.csv(merged_data, "Pre_MPR_filtered_merged_data.csv", row.names = FALSE)
  
  
  ######根据差异基因画多组火山图###----
  library(ggplot2) # Create Elegant Data Visualisations Using the Grammar of Graphics
  library(dplyr) # A Grammar of Data Manipulation
  library(RColorBrewer) # ColorBrewer Palettes
  library(grid) # The Grid Graphics Package
  library(scales)
  
  df <- read.table("Pre_MPR_filtered_merged_data.csv", header = 1, check.names = F, sep = ",")
  df$group <- factor(df$group, levels = c('T cells', 'B cells','Plasma cells', 'Myeloid cells',  
                                          'Mast cells',"Endothelial cells","Fibroblasts","Smooth muscle cells",
                                          "Pericytes","Epithelial cells"))
  
  
  df$group2<-as.factor(ifelse(df$p_val_adj < 0.05 & abs(df$avg_log2FC) >= 0.25, 
                              ifelse(df$avg_log2FC>= 0.25 ,'up','down'),'ns'))
  df$group2 <- factor(df$group2, levels = c("up", "down", "ns"))
  ##确定添加标签的数据
  df$label<-ifelse(df$p_val_adj<0.05&abs(df$avg_log2FC)>=4,"Y","N")
  df$label<-ifelse(df$label == 'Y', as.character(df$X), '')
  write.csv(df,"Pre_MPR_filtered_merged_data_filter_label.csv")
  
  df_bg <- df %>%
    group_by(group) %>%
    summarize(max_log2FC = max(avg_log2FC),min_log2FC = min(avg_log2FC))
  
  df_bg$group<- factor(df_bg$group, levels = c('T cells', 'B cells','Plasma cells', 'Myeloid cells',  
                                               'Mast cells',"Endothelial cells","Fibroblasts","Smooth muscle cells",
                                               "Pericytes","Epithelial cells"))
  
  
  library(ggrepel)
  
  cell_type_color <- c("#FF6347","#ffc556","#1E90FF","#0000FF","#9370DB","#00FFFF","#32cd32","#006400","#FF00ff","#20B2AA")
  
  vol_plot_pre_mpr<-ggplot()+
    ##y轴正半轴的灰色背景
    geom_col(data = df_bg, 
             mapping = aes(group,max_log2FC),
             fill = "grey85", width = 0.8, alpha = 0.5) +
    ##y轴负半轴的灰色背景
    geom_col(data = df_bg, 
             mapping = aes(group, min_log2FC),
             fill = "grey85", width = 0.8, alpha = 0.5)+
    #添加各组的数据点
    geom_jitter(data = df,
                mapping = aes(x = group, y = avg_log2FC, color = group2),
                size= 2,width = 0.4, alpha = 0.7)+
    # 通过在X=0的位置添加方块进行展示分组信息，采用geom_col方法添加
    geom_col(data = df_bg,
             mapping = aes(x= group, y = 0.2, fill = group),
             width = 0.8)+
    geom_col(data = df_bg,
             mapping = aes(x= group, y = -0.2, fill = group),
             width = 0.8)+
    # 在方块中添加分组的文字信息
    geom_text(data=df_bg,
              mapping = aes(x=group, y=0, label=group),
              size = 3.2, color ="white",fontface = "bold")+
    #根据需要决定是添加辅助线
    # geom_hline(yintercept = 4, lty=2, color = '#ae63e4', lwd=0.8)+
    # geom_hline(yintercept = -4, lty=2, color = '#ae63e4', lwd=0.8)+
    #颜色设置
    scale_color_manual(values = c("#e42313", "#0061d5", "#8b8c8d"))+
    scale_fill_manual(values = cell_type_color)+
    #添加显著性标签
    geom_text_repel(data = df,
                    mapping = aes(x = group, y = avg_log2FC, label = label),
                    max.overlaps = 10000,
                    size=4,
                    box.padding=unit(0.3,'lines'),
                    point.padding=unit(0.3, 'lines'),
                    segment.color='black',
                    show.legend=FALSE)+
    # 主题设置
    theme_classic()+
    theme(axis.line.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.line.y = element_line(linewidth = 0.8),
          axis.text.y = element_text(size = 12, color = "black"),
          axis.title = element_text(size = 14, color = "black"),
          axis.ticks.y = element_line(linewidth = 0.8))+
    labs(x = "Pre_MPR Bacteria+ vs Bacteria-", y = "avg_log2FC", fill= NULL, color = NULL)+
    #调整图例
    guides(color=guide_legend(override.aes = list(size=7,alpha=1)))
}
ggsave("Pre_MPR_火山图.jpg",vol_plot_pre_mpr,width = 15,height = 8)


####Post_MPR分组有无微生物基因差异分析----

Idents(MPRNMPR_object_miMPRobe) <- "MPR_Response2"
MPRNMPR_object_miMPRobe@meta.data$MPR_Response2
Post_MPR_object=subset(MPRNMPR_object_miMPRobe,idents = "Post_MPR")

Idents(Post_MPR_object)
unique(Post_MPR_object@meta.data$MPR_Response2)


{cell_type_list <- c('T cells', 'B cells','Plasma cells', 'Myeloid cells',  
                     'Mast cells',"Endothelial cells","Fibroblasts","Smooth muscle cells",
                     "Pericytes","Epithelial cells")
  dir.create("./Post_MPR_deg_all")
  for (cell in cell_type_list){
    print(cell)
    Idents(Post_MPR_object)<-"cell_type_new"
    sub_object=subset(Post_MPR_object,idents = cell)
    Idents(sub_object)<-"group_microbe2"
    deg_all=FindMarkers(sub_object, ident.1 = c("Bacteria+"), ident.2 = "Bacteria-")
    new_filename <- paste(cell, "_deg_all",".csv", sep="")
    print(new_filename)
    write.csv(deg_all,paste("./Post_MPR_deg_all/",new_filename,sep=""))
  }
  
  # 获取所有csv文件列表
  files <- list.files(path = "./Post_MPR_deg_all/", pattern = "\\.csv$", recursive = TRUE)
  
  # 创建一个空数据框用于存储合并后的数据
  merged_data <- data.frame()
  # 循环读取每个csv文件并按照某一列降序排序，进行阈值筛选，并合并到merged_data中
  for (file in files) {
    data <- read.csv(paste("./Post_MPR_deg_all/",file,sep=""))
    sorted_data <- data[order(-data$pct.1), ]
    file<-basename(file)
    print(file)
    file_name <- rep(substr(file, 1, nchar(file)-12), nrow(filtered_data))
    list1<- c("B cells","Plasma cells")
    if(unique(file_name) %in% list1){
      filtered_data <- sorted_data[sorted_data$pct.1 >= 0.25, ]}
    else{
      filtered_data <- sorted_data[sorted_data$pct.1 >= 0.25, ] 
    }
    file_name <- rep(substr(file, 1, nchar(file)-12), nrow(filtered_data))
    filtered_data$group <- file_name
    merged_data <- rbind(merged_data, filtered_data)
  }
  
  # 导出合并后的数据为新的csv文件
  write.csv(merged_data, "Post_MPR_filtered_merged_data.csv", row.names = FALSE)
  
  
  ######根据差异基因画多组火山图###----
  library(ggplot2) # Create Elegant Data Visualisations Using the Grammar of Graphics
  library(dplyr) # A Grammar of Data Manipulation
  library(RColorBrewer) # ColorBrewer Palettes
  library(grid) # The Grid Graphics Package
  library(scales)
  
  df <- read.table("Post_MPR_filtered_merged_data.csv", header = 1, check.names = F, sep = ",")
  df$group <- factor(df$group, levels = c('T cells', 'B cells','Plasma cells', 'Myeloid cells',  
                                          'Mast cells',"Endothelial cells","Fibroblasts","Smooth muscle cells",
                                          "Pericytes","Epithelial cells"))
  
  
  df$group2<-as.factor(ifelse(df$p_val_adj < 0.05 & abs(df$avg_log2FC) >= 0.25, 
                              ifelse(df$avg_log2FC>= 0.25 ,'up','down'),'ns'))
  df$group2 <- factor(df$group2, levels = c("up", "down", "ns"))
  ##确定添加标签的数据
  df$label<-ifelse(df$p_val_adj<0.05&abs(df$avg_log2FC)>=2,"Y","N")
  df$label<-ifelse(df$label == 'Y', as.character(df$X), '')
  write.csv(df,"Post_MPR_filtered_merged_data_filter_label.csv")
  
  df_bg <- df %>%
    group_by(group) %>%
    summarize(max_log2FC = max(avg_log2FC),min_log2FC = min(avg_log2FC))
  
  df_bg$group<- factor(df_bg$group, levels = c('T cells', 'B cells','Plasma cells', 'Myeloid cells',  
                                               'Mast cells',"Endothelial cells","Fibroblasts","Smooth muscle cells",
                                               "Pericytes","Epithelial cells"))
  
  
  library(ggrepel)
  
  cell_type_color <- c("#FF6347","#ffc556","#1E90FF","#0000FF","#9370DB","#00FFFF","#32cd32","#006400","#FF00ff","#20B2AA")
  
  vol_plot_post_mpr<-ggplot()+
    ##y轴正半轴的灰色背景
    geom_col(data = df_bg, 
             mapping = aes(group,max_log2FC),
             fill = "grey85", width = 0.8, alpha = 0.5) +
    ##y轴负半轴的灰色背景
    geom_col(data = df_bg, 
             mapping = aes(group, min_log2FC),
             fill = "grey85", width = 0.8, alpha = 0.5)+
    #添加各组的数据点
    geom_jitter(data = df,
                mapping = aes(x = group, y = avg_log2FC, color = group2),
                size= 2,width = 0.4, alpha = 0.7)+
    # 通过在X=0的位置添加方块进行展示分组信息，采用geom_col方法添加
    geom_col(data = df_bg,
             mapping = aes(x= group, y = 0.1, fill = group),
             width = 0.8)+
    geom_col(data = df_bg,
             mapping = aes(x= group, y = -0.1, fill = group),
             width = 0.8)+
    # 在方块中添加分组的文字信息
    geom_text(data=df_bg,
              mapping = aes(x=group, y=0, label=group),
              size = 3.2, color ="white",fontface = "bold")+
    #根据需要决定是添加辅助线
    # geom_hline(yintercept = 4, lty=2, color = '#ae63e4', lwd=0.8)+
    # geom_hline(yintercept = -4, lty=2, color = '#ae63e4', lwd=0.8)+
    #颜色设置
    scale_color_manual(values = c("#e42313", "#0061d5", "#8b8c8d"))+
    scale_fill_manual(values = cell_type_color)+
    #添加显著性标签
    geom_text_repel(data = df,
                    mapping = aes(x = group, y = avg_log2FC, label = label),
                    max.overlaps = 10000,
                    size=4,
                    box.padding=unit(0.3,'lines'),
                    point.padding=unit(0.3, 'lines'),
                    segment.color='black',
                    show.legend=FALSE)+
    # 主题设置
    theme_classic()+
    theme(axis.line.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.line.y = element_line(linewidth = 0.8),
          axis.text.y = element_text(size = 12, color = "black"),
          axis.title = element_text(size = 14, color = "black"),
          axis.ticks.y = element_line(linewidth = 0.8))+
    labs(x = "Post_MPR Bacteria+ vs Bacteria-", y = "avg_log2FC", fill= NULL, color = NULL)+
    #调整图例
    guides(color=guide_legend(override.aes = list(size=7,alpha=1)))
}
ggsave("Post_MPR_火山图.jpg",vol_plot_post_mpr,width = 15,height = 8)

library(cowplot)

vol_plot_pre_nmpr|vol_plot_pre_mpr+vol_plot_post_mpr
plots_vol<-plot_grid(vol_plot_pre_nmpr,
                  vol_plot_pre_mpr,
                  vol_plot_post_mpr,
                  align = "h",
                  nrow =2,
                  ncol=2)


ggsave("不同分组不同细胞类型有无微生物差异基因火山图.jpg",plots_vol,width = 20,height = 15)

####上调下调基因venn 图----

devtools::install_github("jolars/eulerr") #开发版
library(eulerr)

setwd("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对细胞大类功能影响/不同分组上调下调基因venn图")

dt_up <- read.csv('up_venn_data.csv',header = T)
up <- list(
  Pre_NMPR = unique(dt_up$Pre_NMPR_up),
  Pre_MPR = unique(dt_up$Pre_MPR_up),
  Post_MPR = unique(dt_up$Post_MPR_up)
)

View(up)

library(ggvenn)

ggvenn(up,
       fill_color = c("#0099e5", "#ff4c4c", "#34bf49"), 
       fill_alpha = 0.5,
       show_percentage = T,
       digits = 1,
       set_name_color = "black",
       set_name_size = 5,
       text_color = "black",
       text_size = 4, 
       stroke_color = "white",
       stroke_alpha = 0.5,
       stroke_size = 0.5,
       stroke_linetype = "solid"
)->p1_up
p1_up



plot(
  euler(up,shape = "circle"),
  fills = list(fill = c("#0099e5","#ff4c4c","#34bf49"),alpha=0.6), 
  labels=list(col="black",font=2, fontfamily = "serif",cex=1.5),
  edges = list(col = "white", lwd=2, lty=1),
  quantities = list(type = c("counts","percent"),cex=1),
) ->p_up
p_up



dt_down <- read.csv('down_venn_data.csv',header = T)
down <- list(
  Pre_NMPR = unique(dt_down$Pre_NMPR_down),
  Pre_MPR = unique(dt_down$Pre_MPR_down),
  Post_MPR = unique(dt_down$Post_MPR_down)
)

View(down)

library(ggvenn)

ggvenn(down,
       fill_color = c("#0099e5", "#ff4c4c", "#34bf49"), 
       fill_alpha = 0.5,
       show_percentage = T,
       digits = 1,
       set_name_color = "black",
       set_name_size = 5,
       text_color = "black",
       text_size = 4, 
       stroke_color = "white",
       stroke_alpha = 0.5,
       stroke_size = 0.5,
       stroke_linetype = "solid"
)->p1_down
p1_down



plot(
  euler(down,shape = "circle"),
  fills = list(fill = c("#0099e5","#ff4c4c","#34bf49"),alpha=0.6), 
  labels=list(col="black",font=2, fontfamily = "serif",cex=1.5),
  edges = list(col = "white", lwd=2, lty=1),
  quantities = list(type = c("counts","percent"),cex=1),
) ->p_down
p_down



####组合venn图----

library(patchwork)
plots_venn<-plot_grid(p1_up,
                     p1_down,
                     nrow =1,
                     ncol=2)

plots_venn

ggsave("差异基因venn图.jpg",plots_venn,width = 7,height = 6)



plots_vol_venn<-plot_grid(vol_plot_pre_nmpr,
                     vol_plot_pre_mpr,
                     vol_plot_post_mpr,
                     plots_venn,
                     align = "h",
                     nrow =2,
                     ncol=2)

ggsave("../火山图和差异基因venn图.jpg",plots_vol_venn,width = 20,height = 15)
ggsave("../火山图和差异基因venn图.pdf",plots_vol_venn,width = 20,height = 15)

####venn 交叉集元素提取----
library(VennDiagram)

inter <- get.venn.partitions(up)
View(inter)

for (i in 1:nrow(inter)) inter[i,'values'] <- paste(inter[[i,'..values..']], collapse = ', ')
write.table(inter[-c(4, 5)], 'venn4_inter_up.txt', row.names = FALSE, sep = '\t', quote = FALSE)

inter <- get.venn.partitions(down)
View(inter)

for (i in 1:nrow(inter)) inter[i,'values'] <- paste(inter[[i,'..values..']], collapse = ', ')
write.table(inter[-c(4, 5)], 'venn4_inter_down.txt', row.names = FALSE, sep = '\t', quote = FALSE)

####High Fuso 和total Bacteria- 有无微生物差异分析----
####添加微生物信息

miMPRobe_data<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/PDCD1每个细胞以及微生物分组表达情况/16s_filter.combined.genus.umi.matrix.csv")

metadata<- FetchData(MPRNMPR_object_miMPRobe_remove,"cell_id")
metadata$cell_id <- rownames(metadata)
metadata<- left_join(x=metadata,y=miMPRobe_data,by = join_by(cell_id==barcode))
rownames(metadata)<-metadata$cell_id
MPRNMPR_object_miMPRobe_remove<- AddMetaData(MPRNMPR_object_miMPRobe_remove,metadata = metadata)
table(rownames(MPRNMPR_object_miMPRobe_remove@meta.data) %in% miMPRobe_data$barcode)

MPRNMPR_object_miMPRobe_remove@meta.data$MPR_Response2
MPRNMPR_object_miMPRobe_remove$group_microbe2
MPRNMPR_object_miMPRobe_remove$fuso_group<-MPRNMPR_object_miMPRobe_remove$group_microbe2
Idents(MPRNMPR_object_miMPRobe_remove)<-"cell_type_new"
Immune_cell=subset(MPRNMPR_object_miMPRobe_remove,idents = c("T cells","B cells",
                                                             "Plasma cells","Myeloid cells","Mast cells"))

Immune_cell$fuso_group <- ifelse(Immune_cell$Fusobacterium > 0, 
                                                         "Fuso", NA)
Immune_cell$fuso_group_new <- ifelse(!is.na(Immune_cell$fuso_group), 
                                           "Fuso", 
                                     Immune_cell$group_microbe2)

Idents(Immune_cell)<-"fuso_group_new"
table(Immune_cell$fuso_group_new)

Fuso_Immune_cell_deg_all=FindMarkers(Immune_cell, ident.1 = c("Fuso"), ident.2 = "Bacteria-")

dim(Fuso_Immune_cell_deg_all)
write.csv(Fuso_Immune_cell_deg_all,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/Fuso/Fuso_Immune_cell_deg_all.csv")
Fuso_Immune_cell_deg_all<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/Fuso/Fuso_Immune_cell_deg_all.csv",header = T)

write.csv(Fuso_Immune_cell_deg_all,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/免疫细胞/Fuso_Immune_cell_deg_all.csv")
Fuso_Immune_cell_deg_all<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/Fuso/Fuso_Immune_cell_deg_all.csv",header = T)

####Fuso_Immune_deg----

# gene_name<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/免疫细胞/免疫细胞gene_name.csv")
# 
# intersection <- intersect(gene_name$Gene.symbol, Fuso_Immune_cell_deg_all$X)
# 
# Fuso_Immune_cell_deg_all_inter<- Fuso_Immune_cell_deg_all %>% filter(Fuso_Immune_cell_deg_all$X %in%intersection )

df <-Fuso_Immune_cell_deg_all

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
  ggtitle("Immune cell DEG (Bacteria+ vs Bacteria-)") + #标题
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

###################################
library(org.Hs.eg.db)#物种注释包(Homo sapiens)
library(clusterProfiler)#富集分析
#添加entrez ID列：
##：symbol转entrez ID：
symbol <- Fuso_Immune_cell_deg_all$X
entrez <- bitr(symbol,
               fromType = "SYMBOL",#现有的ID类型
               toType = "ENTREZID",#需转换的ID类型
               OrgDb = "org.Hs.eg.db")
head(entrez)
genelist <- Fuso_Immune_cell_deg_all$avg_log2FC
names(genelist) <-Fuso_Immune_cell_deg_all$X
head(genelist)
#genelist过滤(ID转换中丢失的部分基因)：
genelist <- genelist[names(genelist) %in% entrez[,1]]
names(genelist) <- entrez[match(names(genelist),entrez[,1]),2]
length(genelist)#过滤后剩16161个基因
head(genelist)
#将genelist按照log2FC值从高到低进行排序：
genelist <- sort(genelist,decreasing = T)
head(genelist)

#GSEA_KEGG富集分析：
KEGG_ges <- gseKEGG(
  geneList = genelist,
  organism = "hsa",#不同物种选择可官方查询：https://www.genome.jp/kegg/catalog/org_list.html
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose = FALSE,
  eps = 0
)

#提取结果表：
KEGG_ges_result <- KEGG_ges@result

#显示前10列（第十一列为富集到该基因集的核心基因列表，这里过长不展示）：
KEGG_ges_result[1:3,1:10]
library(ggridges)
library(ggplot2)
library(enrichplot)
nrow(KEGG_ges_result)
p <- ridgeplot(KEGG_ges,
               showCategory = 15,
               fill = "p.adjust",
               decreasing  = T)
p

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/Fuso/Fuso_免疫细胞/Fuso_ridgeplot_KEGG.pdf",p,width=7,height=6)


##GSEA_GO富集分析：
GO_ges <- gseGO(geneList = genelist,
                OrgDb = org.Hs.eg.db,
                ont = "CC", #one of "BP", "MF", and "CC" subontologies, or "ALL" for all three.
                minGSSize = 10,
                maxGSSize = 500,
                pvalueCutoff = 0.05,
                eps = 0,
                verbose = FALSE)
GO_ges[1:3,1:10] #展示同样省略最后一列
nrow(GO_ges)

p <- ridgeplot(GO_ges,
               showCategory = 29,
               fill = "p.adjust",
               decreasing  = T)
p

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/Fuso/Fuso_免疫细胞/Fuso_ridgeplot_GO.pdf",p,width=7,height=9)


##GSEA_Reactome富集分析:
BiocManager::install("ReactomePA")
library(ReactomePA)
ret_ges <- gsePathway(genelist,
                      organism = "human",
                      minGSSize = 10,
                      maxGSSize = 500,
                      pvalueCutoff = 0.05,
                      pAdjustMethod = "BH",
                      verbose = FALSE,
                      eps = 0)
ret_ges_result<-ret_ges[1:34,1:10]
nrow(ret_ges)
p <- ridgeplot(ret_ges,
               showCategory = 34,
               fill = "p.adjust",
               decreasing  = T)
p
?ridgeplot

pp <- dotplot(ret_ges,
              showCategory = 34)
pp

####可视化柱状图
library(stringr)
library(ggplot2)
ret_ges_result <- ret_ges_result[order(ret_ges_result$NES),]
ret_ges_result$yax <- ifelse(ret_ges_result$NES >0, -0.02, 0.02)
ret_ges_result$col <- ifelse(ret_ges_result$NES > 0, "blue","red")
ret_ges_result$NES
names(ret_ges_result)

p <- 
  ggplot(ret_ges_result,aes(y=NES,x=reorder(Description,NES),label = Description))+
  geom_text(aes(y = yax,hjust = NES > 0),size=4)+ 
  geom_bar(stat="identity",aes(fill=col))+
  geom_hline(yintercept = 0)+
  coord_flip(ylim = c(-2,2.5))+
  labs(y = "Normalised enrichment score", x = "",title="")+
  scale_fill_manual(values = c("#FB6542","#375E97"),name="Expression",labels=c("Up","Down"))+
  scale_x_discrete(breaks = NULL)+
  theme_classic(base_size =  12)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.y = element_blank(),
        legend.position = 'bottom')
p


ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/Fuso/Fuso_免疫细胞/Fuso_barplot_ReactomePA.pdf",p,width=7,height=9)


##GSEA_DO(Disease Ontology)富集分析:
library(DOSE)
DO_ges <- gseDO(genelist,
                minGSSize = 10,
                maxGSSize = 500,
                pvalueCutoff = 0.05,
                pAdjustMethod = "BH",
                verbose = FALSE,
                eps = 0)
DO_ges[1:3,1:10]
nrow(DO_ges)

p <- ridgeplot(DO_ges,
               showCategory = 10,
               fill = "p.adjust",
               decreasing  = T)
p


#####GSEA_MSigDB富集分析：
library(gggsea)
library(ggplot2)
library(GSVA)
library(enrichplot)
library(GSEABase)
library(fgsea)
library(clusterProfiler)
gmtfile <- "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/h.all.v2024.1.Hs.symbols.gmt"
hallmark <- read.gmt(gmtfile)
alldiff <- Fuso_Immune_cell_deg_all[order(Fuso_Immune_cell_deg_all$avg_log2FC,decreasing = T),]
id <- alldiff$avg_log2FC
names(id) <- alldiff$X

## 这里去掉了基因集前缀
hallmark$term <- gsub('HALLMARK_','',hallmark$term)
hallmark.list <- hallmark %>% split(.$term) %>% lapply( "[[", 2)

gsea.re1<-clusterProfiler::GSEA(id,TERM2GENE=hallmark,verbose = T)
g1<-as.data.frame(gsea.re1)
g1<-subset(g1,p.adjust<0.05)
g1<-g1[order(g1$NES,decreasing = T),]
type(g1)
head(g1)

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


library(gggsea)

library(ggsci)

names(hallmark.list)
se_hall<-c(head(g2$pathway,9),tail(g1$pathway,9))

#3 所选中的基因集
sig1<-g2[g2$pathway%in%se_hall,]

hallmark.se<-hallmark.list[sig1$pathway]
df.new <- gseaCurve(id, hallmark.se, sig1)

pal_line<-pal_lancet()(3)#各曲线颜色

p_nes<-ggplot() + 
  geom_gsea(df.new,
            prettyGSEA=T,
            tickcolor='grey30',#竖线颜色
            ticksize=0.1,#线条粗细
            #colour='grey',#色块边框
            linecolor=pal_line,
            linesize=2,lty=1,
            #多个图形时使用，设定每行列图形个数
            nrow = 3,
            ncol =1
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

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/Fuso/Fuso_免疫细胞/Fuso_hallmark.pdf",p_nes,width=4,height=6)




####Strep和total bacteria- 有无微生物差异----

miMPRobe_data<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/PDCD1每个细胞以及微生物分组表达情况/16s_filter.combined.genus.umi.matrix.csv")
metadata<- FetchData(MPRNMPR_object_miMPRobe_remove,"cell_id")
metadata$cell_id <- rownames(metadata)
metadata<- left_join(x=metadata,y=miMPRobe_data,by = join_by(cell_id==barcode))
rownames(metadata)<-metadata$cell_id
MPRNMPR_object_miMPRobe_remove<- AddMetaData(MPRNMPR_object_miMPRobe_remove,metadata = metadata)
table(rownames(MPRNMPR_object_miMPRobe_remove@meta.data) %in% miMPRobe_data$barcode)
MPRNMPR_object_miMPRobe_remove@meta.data$MPR_Response2
MPRNMPR_object_miMPRobe_remove$group_microbe2
MPRNMPR_object_miMPRobe_remove$Strep_group<-MPRNMPR_object_miMPRobe_remove$group_microbe2
Idents(MPRNMPR_object_miMPRobe_remove)<-"cell_type_new"

#####strep 免疫细胞----

Immune_cell=subset(MPRNMPR_object_miMPRobe_remove,idents = c("T cells","B cells",
                                                             "Plasma cells","Myeloid cells","Mast cells"))

Immune_cell$Strep_group <- ifelse(Immune_cell$Streptococcus > 0, 
                                 "Strep", NA)
Immune_cell$Strep_group_new <- ifelse(!is.na(Immune_cell$Strep_group), 
                                     "Strep", 
                                     Immune_cell$group_microbe2)

Idents(Immune_cell)<-"Strep_group_new"
table(Immune_cell$Strep_group_new)

Strep_Immune_cell_deg_all=FindMarkers(Immune_cell, ident.1 = c("Strep"), ident.2 = "Bacteria-")

dim(Strep_Immune_cell_deg_all)
write.csv(Strep_Immune_cell_deg_all,"./Strep_Immune_cell_deg_all.csv")
Strep_Immune_cell_deg_all<- read.csv("./Strep_Immune_cell_deg_all.csv",header = T)

write.csv(Strep_Immune_cell_deg_all,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/免疫细胞/Strep_Immune_cell_deg_all.csv")
Strep_Immune_cell_deg_all<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/Strep/Strep_Immune_cell_deg_all.csv",header = T)

gene_name<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/免疫细胞/免疫细胞gene_name.csv")

intersection <- intersect(gene_name$Gene.symbol, Strep_Immune_cell_deg_all$X)
library(dplyr)
Strep_Immune_cell_deg_all_inter<- Strep_Immune_cell_deg_all %>% filter(Strep_Immune_cell_deg_all$X %in%intersection )

df <- Strep_Immune_cell_deg_all

df$group2<-as.factor(ifelse(df$p_val_adj < 0.05 & abs(df$avg_log2FC) >= 0.25, 
                            ifelse(df$avg_log2FC>= 0.25 ,'Up','Down'),'NS'))
df$group2 <- factor(df$group2, levels = c("Up", "Down", "NS"))
unique(df$group2)
View(df)
##确定添加标签的数据

df$label<-ifelse(df$p_val_adj<0.05&abs(df$avg_log2FC)>=7,"Y","N")
df$label<-ifelse(df$label == 'Y', as.character(df$X), '')
write.csv(df,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/Strep/Strep_Immune_cell_deg_all数据.csv")
df<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/Strep/Strep_Immune_cell_deg_all数据.csv")
df$group2 <- factor(df$group2, levels = c("Up", "Down", "NS"))
View(df)
library(tidyverse)
library(ggrepel)

library(ggplot2) ##绘图使用
library(ggprism) ##设置主题私用
library(ggrepel)

degp<- ggplot(df, aes(x =avg_log2FC, y=-log10(p_val_adj), colour=group2)) +
  geom_point(alpha=0.85, size=3) +  #点的透明度和大小
  scale_color_manual(values=c('brown','steelblue','gray')) + #调整点的颜色
  #xlim(c(-1, 1)) +  ##调整x轴的取值范围，max(abs(BRCA_Match_DEG$logFC))最大值是多少，再进行取舍
  geom_vline(xintercept=c(-0.25,0.25),lty=4,col="black",lwd=0.8) + #添加x轴辅助线,lty函数调整线的类型
  geom_hline(yintercept = -log10(0.05), lty=4,col="black",lwd=0.8) +  #添加y轴辅助线
  labs(x="avg_log2FC", y="-log10(p_val_adj)") +  #x、y轴标签
  ggtitle("Immune cell DEG (Strep+ vs Bacteria-)") + #标题
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
ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/Strep/Strep_免疫细胞火山图.pdf",degp,width = 10,height=7)

###################################
Strep_Immune_cell_deg_all<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/Strep/Strep_Immune_cell_deg_all.csv",header = T)

library(org.Hs.eg.db)#物种注释包(Homo sapiens)
library(clusterProfiler)#富集分析
#添加entrez ID列：
##：symbol转entrez ID：
symbol <- Strep_Immune_cell_deg_all$X
entrez <- bitr(symbol,
               fromType = "SYMBOL",#现有的ID类型
               toType = "ENTREZID",#需转换的ID类型
               OrgDb = "org.Hs.eg.db")
head(entrez)
genelist <- Strep_Immune_cell_deg_all$avg_log2FC
names(genelist) <-Strep_Immune_cell_deg_all$X
head(genelist)
#genelist过滤(ID转换中丢失的部分基因)：
genelist <- genelist[names(genelist) %in% entrez[,1]]
names(genelist) <- entrez[match(names(genelist),entrez[,1]),2]
length(genelist)#过滤后剩16161个基因
head(genelist)
#将genelist按照log2FC值从高到低进行排序：
genelist <- sort(genelist,decreasing = T)
head(genelist)

#GSEA_KEGG富集分析：
KEGG_ges <- gseKEGG(
  geneList = genelist,
  organism = "hsa",#不同物种选择可官方查询：https://www.genome.jp/kegg/catalog/org_list.html
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose = FALSE,
  eps = 0
)

####
####no term enriched under specific pvalueCutoff...
####There were 16 warnings (use warnings() to see them)


#提取结果表：
KEGG_ges_result <- KEGG_ges@result

#显示前10列（第十一列为富集到该基因集的核心基因列表，这里过长不展示）：
KEGG_ges_result[1:3,1:10]
library(ggridges)
library(ggplot2)
library(enrichplot)
nrow(KEGG_ges_result)
p <- ridgeplot(KEGG_ges,
               showCategory = 15,
               fill = "p.adjust",
               decreasing  = T)
p

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/Fuso/Fuso_免疫细胞/Fuso_ridgeplot_KEGG.pdf",p,width=7,height=6)

####KEGG 无结果

##GSEA_GO富集分析：
GO_ges <- gseGO(geneList = genelist,
                OrgDb = org.Hs.eg.db,
                ont = "CC", #one of "BP", "MF", and "CC" subontologies, or "ALL" for all three.
                minGSSize = 10,
                maxGSSize = 500,
                pvalueCutoff = 0.05,
                eps = 0,
                verbose = FALSE)
GO_ges[1:3,1:10] #展示同样省略最后一列
nrow(GO_ges)

p <- ridgeplot(GO_ges,
               showCategory = 29,
               fill = "p.adjust",
               decreasing  = T)
p

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/Strep/strep_免疫细胞/Strep_ridgeplot_GO.pdf",p,width=7,height=7)


##GSEA_Reactome富集分析:
BiocManager::install("ReactomePA")
library(ReactomePA)
ret_ges <- gsePathway(genelist,
                      organism = "human",
                      minGSSize = 10,
                      maxGSSize = 500,
                      pvalueCutoff = 0.05,
                      pAdjustMethod = "BH",
                      verbose = FALSE,
                      eps = 0)
ret_ges_result<-ret_ges[1:34,1:10]
nrow(ret_ges)
p <- ridgeplot(ret_ges,
               showCategory = 34,
               fill = "p.adjust",
               decreasing  = T)
p
?ridgeplot

pp <- dotplot(ret_ges,
              showCategory = 34)
pp

####可视化柱状图
library(stringr)
library(ggplot2)
ret_ges_result <- ret_ges_result[order(ret_ges_result$NES),]
ret_ges_result$yax <- ifelse(ret_ges_result$NES >0, -0.02, 0.02)
ret_ges_result$col <- ifelse(ret_ges_result$NES > 0, "blue","red")
ret_ges_result$NES
names(ret_ges_result)

p <- 
  ggplot(ret_ges_result,aes(y=NES,x=reorder(Description,NES),label = Description))+
  geom_text(aes(y = yax,hjust = NES > 0),size=4)+ 
  geom_bar(stat="identity",aes(fill=col))+
  geom_hline(yintercept = 0)+
  coord_flip(ylim = c(-2,2.5))+
  labs(y = "Normalised enrichment score", x = "",title="")+
  scale_fill_manual(values = c("#FB6542","#375E97"),name="Expression",labels=c("Up","Down"))+
  scale_x_discrete(breaks = NULL)+
  theme_classic(base_size =  12)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.y = element_blank(),
        legend.position = 'bottom')
p


ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/Fuso/Fuso_免疫细胞/Fuso_barplot_ReactomePA.pdf",p,width=7,height=9)


##GSEA_DO(Disease Ontology)富集分析:
library(DOSE)
DO_ges <- gseDO(genelist,
                minGSSize = 10,
                maxGSSize = 500,
                pvalueCutoff = 0.05,
                pAdjustMethod = "BH",
                verbose = FALSE,
                eps = 0)
DO_ges[1:3,1:10]
nrow(DO_ges)

p <- ridgeplot(DO_ges,
               showCategory = 10,
               fill = "p.adjust",
               decreasing  = T)
p


#####GSEA_MSigDB富集分析：
library(gggsea)
library(ggplot2)
library(GSVA)
library(enrichplot)
library(GSEABase)
library(fgsea)
library(clusterProfiler)
gmtfile <- "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/h.all.v2024.1.Hs.symbols.gmt"
hallmark <- read.gmt(gmtfile)
alldiff <- Strep_Immune_cell_deg_all[order(Strep_Immune_cell_deg_all$avg_log2FC,decreasing = T),]
id <- alldiff$avg_log2FC
names(id) <- alldiff$X

## 这里去掉了基因集前缀
hallmark$term <- gsub('HALLMARK_','',hallmark$term)
hallmark.list <- hallmark %>% split(.$term) %>% lapply( "[[", 2)

gsea.re1<-clusterProfiler::GSEA(id,TERM2GENE=hallmark,verbose = T)
g1<-as.data.frame(gsea.re1)
g1<-subset(g1,p.adjust<0.05)
g1<-g1[order(g1$NES,decreasing = T),]
type(g1)
head(g1)

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


library(gggsea)

library(ggsci)

names(hallmark.list)
se_hall<-c(head(g2$pathway,9),tail(g1$pathway,9))

#3 所选中的基因集
sig1<-g2[g2$pathway%in%se_hall,]

hallmark.se<-hallmark.list[sig1$pathway]
df.new <- gseaCurve(id, hallmark.se, sig1)

pal_line<-pal_lancet()(3)#各曲线颜色

p_nes<-ggplot() + 
  geom_gsea(df.new,
            prettyGSEA=T,
            tickcolor='grey30',#竖线颜色
            ticksize=0.1,#线条粗细
            #colour='grey',#色块边框
            linecolor=pal_line,
            linesize=2,lty=1,
            #多个图形时使用，设定每行列图形个数
            nrow = 3,
            ncol =1
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

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/Fuso/Fuso_免疫细胞/Fuso_hallmark.pdf",p_nes,width=4,height=6)



####Strep分组火山图----

file_path = "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/Strep/免疫细胞分组/"
cell_type_list <- c("Pre_NMPR","Pre_MPR","Post_MPR")
dir.create(file_path)
for (cell in cell_type_list){
  print(cell)
  Idents(Immune_cell)<-"MPR_Response2"
  sub_object=subset(Immune_cell,idents = cell)
  Idents(sub_object)<-"Strep_group_new"
  deg_all=FindMarkers(sub_object, ident.1 = c("Strep"), ident.2 = "Bacteria-")
  new_filename <- paste(cell, "_deg_all",".csv", sep="")
  print(new_filename)
  write.csv(deg_all,paste("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/Strep/免疫细胞分组/",new_filename,sep=""))
}

# 获取所有csv文件列表

files <- list.files(path = file_path, pattern = "\\.csv$", recursive = TRUE)

files <- list.files(path = file_path, pattern = "\\.csv$", recursive = TRUE)

# 创建一个空数据框用于存储合并后的数据
merged_data <- data.frame()
# 循环读取每个csv文件并按照某一列降序排序，进行阈值筛选，并合并到merged_data中
for (file in files) {
  data <- read.csv(paste(file_path,file,sep=""))
  sorted_data <- data[order(-data$p_val_adj), ]
  file<-basename(file)
  print(file)
  file_name <- rep(substr(file, 1, nchar(file)-12), nrow(sorted_data))
  #filtered_data <- sorted_data[sorted_data$pct.1 >= 0.25, ]
  sorted_data$group <- file_name
  merged_data <- rbind(merged_data, sorted_data)
}

# 导出合并后的数据为新的csv文件
write.csv(merged_data, "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/Strep/免疫细胞分组/immune_cell_merged_data.csv", row.names = FALSE)


###Capnocytophaga和total bacteria- 有无微生物差异----

miMPRobe_data<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/PDCD1每个细胞以及微生物分组表达情况/16s_filter.combined.genus.umi.matrix.csv")
metadata<- FetchData(MPRNMPR_object_miMPRobe_remove,"cell_id")
metadata$cell_id <- rownames(metadata)
metadata<- left_join(x=metadata,y=miMPRobe_data,by = join_by(cell_id==barcode))
rownames(metadata)<-metadata$cell_id
MPRNMPR_object_miMPRobe_remove<- AddMetaData(MPRNMPR_object_miMPRobe_remove,metadata = metadata)
table(rownames(MPRNMPR_object_miMPRobe_remove@meta.data) %in% miMPRobe_data$barcode)
MPRNMPR_object_miMPRobe_remove@meta.data$MPR_Response2
MPRNMPR_object_miMPRobe_remove$group_microbe2
MPRNMPR_object_miMPRobe_remove$Strep_group<-MPRNMPR_object_miMPRobe_remove$group_microbe2
Idents(MPRNMPR_object_miMPRobe_remove)<-"cell_type_new"

#####Capno 免疫细胞----

Immune_cell=subset(MPRNMPR_object_miMPRobe_remove,idents = c("T cells","B cells",
                                                             "Plasma cells","Myeloid cells","Mast cells"))

Immune_cell$Capno_group <- ifelse(Immune_cell$Capnocytophaga > 0, 
                                  "Capno", NA)
Immune_cell$Capno_group_new <- ifelse(!is.na(Immune_cell$Capno_group), 
                                      "Capno", 
                                      Immune_cell$group_microbe2)

Idents(Immune_cell)<-"Capno_group_new"
table(Immune_cell$Capno_group_new)

Capno_Immune_cell_deg_all=FindMarkers(Immune_cell, ident.1 = c("Capno"), ident.2 = "Bacteria-")

dim(Capno_Immune_cell_deg_all)

write.csv(Capno_Immune_cell_deg_all,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/Capno/Capno_Immune_cell_deg_all.csv")
Capno_Immune_cell_deg_all<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/Capno/Capno_Immune_cell_deg_all.csv",header = T)

gene_name<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/免疫细胞/免疫细胞gene_name.csv")

intersection <- intersect(gene_name$Gene.symbol, Capno_Immune_cell_deg_all$X)
library(dplyr)
Capno_Immune_cell_deg_all_inter<- Capno_Immune_cell_deg_all %>% filter(Capno_Immune_cell_deg_all$X %in%intersection )

df<-Capno_Immune_cell_deg_all_inter

df <- Capno_Immune_cell_deg_all

df$group2<-as.factor(ifelse(df$p_val_adj < 0.05 & abs(df$avg_log2FC) >= 0.25, 
                            ifelse(df$avg_log2FC>= 0.25 ,'Up','Down'),'NS'))
df$group2 <- factor(df$group2, levels = c("Up", "Down", "NS"))
unique(df$group2)

##确定添加标签的数据

df$label<-ifelse(df$p_val_adj<0.05&abs(df$avg_log2FC)>=0.25,"Y","N")
df$label<-ifelse(df$label == 'Y', as.character(df$X), '')
write.csv(df,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/Capno/Capno_Immune_cell_deg_all数据.csv")
df<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/Capno/Capno_Immune_cell_deg_all数据.csv")
df$group2 <- factor(df$group2, levels = c("Up", "Down", "NS"))
View(df)
library(tidyverse)
library(ggrepel)
library(ggplot2) ##绘图使用
library(ggprism) ##设置主题私用
library(ggrepel)

degp<- ggplot(df, aes(x =avg_log2FC, y=-log10(p_val_adj), colour=group2)) +
  geom_point(alpha=0.85, size=3) +  #点的透明度和大小
  scale_color_manual(values=c('brown','steelblue','gray')) + #调整点的颜色
  #xlim(c(-1, 1)) +  ##调整x轴的取值范围，max(abs(BRCA_Match_DEG$logFC))最大值是多少，再进行取舍
  geom_vline(xintercept=c(-0.25,0.25),lty=4,col="black",lwd=0.8) + #添加x轴辅助线,lty函数调整线的类型
  geom_hline(yintercept = -log10(0.05), lty=4,col="black",lwd=0.8) +  #添加y轴辅助线
  labs(x="avg_log2FC", y="-log10(p_val_adj)") +  #x、y轴标签
  ggtitle("Immune cell DEG (Capno+ vs Bacteria-)") + #标题
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
ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/Capno/Capno_免疫细胞/Capno_免疫细胞火山图.pdf",degp,width = 7,height=6)


###################################
library(org.Hs.eg.db)#物种注释包(Homo sapiens)
library(clusterProfiler)#富集分析
#添加entrez ID列：
##：symbol转entrez ID：
symbol <- Capno_Immune_cell_deg_all$X
entrez <- bitr(symbol,
               fromType = "SYMBOL",#现有的ID类型
               toType = "ENTREZID",#需转换的ID类型
               OrgDb = "org.Hs.eg.db")
head(entrez)
genelist <- Capno_Immune_cell_deg_all$avg_log2FC
names(genelist) <-Capno_Immune_cell_deg_all$X
head(genelist)
#genelist过滤(ID转换中丢失的部分基因)：
genelist <- genelist[names(genelist) %in% entrez[,1]]
names(genelist) <- entrez[match(names(genelist),entrez[,1]),2]
length(genelist)#过滤后剩16161个基因
head(genelist)
#将genelist按照log2FC值从高到低进行排序：
genelist <- sort(genelist,decreasing = T)
head(genelist)

#GSEA_KEGG富集分析：
KEGG_ges <- gseKEGG(
  geneList = genelist,
  organism = "hsa",#不同物种选择可官方查询：https://www.genome.jp/kegg/catalog/org_list.html
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose = FALSE,
  eps = 0
)

#提取结果表：
KEGG_ges_result <- KEGG_ges@result

#显示前10列（第十一列为富集到该基因集的核心基因列表，这里过长不展示）：
KEGG_ges_result[1:3,1:10]
library(ggridges)
library(ggplot2)
library(enrichplot)
nrow(KEGG_ges_result)
p <- ridgeplot(KEGG_ges,
               showCategory = 15,
               fill = "p.adjust",
               decreasing  = T)
p

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/Capno/Capno_免疫细胞/Capno_ridgeplot_KEGG.pdf",p,width=7,height=6)


##GSEA_GO富集分析：
GO_ges <- gseGO(geneList = genelist,
                OrgDb = org.Hs.eg.db,
                ont = "CC", #one of "BP", "MF", and "CC" subontologies, or "ALL" for all three.
                minGSSize = 10,
                maxGSSize = 500,
                pvalueCutoff = 0.05,
                eps = 0,
                verbose = FALSE)
GO_ges[1:3,1:10] #展示同样省略最后一列
nrow(GO_ges)

p <- ridgeplot(GO_ges,
               showCategory = 36,
               fill = "p.adjust",
               decreasing  = T)
p

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/Capno/Capno_免疫细胞/Capno_ridgeplot_GO.pdf",p,width=7,height=9)


##GSEA_Reactome富集分析:
BiocManager::install("ReactomePA")
library(ReactomePA)
ret_ges <- gsePathway(genelist,
                      organism = "human",
                      minGSSize = 10,
                      maxGSSize = 500,
                      pvalueCutoff = 0.05,
                      pAdjustMethod = "BH",
                      verbose = FALSE,
                      eps = 0)
ret_ges_result<-ret_ges[1:34,1:10]
nrow(ret_ges)
p <- ridgeplot(ret_ges,
               showCategory = 10,
               fill = "p.adjust",
               decreasing  = T)
p
?ridgeplot

pp <- dotplot(ret_ges,
              showCategory = 34)
pp

####可视化柱状图
library(stringr)
library(ggplot2)
ret_ges_result <- ret_ges_result[order(ret_ges_result$NES),]
ret_ges_result$yax <- ifelse(ret_ges_result$NES >0, -0.02, 0.02)
ret_ges_result$col <- ifelse(ret_ges_result$NES > 0, "blue","red")
ret_ges_result$NES
names(ret_ges_result)

p <- 
  ggplot(ret_ges_result,aes(y=NES,x=reorder(Description,NES),label = Description))+
  geom_text(aes(y = yax,hjust = NES > 0),size=4)+ 
  geom_bar(stat="identity",aes(fill=col))+
  geom_hline(yintercept = 0)+
  coord_flip(ylim = c(-2,2.5))+
  labs(y = "Normalised enrichment score", x = "",title="")+
  scale_fill_manual(values = c("#FB6542","#375E97"),name="Expression",labels=c("Up","Down"))+
  scale_x_discrete(breaks = NULL)+
  theme_classic(base_size =  12)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.y = element_blank(),
        legend.position = 'bottom')
p


ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/Capno/Capno_免疫细胞/Capno_barplot_ReactomePA.pdf",p,width=7,height=9)


##GSEA_DO(Disease Ontology)富集分析:
library(DOSE)
DO_ges <- gseDO(genelist,
                minGSSize = 10,
                maxGSSize = 500,
                pvalueCutoff = 0.05,
                pAdjustMethod = "BH",
                verbose = FALSE,
                eps = 0)
DO_ges[1:3,1:10]
nrow(DO_ges)

p <- ridgeplot(DO_ges,
               showCategory = 10,
               fill = "p.adjust",
               decreasing  = T)
p


#####GSEA_MSigDB富集分析：
library(gggsea)
library(ggplot2)
library(GSVA)
library(enrichplot)
library(GSEABase)
library(fgsea)
library(clusterProfiler)
gmtfile <- "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/h.all.v2024.1.Hs.symbols.gmt"
hallmark <- read.gmt(gmtfile)
alldiff <- Capno_Immune_cell_deg_all[order(Capno_Immune_cell_deg_all$avg_log2FC,decreasing = T),]
id <- alldiff$avg_log2FC
names(id) <- alldiff$X

## 这里去掉了基因集前缀
hallmark$term <- gsub('HALLMARK_','',hallmark$term)
hallmark.list <- hallmark %>% split(.$term) %>% lapply( "[[", 2)

gsea.re1<-clusterProfiler::GSEA(id,TERM2GENE=hallmark,verbose = T)
g1<-as.data.frame(gsea.re1)
g1<-subset(g1,p.adjust<0.05)
g1<-g1[order(g1$NES,decreasing = T),]
type(g1)
head(g1)

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


library(gggsea)

library(ggsci)

names(hallmark.list)
se_hall<-c(head(g2$pathway,9),tail(g1$pathway,9))

#3 所选中的基因集
sig1<-g2[g2$pathway%in%se_hall,]

hallmark.se<-hallmark.list[sig1$pathway]
df.new <- gseaCurve(id, hallmark.se, sig1)

pal_line<-pal_lancet()(3)#各曲线颜色

p_nes<-ggplot() + 
  geom_gsea(df.new,
            prettyGSEA=T,
            tickcolor='grey30',#竖线颜色
            ticksize=0.1,#线条粗细
            #colour='grey',#色块边框
            linecolor=pal_line,
            linesize=2,lty=1,
            #多个图形时使用，设定每行列图形个数
            nrow = 2,
            ncol =4
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

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/Capno/Capno_免疫细胞/Capno_hallmark.pdf",p_nes,width=10,height=4)


##上皮细胞有无微生物差异分析----

######Fuso epi----

miMPRobe_data<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/PDCD1每个细胞以及微生物分组表达情况/16s_filter.combined.genus.umi.matrix.csv")

metadata<- FetchData(MPRNMPR_object_miMPRobe_remove,"cell_id")
metadata$cell_id <- rownames(metadata)
metadata<- left_join(x=metadata,y=miMPRobe_data,by = join_by(cell_id==barcode))
rownames(metadata)<-metadata$cell_id
MPRNMPR_object_miMPRobe_remove<- AddMetaData(MPRNMPR_object_miMPRobe_remove,metadata = metadata)
table(rownames(MPRNMPR_object_miMPRobe_remove@meta.data) %in% miMPRobe_data$barcode)

MPRNMPR_object_miMPRobe_remove$fuso_group<-MPRNMPR_object_miMPRobe_remove$group_microbe2
Idents(MPRNMPR_object_miMPRobe_remove)<-"cell_type_new"
epi_cell=subset(MPRNMPR_object_miMPRobe_remove,idents = c("Epithelial cells"))

epi_cell$fuso_group <- ifelse(epi_cell$Fusobacterium > 0, 
                              "Fuso", NA)
epi_cell$fuso_group_new <- ifelse(!is.na(epi_cell$fuso_group), 
                                  "Fuso", 
                                  epi_cell$group_microbe2)

Idents(epi_cell)<-"fuso_group_new"
table(epi_cell$fuso_group_new)

Fuso_epi_cell_deg_all=FindMarkers(epi_cell, ident.1 = c("Fuso"), ident.2 = "Bacteria-")

dim(Fuso_epi_cell_deg_all)
write.csv(Fuso_epi_cell_deg_all,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/上皮细胞/上皮细胞Fuso/Fuso_epi_cell_deg_all.csv")
Fuso_epi_cell_deg_all<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/上皮细胞/上皮细胞Fuso/Fuso_epi_cell_deg_all.csv",header = T)



gene_name<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/上皮细胞/上皮细胞gene_name.csv")

intersection <- intersect(gene_name$Gene.symbol, Fuso_epi_cell_deg_all$X)

Fuso_epi_cell_deg_all_filter_inter<- Fuso_epi_cell_deg_all %>% filter(Fuso_epi_cell_deg_all$X %in%intersection )

df <- Fuso_epi_cell_deg_all_filter_inter
df <-Fuso_epi_cell_deg_all

df$group2<-as.factor(ifelse(df$p_val_adj < 0.05 & abs(df$avg_log2FC) >= 0.25, 
                            ifelse(df$avg_log2FC>= 0.25 ,'Up','Down'),'NS'))
df$group2 <- factor(df$group2, levels = c("Up", "Down", "NS"))


df$group2<-as.factor(ifelse(df$p_val_adj < 0.05 & abs(df$avg_log2FC) >= 0.25, 
                            ifelse(df$avg_log2FC>= 0.25 ,'Up','Down'),'NS'))
df$group2 <- factor(df$group2, levels = c("Up", "Down", "NS"))
unique(df$group2)
View(df)
##确定添加标签的数据

df$label<-ifelse(df$p_val_adj<0.05&abs(df$avg_log2FC)>=0.25,"Y","N")
df$label<-ifelse(df$label == 'Y', as.character(df$X), '')
df$group2 <- factor(df$group2, levels = c("Up", "Down", "NS"))
View(df)
library(tidyverse)
library(ggrepel)

library(ggplot2) ##绘图使用
library(ggprism) ##设置主题私用
library(ggrepel)

degp<- ggplot(df, aes(x =avg_log2FC, y=-log10(p_val_adj), colour=group2)) +
  geom_point(alpha=0.85, size=3) +  #点的透明度和大小
  scale_color_manual(values=c('brown','steelblue','gray')) + #调整点的颜色
  #xlim(c(-1, 1)) +  ##调整x轴的取值范围，max(abs(BRCA_Match_DEG$logFC))最大值是多少，再进行取舍
  geom_vline(xintercept=c(-0.25,0.25),lty=4,col="black",lwd=0.8) + #添加x轴辅助线,lty函数调整线的类型
  geom_hline(yintercept = -log10(0.05), lty=4,col="black",lwd=0.8) +  #添加y轴辅助线
  labs(x="avg_log2FC", y="-log10(p_val_adj)") +  #x、y轴标签
  ggtitle("Epithelial cells DEG (Fuso+ vs Bacteria-)") + #标题
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

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/上皮细胞/上皮细胞Fuso/Fuso上皮细胞火山图.pdf",degp,width=8,height=7)

####Fuso GESA 分析----

library(org.Hs.eg.db)#物种注释包(Homo sapiens)
library(clusterProfiler)#富集分析
#添加entrez ID列：
##：symbol转entrez ID：
symbol <- Fuso_epi_cell_deg_all$X
entrez <- bitr(symbol,
               fromType = "SYMBOL",#现有的ID类型
               toType = "ENTREZID",#需转换的ID类型
               OrgDb = "org.Hs.eg.db")
head(entrez)
genelist <- Fuso_epi_cell_deg_all$avg_log2FC
names(genelist) <-Fuso_epi_cell_deg_all$X
head(genelist)
#genelist过滤(ID转换中丢失的部分基因)：
genelist <- genelist[names(genelist) %in% entrez[,1]]
names(genelist) <- entrez[match(names(genelist),entrez[,1]),2]
length(genelist)#过滤后剩16161个基因
head(genelist)
#将genelist按照log2FC值从高到低进行排序：
genelist <- sort(genelist,decreasing = T)
head(genelist)

#GSEA_KEGG富集分析：
KEGG_ges <- gseKEGG(
  geneList = genelist,
  organism = "hsa",#不同物种选择可官方查询：https://www.genome.jp/kegg/catalog/org_list.html
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose = FALSE,
  eps = 0
)

#提取结果表：
KEGG_ges_result <- KEGG_ges@result

#显示前10列（第十一列为富集到该基因集的核心基因列表，这里过长不展示）：
KEGG_ges_result[1:3,1:10]
library(ggridges)
library(ggplot2)
library(enrichplot)
nrow(KEGG_ges_result)
p <- ridgeplot(KEGG_ges,
               showCategory = 15,
               fill = "p.adjust",
               decreasing  = T)
p

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/上皮细胞/上皮细胞Fuso/Fuso_ridgeplot_KEGG.pdf",p,width=7,height=6)


##GSEA_GO富集分析：
GO_ges <- gseGO(geneList = genelist,
                OrgDb = org.Hs.eg.db,
                ont = "CC", #one of "BP", "MF", and "CC" subontologies, or "ALL" for all three.
                minGSSize = 10,
                maxGSSize = 500,
                pvalueCutoff = 0.05,
                eps = 0,
                verbose = FALSE)
GO_ges[1:3,1:10] #展示同样省略最后一列
nrow(GO_ges)

p <- ridgeplot(GO_ges,
               showCategory = 29,
               fill = "p.adjust",
               decreasing  = T)
p

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/上皮细胞/上皮细胞Fuso/Fuso_ridgeplot_GO.pdf",p,width=7,height=9)


##GSEA_Reactome富集分析:
BiocManager::install("ReactomePA")
library(ReactomePA)
ret_ges <- gsePathway(genelist,
                      organism = "human",
                      minGSSize = 10,
                      maxGSSize = 500,
                      pvalueCutoff = 0.05,
                      pAdjustMethod = "BH",
                      verbose = FALSE,
                      eps = 0)
ret_ges_result<-ret_ges[1:34,1:10]
nrow(ret_ges)
p <- ridgeplot(ret_ges,
               showCategory = 34,
               fill = "p.adjust",
               decreasing  = T)
p
?ridgeplot

pp <- dotplot(ret_ges,
              showCategory = 34)
pp

####可视化柱状图
library(stringr)
library(ggplot2)
ret_ges_result <- ret_ges_result[order(ret_ges_result$NES),]
ret_ges_result$yax <- ifelse(ret_ges_result$NES >0, -0.02, 0.02)
ret_ges_result$col <- ifelse(ret_ges_result$NES > 0, "blue","red")
ret_ges_result$NES
names(ret_ges_result)

p <- 
  ggplot(ret_ges_result,aes(y=NES,x=reorder(Description,NES),label = Description))+
  geom_text(aes(y = yax,hjust = NES > 0),size=4)+ 
  geom_bar(stat="identity",aes(fill=col))+
  geom_hline(yintercept = 0)+
  coord_flip(ylim = c(-2,2.5))+
  labs(y = "Normalised enrichment score", x = "",title="")+
  scale_fill_manual(values = c("#FB6542","#375E97"),name="Expression",labels=c("Up","Down"))+
  scale_x_discrete(breaks = NULL)+
  theme_classic(base_size =  12)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.y = element_blank(),
        legend.position = 'bottom')
p



##GSEA_DO(Disease Ontology)富集分析:
library(DOSE)
DO_ges <- gseDO(genelist,
                minGSSize = 10,
                maxGSSize = 500,
                pvalueCutoff = 0.05,
                pAdjustMethod = "BH",
                verbose = FALSE,
                eps = 0)
DO_ges[1:3,1:10]
nrow(DO_ges)

p <- ridgeplot(DO_ges,
               showCategory = 10,
               fill = "p.adjust",
               decreasing  = T)
p


#####GSEA_MSigDB富集分析：
library(gggsea)
library(ggplot2)
library(GSVA)
library(enrichplot)
library(GSEABase)
library(fgsea)
library(clusterProfiler)
gmtfile <- "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/h.all.v2024.1.Hs.symbols.gmt"
hallmark <- read.gmt(gmtfile)
alldiff <- Fuso_epi_cell_deg_all[order(Fuso_epi_cell_deg_all$avg_log2FC,decreasing = T),]
id <- alldiff$avg_log2FC
names(id) <- alldiff$X

## 这里去掉了基因集前缀
hallmark$term <- gsub('HALLMARK_','',hallmark$term)
hallmark.list <- hallmark %>% split(.$term) %>% lapply( "[[", 2)

gsea.re1<-clusterProfiler::GSEA(id,TERM2GENE=hallmark,verbose = T)
g1<-as.data.frame(gsea.re1)
g1<-subset(g1,p.adjust<0.05)
g1<-g1[order(g1$NES,decreasing = T),]
type(g1)
head(g1)

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


library(gggsea)

library(ggsci)

names(hallmark.list)
se_hall<-c(head(g2$pathway,9),tail(g1$pathway,9))

#3 所选中的基因集
sig1<-g2[g2$pathway%in%se_hall,]

hallmark.se<-hallmark.list[sig1$pathway]
df.new <- gseaCurve(id, hallmark.se, sig1)

pal_line<-pal_lancet()(7)#各曲线颜色

p_nes<-ggplot() + 
  geom_gsea(df.new,
            prettyGSEA=T,
            tickcolor='grey30',#竖线颜色
            ticksize=0.1,#线条粗细
            #colour='grey',#色块边框
            linecolor=pal_line,
            linesize=2,lty=1,
            #多个图形时使用，设定每行列图形个数
            nrow = 2,
            ncol =4
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

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/上皮细胞/上皮细胞Fuso/Fuso_hallmark.pdf",p_nes,width=10,height=4)

####strep epi----

miMPRobe_data<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/PDCD1每个细胞以及微生物分组表达情况/16s_filter.combined.genus.umi.matrix.csv")

metadata<- FetchData(MPRNMPR_object_miMPRobe_remove,"cell_id")
metadata$cell_id <- rownames(metadata)
metadata<- left_join(x=metadata,y=miMPRobe_data,by = join_by(cell_id==barcode))
rownames(metadata)<-metadata$cell_id
MPRNMPR_object_miMPRobe_remove<- AddMetaData(MPRNMPR_object_miMPRobe_remove,metadata = metadata)
table(rownames(MPRNMPR_object_miMPRobe_remove@meta.data) %in% miMPRobe_data$barcode)

MPRNMPR_object_miMPRobe_remove$fuso_group<-MPRNMPR_object_miMPRobe_remove$group_microbe2
Idents(MPRNMPR_object_miMPRobe_remove)<-"cell_type_new"
epi_cell=subset(MPRNMPR_object_miMPRobe_remove,idents = c("Epithelial cells"))

epi_cell$strep_group <- ifelse(epi_cell$Streptococcus > 0, 
                              "strep", NA)
epi_cell$strep_group_new <- ifelse(!is.na(epi_cell$strep_group), 
                                  "Strep", 
                                  epi_cell$group_microbe2)

Idents(epi_cell)<-"strep_group_new"
table(epi_cell$strep_group_new)

Strep_epi_cell_deg_all=FindMarkers(epi_cell, ident.1 = c("Strep"), ident.2 = "Bacteria-")

dim(Strep_epi_cell_deg_all)
write.csv(Strep_epi_cell_deg_all,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/上皮细胞/上皮细胞Strep/Strep_epi_cell_deg_all.csv")
Strep_epi_cell_deg_all<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/上皮细胞/上皮细胞Strep/Strep_epi_cell_deg_all.csv",header = T)


gene_name<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/上皮细胞/上皮细胞gene_name.csv")

intersection <- intersect(gene_name$Gene.symbol, Strep_epi_cell_deg_all$X)

Strep_epi_cell_deg_all_inter<- Strep_epi_cell_deg_all %>% filter(Strep_epi_cell_deg_all$X %in%intersection )

df <- Strep_epi_cell_deg_all_inter

df$group2<-as.factor(ifelse(df$p_val_adj < 0.05 & abs(df$avg_log2FC) >= 0.25, 
                            ifelse(df$avg_log2FC>= 0.25 ,'Up','Down'),'NS'))
df$group2 <- factor(df$group2, levels = c("Up", "Down", "NS"))
unique(df$group2)
head(df)
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
  #xlim(c(-1, 1)) +  ##调整x轴的取值范围，max(abs(BRCA_Match_DEG$logFC))最大值是多少，再进行取舍
  geom_vline(xintercept=c(-0.25,0.25),lty=4,col="black",lwd=0.8) + #添加x轴辅助线,lty函数调整线的类型
  geom_hline(yintercept = -log10(0.05), lty=4,col="black",lwd=0.8) +  #添加y轴辅助线
  labs(x="avg_log2FC", y="-log10(p_val_adj)") +  #x、y轴标签
  ggtitle("Epithelial cells DEG (Strep+ vs Bacteria-)") + #标题
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

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/上皮细胞/上皮细胞Strep/上皮细胞火山图.pdf",degp,width = 7,height=7)


####Strep GESA 分析----
library(org.Hs.eg.db)#物种注释包(Homo sapiens)
library(clusterProfiler)#富集分析
#添加entrez ID列：
##：symbol转entrez ID：
symbol <- Strep_epi_cell_deg_all$X
entrez <- bitr(symbol,
               fromType = "SYMBOL",#现有的ID类型
               toType = "ENTREZID",#需转换的ID类型
               OrgDb = "org.Hs.eg.db")
head(entrez)
genelist <- Strep_epi_cell_deg_all$avg_log2FC
names(genelist) <-Strep_epi_cell_deg_all$X
head(genelist)
#genelist过滤(ID转换中丢失的部分基因)：
genelist <- genelist[names(genelist) %in% entrez[,1]]
names(genelist) <- entrez[match(names(genelist),entrez[,1]),2]
length(genelist)#过滤后剩16161个基因
head(genelist)
#将genelist按照log2FC值从高到低进行排序：
genelist <- sort(genelist,decreasing = T)
head(genelist)

#GSEA_KEGG富集分析：
KEGG_ges <- gseKEGG(
  geneList = genelist,
  organism = "hsa",#不同物种选择可官方查询：https://www.genome.jp/kegg/catalog/org_list.html
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose = FALSE,
  eps = 0
)

#提取结果表：
KEGG_ges_result <- KEGG_ges@result

#显示前10列（第十一列为富集到该基因集的核心基因列表，这里过长不展示）：
KEGG_ges_result[1:3,1:10]
library(ggridges)
library(ggplot2)
library(enrichplot)
nrow(KEGG_ges_result)
p <- ridgeplot(KEGG_ges,
               showCategory = 15,
               fill = "p.adjust",
               decreasing  = T)
p

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/上皮细胞/上皮细胞Strep/Strep_ridgeplot_KEGG.pdf",p,width=7,height=6)


##GSEA_GO富集分析：
GO_ges <- gseGO(geneList = genelist,
                OrgDb = org.Hs.eg.db,
                ont = "CC", #one of "BP", "MF", and "CC" subontologies, or "ALL" for all three.
                minGSSize = 10,
                maxGSSize = 500,
                pvalueCutoff = 0.05,
                eps = 0,
                verbose = FALSE)
GO_ges[1:3,1:10] #展示同样省略最后一列
nrow(GO_ges)

p <- ridgeplot(GO_ges,
               showCategory = 15,
               fill = "p.adjust",
               decreasing  = T)
p

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/上皮细胞/上皮细胞Strep/Strep_ridgeplot_GO.pdf",p,width=7,height=9)


##GSEA_Reactome富集分析:
BiocManager::install("ReactomePA")
library(ReactomePA)
ret_ges <- gsePathway(genelist,
                      organism = "human",
                      minGSSize = 10,
                      maxGSSize = 500,
                      pvalueCutoff = 0.05,
                      pAdjustMethod = "BH",
                      verbose = FALSE,
                      eps = 0)
ret_ges_result<-ret_ges[1:34,1:10]
nrow(ret_ges)
p <- ridgeplot(ret_ges,
               showCategory = 34,
               fill = "p.adjust",
               decreasing  = T)
p
?ridgeplot

pp <- dotplot(ret_ges,
              showCategory = 34)
pp

####可视化柱状图
library(stringr)
library(ggplot2)
ret_ges_result <- ret_ges_result[order(ret_ges_result$NES),]
ret_ges_result$yax <- ifelse(ret_ges_result$NES >0, -0.02, 0.02)
ret_ges_result$col <- ifelse(ret_ges_result$NES > 0, "blue","red")
ret_ges_result$NES
names(ret_ges_result)

p <- 
  ggplot(ret_ges_result,aes(y=NES,x=reorder(Description,NES),label = Description))+
  geom_text(aes(y = yax,hjust = NES > 0),size=4)+ 
  geom_bar(stat="identity",aes(fill=col))+
  geom_hline(yintercept = 0)+
  coord_flip(ylim = c(-2,2.5))+
  labs(y = "Normalised enrichment score", x = "",title="")+
  scale_fill_manual(values = c("#FB6542","#375E97"),name="Expression",labels=c("Up","Down"))+
  scale_x_discrete(breaks = NULL)+
  theme_classic(base_size =  12)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.y = element_blank(),
        legend.position = 'bottom')
p



##GSEA_DO(Disease Ontology)富集分析:
library(DOSE)
DO_ges <- gseDO(genelist,
                minGSSize = 10,
                maxGSSize = 500,
                pvalueCutoff = 0.05,
                pAdjustMethod = "BH",
                verbose = FALSE,
                eps = 0)
DO_ges[1:3,1:10]
nrow(DO_ges)

p <- ridgeplot(DO_ges,
               showCategory = 10,
               fill = "p.adjust",
               decreasing  = T)
p


#####GSEA_MSigDB富集分析：
library(gggsea)
library(ggplot2)
library(GSVA)
library(enrichplot)
library(GSEABase)
library(fgsea)
library(clusterProfiler)
gmtfile <- "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/h.all.v2024.1.Hs.symbols.gmt"
hallmark <- read.gmt(gmtfile)
alldiff <- Fuso_epi_cell_deg_all[order(Fuso_epi_cell_deg_all$avg_log2FC,decreasing = T),]
id <- alldiff$avg_log2FC
names(id) <- alldiff$X

## 这里去掉了基因集前缀
hallmark$term <- gsub('HALLMARK_','',hallmark$term)
hallmark.list <- hallmark %>% split(.$term) %>% lapply( "[[", 2)

gsea.re1<-clusterProfiler::GSEA(id,TERM2GENE=hallmark,verbose = T)
g1<-as.data.frame(gsea.re1)
g1<-subset(g1,p.adjust<0.05)
g1<-g1[order(g1$NES,decreasing = T),]
type(g1)
head(g1)

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


library(gggsea)

library(ggsci)

names(hallmark.list)
se_hall<-c(head(g2$pathway,9),tail(g1$pathway,9))

#3 所选中的基因集
sig1<-g2[g2$pathway%in%se_hall,]

hallmark.se<-hallmark.list[sig1$pathway]
df.new <- gseaCurve(id, hallmark.se, sig1)

pal_line<-pal_lancet()(7)#各曲线颜色

p_nes<-ggplot() + 
  geom_gsea(df.new,
            prettyGSEA=T,
            tickcolor='grey30',#竖线颜色
            ticksize=0.1,#线条粗细
            #colour='grey',#色块边框
            linecolor=pal_line,
            linesize=2,lty=1,
            #多个图形时使用，设定每行列图形个数
            nrow = 2,
            ncol =4
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

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/上皮细胞/上皮细胞Strep/Strep_hallmark.pdf",p_nes,width=10,height=4)


####Capnocytophaga epi----

miMPRobe_data<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/PDCD1每个细胞以及微生物分组表达情况/16s_filter.combined.genus.umi.matrix.csv")

metadata<- FetchData(MPRNMPR_object_miMPRobe_remove,"cell_id")
metadata$cell_id <- rownames(metadata)
metadata<- left_join(x=metadata,y=miMPRobe_data,by = join_by(cell_id==barcode))
rownames(metadata)<-metadata$cell_id
MPRNMPR_object_miMPRobe_remove<- AddMetaData(MPRNMPR_object_miMPRobe_remove,metadata = metadata)
table(rownames(MPRNMPR_object_miMPRobe_remove@meta.data) %in% miMPRobe_data$barcode)

MPRNMPR_object_miMPRobe_remove$fuso_group<-MPRNMPR_object_miMPRobe_remove$group_microbe2
Idents(MPRNMPR_object_miMPRobe_remove)<-"cell_type_new"
epi_cell=subset(MPRNMPR_object_miMPRobe_remove,idents = c("Epithelial cells"))

epi_cell$Capno_group <- ifelse(epi_cell$Capnocytophaga > 0, 
                               "Capno", NA)
epi_cell$Capno_group_new <- ifelse(!is.na(epi_cell$Capno_group), 
                                   "Capno", 
                                   epi_cell$group_microbe2)

Idents(epi_cell)<-"Capno_group_new"
table(epi_cell$Capno_group_new)

Capno_epi_cell_deg_all=FindMarkers(epi_cell, ident.1 = c("Capno"), ident.2 = "Bacteria-")

dim(Capno_epi_cell_deg_all)
write.csv(Capno_epi_cell_deg_all,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/上皮细胞/上皮细胞Capno/Capno_epi_cell_deg_all.csv")
Capno_epi_cell_deg_all<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/上皮细胞/上皮细胞Capno/Capno_epi_cell_deg_all.csv",header = T)


gene_name<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/上皮细胞/上皮细胞gene_name.csv")

intersection <- intersect(gene_name$Gene.symbol, Capno_epi_cell_deg_all$X)

Capno_epi_cell_deg_all_inter<- Capno_epi_cell_deg_all %>% filter(Capno_epi_cell_deg_all$X %in%intersection )

df <- Capno_epi_cell_deg_all_inter

df$group2<-as.factor(ifelse(df$p_val_adj < 0.05 & abs(df$avg_log2FC) >= 0.25, 
                            ifelse(df$avg_log2FC>= 0.25 ,'Up','Down'),'NS'))
df$group2 <- factor(df$group2, levels = c("Up", "Down", "NS"))
unique(df$group2)
head(df)
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
  #xlim(c(-1, 1)) +  ##调整x轴的取值范围，max(abs(BRCA_Match_DEG$logFC))最大值是多少，再进行取舍
  geom_vline(xintercept=c(-0.25,0.25),lty=4,col="black",lwd=0.8) + #添加x轴辅助线,lty函数调整线的类型
  geom_hline(yintercept = -log10(0.05), lty=4,col="black",lwd=0.8) +  #添加y轴辅助线
  labs(x="avg_log2FC", y="-log10(p_val_adj)") +  #x、y轴标签
  ggtitle("Epithelial cells DEG (Strep+ vs Bacteria-)") + #标题
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

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/上皮细胞/上皮细胞Capno/上皮细胞火山图.pdf",degp,width = 7,height=7)


####Strep GESA 分析----
library(org.Hs.eg.db)#物种注释包(Homo sapiens)
library(clusterProfiler)#富集分析
#添加entrez ID列：
##：symbol转entrez ID：
symbol <- Capno_epi_cell_deg_all$X
entrez <- bitr(symbol,
               fromType = "SYMBOL",#现有的ID类型
               toType = "ENTREZID",#需转换的ID类型
               OrgDb = "org.Hs.eg.db")
head(entrez)
genelist <- Capno_epi_cell_deg_all$avg_log2FC
names(genelist) <-Capno_epi_cell_deg_all$X
head(genelist)
#genelist过滤(ID转换中丢失的部分基因)：
genelist <- genelist[names(genelist) %in% entrez[,1]]
names(genelist) <- entrez[match(names(genelist),entrez[,1]),2]
length(genelist)#过滤后剩16161个基因
head(genelist)
#将genelist按照log2FC值从高到低进行排序：
genelist <- sort(genelist,decreasing = T)
head(genelist)

#GSEA_KEGG富集分析：
KEGG_ges <- gseKEGG(
  geneList = genelist,
  organism = "hsa",#不同物种选择可官方查询：https://www.genome.jp/kegg/catalog/org_list.html
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose = FALSE,
  eps = 0
)

#提取结果表：
KEGG_ges_result <- KEGG_ges@result

#显示前10列（第十一列为富集到该基因集的核心基因列表，这里过长不展示）：
KEGG_ges_result[1:3,1:10]
library(ggridges)
library(ggplot2)
library(enrichplot)
nrow(KEGG_ges_result)
p <- ridgeplot(KEGG_ges,
               showCategory = 15,
               fill = "p.adjust",
               decreasing  = T)
p

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/上皮细胞/上皮细胞Capno/Capno_ridgeplot_KEGG.pdf",p,width=7,height=6)


##GSEA_GO富集分析：
GO_ges <- gseGO(geneList = genelist,
                OrgDb = org.Hs.eg.db,
                ont = "CC", #one of "BP", "MF", and "CC" subontologies, or "ALL" for all three.
                minGSSize = 10,
                maxGSSize = 500,
                pvalueCutoff = 0.05,
                eps = 0,
                verbose = FALSE)
GO_ges[1:3,1:10] #展示同样省略最后一列
nrow(GO_ges)

p <- ridgeplot(GO_ges,
               showCategory = 15,
               fill = "p.adjust",
               decreasing  = T)
p

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/上皮细胞/上皮细胞Capno/Capno_ridgeplot_GO.pdf",p,width=7,height=5)


##GSEA_Reactome富集分析:
BiocManager::install("ReactomePA")
library(ReactomePA)
ret_ges <- gsePathway(genelist,
                      organism = "human",
                      minGSSize = 10,
                      maxGSSize = 500,
                      pvalueCutoff = 0.05,
                      pAdjustMethod = "BH",
                      verbose = FALSE,
                      eps = 0)
ret_ges_result<-ret_ges[1:34,1:10]
nrow(ret_ges)
p <- ridgeplot(ret_ges,
               showCategory = 34,
               fill = "p.adjust",
               decreasing  = T)
p
?ridgeplot

pp <- dotplot(ret_ges,
              showCategory = 34)
pp

####可视化柱状图
library(stringr)
library(ggplot2)
ret_ges_result <- ret_ges_result[order(ret_ges_result$NES),]
ret_ges_result$yax <- ifelse(ret_ges_result$NES >0, -0.02, 0.02)
ret_ges_result$col <- ifelse(ret_ges_result$NES > 0, "blue","red")
ret_ges_result$NES
names(ret_ges_result)

p <- 
  ggplot(ret_ges_result,aes(y=NES,x=reorder(Description,NES),label = Description))+
  geom_text(aes(y = yax,hjust = NES > 0),size=4)+ 
  geom_bar(stat="identity",aes(fill=col))+
  geom_hline(yintercept = 0)+
  coord_flip(ylim = c(-2,2.5))+
  labs(y = "Normalised enrichment score", x = "",title="")+
  scale_fill_manual(values = c("#FB6542","#375E97"),name="Expression",labels=c("Up","Down"))+
  scale_x_discrete(breaks = NULL)+
  theme_classic(base_size =  12)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.y = element_blank(),
        legend.position = 'bottom')
p



##GSEA_DO(Disease Ontology)富集分析:
library(DOSE)
DO_ges <- gseDO(genelist,
                minGSSize = 10,
                maxGSSize = 500,
                pvalueCutoff = 0.05,
                pAdjustMethod = "BH",
                verbose = FALSE,
                eps = 0)
DO_ges[1:3,1:10]
nrow(DO_ges)

p <- ridgeplot(DO_ges,
               showCategory = 10,
               fill = "p.adjust",
               decreasing  = T)
p


#####GSEA_MSigDB富集分析：
library(gggsea)
library(ggplot2)
library(GSVA)
library(enrichplot)
library(GSEABase)
library(fgsea)
library(clusterProfiler)
gmtfile <- "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/h.all.v2024.1.Hs.symbols.gmt"
hallmark <- read.gmt(gmtfile)
alldiff <- Capno_epi_cell_deg_all[order(Capno_epi_cell_deg_all$avg_log2FC,decreasing = T),]
id <- alldiff$avg_log2FC
names(id) <- alldiff$X

## 这里去掉了基因集前缀
hallmark$term <- gsub('HALLMARK_','',hallmark$term)
hallmark.list <- hallmark %>% split(.$term) %>% lapply( "[[", 2)

gsea.re1<-clusterProfiler::GSEA(id,TERM2GENE=hallmark,verbose = T)
g1<-as.data.frame(gsea.re1)
g1<-subset(g1,p.adjust<0.05)
g1<-g1[order(g1$NES,decreasing = T),]
type(g1)
head(g1)

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


library(gggsea)

library(ggsci)

names(hallmark.list)
se_hall<-c(head(g2$pathway,9),tail(g1$pathway,9))

#3 所选中的基因集
sig1<-g2[g2$pathway%in%se_hall,]

hallmark.se<-hallmark.list[sig1$pathway]
df.new <- gseaCurve(id, hallmark.se, sig1)

pal_line<-pal_lancet()(7)#各曲线颜色

p_nes<-ggplot() + 
  geom_gsea(df.new,
            prettyGSEA=T,
            tickcolor='grey30',#竖线颜色
            ticksize=0.1,#线条粗细
            #colour='grey',#色块边框
            linecolor=pal_line,
            linesize=2,lty=1,
            #多个图形时使用，设定每行列图形个数
            nrow = 2,
            ncol =4
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

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/上皮细胞/上皮细胞Capno/Capno_hallmark.pdf",p_nes,width=10,height=4)



######根据差异基因画多组火山图###----
library(ggplot2) # Create Elegant Data Visualisations Using the Grammar of Graphics
library(dplyr) # A Grammar of Data Manipulation
library(RColorBrewer) # ColorBrewer Palettes
library(grid) # The Grid Graphics Package
library(scales)

df <- read.table("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/Strep/免疫细胞分组/immune_cell_merged_data.csv", header = 1, check.names = F, sep = ",")
df$group <- factor(df$group, levels = c("Pre_NMPR","Pre_MPR","Post_MPR"))

df$group2<-as.factor(ifelse(df$p_val_adj < 0.05 & abs(df$avg_log2FC) >= 0.25, 
                            ifelse(df$avg_log2FC>=0.25 ,'Up','Down'),'NS'))
df$group2 <- factor(df$group2, levels = c("Up", "Down", "NS"))
##确定添加标签的数据
df$label<-ifelse(df$p_val_adj<0.05&abs(df$avg_log2FC)>=0.25,"Y","N")
df$label<-ifelse(df$label == 'Y', as.character(df$X), '')

write.csv(df,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/Strep/免疫细胞分组/immune_cell_merged_data注释.csv")

df<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/Strep/免疫细胞分组/immune_cell_merged_data注释.csv")

df$label2 <- ifelse(df$avg_log2FC > 20, df$label, "")

df <- df %>%
  group_by(group) %>%
  filter(label != "") %>% 
  arrange(desc(avg_log2FC)) %>%        # 按 logfc 降序排列
  mutate(label2 = ifelse(row_number() <= 5, label, "")) %>%
  ungroup()  

head(df)

gene_name<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/免疫细胞gene_name.csv")
intersection <- intersect(gene_name$Gene.symbol, df$X)

Imerged_data_filter_inter<- df %>% filter(df$X %in%intersection )

write.csv(Imerged_data_filter_inter,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/Strep/免疫细胞分组/merged_data_filter_label.csv")


df$group2<- factor(df$group2,levels=c("Up","NS"))

df_bg <- df %>%
  group_by(group) %>%
  summarize(max_log2FC = max(avg_log2FC),min_log2FC = min(avg_log2FC))

df_bg$group<- factor(df_bg$group, levels = c("Pre_NMPR","Pre_MPR","Post_MPR"))


library(ggrepel)

cell_type_color <- c("#FF6347","#ffc556","#1E90FF","#0000FF","#9370DB","#00FFFF","#32cd32","#006400","#FF00ff","#20B2AA")

treatment_color <- c("#4974a4","#4dae47","#f29600")

vol_plot_pre_mpr<-ggplot()+
  ##y轴正半轴的灰色背景
  geom_col(data = df_bg, 
           mapping = aes(group,max_log2FC),
           fill = "grey85", width = 0.8, alpha = 0.5) +
  ##y轴负半轴的灰色背景
  geom_col(data = df_bg, 
           mapping = aes(group, min_log2FC),
           fill = "grey85", width = 0.8, alpha = 0.5)+
  #添加各组的数据点
  geom_jitter(data = df,
              mapping = aes(x = group, y = avg_log2FC, color = group2),
              size= 2,width = 0.4, alpha = 0.7)+
  # 通过在X=0的位置添加方块进行展示分组信息，采用geom_col方法添加
  geom_col(data = df_bg,
           mapping = aes(x= group, y = 1, fill = group),
           width = 0.8)+
  geom_col(data = df_bg,
           mapping = aes(x= group, y = -1, fill = group),
           width = 0.8)+
  # 在方块中添加分组的文字信息
  geom_text(data=df_bg,
            mapping = aes(x=group, y=0, label=group),
            size = 5, color ="white",fontface = "bold")+
  #根据需要决定是添加辅助线
  # geom_hline(yintercept = 4, lty=2, color = '#ae63e4', lwd=0.8)+
  # geom_hline(yintercept = -4, lty=2, color = '#ae63e4', lwd=0.8)+
  #颜色设置
  scale_color_manual(values = c('brown','steelblue','gray'))+
  scale_fill_manual(values = treatment_color)+
  #添加显著性标签
  geom_text_repel(data = df,
                  mapping = aes(x = group, y = avg_log2FC, label = label2),
                  max.overlaps = 10000,
                  size=4,
                  box.padding=unit(0.3,'lines'),
                  point.padding=unit(0.3, 'lines'),
                  segment.color='black',
                  show.legend=FALSE)+
  # 主题设置
  theme_classic()+
  labs(x = "Srep+ vs Bacteria-", y = "avg_log2FC", fill= NULL, color = NULL)+
  theme(axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.y = element_line(linewidth = 0.8),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 14, color = "black"),
        axis.ticks.y = element_line(linewidth = 0.8),
        legend.text = element_text(size = 12),
        axis.title.x = element_text(size = 18),  # X轴标题字体
        axis.title.y = element_text(size = 18))+
  #调整图例
  guides(color=guide_legend(override.aes = list(size=7,alpha=1)))
vol_plot_pre_mpr

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/Strep/免疫细胞分组/免疫细胞火山图分组.pdf",vol_plot_pre_mpr,width = 8,height=7)

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/Strep/免疫细胞分组/免疫细胞火山图分组没有标签.pdf",vol_plot_pre_mpr,width = 8,height=5)

####组合免疫细胞和上皮细胞有无微生物功能差异----
####免疫细胞整体分析----
Idents(MPRNMPR_object_miMPRobe_remove) <- "cell_type_new"
MPRNMPR_object_miMPRobe_remove@meta.data$MPR_Response2
Immune_cell=subset(MPRNMPR_object_miMPRobe_remove,idents = c("T cells","B cells",
                                                                "Plasma cells","Myeloid cells","Mast cells"))

Immune_cell@meta.data$group_microbe2
Idents(Immune_cell) <- "group_microbe2"
Immune_cell_deg_all=FindMarkers(Immune_cell, ident.1 = c("Bacteria+"), ident.2 = "Bacteria-")
dim(Immune_cell_deg_all)
names(Immune_cell_deg_all)
write.csv(Immune_cell_deg_all,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/免疫细胞/Immune_cell_deg_all.csv")
Immune_cell_deg_all<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/免疫细胞/Immune_cell_deg_all.csv",header = T)

####根据过滤掉pt 小于0.1 的行

names(Immune_cell_deg_all)

Immune_cell_deg_all_filter<- Immune_cell_deg_all %>% filter(pct.1>0.1 &pct.2>0.1)

####读取免疫细胞相关基因，从差异基因集中筛选

gene_name<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/免疫细胞/免疫细胞gene_name.csv")

intersection <- intersect(gene_name$Gene.symbol, Immune_cell_deg_all_filter$X)

Immune_cell_deg_all_filter_inter<- Immune_cell_deg_all_filter %>% filter(Immune_cell_deg_all_filter$X %in%intersection )

df <- Immune_cell_deg_all_filter_inter

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
  ggtitle("Immune cell DEG (Bacteria+ vs Bacteria-)") + #标题
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
ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/免疫细胞/免疫细胞火山图.jpg",degp,width = 10,height=7)
ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/免疫细胞/免疫细胞火山图.pdf",degp,width = 10,height=7)


####免疫细胞分组分析----
# 获取所有csv文件列表
file_path = "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/免疫细胞/"
  cell_type_list <- c("Pre_NMPR","Pre_MPR","Post_MPR")
  dir.create(file_path)
  for (cell in cell_type_list){
    print(cell)
    Idents(Immune_cell)<-"MPR_Response2"
    sub_object=subset(Immune_cell,idents = cell)
    Idents(sub_object)<-"group_microbe2"
    deg_all=FindMarkers(sub_object, ident.1 = c("Bacteria+"), ident.2 = "Bacteria-")
    new_filename <- paste(cell, "_deg_all",".csv", sep="")
    print(new_filename)
    write.csv(deg_all,paste("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/免疫细胞/",new_filename,sep=""))
  }
  
  # 获取所有csv文件列表
  
  files <- list.files(path = file_path, pattern = "\\.csv$", recursive = TRUE)
  
  files <- list.files(path = file_path, pattern = "\\.csv$", recursive = TRUE)
  
  # 创建一个空数据框用于存储合并后的数据
  merged_data <- data.frame()
  # 循环读取每个csv文件并按照某一列降序排序，进行阈值筛选，并合并到merged_data中
  for (file in files) {
    data <- read.csv(paste(file_path,file,sep=""))
    sorted_data <- data[order(-data$p_val_adj), ]
    file<-basename(file)
    print(file)
    file_name <- rep(substr(file, 1, nchar(file)-12), nrow(sorted_data))
    #filtered_data <- sorted_data[sorted_data$pct.1 >= 0.25, ]
    sorted_data$group <- file_name
    merged_data <- rbind(merged_data, sorted_data)
  }
  
  # 导出合并后的数据为新的csv文件
  write.csv(merged_data, "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/免疫细胞/immune_cell_merged_data.csv", row.names = FALSE)
  
  
  ######根据差异基因画多组火山图###----
  library(ggplot2) # Create Elegant Data Visualisations Using the Grammar of Graphics
  library(dplyr) # A Grammar of Data Manipulation
  library(RColorBrewer) # ColorBrewer Palettes
  library(grid) # The Grid Graphics Package
  library(scales)
  
  df <- read.table("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/免疫细胞/immune_cell_merged_data.csv", header = 1, check.names = F, sep = ",")
  df$group <- factor(df$group, levels = c("Pre_NMPR","Pre_MPR","Post_MPR"))
  
  df$group2<-as.factor(ifelse(df$p_val_adj < 0.05 & abs(df$avg_log2FC) >= 0.25, 
                              ifelse(df$avg_log2FC>= 0.25 ,'Up','Down'),'NS'))
  df$group2 <- factor(df$group2, levels = c("Up", "Down", "NS"))
  ##确定添加标签的数据
  df$label<-ifelse(df$p_val_adj<0.05&abs(df$avg_log2FC)>=0.25,"Y","N")
  df$label<-ifelse(df$label == 'Y', as.character(df$X), '')
 
  gene_name<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/免疫细胞gene_name.csv")
  intersection <- intersect(gene_name$Gene.symbol, df$X)
  
  Imerged_data_filter_inter<- df %>% filter(df$X %in%intersection )
  
  write.csv(Imerged_data_filter_inter,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/免疫细胞/merged_data_filter_label.csv")

  
  df <- Imerged_data_filter_inter

  df_bg <- df %>%
    group_by(group) %>%
    summarize(max_log2FC = max(avg_log2FC),min_log2FC = min(avg_log2FC))
  
  df_bg$group<- factor(df_bg$group, levels = c("Pre_NMPR","Pre_MPR","Post_MPR"))
  
  
  library(ggrepel)
  
  cell_type_color <- c("#FF6347","#ffc556","#1E90FF","#0000FF","#9370DB","#00FFFF","#32cd32","#006400","#FF00ff","#20B2AA")
  
  treatment_color <- c("#4974a4","#4dae47","#f29600")
  
  vol_plot_pre_mpr<-ggplot()+
    ##y轴正半轴的灰色背景
    geom_col(data = df_bg, 
             mapping = aes(group,max_log2FC),
             fill = "grey85", width = 0.8, alpha = 0.5) +
    ##y轴负半轴的灰色背景
    geom_col(data = df_bg, 
             mapping = aes(group, min_log2FC),
             fill = "grey85", width = 0.8, alpha = 0.5)+
    #添加各组的数据点
    geom_jitter(data = df,
                mapping = aes(x = group, y = avg_log2FC, color = group2),
                size= 2,width = 0.4, alpha = 0.7)+
    # 通过在X=0的位置添加方块进行展示分组信息，采用geom_col方法添加
    geom_col(data = df_bg,
             mapping = aes(x= group, y = 0.2, fill = group),
             width = 0.8)+
    geom_col(data = df_bg,
             mapping = aes(x= group, y = -0.2, fill = group),
             width = 0.8)+
    # 在方块中添加分组的文字信息
    geom_text(data=df_bg,
              mapping = aes(x=group, y=0, label=group),
              size = 5, color ="white",fontface = "bold")+
    #根据需要决定是添加辅助线
    # geom_hline(yintercept = 4, lty=2, color = '#ae63e4', lwd=0.8)+
    # geom_hline(yintercept = -4, lty=2, color = '#ae63e4', lwd=0.8)+
    #颜色设置
    scale_color_manual(values = c('brown','steelblue','gray'))+
    scale_fill_manual(values = treatment_color)+
    #添加显著性标签
    geom_text_repel(data = df,
                    mapping = aes(x = group, y = avg_log2FC, label = label),
                    max.overlaps = 10000,
                    size=4,
                    box.padding=unit(0.3,'lines'),
                    point.padding=unit(0.3, 'lines'),
                    segment.color='black',
                    show.legend=FALSE)+
    # 主题设置
    theme_classic()+
    labs(x = "Bacteria+ vs Bacteria-", y = "avg_log2FC", fill= NULL, color = NULL)+
    theme(axis.line.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.line.y = element_line(linewidth = 0.8),
          axis.text.y = element_text(size = 12, color = "black"),
          axis.title = element_text(size = 14, color = "black"),
          axis.ticks.y = element_line(linewidth = 0.8),
          legend.text = element_text(size = 12),
          axis.title.x = element_text(size = 18),  # X轴标题字体
          axis.title.y = element_text(size = 18))+
    #调整图例
    guides(color=guide_legend(override.aes = list(size=7,alpha=1)))
  vol_plot_pre_mpr

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/免疫细胞/免疫细胞火山图分组.jpg",vol_plot_pre_mpr,width = 8,height=7)
ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/免疫细胞/免疫细胞火山图分组.pdf",vol_plot_pre_mpr,width = 8,height=7)

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/免疫细胞/免疫细胞火山图分组没有标签.jpg",vol_plot_pre_mpr,width = 8,height=5)
ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/免疫细胞/免疫细胞火山图分组没有标签.pdf",vol_plot_pre_mpr,width = 8,height=5)

####免疫细胞分组bacteria+ 和bacteria-----

######Pre_NMPR_deg_all----

Pre_NMPR_deg_all <- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/免疫细胞/Pre_NMPR_deg_all.csv")

deg<-Pre_NMPR_deg_all 

gene_name<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/免疫细胞/免疫细胞gene_name.csv")

intersection <- intersect(gene_name$Gene.symbol, deg$X)

deg_all_inter<-deg %>% filter(deg$X %in%intersection )

df <-deg_all_inter

df$group2<-as.factor(ifelse(df$p_val_adj < 0.05 & abs(df$avg_log2FC) >= 0.25, 
                            ifelse(df$avg_log2FC>= 0.25 ,'Up','Down'),'NS'))
df$group2 <- factor(df$group2, levels = c("Up", "Down", "NS"))
unique(df$group2)
head(df)
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
  #xlim(c(-1, 1)) +  ##调整x轴的取值范围，max(abs(BRCA_Match_DEG$logFC))最大值是多少，再进行取舍
  geom_vline(xintercept=c(-0.25,0.25),lty=4,col="black",lwd=0.8) + #添加x轴辅助线,lty函数调整线的类型
  geom_hline(yintercept = -log10(0.05), lty=4,col="black",lwd=0.8) +  #添加y轴辅助线
  labs(x="avg_log2FC", y="-log10(p_val_adj)") +  #x、y轴标签
  ggtitle("Pre_NMPR_immune cell DEG (Bacteria+ vs Bacteria-)") + #标题
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

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/免疫细胞/免疫细胞分组/Pre_NMPR_deg_火山图.pdf",degp,width = 7,height=7)

######Pre_NMPR----

library(org.Hs.eg.db)#物种注释包(Homo sapiens)
library(clusterProfiler)#富集分析
#添加entrez ID列：
##：symbol转entrez ID：
symbol <- deg$X
entrez <- bitr(symbol,
               fromType = "SYMBOL",#现有的ID类型
               toType = "ENTREZID",#需转换的ID类型
               OrgDb = "org.Hs.eg.db")
head(entrez)
genelist <- deg$avg_log2FC
names(genelist) <-deg$X
head(genelist)
#genelist过滤(ID转换中丢失的部分基因)：
genelist <- genelist[names(genelist) %in% entrez[,1]]
names(genelist) <- entrez[match(names(genelist),entrez[,1]),2]
length(genelist)#过滤后剩16161个基因
head(genelist)
#将genelist按照log2FC值从高到低进行排序：
genelist <- sort(genelist,decreasing = T)
head(genelist)

#GSEA_KEGG富集分析：
KEGG_ges <- gseKEGG(
  geneList = genelist,
  organism = "hsa",#不同物种选择可官方查询：https://www.genome.jp/kegg/catalog/org_list.html
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose = FALSE,
  eps = 0
)

#提取结果表：
KEGG_ges_result <- KEGG_ges@result

#显示前10列（第十一列为富集到该基因集的核心基因列表，这里过长不展示）：
KEGG_ges_result[1:3,1:10]
library(ggridges)
library(ggplot2)
library(enrichplot)
nrow(KEGG_ges_result)
p <- ridgeplot(KEGG_ges,
               showCategory = 15,
               fill = "p.adjust",
               decreasing  = T)
p

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/免疫细胞/免疫细胞分组/Pre_NMPR_ridgeplot_KEGG.pdf",p,width=7,height=6)


##GSEA_GO富集分析：
GO_ges <- gseGO(geneList = genelist,
                OrgDb = org.Hs.eg.db,
                ont = "CC", #one of "BP", "MF", and "CC" subontologies, or "ALL" for all three.
                minGSSize = 10,
                maxGSSize = 500,
                pvalueCutoff = 0.05,
                eps = 0,
                verbose = FALSE)
GO_ges[1:3,1:10] #展示同样省略最后一列
nrow(GO_ges)

p <- ridgeplot(GO_ges,
               showCategory = 15,
               fill = "p.adjust",
               decreasing  = T)
p

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/免疫细胞/免疫细胞分组/Pre_NMPR_ridgeplot_GO.pdf",p,width=7,height=9)


##GSEA_Reactome富集分析:
BiocManager::install("ReactomePA")
library(ReactomePA)
ret_ges <- gsePathway(genelist,
                      organism = "human",
                      minGSSize = 10,
                      maxGSSize = 500,
                      pvalueCutoff = 0.05,
                      pAdjustMethod = "BH",
                      verbose = FALSE,
                      eps = 0)
ret_ges_result<-ret_ges[1:34,1:10]
nrow(ret_ges)
p <- ridgeplot(ret_ges,
               showCategory = 34,
               fill = "p.adjust",
               decreasing  = T)
p
?ridgeplot

pp <- dotplot(ret_ges,
              showCategory = 34)
pp

####可视化柱状图
library(stringr)
library(ggplot2)
ret_ges_result <- ret_ges_result[order(ret_ges_result$NES),]
ret_ges_result$yax <- ifelse(ret_ges_result$NES >0, -0.02, 0.02)
ret_ges_result$col <- ifelse(ret_ges_result$NES > 0, "blue","red")
ret_ges_result$NES
names(ret_ges_result)

p <- 
  ggplot(ret_ges_result,aes(y=NES,x=reorder(Description,NES),label = Description))+
  geom_text(aes(y = yax,hjust = NES > 0),size=4)+ 
  geom_bar(stat="identity",aes(fill=col))+
  geom_hline(yintercept = 0)+
  coord_flip(ylim = c(-2,2.5))+
  labs(y = "Normalised enrichment score", x = "",title="")+
  scale_fill_manual(values = c("#FB6542","#375E97"),name="Expression",labels=c("Up","Down"))+
  scale_x_discrete(breaks = NULL)+
  theme_classic(base_size =  12)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.y = element_blank(),
        legend.position = 'bottom')
p

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/免疫细胞/免疫细胞分组/Pre_NMPR_barplot_ReactomePA.pdf",p,width=7,height=9)


##GSEA_DO(Disease Ontology)富集分析:
library(DOSE)
DO_ges <- gseDO(genelist,
                minGSSize = 10,
                maxGSSize = 500,
                pvalueCutoff = 0.05,
                pAdjustMethod = "BH",
                verbose = FALSE,
                eps = 0)
DO_ges[1:3,1:10]
nrow(DO_ges)

p <- ridgeplot(DO_ges,
               showCategory = 10,
               fill = "p.adjust",
               decreasing  = T)
p


#####GSEA_MSigDB富集分析：
library(gggsea)
library(ggplot2)
library(GSVA)
library(enrichplot)
library(GSEABase)
library(fgsea)
library(clusterProfiler)
gmtfile <- "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/h.all.v2024.1.Hs.symbols.gmt"
hallmark <- read.gmt(gmtfile)
alldiff <- deg[order(deg$avg_log2FC,decreasing = T),]
id <- alldiff$avg_log2FC
names(id) <- alldiff$X

## 这里去掉了基因集前缀
hallmark$term <- gsub('HALLMARK_','',hallmark$term)
hallmark.list <- hallmark %>% split(.$term) %>% lapply( "[[", 2)

gsea.re1<-clusterProfiler::GSEA(id,TERM2GENE=hallmark,verbose = T)
g1<-as.data.frame(gsea.re1)
g1<-subset(g1,p.adjust<0.05)
g1<-g1[order(g1$NES,decreasing = T),]
type(g1)
head(g1)

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


library(gggsea)

library(ggsci)

names(hallmark.list)
se_hall<-c(head(g2$pathway,9),tail(g1$pathway,9))

#3 所选中的基因集
sig1<-g2[g2$pathway%in%se_hall,]

hallmark.se<-hallmark.list[sig1$pathway]
df.new <- gseaCurve(id, hallmark.se, sig1)

pal_line<-pal_lancet()(7)#各曲线颜色

p_nes<-ggplot() + 
  geom_gsea(df.new,
            prettyGSEA=T,
            tickcolor='grey30',#竖线颜色
            ticksize=0.1,#线条粗细
            #colour='grey',#色块边框
            linecolor=pal_line,
            linesize=2,lty=1,
            #多个图形时使用，设定每行列图形个数
            nrow = 2,
            ncol =4
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

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/免疫细胞/免疫细胞分组/Pre_NMPR_hallmark.pdf",p_nes,width=4,height=4)


####Pre_MPR----

Pre_MPR_deg_all <- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/免疫细胞/Pre_MPR_deg_all.csv")

deg<-Pre_MPR_deg_all 

gene_name<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/免疫细胞/免疫细胞gene_name.csv")

intersection <- intersect(gene_name$Gene.symbol, deg$X)

deg_all_inter<-deg %>% filter(deg$X %in%intersection )

df <-deg_all_inter

df$group2<-as.factor(ifelse(df$p_val_adj < 0.05 & abs(df$avg_log2FC) >= 0.25, 
                            ifelse(df$avg_log2FC>= 0.25 ,'Up','Down'),'NS'))
df$group2 <- factor(df$group2, levels = c("Up", "Down", "NS"))
unique(df$group2)
head(df)
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
  #xlim(c(-1, 1)) +  ##调整x轴的取值范围，max(abs(BRCA_Match_DEG$logFC))最大值是多少，再进行取舍
  geom_vline(xintercept=c(-0.25,0.25),lty=4,col="black",lwd=0.8) + #添加x轴辅助线,lty函数调整线的类型
  geom_hline(yintercept = -log10(0.05), lty=4,col="black",lwd=0.8) +  #添加y轴辅助线
  labs(x="avg_log2FC", y="-log10(p_val_adj)") +  #x、y轴标签
  ggtitle("Pre_NMPR_immune cell DEG (Bacteria+ vs Bacteria-)") + #标题
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

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/免疫细胞/免疫细胞分组/Pre_MPR_deg_火山图.pdf",degp,width = 7,height=7)


library(org.Hs.eg.db)#物种注释包(Homo sapiens)
library(clusterProfiler)#富集分析
#添加entrez ID列：
##：symbol转entrez ID：
symbol <- deg$X
entrez <- bitr(symbol,
               fromType = "SYMBOL",#现有的ID类型
               toType = "ENTREZID",#需转换的ID类型
               OrgDb = "org.Hs.eg.db")
head(entrez)
genelist <- deg$avg_log2FC
names(genelist) <-deg$X
head(genelist)
#genelist过滤(ID转换中丢失的部分基因)：
genelist <- genelist[names(genelist) %in% entrez[,1]]
names(genelist) <- entrez[match(names(genelist),entrez[,1]),2]
length(genelist)#过滤后剩16161个基因
head(genelist)
#将genelist按照log2FC值从高到低进行排序：
genelist <- sort(genelist,decreasing = T)
head(genelist)

#GSEA_KEGG富集分析：
KEGG_ges <- gseKEGG(
  geneList = genelist,
  organism = "hsa",#不同物种选择可官方查询：https://www.genome.jp/kegg/catalog/org_list.html
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose = FALSE,
  eps = 0
)

#提取结果表：
KEGG_ges_result <- KEGG_ges@result

#显示前10列（第十一列为富集到该基因集的核心基因列表，这里过长不展示）：
KEGG_ges_result[1:3,1:10]
library(ggridges)
library(ggplot2)
library(enrichplot)
nrow(KEGG_ges_result)
p <- ridgeplot(KEGG_ges,
               showCategory = 25,
               fill = "p.adjust",
               decreasing  = T)
p

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/免疫细胞/免疫细胞分组/Pre_MPR_ridgeplot_KEGG.pdf",p,width=7,height=6)


##GSEA_GO富集分析：
GO_ges <- gseGO(geneList = genelist,
                OrgDb = org.Hs.eg.db,
                ont = "CC", #one of "BP", "MF", and "CC" subontologies, or "ALL" for all three.
                minGSSize = 10,
                maxGSSize = 500,
                pvalueCutoff = 0.05,
                eps = 0,
                verbose = FALSE)
GO_ges[1:49,1:10] #展示同样省略最后一列
nrow(GO_ges)
GO_ges_result<-GO_ges@result
p <- ridgeplot(GO_ges,
               showCategory = 15,
               fill = "p.adjust",
               decreasing  = T)
p



####可视化柱状图
library(stringr)
library(ggplot2)
GO_ges_result<-as.data.frame(GO_ges_result)
GO_ges_result<-subset(GO_ges_result,p.adjust<0.001)
GO_ges_result<-GO_ges_result[order(GO_ges_result$NES,decreasing = T),]
GO_ges_result$yax <- ifelse(GO_ges_result$NES >0, -0.02, 0.02)
GO_ges_result$col <- ifelse(GO_ges_result$NES > 0, "blue","red")
GO_ges_result$NES
names(GO_ges_result)

p <- 
  ggplot(GO_ges_result,aes(y=NES,x=reorder(Description,NES),label = Description))+
  geom_text(aes(y = yax,hjust = NES > 0),size=4)+ 
  geom_bar(stat="identity",aes(fill=col))+
  geom_hline(yintercept = 0)+
  coord_flip(ylim = c(-2,2.5))+
  labs(y = "Normalised enrichment score", x = "",title="")+
  scale_fill_manual(values = c("#FB6542","#375E97"),name="Expression",labels=c("Up","Down"))+
  scale_x_discrete(breaks = NULL)+
  theme_classic(base_size =  12)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.y = element_blank(),
        legend.position = 'bottom')
p


ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/免疫细胞/免疫细胞分组/Pre_MPR_bar_GO_0.001.pdf",p,width=7,height=7)


##GSEA_Reactome富集分析:
BiocManager::install("ReactomePA")
library(ReactomePA)
ret_ges <- gsePathway(genelist,
                      organism = "human",
                      minGSSize = 10,
                      maxGSSize = 500,
                      pvalueCutoff = 0.05,
                      pAdjustMethod = "BH",
                      verbose = FALSE,
                      eps = 0)
ret_ges_result<-ret_ges[1:34,1:10]
nrow(ret_ges)
p <- ridgeplot(ret_ges,
               showCategory = 34,
               fill = "p.adjust",
               decreasing  = T)
p
?ridgeplot

pp <- dotplot(ret_ges,
              showCategory = 34)
pp

####可视化柱状图
library(stringr)
library(ggplot2)
ret_ges_result <- ret_ges_result[order(ret_ges_result$NES),]
ret_ges_result$yax <- ifelse(ret_ges_result$NES >0, -0.02, 0.02)
ret_ges_result$col <- ifelse(ret_ges_result$NES > 0, "blue","red")
ret_ges_result$NES
names(ret_ges_result)

p <- 
  ggplot(ret_ges_result,aes(y=NES,x=reorder(Description,NES),label = Description))+
  geom_text(aes(y = yax,hjust = NES > 0),size=4)+ 
  geom_bar(stat="identity",aes(fill=col))+
  geom_hline(yintercept = 0)+
  coord_flip(ylim = c(-2,2.5))+
  labs(y = "Normalised enrichment score", x = "",title="")+
  scale_fill_manual(values = c("#FB6542","#375E97"),name="Expression",labels=c("Up","Down"))+
  scale_x_discrete(breaks = NULL)+
  theme_classic(base_size =  12)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.y = element_blank(),
        legend.position = 'bottom')
p

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/免疫细胞/免疫细胞分组/Pre_MPR_barplot_ReactomePA.pdf",p,width=7,height=9)


##GSEA_DO(Disease Ontology)富集分析:
library(DOSE)
DO_ges <- gseDO(genelist,
                minGSSize = 10,
                maxGSSize = 500,
                pvalueCutoff = 0.05,
                pAdjustMethod = "BH",
                verbose = FALSE,
                eps = 0)
DO_ges[1:3,1:10]
nrow(DO_ges)

p <- ridgeplot(DO_ges,
               showCategory = 10,
               fill = "p.adjust",
               decreasing  = T)
p


#####GSEA_MSigDB富集分析：
library(gggsea)
library(ggplot2)
library(GSVA)
library(enrichplot)
library(GSEABase)
library(fgsea)
library(clusterProfiler)
gmtfile <- "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/h.all.v2024.1.Hs.symbols.gmt"
hallmark <- read.gmt(gmtfile)
alldiff <- deg[order(deg$avg_log2FC,decreasing = T),]
id <- alldiff$avg_log2FC
names(id) <- alldiff$X

## 这里去掉了基因集前缀
hallmark$term <- gsub('HALLMARK_','',hallmark$term)
hallmark.list <- hallmark %>% split(.$term) %>% lapply( "[[", 2)

gsea.re1<-clusterProfiler::GSEA(id,TERM2GENE=hallmark,verbose = T)
g1<-as.data.frame(gsea.re1)
g1<-subset(g1,p.adjust<0.05)
g1<-g1[order(g1$NES,decreasing = T),]
type(g1)
head(g1)

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

library(gggsea)
library(ggsci)

names(hallmark.list)
se_hall<-c(head(g2$pathway,11),tail(g1$pathway,9))

#3 所选中的基因集
sig1<-g2[g2$pathway%in%se_hall,]

hallmark.se<-hallmark.list[sig1$pathway]
df.new <- gseaCurve(id, hallmark.se, sig1)

pal_line<-pal_lancet()(11)#各曲线颜色

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
            ncol =3
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

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/免疫细胞/免疫细胞分组/Pre_MPR_hallmark.pdf",p_nes,width=6,height=8)



####Post_MPR----

Post_MPR_deg_all <- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/免疫细胞/Post_MPR_deg_all.csv")

deg<-Post_MPR_deg_all 

gene_name<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/免疫细胞/免疫细胞gene_name.csv")

intersection <- intersect(gene_name$Gene.symbol, deg$X)

deg_all_inter<-deg %>% filter(deg$X %in%intersection )

df <-deg_all_inter

df$group2<-as.factor(ifelse(df$p_val_adj < 0.05 & abs(df$avg_log2FC) >= 0.25, 
                            ifelse(df$avg_log2FC>= 0.25 ,'Up','Down'),'NS'))
df$group2 <- factor(df$group2, levels = c("Up", "Down", "NS"))
unique(df$group2)
head(df)
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
  #xlim(c(-1, 1)) +  ##调整x轴的取值范围，max(abs(BRCA_Match_DEG$logFC))最大值是多少，再进行取舍
  geom_vline(xintercept=c(-0.25,0.25),lty=4,col="black",lwd=0.8) + #添加x轴辅助线,lty函数调整线的类型
  geom_hline(yintercept = -log10(0.05), lty=4,col="black",lwd=0.8) +  #添加y轴辅助线
  labs(x="avg_log2FC", y="-log10(p_val_adj)") +  #x、y轴标签
  ggtitle("Pre_NMPR_immune cell DEG (Bacteria+ vs Bacteria-)") + #标题
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

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/免疫细胞/免疫细胞分组/Pre_MPR_deg_火山图.pdf",degp,width = 7,height=7)


library(org.Hs.eg.db)#物种注释包(Homo sapiens)
library(clusterProfiler)#富集分析
#添加entrez ID列：
##：symbol转entrez ID：
symbol <- deg$X
entrez <- bitr(symbol,
               fromType = "SYMBOL",#现有的ID类型
               toType = "ENTREZID",#需转换的ID类型
               OrgDb = "org.Hs.eg.db")
head(entrez)
genelist <- deg$avg_log2FC
names(genelist) <-deg$X
head(genelist)
#genelist过滤(ID转换中丢失的部分基因)：
genelist <- genelist[names(genelist) %in% entrez[,1]]
names(genelist) <- entrez[match(names(genelist),entrez[,1]),2]
length(genelist)#过滤后剩16161个基因
head(genelist)
#将genelist按照log2FC值从高到低进行排序：
genelist <- sort(genelist,decreasing = T)
head(genelist)

#GSEA_KEGG富集分析：
KEGG_ges <- gseKEGG(
  geneList = genelist,
  organism = "hsa",#不同物种选择可官方查询：https://www.genome.jp/kegg/catalog/org_list.html
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose = FALSE,
  eps = 0
)

#提取结果表：
KEGG_ges_result <- KEGG_ges@result

#显示前10列（第十一列为富集到该基因集的核心基因列表，这里过长不展示）：
KEGG_ges_result[1:3,1:10]
library(ggridges)
library(ggplot2)
library(enrichplot)
nrow(KEGG_ges_result)
p <- ridgeplot(KEGG_ges,
               showCategory = 25,
               fill = "p.adjust",
               decreasing  = T)
p

library(stringr)
library(ggplot2)
KEGG_ges_result <- KEGG_ges_result[order(KEGG_ges_result$NES),]
KEGG_ges_result$yax <- ifelse(KEGG_ges_result$NES >0, -0.02, 0.02)
KEGG_ges_result$col <- ifelse(KEGG_ges_result$NES > 0, "blue","red")
KEGG_ges_result$NES
names(KEGG_ges_result)

p <- 
  ggplot(KEGG_ges_result,aes(y=NES,x=reorder(Description,NES),label = Description))+
  geom_text(aes(y = yax,hjust = NES > 0),size=4)+ 
  geom_bar(stat="identity",aes(fill=col))+
  geom_hline(yintercept = 0)+
  coord_flip(ylim = c(-2,2.5))+
  labs(y = "Normalised enrichment score", x = "",title="")+
  scale_fill_manual(values = c("#FB6542","#375E97"),name="Expression",labels=c("Up","Down"))+
  scale_x_discrete(breaks = NULL)+
  theme_classic(base_size =  12)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.y = element_blank(),
        legend.position = 'bottom')
p


ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/免疫细胞/免疫细胞分组/Post_MPR_barplot_KEGG.pdf",p,width=7,height=8)


##GSEA_GO富集分析：
GO_ges <- gseGO(geneList = genelist,
                OrgDb = org.Hs.eg.db,
                ont = "CC", #one of "BP", "MF", and "CC" subontologies, or "ALL" for all three.
                minGSSize = 10,
                maxGSSize = 500,
                pvalueCutoff = 0.05,
                eps = 0,
                verbose = FALSE)
GO_ges[1:49,1:10] #展示同样省略最后一列
nrow(GO_ges)
GO_ges_result<-GO_ges@result
p <- ridgeplot(GO_ges,
               showCategory = 15,
               fill = "p.adjust",
               decreasing  = T)
p

####可视化柱状图
library(stringr)
library(ggplot2)

GO_ges_result<-subset(GO_ges_result,p.adjust<0.001)
GO_ges_result<-GO_ges_result[order(GO_ges_result$NES,decreasing = T),]


GO_ges_result <- GO_ges_result[order(GO_ges_result$NES),]
GO_ges_result$yax <- ifelse(GO_ges_result$NES >0, -0.02, 0.02)
GO_ges_result$col <- ifelse(GO_ges_result$NES > 0, "blue","red")
GO_ges_result$NES
names(GO_ges_result)

p <- 
  ggplot(GO_ges_result,aes(y=NES,x=reorder(Description,NES),label = Description))+
  geom_text(aes(y = yax,hjust = NES > 0),size=4)+ 
  geom_bar(stat="identity",aes(fill=col))+
  geom_hline(yintercept = 0)+
  coord_flip(ylim = c(-3,3))+
  labs(y = "Normalised enrichment score", x = "",title="")+
  scale_fill_manual(values = c("#FB6542","#375E97"),name="Expression",labels=c("Up","Down"))+
  scale_x_discrete(breaks = NULL)+
  theme_classic(base_size =  12)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.y = element_blank(),
        legend.position = 'bottom')
p

p <- 
  ggplot(GO_ges_result,aes(y=NES,x=reorder(Description,NES),label = Description))+
  geom_text(aes(y = yax,hjust = NES > 0),size=4,angle=90)+ 
  geom_bar(stat="identity",aes(fill=col))+
  geom_hline(yintercept = 0)+
  #coord_flip(xlim = c(-2,2.5))+
  labs(y = "Normalised enrichment score", x = "",title="")+
  scale_fill_manual(values = c("#FB6542","#375E97"),name="Expression",labels=c("Up","Down"))+
  scale_x_discrete(breaks = NULL)+
  theme_classic(base_size =  12)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.x = element_blank(),
        legend.position = 'bottom')
p

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/免疫细胞/免疫细胞分组/Post_MPR_bar_GO_0.001.pdf",p,width=7,height=7)


##GSEA_Reactome富集分析:
BiocManager::install("ReactomePA")
library(ReactomePA)
ret_ges <- gsePathway(genelist,
                      organism = "human",
                      minGSSize = 10,
                      maxGSSize = 500,
                      pvalueCutoff = 0.05,
                      pAdjustMethod = "BH",
                      verbose = FALSE,
                      eps = 0)
ret_ges_result<-ret_ges[1:34,1:10]
nrow(ret_ges)
p <- ridgeplot(ret_ges,
               showCategory = 34,
               fill = "p.adjust",
               decreasing  = T)
p
?ridgeplot

pp <- dotplot(ret_ges,
              showCategory = 34)
pp

####可视化柱状图
library(stringr)
library(ggplot2)
ret_ges_result <- ret_ges_result[order(ret_ges_result$NES),]
ret_ges_result$yax <- ifelse(ret_ges_result$NES >0, -0.02, 0.02)
ret_ges_result$col <- ifelse(ret_ges_result$NES > 0, "blue","red")
ret_ges_result$NES
names(ret_ges_result)

p <- 
  ggplot(ret_ges_result,aes(y=NES,x=reorder(Description,NES),label = Description))+
  geom_text(aes(y = yax,hjust = NES > 0),size=4)+ 
  geom_bar(stat="identity",aes(fill=col))+
  geom_hline(yintercept = 0)+
  coord_flip(ylim = c(-2,2.5))+
  labs(y = "Normalised enrichment score", x = "",title="")+
  scale_fill_manual(values = c("#FB6542","#375E97"),name="Expression",labels=c("Up","Down"))+
  scale_x_discrete(breaks = NULL)+
  theme_classic(base_size =  12)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.y = element_blank(),
        legend.position = 'bottom')
p

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/免疫细胞/免疫细胞分组/Pre_MPR_barplot_ReactomePA.pdf",p,width=7,height=9)


##GSEA_DO(Disease Ontology)富集分析:
library(DOSE)
DO_ges <- gseDO(genelist,
                minGSSize = 10,
                maxGSSize = 500,
                pvalueCutoff = 0.05,
                pAdjustMethod = "BH",
                verbose = FALSE,
                eps = 0)
DO_ges[1:3,1:10]
nrow(DO_ges)

p <- ridgeplot(DO_ges,
               showCategory = 10,
               fill = "p.adjust",
               decreasing  = T)
p


#####GSEA_MSigDB富集分析：
library(gggsea)
library(ggplot2)
library(GSVA)
library(enrichplot)
library(GSEABase)
library(fgsea)
library(clusterProfiler)
gmtfile <- "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/h.all.v2024.1.Hs.symbols.gmt"
hallmark <- read.gmt(gmtfile)
alldiff <- deg[order(deg$avg_log2FC,decreasing = T),]
id <- alldiff$avg_log2FC
names(id) <- alldiff$X

## 这里去掉了基因集前缀
hallmark$term <- gsub('HALLMARK_','',hallmark$term)
hallmark.list <- hallmark %>% split(.$term) %>% lapply( "[[", 2)

gsea.re1<-clusterProfiler::GSEA(id,TERM2GENE=hallmark,verbose = T)
g1<-as.data.frame(gsea.re1)
g1<-subset(g1,p.adjust<0.05)
g1<-g1[order(g1$NES,decreasing = T),]
type(g1)
head(g1)

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

library(gggsea)
library(ggsci)

names(hallmark.list)
se_hall<-c(head(g2$pathway,15),tail(g1$pathway,9))

#3 所选中的基因集
sig1<-g2[g2$pathway%in%se_hall,]

hallmark.se<-hallmark.list[sig1$pathway]
df.new <- gseaCurve(id, hallmark.se, sig1)

pal_line<-pal_lancet()(11)#各曲线颜色

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
            ncol =3
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


library(stringr)
library(ggplot2)

gsea.re2 <- fgseaMultilevel(pathways = hallmark.list,
                            stats = id,
                            minSize=1,
                            maxSize=10000)

gsea.re2 <- gsea.re2[gsea.re2$padj< 0.05,]
gsea.re2 <- gsea.re2[order(gsea.re2$NES),]
gsea.re2$yax <- ifelse(gsea.re2$NES >0, -0.02, 0.02)
gsea.re2$col <- ifelse(gsea.re2$NES > 0, "blue","red")
gsea.re2$NES
names(gsea.re2)

p <- ggplot(gsea.re2,aes(y=NES,x=reorder(pathway,NES),label = pathway))+
  geom_text(aes(y = yax,hjust = NES > 0),size=4)+ 
  geom_bar(stat="identity",aes(fill=col))+
  geom_hline(yintercept = 0)+
  coord_flip(ylim = c(-2,2.5))+
  labs(y = "Normalised enrichment score", x = "",title="")+
  scale_fill_manual(values = c("#FB6542","#375E97"),name="Expression",labels=c("Up","Down"))+
  scale_x_discrete(breaks = NULL)+
  theme_classic(base_size =  12)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.y = element_blank(),
        legend.position = 'bottom')
p


ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/免疫细胞/免疫细胞分组/Post_MPR_hallmark.pdf",p_nes,width=10,height=10)

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/免疫细胞/免疫细胞分组/Post_MPR_hallmark_barplot0.05.pdf",p,width=7,height=7)


####整个上皮细胞有无微生物功能影响----
####获取恶性和非恶性细胞之间的差异基因

Idents(MPRNMPR_object_miMPRobe_remove) <- "cell_type_new"
unique(Idents(MPRNMPR_object_miMPRobe_remove))
MPRNMPR_object_miMPRobe_remove@meta.data$MPR_Response2
Epithelial_cell=subset(MPRNMPR_object_miMPRobe_remove,idents = c("Epithelial cells"))

Epithelial_cell@meta.data$epi_infercnv_group

Idents(Epithelial_cell) <- "epi_infercnv_group"

Malignant_Epithelial_cell_deg_all=FindMarkers(Epithelial_cell, ident.1 = c("Malignant epithelial cells"), ident.2 = "Non-malignant epithelial cells")
Malignant_Epithelial_cell_deg_all_filter<- Malignant_Epithelial_cell_deg_all %>% filter(pct.1>0.1 &pct.2>0.1)
Malignant_Epithelial_cell_deg_all_filter$group2<-as.factor(ifelse(Malignant_Epithelial_cell_deg_all_filter$p_val_adj < 0.05 & abs(Malignant_Epithelial_cell_deg_all_filter$avg_log2FC) >= 0.25, 
                            ifelse(Malignant_Epithelial_cell_deg_all_filter$avg_log2FC>= 0.25 ,'Up','Down'),'NS'))
Malignant_Epithelial_cell_deg_all_filter$group2 <- factor(Malignant_Epithelial_cell_deg_all_filter$group2, levels = c("Up", "Down", "NS"))

write.csv(Malignant_Epithelial_cell_deg_all_filter,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/上皮细胞/恶性和非恶性差异基因.csv")

#+++++++++++++++++
Epithelial_cell@meta.data$group_microbe2
Idents(Epithelial_cell) <- "group_microbe2"
Epithelial_cell_deg_all=FindMarkers(Epithelial_cell, ident.1 = c("Bacteria+"), ident.2 = "Bacteria-")
dim(Epithelial_cell_deg_all)
names(Epithelial_cell_deg_all)
write.csv(Epithelial_cell_deg_all,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/上皮细胞/Epithelial_cell_deg_all.csv")
Epithelial_cell_deg_all<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/上皮细胞/Epithelial_cell_deg_all.csv",header = T)

####根据过滤掉pt 小于0.1 的行

names(Epithelial_cell_deg_all)

Epithelial_cell_deg_all_filter<- Epithelial_cell_deg_all %>% filter(pct.1>0.1 &pct.2>0.1)

####读取免疫细胞相关基因，从差异基因集中筛选

gene_name<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/免疫细胞/免疫细胞gene_name.csv")

gene_name<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/上皮细胞/上皮细胞gene_name.csv")

intersection <- intersect(gene_name$Gene.symbol, Epithelial_cell_deg_all_filter$X)

Epithelial_cell_deg_all_filter_inter<- Epithelial_cell_deg_all_filter %>% filter(Epithelial_cell_deg_all_filter$X %in%intersection )

df <- Epithelial_cell_deg_all_filter_inter

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
  ggtitle("Epithelial_cell DEG (Bacteria+ vs Bacteria-)") + #标题
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
ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/上皮细胞/上皮细胞火山图.jpg",degp,width = 10,height=7)
ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/上皮细胞/上皮细胞火山图.pdf",degp,width = 10,height=7)

####上皮细胞恶性非恶性分组----
Epithelial_cell@meta.data$epi_infercnv_group
file_path = "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/上皮细胞/恶性非恶性差异基因/"
cell_type_list <- c("Non-malignant epithelial cells","Malignant epithelial cells")
dir.create(file_path)
for (cell in cell_type_list){
  print(cell)
  Idents(Epithelial_cell)<-"epi_infercnv_group"
  sub_object=subset(Epithelial_cell,idents = cell)
  Idents(sub_object)<-"group_microbe2"
  deg_all=FindMarkers(sub_object, ident.1 = c("Bacteria+"), ident.2 = "Bacteria-")
  new_filename <- paste(cell, "_deg_all",".csv", sep="")
  print(new_filename)
  write.csv(deg_all,paste("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/上皮细胞/恶性非恶性差异基因",new_filename,sep=""))
}

# 获取所有csv文件列表

files <- list.files(path = file_path, pattern = "\\.csv$", recursive = TRUE)

files <- list.files(path = file_path, pattern = "\\.csv$", recursive = TRUE)

# 创建一个空数据框用于存储合并后的数据
merged_data <- data.frame()
# 循环读取每个csv文件并按照某一列降序排序，进行阈值筛选，并合并到merged_data中
for (file in files) {
  data <- read.csv(paste(file_path,file,sep=""))
  sorted_data <- data[order(-data$p_val_adj), ]
  file<-basename(file)
  print(file)
  file_name <- rep(substr(file, 1, nchar(file)-12), nrow(sorted_data))
  #filtered_data <- sorted_data[sorted_data$pct.1 >= 0.25, ]
  sorted_data$group <- file_name
  merged_data <- rbind(merged_data, sorted_data)
}

# 导出合并后的数据为新的csv文件
write.csv(merged_data, "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/上皮细胞/恶性非恶性差异基因/Epithelial_cell_merged_data.csv", row.names = FALSE)


library(ggplot2) # Create Elegant Data Visualisations Using the Grammar of Graphics
library(dplyr) # A Grammar of Data Manipulation
library(RColorBrewer) # ColorBrewer Palettes
library(grid) # The Grid Graphics Package
library(scales)

df <- read.table("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/上皮细胞/恶性非恶性差异基因/Epithelial_cell_merged_data.csv", header = 1, check.names = F, sep = ",")

df$group <- factor(df$group, levels = c("Non-malignant epi","Malignant epi"))

df$group2<-as.factor(ifelse(df$p_val_adj < 0.05 & abs(df$avg_log2FC) >= 0.25, 
                            ifelse(df$avg_log2FC>= 0.25 ,'Up','Down'),'NS'))
df$group2 <- factor(df$group2, levels = c("Up", "Down", "NS"))
####确定添加标签的数据
df$label<-ifelse(df$p_val_adj<0.05&abs(df$avg_log2FC)>=0.25,"Y","N")
df$label<-ifelse(df$label == 'Y', as.character(df$X), '')

gene_name<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/上皮细胞/上皮细胞gene_name.csv")
intersection <- intersect(gene_name$Gene.symbol, df$X)

Imerged_data_filter_inter<- df %>% filter(df$X %in%intersection )

write.csv(Imerged_data_filter_inter,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/上皮细胞/恶性非恶性差异基因/merged_data_filter_label.csv")

df <- Imerged_data_filter_inter

df_bg <- df %>%
  group_by(group) %>%
  summarize(max_log2FC = max(avg_log2FC),min_log2FC = min(avg_log2FC))

df_bg$group<- factor(df_bg$group, levels = c("Non-malignant epi","Malignant epi"))


library(ggrepel)

cell_type_color <- c("#FF6347","#ffc556","#1E90FF","#0000FF","#9370DB","#00FFFF","#32cd32","#006400","#FF00ff","#20B2AA")

malignant_type_color<-c("#1fa84a","#9c67a4")

vol_plot_pre_mpr<-ggplot()+
  ##y轴正半轴的灰色背景
  geom_col(data = df_bg, 
           mapping = aes(group,max_log2FC),
           fill = "grey85", width = 0.8, alpha = 0.5) +
  ##y轴负半轴的灰色背景
  geom_col(data = df_bg, 
           mapping = aes(group, min_log2FC),
           fill = "grey85", width = 0.8, alpha = 0.5)+
  #添加各组的数据点
  geom_jitter(data = df,
              mapping = aes(x = group, y = avg_log2FC, color = group2),
              size= 2,width = 0.4, alpha = 0.7)+
  # 通过在X=0的位置添加方块进行展示分组信息，采用geom_col方法添加
  geom_col(data = df_bg,
           mapping = aes(x= group, y = 0.2, fill = group),
           width = 0.8)+
  geom_col(data = df_bg,
           mapping = aes(x= group, y = -0.2, fill = group),
           width = 0.8)+
  # 在方块中添加分组的文字信息
  geom_text(data=df_bg,
            mapping = aes(x=group, y=0, label=group),
            size = 5, color ="white",fontface = "bold")+
  #根据需要决定是添加辅助线
  # geom_hline(yintercept = 4, lty=2, color = '#ae63e4', lwd=0.8)+
  # geom_hline(yintercept = -4, lty=2, color = '#ae63e4', lwd=0.8)+
  #颜色设置
  scale_color_manual(values = c('brown','steelblue','gray'))+
  scale_fill_manual(values = malignant_type_color)+
  #添加显著性标签
  geom_text_repel(data = df,
                  mapping = aes(x = group, y = avg_log2FC, label = label),
                  max.overlaps = 10000,
                  size=4,
                  box.padding=unit(0.3,'lines'),
                  point.padding=unit(0.3, 'lines'),
                  segment.color='black',
                  show.legend=FALSE)+
  # 主题设置
  theme_classic()+
  labs(x = "Bacteria+ vs Bacteria-", y = "avg_log2FC", fill= NULL, color = NULL)+
  theme(axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.y = element_line(linewidth = 0.8),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 14, color = "black"),
        axis.ticks.y = element_line(linewidth = 0.8),
  legend.text = element_text(size = 12),
  axis.title.x = element_text(size = 18),  # X轴标题字体
  axis.title.y = element_text(size = 18))+
  #调整图例
  guides(color=guide_legend(override.aes = list(size=7,alpha=1)))
vol_plot_pre_mpr

####恶性基因差异基因的差异基因
m_up<- read.csv("clipboard")
nm_up<-read.csv("clipboard")

up_差异1<-setdiff(m_up[,1],nm_up[,1])
up_差异2<-setdiff(nm_up[,1],m_up[,1])

m_down<- read.csv("clipboard")
nm_down<-read.csv("clipboard")

down_差异1<-setdiff(m_down[,1],nm_down[,1])

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/上皮细胞/恶性非恶性差异基因/上皮细胞火山图分组.jpg",vol_plot_pre_mpr,width = 8,height=7)
ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/上皮细胞/恶性非恶性差异基因/上皮细胞火山图分组.pdf",vol_plot_pre_mpr,width = 8,height=7)


####上皮细胞恶性非恶性GSEA----


####恶性细胞----

mepi_deg_all <- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/上皮细胞/恶性非恶性差异基因/Malignant epithelial cells_deg_all.csv")

deg<-mepi_deg_all

gene_name<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/上皮细胞/上皮细胞gene_name.csv")

intersection <- intersect(gene_name$Gene.symbol, deg$X)

deg_all_inter<-deg %>% filter(deg$X %in%intersection )

df <-deg_all_inter

df$group2<-as.factor(ifelse(df$p_val_adj < 0.05 & abs(df$avg_log2FC) >= 0.25, 
                            ifelse(df$avg_log2FC>= 0.25 ,'Up','Down'),'NS'))
df$group2 <- factor(df$group2, levels = c("Up", "Down", "NS"))
unique(df$group2)
head(df)
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
  #xlim(c(-1, 1)) +  ##调整x轴的取值范围，max(abs(BRCA_Match_DEG$logFC))最大值是多少，再进行取舍
  geom_vline(xintercept=c(-0.25,0.25),lty=4,col="black",lwd=0.8) + #添加x轴辅助线,lty函数调整线的类型
  geom_hline(yintercept = -log10(0.05), lty=4,col="black",lwd=0.8) +  #添加y轴辅助线
  labs(x="avg_log2FC", y="-log10(p_val_adj)") +  #x、y轴标签
  ggtitle("Malignant epithelial cells DEG (Bacteria+ vs Bacteria-)") + #标题
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

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/上皮细胞/恶性非恶性差异基因/Malignant_epithelial_cells_deg_火山图.pdf",degp,width = 7,height=6)



library(org.Hs.eg.db)#物种注释包(Homo sapiens)
library(clusterProfiler)#富集分析
#添加entrez ID列：
##：symbol转entrez ID：
symbol <- deg$X
entrez <- bitr(symbol,
               fromType = "SYMBOL",#现有的ID类型
               toType = "ENTREZID",#需转换的ID类型
               OrgDb = "org.Hs.eg.db")
head(entrez)
genelist <- deg$avg_log2FC
names(genelist) <-deg$X
head(genelist)
#genelist过滤(ID转换中丢失的部分基因)：
genelist <- genelist[names(genelist) %in% entrez[,1]]
names(genelist) <- entrez[match(names(genelist),entrez[,1]),2]
length(genelist)#过滤后剩16161个基因
head(genelist)
#将genelist按照log2FC值从高到低进行排序：
genelist <- sort(genelist,decreasing = T)
head(genelist)

#GSEA_KEGG富集分析：
KEGG_ges <- gseKEGG(
  geneList = genelist,
  organism = "hsa",#不同物种选择可官方查询：https://www.genome.jp/kegg/catalog/org_list.html
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose = FALSE,
  eps = 0
)

#提取结果表：
KEGG_ges_result <- KEGG_ges@result

#显示前10列（第十一列为富集到该基因集的核心基因列表，这里过长不展示）：
KEGG_ges_result[1:3,1:10]
library(ggridges)
library(ggplot2)
library(enrichplot)
nrow(KEGG_ges_result)
p <- ridgeplot(KEGG_ges,
               showCategory = 25,
               fill = "p.adjust",
               decreasing  = T)
p

library(stringr)
library(ggplot2)
KEGG_ges_result <- KEGG_ges_result[KEGG_ges_result$padj< 0.05,]
KEGG_ges_result <- KEGG_ges_result[order(KEGG_ges_result$NES),]
KEGG_ges_result$yax <- ifelse(KEGG_ges_result$NES >0, -0.02, 0.02)
KEGG_ges_result$col <- ifelse(KEGG_ges_result$NES > 0, "blue","red")
KEGG_ges_result$NES
names(KEGG_ges_result)

p <- 
  ggplot(KEGG_ges_result,aes(y=NES,x=reorder(Description,NES),label = Description))+
  geom_text(aes(y = yax,hjust = NES > 0),size=4)+ 
  geom_bar(stat="identity",aes(fill=col))+
  geom_hline(yintercept = 0)+
  coord_flip(ylim = c(-2,2.5))+
  labs(y = "Normalised enrichment score", x = "",title="")+
  scale_fill_manual(values = c("#FB6542","#375E97"),name="Expression",labels=c("Up","Down"))+
  scale_x_discrete(breaks = NULL)+
  theme_classic(base_size =  12)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.y = element_blank(),
        legend.position = 'bottom')
p


ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/上皮细胞/恶性非恶性差异基因/Malignant_epithelial_cells_barplot_KEGG.pdf",p,width=7,height=12)


##GSEA_GO富集分析：
GO_ges <- gseGO(geneList = genelist,
                OrgDb = org.Hs.eg.db,
                ont = "ALL", #one of "BP", "MF", and "CC" subontologies, or "ALL" for all three.
                minGSSize = 10,
                maxGSSize = 500,
                pvalueCutoff = 0.05,
                eps = 0,
                verbose = FALSE)
GO_ges[1:49,1:10] #展示同样省略最后一列
nrow(GO_ges)
GO_ges_result<-GO_ges@result
p <- ridgeplot(GO_ges,
               showCategory = 15,
               fill = "p.adjust",
               decreasing  = T)
p

####可视化柱状图
library(stringr)
library(ggplot2)
GO_ges_result <- GO_ges_result[GO_ges_result$padj< 0.05,]
GO_ges_result <- GO_ges_result[order(GO_ges_result$NES),]
GO_ges_result$yax <- ifelse(GO_ges_result$NES >0, -0.02, 0.02)
GO_ges_result$col <- ifelse(GO_ges_result$NES > 0, "blue","red")
GO_ges_result$NES
names(GO_ges_result)

p <- 
  ggplot(GO_ges_result,aes(y=NES,x=reorder(Description,NES),label = Description))+
  geom_text(aes(y = yax,hjust = NES > 0),size=4)+ 
  geom_bar(stat="identity",aes(fill=col))+
  geom_hline(yintercept = 0)+
  coord_flip(ylim = c(-2,2.5))+
  labs(y = "Normalised enrichment score", x = "",title="")+
  scale_fill_manual(values = c("#FB6542","#375E97"),name="Expression",labels=c("Up","Down"))+
  scale_x_discrete(breaks = NULL)+
  theme_classic(base_size =  12)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.y = element_blank(),
        legend.position = 'bottom')
p


##GSEA_Reactome富集分析:
BiocManager::install("ReactomePA")
library(ReactomePA)
ret_ges <- gsePathway(genelist,
                      organism = "human",
                      minGSSize = 10,
                      maxGSSize = 500,
                      pvalueCutoff = 0.05,
                      pAdjustMethod = "BH",
                      verbose = FALSE,
                      eps = 0)
ret_ges_result<-ret_ges[1:34,1:10]
nrow(ret_ges)
p <- ridgeplot(ret_ges,
               showCategory = 34,
               fill = "p.adjust",
               decreasing  = T)
p
?ridgeplot

pp <- dotplot(ret_ges,
              showCategory = 34)
pp

####可视化柱状图
library(stringr)
library(ggplot2)
ret_ges_result <- ret_ges_result[order(ret_ges_result$NES),]
ret_ges_result$yax <- ifelse(ret_ges_result$NES >0, -0.02, 0.02)
ret_ges_result$col <- ifelse(ret_ges_result$NES > 0, "blue","red")
ret_ges_result$NES
names(ret_ges_result)

p <- 
  ggplot(ret_ges_result,aes(y=NES,x=reorder(Description,NES),label = Description))+
  geom_text(aes(y = yax,hjust = NES > 0),size=4)+ 
  geom_bar(stat="identity",aes(fill=col))+
  geom_hline(yintercept = 0)+
  coord_flip(ylim = c(-2,2.5))+
  labs(y = "Normalised enrichment score", x = "",title="")+
  scale_fill_manual(values = c("#FB6542","#375E97"),name="Expression",labels=c("Up","Down"))+
  scale_x_discrete(breaks = NULL)+
  theme_classic(base_size =  12)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.y = element_blank(),
        legend.position = 'bottom')
p


##GSEA_DO(Disease Ontology)富集分析:
library(DOSE)
DO_ges <- gseDO(genelist,
                minGSSize = 10,
                maxGSSize = 500,
                pvalueCutoff = 0.05,
                pAdjustMethod = "BH",
                verbose = FALSE,
                eps = 0)
DO_ges[1:3,1:10]
nrow(DO_ges)

p <- ridgeplot(DO_ges,
               showCategory = 10,
               fill = "p.adjust",
               decreasing  = T)
p


#####GSEA_MSigDB富集分析：
library(gggsea)
library(ggplot2)
library(GSVA)
library(enrichplot)
library(GSEABase)
library(fgsea)
library(clusterProfiler)
gmtfile <- "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/h.all.v2024.1.Hs.symbols.gmt"
hallmark <- read.gmt(gmtfile)
alldiff <- deg[order(deg$avg_log2FC,decreasing = T),]
id <- alldiff$avg_log2FC
names(id) <- alldiff$X

## 这里去掉了基因集前缀
hallmark$term <- gsub('HALLMARK_','',hallmark$term)
hallmark.list <- hallmark %>% split(.$term) %>% lapply( "[[", 2)

gsea.re1<-clusterProfiler::GSEA(id,TERM2GENE=hallmark,verbose = T)
g1<-as.data.frame(gsea.re1)
g1<-subset(g1,p.adjust<0.05)
g1<-g1[order(g1$NES,decreasing = T),]
type(g1)
head(g1)

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

library(gggsea)
library(ggsci)

names(hallmark.list)
se_hall<-c(head(g2$pathway,9),tail(g1$pathway,9))

#3 所选中的基因集
sig1<-g2[g2$pathway%in%se_hall,]

hallmark.se<-hallmark.list[sig1$pathway]
df.new <- gseaCurve(id, hallmark.se, sig1)

pal_line<-pal_lancet()(11)#各曲线颜色

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
            ncol =3
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


library(stringr)
library(ggplot2)
gsea.re2 <- gsea.re2[order(gsea.re2$NES),]
gsea.re2$yax <- ifelse(gsea.re2$NES >0, -0.02, 0.02)
gsea.re2$col <- ifelse(gsea.re2$NES > 0, "blue","red")
gsea.re2$NES
names(gsea.re2)

p <- ggplot(gsea.re2,aes(y=NES,x=reorder(pathway,NES),label = pathway))+
  geom_text(aes(y = yax,hjust = NES > 0),size=4)+ 
  geom_bar(stat="identity",aes(fill=col))+
  geom_hline(yintercept = 0)+
  coord_flip(ylim = c(-2,2.5))+
  labs(y = "Normalised enrichment score", x = "",title="")+
  scale_fill_manual(values = c("#FB6542","#375E97"),name="Expression",labels=c("Up","Down"))+
  scale_x_discrete(breaks = NULL)+
  theme_classic(base_size =  12)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.y = element_blank(),
        legend.position = 'bottom')
p


ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/上皮细胞/恶性非恶性差异基因/Malignant_epithelial_cells_hallmark.pdf",p_nes,width=10,height=10)

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/上皮细胞/恶性非恶性差异基因/Malignant_epithelial_cells_barplot_hallmark.pdf",p,width=10,height=10)


####非恶性细胞----

non_mepi_deg_all <- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/上皮细胞/恶性非恶性差异基因/Non-malignant epithelial cells_deg_all.csv")

deg<-non_mepi_deg_all

gene_name<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/上皮细胞/上皮细胞gene_name.csv")

intersection <- intersect(gene_name$Gene.symbol, deg$X)

deg_all_inter<-deg %>% filter(deg$X %in%intersection )

df <-deg_all_inter

df$group2<-as.factor(ifelse(df$p_val_adj < 0.05 & abs(df$avg_log2FC) >= 0.25, 
                            ifelse(df$avg_log2FC>= 0.25 ,'Up','Down'),'NS'))
df$group2 <- factor(df$group2, levels = c("Up", "Down", "NS"))
unique(df$group2)
head(df)
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
  #xlim(c(-1, 1)) +  ##调整x轴的取值范围，max(abs(BRCA_Match_DEG$logFC))最大值是多少，再进行取舍
  geom_vline(xintercept=c(-0.25,0.25),lty=4,col="black",lwd=0.8) + #添加x轴辅助线,lty函数调整线的类型
  geom_hline(yintercept = -log10(0.05), lty=4,col="black",lwd=0.8) +  #添加y轴辅助线
  labs(x="avg_log2FC", y="-log10(p_val_adj)") +  #x、y轴标签
  ggtitle("Malignant epithelial cells DEG (Bacteria+ vs Bacteria-)") + #标题
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

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/上皮细胞/恶性非恶性差异基因/non-Malignant_epithelial_cells_deg_火山图.pdf",degp,width = 7,height=6)



library(org.Hs.eg.db)#物种注释包(Homo sapiens)
library(clusterProfiler)#富集分析
#添加entrez ID列：
##：symbol转entrez ID：
symbol <- deg$X
entrez <- bitr(symbol,
               fromType = "SYMBOL",#现有的ID类型
               toType = "ENTREZID",#需转换的ID类型
               OrgDb = "org.Hs.eg.db")
head(entrez)
genelist <- deg$avg_log2FC
names(genelist) <-deg$X
head(genelist)
#genelist过滤(ID转换中丢失的部分基因)：
genelist <- genelist[names(genelist) %in% entrez[,1]]
names(genelist) <- entrez[match(names(genelist),entrez[,1]),2]
length(genelist)#过滤后剩16161个基因
head(genelist)
#将genelist按照log2FC值从高到低进行排序：
genelist <- sort(genelist,decreasing = T)
head(genelist)

#GSEA_KEGG富集分析：
KEGG_ges <- gseKEGG(
  geneList = genelist,
  organism = "hsa",#不同物种选择可官方查询：https://www.genome.jp/kegg/catalog/org_list.html
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose = FALSE,
  eps = 0
)

#提取结果表：
KEGG_ges_result <- KEGG_ges@result

#显示前10列（第十一列为富集到该基因集的核心基因列表，这里过长不展示）：
KEGG_ges_result[1:3,1:10]
library(ggridges)
library(ggplot2)
library(enrichplot)
nrow(KEGG_ges_result)
p <- ridgeplot(KEGG_ges,
               showCategory = 25,
               fill = "p.adjust",
               decreasing  = T)
p

library(stringr)
library(ggplot2)
KEGG_ges_result <- KEGG_ges_result[order(KEGG_ges_result$NES),]
KEGG_ges_result$yax <- ifelse(KEGG_ges_result$NES >0, -0.02, 0.02)
KEGG_ges_result$col <- ifelse(KEGG_ges_result$NES > 0, "blue","red")
KEGG_ges_result$NES
names(KEGG_ges_result)

p <- 
  ggplot(KEGG_ges_result,aes(y=NES,x=reorder(Description,NES),label = Description))+
  geom_text(aes(y = yax,hjust = NES > 0),size=4)+ 
  geom_bar(stat="identity",aes(fill=col))+
  geom_hline(yintercept = 0)+
  coord_flip(ylim = c(-2,2.5))+
  labs(y = "Normalised enrichment score", x = "",title="")+
  scale_fill_manual(values = c("#FB6542","#375E97"),name="Expression",labels=c("Up","Down"))+
  scale_x_discrete(breaks = NULL)+
  theme_classic(base_size =  12)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.y = element_blank(),
        legend.position = 'bottom')
p


ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/上皮细胞/恶性非恶性差异基因/non-Malignant_epithelial_cells_barplot_KEGG.pdf",p,width=7,height=12)


##GSEA_GO富集分析：
GO_ges <- gseGO(geneList = genelist,
                OrgDb = org.Hs.eg.db,
                ont = "BP", #one of "BP", "MF", and "CC" subontologies, or "ALL" for all three.
                minGSSize = 10,
                maxGSSize = 500,
                pvalueCutoff = 0.05,
                eps = 0,
                verbose = FALSE)
GO_ges[1:49,1:10] #展示同样省略最后一列
nrow(GO_ges)
GO_ges_result<-GO_ges@result
p <- ridgeplot(GO_ges,
               showCategory = 15,
               fill = "p.adjust",
               decreasing  = T)
p

####可视化柱状图
library(stringr)
library(ggplot2)
GO_ges_result <- GO_ges_result[order(GO_ges_result$NES),]
GO_ges_result$yax <- ifelse(GO_ges_result$NES >0, -0.02, 0.02)
GO_ges_result$col <- ifelse(GO_ges_result$NES > 0, "blue","red")
GO_ges_result$NES
names(GO_ges_result)

p <- 
  ggplot(GO_ges_result,aes(y=NES,x=reorder(Description,NES),label = Description))+
  geom_text(aes(y = yax,hjust = NES > 0),size=4)+ 
  geom_bar(stat="identity",aes(fill=col))+
  geom_hline(yintercept = 0)+
  coord_flip(ylim = c(-2,2.5))+
  labs(y = "Normalised enrichment score", x = "",title="")+
  scale_fill_manual(values = c("#FB6542","#375E97"),name="Expression",labels=c("Up","Down"))+
  scale_x_discrete(breaks = NULL)+
  theme_classic(base_size =  12)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.y = element_blank(),
        legend.position = 'bottom')
p


##GSEA_Reactome富集分析:
BiocManager::install("ReactomePA")
library(ReactomePA)
ret_ges <- gsePathway(genelist,
                      organism = "human",
                      minGSSize = 10,
                      maxGSSize = 500,
                      pvalueCutoff = 0.05,
                      pAdjustMethod = "BH",
                      verbose = FALSE,
                      eps = 0)
ret_ges_result<-ret_ges[1:34,1:10]
nrow(ret_ges)
p <- ridgeplot(ret_ges,
               showCategory = 34,
               fill = "p.adjust",
               decreasing  = T)
p
?ridgeplot

pp <- dotplot(ret_ges,
              showCategory = 34)
pp

####可视化柱状图
library(stringr)
library(ggplot2)
ret_ges_result <- ret_ges_result[order(ret_ges_result$NES),]
ret_ges_result$yax <- ifelse(ret_ges_result$NES >0, -0.02, 0.02)
ret_ges_result$col <- ifelse(ret_ges_result$NES > 0, "blue","red")
ret_ges_result$NES
names(ret_ges_result)

p <- 
  ggplot(ret_ges_result,aes(y=NES,x=reorder(Description,NES),label = Description))+
  geom_text(aes(y = yax,hjust = NES > 0),size=4)+ 
  geom_bar(stat="identity",aes(fill=col))+
  geom_hline(yintercept = 0)+
  coord_flip(ylim = c(-2,2.5))+
  labs(y = "Normalised enrichment score", x = "",title="")+
  scale_fill_manual(values = c("#FB6542","#375E97"),name="Expression",labels=c("Up","Down"))+
  scale_x_discrete(breaks = NULL)+
  theme_classic(base_size =  12)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.y = element_blank(),
        legend.position = 'bottom')
p


##GSEA_DO(Disease Ontology)富集分析:
library(DOSE)
DO_ges <- gseDO(genelist,
                minGSSize = 10,
                maxGSSize = 500,
                pvalueCutoff = 0.05,
                pAdjustMethod = "BH",
                verbose = FALSE,
                eps = 0)
DO_ges[1:3,1:10]
nrow(DO_ges)

p <- ridgeplot(DO_ges,
               showCategory = 10,
               fill = "p.adjust",
               decreasing  = T)
p


#####GSEA_MSigDB富集分析：
library(gggsea)
library(ggplot2)
library(GSVA)
library(enrichplot)
library(GSEABase)
library(fgsea)
library(clusterProfiler)
gmtfile <- "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/h.all.v2024.1.Hs.symbols.gmt"
hallmark <- read.gmt(gmtfile)
alldiff <- deg[order(deg$avg_log2FC,decreasing = T),]
id <- alldiff$avg_log2FC
names(id) <- alldiff$X

## 这里去掉了基因集前缀
hallmark$term <- gsub('HALLMARK_','',hallmark$term)
hallmark.list <- hallmark %>% split(.$term) %>% lapply( "[[", 2)

gsea.re1<-clusterProfiler::GSEA(id,TERM2GENE=hallmark,verbose = T)
g1<-as.data.frame(gsea.re1)
g1<-subset(g1,p.adjust<0.05)
g1<-g1[order(g1$NES,decreasing = T),]
type(g1)
head(g1)

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

library(gggsea)
library(ggsci)

names(hallmark.list)
se_hall<-c(head(g2$pathway,12),tail(g1$pathway,9))

#3 所选中的基因集
sig1<-g2[g2$pathway%in%se_hall,]

hallmark.se<-hallmark.list[sig1$pathway]
df.new <- gseaCurve(id, hallmark.se, sig1)

pal_line<-pal_lancet()(11)#各曲线颜色

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
            ncol =3
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


library(stringr)
library(ggplot2)
gsea.re2 <- gsea.re2[order(gsea.re2$NES),]
gsea.re2$yax <- ifelse(gsea.re2$NES >0, -0.02, 0.02)
gsea.re2$col <- ifelse(gsea.re2$NES > 0, "blue","red")
gsea.re2$NES
names(gsea.re2)

p <- ggplot(gsea.re2,aes(y=NES,x=reorder(pathway,NES),label = pathway))+
  geom_text(aes(y = yax,hjust = NES > 0),size=4)+ 
  geom_bar(stat="identity",aes(fill=col))+
  geom_hline(yintercept = 0)+
  coord_flip(ylim = c(-2,2.5))+
  labs(y = "Normalised enrichment score", x = "",title="")+
  scale_fill_manual(values = c("#FB6542","#375E97"),name="Expression",labels=c("Up","Down"))+
  scale_x_discrete(breaks = NULL)+
  theme_classic(base_size =  12)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.y = element_blank(),
        legend.position = 'bottom')
p


ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/上皮细胞/恶性非恶性差异基因/non-Malignant_epithelial_cells_hallmark.pdf",p_nes,width=10,height=10)

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/上皮细胞/恶性非恶性差异基因/non-Malignant_epithelial_cells_barplot_hallmark.pdf",p,width=10,height=10)



####上皮细胞响应分组分析----
# 获取所有csv文件列表
file_path = "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/上皮细胞/响应分组/"
cell_type_list <- c("Pre_NMPR","Pre_MPR","Post_MPR")
dir.create(file_path)
for (cell in cell_type_list){
  print(cell)
  Idents(Epithelial_cell)<-"MPR_Response2"
  sub_object=subset(Epithelial_cell,idents = cell)
  Idents(sub_object)<-"group_microbe2"
  deg_all=FindMarkers(sub_object, ident.1 = c("Bacteria+"), ident.2 = "Bacteria-")
  new_filename <- paste(cell, "_deg_all",".csv", sep="")
  print(new_filename)
  write.csv(deg_all,paste("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/上皮细胞/响应分组/",new_filename,sep=""))
}

# 获取所有csv文件列表

files <- list.files(path = file_path, pattern = "\\.csv$", recursive = TRUE)


# 创建一个空数据框用于存储合并后的数据
merged_data <- data.frame()
# 循环读取每个csv文件并按照某一列降序排序，进行阈值筛选，并合并到merged_data中
for (file in files) {
  data <- read.csv(paste(file_path,file,sep=""))
  sorted_data <- data[order(-data$p_val_adj), ]
  file<-basename(file)
  print(file)
  file_name <- rep(substr(file, 1, nchar(file)-12), nrow(sorted_data))
  #filtered_data <- sorted_data[sorted_data$pct.1 >= 0.25, ]
  sorted_data$group <- file_name
  merged_data <- rbind(merged_data, sorted_data)
}

# 导出合并后的数据为新的csv文件
write.csv(merged_data, "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/上皮细胞/响应分组/Epithelial_cell_merged_data.csv", row.names = FALSE)


######根据差异基因画多组火山图###----

library(ggplot2) # Create Elegant Data Visualisations Using the Grammar of Graphics
library(dplyr) # A Grammar of Data Manipulation
library(RColorBrewer) # ColorBrewer Palettes
library(grid) # The Grid Graphics Package
library(scales)

df <- read.table("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/上皮细胞/响应分组/Epithelial_cell_merged_data.csv", header = 1, check.names = F, sep = ",")
df$group <- factor(df$group, levels = c("Pre_NMPR","Pre_MPR","Post_MPR"))

df$group2<-as.factor(ifelse(df$p_val_adj < 0.05 & abs(df$avg_log2FC) >= 0.25, 
                            ifelse(df$avg_log2FC>= 0.25 ,'Up','Down'),'NS'))
df$group2 <- factor(df$group2, levels = c("Up", "Down", "NS"))
##确定添加标签的数据
df$label<-ifelse(df$p_val_adj<0.05&abs(df$avg_log2FC)>=0.25,"Y","N")
df$label<-ifelse(df$label == 'Y', as.character(df$X), '')

gene_name<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/上皮细胞/上皮细胞gene_name.csv")
intersection <- intersect(gene_name$Gene.symbol, df$X)

Imerged_data_filter_inter<- df %>% filter(df$X %in%intersection )

write.csv(Imerged_data_filter_inter,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/上皮细胞/响应分组/merged_data_filter_label.csv")


df <- Imerged_data_filter_inter

df_bg <- df %>%
  group_by(group) %>%
  summarize(max_log2FC = max(avg_log2FC),min_log2FC = min(avg_log2FC))

df_bg$group<- factor(df_bg$group, levels = c("Pre_NMPR","Pre_MPR","Post_MPR"))


library(ggrepel)

cell_type_color <- c("#FF6347","#ffc556","#1E90FF","#0000FF","#9370DB","#00FFFF","#32cd32","#006400","#FF00ff","#20B2AA")

treatment_color <- c("#4974a4","#4dae47","#f29600")

vol_plot_pre_mpr<-ggplot()+
  ##y轴正半轴的灰色背景
  geom_col(data = df_bg, 
           mapping = aes(group,max_log2FC),
           fill = "grey85", width = 0.8, alpha = 0.5) +
  ##y轴负半轴的灰色背景
  geom_col(data = df_bg, 
           mapping = aes(group, min_log2FC),
           fill = "grey85", width = 0.8, alpha = 0.5)+
  #添加各组的数据点
  geom_jitter(data = df,
              mapping = aes(x = group, y = avg_log2FC, color = group2),
              size= 2,width = 0.4, alpha = 0.7)+
  # 通过在X=0的位置添加方块进行展示分组信息，采用geom_col方法添加
  geom_col(data = df_bg,
           mapping = aes(x= group, y = 0.2, fill = group),
           width = 0.8)+
  geom_col(data = df_bg,
           mapping = aes(x= group, y = -0.2, fill = group),
           width = 0.8)+
  # 在方块中添加分组的文字信息
  geom_text(data=df_bg,
            mapping = aes(x=group, y=0, label=group),
            size = 5, color ="white",fontface = "bold")+
  #根据需要决定是添加辅助线
  # geom_hline(yintercept = 4, lty=2, color = '#ae63e4', lwd=0.8)+
  # geom_hline(yintercept = -4, lty=2, color = '#ae63e4', lwd=0.8)+
  #颜色设置
  scale_color_manual(values = c('brown','steelblue','gray'))+
  scale_fill_manual(values = treatment_color)+
  #添加显著性标签
  geom_text_repel(data = df,
                  mapping = aes(x = group, y = avg_log2FC, label = label),
                  max.overlaps = 10000,
                  size=4,
                  box.padding=unit(0.3,'lines'),
                  point.padding=unit(0.3, 'lines'),
                  segment.color='black',
                  show.legend=FALSE)+
  # 主题设置
  theme_classic()+
  labs(x = "Bacteria+ vs Bacteria-", y = "avg_log2FC", fill= NULL, color = NULL)+
  theme(axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.y = element_line(linewidth = 0.8),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 14, color = "black"),
        axis.ticks.y = element_line(linewidth = 0.8),
        legend.text = element_text(size = 12),
        axis.title.x = element_text(size = 18),  # X轴标题字体
        axis.title.y = element_text(size = 18))+
  #调整图例
  guides(color=guide_legend(override.aes = list(size=7,alpha=1)))
vol_plot_pre_mpr

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/上皮细胞/响应分组/上皮细胞火山图响应分组.jpg",vol_plot_pre_mpr,width = 8,height=7)
ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/上皮细胞/响应分组/上皮细胞火山图响应分组.pdf",vol_plot_pre_mpr,width = 8,height=7)



####GSEA 免疫细胞通路分析----

options(BioC_mirror="https://mirrors.westlake.edu.cn/bioconductor")
BiocManager::install('clusterProfiler')
BiocManager::install("UCSC.utils")
BiocManager::install("GenomeInfoDbData")
BiocManager::install("enrichplot",force =TRUE)  #画图需要
BiocManager::install("org.Hs.eg.db") #基因注释需要
BiocManager::install("GO.db")
BiocManager::install("HDO.db")
nBiocManager::install("Biostrings")
BiocManager::install("AnnotationDbi",force = TRUE)
BiocManager::install("GSVA")
BiocManager::install("fgsea")
BiocManager::install("enrichplot")
devtools::install_github("nicolash2/gggsea")

library(gggsea)
library(ggplot2)
library(GSVA)
library(enrichplot)
library(GSEABase)
library(fgsea)
library(clusterProfiler)
######GSEA_h.all.v2024.1.Hs.symbols.gmt####
##ID转换（将gene_symbol转为ENTREZID）
deg<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/免疫细胞/Immune_cell_deg_all.csv",header = T)
names(deg)

colnames(deg)[1] <- "gene_name"
head(deg)
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
type(g1)
head(g1)
View(g1$Description)

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
write.csv(g2$pathway,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/免疫细胞/免疫细胞GSEA/g2_pathway_fgseaMultilevel.csv")
write.csv(g1$Description,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/免疫细胞/免疫细胞GSEA/g2_pathway_clusterProfiler.csv")
g2$pathway
save(gsea.re1,g1,gsea.re2,g2,file = 'D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/免疫细胞/免疫细胞GSEA/gsea.RData')

####后面使用ges 画图
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
num2=6
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
se_hall<-c(head(g2$pathway,8),tail(g1$Description,10))

#3 所选中的基因集
sig1<-g2[g2$pathway%in%se_hall,]

hallmark.se<-hallmark.list[sig1$pathway]
df.new <- gseaCurve(id, hallmark.se, sig1)

#主要由三部分组成
# ggplot() + geom_gseaLine(df.new) + theme_gsea() #曲线
# ggplot() + geom_gseaTicks (df.new) + theme_gsea()#竖线
# ggplot() + geom_gseaGradient(df.new) + theme_gsea()#色块

pal_line<-pal_lancet()(10)#各曲线颜色

p_nes<-ggplot() + 
  geom_gsea(df.new,
            prettyGSEA=T,
            tickcolor='grey30',#竖线颜色
            ticksize=0.1,#线条粗细
            #colour='grey',#色块边框
            linecolor=pal_line,
            linesize=2,lty=1,
            #多个图形时使用，设定每行列图形个数
            nrow =2,
            ncol =4
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

p_nes <- ggplot(df.new, aes(x = x, y = y, color = set, group = set)) + 
  geom_gsea(df.new,
    prettyGSEA = TRUE,
    tickcolor = 'grey30',
    ticksize = 0.1,
    linecolor = pal_line,
    linesize = 2,
    lty = 1
  ) +
  theme_bw() +
  xlab(bquote(italic('Rank'))) +
  ylab(bquote(italic('Enrichment Score'))) +
  theme(
    axis.text.x = element_text(size = 10, angle = -30, face = 'plain', hjust = 0.5),
    axis.text.y = element_text(size = 10, face = 'plain', vjust = 0.5),
    axis.title.x = element_text(size = 12, face = 'plain', hjust = 0.5),
    axis.title.y = element_text(size = 12, face = 12, vjust = 0.5),
    legend.title = element_blank(),
    legend.position = "right"
  ) +
  scale_color_manual(values = pal_line)  # Ensure custom colors are applied

# Display the plot
print(p_nes)

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/免疫细胞/免疫细胞hallmark/免疫细胞hallmark_GSEA.pdf",p_nes,width = 10,height=4)

dotplot(gsea.re1, showCategory = 10, split = ".sign") + facet_grid(~.sign) +
  theme(plot.title = element_text(size = 10, color = "black", hjust = 0.5),
        axis.title = element_text(size = 10,color = "black"), 
        axis.text = element_text(size = 10,color = "black"),
        axis.text.x = element_text(angle = 0, hjust = 1 ),
        legend.position = "right",
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10))


##########################################################################
deg<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/免疫细胞/Immune_cell_deg_all.csv",header = T)
colnames(deg)[1] <- "gene_name"
geneList <- deg$avg_log2FC                 # 获取GeneList
names(geneList) <- deg$gene_name      # 对GeneList命名
geneList <- sort(geneList, decreasing = T)  
geneSet <- read.gmt("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/h.all.v2024.1.Hs.symbols.gmt")
head(geneSet)
GSEA_enrichment <- GSEA(geneList,                 # 排序后的gene
                        TERM2GENE = geneSet, # 基因集
                        pvalueCutoff = 0.05,      # P值阈值
                        minGSSize = 20,           # 最小基因数量
                        maxGSSize = 1000,         # 最大基因数量
                        eps = 0,                  # P值边界
                        pAdjustMethod = "BH")     # 校正P值的计算方法


result <- data.frame(GSEA_enrichment)
dim(GSEA_enrichment@result)
library(enrichplot) 
gseaplot2(GSEA_enrichment, "HALLMARK_MYOGENESIS", color = "red3", pvalue_table = T)
gseaplot2(GSEA_enrichment, c("VEGF_A_UP.V1_UP", "VEGF_A_UP.V1_DN"), color = c("red3", "blue4"), pvalue_table = T)
dotplot(GSEA_enrichment, showCategory = 15, color = "p.adjust")
library(ggplot2)     # 画图图
dotplot(GSEA_enrichment, showCategory = 10, split = ".sign") + facet_grid(~.sign) +
  theme(plot.title = element_text(size = 10, color = "black", hjust = 0.5),
        axis.title = element_text(size = 10,color = "black"), 
        axis.text = element_text(size = 10,color = "black"),
        axis.text.x = element_text(angle = 0, hjust = 1 ),
        legend.position = "right",
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10))

####################################################################################
####可视化柱状图
library(stringr)
library(ggplot2)
g2$process <- gsub("HALLMARK_","",g2$pathway)
g2$process <- gsub("_"," ",g2$process)
g2$process <- str_to_sentence(g2$process)
g2$process <- gsub("E2f","E2F",g2$process)
g2$process <- gsub("Ifn","IFN",g2$process)
g2$process <- gsub("G2m","G2M",g2$process)
g2$process <- gsub("Il6 jak stat3","IL6 JAK STAT3",g2$process)
g2$process <- gsub("Myc","MYC",g2$process)
g2$process <- gsub("Mtorc1","MTORC1",g2$process)
g2$process <- gsub("Tnfa signaling via nfkb","TNFA signaling via NFKB",g2$process)
g2$process <- gsub("Il2 stat","IL2 STAT",g2$process)
g2 <- g2[order(g2$NES),]
g2$yax <- ifelse(g2$NES >0, -0.02, 0.02)
g2$col <- ifelse(g2$NES > 0, "blue","red")
g2$NES
p <- 
  ggplot(g2,aes(y=NES,x=reorder(process,NES),label = process))+
  geom_text(aes(y = yax,hjust = NES > 0),size=4)+ 
  geom_bar(stat="identity",aes(fill=col))+
  geom_hline(yintercept = 0)+
  coord_flip(ylim = c(-2,2.5))+
  labs(y = "Normalised enrichment score", x = "",title="")+
  scale_fill_manual(values = c("#FB6542","#375E97"),name="Expression",labels=c("Up","Down"))+
  scale_x_discrete(breaks = NULL)+
  theme_classic(base_size =  12)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.y = element_blank(),
        legend.position = 'bottom')
p
ggsave(p,file='gsea.pdf',width = 5,height = 5)

#####GO数据库##################################################################################
library(dplyr)
library(org.Hs.eg.db) #物种注释包(Homo sapiens)
library(clusterProfiler) #富集分析
library(enrichplot)
library(ComplexHeatmap)
library(circlize)
gene <- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/免疫细胞/Immune_cell_deg_all.csv",header = T)
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
save(Go_ges,file = 'D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/免疫细胞/免疫细胞GSEA/gsea_GO.RData')
write.csv(Go_ges_result, file = c('D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/免疫细胞/免疫细胞GSEA/gsea(Go).csv'))

gseaplot2(Go_gseresult,1:7,pvalue_table = TRUE)
View(Go_ges)
my_colors <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b")
####
selected_gene_sets <- c(5, 8, 17, 25, 57, 66)
selected_gene_sets <- c("blood vessel development", "locomotion", "cell motility",
                        "T cell receptor complex", "antigen binding", "adaptive immune response")

p1 <- gseaplot2(Go_ges,
                geneSetID = selected_gene_sets , #目标基因集，也可用ID，如'hsa00830'
                color = my_colors,
                rel_heights = c(1.5, 0.5, 1), #子图高度
                subplots = 1:3, #显示哪些子图
                pvalue_table = TRUE, #是否显示pvalue表
                #title = Go_ges$Description[18],
                ES_geom = 'line')+ #'dot'将线转换为点
p1

p2 <- gseaplot2(Go_ges,
                geneSetID = selected_gene_sets,
                color = my_colors,
                rel_heights = c(1.5, 0.5, 1),
                subplots = 1:2,
                pvalue_table = T,
                title = 'Immune cell bacteria+ vs bacteria+',
                ES_geom = 'line')
p2
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
library(ggplot2)
ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/免疫细胞/免疫细胞GSEA/上调图_GO注释.pdf",p8,width =9,height=7)
ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/免疫细胞/免疫细胞GSEA/下调图_GO注释.pdf",p9,width =9,height=7)


KEGG_gseresult <- gseKEGG(gene_fc, 
                          geneList = genelist,
                          organism = "hsa",
                          minGSSize = 10,
                          maxGSSize = 500,
                          pvalueCutoff = 0.05,
                          pAdjustMethod = "BH",
                          verbose = FALSE,
                          eps = 0) #使用GSEA进行KEGG富集分析
#保存富集分析结果
gseaplot2(KEGG_gseresult,1:10,pvalue_table = TRUE)

kegg_results<-as.data.frame(KEGG_gseresult)

write.csv(kegg_results, file ="D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/免疫细胞/免疫细胞GSEA/KEGG_gseresult.csv")

####免疫细胞整体分析GO、KEGG----

library(org.Hs.eg.db)#物种注释包(Homo sapiens)
library(clusterProfiler)#富集分析

load("TCGA_CHOL_DESeq2.Rdata")
head(DESeq2)
length(rownames(DESeq2))#共16345个基因

Immune_cell_deg_all<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/免疫细胞/Immune_cell_deg_all.csv",header = T)
names(Immune_cell_deg_all)

symbol <- Immune_cell_deg_all$X
entrez <- bitr(symbol,
               fromType = "SYMBOL",#现有的ID类型
               toType = "ENTREZID",#需转换的ID类型
               OrgDb = "org.Hs.eg.db")
head(entrez)

#准备genelist文件：
##需要的genelist格式：entrez ID+log2fc

genelist <- Immune_cell_deg_all$avg_log2FC

names(genelist) <- Immune_cell_deg_all$X

head(genelist)

#genelist过滤(ID转换中丢失的部分基因)：
genelist <- genelist[names(genelist) %in% entrez[,1]]
names(genelist) <- entrez[match(names(genelist),entrez[,1]),2]
length(genelist)#过滤后剩16161个基因
head(genelist)

#将genelist按照log2FC值从高到低进行排序：
genelist <- sort(genelist,decreasing = T)
head(genelist)

#####GSEA_KEGG富集分析----：
KEGG_ges <- gseKEGG(
  geneList = genelist,
  organism = "hsa",#不同物种选择可官方查询：https://www.genome.jp/kegg/catalog/org_list.html
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose = FALSE,
  eps = 0
)

#提取结果表：
KEGG_ges_result <- KEGG_ges@result

#显示前10列（第十一列为富集到该基因集的核心基因列表，这里过长不展示）：
KEGG_ges_result[1:21,1:10]
nrow(KEGG_ges_result)

#GSEA可视化：
##山峦图：
library(ggridges)
library(ggplot2)
library(enrichplot)

p <- ridgeplot(KEGG_ges,
               showCategory = 21,
               fill = "p.adjust",
               decreasing  = T)
p

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/免疫细胞/免疫细胞KEGG/免疫细胞_ridgeplot_KEGG.pdf",p,width=7,height=7)

####GSEA_GO富集分析---
GO_ges <- gseGO(geneList = genelist,
                OrgDb = org.Hs.eg.db,
                ont = "CC", #one of "BP", "MF", and "CC" subontologies, or "ALL" for all three.
                minGSSize = 10,
                maxGSSize = 500,
                pvalueCutoff = 0.05,
                eps = 0,
                verbose = FALSE)
GO_ges[1:3,1:10] #展示同样省略最后一列
nrow(GO_ges)
p <- ridgeplot(GO_ges,
               showCategory = 15,
               fill = "p.adjust",
               decreasing  = T)
p


ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/免疫细胞/免疫细胞GO/免疫细胞_ridgeplot_GO.pdf",p,width=7,height=7)


####GSEA_Reactome富集分析----
BiocManager::install("ReactomePA")
library(ReactomePA)
ret_ges <- gsePathway(genelist,
                      organism = "human",
                      minGSSize = 10,
                      maxGSSize = 500,
                      pvalueCutoff = 0.05,
                      pAdjustMethod = "BH",
                      verbose = FALSE,
                      eps = 0)
ret_ges[1:3,1:10]
nrow(ret_ges)

p <- ridgeplot(ret_ges,
               showCategory = 15,
               fill = "p.adjust",
               decreasing  = T)
p

####GSEA_DO(Disease Ontology)富集分析----
library(DOSE)
DO_ges <- gseDO(genelist,
                minGSSize = 10,
                maxGSSize = 500,
                pvalueCutoff = 0.05,
                pAdjustMethod = "BH",
                verbose = FALSE,
                eps = 0)
DO_ges[1:3,1:10]


####GSEA_MSigDB富集分析----
#包的下载和载入：
devtools::install_github("ToledoEM/msigdf")
library(msigdf)
#提取C2注释(human)：
library(dplyr)
h <- msigdf.human %>%
  filter(category_code == "h") %>% select(geneset, symbol) %>% as.data.frame

head(h)

#这里genelist需要的是symbol：
#重新准备genelist文件：
genelist2 <- Immune_cell_deg_all$avg_log2FC
names(genelist2) <-Immune_cell_deg_all$X

#genelist过滤(ID转换中丢失的部分基因)：
genelist2 <- genelist2[names(genelist2) %in% entrez[,1]]

#将genelist按照log2FC值从高到低进行排序：
genelist2 <- sort(genelist2,decreasing = T)
head(genelist2)

alldiff <- Immune_cell_deg_all[order(Immune_cell_deg_all$avg_log2FC,decreasing = T),]
id <- alldiff$avg_log2FC
names(id) <- alldiff$gene_name

h_ges <- GSEA(genelist2,
               TERM2GENE = h,
               minGSSize = 10,
               maxGSSize = 500,
               pvalueCutoff = 0.05,
               pAdjustMethod = "BH",
               verbose = FALSE,
               eps = 0)
h_ges_result <- h_ges@result

h_ges_result[,1:10]
View(h_ges_result[,1:10])


nrow(h_ges_result)

p <- ridgeplot(h_ges,
               showCategory = 9,
               fill = "p.adjust",
               decreasing  = T)
p

library(enrichplot) 
gseaplot2(h_ges, "HOEK_NK_CELL_2011_2012_TIV_3D_VS_0DY_ADULT_3D_DN", color = "red3", pvalue_table = T)
gseaplot2(h_ges, c("VEGF_A_UP.V1_UP", "VEGF_A_UP.V1_DN"), color = c("red3", "blue4"), pvalue_table = T)
dotplot(c7_ges, showCategory = 9, color = "p.adjust")

####ggsea 可视化
library(ggsci)
col_gsea1<-pal_simpsons()(16)
num2=9
gseaplot2(h_ges,geneSetID = rownames(h_ges_result)[1:num2],
          title = "",#标题
          color = col_gsea1[1:num2],#颜色
          base_size = 14,#基础大小
          rel_heights = c(1, 0.2, 0.4),#小图相对高度
          subplots = 1:3,#展示小图
          pvalue_table = FALSE,#p值表格
          ES_geom = "line"#line or dot
)

####免疫细胞Fusobacterium----





####免疫细胞Streptococcus----





####成纤维GESA细胞有无微生物功能分析----

Idents(MPRNMPR_object_miMPRobe_remove) <- "cell_type_new"

unique(Idents(MPRNMPR_object_miMPRobe_remove))

MPRNMPR_object_miMPRobe_remove@meta.data$MPR_Response2

Fib_cell=subset(MPRNMPR_object_miMPRobe_remove,idents = c("Fibroblasts"))

Fib_cell@meta.data$group_microbe2
Idents(Fib_cell) <- "group_microbe2"
Fib_cell_deg_all=FindMarkers(Fib_cell, ident.1 = c("Bacteria+"), ident.2 = "Bacteria-")

write.csv(Fib_cell_deg_all,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/成纤维细胞/Fib_cell_deg_all.csv")

####成纤维细胞DEG----
names(df)

df <- Fib_cell_deg_all

df$group2<-as.factor(ifelse(df$p_val_adj < 0.05 & abs(df$avg_log2FC) >= 0.25, 
                            ifelse(df$avg_log2FC>= 0.25 ,'Up','Down'),'NS'))
df$group2 <- factor(df$group2, levels = c("Up", "Down", "NS"))

##确定添加标签的数据

df$label<-ifelse(df$p_val_adj<0.05&abs(df$avg_log2FC)>=0.25,"Y","N")
df$label<-ifelse(df$label == 'Y', as.character(rownames(df)), '')

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
  ggtitle("Immune cell DEG (Bacteria+ vs Bacteria-)") + #标题
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
ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/免疫细胞/免疫细胞火山图.jpg",degp,width = 10,height=7)
ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/免疫细胞/免疫细胞火山图.pdf",degp,width = 10,height=7)


####成纤维GESA hallmarker----

library(gggsea)
library(ggplot2)
library(GSVA)
library(enrichplot)
library(GSEABase)
library(fgsea)
library(clusterProfiler)
######GSEA_h.all.v2024.1.Hs.symbols.gmt####
##ID转换（将gene_symbol转为ENTREZID）
deg<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/成纤维细胞/Fib_cell_deg_all.csv",header = T)
names(deg)
colnames(deg)[1] <- "gene_name"
head(deg)
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
type(g1)
head(g1)
View(g1$Description)

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
write.csv(g2$pathway,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/成纤维细胞/g2_pathway_fgseaMultilevel.csv")
write.csv(g1$Description,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/成纤维细胞/g2_pathway_clusterProfiler.csv")

####后面使用ges 画图
####可视化
library(ggsci)
col_gsea1<-pal_simpsons()(16)

####设置多个名称
num2=6
gseaplot2(gsea.re1,geneSetID = rownames(g1)[1:num2],
          title = "",#标题
          color = col_gsea1[1:num2],#颜色
          base_size = 14,#基础大小
          rel_heights = c(1, 0.2, 0.4),#小图相对高度
          subplots = 1:3,#展示小图
          pvalue_table = FALSE,#p值表格
          ES_geom = "line"#line or dot
)


####ggsea 可视化
names(hallmark.list)
se_hall<-c(head(g2$pathway,10),tail(g1$pathway,9))

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

print(p_nes)

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/免疫细胞/免疫细胞GSEA/免疫细胞hallmark_GSEA.pdf",p_nes,width = 8,height=10)






####免疫细胞GESA分组----

library(enrichplot)
library(ComplexHeatmap)
library(circlize)
gene <- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/免疫细胞/merged_data_filter_label.csv",header = T)
gene_Post_MPR<-gene[gene$group =="Post_MPR",] 
colnames(gene_Post_MPR)[2] <- "gene_name"
gene_Post_MPR$regulate <-ifelse(gene_Post_MPR$avg_log2FC>0,"up","down")
gene_ID=bitr(gene_Post_MPR$gene_name,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db")
colnames(gene_ID)
View(gene_ID)
table(is.na(gene_ID))
gene_df <- base::merge(gene_ID,gene_Post_MPR,by.x="SYMBOL",by.y="gene_name")

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
saveRDS(Go_ges,file = 'D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/免疫细胞/免疫细胞GSEA/gsea_GO_Post_MPR.rds')
write.csv(Go_ges_result, file = c('D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/免疫细胞/免疫细胞GSEA/gsea(Go)_Post_MPR.csv'))

readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/免疫细胞/免疫细胞GSEA/gsea_GO_Post_MPR.rds")

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/免疫细胞/免疫细胞GSEA/免疫细胞hallmark_GSEA.pdf",p_nes,width = 8,height=10)

selected_gene_sets<-c(1,6,13,15,26,52)
post<-c("inflammatory response","tumor necrosis factor superfamily cytokine production",
  "regulation of tumor necrosis factor superfamily cytokine production",
  "B cell receptor signaling pathway","antigen receptor-mediated signaling pathway",
  "immune response-activating cell surface receptor signaling pathway")
my_colors <- c("#FF6F61", "#6B5B95", "#88B04B", "#F7CAC9", "#92A8D1", "#955251")

library(GseaVis)
#指定名称：图纸
selected_gene_sets <- c(1,6,13)
p1<-gseaNb(object = Go_ges,
           geneSetID =selected_gene_sets)
p1

library(GseaVis)
#指定名称：图纸

# clsaasic with pvalue
p8<-gseaNb(object = Go_ges,
           geneSetID = selected_gene_sets,
           curveCol=c("#FF6F61", "#6B5B95", "#88B04B"),
           addPval = F,
           pvalX = 0.75,pvalY = 0.8,
           pCol = 'black',
           pHjust = 0,
           subPlot = 2)
p8
selected_gene_sets=c(15,26,52)
p9<-gseaNb(object = Go_ges,
           geneSetID = selected_gene_sets,
           curveCol=c("#F7CAC9", "#92A8D1", "#955251"),
           addPval = F,
           pvalX = 0.75,pvalY = 0.8,
           pCol = 'black',
           pHjust = 0,
           subPlot = 2)
p9
library(ggplot2)
ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/免疫细胞/免疫细胞GSEA/分组上调图_GO注释.pdf",p8,width =9,height=7)
ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/免疫细胞/免疫细胞GSEA/分组下调图_GO注释.pdf",p9,width =9,height=7)


######################################################
gene <- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/免疫细胞/merged_data_filter_label.csv",header = T)
gene_Pre_MPR<-gene[gene$group =="Pre_MPR",] 
colnames(gene_Pre_MPR)[2] <- "gene_name"
gene_Pre_MPR$regulate <-ifelse(gene_Pre_MPR$avg_log2FC>0,"up","down")
gene_ID=bitr(gene_Pre_MPR$gene_name,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db")
colnames(gene_ID)
View(gene_ID)
table(is.na(gene_ID))
gene_df <- base::merge(gene_ID,gene_Pre_MPR,by.x="SYMBOL",by.y="gene_name")

#按照LOG_FOIDCHANGE进行排序
colnames(gene_df)
gene_df_sort <- gene_df[order(gene_df$avg_log2FC, decreasing = T),]
gene_fc <- gene_df_sort$avg_log2FC
names(gene_fc) <- gene_df_sort$ENTREZID

#GSEA分析——以GO数据库为例
Go_gseresult <- gseGO(gene_fc, 'org.Hs.eg.db', keyType = "ENTREZID", ont="all", 
                      minGSSize = 10, maxGSSize = 1000, pvalueCutoff=1,
                      pAdjustMethod = "BH",
                      verbose = FALSE,
                      eps = 0)

Go_ges <- setReadable(Go_gseresult,
                      OrgDb = org.Hs.eg.db,
                      keyType = "ENTREZID")

Go_ges_result <- Go_ges@result
save(Go_ges,file = 'D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/免疫细胞/免疫细胞GSEA/gsea_GO_Pre_MPR.RData')
write.csv(Go_ges_result, file = c('D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/免疫细胞/免疫细胞GSEA/gsea(Go)_Pre_MPR.csv'))


gene <- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/免疫细胞/merged_data_filter_label.csv",header = T)
gene_Pre_NMPR<-gene[gene$group =="Pre_NMPR",] 
colnames(gene_Pre_NMPR)[2] <- "gene_name"
gene_Pre_NMPR$regulate <-ifelse(gene_Pre_NMPR$avg_log2FC>0,"up","down")
gene_ID=bitr(gene_Pre_NMPR$gene_name,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db")
colnames(gene_ID)
View(gene_ID)
table(is.na(gene_ID))
gene_df <- base::merge(gene_ID,gene_Pre_NMPR,by.x="SYMBOL",by.y="gene_name")

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
save(Go_ges,file = 'D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/免疫细胞/免疫细胞GSEA/gsea_GO_Pre_NMPR.RData')
write.csv(Go_ges_result, file = c('D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/免疫细胞/免疫细胞GSEA/gsea(Go)_Pre_NMPR.csv'))

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
library(ggplot2)
ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/免疫细胞/免疫细胞GSEA/上调图_GO注释.pdf",p8,width =9,height=7)
ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/免疫细胞/免疫细胞GSEA/下调图_GO注释.pdf",p9,width =9,height=7)


####表皮细胞GSEA分析----

library(dplyr)
library(org.Hs.eg.db) #物种注释包(Homo sapiens)
library(clusterProfiler) #富集分析
library(enrichplot)
library(ComplexHeatmap)
library(circlize)

gene <- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/上皮细胞/Epithelial_cell_deg_all.csv",header = T)
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
save(Go_ges,file = 'D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/上皮细胞/上皮细胞GSEA/上皮gsea_GO.RData')
write.csv(Go_ges_result, file = c('D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/上皮细胞/上皮细胞GSEA/上皮gsea(Go).csv'))

gseaplot2(Go_gseresult,1:7,pvalue_table = TRUE)
View(Go_ges)
my_colors <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b")
####
selected_gene_sets <- c(5, 8, 17, 25, 57, 66)
selected_gene_sets <- c("blood vessel development", "locomotion", "cell motility",
                        "T cell receptor complex", "antigen binding", "adaptive immune response")

library(GseaVis)
#指定名称：图纸
selected_gene_sets <- c(5, 8, 17)
p1<-gseaNb(object = Go_ges,
           geneSetID =selected_gene_sets)

p1
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
library(ggplot2)
ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/上皮细胞/上皮细胞GSEA/上调图_GO注释.pdf",p8,width =9,height=7)
ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/上皮细胞/上皮细胞GSEA/下调图_GO注释.pdf",p9,width =9,height=7)


####上皮细胞hallmark----

library(gggsea)
library(ggplot2)
library(GSVA)
library(enrichplot)
library(GSEABase)
library(fgsea)
library(clusterProfiler)
##ID转换（将gene_symbol转为ENTREZID）

deg<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/上皮细胞/Epithelial_cell_deg_all.csv",header = T)
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
nrow(g1)

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
g2$padj
write.csv(gsea.re1@result,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/上皮细胞/上皮细胞GSEA/上皮_g2_pathway_clusterProfiler.csv")
gsea.re2_df<- as.data.frame(gsea.re2)
sapply(gsea.re2_df, class)
gsea.re2_df$leadingEdge <- sapply(gsea.re2_df$leadingEdge, toString)
write.csv(gsea.re2_df,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/上皮细胞/上皮细胞GSEA/上皮_g2_pathway_fgseaMultilevel.csv")
View(gsea.re2)
g2$pathway
save(gsea.re1,g1,gsea.re2,g2,file = 'D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/上皮细胞/上皮细胞GSEA/上皮_gsea.RData')

####后面使用ges 画图
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
num2=6
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
se_hall<-c(head(g2$pathway,11),tail(g1$pathway,11))

#3 所选中的基因集
sig1<-g2[g2$pathway%in%se_hall,]

hallmark.se<-hallmark.list[sig1$pathway]
df.new <- gseaCurve(id, hallmark.se, sig1)

#主要由三部分组成
# ggplot() + geom_gseaLine(df.new) + theme_gsea() #曲线
# ggplot() + geom_gseaTicks (df.new) + theme_gsea()#竖线
# ggplot() + geom_gseaGradient(df.new) + theme_gsea()#色块

pal_line<-pal_lancet()(11)#各曲线颜色

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
ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/上皮细胞/上皮细胞GSEA/上皮细胞hallmark_GSEA.pdf",p_nes,width = 8,height=10)


dotplot(gsea.re1, showCategory = 10, split = ".sign") + facet_grid(~.sign) +
  theme(plot.title = element_text(size = 10, color = "black", hjust = 0.5),
        axis.title = element_text(size = 10,color = "black"), 
        axis.text = element_text(size = 10,color = "black"),
        axis.text.x = element_text(angle = 0, hjust = 1 ),
        legend.position = "right",
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10))

####柱状图
library(stringr)
library(ggplot2)

View(gsea.re1@result)

q<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/上皮细胞/上皮细胞GSEA/上皮_g2_pathway_clusterProfiler.csv")
names(q)
q <- q[order(q$NES),]
q$yax <- ifelse(q$NES >0,-0.02,0.02)
q$col<-ifelse(q$NES>0,"blue","red")
p<-ggplot(q,aes(y=NES,x=reorder(Description,NES),label = Description))+  
  geom_text(aes(y = yax,hjust = NES > 0),size=4)+   
  geom_bar(stat="identity",aes(fill=col),width = 0.8)+  geom_hline(yintercept = 0)+  
  coord_flip(ylim = c(-2.5,2.5))+  
  labs(y = "Normalised enrichment score", x = "",title="")+  
  scale_fill_manual(values = c("#EBB3BE","#89ad9f"),name="Expression",labels=c("Up","Down"))+  
  scale_x_discrete(breaks = NULL)+  
  theme_classic(base_size =  12)+  
  theme(panel.grid.major = element_blank(),        
        panel.grid.minor = element_blank(),        
        axis.line.y = element_blank(),        
        legend.position = 'bottom')
p
ggsave(p,file='D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/上皮细胞/上皮细胞GSEA/上皮细胞gsea_hallmark柱状图.pdf',width = 6,height = 6)

#####上皮细胞整体分析
library(org.Hs.eg.db)#物种注释包(Homo sapiens)
library(clusterProfiler)#富集分析
#添加entrez ID列：
##：symbol转entrez ID：
symbol <- deg$X
entrez <- bitr(symbol,
               fromType = "SYMBOL",#现有的ID类型
               toType = "ENTREZID",#需转换的ID类型
               OrgDb = "org.Hs.eg.db")
head(entrez)
genelist <- deg$avg_log2FC
names(genelist) <-deg$X
head(genelist)
#genelist过滤(ID转换中丢失的部分基因)：
genelist <- genelist[names(genelist) %in% entrez[,1]]
names(genelist) <- entrez[match(names(genelist),entrez[,1]),2]
length(genelist)#过滤后剩16161个基因
head(genelist)
#将genelist按照log2FC值从高到低进行排序：
genelist <- sort(genelist,decreasing = T)
head(genelist)

#GSEA_KEGG富集分析：
KEGG_ges <- gseKEGG(
  geneList = genelist,
  organism = "hsa",#不同物种选择可官方查询：https://www.genome.jp/kegg/catalog/org_list.html
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose = FALSE,
  eps = 0
)

#提取结果表：
KEGG_ges_result <- KEGG_ges@result

#显示前10列（第十一列为富集到该基因集的核心基因列表，这里过长不展示）：
KEGG_ges_result[1:3,1:10]
library(ggridges)
library(ggplot2)
library(enrichplot)
nrow(KEGG_ges_result)
p <- ridgeplot(KEGG_ges,
               showCategory = 15,
               fill = "p.adjust",
               decreasing  = T)
p

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/上皮细胞/上皮细胞整体_ridgeplot_KEGG.pdf",p,width=7,height=6)


##GSEA_GO富集分析：
GO_ges <- gseGO(geneList = genelist,
                OrgDb = org.Hs.eg.db,
                ont = "CC", #one of "BP", "MF", and "CC" subontologies, or "ALL" for all three.
                minGSSize = 10,
                maxGSSize = 500,
                pvalueCutoff = 0.05,
                eps = 0,
                verbose = FALSE)
GO_ges[1:3,1:10] #展示同样省略最后一列
nrow(GO_ges)

p <- ridgeplot(GO_ges,
               showCategory =15,
               fill = "p.adjust",
               decreasing  = T)
p

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/上皮细胞/上皮细胞整体_ridgeplot_GO.pdf",p,width=7,height=9)


##GSEA_Reactome富集分析:
BiocManager::install("ReactomePA")
library(ReactomePA)
ret_ges <- gsePathway(genelist,
                      organism = "human",
                      minGSSize = 10,
                      maxGSSize = 500,
                      pvalueCutoff = 0.05,
                      pAdjustMethod = "BH",
                      verbose = FALSE,
                      eps = 0)
ret_ges_result<-ret_ges[1:34,1:10]
nrow(ret_ges)
p <- ridgeplot(ret_ges,
               showCategory = 15,
               fill = "p.adjust",
               decreasing  = T)
p
?ridgeplot

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/上皮细胞/上皮细胞整体_ridgeplot_ReactomePA.pdf",p,width=7,height=9)

pp <- dotplot(ret_ges,
              showCategory = 34)
pp

####可视化柱状图
library(stringr)
library(ggplot2)
ret_ges_result <- ret_ges_result[order(ret_ges_result$NES),]
ret_ges_result$yax <- ifelse(ret_ges_result$NES >0, -0.02, 0.02)
ret_ges_result$col <- ifelse(ret_ges_result$NES > 0, "blue","red")
ret_ges_result$NES
names(ret_ges_result)

p <- 
  ggplot(ret_ges_result,aes(y=NES,x=reorder(Description,NES),label = Description))+
  geom_text(aes(y = yax,hjust = NES > 0),size=4)+ 
  geom_bar(stat="identity",aes(fill=col))+
  geom_hline(yintercept = 0)+
  coord_flip(ylim = c(-2,2.5))+
  labs(y = "Normalised enrichment score", x = "",title="")+
  scale_fill_manual(values = c("#FB6542","#375E97"),name="Expression",labels=c("Up","Down"))+
  scale_x_discrete(breaks = NULL)+
  theme_classic(base_size =  12)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.y = element_blank(),
        legend.position = 'bottom')
p

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/上皮细胞/上皮细胞整体_柱状图_ReactomePA.pdf",p,width=7,height=9)


##GSEA_DO(Disease Ontology)富集分析:
library(DOSE)
DO_ges <- gseDO(genelist,
                minGSSize = 10,
                maxGSSize = 500,
                pvalueCutoff = 0.05,
                pAdjustMethod = "BH",
                verbose = FALSE,
                eps = 0)
DO_ges[1:3,1:10]
nrow(DO_ges)

p <- ridgeplot(DO_ges,
               showCategory = 10,
               fill = "p.adjust",
               decreasing  = T)
p


#####GSEA_MSigDB富集分析：
library(gggsea)
library(ggplot2)
library(GSVA)
library(enrichplot)
library(GSEABase)
library(fgsea)
library(clusterProfiler)
gmtfile <- "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/h.all.v2024.1.Hs.symbols.gmt"
hallmark <- read.gmt(gmtfile)
alldiff <- deg[order(deg$avg_log2FC,decreasing = T),]
id <- alldiff$avg_log2FC
names(id) <- alldiff$X

## 这里去掉了基因集前缀
hallmark$term <- gsub('HALLMARK_','',hallmark$term)
hallmark.list <- hallmark %>% split(.$term) %>% lapply( "[[", 2)

gsea.re1<-clusterProfiler::GSEA(id,TERM2GENE=hallmark,verbose = T)
g1<-as.data.frame(gsea.re1)
g1<-subset(g1,p.adjust<0.05)
g1<-g1[order(g1$NES,decreasing = T),]
type(g1)
head(g1)
g1$Description

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


library(gggsea)

library(ggsci)

names(hallmark.list)
se_hall<-c(head(g2$pathway,13),tail(g1$pathway,9))

#3 所选中的基因集
sig1<-g2[g2$pathway%in%se_hall,]

hallmark.se<-hallmark.list[sig1$pathway]
df.new <- gseaCurve(id, hallmark.se, sig1)

pal_line<-pal_lancet()(13)#各曲线颜色

p_nes<-ggplot() + 
  geom_gsea(df.new,
            prettyGSEA=T,
            tickcolor='grey30',#竖线颜色
            ticksize=0.1,#线条粗细
            #colour='grey',#色块边框
            linecolor=pal_line,
            linesize=2,lty=1,
            #多个图形时使用，设定每行列图形个数
            nrow = 4,
            ncol =4
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

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/上皮细胞/上皮细胞Fuso/Fuso_hallmark.pdf",p_nes,width=10,height=4)



####上皮细胞hallmark分组分析----

library(enrichplot)
library(ComplexHeatmap)
library(circlize)
library(clusterProfiler)

gene <- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/上皮细胞/响应分组/merged_data_filter_label.csv",header = T)

gmtfile <- "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/h.all.v2024.1.Hs.symbols.gmt"
hallmark <- read.gmt(gmtfile)

gene_Post_MPR<-gene[gene$group =="Post_MPR",] 

alldiff <- gene_Post_MPR[order(gene_Post_MPR$avg_log2FC,decreasing = T),]
id <- alldiff$avg_log2FC
names(id) <- alldiff$gene_name

## 这里去掉了基因集前缀
hallmark$term <- gsub('HALLMARK_','',hallmark$term)
hallmark.list <- hallmark %>% split(.$term) %>% lapply( "[[", 2)

gsea.re1<-clusterProfiler::GSEA(id,TERM2GENE=hallmark,verbose = T)
g1<-as.data.frame(gsea.re1)
g1<-subset(g1,p.adjust<0.05)
g1<-g1[order(g1$NES,decreasing = T),]

write.csv(gsea.re1@result,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/上皮细胞/响应分组/上皮_g2_pathway_clusterProfiler_Post_MPR.csv")

saveRDS(gsea.re1,file = 'D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/上皮细胞/响应分组/上皮_gsea_Post_MPR.rds')

##############################################
####Pre_MPR
gene_pre_MPR<-gene[gene$group =="Pre_MPR",] 

alldiff <- gene_pre_MPR[order(gene_pre_MPR$avg_log2FC,decreasing = T),]
id <- alldiff$avg_log2FC
names(id) <- alldiff$gene_name

## 这里去掉了基因集前缀
hallmark$term <- gsub('HALLMARK_','',hallmark$term)
hallmark.list <- hallmark %>% split(.$term) %>% lapply( "[[", 2)

gsea.re1<-clusterProfiler::GSEA(id,TERM2GENE=hallmark,verbose = T)
g1<-as.data.frame(gsea.re1)
g1<-subset(g1,p.adjust<0.05)
g1<-g1[order(g1$NES,decreasing = T),]

write.csv(gsea.re1@result,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/上皮细胞/响应分组/上皮_g2_pathway_clusterProfiler_Post_MPR.csv")

saveRDS(gsea.re1,file = 'D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/上皮细胞/响应分组/上皮_gsea_Post_MPR.rds')

#############################################

gene_pre_NMPR<-gene[gene$group =="Pre_NMPR",] 

alldiff <- gene_pre_NMPR[order(gene_pre_NMPR$avg_log2FC,decreasing = T),]
id <- alldiff$avg_log2FC
names(id) <- alldiff$gene_name

## 这里去掉了基因集前缀
hallmark$term <- gsub('HALLMARK_','',hallmark$term)
hallmark.list <- hallmark %>% split(.$term) %>% lapply( "[[", 2)

gsea.re1<-clusterProfiler::GSEA(id,TERM2GENE=hallmark,verbose = T)
g1<-as.data.frame(gsea.re1)
g1<-subset(g1,p.adjust<0.05)
g1<-g1[order(g1$NES,decreasing = T),]

write.csv(gsea.re1@result,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/上皮细胞/响应分组/上皮_g2_pathway_clusterProfiler_pre_NMPR.csv")

saveRDS(gsea.re1,file = 'D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/上皮细胞/响应分组/上皮_gsea_pre_NMPR.rds')

####后面使用ges 画图
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
num2=6
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

plotGseaTable(hallmark.list[g1$pathway],
              id, 
              gsea.re1,gseaParam = 0.5,
              colwidths = c(0.5,0.2,0.1,0.1,0.1)
)

####ggsea 可视化
names(hallmark.list)
se_hall<-c(head(g1$pathway,11),tail(g1$pathway,11))

#3 所选中的基因集
sig1<-g1[g1$pathway%in%se_hall,]

hallmark.se<-hallmark.list[sig1$pathway]
df.new <- gseaCurve(id, hallmark.se, sig1)

#主要由三部分组成
# ggplot() + geom_gseaLine(df.new) + theme_gsea() #曲线
# ggplot() + geom_gseaTicks (df.new) + theme_gsea()#竖线
# ggplot() + geom_gseaGradient(df.new) + theme_gsea()#色块

pal_line<-pal_lancet()(11)#各曲线颜色

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
ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/上皮细胞/响应分组/上皮细胞hallmark_GSEA.pdf",p_nes,width = 8,height=10)

####上皮细胞Go分组分析----
library(enrichplot)
library(ComplexHeatmap)
library(circlize)
gene <- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/上皮细胞/响应分组/merged_data_filter_label.csv",header = T)
gene_Post_MPR<-gene[gene$group =="Post_MPR",] 
colnames(gene_Post_MPR)[2] <- "gene_name"
gene_Post_MPR$regulate <-ifelse(gene_Post_MPR$avg_log2FC>0,"up","down")
gene_ID=bitr(gene_Post_MPR$gene_name,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db")
colnames(gene_ID)
View(gene_ID)
table(is.na(gene_ID))
gene_df <- base::merge(gene_ID,gene_Post_MPR,by.x="SYMBOL",by.y="gene_name")

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
saveRDS(Go_ges,file = 'D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/上皮细胞/响应分组/gsea_GO_Post_MPR.rds')
write.csv(Go_ges_result, file = c('D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/上皮细胞/响应分组/gsea(Go)_Post_MPR.csv'))

######################################
gene_pre_MPR<-gene[gene$group =="Pre_MPR",] 
colnames(gene_pre_MPR)[2] <- "gene_name"
gene_pre_MPR$regulate <-ifelse(gene_pre_MPR$avg_log2FC>0,"up","down")
gene_ID=bitr(gene_pre_MPR$gene_name,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db")
colnames(gene_ID)
View(gene_ID)
table(is.na(gene_ID))
gene_df <- base::merge(gene_ID,gene_pre_MPR,by.x="SYMBOL",by.y="gene_name")

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
saveRDS(Go_ges,file = 'D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/上皮细胞/响应分组/gsea_GO_Pre_MPR.rds')
write.csv(Go_ges_result, file = c('D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/上皮细胞/响应分组/gsea(Go)_Pre_MPR.csv'))

#############################################
gene_pre_NMPR<-gene[gene$group =="Pre_NMPR",] 
colnames(gene_pre_NMPR)[2] <- "gene_name"
gene_pre_NMPR$regulate <-ifelse(gene_pre_NMPR$avg_log2FC>0,"up","down")
gene_ID=bitr(gene_pre_NMPR$gene_name,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db")
colnames(gene_ID)
View(gene_ID)
table(is.na(gene_ID))
gene_df <- base::merge(gene_ID,gene_pre_NMPR,by.x="SYMBOL",by.y="gene_name")

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
saveRDS(Go_ges,file = 'D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/上皮细胞/响应分组/gsea_GO_pre_NMPR.rds')
write.csv(Go_ges_result, file = c('D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/上皮细胞/响应分组/gsea(Go)_pre_NMPR.csv'))




readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/免疫细胞/免疫细胞GSEA/gsea_GO_Post_MPR.rds")

selected_gene_sets<-c(1,6,13,15,26,52)
post<-c("inflammatory response","tumor necrosis factor superfamily cytokine production",
        "regulation of tumor necrosis factor superfamily cytokine production",
        "B cell receptor signaling pathway","antigen receptor-mediated signaling pathway",
        "immune response-activating cell surface receptor signaling pathway")
my_colors <- c("#FF6F61", "#6B5B95", "#88B04B", "#F7CAC9", "#92A8D1", "#955251")

library(GseaVis)
#指定名称：图纸
selected_gene_sets <- c(1,6,13)
p1<-gseaNb(object = Go_ges,
           geneSetID =selected_gene_sets)
p1

library(GseaVis)
#指定名称：图纸

# clsaasic with pvalue
p8<-gseaNb(object = Go_ges,
           geneSetID = selected_gene_sets,
           curveCol=c("#FF6F61", "#6B5B95", "#88B04B"),
           addPval = F,
           pvalX = 0.75,pvalY = 0.8,
           pCol = 'black',
           pHjust = 0,
           subPlot = 2)
p8
selected_gene_sets=c(15,26,52)
p9<-gseaNb(object = Go_ges,
           geneSetID = selected_gene_sets,
           curveCol=c("#F7CAC9", "#92A8D1", "#955251"),
           addPval = F,
           pvalX = 0.75,pvalY = 0.8,
           pCol = 'black',
           pHjust = 0,
           subPlot = 2)
p9
library(ggplot2)
ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/免疫细胞/免疫细胞GSEA/分组上调图_GO注释.pdf",p8,width =9,height=7)
ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/免疫细胞/免疫细胞GSEA/分组下调图_GO注释.pdf",p9,width =9,height=7)



####ReactomePA----
BiocManager::install("ReactomePA")
library(ReactomePA)
library(clusterProfiler)

gene<-read.csv("Post_MPR_filtered_merged_data.csv")
colnames(gene)[1] <- "gene_name"
gene$regulate <-ifelse(gene$avg_log2FC>0,"up","down")
gene_up<- gene[gene$group=="T cells",]

rt1 = as.data.frame(gene_up[!duplicated(gene_up$gene_name),])
colnames(rt1)=c("symbol")
rt2 <- bitr(rt1$symbol, fromType="SYMBOL",
            toType=c("ENTREZID","ENSEMBL"),
            OrgDb="org.Hs.eg.db",
            drop = FALSE)
rt3=rt2[is.na(rt2[,"ENTREZID"])==F,]

eReac <- enrichPathway(gene = rt3$ENTREZID,
                       organism = 'human',
                       pvalueCutoff = 0.05)

barplot(eReac , color = "p.adjust", font.size = 12,
        title = "Reactome pathway", showCategory = 20)


####heatmap 验证差异基因##----

markers <- c("KRT16")
markers <- as.data.frame(markers)

markerdata <- ScaleData(MPRNMPR_object_miMPRobe, features = as.character(unique(markers$markers)), assay = "RNA")

aver_dt<-AverageExpression(markerdata,
                           features = "KRT16",
                           group.by =c('group_microbe2'))


aver_dt <- as.data.frame(aver_dt$RNA)

aver_dt<-as.data.frame(t(aver_dt))
colnames(aver_dt)<- "Average_Expression"


PDCD1_plot2<-ggplot(df2, aes(cell_type,response)) + 
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
        axis.text.y = element_text(size = 12,color = 'black',
        ),
        axis.text.x = element_text(size = 15,color = 'black',angle =40,hjust =1,vjust=1),
        legend.title = element_text(size=12, color = "black"), 
        legend.text = element_text(size=10,color = "black",angle =45),
        legend.position = "top") + 
  scale_y_discrete(position = "right")

PDCD1_plot2



#################################################################################################
#单细胞微生物细胞大类#########################################################
################################################################################################
#################################################################################
####单细胞微生物数据###
dir()
setwd("./OSCC_data/MPR_NMPR_group")
MPRNMPR_object_miMPRobe<- readRDS("./INVADE-seq/MPRNMPR_object_miMPRobe_v1_harmony.rds")
MPRNMPR_object_miMPRobe = UpdateSeuratObject(MPRNMPR_object_miMPRobe)

save(MPRNMPR_object_miMPRobe, file= "./INVADE-seq/MPRNMPR_object_miMPRobe.RData")

####umi加入胞内菌到metadata#######----

miMPRobe_data<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/PDCD1每个细胞以及微生物分组表达情况/16s_filter.combined.genus.umi.matrix.csv")

View(miMPRobe_data)

metadata<- FetchData(MPRNMPR_object_miMPRobe_remove,"cell_id")
metadata$cell_id <- rownames(metadata)

View(metadata)

metadata<- left_join(x=metadata,y=miMPRobe_data,by = join_by(cell_id==barcode))

rownames(metadata)<-metadata$cell_id

MPRNMPR_object_miMPRobe_remove<- AddMetaData(MPRNMPR_object_miMPRobe_remove,metadata = metadata)

View(MPRNMPR_object_miMPRobe_remove)

table(rownames(MPRNMPR_object_miMPRobe_remove@meta.data) %in% miMPRobe_data$barcode)
####查看微生物细胞数目
num_cells <- ncol(MPRNMPR_object_miMPRobe_remove)

table(rownames(MPRNMPR_object_miMPRobe_remove@meta.data)%in% miMPRobe_data$barcode)

####添加微生物有无分组
x <- ifelse(rownames(MPRNMPR_object_miMPRobe@meta.data) %in% miMPRobe_data$barcode,"Bacteria detected","No bacteria detected")

x <- ifelse(rownames(MPRNMPR_object_miMPRobe@meta.data) %in% miMPRobe_data$barcode,"Bacteria+","Bacteria-")

MPRNMPR_object_miMPRobe@meta.data$group_microbe2 <- x


#rownames(miMPRobe_data)<-miMPRobe_data$barcode

#MPRNMPR_object_miMPRobe<- subset(MPRNMPR_object_miMPRobe,cells =rownames(miMPRobe_data))

####微生物featureplot  提取降维信息用ggplot 画图####----
library(tidyverse)
MPRNMPR.umap.harmonyredu<- MPRNMPR_object_miMPRobe_remove@reductions$umap.harmony@cell.embeddings%>%
  as.data.frame() %>% 
  cbind(cell_type = MPRNMPR_object_miMPRobe_remove@meta.data$cell_type_new)%>% 
  cbind(Patients = MPRNMPR_object_miMPRobe_remove@meta.data$Patients) %>% 
  cbind(Treatments = MPRNMPR_object_miMPRobe_remove@meta.data$Treatments) %>% 
  cbind(MPR_Response2 = MPRNMPR_object_miMPRobe_remove@meta.data$MPR_Response2)%>% 
  cbind(microbe_detected = MPRNMPR_object_miMPRobe_remove@meta.data$group_microbe)%>%
  cbind(microbe_detected2 = MPRNMPR_object_miMPRobe_remove@meta.data$group_microbe2)
  
names(MPRNMPR.umap.harmonyredu)

###颜色设置###画小坐标轴####https://cloud.tencent.com/developer/article/1924260

cell_type_color <- c("#FF6347","#ffc556","#1E90FF","#0000FF","#9370DB","#00FFFF","#32cd32","#006400","#FF00ff","#20B2AA")

DimPlot(MPRNMPR_object_miMPRobe_remove, group.by = "group_microbe2", reduction = "umap.harmony",raster=FALSE)


MPRNMPR.umap.harmonyredu$microbe_detected2 <- factor(MPRNMPR.umap.harmonyredu$microbe_detected2,levels=c("Bacteria+","Bacteria-"))



MPRNMPR_cell_type <- ggplot(MPRNMPR.umap.harmonyredu,aes(x= umapharmony_1, y = umapharmony_2,color = microbe_detected2)) +  
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
  geom_segment(aes(x = min(MPRNMPR.umap.harmonyredu$umapharmony_1) , y = min(MPRNMPR.umap.harmonyredu$umapharmony_2) ,
                   xend = min(MPRNMPR.umap.harmonyredu$umapharmony_1) +5, yend = min(MPRNMPR.umap.harmonyredu$umapharmony_2) ),
               colour = "black", size=1,arrow = arrow(length = unit(0.5,"cm")))+ 
  geom_segment(aes(x = min(MPRNMPR.umap.harmonyredu$umapharmony_1)  , y = min(MPRNMPR.umap.harmonyredu$umapharmony_2)  ,
                   xend = min(MPRNMPR.umap.harmonyredu$umapharmony_1) , yend = min(MPRNMPR.umap.harmonyredu$umapharmony_2) + 5),
               colour = "black", size=1,arrow = arrow(length = unit(0.5,"cm"))) +
  annotate("text", x = min(MPRNMPR.umap.harmonyredu$umapharmony_1) +1.5, y = min(MPRNMPR.umap.harmonyredu$umapharmony_2) -1, label = "UMAP_1",
           color="black",size = 4, fontface="bold" ) + 
  annotate("text", x = min(MPRNMPR.umap.harmonyredu$umapharmony_1) -1, y = min(MPRNMPR.umap.harmonyredu$umapharmony_2) + 1.5, label = "UMAP_2",
           color="black",size = 4, fontface="bold" ,angle=90)

MPRNMPR_cell_type

MPRNMPR_cell_type_med <- MPRNMPR.umap.harmonyredu %>%
  group_by(cell_type) %>%
  summarise(
    umapharmony_1 = median(umapharmony_1),
    umapharmony_2 = median(umapharmony_2)
  )

library(ggrepel)
MPRNMPR_cell_type_med

MPRNMPR_cell_type2<-MPRNMPR_cell_type +geom_label_repel(aes(label=cell_type),size=4,color="black",fontface="bold",data = MPRNMPR_cell_type_med,
                                                        point.padding = NA,label.size = NA, fill = rgb(1, 1, 1, alpha = 0.2),
                                                        segment.size=0.5,force = 1,nudge_x=0.5, nudge_y = 0,direction="y",max.overlaps=50)+
  theme(legend.position = "none")

MPRNMPR_cell_type2<-MPRNMPR_cell_type +geom_label_repel(aes(label=cell_type),size=4,color="black",fontface="bold",data = MPRNMPR_cell_type_med,
                                                        point.padding = NA,label.size = NA, fill = rgb(1, 1, 1, alpha = 0.2),
                                                        segment.size=0.5,force = 1,nudge_x=0.5, nudge_y = 0,direction="y",max.overlaps=50)+
  theme(legend.position = "top")

#MPRNMPR_cell_type2

ggsave('D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/INVADE-seq/微生物featureplot/Microbe_MPRNMPR_cell_type2.jpg',MPRNMPR_cell_type2,width =8,height =8)

ggsave('./INVADE-seq/Microbe_MPRNMPR_cell_type2.pdf',MPRNMPR_cell_type2,width =8,height = 7)

####微生物featureplot----

####umi加入胞内菌到metadata#######----

MPRNMPR_object_miMPRobe<- readRDS("./INVADE-seq/MPRNMPR_object_miMPRobe_v1_harmony.rds")
MPRNMPR_object_miMPRobe = UpdateSeuratObject(MPRNMPR_object_miMPRobe)

saveRDS(MPRNMPR_object_miMPRobe, file = "./INVADE-seq/MPRNMPR_object_miMPRobe_v1_harmony.rds")

####不同分组检测到微生物umi数量百分比----

Cellratio <- table(MPRNMPR_object_miMPRobe_remove@meta.data$group_microbe2,MPRNMPR_object_miMPRobe_remove@meta.data$orig.ident)
Cellratio <- as.data.frame(Cellratio)
colnames(Cellratio) <- c("Bacteria_group","Samples","Freq")
library(reshape2)
cellper <- dcast(Cellratio,Samples~Bacteria_group, value.var = "Freq")#长数据转为宽数据
rownames(cellper) <- cellper[,1]
cellper <- cellper[,-1]

View(cellper)
rownames(cellper)

####微生物细胞数目统计boxplot 添加分组信息###
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

cellperlong<-melt(cellper, id.vars = c("MPR_Response","Treatments","MPR_Response2"), #需保留的不参与聚合的变量列名
                  measure.vars = c('Bacteria+','Bacteria-'),#需要聚合的变量s1-s10
                  variable.name = c('Bacteria_group'),#聚合变量的新列名
                  value.name = 'value')#聚合值的新列名
cellperlong$Bacteria_group

Bacteria_group_color <- c("#f4a460","#6ca6cd")
cellperlong$MPR_Response2 <- factor(cellperlong$MPR_Response2,levels=c("Pre_NMPR","Pre_MPR","Post_MPR"))
library(cowplot)
library(ggpubr)
library(ggsignif)
Bacteria_group<- cellperlong$Bacteria_group
cell_number_total<-ggplot(cellperlong,aes(x=MPR_Response2,y=value)) + #ggplot作图
  geom_boxplot(aes(fill=Bacteria_group),outlier.shape = NA,lwd= 1)+
  #geom_boxplot(outlier.colour="red", outlier.shape=7,outlier.size=1) +
  scale_fill_manual(values = Bacteria_group_color)+
  #geom_jitter(shape = 21,aes(fill=Bacteria_group),width = 0.25) + 
  #stat_summary(fun=mean, geom="point", color="grey60") +
  theme_cowplot() +
  theme(axis.text = element_text(size = 25),
        axis.line.x=element_line(linetype=1,color="black",size=1),
        axis.ticks.x=element_line(color="black",size=2,lineend = 1),
        axis.line.y=element_line(linetype=1,color="black",size=1),
        axis.ticks.y=element_line(color="black",size=2,lineend = 1),
        axis.title = element_text(size = 25),
        legend.text = element_text(size = 25),
        axis.text.x = element_text(vjust=1,hjust=1,angle=30,size=25),
        axis.text.y = element_text(vjust=1,hjust=1,angle=0,size=20),
        legend.title = element_text(size = 25),
        plot.title = element_text(size = 25),
        legend.position = 'top') + 
  labs(title = "",y='Total cell counts',x="") +
  guides(fill = guide_legend(title = NULL))

compare_means(value ~ Bacteria_group, data =cellperlong, group.by = "MPR_Response2")
  

cell_number_total2<-cell_number_total+stat_compare_means(aes(comparisons=Bacteria_group),data=cellperlong, method="wilcox.test",
                                     label="p.signif",show.legend = F,size=12,bracket.size = 1)



ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/INVADE-seq/微生物细胞数目boxplot_棒棒图/bacterial_cell_number_total.jpg",cell_number_total,height = 14,width = 11)

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/INVADE-seq/微生物细胞数目boxplot_棒棒图/bacterial_cell_number_total.pdf",cell_number_total,height = 14,width = 11)


cellperlong_bacteria_d<- cellperlong[cellperlong$Bacteria_group %in% "Bacteria+",]

my_comparisons <- list(c("Pre_MPR", "Pre_NMPR"),c("Pre_MPR", "Post_MPR"),c("Pre_NMPR", "Post_MPR"))

cellperlong_bacteria_d$MPR_Response2 <- factor(cellperlong_bacteria_d$MPR_response2,levels = c('Pre_NMPR','Pre_MPR','Post_MPR'))


cell_number_bacteria<-ggplot(cellperlong_bacteria_d,aes(x=MPR_Response2,y=value)) + #ggplot作图
  geom_boxplot(aes(fill=MPR_Response2),outlier.shape = NA,lwd= 1)+
  #geom_boxplot(outlier.colour="red", outlier.shape=7,outlier.size=1) +
  scale_fill_manual(values =c("blue","red","green"))+
  #geom_jitter(shape = 21,aes(fill=Bacteria_group),width = 0.25) + 
  #stat_summary(fun=mean, geom="point", color="grey60") +
  theme_cowplot() +
  theme(axis.text = element_text(size = 25),
        axis.line.x=element_line(linetype=1,color="black",size=1),
        axis.ticks.x=element_line(color="black",size=2,lineend = 1),
        axis.line.y=element_line(linetype=1,color="black",size=1),
        axis.ticks.y=element_line(color="black",size=2,lineend = 1),
        axis.title = element_text(size = 25),
        legend.text = element_text(size = 25),
        axis.text.x = element_text(vjust=1,hjust=1,angle=30,size=25),
        axis.text.y = element_text(vjust=1,hjust=1,angle=0,size=20),
        legend.title = element_text(size = 25),
        plot.title = element_text(size = 25),
        legend.position = 'top') + 
  labs(title = "",y='Bacteria+ cell counts',x="") +
  guides(fill = guide_legend(title = NULL))+
  stat_compare_means(comparisons=my_comparisons,method="wilcox.test",
                     show.legend = F,size=10,bracket.size = 1)


cell_number_bacteria


ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/INVADE-seq/微生物细胞数目boxplot_棒棒图/bacterial_cell_number.jpg",cell_number_bacteria,height = 11,width = 9)

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/INVADE-seq/微生物细胞数目boxplot_棒棒图/bacterial_cell_number.pdf",cell_number_bacteria,height = 11,width = 9)



####微生物counts 加入到metadata----

reads_umi_count<- read.csv("./INVADE-seq/16s_read_umi_count.csv")
rownames(reads_umi_count) <- reads_umi_count$barcode
View(reads_umi_count)


metadata<- FetchData(MPRNMPR_object_miMPRobe_remove,"cell_id")
metadata$cell_id <- rownames(metadata)

View(metadata)

metadata<- left_join(x=metadata,y=reads_umi_count,by = join_by(cell_id==barcode))

rownames(metadata)<-metadata$cell_id

MPRNMPR_object_miMPRobe_remove<- AddMetaData(MPRNMPR_object_miMPRobe_remove,metadata = metadata)

View(MPRNMPR_object_miMPRobe_remove)

table(MPRNMPR_object_miMPRobe_remove@meta.data$Reads_counts)
Idents(MPRNMPR_object_miMPRobe_remove)


Cellratio <- cbind(MPRNMPR_object_miMPRobe_remove@meta.data$orig.ident,
                   MPRNMPR_object_miMPRobe_remove@meta.data$cell_type_new,
                   MPRNMPR_object_miMPRobe_remove@meta.data$group_microbe2,
                   MPRNMPR_object_miMPRobe_remove@meta.data$MPR_Response2,
                   MPRNMPR_object_miMPRobe_remove@meta.data$CRNCR_Response2,
                   MPRNMPR_object_miMPRobe_remove@meta.data$Reads_counts,
                   MPRNMPR_object_miMPRobe_remove@meta.data$UMIs_counts)
df<- Cellratio %>% 
   as.data.frame() %>%
  drop_na()

colnames(df) <- c("Samples","cell_type_new","group_microbe2","MPR_Response2","CRNCR_Response2","Reads_counts","UMIs_counts")


####微生物棒棒图根据分组信息reads_umis_counts----

# 数据读取
#data <- read.xlsx("data.xlsx",check.names = F)

# 细胞类型因子化，固定排序
df$cell_type_new<- factor(df$cell_type_new, levels = c('T cells', 'B cells','Plasma cells', 'Myeloid cells',  
                                                               'Mast cells',"Endothelial cells","Fibroblasts","Smooth muscle cells",
                                                               "Pericytes","Epithelial cells"))

# 颜色配置
cell_type_color <- c("#FF6347","#ffc556","#1E90FF","#0000FF","#9370DB","#00FFFF","#32cd32","#006400","#FF00ff","#20B2AA")

names(df)

ggplot(df,aes(x=cell_type_new,y=round(as.numeric(UMIs_counts),1)))+
  geom_segment(aes(x=cell_type_new,
                   xend=cell_type_new,
                   y=0,
                   yend=round(as.numeric(UMIs_counts),1)))+
  geom_point(aes(x=cell_type_new,y=round(as.numeric(UMIs_counts),1),size=5))

df_cell_type <- df %>% group_by(cell_type_new) %>% summarise(UMIs_counts = sum(as.numeric(UMIs_counts)),Reads_counts=sum(sum(as.numeric(Reads_counts))))


library(ggbreak)
df_cell_type$cell_type_new<- factor(df_cell_type$cell_type_new, levels = c('T cells', 'B cells','Plasma cells', 'Myeloid cells',  
                                                       'Mast cells',"Endothelial cells","Fibroblasts","Smooth muscle cells",
                                                       "Pericytes","Epithelial cells"))

pop_cell_type<- ggplot(df_cell_type, aes(x = cell_type_new, y = as.numeric(UMIs_counts))) +
  geom_segment(aes(x = cell_type_new, xend = cell_type_new, y = 0, yend = as.numeric(UMIs_counts) - 0.5), color = "gray",size=2) +
  geom_point(aes(x = cell_type_new, y = as.numeric(UMIs_counts), color = cell_type_new), size = 8) +
  #geom_point(aes(x = cell_type_new, y = as.numeric(UMIs_counts), color = cell_type_new), shape = 1, size = 6) +
  geom_text(aes(label = as.numeric(UMIs_counts), x = cell_type_new, y =as.numeric(UMIs_counts) + 0.7), hjust = 0.4,vjust=-0.7, size = 5) + 
  scale_color_manual(values = cell_type_color) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 80000)) +
  #coord_flip()+
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(size = 0.7),
        axis.ticks = element_line(size = 0.7),
        axis.ticks.length = unit(0.1, "cm"),
        axis.text.y = element_text(size = 15, colour = "black"),
        axis.text.x = element_text(size = 10, colour = "black",angle = 0),
        axis.title.x = element_text(size = 14, colour = "black"), 
        axis.title.y = element_blank(),  
        legend.position = "none")+
  labs(y = "Bacterial UMI counts")


pop_cell_type2 <-pop_cell_type+scale_y_break(c(15000,75000),#截断位置及范围
              space = 0.3,#间距大小
             scales = 0.3)+geom_hline( 
  yintercept = 0,  # 水平线位置
  linetype = "dashed",  # 虚线
  size = 1, # 连线的粗细 
  colour="gray40" # 连线的颜色
) +coord_flip()

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/INVADE-seq/微生物细胞数目boxplot_棒棒图/microbe_cell_type.jpg",pop_cell_type2,width=14,height = 7)
ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/INVADE-seq/微生物细胞数目boxplot_棒棒图/microbe_cell_type.pdf",pop_cell_type2,width=14,height = 7)

####微生物counts 棒棒糖图----
####不同分组microbe_counts_boxplot----

microbe_group <- cbind(MPRNMPR_object_miMPRobe_remove@meta.data$cell_id,
                       MPRNMPR_object_miMPRobe_remove@meta.data$orig.ident,
                       MPRNMPR_object_miMPRobe_remove@meta.data$cell_type_new,
                       MPRNMPR_object_miMPRobe_remove@meta.data$group_microbe2,
                       MPRNMPR_object_miMPRobe_remove@meta.data$MPR_Response2,
                       MPRNMPR_object_miMPRobe_remove@meta.data$CRNCR_Response2,
                       MPRNMPR_object_miMPRobe_remove@meta.data$Reads_counts,
                       MPRNMPR_object_miMPRobe_remove@meta.data$UMIs_counts) %>%
  as.data.frame()




colnames(microbe_group) <-c("cell_id","samples","cell_type_new","group_microbe2","MPR_Response2","CRNCR_Response2","Reads_counts","UMIs_counts")

ncol(microbe_group)
nrow(microbe_group)
names(microbe_group)
View(microbe_group)
rownames(microbe_group) <- microbe_group[,1]
microbe_group2<-microbe_group %>% 
  as.data.frame() %>% 
  drop_na() %>%
  group_by(MPR_Response2,samples) %>% 
  summarise(umi_counts_sum=sum(as.numeric(UMIs_counts)))%>% as.data.frame()

microbe_group2$MPR_Response2<-factor(microbe_group2$MPR_Response2,levels = c('Pre_NMPR','Pre_MPR','Post_MPR'))

color<-  c("#4974a4","#4dae47","#f29600")

View()
library(ggpubr)
names(microbe_group2)


ggplot(microbe_group2[-c(11,8),],aes(x=MPR_Response2,umi_counts_sum))+
  geom_boxplot(aes(fill = MPR_Response2))


microbe_MPR_response_counts<-ggplot(microbe_group2[-c(11,8),],aes(x=MPR_Response2,umi_counts_sum,fill= MPR_Response2)) +
  #stat_boxplot(geom = "errorbar", linewidth=1)+
  geom_boxplot(linewidth=0.5)+
  scale_fill_manual(values =c("#4974a4","#4dae47","#f29600"))+
  #geom_jitter(width = 0.1,alpha =0.5,size=2)+
  #geom_text(data=microbe_group2,aes(x=MPR_Response2,y=mean+1.5*std,label=Letters))+
  labs(title = "", x = "", y = "Bacterial UMI counts")+
  theme_bw()+
  theme(panel.background = element_blank(),
        panel.grid = element_blank(), 
        panel.border = element_rect(color="black", linewidth=1, linetype="solid"),
        #strip.text.y = element_text(size = 12, colour = "black",face = "bold"), #分面标题文字设置
        plot.title = element_text(hjust = 0.5,size=13,face = "bold"),
        axis.text.x = element_text(size=14,colour ="black",face = "bold",angle = 0,vjust=-0.001),
        axis.text.y = element_text(size=15,colour ="black",face = "bold"),
        axis.title.y= element_text(size=15,colour ="black",face = "bold"),
        legend.position = "none")+
  geom_signif(mapping=aes(x=MPR_Response2,y=umi_counts_sum), # 不同组别的显著性
              comparisons = list(c("Pre_NMPR", "Pre_MPR"), # 哪些组进行比较
                                 c("Pre_MPR", "Post_MPR"),
                                 c("Pre_NMPR", "Post_MPR")),
              map_signif_level=T, # T显示显著性，F显示p value
              tip_length=c(0,0,0,0,0,0,0,0,0,0,0,0), # 修改显著性线两端的长短
              y_position = c(1500,1800,2000), # 设置显著性线的位置高度
              size=1, # 修改线的粗细
              textsize = 4, # 修改显著性标记的大小
               test = "wilcox.test") # 检验的类型,可以更改
 
ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/INVADE-seq/微生物细胞数目boxplot_棒棒图/microbe_MPR_response_counts.jpg",microbe_MPR_response_counts,height = 7,width = 7)
ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/INVADE-seq/微生物细胞数目boxplot_棒棒图/microbe_MPR_response_counts.pdf",microbe_MPR_response_counts,height = 7,width = 7)

compare_means(umi_counts_sum ~ MPR_Response2, data=microbe_group2[-c(11,8),])


#### 不同分组不同细胞类型微生物counts 环状条形图----
###准备数据
library(tidyverse)

microbe_counts <- cbind(MPRNMPR_object_miMPRobe_remove@meta.data$MPR_Response2,
                       MPRNMPR_object_miMPRobe_remove@meta.data$cell_type_new,
                       MPRNMPR_object_miMPRobe_remove@meta.data$UMIs_counts) %>%
  as.data.frame() %>% drop_na()


head(microbe_counts)


colnames(microbe_counts)<-c("MPR_Response2","cell_type_new","UMIs_counts")


head(microbe_counts)

microbe_counts$UMIs_counts<- as.numeric(microbe_counts$UMIs_counts)

microbe_counts2<-aggregate(UMIs_counts~MPR_Response2+cell_type_new,microbe_counts,sum)

View(microbe_counts2)


View(microbe_counts2)

library(tidyverse)
library(ggtext)
library(plotrix)
library(gg.gap)

df <-microbe_counts2 %>% 
  select(cell_type_new,UMIs_counts,MPR_Response2,id)

View(df)

group1_order <- c("Pre_NMPR", "Pre-MPR", "Post-MPR")
group2_order <- c('T cells', 'B cells','Plasma cells', 'Myeloid cells',  
                  'Mast cells',"Endothelial cells","Fibroblasts","Smooth muscle cells",
                  "Pericytes","Epithelial cells")

df$MPR_Response2 <- factor(df$MPR_Response2, levels = group1_order)
df$cell_type_new <- factor(df$cell_type_new, levels = group2_order)

library(dplyr)

sorted_df <- df %>%
  arrange(MPR_Response2,cell_type_new)

print(sorted_df)

write.csv(sorted_df,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/INVADE-seq/微生物数目分组不同细胞柱状图/微生物丰度环状条形图数据.csv")

df <- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/INVADE-seq/微生物数目分组不同细胞柱状图/微生物丰度环状条形图数据.csv",header=T)



adjust_label_angles <- function(df, id_column) {
  nBar <- nrow(df)
  angle <- 90 - 360 * (df[[id_column]] - 0.5) / nBar
  
  df$hjust <- as.numeric(angle < -90)
  df$angle <- (angle + 180) * (angle < -90) + angle * (angle >= -90)
  
  return(df)
}


label_data <- adjust_label_angles(df, "id")

col <- c("#81A88D","#0A9F9D","#E6A0C4")


# p<- ggplot(df) +
#   geom_bar(aes(x = id, y = UMIs_counts,fill=MPR_Response2),
#            stat = "identity", position = "dodge",
#            show.legend = T,alpha = .9,linewidth = 0.5) 
  
# p2<- gg.gap(plot = p,
#          segments = c(1000, 69000),
#          tick_width = 3,
#          rel_heights = c(4, 0, 2),# 设置分隔为的三个部分的宽度
#          ylim = c(0, 70000)
#   ) 
# 
# p2
  
p3<- ggplot(df) +
  geom_bar(aes(x = id, y = scale(UMIs_counts),fill=MPR_Response2),
           stat = "identity", position = "dodge",
           show.legend = T,alpha = .9,linewidth =1)+
  coord_cartesian(ylim = c(0, 75000))+
  coord_polar()+
  scale_y_continuous(limits = c(-7500,7500),expand = c(10,10)) +
  geom_point(data=microbe_counts2,aes(x = id, y = UMIs_counts,size=UMIs_counts),shape=19,alpha=0.7,show.legend = T)+
  geom_text(data=label_data, 
            aes(x=id, y=UMIs_counts, label=paste(cell_type_new,UMIs_counts,sep="\t"), hjust=hjust),
            fontface="bold",
            alpha=1, size=3, angle= label_data$angle, inherit.aes = FALSE )+
  scale_x_continuous(expand = c(0,0))+
  scale_fill_manual(labels = paste("<span style='color:",col,
                                   "'>",unique(df$MPR_Response2),
                                   "</span>"),values = col,
                    name=" ")+
  guides(size="none")+
  labs( x=NULL,y=NULL) +
  theme( panel.background = element_blank(),
         plot.background = element_blank(),
         axis.text = element_blank(),
         legend.text = element_markdown(size = 8, face = "bold"),
         legend.title = element_text(vjust = 0.5,hjust=0.5,face="bold"),
         legend.position = c(0.5,0.5),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         axis.ticks = element_blank(),
         plot.margin = margin(0,0,0,0))

p3

####微生物counts 分面柱状图----
library(tidyverse)

microbe_counts <- cbind(MPRNMPR_object_miMPRobe_remove@meta.data$MPR_Response2,
                        MPRNMPR_object_miMPRobe_remove@meta.data$cell_type_new,
                        MPRNMPR_object_miMPRobe_remove@meta.data$UMIs_counts) %>%
  as.data.frame() %>% drop_na()


head(microbe_counts)


colnames(microbe_counts)<-c("MPR_Response2","cell_type_new","UMIs_counts")


head(microbe_counts)

microbe_counts$UMIs_counts<- as.numeric(microbe_counts$UMIs_counts)

microbe_counts2<-aggregate(UMIs_counts~MPR_Response2+cell_type_new,microbe_counts,sum)

View(microbe_counts2)


View(microbe_counts2)

library(tidyverse)
library(ggtext)
library(plotrix)
library(gg.gap)

df <-microbe_counts2 %>% 
  select(cell_type_new,UMIs_counts,MPR_Response2,id)

View(df)

group1_order <- c("Pre_NMPR", "Pre_MPR", "Post_MPR")
group2_order <- c('T cells', 'B cells','Plasma cells', 'Myeloid cells',  
                  'Mast cells',"Endothelial cells","Fibroblasts","Smooth muscle cells",
                  "Pericytes","Epithelial cells")

df$MPR_Response2 <- factor(df$MPR_Response2, levels = group1_order)
df$cell_type_new <- factor(df$cell_type_new, levels = group2_order)

library(dplyr)

sorted_df <- df %>%
  arrange(MPR_Response2,cell_type_new)

print(sorted_df)

write.csv(sorted_df,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/INVADE-seq/微生物数目分组不同细胞柱状图/微生物丰度环状条形图数据.csv")

df <- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/INVADE-seq/微生物数目分组不同细胞柱状图/微生物丰度环状条形图数据.csv",header=T)



cell_type_color <- c("#FF6347","#ffc556","#1E90FF","#0000FF","#9370DB","#00FFFF","#32cd32","#006400","#FF00ff","#20B2AA")

df$cell_type_new<- factor(df$cell_type_new, levels = c('T cells', 'B cells','Plasma cells', 'Myeloid cells',  
                                                               'Mast cells',"Endothelial cells","Fibroblasts","Smooth muscle cells",
                                                               "Pericytes","Epithelial cells"))

df$MPR_Response2<- factor(df$MPR_Response2, levels = c('Pre-NMPR', 'Pre-MPR','Post-MPR'))

colnames(df)
p1_down<-ggplot(df) +
  geom_bar(aes(x =id, y = UMIs_counts,fill=cell_type_new),
           stat = "identity", position = "dodge",width=0.9,
           show.legend = T,alpha = .9,linewidth =1)+
  scale_fill_manual(values =cell_type_color)+
  #geom_jitter(width = 0.1,alpha =0.5,size=2)+
  geom_text(data=df,aes(x=id,y= UMIs_counts+200,label=UMIs_counts),size=3.5)+
  #facet_grid(MPR_Response2~., scales = "free_y") + #按variable纵向分面
  labs(title = "", x = "", y = "")+
  geom_vline(xintercept=c(11.5,23.5),linetype="dotted",size=1)+
  theme_bw()+
  coord_cartesian(ylim = c(0, 7500))+
  scale_y_continuous(breaks=c(0,2000,6000), labels =c(0,2000,6000))+
  theme(panel.background = element_blank(),
        panel.grid = element_blank(), 
        panel.border = element_rect(linewidth=1),
        strip.text.y = element_text(size = 12, colour = "black",face = "bold"), #分面标题文字设置
        plot.title = element_text(hjust = 0.5,size=13,face = "bold"),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=15,colour ="black"),
        legend.position = "none")

p1_down

p1_upper<-ggplot(df) +
  geom_bar(aes(x =id, y = UMIs_counts,fill=cell_type_new),
           stat = "identity", position = "dodge",width=0.9,
           show.legend = T,alpha = .9,linewidth =1)+
  scale_fill_manual(values =cell_type_color)+
  geom_text(data=df,aes(x=id,y= UMIs_counts+200,label=UMIs_counts),size=3.5)+
  #geom_jitter(width = 0.1,alpha =0.5,size=2)+
  #geom_text(data=df,aes(x=MPR_Response2,label=UMIs_counts))+
  #facet_grid(MPR_Response2~., scales = "free_y") + #按variable纵向分面
  labs(title = "", x = "", y = "")+
  geom_vline(xintercept=c(11.5,23.5),linetype="dotted",size=1)+
  theme_bw()+
  coord_cartesian(ylim = c(60000, 75000))+
  scale_y_continuous(breaks=c(60000, 75000), labels =c(60000, 75000))+
  theme(panel.background = element_blank(),
        panel.grid = element_blank(), 
        panel.border = element_rect(linewidth=1),
        strip.text.y = element_text(size = 12, colour = "black",face = "bold"), #分面标题文字设置
        plot.title = element_text(hjust = 0.5,size=13,face = "bold"),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=15,colour ="black"),
        legend.position = "none")

p1_upper

library(ggpubr)
p_upper_down<- ggarrange(p1_upper, 
          p1_down, 
          heights = c(2, 3), 
          widths = c(2, 0.5), 
          ncol = 1, 
          nrow = 2, 
          align = c("hv"),
          common.legend = T, 
          legend="bottom")+guides(color = guide_legend(ncol = 2))

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/INVADE-seq/微生物数目分组不同细胞柱状图/microbe_counts_bar.jpg",p_upper_down,width = 10,height = 9)
ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/INVADE-seq/微生物数目分组不同细胞柱状图/microbe_counts_bar.pdf",p_upper_down,width = 10,height = 9)

# p2<- gg.gap(plot = p1,
#            segments = c(1000, 69000),
#            tick_width = 3,
#            rel_heights = c(4, 0,2)# 设置分隔为的三个部分的宽度
#            )
# 
# p2
# 
# ?gg.gap


####胞内菌丰度堆砌图----

microbe_group3<-microbe_group %>% 
  as.data.frame() %>% 
  drop_na() 

umi_16s<- read.csv("./INVADE-seq/16s_filter.combined.genus.umi.matrix.csv",header = T,row.names =1)
ncol(umi_16s)
nrow(umi_16s)
names(umi_16s)
umi_16s_MPR<-umi_16s[match(microbe_group$cell_id,rownames(umi_16s)),] %>% drop_na()
ncol(umi_16s_MPR)
nrow(umi_16s_MPR)
View(umi_16s_MPR)


write.csv(umi_16s_MPR,"./INVADE-seq/umi_16s_MPR.csv")

View(umi_16s_MPR)


# row=as.numeric(length(row.names(umi_16s_MPR))) #行数
# col=as.numeric(length(colnames(umi_16s_MPR)))#列数
# 
# row_sum=rep(rowSums(umi_16s_MPR), col)#行求和
# row_sum_mat=matrix(row_sum, nrow = row, ncol = col)#因为前面单次重复的个数等于列数，因此需要的转化
# ncol(row_sum)
# nrow(row_sum)
# View(row_sum)
# 
# umi_16s_MPR_relative=umi_16s_MPR/t(row_sum)
# 
# rowSums(umi_16s_MPR_relative)


umi_16s_MPR_t<- t(umi_16s_MPR) %>% as.data.frame()
View(umi_16s_MPR_t)

###计算相对丰度
relative_abundance <- umi_16s_MPR_t %>%
  mutate(across(everything(),~ round(.x / sum(.x), 10)))

View(relative_abundance)
colSums(relative_abundance)

relative_abundance$sum <- rowSums(relative_abundance)

relative_abundance<- relative_abundance[order(relative_abundance$sum,decreasing = TRUE),]

write.csv(relative_abundance,"./INVADE-seq/umi_16s_MPR_relative_order.csv")

genus_top50<- relative_abundance[1:50,-ncol(relative_abundance)]
View(genus_top50)
genus_top50["Other",] <- 1-colSums(genus_top50)


library(reshape2)
genus_top50$genus <- factor(rownames(genus_top50),levels = rev(rownames(genus_top50)))
genus_top50<- melt(genus_top50,id="genus")
names(genus_top50)[2]<- "cell_id"

View(genus_top50)
names(microbe_group3)
rownames(microbe_group3)

genus_top50$samples<-microbe_group3[match(genus_top50$cell_id, rownames(microbe_group3)),2]
genus_top50$cell_type_new<-microbe_group3[match(genus_top50$cell_id, rownames(microbe_group3)),3]
#genus_top50$group_microbe2<- microbe_group3[match(genus_top50$cell_id, rownames(microbe_group3)),4]

genus_top50$MPR_Response2<- microbe_group3[match(genus_top50$cell_id, rownames(microbe_group3)),5]

genus_top50$CRNCR_Response2<- microbe_group3[match(genus_top50$cell_id, rownames(microbe_group3)),6]

library(randomcoloR)
palette <- randomColor(count = 60)  #随机生成60种颜色，其实里面有重复的
palette <- distinctColorPalette(51)

genus_top50$MPR_Response2 <- factor(genus_top50$MPR_Response2,levels=c("Pre_NMPR","Pre_MPR","Post_MPR"))

##物种丰度堆砌图
genus_inveda<-ggplot()+geom_bar(data=genus_top50,
                                  aes(x=MPR_Response2,
                                      weight=value,
                                      fill=reorder(genus,-value)),
                                  position = "fill", # ggplot2会自行计算相对丰度，无需提前计算。
                                  width=0.5)+
  #facet_grid(.~depth)+
  scale_fill_manual(values = palette)+ # 颜色与绝对丰度堆叠柱形图保持一致
  #scale_fill_manual(values = palette)+
  scale_y_continuous(#expand=c(0,0), #设置横坐标轴紧挨柱形图
    name = "Relative abundance (%)",
    limits = c(0,1),
    breaks = seq(0,1,0.25),
    labels = paste(seq(0,100,25),"%")
  )+
  theme_bw()+
  guides(fill=guide_legend(title = "Cell-associated genera",ncol = 2))+
  labs(x=NULL)+
  theme(panel.border = element_rect(color = "black", size = 1, fill = NA),
        legend.position="right",
        axis.title = element_text(face = "bold", 
                                  size = 15,colour = "black"))+
  theme(axis.text = element_text(face = "bold", 
                                 size = 15,color="black"),
        strip.text.x = element_text(face = "bold", 
                                    size =12,color="black"))+
  theme(panel.grid=element_blank())+
  theme(legend.title = element_text(face = "bold", 
                                    size =12,color="black"),
        legend.text = element_text(face = "bold", 
                                   size =10,color="black"))
genus_inveda

ggsave("./INVADE-seq/genus_inveda_top50.jpg",genus_inveda,width =12,height =7)
ggsave("./INVADE-seq/genus_inveda_top50.pdf",genus_inveda,width =9,height =7)

#####top 30

genus_top30<- relative_abundance[1:30,-ncol(relative_abundance)]
View(genus_top30)
genus_top30["Other",] <- 1-colSums(genus_top30)


library(reshape2)
genus_top30$genus <- factor(rownames(genus_top30),levels = rev(rownames(genus_top30)))
genus_top30<- melt(genus_top30,id="genus")
names(genus_top30)[2]<- "cell_id"

View(genus_top30)
names(microbe_group3)
rownames(microbe_group3)

genus_top30$samples<-microbe_group3[match(genus_top30$cell_id, rownames(microbe_group3)),2]
genus_top30$cell_type_new<-microbe_group3[match(genus_top30$cell_id, rownames(microbe_group3)),3]
#genus_top50$group_microbe2<- microbe_group3[match(genus_top50$cell_id, rownames(microbe_group3)),4]

genus_top30$MPR_Response2<- microbe_group3[match(genus_top30$cell_id, rownames(microbe_group3)),5]

genus_top30$CRNCR_Response2<- microbe_group3[match(genus_top30$cell_id, rownames(microbe_group3)),6]

library(randomcoloR)
palette <- randomColor(count = 60)  #随机生成60种颜色，其实里面有重复的
palette <- distinctColorPalette(51)

genus_top30$MPR_Response2 <- factor(genus_top30$MPR_Response2,levels=c("Pre_NMPR","Pre_MPR","Post_MPR"))

##物种丰度堆砌图
genus_inveda<-ggplot()+geom_bar(data=genus_top30,
                                aes(x=MPR_Response2,
                                    weight=value,
                                    fill=reorder(genus,-value)),
                                position = "fill", # ggplot2会自行计算相对丰度，无需提前计算。
                                width=0.5)+
  #facet_grid(.~depth)+
  scale_fill_manual(values = palette)+ # 颜色与绝对丰度堆叠柱形图保持一致
  #scale_fill_manual(values = palette)+
  scale_y_continuous(#expand=c(0,0), #设置横坐标轴紧挨柱形图
    name = "Relative abundance (%)",
    limits = c(0,1),
    breaks = seq(0,1,0.25),
    labels = paste(seq(0,100,25),"%")
  )+
  theme_bw()+
  guides(fill=guide_legend(title = "Cell-associated genera",ncol = 2))+
  labs(x=NULL)+
  theme(panel.border = element_rect(color = "black", size = 1, fill = NA),
        legend.position="right",
        axis.title = element_text(face = "bold", 
                                  size = 15,colour = "black"))+
  theme(axis.text = element_text(face = "bold", 
                                 size = 15,color="black"),
        strip.text.x = element_text(face = "bold", 
                                    size =12,color="black"))+
  theme(panel.grid=element_blank())+
  theme(legend.title = element_text(face = "bold", 
                                    size =12,color="black"),
        legend.text = element_text(face = "bold", 
                                   size =10,color="black"))
genus_inveda

ggsave("./INVADE-seq/genus_inveda_top30.jpg",genus_inveda,width =12,height =7)
ggsave("./INVADE-seq/genus_inveda_top30.pdf",genus_inveda,width =9,height =7)


####物种丰度桑葚图----
library(ggalluvial) #加载ggalluvial包
library(reshape2)
#install.packages("ggalluvial")
genus_top10<- relative_abundance[1:10,-ncol(relative_abundance)]
genus_top10["Other",] <- 1-colSums(genus_top10)
genus_top10$genus <- factor(rownames(genus_top10),levels = rev(rownames(genus_top10)))
genus_top10<- melt(genus_top10,id="genus")
names(genus_top10)[2]<- "cell_id"
names(microbe_group3)
rownames(microbe_group3)

genus_top10$samples<-microbe_group3[match(genus_top10$cell_id, rownames(microbe_group3)),2]
genus_top10$cell_type_new<-microbe_group3[match(genus_top10$cell_id, rownames(microbe_group3)),3]
genus_top10$MPR_Response2<- microbe_group3[match(genus_top10$cell_id, rownames(microbe_group3)),5]
genus_top10$CRNCR_Response2<- microbe_group3[match(genus_top10$cell_id, rownames(microbe_group3)),6]

library(randomcoloR)
palette <- randomColor(count = 60)  #随机生成60种颜色，其实里面有重复的

genus_top10_2<-genus_top10 %>% group_by(MPR_Response2,cell_type_new,genus) %>% summarise(value2=sum(value))%>% as.data.frame()
View(genus_top10_2)
# 挑选出列名为col中元素为delete的行数
del <- which(genus_top10_2$genus=="Other")
# 删除这些行
genus_top10_2 <- genus_top10_2[-del,]

genus_top10_2$MPR_Response2 <- factor(genus_top10_2$MPR_Response2,levels=c("Pre_NMPR","Pre_MPR","Post_MPR"))
genus_top10_2$cell_type_new <- factor(genus_top10_2$cell_type_new,levels=c('T cells', 'B cells','Plasma cells', 'Myeloid cells',  
                                                                           'Mast cells',"Endothelial cells","Fibroblasts","Smooth muscle cells",
                                                                           "Pericytes","Epithelial cells"))
palette <- distinctColorPalette(23)
Mulberry<-ggplot(genus_top10_2,
       aes(y =as.numeric(value2),
           axis1 = MPR_Response2, axis2 = genus, axis3 = cell_type_new))+#定义图形绘制
  geom_alluvium(aes(fill = genus),width = 0, reverse = FALSE,discern = TRUE)+#控制线条流向
  scale_fill_manual(values = palette)+
  geom_stratum(width = 1/4,fill="#F0fff0",color="blue",alpha=0.8,reverse = FALSE,discern = TRUE) +#控制中间框的宽度 
  #geom_stratum(fill ="white",color="skyblue",alpha=1,width =1/5)+
  geom_text(stat = "stratum", aes(label = after_stat(stratum)),reverse = FALSE, size = 4,angle=0,discern = TRUE)+ #定义中间的文字
  scale_x_continuous(breaks = 1:3, labels = c("Response", "Genus", "Cell_type"))+#定义X轴上图标排序
  theme(legend.position = "none") +#定义图例有无
  #ggtitle("Cell-associated genera") +#定义图名
  theme_bw()+
  labs(x=NULL,y=NULL)+
  theme(legend.position="none",
        panel.grid.major =element_blank(),
        panel.border = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_text(face = "bold", 
                                  size = 15,colour = "black"))+
  theme(axis.text = element_text(face = "bold", 
                                 size = 15,color="black"),
        strip.text.x = element_text(face = "bold", 
                                    size =12,color="black"))+
  theme(panel.grid=element_blank())+
  theme(legend.title = element_text(face = "bold", 
                                    size =12,color="black"),
        legend.text = element_text(face = "bold", 
                                   size =10,color="black"))

Mulberry
ggsave("./INVADE-seq/microbe_mulberry.jpg",Mulberry,height = 8,width = 15)

ggsave("./INVADE-seq/microbe_mulberry.pdf",Mulberry,height = 8,width = 15)

####物种丰度桑葚图----

####不同分组不同细胞类型物种丰度热图####----

#group 数据准备

microbe_group <- cbind(MPRNMPR_object_miMPRobe@meta.data$cell_id,
                       MPRNMPR_object_miMPRobe@meta.data$orig.ident,
                       MPRNMPR_object_miMPRobe@meta.data$cell_type_new,
                       MPRNMPR_object_miMPRobe@meta.data$group_microbe2,
                       MPRNMPR_object_miMPRobe@meta.data$MPR_Response2,
                       MPRNMPR_object_miMPRobe@meta.data$CRNCR_Response2,
                       MPRNMPR_object_miMPRobe@meta.data$Reads_counts,
                       MPRNMPR_object_miMPRobe@meta.data$UMIs_counts,
                       MPRNMPR_object_miMPRobe@meta.data$Patients) %>%
  as.data.frame()

ncol(microbe_group)
nrow(microbe_group)

colnames(microbe_group)<-c("cell_id","orig.ident","cell_type_new","group_microbe2",
                           "MPR_Response2","CRNCR_Response2","Reads_counts","UMIs_counts","Patients")

View(microbe_group)

####胞内菌丰度数据准备heatmap----
setwd("./MPR_NMPR_group")

umi_16s<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/INVADE-seq/16s_filter.combined.genus.umi.matrix.csv",header = T,row.names =1)
ncol(umi_16s)
nrow(umi_16s)
names(umi_16s)

umi_16s_MPR<-umi_16s[match(microbe_group$cell_id,rownames(umi_16s)),] %>% drop_na()
ncol(umi_16s_MPR)
nrow(umi_16s_MPR)
names(umi_16s)
microbe_group_heatmap<- microbe_group[match(rownames(umi_16s_MPR),microbe_group$cell_id),]
ncol(microbe_group_heatmap)
nrow(microbe_group_heatmap)
View(microbe_group_heatmap)

###计算相对丰度
umi_16s_MPR_t<- t(umi_16s_MPR) %>% as.data.frame()
relative_abundance <- umi_16s_MPR_t %>%
  mutate(across(everything(),~ round(.x / sum(.x), 10)))

View(relative_abundance)
colSums(relative_abundance)

##整理数据计算Occurrence frequency
relative_abundance_t<-t(relative_abundance) %>% as.data.frame()
View(relative_abundance_t)
colnames(microbe_group_heatmap)
relative_abundance_t$Patients <-microbe_group_heatmap[match(rownames(relative_abundance_t),microbe_group_heatmap$cell_id),9]
relative_abundance_t$cell_type_new<-microbe_group_heatmap[match(rownames(relative_abundance_t),microbe_group_heatmap$cell_id),3]

relative_abundance_t$MPR_Response2<-microbe_group_heatmap[match(rownames(relative_abundance_t),microbe_group_heatmap$cell_id),5]

genus_list<- colnames(relative_abundance_t)[c(-514,-515,-516)]
result <- data.frame()
View(relative_abundance_t)

for (i in 1:length(genus_list)){
  col_data<- aggregate(relative_abundance_t[genus_list[i]], by=list(type=relative_abundance_t$MPR_Response2,relative_abundance_t$Patients),sum)
  col_name <- paste0(genus_list[i])
  if(i == 1) {
    result <- data.frame(col_data)
  } else {
    result <- cbind(result, col_data[[genus_list[i]]])
  }
  colnames(result)[i+2] <- col_name
}

View(result)

write.csv(result,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/INVADE-seq/微生物丰度热图/每个病人物种丰度.csv")

result_Occurrence <- result %>%
  group_by(type,Group.2) %>%
  summarise(across(everything(), ~sum(. > 0, na.rm = TRUE)))

View(result_Occurrence)
genus_list<- colnames(result_Occurrence)[c(-1,-2)]
result_Occurrence_frequency <- result_Occurrence %>%select(1,2,where(is.numeric)) %>% 
  group_by(type,Group.2) %>%
  summarize_at(vars(genus_list), sum) %>% as.data.frame()

View(result_Occurrence_frequency)


rownames(result_Occurrence_frequency)<-paste(result_Occurrence_frequency$type,result_Occurrence_frequency$Group.2,sep="_")
rownames(result_Occurrence_frequency)<-c("Post_MPR_occurrence_frequency","Pre_MPR_occurrence_frequency","Pre_NMPR_occurrence_frequency")
result_Occurrence_frequency <- result_Occurrence_frequency[,-c(1,2)]

# 使用 apply() 函数将每行的元素逐个除以相应的数字

#分组个数
divisor <- c(4, 6, 4)
#每个病人除以每个病人divisor
divisor <-c(rep(1,14))
result_Occurrence_frequency2 <- apply(t(result_Occurrence_frequency), 1, function(x) x / divisor)
View(result_Occurrence_frequency2)

result_Occurrence_frequency3 <- t(result_Occurrence_frequency2 *100) %>% as.data.frame()
View(result_Occurrence_frequency3)

###物种总数##

relative_abundance$sum <- rowSums(relative_abundance)
# 合并Occurrence
View(relative_abundance)

relative_abundance$Post_MPR_occurrence_frequency <- result_Occurrence_frequency3[match(rownames(relative_abundance),rownames(result_Occurrence_frequency3)),1]
relative_abundance$Pre_MPR_occurrence_frequency <- result_Occurrence_frequency3[match(rownames(relative_abundance),rownames(result_Occurrence_frequency3)),2]
relative_abundance$Pre_NMPR_occurrence_frequency <- result_Occurrence_frequency3[match(rownames(relative_abundance),rownames(result_Occurrence_frequency3)),3]

colnames(result_Occurrence_frequency3)

###根据多个分组降序
relative_abundance_sorted <- relative_abundance %>%
  arrange(desc(sum), desc(Post_MPR_occurrence_frequency), desc(Pre_MPR_occurrence_frequency),desc(Pre_NMPR_occurrence_frequency))

View(relative_abundance_sorted)

write.csv(relative_abundance_sorted,"relative_abundance_sorted.csv")

##只根据sum 降序
relative_abundance_sorted_sum<- relative_abundance[order(relative_abundance$sum,decreasing = TRUE),]
View(relative_abundance_sorted_sum)
write.csv(relative_abundance_sorted_sum,"relative_abundance_sorted_sum.csv")

genus_top30<- relative_abundance[1:30,-ncol(relative_abundance)]
View(genus_top30)


###提取分组信息##
#手动标准化

umi_16s_MPR_top<- umi_16s_MPR[match(rownames(genus_top30),colnames(umi_16s_MPR),)]
umi_16s_MPR_top$cell_type_new <- microbe_group_heatmap[match(rownames(umi_16s_MPR_top), microbe_group_heatmap$cell_id),3]
umi_16s_MPR_top$MPR_Response2 <- microbe_group_heatmap[match(rownames(umi_16s_MPR_top), microbe_group_heatmap$cell_id),5]

umi_16s_MPR_heatmap<-umi_16s_MPR_top %>% group_by(MPR_Response2,cell_type_new) %>% 
  summarise_all(sum, na.rm=TRUE)
View(umi_16s_MPR_heatmap)


normalize <- function(x) {
  return ((x - min(x)) / (max(x) - min(x)))
}
umi_16s_MPR_heatmap=as.data.frame(lapply(umi_16s_MPR_heatmap[,3:30],normalize))%>% cbind(umi_16s_MPR_heatmap[,c(1,2)])
View(umi_16s_MPR_heatmap)


umi_16s_MPR_heatmap<-scale(umi_16s_MPR_heatmap[,c(-1,-2)], center = TRUE, scale = TRUE)%>% cbind(umi_16s_MPR_heatmap[,c(1,2)])
View(umi_16s_MPR_heatmap)

#for (i in 1:nrow(umi_16s_MPR_heatmap)) umi_16s_MPR_heatmap[i, ] <- scale(log(unlist(umi_16s_MPR_heatmap[, ]) + 1, 2))

#for (i in 3:ncol(umi_16s_MPR_heatmap)) umi_16s_MPR_heatmap[,i] <- scale(log(unlist(umi_16s_MPR_heatmap[,i]) + 1, 2))

View(umi_16s_MPR_heatmap)


umi_16s_MPR_heatmap_scaled <- as.matrix(umi_16s_MPR_heatmap)


write.csv(umi_16s_MPR_heatmap,"./umi_16s_MPR_heatmap_top30_scaled.csv")


####微生物丰度热图top30----

setwd("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/INVADE-seq/微生物丰度热图")
library(tidyverse)
library(ggnewscale)
library(MetBrewer)
library(patchwork)
library(ggtext)
df <- read_tsv("./INVADE-seq/微生物丰度热图/data.xls") %>% 
  select(where(is.character)) %>%
  pivot_longer(-sample.ID)


df_top <- read.csv("./umi_16s_MPR_heatmap_top30_scaled.csv",row.names = 1,header = T) 

df <- read.csv("./umi_16s_MPR_heatmap_top30_scaled.csv") %>% 
  select(where(is.character)) %>%
  pivot_longer(-sample_id)

View(df)
ncol(df)
nrow(df)

df$sample <- factor(df$sample_id,levels = df$sample_id %>% unique())
df$sample <- factor(df$sample_id,levels = c("S1","S2","S3","S4","S5","S6","S7","S8","S9","S10",
                                            "S11","S12","S13","S14","S15","S16","S17","S18","S19","S20",
                                            "S21","S22","S23","S24","S25","S26","S27","S28","S29","S30"))
df$name <- factor(df$name,levels = df$name %>% unique() %>% rev())
View(df)
df1 <- df %>% filter(name=="MPR_Response2") %>% dplyr::rename("MPR_Response2"="value")
df2 <- df %>% filter(name=="cell_type_new") %>% 
  dplyr::rename("cell_type_new"="value")
View(df2)

df5 <- read.csv("./umi_16s_MPR_heatmap_top30_scaled.csv") %>% select(1,where(is.numeric)) %>% 
  pivot_longer(-sample_id) 
df5$value<- scale(df5$value)

#mutate(value=log(value)) %>% 
#dplyr::rename("median centered<br>log-transformed FPKM"="value")
nrow(df5)
ncol(df5)
View(df1)
View(df2)
View(df)
cell_type_color <- c("#FF6347","#ffc556","#1E90FF","#0000FF","#9370DB","#00FFFF","#32cd32","#006400","#FF00ff","#20B2AA")
response_color <- c("#4974a4","#4dae47","#f29600")

df1$MPR_Response2<- factor(df1$MPR_Response2,levels=c("Pre_NMPR","Pre_MPR","Post_MPR"))

df2$cell_type_new<- factor(df2$cell_type_new,levels = c('T cells', 'B cells','Plasma cells', 'Myeloid cells',  
                                                        'Mast cells',"Endothelial cells","Fibroblasts","Smooth muscle cells",
                                                        "Pericytes","Epithelial cells"))
df2$sample_id<- factor(df2$sample_id,levels = c("S1","S2","S3","S4","S5","S6","S7","S8","S9","S10",
                                               "S11","S12","S13","S14","S15","S16","S17","S18","S19","S20",
                                               "S21","S22","S23","S24","S25","S26","S27","S28","S29","S30"))

df1$sample_id<- factor(df1$sample_id,levels = c("S1","S2","S3","S4","S5","S6","S7","S8","S9","S10",
                                                "S11","S12","S13","S14","S15","S16","S17","S18","S19","S20",
                                                "S21","S22","S23","S24","S25","S26","S27","S28","S29","S30"))
df5$sample_id<- factor(df5$sample_id,levels = c("S1","S2","S3","S4","S5","S6","S7","S8","S9","S10",
                                                "S11","S12","S13","S14","S15","S16","S17","S18","S19","S20",
                                                "S21","S22","S23","S24","S25","S26","S27","S28","S29","S30"))

df5$name<- factor(df5$name,levels=rev(c("Massilia", "Herbaspirillum","Corynebacterium","Sphingomonas",               
                                        "Treponema","Acinetobacter","Achromobacter",              
                                        "Fusobacterium","Arthrobacter","Rothia",                     
                                        "Lawsonella","Rhodococcus","Blattabacterium",            
                                        "Methylobacterium","Duganella" ,"Dietzia",                    
                                        "unclassified.Prevotellaceae","Streptococcus","Prevotella" ,                
                                        "Aerococcus","Planococcus","Staphylococcus",             
                                        "Campylobacter","Capnocytophaga","Bradyrhizobium",            
                                        "Propionibacterium","Lactobacillus","Actinomyces",                
                                        "Blastococcus","Belnapia")))

View(df5)
p1 <- df %>% ggplot(aes(x=sample,y=name))+
  geom_tile(data=df2,
            color="black",aes(fill=cell_type_new))+
  scale_fill_manual(values=cell_type_color)+
  new_scale_fill()+
  geom_tile(data=df1,color="black",aes(fill=MPR_Response2))+
  scale_fill_manual(values=response_color)+
  new_scale_fill()+
  labs(x=NULL,y=NULL)+
  theme(axis.text.x = element_blank(),
        axis.text.y=element_text(color="black"),
        #legend.position = c(1.6,-2),
        legend.location = "left",
        legend.direction = "vertical",
        legend.spacing.x =unit(0.1,"cm"),
        legend.spacing.y = unit(0.1,"cm"),
        legend.key.height = unit(0.5,"cm"),
        legend.key.width = unit(0.5,"cm"),
        plot.margin=unit(c(0,1,-20,0.5),unit="cm"),
        plot.background = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank(),
        legend.text = element_text(color="black"))

p1
col = rainbow(4)
names(col) <- letters[1:4]
#met.brewer("Cassatt1")
##00Bfff
#colorRampPalette(colors = c("blue","white","red"))(100)
p2 <- df %>% ggplot(aes(x=sample_id,y=name))+
  geom_tile(data=df5,color="black",
            aes(fill=value))+
  scale_fill_gradientn(colors = c('#00Bfff', "yellow", 'red'),na.value = NA)+
  new_scale_fill()+
  labs(x= NULL,y=NULL)+
  theme(axis.text.x = element_blank(),
        axis.text.y=element_text(color="black"),
        #legend.position = c(1.6,-2),
        legend.location = "left",
        legend.direction = "vertical",
        plot.margin=unit(c(-20,1,0.5,0.5),unit="cm"),
        plot.background = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank(),
        legend.title =element_markdown(),
        legend.text = element_text(color="black"))

p2


plot <- p1/p2+plot_layout(heights = c(0.1,2))

plot
ggsave(plot,file="D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/INVADE-seq/微生物丰度热图/heatmap3分组颜色改变.pdf",dpi=300,width =9,height=10)


####每个病人Top30物种丰度堆砌图----

patient_microbe<-read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/INVADE-seq/微生物丰度热图/每个病人物种丰度.csv")
rownames(patient_microbe)<- paste(patient_microbe$type,patient_microbe$Group.2,sep="-")
colnames(patient_microbe)
patient_microbe<- patient_microbe[,-c(1,2)]
relative_abundance<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/INVADE-seq/物种occurance frequency/relative_abundance_sorted_sum只根据sum降序丰度.csv")
colnames(relative_abundance)[1]
genus_top30<- relative_abundance$X[1:30]


patient_microbe_abundance<-patient_microbe[,which(colnames(patient_microbe)%in%genus_top30)]

patient_microbe_abundance<-patient_microbe_abundance[,match(genus_top30,colnames(patient_microbe_abundance))]

#patient_microbe_abundance[,"Other"] <- rowSums(patient_microbe)-rowSums(patient_microbe_abundance)

patient_microbe_abundance <- patient_microbe_abundance %>%
  mutate(across(everything(), ~ . / sum(.) * 100))

View(patient_microbe_abundance)
library(reshape2)
patient_microbe_abundance$patient <- factor(rownames(patient_microbe_abundance),levels = rev(rownames(patient_microbe_abundance)) )
View(patient_microbe_abundance)
patient_microbe_abundance<- melt(patient_microbe_abundance,id="patient")
names(patient_microbe_abundance)
# df<- patient_microbe_abundance %>%
#   separate(full_name, into = c("first_name", "last_name"), sep = " ")
patient_microbe_abundance$patient <- factor(patient_microbe_abundance$patient,levels=c("Pre_NMPR-P3","Pre_NMPR-P5",
                                                                                       "Pre_NMPR-P6","Pre_NMPR-P10",
                                                                                       "Pre_MPR-P1","Pre_MPR-P2",
                                                                                       "Pre_MPR-P4","Pre_MPR-P7",
                                                                                       "Pre_MPR-P8","Pre_MPR-P9",
                                                                                       "Post_MPR-P1","Post_MPR-P4",
                                                                                       "Post_MPR-P7","Post_MPR-P8"))


patient_microbe_abundance$variable<-factor(patient_microbe_abundance$variable,levels=rev(unique(patient_microbe_abundance$variable)))
####绘制堆积柱状图
library(ggplot2)
library(ggprism)
library(reshape)
library(ggalluvial)
library(randomcoloR)
palette <- distinctColorPalette(14)

#[1] "#7BABD6" "#954DE9" "#A782D8" "#C9E490" "#6DDCA8" "#D74DCB" "#90DFDD" "#DBBB56" "#D9B6D0" "#D3EB53" "#DA8A68" "#D3D3B5"
#[13] "#68E563" "#E46792"

p1 <- ggplot(data=patient_microbe_abundance,aes(variable,value,fill=patient))+
  geom_bar(stat="identity", position="fill",color="black", width=0.6,size=0.4)+
  scale_fill_manual(values=palette)+
  theme(
    axis.title=element_text(size=15,face="plain",color="black"),
    axis.text = element_text(size=15,face="plain",color="black"),
    legend.title=element_text(size=15,face="plain",color="black"),
    legend.position = "right",
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black", size = 0.4))+theme_bw()+
  theme(text=element_text(family="B",size=20))+
  theme(axis.ticks.length=unit(-0.25, "cm"), 
        axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
        axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")) )
p1
windowsFonts(A=windowsFont("Times New Roman"),
             B=windowsFont("Arial"))

library(ggalluvial)

p2 <- ggplot(data=patient_microbe_abundance,aes(variable,value,fill=patient,stratum = patient, alluvium = patient)) +
  geom_stratum(color="black",width=0.6,size=0.5)+
  geom_flow(alpha = 0.5) +  #绘制同类别之间的连接线
  scale_fill_manual(values=palette) +
  scale_y_continuous(name = "Relative abundance (%)")+
  theme(
    axis.title=element_text(size=15,face="Arial",color="black"),
    axis.text = element_text(size=15,face="Arial",color="black"),
    legend.title=element_text(size=15,face="Arial",color="black"),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black", size = 0.5))+theme_bw()+
  theme(text=element_text(family="B",size=20))+
  theme(axis.ticks.length=unit(-0.25, "cm"),
        axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")),
        axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")) )+
  guides(fill=guide_legend(keywidth = 1, keyheight = 1)) +
  theme_prism(base_fontface = "plain", 
              base_size = 16, 
              base_line_size = 0.8, 
              axis_text_angle = 0)+
  theme(legend.position = 'right')+coord_flip()
p2

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/INVADE-seq/微生物丰度热图/热图旁边堆积柱状图.jpg",
       p2, width = 12,height = 15)
ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/INVADE-seq/微生物丰度热图/热图旁边堆积柱状图.pdf",
       p2, width = 12,height = 15)

#拼图
library(patchwork)
p1 + p2 + plot_layout(ncol = 2, widths = c(2, 2))




####每个病人热图----

library(tidyverse)
library(ggnewscale)
library(MetBrewer)
library(patchwork)
library(ggtext)

patient_microbe<-read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/INVADE-seq/微生物丰度热图/每个病人物种丰度.csv")
rownames(patient_microbe)<- paste(patient_microbe$type,patient_microbe$Group.2,sep="-")
patient_microbe<- patient_microbe[,-c(1,2)]
relative_abundance<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/INVADE-seq/物种occurance frequency/relative_abundance_sorted_sum只根据sum降序丰度.csv")
colnames(relative_abundance)[1]
genus_top30<- relative_abundance$X[1:30]

patient_microbe_abundance<-patient_microbe[,which(colnames(patient_microbe)%in%genus_top30)]
patient_microbe_abundance<-patient_microbe_abundance[,match(genus_top30,colnames(patient_microbe_abundance))]
View(patient_microbe_abundance)

####添加分组信息
patient_microbe_abundance$group <- rownames(patient_microbe_abundance)

df <- patient_microbe_abundance %>% 
  pivot_longer(-group)

df$group <- factor(df$group,levels=c("Pre_NMPR-P3","Pre_NMPR-P5",
                                                                                       "Pre_NMPR-P6","Pre_NMPR-P10",
                                                                                       "Pre_MPR-P1","Pre_MPR-P2",
                                                                                       "Pre_MPR-P4","Pre_MPR-P7",
                                                                                       "Pre_MPR-P8","Pre_MPR-P9",
                                                                                       "Post_MPR-P1","Post_MPR-P4",
                                                                                       "Post_MPR-P7","Post_MPR-P8"))

df$name<- factor(df$name,levels=rev(c("Massilia", "Herbaspirillum","Corynebacterium","Sphingomonas",               
                                        "Treponema","Acinetobacter","Achromobacter",              
                                        "Fusobacterium","Arthrobacter","Rothia",                     
                                        "Lawsonella","Rhodococcus","Blattabacterium",            
                                        "Methylobacterium","Duganella" ,"Dietzia",                    
                                        "unclassified.Prevotellaceae","Streptococcus","Prevotella" ,                
                                        "Aerococcus","Planococcus","Staphylococcus",             
                                        "Campylobacter","Capnocytophaga","Bradyrhizobium",            
                                        "Propionibacterium","Lactobacillus","Actinomyces",                
                                        "Blastococcus","Belnapia")))

df$value<- scale(df$value)
df5 <-patient_microbe_abundance%>% select(31,where(is.numeric)) %>% 
  pivot_longer(-group)  

cell_type_color <- c("#FF6347","#ffc556","#1E90FF","#0000FF","#9370DB","#00FFFF","#32cd32","#006400","#FF00ff","#20B2AA")
response_color <-c("#4974a4","#4dae47","#f29600")


p1 <- df %>% ggplot(aes(x=sample,y=name))+
  geom_tile(data=df2,
            color="black",aes(fill=cell_type_new))+
  scale_fill_manual(values=cell_type_color)+
  new_scale_fill()+
  geom_tile(data=df1,color="black",aes(fill=MPR_Response2))+
  scale_fill_manual(values=response_color)+
  new_scale_fill()+
  labs(x=NULL,y=NULL)+
  theme(axis.text.x = element_blank(),
        axis.text.y=element_text(color="black"),
        #legend.position = c(1.6,-2),
        legend.location = "left",
        legend.direction = "vertical",
        legend.spacing.x =unit(0.1,"cm"),
        legend.spacing.y = unit(0.1,"cm"),
        legend.key.height = unit(0.5,"cm"),
        legend.key.width = unit(0.5,"cm"),
        plot.margin=unit(c(0,1,-20,0.5),unit="cm"),
        plot.background = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank(),
        legend.text = element_text(color="black"))

p1
col = rainbow(4)
names(col) <- letters[1:4]
#met.brewer("Cassatt1")
##00Bfff
#colorRampPalette(colors = c("blue","white","red"))(100)
p2 <- df %>% ggplot(aes(x=name,y=group))+
  geom_tile(data=df,color="black",
            aes(fill=value))+
  scale_fill_gradientn(colors = c('#00Bfff', "yellow", 'red'),na.value = NA)+
  new_scale_fill()+
  labs(x= NULL,y=NULL)+
  theme(axis.text.x = element_blank(),
        axis.text.y=element_text(color="black"),
        #legend.position = c(1.6,-2),
        legend.location = "left",
        legend.direction = "vertical",
        plot.margin=unit(c(-20,1,0.5,0.5),unit="cm"),
        plot.background = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank(),
        legend.title =element_markdown(),
        legend.text = element_text(color="black"))

p2


plot <- p1/p2+plot_layout(heights = c(0.1,2))

plot
ggsave(plot,file="./heatmap.pdf",dpi=300,width =9,height=10)


####画热图旁边的柱状图#####
library(tidyverse)

occurrence_freq<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/INVADE-seq/微生物丰度热图/relative_abundance_sortedsum和occur一起降序.csv") %>%
  select(X,Post_MPR_occurrence_frequency,Pre_MPR_occurrence_frequency,Pre_NMPR_occurrence_frequency) %>% slice_head(n=30)

colnames(occurrence_freq)<-c("genus","Post_MPR","Pre_MPR","Pre_NMPR")
View(occurrence_freq)

occurrence_freq_long<- occurrence_freq %>% pivot_longer(-genus)

library(ggtext)
library(ggpubr)
library(rstatix)

occurrence_freq_long$name<-factor(occurrence_freq_long$name,levels=c("Pre_NMPR","Pre_MPR","Post_MPR"))

occurrence_freq_long$genus<-factor(occurrence_freq_long$genus,levels=rev(c("Massilia", "Herbaspirillum","Corynebacterium","Sphingomonas",               
                                                                           "Treponema","Acinetobacter","Achromobacter",              
                                                                           "Fusobacterium","Arthrobacter","Rothia",                     
                                                                           "Lawsonella","Rhodococcus","Blattabacterium",            
                                                                           "Methylobacterium","Duganella" ,"Dietzia",                    
                                                                           "unclassified.Prevotellaceae","Streptococcus","Prevotella" ,                
                                                                           "Aerococcus","Planococcus","Staphylococcus",             
                                                                           "Campylobacter","Capnocytophaga","Bradyrhizobium",            
                                                                           "Propionibacterium","Lactobacillus","Actinomyces",                
                                                                           "Blastococcus","Belnapia")))

occurrence_freq_plot<- ggplot(data = occurrence_freq_long,mapping= aes(x = genus, y = value, fill = name))+
  geom_bar( stat="identity" ,position=position_dodge(1))+theme_test()+
  geom_text(aes(label =paste(round(value,0),"%",sep = ""), x = genus, y =value),hjust =0.7,vjust=0.2,angle=0, size = 5) +
  scale_fill_brewer(palette = 'Greens')+
  labs(x="",y = " ")+
  theme(strip.text.x =element_text(size = 10, colour = "black"),
        strip.text.y = element_text(size = 10, colour = "black"),
        strip.background.x = element_rect(fill = "white", colour = "black",linewidth = 0.8),
        strip.background.y = element_rect(fill = "white", colour = "black",linewidth = 0.8),
        axis.title.y = element_text(size=17, color="black", face="bold"),
        axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_text(size=12, color="black", face="bold"),
        legend.position = "none")+
  coord_flip() + #绘制水平柱状图
  facet_wrap(.~name)


ggsave(occurrence_freq_plot,file="./热图旁边的occurrence_freq.jpg",dpi=300,width =7,height=10)
ggsave(occurrence_freq_plot,file="./热图旁边的occurrence_freq.pdf",dpi=300,width =7,height=10)


#########每个分组的不同细胞类型top 的属##----
##提取每个分组不同细胞类型top10的属


microbe_group <- cbind(MPRNMPR_object_miMPRobe@meta.data$cell_id,
                       MPRNMPR_object_miMPRobe@meta.data$orig.ident,
                       MPRNMPR_object_miMPRobe@meta.data$cell_type_new,
                       MPRNMPR_object_miMPRobe@meta.data$group_microbe2,
                       MPRNMPR_object_miMPRobe@meta.data$MPR_Response2,
                       MPRNMPR_object_miMPRobe@meta.data$CRNCR_Response2,
                       MPRNMPR_object_miMPRobe@meta.data$Reads_counts,
                       MPRNMPR_object_miMPRobe@meta.data$UMIs_counts,
                       MPRNMPR_object_miMPRobe@meta.data$Patients) %>%
  as.data.frame()

colnames(microbe_group)<-c("cell_id","orig.ident","cell_type_new","group_microbe2",
                           "MPR_Response2","CRNCR_Response2","Reads_counts","UMIs_counts","Patients")

umi_16s<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/INVADE-seq/16s_filter.combined.genus.umi.matrix.csv",header = T,row.names =1)
umi_16s_MPR<-umi_16s[match(microbe_group$cell_id,rownames(umi_16s)),] %>% drop_na()
microbe_group_heatmap<- microbe_group[match(rownames(umi_16s_MPR),microbe_group$cell_id),]
names(microbe_group_heatmap)

###计算相对丰度
umi_16s_MPR_t<- t(umi_16s_MPR) %>% as.data.frame()
relative_abundance <- umi_16s_MPR_t %>%
  mutate(across(everything(),~ round(.x / sum(.x), 10)))

# relative_abundance$sum <- rowSums(relative_abundance)
# relative_abundance<- relative_abundance[order(relative_abundance$sum,decreasing = TRUE),]
# genus_rank<- as.data.frame(cbind(rownames(relative_abundance),relative_abundance$sum))
# View(genus_rank)
# colnames(genus_rank)<-c("Genus","sum")


relative_abundance_t<-t(relative_abundance) %>% as.data.frame()
View(relative_abundance_t)
colnames(microbe_group_heatmap)
relative_abundance_t$cell_type_new<-microbe_group_heatmap[match(rownames(relative_abundance_t),microbe_group_heatmap$cell_id),3]

relative_abundance_t$MPR_Response2<-microbe_group_heatmap[match(rownames(relative_abundance_t),microbe_group_heatmap$cell_id),5]

relative_abundance_t$Patients<-microbe_group_heatmap[match(rownames(relative_abundance_t),microbe_group_heatmap$cell_id),9]

genus_list<- colnames(relative_abundance_t)[c(-514,-515)]
result <- data.frame()
View(result)

for (i in 1:length(genus_list)){
  col_data<- aggregate(relative_abundance_t[genus_list[i]], by=list(type=relative_abundance_t$cell_type_new,relative_abundance_t$Patients),sum)
  col_name <- paste0(genus_list[i])
  if(i == 1) {
    result <- data.frame(col_data)
  } else {
    result <- cbind(result, col_data[[genus_list[i]]])
  }
  colnames(result)[i+2] <- col_name
}

View(result)
colnames(result)
write.csv(result,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/INVADE-seq/不同分组不同细胞类型丰度条形图/不同分组不同细胞类型微生物丰度数据.csv")
write.csv(result,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/INVADE-seq/不同分组不同细胞类型丰度条形图/不同病人不同细胞类型微生物丰度数据.csv")

result<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/INVADE-seq/不同分组不同细胞类型丰度条形图/不同分组不同细胞类型微生物丰度数据.csv",
                  row.names = 1,header = T)

###根据数值获得top10
df <- result %>%
  select(Group.2,type,where(is.numeric)) %>%
  pivot_longer(-c(Group.2,type)) %>%
  arrange(Group.2,type, desc(value)) %>%
  group_by(Group.2,type) %>%
  top_n(10, value)

###根据数值获得前10行
df <- result %>%
  select(Group.2,type,where(is.numeric)) %>%
  pivot_longer(-c(Group.2,type)) %>%
  arrange(Group.2,type, desc(value)) %>%
  group_by(Group.2,type) %>%
  slice_head(n = 10)
View(df)

write.csv(df,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/INVADE-seq/不同分组不同细胞类型丰度条形图/不同分组不同细胞类型微生物丰度top数据.csv")

write.csv(df,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/INVADE-seq/不同分组不同细胞类型丰度条形图/不同病人不同细胞类型微生物丰度top数据.csv")

df2 <- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/INVADE-seq/不同分组不同细胞类型丰度条形图/不同分组不同细胞类型微生物丰度top数据.csv")

df2 <- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/INVADE-seq/不同分组不同细胞类型丰度条形图/不同病人不同细胞类型微生物丰度top数据.csv")

library(randomcoloR)
color <- distinctColorPalette(56)
df2$Group.2<-factor(df2$Group.2,levels = c("Pre_NMPR","Pre_MPR","Post_MPR"))
df2$Group.2<-factor(df2$Group.2,levels =c("P1","P2","P3","P4","P5","P6","P7","P8","P9","P10") )
df2$type <- factor(df2$type,levels = c('T cells', 'B cells','Plasma cells', 'Myeloid cells',  
                                       'Mast cells',"Endothelial cells","Fibroblasts","Smooth muscle cells",
                                       "Pericytes","Epithelial cells")) %>% rev()
df2$name <- factor(df2$name,levels = df2$name %>% unique() %>% rev())

library(tidyverse)
library(RColorBrewer)
library(MetBrewer)
library(ggnewscale)

# df <- read_tsv("DEG.Ko.enrich.xls") %>% 
#   select(3,5,6,8,9) %>% 
#   mutate(count=Input.number/Background.number) %>% 
#   arrange(desc(count)) %>% head(20)
# 
# View(df)

####循环画图不同分组不同细胞微生物丰度##----

levels = c('T cells', 'B cells','Plasma cells', 'Myeloid cells',  
           'Mast cells',"Endothelial cells","Fibroblasts","Smooth muscle cells",
           "Pericytes","Epithelial cells")

#设置物种对应色板
genus_56<- unique(df2$name) %>% as.data.frame()
genus_color<-distinctColorPalette(56) 
genus_color <-c("#B29EE5","#DECEDB","#63EBC7","#E079AA","#A49FBC","#A17F65", "#9BB479", "#64B496", "#E9AE74", "#4E5EE4", "#DCE4BB",
"#55D1DE", "#ECACE2", "#5DEC90", "#B843ED", "#AACDE7", "#DD843C", "#C395B5", "#EE523D", "#A24796", "#E2BC3D", "#667CAE",
"#CEEDAA", "#E16ADE", "#846579", "#DB35D7", "#95E8A2", "#5D39EF", "#7EAAEA", "#8867AD", "#B0EDD0", "#9550CF", "#4EB24F",
"#7AEC49", "#E94EB1", "#E08DE2" ,"#CAE5EE", "#E7C5AA", "#D6E841", "#587170", "#EEB3C2", "#B27DE7", "#E9F0E4", "#5C83E5",
"#E35178", "#60C0EB", "#AEE971", "#DE8781", "#99EBEA", "#E8D78F", "#B4BDB4", "#938839", "#719DB7", "#E2EA83", "#DCC7F0",
 "#A0BC4D")
genus_color_df <- cbind(genus_56,genus_color) %>% as.data.frame()
colnames(genus_color_df)<- c("genus","color")
View(genus_color_df)



cell_type_list<-c('T cells', 'B cells','Plasma cells', 'Myeloid cells',  
                  'Mast cells',"Endothelial cells","Fibroblasts","Smooth muscle cells",
                  "Pericytes","Epithelial cells")

plots=list()

cell<-'T cells'
df3<- df2 %>% filter(grepl("Pre_NMPR","Pre_MPR","Post_MPR", Group.2) | grepl("B cells", type)) %>% filter_if(is.numeric, all_vars((.) != 0 ))
cell_color <-genus_color_df %>% filter(genus %in% df3$name)
cell_color <- cell_color[match(df3$name,cell_color$genus),]
cell_color$color
cell_color$genus


for (cell in cell_type_list){
  df3<- df2 %>% filter(grepl("Pre_NMPR","Pre_MPR","Post_MPR", Group.2) | grepl(!!cell, type)) %>% filter_if(is.numeric, all_vars((.) != 0 ))
  cell_color <-cell_color <- eval(parse(text = paste0("c(", paste0("'", genus_color_df$genus, "'", " = '", genus_color_df$color, "'", collapse = ", "), ")")))
  df3$Group.2<- factor(df3$Group.2,levels=c("Pre_NMPR","Pre_MPR","Post_MPR"))
  df3$name <- factor(df3$name,levels = unique(df3$name))
  list1<- c("T cells","Endothelial cells",'Myeloid cells') ##75
  list2<- c('Plasma cells','Mast cells',"Smooth muscle cells","B cells") ##15
  list3<- c("Fibroblasts","Epithelial cells") ##130
  if (cell %in% list1){plots[[cell]]=df3 %>% 
    ggplot(aes(type,name))+
    geom_segment(aes(x=0,xend=value,y=name,yend=name,color=name),
                 size=2,show.legend = F,alpha=1)+
    geom_text(aes(label =name, x = value, y =name),hjust =-0.2,vjust=0.5,angle=0, size = 4) + 
    geom_point(aes(x = value, y = name,size=1,color=name,fill=name),shape=19)+
    scale_fill_manual(values=cell_color)+
    scale_color_manual(values=cell_color)+
    facet_grid(rows = vars(Group.2),cols = vars(type),scale="free_y")+
    coord_cartesian(xlim = c(0, 75))+
    geom_hline(yintercept=3,linetype='dashed', size =0.8)+
    labs(x="Microbial_relative_abundance",y=NULL)+
    theme(panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.text.x =element_text(size = 12, colour = "black"),
          strip.text.y = element_text(size = 12, colour = "black"),
          strip.background.x = element_rect(fill = "white", colour = "black",linewidth = 1.5),
          strip.background.y = element_rect(fill = "white", colour = "black",linewidth = 1.5),
          axis.line = element_line(size = 0.7),
          axis.ticks = element_line(size = 0.7),
          axis.ticks.length = unit(0.1, "cm"),
          axis.text.x = element_text(size = 15, colour = "black"),
          axis.text.y = element_blank(),
          axis.title.y = element_blank(), 
          axis.title.x = element_text(size = 14, colour = "black"),  
          axis.ticks.y = element_blank(), 
          legend.key=element_blank(),
          legend.title = element_text(color="black",size=10), 
          legend.text = element_text(color="black",size=10), 
          legend.key.size=unit(0.5,'cm'),
          legend.background=element_blank(),
          legend.position ="none")+
    guides(color = guide_legend(override.aes=list(shape = 15)))
  }else if(cell %in% list2){plots[[cell]]=df3 %>% 
    ggplot(aes(type,name))+
    geom_segment(aes(x=0,xend=value,y=name,yend=name,color=name),
                 size=2,show.legend = F,alpha=1)+
    geom_text(aes(label =name, x = value, y =name),hjust =-0.2,vjust=0.5,angle=0, size = 4) + 
    geom_point(aes(x = value, y = name,size=1,color=name,fill=name),shape=19)+
    scale_fill_manual(values=cell_color)+
    scale_color_manual(values=cell_color)+
    facet_grid(rows = vars(Group.2),cols = vars(type),scale="free_y")+
    coord_cartesian(xlim = c(0, 15))+
    geom_hline(yintercept=3,linetype='dashed', size =0.8)+
    labs(x="Microbial_relative_abundance",y=NULL)+
    theme(panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.text.x =element_text(size = 12, colour = "black"),
          strip.text.y = element_text(size = 12, colour = "black"),
          strip.background.x = element_rect(fill = "white", colour = "black",linewidth = 1.5),
          strip.background.y = element_rect(fill = "white", colour = "black",linewidth = 1.5),
          axis.line = element_line(size = 0.7),
          axis.ticks = element_line(size = 0.7),
          axis.ticks.length = unit(0.1, "cm"),
          axis.text.x = element_text(size = 15, colour = "black"),
          axis.text.y = element_blank(),
          axis.title.y = element_blank(), 
          axis.title.x = element_text(size = 14, colour = "black"),  
          axis.ticks.y = element_blank(), 
          legend.key=element_blank(),
          legend.title = element_text(color="black",size=10), 
          legend.text = element_text(color="black",size=10), 
          legend.key.size=unit(0.5,'cm'),
          legend.background=element_blank(),
          legend.position ="none")+
    guides(color = guide_legend(override.aes=list(shape = 15)))}
  else if(cell %in% list3){plots[[cell]]=df3 %>% 
      ggplot(aes(type,name))+
      geom_segment(aes(x=0,xend=value,y=name,yend=name,color=name),
                   size=2,show.legend = F,alpha=1)+
      geom_text(aes(label =name, x = value, y =name),hjust =-0.2,vjust=0.5,angle=0, size = 4) + 
      geom_point(aes(x = value, y = name,size=1,color=name,fill=name),shape=19)+
      scale_fill_manual(values=cell_color)+
      scale_color_manual(values=cell_color)+
      facet_grid(rows = vars(Group.2),cols = vars(type),scale="free_y")+
      coord_cartesian(xlim = c(0, 130))+
      geom_hline(yintercept=3,linetype='dashed', size =0.8)+
      labs(x="Microbial_relative_abundance",y=NULL)+
      theme(panel.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            strip.text.x =element_text(size = 12, colour = "black"),
            strip.text.y = element_text(size = 12, colour = "black"),
            strip.background.x = element_rect(fill = "white", colour = "black",linewidth = 1.5),
            strip.background.y = element_rect(fill = "white", colour = "black",linewidth = 1.5),
            axis.line = element_line(size = 0.7),
            axis.ticks = element_line(size = 0.7),
            axis.ticks.length = unit(0.1, "cm"),
            axis.text.x = element_text(size = 15, colour = "black"),
            axis.text.y = element_blank(),
            axis.title.y = element_blank(), 
            axis.title.x = element_text(size = 14, colour = "black"),  
            axis.ticks.y = element_blank(), 
            legend.key=element_blank(),
            legend.title = element_text(color="black",size=10), 
            legend.text = element_text(color="black",size=10), 
            legend.key.size=unit(0.5,'cm'),
            legend.background=element_blank(),
            legend.position ="none")+
      guides(color = guide_legend(override.aes=list(shape = 15)))}
  else{plots[[cell]]=df3 %>% 
    ggplot(aes(type,name))+
    geom_segment(aes(x=0,xend=value,y=name,yend=name,color=name),
                 size=2,show.legend = F,alpha=1)+
    geom_text(aes(label =name, x = value, y =name),hjust =-0.2,vjust=0.5,angle=0, size = 4) + 
    geom_point(aes(x = value, y = name,size=1,color=name,fill=name),shape=19)+
    scale_fill_manual(values=cell_color)+
    scale_color_manual(values=cell_color)+
    facet_grid(rows = vars(Group.2),cols = vars(type),scale="free_y")+
    coord_cartesian(xlim = c(0, 20))+
    geom_hline(yintercept=3,linetype='dashed', size =0.8)+
    labs(x="Microbial_relative_abundance",y=NULL)+
    theme(panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.text.x =element_text(size = 12, colour = "black"),
          strip.text.y = element_text(size = 12, colour = "black"),
          strip.background.x = element_rect(fill = "white", colour = "black",linewidth = 1.5),
          strip.background.y = element_rect(fill = "white", colour = "black",linewidth = 1.5),
          axis.line = element_line(size = 0.7),
          axis.ticks = element_line(size = 0.7),
          axis.ticks.length = unit(0.1, "cm"),
          axis.text.x = element_text(size = 15, colour = "black"),
          axis.text.y = element_blank(),
          axis.title.y = element_blank(), 
          axis.title.x = element_text(size = 14, colour = "black"),  
          axis.ticks.y = element_blank(), 
          legend.key=element_blank(),
          legend.title = element_text(color="black",size=10), 
          legend.text = element_text(color="black",size=10), 
          legend.key.size=unit(0.5,'cm'),
          legend.background=element_blank(),
          legend.position ="none")+
    guides(color = guide_legend(override.aes=list(shape = 15)))}
  
}



library(cowplot)
library(patchwork)

plots_<-plot_grid(plots[['T cells']],
                  plots[['B cells']],
                  plots[['Plasma cells']],
                  plots[['Myeloid cells']],
                  plots[['Mast cells']],
                  plots[['Endothelial cells']],
                  plots[['Fibroblasts']],
                  plots[['Smooth muscle cells']],
                  plots[['Pericytes']],
                  plots[['Epithelial cells']],
                  align = "h",  #axis = 'l',
                  nrow =3,
                  ncol=3)
setwd("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/INVADE-seq/不同分组不同细胞类型丰度条形图/")
ggsave(plots_,file="MPR_microbial_relative_abundance_top10.jpg",width =12,height = 20)
ggsave(plots_,file="MPR_microbial_relative_abundance_top10.pdf",width =12,height = 20)



#######top5不同分组不同细胞类型丰度----
{df <- result %>%
  select(Group.2,type,where(is.numeric)) %>%
  pivot_longer(-c(Group.2,type)) %>%
  arrange(Group.2,type, desc(value)) %>%
  group_by(Group.2,type) %>%
  slice_head(n = 5)
View(df)

write.csv(df,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/INVADE-seq/不同分组不同细胞类型丰度条形图/不同分组不同细胞类型微生物丰度top_5数据.csv")

df_top5 <- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/INVADE-seq/不同分组不同细胞类型丰度条形图/不同分组不同细胞类型微生物丰度top_5数据.csv")

genus_color <-c("#B29EE5","#DECEDB","#63EBC7","#E079AA","#A49FBC","#A17F65", "#9BB479", "#64B496", "#E9AE74", "#4E5EE4", "#DCE4BB",
                "#55D1DE", "#ECACE2", "#5DEC90", "#B843ED", "#AACDE7", "#DD843C", "#C395B5", "#EE523D", "#A24796", "#E2BC3D", "#667CAE",
                "#CEEDAA", "#E16ADE", "#846579", "#DB35D7", "#95E8A2", "#5D39EF", "#7EAAEA", "#8867AD", "#B0EDD0", "#9550CF", "#4EB24F",
                "#7AEC49", "#E94EB1", "#E08DE2" ,"#CAE5EE", "#E7C5AA", "#D6E841", "#587170", "#EEB3C2", "#B27DE7", "#E9F0E4", "#5C83E5",
                "#E35178", "#60C0EB", "#AEE971", "#DE8781", "#99EBEA", "#E8D78F", "#B4BDB4", "#938839", "#719DB7", "#E2EA83", "#DCC7F0",
                "#A0BC4D")
genus_color_df <- cbind(genus_56,genus_color) %>% as.data.frame()
colnames(genus_color_df)<- c("genus","color")
View(genus_color_df)

cell_type_list<-c('T cells', 'B cells','Plasma cells', 'Myeloid cells',  
                  'Mast cells',"Endothelial cells","Fibroblasts","Smooth muscle cells",
                  "Pericytes","Epithelial cells")

plots=list()
for (cell in cell_type_list){
  df3<- df_top5 %>% filter(grepl("Pre_NMPR","Pre_MPR","Post_MPR", Group.2) | grepl(!!cell, type)) %>% filter_if(is.numeric, all_vars((.) != 0 ))
  cell_color <-cell_color <- eval(parse(text = paste0("c(", paste0("'", genus_color_df$genus, "'", " = '", genus_color_df$color, "'", collapse = ", "), ")")))
  df3$Group.2<- factor(df3$Group.2,levels=c("Pre_NMPR","Pre_MPR","Post_MPR"))
  df3$name <- factor(df3$name,levels = unique(df3$name))
  list1<- c("T cells","Endothelial cells",'Myeloid cells') ##75
  list2<- c('Plasma cells','Mast cells',"Smooth muscle cells","B cells") ##15
  list3<- c("Fibroblasts","Epithelial cells") ##130
  if (cell %in% list1){plots[[cell]]=df3 %>% 
    ggplot(aes(type,name))+
    geom_segment(aes(x=0,xend=value,y=name,yend=name,color=name),
                 size=2,show.legend = F,alpha=1)+
    geom_text(aes(label =name, x = value, y =name),hjust =-0.2,vjust=0.5,angle=0, size = 4) + 
    geom_point(aes(x = value, y = name,size=1,color=name,fill=name),shape=19)+
    scale_fill_manual(values=cell_color)+
    scale_color_manual(values=cell_color)+
    facet_grid(rows = vars(Group.2),cols = vars(type),scale="free_y")+
    coord_cartesian(xlim = c(0, 75))+
    geom_hline(yintercept=3,linetype='dashed', size =0.8)+
    labs(x="Microbial_relative_abundance",y=NULL)+
    theme(panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.text.x =element_text(size = 12, colour = "black"),
          strip.text.y = element_text(size = 12, colour = "black"),
          strip.background.x = element_rect(fill = "white", colour = "black",linewidth = 1.5),
          strip.background.y = element_rect(fill = "white", colour = "black",linewidth = 1.5),
          axis.line = element_line(size = 0.7),
          axis.ticks = element_line(size = 0.7),
          axis.ticks.length = unit(0.1, "cm"),
          axis.text.x = element_text(size = 15, colour = "black"),
          axis.text.y = element_blank(),
          axis.title.y = element_blank(), 
          axis.title.x = element_text(size = 14, colour = "black"),  
          axis.ticks.y = element_blank(), 
          legend.key=element_blank(),
          legend.title = element_text(color="black",size=10), 
          legend.text = element_text(color="black",size=10), 
          legend.key.size=unit(0.5,'cm'),
          legend.background=element_blank(),
          legend.position ="none")+
    guides(color = guide_legend(override.aes=list(shape = 15)))
  }else if(cell %in% list2){plots[[cell]]=df3 %>% 
    ggplot(aes(type,name))+
    geom_segment(aes(x=0,xend=value,y=name,yend=name,color=name),
                 size=2,show.legend = F,alpha=1)+
    geom_text(aes(label =name, x = value, y =name),hjust =-0.2,vjust=0.5,angle=0, size = 4) + 
    geom_point(aes(x = value, y = name,size=1,color=name,fill=name),shape=19)+
    scale_fill_manual(values=cell_color)+
    scale_color_manual(values=cell_color)+
    facet_grid(rows = vars(Group.2),cols = vars(type),scale="free_y")+
    coord_cartesian(xlim = c(0, 15))+
    geom_hline(yintercept=3,linetype='dashed', size =0.8)+
    labs(x="Microbial_relative_abundance",y=NULL)+
    theme(panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.text.x =element_text(size = 12, colour = "black"),
          strip.text.y = element_text(size = 12, colour = "black"),
          strip.background.x = element_rect(fill = "white", colour = "black",linewidth = 1.5),
          strip.background.y = element_rect(fill = "white", colour = "black",linewidth = 1.5),
          axis.line = element_line(size = 0.7),
          axis.ticks = element_line(size = 0.7),
          axis.ticks.length = unit(0.1, "cm"),
          axis.text.x = element_text(size = 15, colour = "black"),
          axis.text.y = element_blank(),
          axis.title.y = element_blank(), 
          axis.title.x = element_text(size = 14, colour = "black"),  
          axis.ticks.y = element_blank(), 
          legend.key=element_blank(),
          legend.title = element_text(color="black",size=10), 
          legend.text = element_text(color="black",size=10), 
          legend.key.size=unit(0.5,'cm'),
          legend.background=element_blank(),
          legend.position ="none")+
    guides(color = guide_legend(override.aes=list(shape = 15)))}
  else if(cell %in% list3){plots[[cell]]=df3 %>% 
    ggplot(aes(type,name))+
    geom_segment(aes(x=0,xend=value,y=name,yend=name,color=name),
                 size=2,show.legend = F,alpha=1)+
    geom_text(aes(label =name, x = value, y =name),hjust =-0.2,vjust=0.5,angle=0, size = 4) + 
    geom_point(aes(x = value, y = name,size=1,color=name,fill=name),shape=19)+
    scale_fill_manual(values=cell_color)+
    scale_color_manual(values=cell_color)+
    facet_grid(rows = vars(Group.2),cols = vars(type),scale="free_y")+
    coord_cartesian(xlim = c(0, 130))+
    geom_hline(yintercept=3,linetype='dashed', size =0.8)+
    labs(x="Microbial_relative_abundance",y=NULL)+
    theme(panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.text.x =element_text(size = 12, colour = "black"),
          strip.text.y = element_text(size = 12, colour = "black"),
          strip.background.x = element_rect(fill = "white", colour = "black",linewidth = 1.5),
          strip.background.y = element_rect(fill = "white", colour = "black",linewidth = 1.5),
          axis.line = element_line(size = 0.7),
          axis.ticks = element_line(size = 0.7),
          axis.ticks.length = unit(0.1, "cm"),
          axis.text.x = element_text(size = 15, colour = "black"),
          axis.text.y = element_blank(),
          axis.title.y = element_blank(), 
          axis.title.x = element_text(size = 14, colour = "black"),  
          axis.ticks.y = element_blank(), 
          legend.key=element_blank(),
          legend.title = element_text(color="black",size=10), 
          legend.text = element_text(color="black",size=10), 
          legend.key.size=unit(0.5,'cm'),
          legend.background=element_blank(),
          legend.position ="none")+
    guides(color = guide_legend(override.aes=list(shape = 15)))}
  else{plots[[cell]]=df3 %>% 
    ggplot(aes(type,name))+
    geom_segment(aes(x=0,xend=value,y=name,yend=name,color=name),
                 size=2,show.legend = F,alpha=1)+
    geom_text(aes(label =name, x = value, y =name),hjust =-0.2,vjust=0.5,angle=0, size = 4) + 
    geom_point(aes(x = value, y = name,size=1,color=name,fill=name),shape=19)+
    scale_fill_manual(values=cell_color)+
    scale_color_manual(values=cell_color)+
    facet_grid(rows = vars(Group.2),cols = vars(type),scale="free_y")+
    coord_cartesian(xlim = c(0, 20))+
    geom_hline(yintercept=3,linetype='dashed', size =0.8)+
    labs(x="Microbial_relative_abundance",y=NULL)+
    theme(panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.text.x =element_text(size = 12, colour = "black"),
          strip.text.y = element_text(size = 12, colour = "black"),
          strip.background.x = element_rect(fill = "white", colour = "black",linewidth = 1.5),
          strip.background.y = element_rect(fill = "white", colour = "black",linewidth = 1.5),
          axis.line = element_line(size = 0.7),
          axis.ticks = element_line(size = 0.7),
          axis.ticks.length = unit(0.1, "cm"),
          axis.text.x = element_text(size = 15, colour = "black"),
          axis.text.y = element_blank(),
          axis.title.y = element_blank(), 
          axis.title.x = element_text(size = 14, colour = "black"),  
          axis.ticks.y = element_blank(), 
          legend.key=element_blank(),
          legend.title = element_text(color="black",size=10), 
          legend.text = element_text(color="black",size=10), 
          legend.key.size=unit(0.5,'cm'),
          legend.background=element_blank(),
          legend.position ="none")+
    guides(color = guide_legend(override.aes=list(shape = 15)))}
  
}



library(cowplot)
library(patchwork)

plots_<-plot_grid(plots[['T cells']],
                  plots[['B cells']],
                  plots[['Plasma cells']],
                  plots[['Myeloid cells']],
                  plots[['Mast cells']],
                  plots[['Endothelial cells']],
                  plots[['Fibroblasts']],
                  plots[['Smooth muscle cells']],
                  plots[['Pericytes']],
                  plots[['Epithelial cells']],
                  align = "h",  #axis = 'l',
                  nrow =3,
                  ncol=3)
setwd("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/INVADE-seq/不同分组不同细胞类型丰度条形图/")
ggsave(plots_,file="MPR_microbial_relative_abundance_top5.jpg",width =15,height = 15)
ggsave(plots_,file="MPR_microbial_relative_abundance_top5.pdf",width =15,height = 15)
}




#tcell_color_group<-c("blue","red","green") 

# tell_pre_nmpr_plot<-tell_pre_nmpr %>% 
#   ggplot(aes(name,type))+
#   geom_segment(aes(x=name,xend=name,y=0,yend=value,color=name),
#                size=5,show.legend = F,alpha=1)+
#   geom_text(aes(label =name, x = name, y =value),hjust =-0.1,vjust=0.5,angle=90, size = 4) + 
#   geom_point(aes(x = name, y = value,size=value,color=name,fill=name),shape=19)+
#   scale_fill_manual(values=tcell_color)+
#   scale_color_manual(values=tcell_color)+
#   facet_grid(rows = vars(type),cols = vars(Group.2),scale="free_y")+
#   coord_cartesian(ylim = c(0, 20))+
#   geom_vline(xintercept=8.6,linetype='dashed', linewidth =0.8)+
#   labs(x=NULL,y="Microbe_relative_abundance")+
#   theme(panel.background = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         strip.text.x =element_text(size = 10, colour = "black"),
#         strip.text.y = element_text(size = 12, colour = "black"),
#         strip.background.x = element_rect(fill = "white", colour = "black",linewidth = 2),
#         strip.background.y = element_rect(fill = "white", colour = "black",linewidth = 2),
#         axis.line = element_line(size = 0.7),
#         axis.ticks = element_line(size = 0.7),
#         axis.ticks.length = unit(0.1, "cm"),
#         axis.text.y = element_text(size = 15, colour = "black"),
#         axis.text.x = element_blank(),
#         axis.title.x = element_blank(), 
#         axis.title.y = element_text(size = 14, colour = "black"),  
#         axis.ticks.x = element_blank(), 
#         legend.key=element_blank(),
#         legend.title = element_text(color="black",size=10), 
#         legend.text = element_text(color="black",size=10), 
#         legend.key.size=unit(0.5,'cm'),
#         legend.background=element_blank(),
#         legend.position ="right")+
#   guides(color = guide_legend(
#     override.aes=list(shape = 15)))

# "#B29EE5"


####2bRAD与胞内菌交叉的属----

setwd("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/INVADE-seq/与2bRAD共有属")

invaseq_microbe<-read.csv("16s_filter.combined.genus.umi.matrix.csv",row.names = 1,header = T)
rad_microbe<- read_tsv("Abundance.filtered.anno.xls")
View(rad_microbe)
rad_genus<- rad_microbe$Genus
unique(rad_genus)
unique(colnames(invaseq_microbe))
invaseq_microbe_filter<- invaseq_microbe[,colnames(invaseq_microbe) %in% rad_genus]

unique(colnames(invaseq_microbe_filter))
write.csv(invaseq_microbe_filter,"16s_filter.combined.genus.umi.matrix_2brad共有过滤.csv")

######重新提取top30的物种----
microbe_group <- cbind(MPRNMPR_object_miMPRobe@meta.data$cell_id,
                       MPRNMPR_object_miMPRobe@meta.data$orig.ident,
                       MPRNMPR_object_miMPRobe@meta.data$cell_type_new,
                       MPRNMPR_object_miMPRobe@meta.data$group_microbe2,
                       MPRNMPR_object_miMPRobe@meta.data$MPR_Response2,
                       MPRNMPR_object_miMPRobe@meta.data$CRNCR_Response2,
                       MPRNMPR_object_miMPRobe@meta.data$Reads_counts,
                       MPRNMPR_object_miMPRobe@meta.data$UMIs_counts) %>%
  as.data.frame()

colnames(microbe_group)<-c("cell_id","orig.ident","cell_type_new","group_microbe2",
                           "MPR_Response2","CRNCR_Response2","Reads_counts","UMIs_counts")

umi_16s<- read.csv("16s_filter.combined.genus.umi.matrix_2brad共有过滤.csv",header = T,row.names =1)
umi_16s_MPR<-umi_16s[match(microbe_group$cell_id,rownames(umi_16s)),] %>% drop_na()
microbe_group_heatmap<- microbe_group[match(rownames(umi_16s_MPR),microbe_group$cell_id),]
names(microbe_group_heatmap)

umi_16s_MPR<-umi_16s[match(microbe_group$cell_id,rownames(umi_16s)),] %>% drop_na()
ncol(umi_16s_MPR)
nrow(umi_16s_MPR)
names(umi_16s)
microbe_group_heatmap<- microbe_group[match(rownames(umi_16s_MPR),microbe_group$cell_id),]
ncol(microbe_group_heatmap)
nrow(microbe_group_heatmap)
View(microbe_group_heatmap)

###计算相对丰度
umi_16s_MPR_t<- t(umi_16s_MPR) %>% as.data.frame()
relative_abundance <- umi_16s_MPR_t %>%
  mutate(across(everything(),~ round(.x / sum(.x), 10)))%>%select_if(~ all(!is.na(.)))

View(relative_abundance)
colSums(relative_abundance)

##整理数据计算Occurrence frequency
relative_abundance_t<-t(relative_abundance) %>% as.data.frame()
View(relative_abundance_t)

relative_abundance_t$Patients <-microbe_group_heatmap[match(rownames(relative_abundance_t),microbe_group_heatmap$cell_id),3]
relative_abundance_t$cell_type_new<-microbe_group_heatmap[match(rownames(relative_abundance_t),microbe_group_heatmap$cell_id),4]

relative_abundance_t$MPR_Response2<-microbe_group_heatmap[match(rownames(relative_abundance_t),microbe_group_heatmap$cell_id),6]

genus_list<- colnames(relative_abundance_t)[c(-177,-178,-179)]
result <- data.frame()
View(relative_abundance_t)

for (i in 1:length(genus_list)){
  col_data<- aggregate(relative_abundance_t[genus_list[i]], by=list(type=relative_abundance_t$MPR_Response2,relative_abundance_t$Patients),sum)
  col_name <- paste0(genus_list[i])
  if(i == 1) {
    result <- data.frame(col_data)
  } else {
    result <- cbind(result, col_data[[genus_list[i]]])
  }
  colnames(result)[i+2] <- col_name
}

View(result)

result_Occurrence <- result %>%
  group_by(type,Group.2) %>%
  summarise(across(everything(), ~sum(. > 0, na.rm = TRUE)))

View(result_Occurrence)
genus_list<- colnames(result_Occurrence)[c(-1,-2)]
result_Occurrence_frequency <- result_Occurrence %>%select(1,where(is.numeric)) %>% 
  group_by(type) %>%
  summarize_at(vars(genus_list), sum) %>% as.data.frame()

View(result_Occurrence_frequency)

divisor <- c(4, 6, 4)

rownames(result_Occurrence_frequency)<-c("Post_MPR_occurrence_frequency","Pre_MPR_occurrence_frequency","Pre_NMPR_occurrence_frequency")
result_Occurrence_frequency <- result_Occurrence_frequency[,-1]
# 使用 apply() 函数将每行的元素逐个除以相应的数字
result_Occurrence_frequency2 <- apply(t(result_Occurrence_frequency), 1, function(x) x / divisor)
View(result_Occurrence_frequency2)

result_Occurrence_frequency3 <- t(result_Occurrence_frequency2 *100) %>% as.data.frame()
View(result_Occurrence_frequency3)

###物种总数##

relative_abundance$sum <- rowSums(relative_abundance)
# 合并Occurrence
View(relative_abundance)

relative_abundance$Post_MPR_occurrence_frequency <- result_Occurrence_frequency3[match(rownames(relative_abundance),rownames(result_Occurrence_frequency3)),1]
relative_abundance$Pre_MPR_occurrence_frequency <- result_Occurrence_frequency3[match(rownames(relative_abundance),rownames(result_Occurrence_frequency3)),2]
relative_abundance$Pre_NMPR_occurrence_frequency <- result_Occurrence_frequency3[match(rownames(relative_abundance),rownames(result_Occurrence_frequency3)),3]

###根据多个分组降序
relative_abundance_sorted <- relative_abundance %>%
  arrange(desc(sum), desc(Post_MPR_occurrence_frequency), desc(Pre_MPR_occurrence_frequency),desc(Pre_NMPR_occurrence_frequency))

View(relative_abundance_sorted)

write.csv(relative_abundance_sorted,"relative_abundance_sorted.csv")

##只根据sum 降序
relative_abundance_sorted_sum<- relative_abundance[order(relative_abundance$sum,decreasing = TRUE),]
View(relative_abundance_sorted_sum)
write.csv(relative_abundance_sorted_sum,"relative_abundance_sorted_sum.csv")

genus_top30<- relative_abundance[1:30,-ncol(relative_abundance)]
View(genus_top30)


###提取分组信息##
#手动标准化

umi_16s_MPR_top<- umi_16s_MPR[match(rownames(genus_top30),colnames(umi_16s_MPR),)]
umi_16s_MPR_top$cell_type_new <- microbe_group_heatmap[match(rownames(umi_16s_MPR_top), microbe_group_heatmap$cell_id),3]
umi_16s_MPR_top$MPR_Response2 <- microbe_group_heatmap[match(rownames(umi_16s_MPR_top), microbe_group_heatmap$cell_id),5]

umi_16s_MPR_heatmap<-umi_16s_MPR_top %>% group_by(MPR_Response2,cell_type_new) %>% 
  summarise_all(sum, na.rm=TRUE)
View(umi_16s_MPR_heatmap)


normalize <- function(x) {
  return ((x - min(x)) / (max(x) - min(x)))
}
umi_16s_MPR_heatmap=as.data.frame(lapply(umi_16s_MPR_heatmap[,3:30],normalize))%>% cbind(umi_16s_MPR_heatmap[,c(1,2)])
View(umi_16s_MPR_heatmap)

umi_16s_MPR_heatmap<-scale(umi_16s_MPR_heatmap[,c(-1,-2)], center = TRUE, scale = TRUE)%>% cbind(umi_16s_MPR_heatmap[,c(1,2)])
View(umi_16s_MPR_heatmap)

#for (i in 1:nrow(umi_16s_MPR_heatmap)) umi_16s_MPR_heatmap[i, ] <- scale(log(unlist(umi_16s_MPR_heatmap[, ]) + 1, 2))

#for (i in 3:ncol(umi_16s_MPR_heatmap)) umi_16s_MPR_heatmap[,i] <- scale(log(unlist(umi_16s_MPR_heatmap[,i]) + 1, 2))

View(umi_16s_MPR_heatmap)

umi_16s_MPR_heatmap_scaled <- as.matrix(umi_16s_MPR_heatmap)
write.csv(umi_16s_MPR_heatmap,"./umi_16s_MPR_heatmap_top30_scaled.csv")

####各种亚型合并object----
library(Seurat)
library(tidyverse)
library(harmony)
library(dplyr)
library(patchwork)
library(SeuratWrappers)

MPRNMPR_object_miMPRobe_remove <- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/MPRNMPR_object_miMPRobe_remove.rds")
MPRNMPR_object_miMPRobe_remove = UpdateSeuratObject(MPRNMPR_object_miMPRobe_remove)
MPRNMPR_object_miMPRobe_remove
View(MPRNMPR_object_miMPRobe_remove)

####读取Tcellobject
T_cells_object <- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/T_cells_object.rds")
####B_cell_object
B_cells_object <- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/B_cells_object_combined使用.rds")
###髓系细胞
myeloid_cells_object <- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/myeloid_cells_object_filter使用.rds")
####成纤维细胞
CAFs_cells_object<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/CAFs_cells_object.rds")
####内皮细胞
edo_cells_object<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/edo_cells_object.rds")
###上皮细胞
epi_cells_object<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/epi_cells_object_nomerged.rds")

####Pericytes
Pericytes_object <- subset(MPRNMPR_object_miMPRobe_remove,idents = c("Pericytes"))
table(Pericytes_object$cell_type_new)

####Smooth muscle cells
SMC_object <- subset(MPRNMPR_object_miMPRobe_remove,idents = c("Smooth muscle cells"))
table(SMC_object$cell_type_new)

levels(Idents(MPRNMPR_object_miMPRobe_remove))

####提取cell_id origt.----

T_cells_id<- cbind(cell_id = T_cells_object@meta.data$cell_id) %>% 
  as.data.frame() %>% 
  cbind(orig.ident = T_cells_object@meta.data$orig.ident) %>% 
  cbind(cell_type =T_cells_object@meta.data$cd4_cd8_group)

B_cells_id<- cbind(cell_id = B_cells_object@meta.data$cell_id) %>% 
  as.data.frame() %>% 
  cbind(orig.ident = B_cells_object@meta.data$orig.ident) %>% 
  cbind(cell_type =B_cells_object@meta.data$B_cells_type)

myeloid_cells_id<- cbind(cell_id = myeloid_cells_object@meta.data$cell_id) %>% 
  as.data.frame() %>% 
  cbind(orig.ident = myeloid_cells_object@meta.data$orig.ident) %>% 
  cbind(cell_type =myeloid_cells_object@meta.data$myeloid_cells_type)

edo_cells_id<- cbind(cell_id = edo_cells_object@meta.data$cell_id) %>% 
  as.data.frame() %>% 
  cbind(orig.ident =edo_cells_object@meta.data$orig.ident) %>% 
  cbind(cell_type =edo_cells_object@meta.data$edo_cells_type)

epi_cells_id<- cbind(cell_id = epi_cells_object@meta.data$cell_id) %>% 
  as.data.frame() %>% 
  cbind(orig.ident =epi_cells_object@meta.data$orig.ident) %>% 
  cbind(cell_type =epi_cells_object@meta.data$epi_infercnv_group)
table(epi_cells_object@meta.data$epi_infercnv_group)

Pericytes_id<- cbind(cell_id = Pericytes_object@meta.data$cell_id) %>% 
  as.data.frame() %>% 
  cbind(orig.ident =Pericytes_object@meta.data$orig.ident) %>% 
  cbind(cell_type =Pericytes_object@meta.data$cell_type_new)

SMC_id<- cbind(cell_id = SMC_object@meta.data$cell_id) %>% 
  as.data.frame() %>% 
  cbind(orig.ident =SMC_object@meta.data$orig.ident) %>% 
  cbind(cell_type =SMC_object@meta.data$cell_type_new)

CAFs_cells_id<- cbind(cell_id = CAFs_cells_object@meta.data$cell_id) %>% 
  as.data.frame() %>% 
  cbind(orig.ident = CAFs_cells_object@meta.data$orig.ident) %>% 
  cbind(cell_type =CAFs_cells_object@meta.data$CAFs_cells_type)

merged_df <- rbind(T_cells_id, B_cells_id,myeloid_cells_id,
                   edo_cells_id,epi_cells_id,Pericytes_id,SMC_id,CAFs_cells_id)

unique(merged_df$cell_type)
cells_id<- cbind(cell_id = MPRNMPR_object_miMPRobe_remove@meta.data$cell_id) %>% 
  as.data.frame() %>% 
  cbind(orig.ident =MPRNMPR_object_miMPRobe_remove@meta.data$orig.ident) %>% 
  cbind(cell_type =MPRNMPR_object_miMPRobe_remove@meta.data$cell_type_new)

# cells_to_remove<-setdiff(cells_id$cell_id, merged_df$cell_id)
# 
# MPRNMPR_object_miMPRobe_remove <- subset(MPRNMPR_object_miMPRobe, cells = Cells(MPRNMPR_object_miMPRobe)[!Cells(MPRNMPR_object_miMPRobe) %in% cells_to_remove])
# 
metadata<- FetchData(MPRNMPR_object_miMPRobe_remove,"cell_id")
metadata$cell_id <- rownames(metadata)

View(metadata)

names(merged_df)<-c("cell_id","orig.ident","cell_subtype")
names(metadata)

metadata<- left_join(x=metadata,y=merged_df,by = join_by(cell_id==cell_id))

rownames(metadata)<-metadata$cell_id

MPRNMPR_object_miMPRobe_remove<- AddMetaData(MPRNMPR_object_miMPRobe_remove,metadata = metadata)

View(MPRNMPR_object_miMPRobe_remove)

table(rownames(MPRNMPR_object_miMPRobe_remove@meta.data) %in% merged_df$cell_id)
table(MPRNMPR_object_miMPRobe_remove@meta.data$cell_subtype)

saveRDS(MPRNMPR_object_miMPRobe_remove, file = "D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/MPRNMPR_object_miMPRobe_remove.rds")

MPRNMPR_object_miMPRobe_remove$cell_subtype
########################################
########################################

library(tidyverse)
MPRNMPR.umap.harmonyredu<- MPRNMPR_object_miMPRobe_remove@reductions$umap.harmony@cell.embeddings%>%
  as.data.frame() %>% 
  cbind(cell_type = MPRNMPR_object_miMPRobe_remove@meta.data$cell_type_new)%>% 
  cbind(Patients = MPRNMPR_object_miMPRobe_remove@meta.data$Patients) %>% 
  cbind(Treatments = MPRNMPR_object_miMPRobe_remove@meta.data$Treatments) %>% 
  cbind(MPR_Response = MPRNMPR_object_miMPRobe_remove@meta.data$MPR_Response)%>%
  cbind(cell_subtype = MPRNMPR_object_miMPRobe_remove@meta.data$cell_subtype)

names(MPRNMPR.umap.harmonyredu)

####重新画cell_subtype图----
# cell_type_color <- c("#FF6347","#ffc556","#1E90FF","#0000FF","#9370DB","#00FFFF","#32cd32","#006400","#FF00ff","#20B2AA")
library(randomcoloR)
cell_subtype_color<- distinctColorPalette(43)

unique(MPRNMPR.umap.harmonyredu$cell_subtype)
subtype <- c("CD4_C1-PDCD1","CD4_C2-IL7R","CD4_C3-FOXP3", "CD4_C4-FOSB","CD8_C1-GNLY","CD8_C2-CTSW","CD8_C3-IL7R","CD8_C4-GZMK","CD8_C5-STMN1", "CD3+CD4-CD8-","NK", #T cells
             "Pre-B cells","Primarily mature naive B cells","Naive B cells","Active B cells","Naive plasma cell", "Plasma cells",#B cells
             "Mast cells","SPP1+TAMs","FOLR2+TAMs", "ZBTB16+TAMs", "CXCL10+TAMs","Neutrophil","Monocyte-derived dendritic cells","cDC1","cDC2","Migratory cDCs", #myeloid cells
             "CXCL2+iCAFs","CXCR4+iCAFs","TPSAB1+iCAFs","MyCAFs","apCAFs","eCAFs", #Fibroblast
             "Non-malignant epithelial cells","Malignant epithelial cells", #上皮细胞恶性非恶性
             "Lymphatic ECs","Capillary-arterial ECs", "ACKR1+Blood ECs", "PTPRC+Blood ECs","MKI67+Blood ECs", "SPP1+ECs", #Endothelial cells
             "Smooth muscle cells",
             "Pericytes")

MPRNMPR.umap.harmonyredu$cell_subtype <- factor(MPRNMPR.umap.harmonyredu$cell_subtype,levels=subtype)

MPRNMPR_cell_type <- ggplot(MPRNMPR.umap.harmonyredu,aes(x= umapharmony_1, y = umapharmony_2,color = cell_subtype)) +  
  geom_point(size = 0.5, alpha =0.5)  +  
  scale_color_manual(values = cell_subtype_color)+
  xlab("UMAP1")+
  ylab("UMAP2")+
  theme(panel.grid.major = element_blank(), #主网格线
        panel.grid.minor = element_blank(), #次网格线
        panel.border = element_rect(color = "black", linewidth=2, fill = NA), #边框
        axis.title = element_blank(),  #轴标题
        axis.text = element_blank(), # 文本
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = 'white'), #背景色
        plot.background=element_rect(fill="white"))+
  theme(legend.title = element_blank(), #去掉legend.title 
        #legend.key=element_rect(fill='white'), #
        legend.key = element_blank(),
        legend.text = element_text(size=18), #设置legend标签的大小
        legend.key.size=unit(1,'cm') ) +  # 设置legend标签之间的大小
  guides(color = guide_legend(override.aes = list(size=8),ncol = 4))+ #设置legend中 点的大小 
  geom_segment(aes(x = min(MPRNMPR.umap.harmonyredu$umapharmony_1) , y = min(MPRNMPR.umap.harmonyredu$umapharmony_2) ,
                   xend = min(MPRNMPR.umap.harmonyredu$umapharmony_1) +3, yend = min(MPRNMPR.umap.harmonyredu$umapharmony_2) ),
               colour = "black", size=1,arrow = arrow(length = unit(0.5,"cm")))+ 
  geom_segment(aes(x = min(MPRNMPR.umap.harmonyredu$umapharmony_1)  , y = min(MPRNMPR.umap.harmonyredu$umapharmony_2)  ,
                   xend = min(MPRNMPR.umap.harmonyredu$umapharmony_1) , yend = min(MPRNMPR.umap.harmonyredu$umapharmony_2) + 3),
               colour = "black", size=1,arrow = arrow(length = unit(0.5,"cm"))) +
  annotate("text", x = min(MPRNMPR.umap.harmonyredu$umapharmony_1) +1.5, y = min(MPRNMPR.umap.harmonyredu$umapharmony_2) -1, label = "UMAP_1",
           color="black",size = 5, fontface="bold" ) + 
  annotate("text", x = min(MPRNMPR.umap.harmonyredu$umapharmony_1) -1, y = min(MPRNMPR.umap.harmonyredu$umapharmony_2) + 1.5, label = "UMAP_2",
           color="black",size = 5, fontface="bold" ,angle=90)+
  theme(legend.position = "bottom")

MPRNMPR_cell_type

ggsave('D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞注释_dotplot_featureplot/第二次umap图/MPRNMPR_cell_subtype.jpg',MPRNMPR_cell_type,width =15,height = 16)

ggsave('D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞注释_dotplot_featureplot/第二次umap图/MPRNMPR_cell_subtype.pdf',MPRNMPR_cell_type,width =15,height = 16)


###颜色设置###画小坐标轴####https://cloud.tencent.com/developer/article/1924260
####重新画cell_type umap----

cell_type_color <- c("#FF6347","#ffc556","#1E90FF","#0000FF","#9370DB","#00FFFF","#32cd32","#006400","#FF00ff","#20B2AA")

MPRNMPR.umap.harmonyredu$cell_type<- factor(MPRNMPR.umap.harmonyredu$cell_type, levels = c('T cells', 'B cells','Plasma cells', 'Myeloid cells',  
                                                                                           'Mast cells',"Endothelial cells","Fibroblasts","Smooth muscle cells",
                                                                                           "Pericytes","Epithelial cells"))
# library(RColorBrewer)
# library(ggplot2)
# par(mar=c(3,4,2,2))
# display.brewer.all()
# br_pal <- brewer.pal(10,"Paired")



MPRNMPR_cell_type <- ggplot(MPRNMPR.umap.harmonyredu,aes(x= umapharmony_1, y = umapharmony_2,color = cell_type)) +  
  geom_point(size = 0.5, alpha =0.5)  +  
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
        #legend.key=element_rect(fill='white'), #
        legend.key = element_blank(),
        legend.text = element_text(size=20), #设置legend标签的大小
        legend.key.size=unit(1,'cm') ) +  # 设置legend标签之间的大小
  guides(color = guide_legend(override.aes = list(size=5)))+ #设置legend中 点的大小 
  geom_segment(aes(x = min(MPRNMPR.umap.harmonyredu$umapharmony_1) , y = min(MPRNMPR.umap.harmonyredu$umapharmony_2) ,
                   xend = min(MPRNMPR.umap.harmonyredu$umapharmony_1) +3, yend = min(MPRNMPR.umap.harmonyredu$umapharmony_2) ),
               colour = "black", size=1,arrow = arrow(length = unit(0.5,"cm")))+ 
  geom_segment(aes(x = min(MPRNMPR.umap.harmonyredu$umapharmony_1)  , y = min(MPRNMPR.umap.harmonyredu$umapharmony_2)  ,
                   xend = min(MPRNMPR.umap.harmonyredu$umapharmony_1) , yend = min(MPRNMPR.umap.harmonyredu$umapharmony_2) + 3),
               colour = "black", size=1,arrow = arrow(length = unit(0.5,"cm"))) +
  annotate("text", x = min(MPRNMPR.umap.harmonyredu$umapharmony_1) +1.5, y = min(MPRNMPR.umap.harmonyredu$umapharmony_2) -1, label = "UMAP_1",
           color="black",size = 7, fontface="bold" ) + 
  annotate("text", x = min(MPRNMPR.umap.harmonyredu$umapharmony_1) -1, y = min(MPRNMPR.umap.harmonyredu$umapharmony_2) + 1.5, label = "UMAP_2",
           color="black",size = 7, fontface="bold" ,angle=90)

MPRNMPR_cell_type

MPRNMPR_cell_type_med <- MPRNMPR.umap.harmonyredu %>%
  group_by(cell_type) %>%
  summarise(
    umapharmony_1 = median(umapharmony_1),
    umapharmony_2 = median(umapharmony_2)
  )

library(ggrepel)
MPRNMPR_cell_type_med

MPRNMPR_cell_type2<-MPRNMPR_cell_type +geom_label_repel(aes(label=cell_type),size=7,color="black",fontface="bold",data = MPRNMPR_cell_type_med,
                                                        point.padding=unit(0.5, "lines"),fill = alpha(c("white"),0.5),
                                                        segment.size=0.5, nudge_x=0.5, nudge_y = 0,direction="y",max.overlaps=50)+
  theme(legend.position = "none")


ggsave('D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞注释_dotplot_featureplot/第二次umap图/MPRNMPR_cell_type2.jpg',MPRNMPR_cell_type2,width =7,height = 6)



#####Patients

Patients_color <- c("#40e0d0","#ee82ee","#7ccd7c","#551a8b","#ffc0cb","#ff6347",
                    "#f4a460","#b2dfee","#00ff7f","#63b8ff")
##调整legend 顺序

MPRNMPR.umap.harmonyredu$Patients <- factor(MPRNMPR.umap.harmonyredu$Patients, levels=c('P1', 'P2', 'P3', 'P4', 
                                                                                        'P5',"P6","P7","P8","P9","P10"))

MPRNMPR_Patients <- ggplot(MPRNMPR.umap.harmonyredu,aes(x= umapharmony_1, y = umapharmony_2,color = Patients)) +  
  geom_point(size = 0.2, alpha =1)  +  
  scale_color_manual(values = Patients_color)+
  xlab("UMAP1")+
  ylab("UMAP2")+
  theme(legend.position = "bottom")+
  guides(color = guide_legend(nrow = 4))+
  theme(panel.grid.major = element_blank(), #主网格线
        panel.grid.minor = element_blank(), #次网格线
        #panel.border = element_blank(), #边框
        axis.title = element_blank(),  #轴标题
        axis.text = element_blank(), # 文本
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = 'white', color="black",linewidth=1), #背景色
        plot.background=element_rect(fill="white"))+
  theme(legend.title = element_blank(), #去掉legend.title 
        #legend.key=element_rect(fill='white'), #
        legend.key = element_blank(),
        legend.text = element_text(size=20), #设置legend标签的大小
        legend.key.size=unit(1,'cm') ) +  # 设置legend标签之间的大小
  guides(color = guide_legend(override.aes = list(size=5)))+ #设置legend中 点的大小 
  geom_segment(aes(x = min(MPRNMPR.umap.harmonyredu$umapharmony_1) , y = min(MPRNMPR.umap.harmonyredu$umapharmony_2) ,
                   xend = min(MPRNMPR.umap.harmonyredu$umapharmony_1) +3, yend = min(MPRNMPR.umap.harmonyredu$umapharmony_2) ),
               colour = "black", size=0.5,arrow = arrow(length = unit(0.3,"cm")))+ 
  geom_segment(aes(x = min(MPRNMPR.umap.harmonyredu$umapharmony_1)  , y = min(MPRNMPR.umap.harmonyredu$umapharmony_2)  ,
                   xend = min(MPRNMPR.umap.harmonyredu$umapharmony_1) , yend = min(MPRNMPR.umap.harmonyredu$umapharmony_2) + 3),
               colour = "black", size=0.5,arrow = arrow(length = unit(0.3,"cm"))) +
  annotate("text", x = min(MPRNMPR.umap.harmonyredu$umapharmony_1) +1.5, y = min(MPRNMPR.umap.harmonyredu$umapharmony_2) -1, label = "UMAP_1",
           color="black",size = 3, fontface="bold" ) + 
  annotate("text", x = min(MPRNMPR.umap.harmonyredu$umapharmony_1) -1, y = min(MPRNMPR.umap.harmonyredu$umapharmony_2) + 1.5, label = "UMAP_2",
           color="black",size = 3, fontface="bold" ,angle=90)

MPRNMPR_Patients



Treatments_color <- c("#32cd32","#b9d3ee")

MPRNMPR_Treatments <- ggplot(MPRNMPR.umap.harmonyredu,aes(x= umapharmony_1, y = umapharmony_2,color = Treatments)) +  
  geom_point(size = 0.2, alpha =1)  +  
  scale_color_manual(values = Treatments_color)+
  #xlab("UMAP1")+
  #ylab("UMAP2")+
  theme(legend.position = "bottom")+
  guides(color = guide_legend(nrow = 1))+
  theme(panel.grid.major = element_blank(), #主网格线
        panel.grid.minor = element_blank(), #次网格线
        #panel.border = element_blank(), #边框
        axis.title = element_blank(),  #轴标题
        axis.text = element_blank(), # 文本
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = 'white', color="black",linewidth=1), #背景色
        plot.background=element_rect(fill="white"))+
  theme(legend.title = element_blank(), #去掉legend.title 
        #legend.key=element_rect(fill='white'),
        legend.key = element_blank(),#
        legend.text = element_text(size=20), #设置legend标签的大小
        legend.key.size=unit(1,'cm')) +  # 设置legend标签之间的大小
  guides(color = guide_legend(override.aes = list(size=5)))+ #设置legend中 点的大小 
  geom_segment(aes(x = min(MPRNMPR.umap.harmonyredu$umapharmony_1) , y = min(MPRNMPR.umap.harmonyredu$umapharmony_2) ,
                   xend = min(MPRNMPR.umap.harmonyredu$umapharmony_1) +3, yend = min(MPRNMPR.umap.harmonyredu$umapharmony_2) ),
               colour = "black", size=0.5,arrow = arrow(length = unit(0.3,"cm")))+ 
  geom_segment(aes(x = min(MPRNMPR.umap.harmonyredu$umapharmony_1)  , y = min(MPRNMPR.umap.harmonyredu$umapharmony_2)  ,
                   xend = min(MPRNMPR.umap.harmonyredu$umapharmony_1) , yend = min(MPRNMPR.umap.harmonyredu$umapharmony_2) + 3),
               colour = "black", size=0.5,arrow = arrow(length = unit(0.3,"cm"))) +
  annotate("text", x = min(MPRNMPR.umap.harmonyredu$umapharmony_1) +1.5, y = min(MPRNMPR.umap.harmonyredu$umapharmony_2) -1, label = "UMAP_1",
           color="black",size = 3, fontface="bold" ) + 
  annotate("text", x = min(MPRNMPR.umap.harmonyredu$umapharmony_1) -1, y = min(MPRNMPR.umap.harmonyredu$umapharmony_2) + 1.5, label = "UMAP_2",
           color="black",size = 3, fontface="bold" ,angle=90)

MPRNMPR_Treatments

library(ggrepel)

MPR_Response_color <- c("#ff8247","#a020f0")


MPRNMPR_MPR_Response <- ggplot(MPRNMPR.umap.harmonyredu,aes(x= umapharmony_1, y = umapharmony_2,color = MPR_Response)) +  
  geom_point(size = 0.2, alpha =1)  +  
  scale_color_manual(values = MPR_Response_color)+
  #xlab("UMAP1")+
  #ylab("UMAP2")+
  theme(legend.position = "bottom")+
  guides(color = guide_legend(nrow = 1))+
  theme(panel.grid.major = element_blank(), #主网格线
        panel.grid.minor = element_blank(), #次网格线
        #panel.border = element_blank(), #边框
        axis.title = element_blank(),  #轴标题
        axis.text = element_blank(), # 文本
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = 'white', color="black",linewidth=1), #背景色
        plot.background=element_rect(fill="white"))+
  theme(legend.title = element_blank(), #去掉legend.title 
        #legend.key=element_rect(fill='white'), #
        legend.key = element_blank(),
        legend.text = element_text(size=20), #设置legend标签的大小
        legend.key.size=unit(1,'cm') ) +  # 设置legend标签之间的大小
  guides(color = guide_legend(override.aes = list(size=5)))+ #设置legend中 点的大小 
  geom_segment(aes(x = min(MPRNMPR.umap.harmonyredu$umapharmony_1) , y = min(MPRNMPR.umap.harmonyredu$umapharmony_2) ,
                   xend = min(MPRNMPR.umap.harmonyredu$umapharmony_1) +3, yend = min(MPRNMPR.umap.harmonyredu$umapharmony_2) ),
               colour = "black", size=0.5,arrow = arrow(length = unit(0.3,"cm")))+ 
  geom_segment(aes(x = min(MPRNMPR.umap.harmonyredu$umapharmony_1)  , y = min(MPRNMPR.umap.harmonyredu$umapharmony_2)  ,
                   xend = min(MPRNMPR.umap.harmonyredu$umapharmony_1) , yend = min(MPRNMPR.umap.harmonyredu$umapharmony_2) + 3),
               colour = "black", size=0.5,arrow = arrow(length = unit(0.3,"cm"))) +
  annotate("text", x = min(MPRNMPR.umap.harmonyredu$umapharmony_1) +1.5, y = min(MPRNMPR.umap.harmonyredu$umapharmony_2) -1, label = "UMAP_1",
           color="black",size = 3, fontface="bold" ) + 
  annotate("text", x = min(MPRNMPR.umap.harmonyredu$umapharmony_1) -1, y = min(MPRNMPR.umap.harmonyredu$umapharmony_2) + 1.5, label = "UMAP_2",
           color="black",size = 3, fontface="bold" ,angle=90)

MPRNMPR_MPR_Response


library(cowplot)

bottom_row <- plot_grid(MPRNMPR_Patients,MPRNMPR_Treatments,MPRNMPR_MPR_Response, nrow =1, align = 'h', rel_widths = c(1,1,1),rel_heights = c(1,1,1))

MPRNMPR_cell_type_merged <-plot_grid(MPRNMPR_cell_type2,bottom_row,
                                     ncol = 1,
                                     #labels = c('A','空图','B'),
                                     rel_widths = c(1,1.3),
                                     rel_heights = c(1,0.7))

ggsave('D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞注释_dotplot_featureplot/第二次umap图/MPRNMPR_cell_type_merged.jpg',MPRNMPR_cell_type_merged,width =15,height = 15)
ggsave('D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞注释_dotplot_featureplot/第二次umap图/MPRNMPR_cell_type_merged.pdf',MPRNMPR_cell_type_merged,width =15,height = 15)


####
####细胞比例九宫格图----
#画图使用数据为长数据，第一列group，第二列物种或者细胞类型，第三列为value
Cellratio <- prop.table(table(MPRNMPR_object_miMPRobe_remove@meta.data$cell_type_new,MPRNMPR_object_miMPRobe_remove@meta.data$MPR_Response2), margin = 2)
Cellratio <- as.data.frame(Cellratio)
View(Cellratio)
colnames(Cellratio) <- c("Cluster","Response","Freq")

write.csv(Cellratio,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞注释_dotplot_featureplot/第二次umap图/Cellratio_细胞比例饼图数据.csv")
###加载
library(tidyverse)
library(ggtext)
library(glue)
# devtools::install_github("AllanCameron/VoronoiPlus")
library(VoronoiPlus) 

au_vor1 <- voronoi_treemap(Freq ~ Cluster,
                           data = Cellratio %>% filter(Response=="Pre_NMPR")) %>% 
  mutate(type="Pre_NMPR")

View(au_vor1)

au_vor2 <- voronoi_treemap(Freq ~ Cluster,
                           data = Cellratio %>% filter(Response=="Pre_MPR")) %>% 
  mutate(type="Pre_MPR")

au_vor3 <- voronoi_treemap(Freq ~ Cluster,
                           data = Cellratio %>% filter(Response=="Post_MPR")) %>% 
  mutate(type="Post_MPR")

groups <- au_vor1 %>% bind_rows(au_vor2,au_vor3)

groups$type<- factor(groups$type, levels = c('Pre_NMPR', 'Pre_MPR','Post_MPR'))
groups$group<- factor(groups$group, levels = c('T cells', 'B cells','Plasma cells', 'Myeloid cells',  
                                               'Mast cells',"Endothelial cells","Fibroblasts","Smooth muscle cells",
                                               "Pericytes","Epithelial cells"))
cell_type_color <- c("#FF6347","#ffc556","#1E90FF","#0000FF","#9370DB","#00FFFF","#32cd32","#006400","#FF00ff","#20B2AA")

label_group<-groups%>% group_by(group,type,value) %>% summarise(x_median=median(x),y_median=median(y))
nrow(value_group)
nrow(label_group)
View(label_group)

label_group$value <- paste(round(label_group$value *100,2),sep="","%")
pie_plot<- ggplot() +
  geom_polygon(data = groups,mapping = aes(x = x, y = y, group = group, fill = group),
               colour = "white",linewidth =0.8) +
  facet_wrap(.~type,scale="free")+
  scale_fill_manual(values =cell_type_color)+
  geom_text(aes(x_median,y_median,label=value),size=7,color="black",fontface="bold",data = label_group)+
  labs(x=NULL,y=NULL)+
  theme_test()+
  theme(axis.text=element_blank(),
        axis.ticks = element_blank(),
        plot.background = element_rect(fill = "white", colour = "white"),
        panel.background = element_rect(fill = "white", colour = "white"),
        panel.spacing.x=unit(-0.1,"cm"),
        panel.spacing.y=unit(-0.1,"cm"),
        plot.margin = margin(2,2,0,0),
        legend.key.height = unit(0.2,"in"), 
        legend.key.width = unit(0.2,"in"),
        legend.title = element_blank(),
        legend.text = element_text(margin = margin(l=0,unit="cm"),size=15),
        strip.background = element_rect(fill="grey96"),
        strip.text = element_text(face="bold",size=19))+
  theme(legend.position = "bottom")
  

pie_plot


ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞注释_dotplot_featureplot/第二次umap图/MPR_pie_细胞比例2.jpg",pie_plot,width=16.5,height=6.5)

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞注释_dotplot_featureplot/第二次umap图/MPR_pie_细胞比例2.pdf",pie_plot,width=16.5,height=6.5)

write.csv(groups,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞注释_dotplot_featureplot/第二次umap图/细胞比例饼图长数据groups.csv")


######重新featureplot----

setwd("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞注释_dotplot_featureplot")

library(Seurat)
library(ggplot2)
library(cowplot)
library(dplyr)
library(tidydr)
library(stringr)
library(viridis)

markers <- c("TRAC", # Tcells
             "MS4A1", #B cells
             "IGLC2",#Plasma cell
             "C1QB", # Myeloid cells
             "TPSB2", # Mast cells
             "PROX1","VEGFC",#Endothelial cells
             "COL12A1", # Fibroblasts
             "ACTA1", #Smooth muscle cells
             "HIGD1B",#Pericytes
             "KRT15",
             "PDCD1") # Endothelial cells

pal <- viridis(n = 10, option = "C")
# "#0D0887FF" "#47039FFF" "#7301A8FF" "#9C179EFF" "#BD3786FF" "#D8576BFF" "#ED7953FF" "#FA9E3BFF" "#FDC926FF" "#F0F921FF"
pal <- viridis(n = 15, option = "D", direction = -1)


p<-FeaturePlot(object =MPRNMPR_object_miMPRobe_remove,features = markers,reduction = "umap.harmony",raster=FALSE, cols = pal, order = T)

ggsave(file="D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞注释_dotplot_featureplot/第二次umap图/MPR细胞大类markerfeatureplot.jpg",p,width=15,height=9)
ggsave(file="D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞注释_dotplot_featureplot/第二次umap图/MPR细胞大类markerfeatureplot.pdf",p,width=15,height=9)


i=1
plots=list()
for (i in 1:length(markers)){
  plots[[i]]=FeaturePlot(object=MPRNMPR_object_miMPRobe,features = markers[i],reduction = "umap.harmony",raster=FALSE, cols = pal, order = T)+
    xlab("UMAP_1")+
    ylab("UMAP_1")+
    theme(legend.position = "none")
}
library(patchwork)
p<-wrap_plots(plots, ncol = 4)+plot_annotation(tag_levels = "A");p

p_legend<-FeaturePlot(object=MPRNMPR_object_miMPRobe_remove,features = "PDCD1",reduction = "umap.harmony",raster=FALSE, cols = pal, order = T)+
  xlab("UMAP_1")+
  ylab("UMAP_1")+
  theme(legend.position = "right")

legend <- get_legend(p_legend+
                       #guides(color = guide_legend(nrow = 1)) +
                       theme(legend.position = "right"))
plt_featureplot2<- plot_grid(p,legend, ncol=2,rel_widths = c(1,0.1))

ggsave(file="D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞注释_dotplot_featureplot/第二次umap图/MPR细胞大类markerfeatureplot_legend.jpg",plt_featureplot2,width=15,height=9)
ggsave(file="D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞注释_dotplot_featureplot/第二次umap图/MPR细胞大类markerfeatureplot_legend.pdf",plt_featureplot2,width=15,height=9)

#ggsave("T_cells_CD3E.jpg",device = "pdf",width = 10,height = 10.5,units = "cm")


####重新画细胞比例差异统计图######################
Cellratio <- prop.table(table(MPRNMPR_object_miMPRobe_remove@meta.data$cell_type_new,MPRNMPR_object_miMPRobe_remove@meta.data$orig.ident), margin = 2)

write.csv(Cellratio,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞比例_OR值计算/每个病人Cellratio.csv")

Cellratio <- as.data.frame(Cellratio)
colnames(Cellratio) <- c("Cluster","Group","Freq")
library(reshape2)
cellper <- dcast(Cellratio,Group~Cluster, value.var = "Freq")#长数据转为宽数据
write.csv(cellper,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞比例scCODA/细胞大类比例.csv")
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


cell_type_groups = c('T cells', 'B cells','Plasma cells', 'Myeloid cells',  
                     'Mast cells',"Endothelial cells","Fibroblasts","Smooth muscle cells",
                     "Pericytes","Epithelial cells")
library(ggplot2)
library(dplyr)
library(ggpubr)
library(cowplot)
######MPR_NMPR 颜色----

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
    geom_boxplot(aes(fill=MPR_Response2),outlier.shape = NA,lwd= 1)+
    #geom_boxplot(outlier.colour="red", outlier.shape=7,outlier.size=1) +
    scale_fill_manual(values = treatment_color)+
    #geom_jitter(shape = 21,aes(fill=Treatments$Treatments),width = 0.25) + 
    #stat_summary(fun=mean, geom="point", color="grey60") +
    theme_cowplot() +
    theme(axis.text = element_text(size = 25),
          axis.line.x=element_line(linetype=1,color="black",size=1),
          axis.ticks.x=element_line(color="black",size=2,lineend = 1),
          axis.line.y=element_line(linetype=1,color="black",size=1),
          axis.ticks.y=element_line(color="black",size=2,lineend = 1),
          axis.title = element_text(size = 25),
          legend.text = element_text(size = 25),
          axis.text.x = element_text(vjust=1,hjust=1,angle=30,size=25),
          axis.text.y = element_text(vjust=1,hjust=1,angle=0,size=20),
          legend.title = element_text(size = 25),
          plot.title = element_text(size = 25),
          legend.position = 'none') + 
    labs(title = cell_type,y='Percentage',x="") +
    guides(fill = guide_legend(title = NULL))
  #geom_errorbar(aes(ymin = lower, ymax = upper),col = "grey60",width =  1)
  
  ###组间t检验分析
  my_comparisons <- list(c("Pre_MPR", "Pre_NMPR"),c("Pre_MPR", "Post_MPR"),c("Pre_NMPR", "Post_MPR"))
  pp2<-pp1+stat_compare_means(comparisons=my_comparisons,
                              method="wilcox.test",size=10,bracket.size = 1)
  pplist[[cell_type]] = pp2
}

# library(cowplot)
cells_<-plot_grid(pplist[['T cells']],
                  pplist[['B cells']],
                  pplist[['Plasma cells']],
                  pplist[['Myeloid cells']],
                  pplist[['Mast cells']],
                  pplist[['Endothelial cells']],
                  pplist[['Fibroblasts']],
                  pplist[['Pericytes']],
                  pplist[['Smooth muscle cells']],
                  pplist[['Epithelial cells']],
                  align = "h",  #axis = 'l',
                  nrow =2)

plt1<- ggplot(cellper_cell_type,aes(x=MPR_Response,y=percent)) + #ggplot作图
  geom_boxplot(aes(fill=Treatments),outlier.shape = NA,lwd= 1)+
  #geom_boxplot(outlier.colour="red", outlier.shape=7,outlier.size=1) +
  scale_fill_manual(values = treatment_color)+
  #geom_jitter(shape = 21,aes(fill=Treatments$Treatments),width = 0.25) + 
  #stat_summary(fun=mean, geom="point", color="grey60") +
  theme_cowplot() +
  theme(axis.text = element_text(size = 25),
        legend.key.size = unit(50, "pt"),
        axis.title = element_text(size = 25,face="bold"),
        legend.text = element_text(size = 20,face="bold"),
        axis.text.x = element_text(vjust=1,hjust=1,angle=30,size=25,face="bold"),
        axis.text.y = element_text(vjust=1,hjust=1,angle=0,size=20,face="bold"),
        legend.title = element_text(size = 20),
        plot.title = element_text(size = 25,face="bold"),
        legend.position = "right") + 
  labs(title = cell_type,y='Percentage',x="") +
  guides(fill = guide_legend(title = NULL))

legend <- get_legend(plt1+
                       #guides(color = guide_legend(nrow = 1)) +
                       theme(legend.position = "right"))

cells_pre2<- plot_grid(cells_,legend,nrow=2,rel_widths=c(1,0.5))

ggsave('D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞比例_OR值计算/MPR_cells_percentage.jpg',cells_,width =30,height =25,limitsize = FALSE)

ggsave('D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/细胞比例_OR值计算/MPR_cells_percentage.pdf',cells_,width =30,height =25,limitsize = FALSE)

####infercnv分析----
###在上皮细胞

####微生物单细胞中分布情况---
library(dplyr)
invade_16s<-read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/有无微生物对免疫细胞和上皮细胞功能影响/16s_filter.combined.genus.umi.matrix.csv")
invade_16s$Massilia
MPRNMPR_object_miMPRobe_remove@meta.data$Massilia <- invade_16s[match(MPRNMPR_object_miMPRobe_remove@meta.data$cell_id, invade_16s$barcode),"Massilia"]
min(is.na(MPRNMPR_object_miMPRobe_remove@meta.data$Massilia))
sum(is.na(MPRNMPR_object_miMPRobe_remove@meta.data$Massilia))
MPRNMPR_object_miMPRobe_remove@meta.data$Massilia[is.na(MPRNMPR_object_miMPRobe_remove@meta.data$Massilia)] <- 0

MPRNMPR_object_miMPRobe_remove@meta.data$MPR_patient_group <- paste(MPRNMPR_object_miMPRobe_remove@meta.data$MPR_Response2,MPRNMPR_object_miMPRobe_remove@meta.data$Patients,sep="-")

library(tidyverse)
MPRNMPR.umap.harmonyredu<- MPRNMPR_object_miMPRobe_remove@reductions$umap.harmony@cell.embeddings%>%
  as.data.frame() %>% 
  cbind(cell_id = MPRNMPR_object_miMPRobe_remove@meta.data$cell_id)%>%
  cbind(MPR_response2 = MPRNMPR_object_miMPRobe_remove@meta.data$MPR_Response2) %>%
  cbind(Patients = MPRNMPR_object_miMPRobe_remove@meta.data$Patients) %>%
  cbind(MPR_patient_group = MPRNMPR_object_miMPRobe_remove@meta.data$MPR_patient_group)

MPRNMPR.umap.harmonyredu$Massilia <- invade_16s[match(MPRNMPR.umap.harmonyredu$cell_id, invade_16s$barcode),"Massilia"]
MPRNMPR.umap.harmonyredu$Massilia[is.na(MPRNMPR.umap.harmonyredu$Massilia)] <- 0
MPRNMPR.umap.harmonyredu_filter<-MPRNMPR.umap.harmonyredu[!MPRNMPR.umap.harmonyredu$Massilia=="0",]
max(MPRNMPR.umap.harmonyredu_filter$Massilia)
MPRNMPR.umap.harmonyredu$MPR_response2<- factor(MPRNMPR.umap.harmonyredu$MPR_response2,levels=c("Pre_NMPR","Pre_MPR","Post_MPR"))
MPRNMPR.umap.harmonyredu_filter$MPR_patient_group<-factor(MPRNMPR.umap.harmonyredu_filter$MPR_patient_group,
                                                   levels=c("Pre_NMPR-P3","Pre_NMPR-P5","Pre_NMPR-P6","Pre_NMPR-P10",
                                                             "Pre_MPR-P1","Pre_MPR-P2","Pre_MPR-P4","Pre_MPR-P7",
                                                           "Pre_MPR-P8","Pre_MPR-P9","Post_MPR-P1", "Post_MPR-P4",
                                                           "Post_MPR-P7","Post_MPR-P8"))


MPRNMPR_cell_type <- ggplot()+ 
  geom_point(data=MPRNMPR.umap.harmonyredu,aes(x=umapharmony_1, y = umapharmony_2),
             size = 0.5, alpha =1,fill="#c1cdc1",color="#c1cdc1")+
  theme_bw()+
  geom_point(data=MPRNMPR.umap.harmonyredu_filter,aes(x= umapharmony_1, y = umapharmony_2,color = Massilia),
             size =0.5, alpha =1)+
  #scale_colour_gradient2(low = "yellow",mid="blue", high = "red")+
  scale_color_gradientn(colors = c("blue","yellow", "red"), values = c(0,0.3, 0.5, 1), limits = c(1, 1650),breaks = c(1, seq(200,1650, by = 400)))+
  xlab("UMAP_1")+
  ylab("UMAP_2")+
  ggtitle("Massilia")+
  facet_wrap(~ MPR_patient_group) +
  labs(color = "UMI counts") +
  theme(panel.grid.major = element_blank(), #主网格线
        panel.grid.minor = element_blank(), #次网格线
        panel.border = element_blank(), #边框
        panel.background = element_rect(fill = 'white'), #背景色
        plot.background=element_rect(fill="white"),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold"), # Customize X-axis title
        axis.title.y = element_text(size = 12, face = "bold"), # Customize Y-axis title
        axis.text.x = element_text(size = 12), # Customize X-axis text
        axis.text.y = element_text(size = 12),
        axis.line = element_line(color = "black", linewidth = 0.8),  # 设置轴线的颜色和大小
        axis.line.x.bottom = element_line(color = "black", linewidth = 0.8),  # 设置底部 X 轴线条
        axis.line.y.left = element_line(color = "black", linewidth = 0.8),
        strip.text = element_text(size = 10, face = "bold", color = "black")) +# 设置左边 Y 轴线条)+
  theme(#legend.title = element_blank(), #去掉legend.title 
    legend.key=element_rect(fill='white'), #
    legend.text = element_text(size=10), #设置legend标签的大小
    legend.key.size=unit(0.5,'cm') ) +  # 设置legend标签之间的大小
  theme(legend.position = "right")
MPRNMPR_cell_type
ggsave('D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/INVADE-seq/微生物featureplot/Massilia.jpg',MPRNMPR_cell_type,width =8,height =7,limitsize = FALSE)
ggsave('D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/INVADE-seq/微生物featureplot/Massilia分组.jpg',MPRNMPR_cell_type,width =8,height =7,limitsize = FALSE)

ggsave('D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/INVADE-seq/微生物featureplot/Massilia_病人分面.jpg',MPRNMPR_cell_type,width =8,height =8,limitsize = FALSE)
ggsave('D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/INVADE-seq/微生物featureplot/Massilia_响应分面.jpg',MPRNMPR_cell_type,width =12,height =5,limitsize = FALSE)


####Herbaspirillum----
invade_16s$Herbaspirillum

MPRNMPR.umap.harmonyredu$Herbaspirillum <- invade_16s[match(MPRNMPR.umap.harmonyredu$cell_id, invade_16s$barcode),"Herbaspirillum"]
MPRNMPR.umap.harmonyredu$Herbaspirillum[is.na(MPRNMPR.umap.harmonyredu$Herbaspirillum)] <- 0

MPRNMPR.umap.harmonyredu_filter<-MPRNMPR.umap.harmonyredu[!MPRNMPR.umap.harmonyredu$Herbaspirillum=="0",]
max(MPRNMPR.umap.harmonyredu_filter$Herbaspirillum)
min(MPRNMPR.umap.harmonyredu_filter$Herbaspirillum)
MPRNMPR.umap.harmonyredu$MPR_response2<- factor(MPRNMPR.umap.harmonyredu$MPR_response2,levels=c("Pre_NMPR","Pre_MPR","Post_MPR"))

MPRNMPR.umap.harmonyredu_filter$MPR_patient_group<-factor(MPRNMPR.umap.harmonyredu_filter$MPR_patient_group,
                                                          levels=c("Pre_NMPR-P3","Pre_NMPR-P5","Pre_NMPR-P6","Pre_NMPR-P10",
                                                                   "Pre_MPR-P1","Pre_MPR-P2","Pre_MPR-P4","Pre_MPR-P7",
                                                                   "Pre_MPR-P8","Pre_MPR-P9","Post_MPR-P1", "Post_MPR-P4",
                                                                   "Post_MPR-P7","Post_MPR-P8"))
MPRNMPR_cell_type <- ggplot()+ 
  geom_point(data=MPRNMPR.umap.harmonyredu,aes(x=umapharmony_1, y = umapharmony_2),
             size = 0.5, alpha =1,fill="#c1cdc1",color="#c1cdc1")+
  theme_bw()+
  geom_point(data=MPRNMPR.umap.harmonyredu_filter,aes(x= umapharmony_1, y = umapharmony_2,color = Herbaspirillum),
             size =0.5, alpha =1)+
  #scale_colour_gradient2(low = "yellow",mid="blue", high = "red")+
  scale_color_gradientn(colors = c("blue","yellow", "red"), values = c(0,0.3, 0.5, 1), limits = c(1, 77),breaks = c(1, seq(20, 77, by = 20)))+
  xlab("UMAP_1")+
  ylab("UMAP_2")+
  ggtitle("Herbaspirillum")+
  labs(color = "UMI counts") +
  facet_wrap(~ MPR_patient_group) +
  theme(panel.grid.major = element_blank(), #主网格线
        panel.grid.minor = element_blank(), #次网格线
        panel.border = element_blank(), #边框
        panel.background = element_rect(fill = 'white'), #背景色
        plot.background=element_rect(fill="white"),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold"), # Customize X-axis title
        axis.title.y = element_text(size = 12, face = "bold"), # Customize Y-axis title
        axis.text.x = element_text(size = 12), # Customize X-axis text
        axis.text.y = element_text(size = 12),
        axis.line = element_line(color = "black", linewidth = 0.8),  # 设置轴线的颜色和大小
        axis.line.x.bottom = element_line(color = "black", linewidth = 0.8),  # 设置底部 X 轴线条
        axis.line.y.left = element_line(color = "black", linewidth = 0.8),
        strip.text = element_text(size = 10, face = "bold", color = "black")) +# 设置左边 Y 轴线条)+
  theme(#legend.title = element_blank(), #去掉legend.title 
    legend.key=element_rect(fill='white'), #
    legend.text = element_text(size=10), #设置legend标签的大小
    legend.key.size=unit(0.5,'cm') ) +  # 设置legend标签之间的大小
  theme(legend.position = "right")
MPRNMPR_cell_type

ggsave('D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/INVADE-seq/微生物featureplot/Herbaspirillum.jpg',MPRNMPR_cell_type,width =8,height =7,limitsize = FALSE)
ggsave('D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/INVADE-seq/微生物featureplot/Herbaspirillum_分组.jpg',MPRNMPR_cell_type,width =8,height =7,limitsize = FALSE)

ggsave('D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/INVADE-seq/微生物featureplot/Herbaspirillum_病人分面.jpg',MPRNMPR_cell_type,width =8,height =8,limitsize = FALSE)
ggsave('D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/INVADE-seq/微生物featureplot/Herbaspirillum_响应分面.jpg',MPRNMPR_cell_type,width =12,height =5,limitsize = FALSE)

####Streptococcus----
MPRNMPR.umap.harmonyredu$Streptococcus <- invade_16s[match(MPRNMPR.umap.harmonyredu$cell_id, invade_16s$barcode),"Streptococcus"]
MPRNMPR.umap.harmonyredu$Streptococcus[is.na(MPRNMPR.umap.harmonyredu$Streptococcus)] <- 0

MPRNMPR.umap.harmonyredu_filter<-MPRNMPR.umap.harmonyredu[!MPRNMPR.umap.harmonyredu$Streptococcus=="0",]
max(MPRNMPR.umap.harmonyredu_filter$Streptococcus)
min(MPRNMPR.umap.harmonyredu_filter$Streptococcus)
MPRNMPR.umap.harmonyredu$MPR_response2<- factor(MPRNMPR.umap.harmonyredu$MPR_response2,levels=c("Pre_NMPR","Pre_MPR","Post_MPR"))
MPRNMPR.umap.harmonyredu_filter$MPR_patient_group<-factor(MPRNMPR.umap.harmonyredu_filter$MPR_patient_group,
                                                          levels=c("Pre_NMPR-P3","Pre_NMPR-P5","Pre_NMPR-P6","Pre_NMPR-P10",
                                                                   "Pre_MPR-P1","Pre_MPR-P2","Pre_MPR-P4","Pre_MPR-P7",
                                                                   "Pre_MPR-P8","Pre_MPR-P9","Post_MPR-P1", "Post_MPR-P4",
                                                                   "Post_MPR-P7","Post_MPR-P8"))
MPRNMPR_cell_type <- ggplot()+ 
  geom_point(data=MPRNMPR.umap.harmonyredu,aes(x=umapharmony_1, y = umapharmony_2),
             size = 0.5, alpha =1,fill="#c1cdc1",color="#c1cdc1")+
  theme_bw()+
  geom_point(data=MPRNMPR.umap.harmonyredu_filter,aes(x= umapharmony_1, y = umapharmony_2,color = Streptococcus),
             size =0.7, alpha =1)+
  #scale_colour_gradient2(low = "yellow",mid="blue", high = "red")+
  scale_color_gradientn(colors = c("blue","yellow", "red"), values = c(0,0.3, 0.5, 1), limits = c(1, 7),breaks = c(1, seq(2, 7, by = 3)))+
  xlab("UMAP_1")+
  ylab("UMAP_2")+
  ggtitle("Streptococcus")+
  labs(color = "UMI counts") +
  facet_wrap(~ MPR_patient_group) +
  theme(panel.grid.major = element_blank(), #主网格线
        panel.grid.minor = element_blank(), #次网格线
        panel.border = element_blank(), #边框
        panel.background = element_rect(fill = 'white'), #背景色
        plot.background=element_rect(fill="white"),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold"), # Customize X-axis title
        axis.title.y = element_text(size = 12, face = "bold"), # Customize Y-axis title
        axis.text.x = element_text(size = 12), # Customize X-axis text
        axis.text.y = element_text(size = 12),
        axis.line = element_line(color = "black", linewidth = 0.8),  # 设置轴线的颜色和大小
        axis.line.x.bottom = element_line(color = "black", linewidth = 0.8),  # 设置底部 X 轴线条
        axis.line.y.left = element_line(color = "black", linewidth = 0.8),
        strip.text = element_text(size = 10, face = "bold", color = "black")) +# 设置左边 Y 轴线条)+
  theme(#legend.title = element_blank(), #去掉legend.title 
    legend.key=element_rect(fill='white'), #
    legend.text = element_text(size=10), #设置legend标签的大小
    legend.key.size=unit(0.5,'cm') ) +  # 设置legend标签之间的大小
  theme(legend.position = "right")
MPRNMPR_cell_type

ggsave('D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/INVADE-seq/微生物featureplot/Streptococcus.jpg',MPRNMPR_cell_type,width =8,height =7,limitsize = FALSE)
ggsave('D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/INVADE-seq/微生物featureplot/Streptococcus_分组.jpg',MPRNMPR_cell_type,width =8,height =7,limitsize = FALSE)
ggsave('D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/INVADE-seq/微生物featureplot/Streptococcus_病人分面.jpg',MPRNMPR_cell_type,width =8,height =8,limitsize = FALSE)
ggsave('D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/INVADE-seq/微生物featureplot/HStreptococcus_响应分面.jpg',MPRNMPR_cell_type,width =12,height =5,limitsize = FALSE)

####Corynebacterium,Fusobacterium,Rothia, Treponema,sphingomonas

MPRNMPR.umap.harmonyredu$Corynebacterium <- invade_16s[match(MPRNMPR.umap.harmonyredu$cell_id, invade_16s$barcode),"Corynebacterium"]
MPRNMPR.umap.harmonyredu$Corynebacterium[is.na(MPRNMPR.umap.harmonyredu$Corynebacterium)] <- 0

MPRNMPR.umap.harmonyredu_filter<-MPRNMPR.umap.harmonyredu[!MPRNMPR.umap.harmonyredu$Corynebacterium=="0",]
max(MPRNMPR.umap.harmonyredu_filter$Corynebacterium)
min(MPRNMPR.umap.harmonyredu_filter$Corynebacterium)
MPRNMPR.umap.harmonyredu$MPR_response2<- factor(MPRNMPR.umap.harmonyredu$MPR_response2,levels=c("Pre_NMPR","Pre_MPR","Post_MPR"))
MPRNMPR.umap.harmonyredu_filter$MPR_patient_group<-factor(MPRNMPR.umap.harmonyredu_filter$MPR_patient_group,
                                                          levels=c("Pre_NMPR-P3","Pre_NMPR-P5","Pre_NMPR-P6","Pre_NMPR-P10",
                                                                   "Pre_MPR-P1","Pre_MPR-P2","Pre_MPR-P4","Pre_MPR-P7",
                                                                   "Pre_MPR-P8","Pre_MPR-P9","Post_MPR-P1", "Post_MPR-P4",
                                                                   "Post_MPR-P7","Post_MPR-P8"))
MPRNMPR_cell_type <- ggplot()+ 
  geom_point(data=MPRNMPR.umap.harmonyredu,aes(x=umapharmony_1, y = umapharmony_2),
             size = 0.5, alpha =1,fill="#c1cdc1",color="#c1cdc1")+
  theme_bw()+
  geom_point(data=MPRNMPR.umap.harmonyredu_filter,aes(x= umapharmony_1, y = umapharmony_2,color = Corynebacterium),
             size =0.7, alpha =1)+
  #scale_colour_gradient2(low = "yellow",mid="blue", high = "red")+
  scale_color_gradientn(colors = c("blue","yellow", "red"), values = c(0,0.3, 0.5, 1), limits = c(1, 3),breaks = c(1, seq(2, 3, by = 3)))+
  xlab("UMAP_1")+
  ylab("UMAP_2")+
  ggtitle("Corynebacterium")+
  labs(color = "UMI counts") +
  #facet_wrap(~ MPR_patient_group) +
  theme(panel.grid.major = element_blank(), #主网格线
        panel.grid.minor = element_blank(), #次网格线
        panel.border = element_blank(), #边框
        panel.background = element_rect(fill = 'white'), #背景色
        plot.background=element_rect(fill="white"),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold"), # Customize X-axis title
        axis.title.y = element_text(size = 12, face = "bold"), # Customize Y-axis title
        axis.text.x = element_text(size = 12), # Customize X-axis text
        axis.text.y = element_text(size = 12),
        axis.line = element_line(color = "black", linewidth = 0.8),  # 设置轴线的颜色和大小
        axis.line.x.bottom = element_line(color = "black", linewidth = 0.8),  # 设置底部 X 轴线条
        axis.line.y.left = element_line(color = "black", linewidth = 0.8),
        strip.text = element_text(size = 10, face = "bold", color = "black")) +# 设置左边 Y 轴线条)+
  theme(#legend.title = element_blank(), #去掉legend.title 
    legend.key=element_rect(fill='white'), #
    legend.text = element_text(size=10), #设置legend标签的大小
    legend.key.size=unit(0.5,'cm') ) +  # 设置legend标签之间的大小
  theme(legend.position = "right")
MPRNMPR_cell_type

ggsave('D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/INVADE-seq/微生物featureplot/Corynebacterium.jpg',MPRNMPR_cell_type,width =8,height =7,limitsize = FALSE)
ggsave('D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/INVADE-seq/微生物featureplot/Corynebacterium_分组.jpg',MPRNMPR_cell_type,width =8,height =7,limitsize = FALSE)
ggsave('D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/INVADE-seq/微生物featureplot/Corynebacterium_病人分面.jpg',MPRNMPR_cell_type,width =8,height =8,limitsize = FALSE)
ggsave('D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/INVADE-seq/微生物featureplot/Corynebacterium_响应分面.jpg',MPRNMPR_cell_type,width =12,height =5,limitsize = FALSE)



####Fusobacterium----

MPRNMPR.umap.harmonyredu$Fusobacterium <- invade_16s[match(MPRNMPR.umap.harmonyredu$cell_id, invade_16s$barcode),"Fusobacterium"]
MPRNMPR.umap.harmonyredu$Fusobacterium[is.na(MPRNMPR.umap.harmonyredu$Fusobacterium)] <- 0

MPRNMPR.umap.harmonyredu_filter<-MPRNMPR.umap.harmonyredu[!MPRNMPR.umap.harmonyredu$Fusobacterium=="0",]
max(MPRNMPR.umap.harmonyredu_filter$Fusobacterium)
min(MPRNMPR.umap.harmonyredu_filter$Fusobacterium)
MPRNMPR.umap.harmonyredu$MPR_response2<- factor(MPRNMPR.umap.harmonyredu$MPR_response2,levels=c("Pre_NMPR","Pre_MPR","Post_MPR"))
MPRNMPR.umap.harmonyredu_filter$MPR_patient_group<-factor(MPRNMPR.umap.harmonyredu_filter$MPR_patient_group,
                                                          levels=c("Pre_NMPR-P3","Pre_NMPR-P5","Pre_NMPR-P6","Pre_NMPR-P10",
                                                                   "Pre_MPR-P1","Pre_MPR-P2","Pre_MPR-P4","Pre_MPR-P7",
                                                                   "Pre_MPR-P8","Pre_MPR-P9","Post_MPR-P1", "Post_MPR-P4",
                                                                   "Post_MPR-P7","Post_MPR-P8"))
MPRNMPR_cell_type <- ggplot()+ 
  geom_point(data=MPRNMPR.umap.harmonyredu,aes(x=umapharmony_1, y = umapharmony_2),
             size = 0.5, alpha =1,fill="#c1cdc1",color="#c1cdc1")+
  theme_bw()+
  geom_point(data=MPRNMPR.umap.harmonyredu_filter,aes(x= umapharmony_1, y = umapharmony_2,color = Fusobacterium),
             size =0.7, alpha =1)+
  #scale_colour_gradient2(low = "yellow",mid="blue", high = "red")+
  scale_color_gradientn(colors = c("blue","yellow", "red"), values = c(0,0.3, 0.5, 1), limits = c(1, 13),breaks = c(1, seq(2, 13, by = 3)))+
  xlab("UMAP_1")+
  ylab("UMAP_2")+
  ggtitle("Fusobacterium")+
  labs(color = "UMI counts") +
  #facet_wrap(~ MPR_patient_group) +
  theme(panel.grid.major = element_blank(), #主网格线
        panel.grid.minor = element_blank(), #次网格线
        panel.border = element_blank(), #边框
        panel.background = element_rect(fill = 'white'), #背景色
        plot.background=element_rect(fill="white"),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold"), # Customize X-axis title
        axis.title.y = element_text(size = 12, face = "bold"), # Customize Y-axis title
        axis.text.x = element_text(size = 12), # Customize X-axis text
        axis.text.y = element_text(size = 12),
        axis.line = element_line(color = "black", linewidth = 0.8),  # 设置轴线的颜色和大小
        axis.line.x.bottom = element_line(color = "black", linewidth = 0.8),  # 设置底部 X 轴线条
        axis.line.y.left = element_line(color = "black", linewidth = 0.8),
        strip.text = element_text(size = 10, face = "bold", color = "black")) +# 设置左边 Y 轴线条)+
  theme(#legend.title = element_blank(), #去掉legend.title 
    legend.key=element_rect(fill='white'), #
    legend.text = element_text(size=10), #设置legend标签的大小
    legend.key.size=unit(0.5,'cm') ) +  # 设置legend标签之间的大小
  theme(legend.position = "right")
MPRNMPR_cell_type

ggsave('D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/INVADE-seq/微生物featureplot/Fusobacterium.jpg',MPRNMPR_cell_type,width =8,height =7,limitsize = FALSE)
ggsave('D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/INVADE-seq/微生物featureplot/Fusobacterium_分组.jpg',MPRNMPR_cell_type,width =8,height =7,limitsize = FALSE)
ggsave('D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/INVADE-seq/微生物featureplot/Fusobacterium_病人分面.jpg',MPRNMPR_cell_type,width =8,height =8,limitsize = FALSE)
ggsave('D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/INVADE-seq/微生物featureplot/Fusobacterium_响应分面.jpg',MPRNMPR_cell_type,width =12,height =5,limitsize = FALSE)

####T 细胞拟时序分析----
BiocManager::install("monocle")
devtools::install_github("cole-trapnell-lab/monocle-release@develop")
library(monocle)
# 读入文件
T_cells_object <- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/T_cells_object.rds")
# 提取用于轨迹推断的亚型
table(T_cells_object$cd4_cd8_group)
Idents(T_cells_object) <-T_cells_object$cd4_cd8_group
#for_monocle <- subset(T_cells_object, idents = c("Malignant cell","HPC"))

#创建CellDataSet对象
exp <- as(as.matrix(T_cells_object[["RNA"]]$counts),'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = T_cells_object@meta.data)
fData <- data.frame(gene_short_name = row.names(exp), row.names = row.names(exp))
fd <- new('AnnotatedDataFrame', data = fData)

cds <- newCellDataSet(exp,
                      phenoData = pd,
                      featureData = fd,
                      expressionFamily=VGAM::negbinomial.size(),
                      lowerDetectionLimit=1)

#### 数据预处理
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

#### 细胞过滤
cds <- detectGenes(cds, min_expr = 0.1)
print(head(fData(cds)))
expressed_genes <- row.names(subset(fData(cds),num_cells_expressed >= 10)) #保存了数据集中至少10个细胞中表达的基因
expressed_genes

####关键基因筛选
pData(cds)$Cluster  <- as.factor(pData(cds)$cd4_cd8_group)
diff_test_res <- differentialGeneTest(cds[expressed_genes,],fullModelFormulaStr = "~Cluster")
head(diff_test_res)
# 选择前1000个高变基因
ordering_genes <- row.names(diff_test_res)[order(diff_test_res$qval)][1:1000]
cds <- setOrderingFilter(cds,ordering_genes = ordering_genes)
plot_ordering_genes(cds)

####基于关键基因降维
# 默认使用DDRTree方法进行降维
cds <- reduceDimension(cds, max_components = 2, method = 'DDRTree')

####沿时间轨迹排序细胞
# 根据关键基因的表达模式及其变化趋势，对细胞进行排序，并据此构建出细胞发育的轨迹。
cds <- orderCells(cds)


####单细胞微生物组成----

####umi加入胞内菌到metadata#######----

miMPRobe_data<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/INVADE-seq/16s_filter.combined.genus.umi.matrix.csv")
names(miMPRobe_data)
View(miMPRobe_data)

miMPRobe_data$barcode

metadata<- FetchData(MPRNMPR_object_miMPRobe_remove,"cell_id")
metadata$cell_id <- rownames(metadata)

View(metadata)

metadata<- left_join(x=metadata,y=miMPRobe_data,by = join_by(cell_id==barcode))

rownames(metadata)<-metadata$cell_id

MPRNMPR_object_miMPRobe_remove<- AddMetaData(MPRNMPR_object_miMPRobe_remove,metadata = metadata)

View(MPRNMPR_object_miMPRobe_remove)

table(rownames(MPRNMPR_object_miMPRobe_remove@meta.data) %in% miMPRobe_data$barcode)

scrna_metadata<- MPRNMPR_object_miMPRobe_remove@meta.data

scrna_metadata <- scrna_metadata[complete.cases(scrna_metadata), ] 

View(scrna_metadata)
######单细胞微生物组成----
library(microeco)
library(ggplot2)
library(dplyr)
names(scrna_metadata)
scrna_metadata$sampleid
genus0<- scrna_metadata[,c(5,66:578)] %>% as.data.frame()
names(genus0)
View(genus0)
genus1 <- genus0 %>%
  group_by(!!sym("sampleid")) %>%  # 将 V53 替换为你的实际列名
  summarise(across(1:513, sum, na.rm = TRUE))  %>% as.data.frame()
names(genus1)
View(genus1)
genus1$sampleid
genus2<- t(genus1)
rownames(genus2)
View(genus2)
colnames(genus2)<- genus2[1,]
genus2<- genus2[-1,] %>% as.data.frame
View(genus2)
genus2[] <- lapply(genus2, function(x) as.numeric(as.character(x)))

###分组信息
names(scrna_metadata)
sample <- scrna_metadata[,c(5,53)] %>% as.data.frame()
View(sample)
otu_cols <- colnames(genus2)
# 去除重复的 sampleid，并保持顺序
filtered_df <- sample %>%
  distinct(sampleid, .keep_all = TRUE) %>%
  arrange(match(sampleid, otu_cols))
# 确保过滤后的数据框与 otu_table 的列名顺序一致
final_df <- filtered_df[filtered_df$sampleid %in% otu_cols, ]
# 如果需要调整列名顺序
final_df <- final_df[match(otu_cols, final_df$sampleid), ]
sample1<- final_df
View(sample1)
rownames(sample1)<- sample1$sampleid
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

genus_df<-trans_abund$new(dataset=df,taxrank="genus",groupmean = "MPR_Response2",input_taxaname=c("Massilia",
                                                                                                  "Herbaspirillum",
                                                                                                  "Corynebacterium",
                                                                                                  "Sphingomonas",
                                                                                                  "Treponema",
                                                                                                  "Acinetobacter",
                                                                                                  "Achromobacter",
                                                                                                  "Fusobacterium",
                                                                                                  "Arthrobacter",
                                                                                                  "Rothia",
                                                                                                  "Lawsonella",
                                                                                                  "Rhodococcus",
                                                                                                  "Blattabacterium",
                                                                                                  "Methylobacterium",
                                                                                                  "Duganella",
                                                                                                  "Dietzia",
                                                                                                  "unclassified.Prevotellaceae",
                                                                                                  "Streptococcus",
                                                                                                  "Prevotella",
                                                                                                  "Aerococcus",
                                                                                                  "Planococcus",
                                                                                                  "Staphylococcus",
                                                                                                  "Campylobacter",
                                                                                                  "Capnocytophaga",
                                                                                                  "Bradyrhizobium",
                                                                                                  "Propionibacterium",
                                                                                                  "Lactobacillus",
                                                                                                  "Actinomyces",
                                                                                                  "Blastococcus",
                                                                                                  "Belnapia"))

genus_df<-trans_abund$new(dataset=df,taxrank="genus",groupmean = "MPR_Response2",input_taxaname=c("Achromobacter"))

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
ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/INVADE-seq/microeco/INVADE_MPR_Response2.pdf",p1,width = 9,height=7)


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
genus_df1$Taxonomy<- factor(genus_df1$Taxonomy,levels = c("Streptococcus", "Fusobacterium"))
head(genus_df1)
p3<- ggplot(genus_df1, aes(x =Taxonomy, y = Abundance,fill=MPR_Response2)) +
  geom_boxplot(outlier.shape = NA) +  
  labs(title = "", x = "", y = "Relative abundance(%)") +  # 添加标题和轴标签
  theme_minimal() +  # 使用简单主题
  scale_fill_manual(values = c("#4974a4","#4dae47","#f29600"))+
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

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/INVADE-seq/microeco/INVADE_MPR_boxplot.pdf",p3,width = 6,height=5)

genus_df<-trans_abund$new(dataset=df,taxrank="genus",input_taxaname=c("Streptococcus","Fusobacterium"),groupmean = "MPR_Response2")
genus_df<-trans_abund$new(dataset=df,taxrank="genus",ntaxa=20,groupmean = "MPR_Response2")
genus_df<-trans_abund$new(dataset=df,taxrank="genus",ntaxa=20)
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

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/INVADE-seq/microeco/INVADE_MPR_heatmap_top30.pdf",p4,width =6,height=7)

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/INVADE-seq/microeco/INVADE_MPR_heatmap_top10每个患者.pdf",p4,width =8,height=8)


####(aMiAD) scores----

install.packages("devtools")
source("https://bioconductor.org/biocLite.R")
biocLite("phyloseq")
BiocManager::install("phyloseq")
install.packages("picante")

install.packages("entropart")

install.packages("vegan")

library(devtools)
install_github("hk1785/aMiAD", force=T)

library(phyloseq)
#URL: https://joey711.github.io/phyloseq/

library(aMiAD)
library(phyloseq)
library(picante)
library(entropart)
library(vegan)

set.seed(100)
data(sim.biom)
rare.biom <- rarefy_even_depth(sim.biom, rngseed=TRUE)
  
Alpha.Diversity(phylo, metrics=c("Observed","Shannon","Simpson","PD","PE","PQE"), Normalize=TRUE)
Alpha.Diversity(sim.biom, Normalize=FALSE)

####读取微生物表格----
metadata<- MPRNMPR_object_miMPRobe_remove@meta.data
head(metadata)
metadata$cell_id
names(metadata)
invade <-read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/INVADE-seq/16s_filter.combined.genus.umi.matrix.csv")
#invade <-read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/INVADE-seq/16s_filter.combined.genus.read.matrix.csv")
invade$barcode <- gsub("-", "_", invade$barcode)

merged_df <- merge(invade, metadata, by.x="barcode",by.y="cell_id")
merged_df$barcode
names(merged_df)
library(dplyr)
result <- merged_df %>%
  group_by(orig.ident) %>%
  summarise(across(2:514, sum, na.rm = TRUE))

col_range <- 548:578

names(result)
result <- as.data.frame(result)
rownames(result)<-result$orig.ident
result<- result[,-1]
View(result)

genus_names <- colnames(result)
tax_table_data <- data.frame(genus = genus_names)
rownames(tax_table_data)<-tax_table_data$genus


####计算betadiversity----

library(vegan)
library(ape)
otu <- result
rownames(otu)
abundnce_dist<- vegdist(otu, method="bray")
pcoa <- cmdscale(abundnce_dist, k=3, eig=T)#高维度数据降低到低维度数据的操作
pcoa_points <- as.data.frame(pcoa$points)#提取三个维度的PCoA值
sum_eig <- sum(pcoa$eig)#提取特征值作为
eig_percent <- round(pcoa$eig/sum_eig*100,1)
colnames(pcoa_points) <- paste0("PCoA", 1:2:3)
pcoa_result <- cbind(pcoa_points,otu)#PCoA值和OTU合并
View(pcoa_result)
group <- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/MPR_NMPR.csv", header=T)
group$orig.ident
rownames(otu)
MPR_response2 <- group$MPR_Response2#提取Moisture分组信息

div <- adonis2(otu ~ MPR_response2, data = group, permutations = 999, method="bray")#adonis
adonis <- paste0("adonis R2: ",round(div$R2,2), "; P-value: ", div$`Pr(>F)`)#提取p值，组成一个字符串

group$MPR_Response2<-factor(group$MPR_Response2,levels = c('Pre_NMPR','Pre_MPR','Post_MPR'))


p<-ggplot(pcoa_result, aes(x=PCoA1, y=PCoA2, color=group$MPR_Response2)) +
  scale_color_manual(values =c("#4974a4","#4dae47","#f29600"))+
  labs(x=paste("PCoA 1 (", eig_percent[1], "%)", sep=""),
       y=paste("PCoA 2 (", eig_percent[2], "%)", sep=""),
       title=adonis) +#x，y轴名以及标题
  theme_bw() +
  theme(
    panel.grid=element_blank(),#去除网格线
    panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
    #legend.justification = c("right", "top"),
    legend.position = c(0.85, 0.85),
    text=element_text(size=16,  family=""),
    plot.title = element_text(hjust = 0.5,size=20),
    axis.text.x = element_text(vjust = 0.5, hjust=0.5,size=18,face="bold",family=""),
    axis.text.y = element_text(vjust = 0.5, hjust=0.5,size=18,face="bold",family=""),
    axis.title.y = element_text(angle = 90, vjust = 0.5, hjust=0.5,size=18,face = "bold"),
    axis.title.x = element_text(angle = 0, vjust = 0.5, hjust=0.5,size=18,face = "bold"),
  )+#调整x，y轴，标题的位置和字体
  geom_point(size=4)+
  guides(colour=guide_legend(title=NULL))+
  stat_ellipse(aes(fill=MPR_response2), linetype = 2,lwd=1, show.legend = F) #代替绘制置信椭圆

p

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/INVADE-seq/betadiveristy/invade_MPR_beta_diversity2.pdf",p,height =5,width = 6)


####计算alphadiversity

library(vegan)
library(picante)

otu <- result


alpha_diversity <- function(x, tree = NULL) {
  observed_species <- estimateR(x)[1, ]
  Chao1 <- estimateR(x)[2, ]
  ACE <- estimateR(x)[4, ]
  Shannon <- diversity(x, index = 'shannon',base = 2)
  Simpson <- diversity(x, index = 'simpson')    #注意，这里是Gini-Simpson 指数
  goods_coverage <- 1 - rowSums(x == 1) / rowSums(x)
  
  #保留四位小数
  Shannon <- sprintf("%0.4f", Shannon)
  Simpson <- sprintf("%0.4f", Simpson)
  goods_coverage <- sprintf("%0.4f", goods_coverage)
  
  
  result <- data.frame(observed_species,Chao1, Shannon, Simpson)
  
  if (!is.null(tree)) {
    PD_whole_tree <- pd(x, tree, include.root = FALSE)[1]
    names(PD_whole_tree) <- 'PD_whole_tree'
    result <- cbind(result, PD_whole_tree)
    
    result <- data.frame(observed_species, ACE,Chao1, Shannon, Simpson,
                         PD_whole_tree ,goods_coverage)
  }
  
  
  result
}
alpha <- alpha_diversity (otu)
# alpha1 <- alpha1 %>% 
#   mutate(across(everything(), as.numeric)) %>% as.data.frame()
type(alpha1)


write.csv(alpha,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/INVADE-seq/aMiAD_score/alpha_sample.csv")

alpha1<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/INVADE-seq/aMiAD_score/alpha11.csv",header = T)
alpha1<- sapply(alpha1, as.numeric)

type(alpha1)

alpha<- alpha %>% as.data.frame()
data(sim.biom)
alpha <- Alpha.Diversity(sim.biom)
write.csv(alpha,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/INVADE-seq/aMiAD_score/alpha.csv")
type(alpha)

col_range <- 548:578
df_subset <- merged_df %>%
  select(orig.ident, all_of(col_range))

# 根据 orig.ident 分组并去除重复行，保留每组的第一行
sample_data <-df_subset

otu_order <- rownames(result)

# 使用 otu_order 对 sample_data 进行排序
sample_data <- sample_data %>%
  slice(match(otu_order, orig.ident))

y=progeny_scores_mean[["Androgen"]]
y <- sample_data(sim.biom)$y
aMiAD_score <- aMiAD(alpha1, y, cov=NULL, model=c("gaussian","binomial"), n.perm=5000)


#fit <- aMiAD(alpha, y, cov=cbind(x1,x2), model="binomial")

aMiAD.plot(aMiAD_score, filename="D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/INVADE-seq/aMiAD_score/Figure1.pdf", fig.title="")
?aMiAD

####alpha箱装图----
sample_data$orig.ident
library(reshape2)
alpha1$MPR_response<-sample_data$MPR_Response2
dat <- melt(alpha1,id.vars = c(7),variable.name = "Alpha")


library(rstatix)
library(ggpubr)
dat$value<- as.numeric(dat$value)
df_p_val1 <- dat %>%group_by(Alpha)%>%
  wilcox_test(value~ MPR_response) %>%
  adjust_pvalue(p.col = "p", method = "bonferroni") %>%
  add_significance(p.col = "p.adj") %>%
  add_xy_position(x = "MPR_response", dodge = 0.8)

response_color <- c("Pre_NMPR"="#4974a4",
                    "Pre_MPR"="#4dae47",
                    "Post_MPR"="#f29600")

levels(dat$MPR_response)
dat$MPR_response <- factor(dat$MPR_response,levels = c("Pre_NMPR","Pre_MPR","Post_MPR"))

PMCMR_compare1 <- function(data,group,compare,value){
  library(multcompView)
  library(multcomp)
  library(PMCMRplus)
  #library(PMCMR)##多组两两比较函数用到的R包
  a <- data.frame(stringsAsFactors = F)
  type <- unique(data[,group])
  for (i in type)
  {
    g1=compare
    sub_dat <- data[data[,group]==i,]
    names(sub_dat)[names(sub_dat)==compare] <- 'g1'
    names(sub_dat)[names(sub_dat)==value] <- 'value'
    sub_dat$g1 <- factor(sub_dat$g1)
    options(warn = -1)
    
    k <- PMCMRplus::kwAllPairsNemenyiTest(value ~ g1,data=sub_dat)
    n <- as.data.frame(k$p.value)
    h <- n %>%
      mutate(compare=rownames(n)) %>%
      gather(group,p,-compare,na.rm = TRUE) %>%
      unite(compare,group,col="G",sep="-")
    dif <- h$p
    names(dif) <- h$G
    difL <- multcompLetters(dif)
    K.labels <- data.frame(difL['Letters'], stringsAsFactors = FALSE)
    K.labels$compare = rownames(K.labels)
    K.labels$type <- i
    
    mean_sd <- merge(aggregate(sub_dat[['value']],by=list(sub_dat[,'g1']),FUN=sd),
                     aggregate(sub_dat[['value']],by=list(sub_dat[,'g1']),FUN=mean),by="Group.1"
    )
    names(mean_sd) <- c('compare','std','mean')
    a <- rbind(a,merge(mean_sd,K.labels,by='compare'))
  }
  names(a) <- c(compare,'std','mean','Letters',group)
  return(a)
}

df1 <- PMCMR_compare1(dat,'Alpha','MPR_response','value')

names(dat)


alpha_plot<- ggplot(dat, aes(x = MPR_response, y = value, fill = MPR_response)) + 
  stat_boxplot(geom = "errorbar", linewidth=1)+
  geom_boxplot(linewidth=1)+
  scale_fill_manual(values = response_color) +
  geom_text(data=df1,aes(x=MPR_response,y=mean+1.5*std,label=Letters))+
  labs(title = "Box plot of alpha diversity",
       x = "",
       y = "Alpha diversity value") +
  facet_wrap(.~Alpha, scales = "free_y",nrow = 1)+
  theme_minimal()+
  theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),units=,"cm"),
        plot.title = element_text(hjust = 0.5, size = 13),
        panel.background = element_blank(),
        plot.background = element_blank(),
        strip.text = element_text(size=12),
        axis.title.y = element_text(size = 14),
        axis.line = element_line(color = "black",size = 0.4),
        axis.text.y = element_text(color="black",size=12),
        axis.text.x = element_text(size=12,colour ="black",angle = 30,vjust=0.5,hjust=0.6),
        legend.position = "none",
        panel.spacing = unit(0,"lines"),
        axis.ticks.x=element_line(color="black",size=0.5,lineend = 1),
        axis.line.y=element_line(linetype=1,color="black",size=0.7),
        axis.ticks.y=element_line(color="black",size=0.5,lineend = 1),
        axis.line.x=element_line(linetype=1,color="black",size=0.7))+
  coord_cartesian()

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/INVADE-seq/aMiAD_score/invade_seq_alpha_diversity.pdf",alpha_plot,width=12,height=4)

+++++++++++++++++++++++++++++++++++
  ++++++++++++++++++++++++++++++++++++++
####progeny计算----  
  
library(Seurat)
devtools::install_github('immunogenomics/presto')
library(presto)
Idents(MPRNMPR_object_miMPRobe_remove)<-"MPR_Response2"
group.markers <- FindAllMarkers(MPRNMPR_object_miMPRobe_remove, only.pos = TRUE, min.pct = 0.25, 
                               logfc.threshold = 0.25, verbose = FALSE)

BiocManager::install("progeny")
library(progeny)
library(dplyr)
library(Seurat)
library(ggplot2)
library(tidyr)
library(readr)
library(pheatmap)
library(tibble)


pbmc <- progeny(MPRNMPR_object_miMPRobe_remove, scale=FALSE, organism="Human", top=500, perm=1, 
                return_assay = TRUE)

pbmc@assays$progeny

pbmc <- Seurat::ScaleData(pbmc, assay = "progeny") 

progeny_scores_df <- t(pbmc@assays[["progeny"]]@scale.data)
write.csv(progeny_scores_df,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/INVADE-seq/aMiAD_score/progeny_scores_df")
head(progeny_scores_df)
progeny_scores_df <- t(pbmc@assays[["progeny"]]@scale.data) %>% cbind("orig.ident"=pbmc$orig.ident,"cell_type_new"=pbmc$cell_type_new,"MPR_Response2"=pbmc$MPR_Response2)%>% as.data.frame()
progeny_scores_df$Cell<- rownames(progeny_scores_df)
write.csv(progeny_scores_df,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/INVADE-seq/aMiAD_score/progeny_scores_df.csv")

progeny_scores_df[1:14] <- progeny_scores_df[1:14] %>% 
  mutate(across(everything(), as.numeric))

progeny_scores_mean <- progeny_scores_df %>%
  group_by(orig.ident) %>%
  summarise(across(1:14, mean, na.rm = TRUE))

progeny_scores_mean$orig.ident2 <- sample_data$orig.ident
names(sample_data)
progeny_scores_mean$MPR_Response2 <- sample_data$MPR_Response2
View(progeny_scores_mean)

write.csv(progeny_scores_mean,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/INVADE-seq/aMiAD_score/progeny_scores_df_mean.csv")

library(aMiAD)
library(phyloseq)
library(picante)
library(entropart)
library(vegan)

progeny_scores_mean1<- progeny_scores_mean[,-1]
progeny_scores_mean1 <-progeny_scores_mean1 %>% 
  mutate(across(everything(), as.numeric))
names(progeny_scores_mean)

aMiAD_scores <- list()

columns_to_analyze <- c("Androgen","EGFR","Estrogen","Hypoxia","JAK-STAT","MAPK",      
                        "NFkB","p53","PI3K","TGFb","TNFa","Trail","VEGF",      
                        "WNT")

for (col_name in columns_to_analyze) {
  y <- progeny_scores_mean[[col_name]]  # 提取当前列的数据
  # 运行 aMiAD 函数
  aMiAD_scores[[col_name]] <- aMiAD(alpha1, y, cov = NULL, model = c("gaussian", "binomial"), n.perm = 5000)
  # 提取结果
}

all_item_by_item <- data.frame()
all_a_miad <- data.frame()

for (col_name in columns_to_analyze) {
  # 提取当前列的结果
  item_by_item_out <- aMiAD_scores[[col_name]]$ItembyItem.out
  a_miad_out <- aMiAD_scores[[col_name]]$aMiAD.out
  
  # 将结果转换为数据框
  df_item_by_item <- as.data.frame(item_by_item_out)
  df_a_miad <- as.data.frame(a_miad_out)
  
  # 添加列名标识
  df_item_by_item$col_name <- col_name
  df_a_miad$col_name <- col_name
  
  # 将当前迭代的结果添加到总的合并数据框中
  all_item_by_item <- bind_rows(all_item_by_item, df_item_by_item)
  all_a_miad <- bind_rows(all_a_miad, df_a_miad)
}


# 打印结果
print(all_item_by_item)
print(all_a_miad)

write.csv(all_item_by_item,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/INVADE-seq/aMiAD_score/aMiAD_progeny_all_item_by_item2.csv")

write.csv(all_a_miad,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/INVADE-seq/aMiAD_score/aMiAD_progeny_all_a_miad2.csv")
rownames(all_item_by_item)
all_item_by_item_filter <- all_item_by_item[grep("observed_species", rownames(all_item_by_item)), ]
names(all_item_by_item)
all_item_by_item_filter$col_name<- factor(all_item_by_item_filter$col_name,levels=rev(c("Trail","Androgen","p53","WNT", "TGFb" ,"Estrogen","Hypoxia","NFkB","TNFa", "JAK-STAT","VEGF","PI3K","EGFR","MAPK")))
names(all_a_miad)
all_a_miad$aMiDivES <- rownames(all_a_miad)
all_a_miad_filter <- all_a_miad %>% filter(str_detect(aMiDivES, "aMiDivES"))
all_a_miad_filter$a_miad_out <- as.numeric(all_a_miad_filter$a_miad_out)
View(wide_df)
all_a_miad_filter$col_name<- factor(all_a_miad_filter$col_name,levels=rev(c("Trail","Androgen","p53","WNT", "TGFb" ,"Estrogen","Hypoxia","NFkB","TNFa", "JAK-STAT","VEGF","PI3K","EGFR","MAPK")))
p<- ggplot(all_a_miad_filter, aes(x = a_miad_out, y = col_name, fill = a_miad_out)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient2(low = "blue", mid = "grey", high = "red", midpoint = 0) +
  labs(title = "",
       x = "aMiDivES",
       y = "PROGENy pathway",
       fill = "aMiDivES") +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = NA, color = "black"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        strip.text.y.left = element_text(angle = 0),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(angle = 90, size=8, hjust = 1, vjust = 0.5, color = "black"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7)
  ) 
p
ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/INVADE-seq/aMiAD_score/PROGENy_pathway_aMiDivES.pdf",p,width =3,height=3)

# Idents(pbmc)<-"cell_type_new"
# CellsClusters <- data.frame(Cell = names(Idents(pbmc)), 
#                             CellType = as.character(Idents(pbmc)),
#                             stringsAsFactors = FALSE)
# 
# progeny_scores_df <- inner_join(progeny_scores_df, CellsClusters)
# head(progeny_scores_df)
# 
# pbmc@meta.data$orig.ident
# Idents(pbmc)<-"orig.ident"
# orig.identClusters <- data.frame(Cell = names(Idents(pbmc)), 
#                             sampleid = as.character(Idents(pbmc)),
#                             stringsAsFactors = FALSE)
# 
# progeny_scores_df <- inner_join(progeny_scores_df, orig.identClusters)
# head(progeny_scores_df)
# 
# Idents(pbmc)<-"MPR_Response2"
# MPR_Response2Clusters <- data.frame(Cell = names(Idents(pbmc)), 
#                                     MPR_Response2 = as.character(Idents(pbmc)),
#                                  stringsAsFactors = FALSE)
# 
# progeny_scores_df <- inner_join(progeny_scores_df, MPR_Response2Clusters)
# 
# head(progeny_scores_df)

## We summarize the Progeny scores by cellpopulation
progeny_scores_df_long<- progeny_scores_df %>%
  pivot_longer(
    cols = -c(orig.ident,cell_type_new, MPR_Response2, CellType, Cell),  # 保留这些列不变
    names_to = "Measurement",        # 新的列名，用于存储原列的列名
    values_to = "Value"              # 新的列名，用于存储原列的值
  )

summarized_progeny_scores1 <- progeny_scores_df_long %>% 
  group_by(Measurement, MPR_Response2) %>%
  summarise(avg = mean(Value), std = sd(Value))
dim(summarized_progeny_scores1)

head(summarized_progeny_scores1)
wide_data <- summarized_progeny_scores1 %>%
  pivot_wider(names_from = Pathway,
              values_from = c(avg, std))
write.csv(wide_data,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/INVADE-seq/aMiAD_score/summarized_progeny_scores1_sampleid.csv")

summarized_progeny_scores_df <- summarized_progeny_scores1 %>%
  dplyr::select(-std) %>%   
  spread(Measurement, avg) %>%
  data.frame(row.names = 1, check.names = FALSE, stringsAsFactors = FALSE) 

summarized_progeny_scores_df  <- summarized_progeny_scores_df [c("Pre_NMPR", "Pre_MPR", "Post_MPR"), ]
summarized_progeny_scores_df_t <- t(summarized_progeny_scores_df)
summarized_progeny_scores_df_t  <- summarized_progeny_scores_df_t[,c("Pre_NMPR", "Pre_MPR", "Post_MPR")]
colnames(summarized_progeny_scores_df_t) <- factor(colnames(summarized_progeny_scores_df_t),levels=c("Pre_NMPR", "Pre_MPR", "Post_MPR"))


paletteLength = 100
myColor = colorRampPalette(c("Darkblue", "white","red"))(paletteLength)


progenyBreaks = c(seq(min(summarized_progeny_scores_df), 0, 
                      length.out=ceiling(paletteLength/2) + 1),
                  seq(max(summarized_progeny_scores_df)/paletteLength, 
                      max(summarized_progeny_scores_df), 
                      length.out=floor(paletteLength/2)))
progeny_hmap = pheatmap(summarized_progeny_scores_df_t,fontsize=12, 
                        fontsize_row = 12, 
                        color=myColor, breaks = progenyBreaks, 
                        main = "PROGENy (500)", angle_col = 45,
                        treeheight_col = 0,  border_color = NA)

rownames(summarized_progeny_scores_df_t)
ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/INVADE-seq/aMiAD_score/progeny_hmap_group.pdf",progeny_hmap,width = 5,height=6)


####KEGG、GO、Reactome、Do、MSigDB进行GSEA富集----

library(org.Hs.eg.db)#物种注释包(Homo sapiens)
library(clusterProfiler)#富集分析

Idents(MPRNMPR_object_miMPRobe_remove)<-"MPR_Response2"

group.markers <- FindAllMarkers(MPRNMPR_object_miMPRobe_remove, only.pos = TRUE, min.pct = 0.25, 
                                logfc.threshold = 0.25, verbose = FALSE)


symbol <- rownames(group.markers)
entrez <- bitr(symbol,
               fromType = "SYMBOL",#现有的ID类型
               toType = "ENTREZID",#需转换的ID类型
               OrgDb = "org.Hs.eg.db")
head(entrez)

names(group.markers)
genelist <- group.markers$avg_log2FC
names(genelist) <- rownames(group.markers)
head(genelist)

genelist <- genelist[names(genelist) %in% entrez[,1]]
names(genelist) <- entrez[match(names(genelist),entrez[,1]),2]
length(genelist)#过滤后剩16161个基因
head(genelist)

#将genelist按照log2FC值从高到低进行排序：
genelist <- sort(genelist,decreasing = T)
head(genelist)

KEGG_ges <- gseKEGG(
  geneList = genelist,
  organism = "hsa",#不同物种选择可官方查询：https://www.genome.jp/kegg/catalog/org_list.html
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose = FALSE,
  eps = 0
)

#提取结果表：
KEGG_ges_result <- KEGG_ges@result

#显示前10列（第十一列为富集到该基因集的核心基因列表，这里过长不展示）：
KEGG_ges_result[1:3,1:10]


#GSEA可视化：
##山峦图：
library(ggridges)
library(ggplot2)
library(enrichplot)

p <- ridgeplot(KEGG_ges,
               showCategory = 15,
               fill = "p.adjust",
               decreasing  = T)
p

##GSEA_GO富集分析：
GO_ges <- gseGO(geneList = genelist,
                OrgDb = org.Hs.eg.db,
                ont = "CC", #one of "BP", "MF", and "CC" subontologies, or "ALL" for all three.
                minGSSize = 10,
                maxGSSize = 500,
                pvalueCutoff = 0.05,
                eps = 0,
                verbose = FALSE)
GO_ges[1:3,1:10] #展示同样省略最后一列


p_pre_nmpr <- ridgeplot(GO_ges,
               showCategory = 15,
               fill = "p.adjust",
               decreasing  = T)
p_pre_nmpr



####分组GSEA_Reactome富集分析----
####pre_NMPR
group.markers <- FindAllMarkers(MPRNMPR_object_miMPRobe_remove, only.pos = TRUE, min.pct = 0.25, 
                                logfc.threshold = 0.25, verbose = FALSE)

names(group.markers)
group.marker_pre_nmpr<- group.markers %>% filter(cluster %in% c("Pre_NMPR"))

unique(group.marker_pre_nmpr$cluster)

symbol <- rownames(group.marker_pre_nmpr)
entrez <- bitr(symbol,
               fromType = "SYMBOL",#现有的ID类型
               toType = "ENTREZID",#需转换的ID类型
               OrgDb = "org.Hs.eg.db")
head(entrez)

genelist <- group.marker_pre_nmpr$avg_log2FC
names(genelist) <- rownames(group.marker_pre_nmpr)
head(genelist)

genelist <- genelist[names(genelist) %in% entrez[,1]]
names(genelist) <- entrez[match(names(genelist),entrez[,1]),2]
length(genelist)#过滤后剩16161个基因
head(genelist)

#将genelist按照log2FC值从高到低进行排序：
genelist <- sort(genelist,decreasing = T)
head(genelist)

BiocManager::install("ReactomePA")
library(ReactomePA)
ret_ges <- gsePathway(genelist,
                      organism = "human",
                      minGSSize = 10,
                      maxGSSize = 500,
                      pvalueCutoff = 0.05,
                      pAdjustMethod = "BH",
                      verbose = FALSE,
                      eps = 0)
ret_ges[1:3,1:10]

p <- ridgeplot(ret_ges,
               showCategory = 15,
               fill = "p.adjust",
               decreasing  = T)
p


####pre_MPR
group.markers <- FindAllMarkers(MPRNMPR_object_miMPRobe_remove, only.pos = TRUE, min.pct = 0.25, 
                                logfc.threshold = 0.25, verbose = FALSE)

names(group.markers)
group.marker_pre_mpr<- group.markers %>% filter(cluster %in% c("Pre_MPR"))

unique(group.marker_pre_mpr$cluster)

symbol <- rownames(group.marker_pre_mpr)
entrez <- bitr(symbol,
               fromType = "SYMBOL",#现有的ID类型
               toType = "ENTREZID",#需转换的ID类型
               OrgDb = "org.Hs.eg.db")
head(entrez)

genelist <- group.marker_pre_mpr$avg_log2FC
names(genelist) <- rownames(group.marker_pre_mpr)
head(genelist)

genelist <- genelist[names(genelist) %in% entrez[,1]]
names(genelist) <- entrez[match(names(genelist),entrez[,1]),2]
length(genelist)#过滤后剩16161个基因
head(genelist)

#将genelist按照log2FC值从高到低进行排序：
genelist <- sort(genelist,decreasing = T)
head(genelist)

BiocManager::install("ReactomePA")
library(ReactomePA)
ret_ges <- gsePathway(genelist,
                      organism = "human",
                      minGSSize = 10,
                      maxGSSize = 500,
                      pvalueCutoff = 0.05,
                      pAdjustMethod = "BH",
                      verbose = FALSE,
                      eps = 0)
ret_ges[1:3,1:10]

p_pre_mpr <- ridgeplot(ret_ges,
               showCategory = 15,
               fill = "p.adjust",
               decreasing  = T)
p_pre_mpr


####post_MPR
group.markers <- FindAllMarkers(MPRNMPR_object_miMPRobe_remove, only.pos = TRUE, min.pct = 0.25, 
                                logfc.threshold = 0.25, verbose = FALSE)

names(group.markers)
group.marker_post_mpr<- group.markers %>% filter(cluster %in% c("Post_MPR"))

unique(group.marker_post_mpr$cluster)

symbol <- rownames(group.marker_post_mpr)
entrez <- bitr(symbol,
               fromType = "SYMBOL",#现有的ID类型
               toType = "ENTREZID",#需转换的ID类型
               OrgDb = "org.Hs.eg.db")
head(entrez)

genelist <- group.marker_post_mpr$avg_log2FC
names(genelist) <- rownames(group.marker_post_mpr)
head(genelist)

genelist <- genelist[names(genelist) %in% entrez[,1]]
names(genelist) <- entrez[match(names(genelist),entrez[,1]),2]
length(genelist)#过滤后剩16161个基因
head(genelist)

#将genelist按照log2FC值从高到低进行排序：
genelist <- sort(genelist,decreasing = T)
head(genelist)

BiocManager::install("ReactomePA")
library(ReactomePA)
ret_ges <- gsePathway(genelist,
                      organism = "human",
                      minGSSize = 10,
                      maxGSSize = 500,
                      pvalueCutoff = 0.05,
                      pAdjustMethod = "BH",
                      verbose = FALSE,
                      eps = 0)
ret_ges[1:3,1:10]

p_post_mpr <- ridgeplot(ret_ges,
                       showCategory = 15,
                       fill = "p.adjust",
                       decreasing  = T)
p_post_mpr



###p53: 




##GSEA_DO(Disease Ontology)富集分析:
library(DOSE)
DO_ges <- gseDO(genelist,
                minGSSize = 10,
                maxGSSize = 500,
                pvalueCutoff = 0.05,
                pAdjustMethod = "BH",
                verbose = FALSE,
                eps = 0)
DO_ges[1:3,1:10]

p <- ridgeplot(DO_ges,
               showCategory = 15,
               fill = "p.adjust",
               decreasing  = T)
p


##GSEA_MSigDB富集分析：
#包的下载和载入：
devtools::install_github("ToledoEM/msigdf")
library(msigdf)
#提取C2注释(human)：
library(dplyr)
c2 <- msigdf.human %>%
  filter(category_code == "c2") %>% select(geneset, symbol) %>% as.data.frame

head(c2)

#这里genelist需要的是symbol：
#重新准备genelist文件：
genelist2 <-group.markers$avg_log2FC 
names(genelist2) <- rownames(group.markers)

#genelist过滤(ID转换中丢失的部分基因)：
genelist2 <- genelist2[names(genelist2) %in% entrez[,1]]

#将genelist按照log2FC值从高到低进行排序：
genelist2 <- sort(genelist2,decreasing = T)
head(genelist2)

c2_ges <- GSEA(genelist2,
               TERM2GENE = c2,
               minGSSize = 10,
               maxGSSize = 500,
               pvalueCutoff = 0.05,
               pAdjustMethod = "BH",
               verbose = FALSE,
               eps = 0)
c2_ges_result <- c2_ges@result
View(c2_ges_result)


#在富集结果表中将entrez转换为symbol：
#转换前：
KEGG_ges@result$core_enrichment[1]
#ID转换：
KEGG_ges_set <- setReadable(KEGG_ges,
                            OrgDb = org.Hs.eg.db,
                            keyType="ENTREZID")

#转换后：
KEGG_ges_set@result$core_enrichment[1]




###获取了TIDE_signatures(Signatures of T cell dysfunction and exclusion predict cancer immunotherapy response)----

#不同病人基因表达量
avg_expression <- AverageExpression(MPRNMPR_object_miMPRobe_remove, group.by = "orig.ident", assays = "RNA")

# 获取平均表达矩阵
avg_expr_matrix <- avg_expression$RNA
avg_expr_matrix

normalized_matrix <- scale(avg_expr_matrix)


write.table(normalized_matrix,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/INVADE-seq/aMiAD_score/dysfunction_interaction_test-1.0 (1)/每个病人的基因表达量.tsv",sep = "\t")


metadata<- MPRNMPR_object_miMPRobe_remove@meta.data
head(metadata)
metadata$cell_id
names(metadata)
invade <-read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/INVADE-seq/16s_filter.combined.genus.umi.matrix.csv")
#invade <-read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/INVADE-seq/16s_filter.combined.genus.read.matrix.csv")
invade$barcode <- gsub("-", "_", invade$barcode)

merged_df <- merge(invade, metadata, by.x="barcode",by.y="cell_id")
merged_df$barcode
names(merged_df)
library(dplyr)
result <- merged_df %>%
  group_by(orig.ident) %>%
  summarise(across(2:514, sum, na.rm = TRUE))

col_range <- 548:578

names(result)
result <- as.data.frame(result)
rownames(result)<-result$orig.ident
result<- result %>% select(-1)
View(result)

genus_names <- colnames(result)
tax_table_data <- data.frame(genus = genus_names)
rownames(tax_table_data)<-tax_table_data$genus

library(vegan)
library(picante)
otu <- result
alpha_diversity <- function(x, tree = NULL) {
  observed_species <- estimateR(x)[1, ]
  Chao1 <- estimateR(x)[2, ]
  ACE <- estimateR(x)[4, ]
  Shannon <- diversity(x, index = 'shannon',base = 2)
  Simpson <- diversity(x, index = 'simpson')    #注意，这里是Gini-Simpson 指数
  goods_coverage <- 1 - rowSums(x == 1) / rowSums(x)
  
  #保留四位小数
  Shannon <- sprintf("%0.4f", Shannon)
  Simpson <- sprintf("%0.4f", Simpson)
  goods_coverage <- sprintf("%0.4f", goods_coverage)
  
  
  result <- data.frame(observed_species,Chao1, Shannon, Simpson)
  
  if (!is.null(tree)) {
    PD_whole_tree <- pd(x, tree, include.root = FALSE)[1]
    names(PD_whole_tree) <- 'PD_whole_tree'
    result <- cbind(result, PD_whole_tree)
    
    result <- data.frame(observed_species, ACE,Chao1, Shannon, Simpson,
                         PD_whole_tree ,goods_coverage)
  }
  
  
  result
}
alpha <- alpha_diversity (otu)

order_alpha<-rownames(alpha)
type(alpha1)
alpha1<- sapply(alpha, as.numeric)

aMiAD_scores <- list()

TIDE_signatures <- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/INVADE-seq/aMiAD_score/dysfunction_interaction_test-1.0 (1)/TIDE_signatures2.csv",row.names = 1)

names(TIDE_signatures)=c("TIDE","IFNG","MSI-SCORE","MERCK18","CD274","CD8","CTL",    
                         "Dysfunction","Exclusion" , "MDSC" ,"CAF" ,"TAM-M2")

columns_to_analyze <- c("TIDE","IFNG","MSI-SCORE","MERCK18","CD274","CD8","CTL",    
                        "Dysfunction","Exclusion" , "MDSC" ,"CAF" ,"TAM-M2")

library(aMiAD)
library(phyloseq)
library(picante)
library(entropart)
library(vegan)

for (col_name in columns_to_analyze) {
  y <- TIDE_signatures[[col_name]]  # 提取当前列的数据
  # 运行 aMiAD 函数
  aMiAD_scores[[col_name]] <- aMiAD(alpha1, y, cov = NULL, model = c("gaussian", "binomial"), n.perm = 5000)
  # 提取结果
}

all_item_by_item <- data.frame()
all_a_miad <- data.frame()

for (col_name in columns_to_analyze) {
  # 提取当前列的结果
  item_by_item_out <- aMiAD_scores[[col_name]]$ItembyItem.out
  a_miad_out <- aMiAD_scores[[col_name]]$aMiAD.out
  
  # 将结果转换为数据框
  df_item_by_item <- as.data.frame(item_by_item_out)
  df_a_miad <- as.data.frame(a_miad_out)
  
  # 添加列名标识
  df_item_by_item$col_name <- col_name
  df_a_miad$col_name <- col_name
  
  # 将当前迭代的结果添加到总的合并数据框中
  all_item_by_item <- bind_rows(all_item_by_item, df_item_by_item)
  all_a_miad <- bind_rows(all_a_miad, df_a_miad)
}


# 打印结果
print(all_item_by_item)
print(all_a_miad)

write.csv(all_item_by_item,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/INVADE-seq/aMiAD_score/dysfunction_interaction_test-1.0 (1)/aMiAD_TIDE_all_item_by_item.csv")
write.csv(all_a_miad,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/INVADE-seq/aMiAD_score/dysfunction_interaction_test-1.0 (1)/aMiAD_TIDE_all_a_miad.csv")

library(stringr)
library(ggplot2)
all_a_miad$aMiDivES <- rownames(all_a_miad)
all_a_miad_filter <- all_a_miad %>% filter(str_detect(aMiDivES, "aMiDivES"))
all_a_miad_filter$a_miad_out <- as.numeric(all_a_miad_filter$a_miad_out)
all_a_miad_filter$col_name<- factor(all_a_miad_filter$col_name,levels=rev(c("TIDE","IFNG","MSI-SCORE","MERCK18","CD274","CD8","CTL","Dysfunction",
                                                                            "Exclusion","MDSC","CAF","TAM-M2")))
p<- ggplot(all_a_miad_filter, aes(x = a_miad_out, y = col_name, fill = a_miad_out)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient2(low = "blue", mid = "grey", high = "red", midpoint = 0) +
  labs(title = "",
       x = "aMiDivES",
       y = "PROGENy pathway",
       fill = "aMiDivES") +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = NA, color = "black"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        strip.text.y.left = element_text(angle = 0),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(angle = 90, size=8, hjust = 1, vjust = 0.5, color = "black"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7)
  ) 
p
ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/INVADE-seq/aMiAD_score/dysfunction_interaction_test-1.0 (1)/aMiAD_TIDE_aMiDivES.pdf",p,width =3,height=3)

library(tidyr)

TIDE_signatures <- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/INVADE-seq/aMiAD_score/dysfunction_interaction_test-1.0 (1)/TIDE_signatures2.csv",row.names = 1)

names(TIDE_signatures)=c("TIDE","IFNG","MSI-SCORE","MERCK18","CD274","CD8","CTL",    
                         "Dysfunction","Exclusion" , "MDSC" ,"CAF" ,"TAM-M2","MPR_Response2")

TIDE_signatures_long<- TIDE_signatures %>%
  pivot_longer(
    cols = -c(MPR_Response2),  # 保留这些列不变
    names_to = "Measurement",        # 新的列名，用于存储原列的列名
    values_to = "Value"              # 新的列名，用于存储原列的值
  )

TIDE_signatures_scores1 <-TIDE_signatures_long%>% 
  group_by(Measurement, MPR_Response2) %>%
  summarise(avg = mean(Value), std = sd(Value))
dim(TIDE_signatures_scores1)

head(TIDE_signatures_scores1)
wide_data <- TIDE_signatures_scores1 %>%
  pivot_wider(names_from = Measurement,
              values_from = c(avg, std))
write.csv(wide_data,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/INVADE-seq/aMiAD_score/dysfunction_interaction_test-1.0 (1)/TIDE_signatures2_heatmap数据.csv")

TIDE_signatures_scores1_df <- TIDE_signatures_scores1 %>%
  dplyr::select(-std) %>%   
  spread(Measurement, avg) %>%
  data.frame(row.names = 1, check.names = FALSE, stringsAsFactors = FALSE) 

TIDE_signatures_scores1_df  <- TIDE_signatures_scores1_df [c("Pre_NMPR", "Pre_MPR", "Post_MPR"), ]
TIDE_signatures_scores1_df_t <- t(TIDE_signatures_scores1_df)
TIDE_signatures_scores1_df_t  <- TIDE_signatures_scores1_df_t[,c("Pre_NMPR", "Pre_MPR", "Post_MPR")]
colnames(TIDE_signatures_scores1_df_t) <- factor(colnames(TIDE_signatures_scores1_df_t),levels=c("Pre_NMPR", "Pre_MPR", "Post_MPR"))



paletteLength = 100
myColor = colorRampPalette(c("Darkblue", "white","red"))(paletteLength)

myColor = colorRampPalette(c("blue", "cyan", "yellow", "red"))(paletteLength)

library(pheatmap)
progenyBreaks = c(seq(min(TIDE_signatures_scores1_df), 0, 
                      length.out=ceiling(paletteLength/2) + 1),
                  seq(max(TIDE_signatures_scores1_df)/paletteLength, 
                      max(TIDE_signatures_scores1_df), 
                      length.out=floor(paletteLength/2)))
progeny_hmap = pheatmap(TIDE_signatures_scores1_df_t,fontsize=12, 
                        fontsize_row = 12, 
                        color=myColor, breaks = progenyBreaks, 
                        main = "TIDE_signatures", angle_col = 45,
                        treeheight_col = 0,  border_color = NA)

rownames(TIDE_signatures_scores1_df_t)
ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/INVADE-seq/aMiAD_score/dysfunction_interaction_test-1.0 (1)/TIDE_signatures_hmap_group.pdf",progeny_hmap,width = 5,height=6)


####免疫细胞类型 aMiAD----

metadata<- MPRNMPR_object_miMPRobe_remove@meta.data
head(metadata)
metadata$cell_id
names(metadata)
invade <-read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/INVADE-seq/16s_filter.combined.genus.umi.matrix.csv")
#invade <-read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/INVADE-seq/16s_filter.combined.genus.read.matrix.csv")
invade$barcode <- gsub("-", "_", invade$barcode)

merged_df <- merge(invade, metadata, by.x="barcode",by.y="cell_id")
merged_df$barcode
names(merged_df)
library(dplyr)
merged_df$merge_cell_type2 <- paste0(merged_df$orig.ident,"_",merged_df$cell_type_new)
result <- merged_df %>%
  group_by(merge_cell_type2) %>%
  summarise(across(2:514, sum, na.rm = TRUE))

names(result)
#result$group<- paste(result$MPR_Response2,result$cell_type_new,sep = "_")

result <- as.data.frame(result)
rownames(result)<-result$merge_cell_type
names(result)
result<- result[, -1]

View(result)

genus_names <- colnames(result)
tax_table_data <- data.frame(genus = genus_names)
rownames(tax_table_data)<-tax_table_data$genus

####计算每个细胞alphadiversity----

library(vegan)
library(picante)

otu <- result

alpha_diversity <- function(x, tree = NULL) {
  observed_species <- estimateR(x)[1, ]
  Chao1 <- estimateR(x)[2, ]
  ACE <- estimateR(x)[4, ]
  Shannon <- diversity(x, index = 'shannon',base = 2)
  Simpson <- diversity(x, index = 'simpson')    #注意，这里是Gini-Simpson 指数
  goods_coverage <- 1 - rowSums(x == 1) / rowSums(x)
  
  #保留四位小数
  Shannon <- sprintf("%0.4f", Shannon)
  Simpson <- sprintf("%0.4f", Simpson)
  goods_coverage <- sprintf("%0.4f", goods_coverage)
  
  
  result <- data.frame(observed_species, ACE,Chao1, Shannon, Simpson, goods_coverage)
  
  if (!is.null(tree)) {
    PD_whole_tree <- pd(x, tree, include.root = FALSE)[1]
    names(PD_whole_tree) <- 'PD_whole_tree'
    result <- cbind(result, PD_whole_tree)
    
    result <- data.frame(observed_species, ACE,Chao1, Shannon, Simpson,
                         PD_whole_tree ,goods_coverage)
  }
  
  
  result
}

alpha <- alpha_diversity (otu)
View(alpha)
rownames(alpha)
new_order <- c("Pre_NMPR_T cells","Pre_NMPR_B cells","Pre_NMPR_Plasma cells","Pre_NMPR_Myeloid cells","Pre_NMPR_Mast cells","Pre_NMPR_Endothelial cells","Pre_NMPR_Fibroblasts","Pre_NMPR_Smooth muscle cells","Pre_NMPR_Pericytes","Pre_NMPR_Epithelial cells",
               
               "Pre_MPR_T cells","Pre_MPR_B cells","Pre_MPR_Plasma cells","Pre_MPR_Myeloid cells","Pre_MPR_Mast cells","Pre_MPR_Endothelial cells","Pre_MPR_Fibroblasts","Pre_MPR_Smooth muscle cells", "Pre_MPR_Pericytes","Pre_MPR_Epithelial cells",
               
               "Post_MPR_T cells","Post_MPR_B cells","Post_MPR_Plasma cells","Post_MPR_Myeloid cells","Post_MPR_Mast cells","Post_MPR_Endothelial cells","Post_MPR_Fibroblasts", "Post_MPR_Smooth muscle cells","Post_MPR_Pericytes","Post_MPR_Epithelial cells")

new_order <- new_order[new_order %in% rownames(alpha)]

alpha <- alpha[new_order, ]
rownames(alpha)

#"Post_MPR_B cells"             "Pre_MPR_B cells"              "Pre_NMPR_B cells"            
# "Post_MPR_Endothelial cells"   "Pre_MPR_Endothelial cells"    "Pre_NMPR_Endothelial cells"  
# "Post_MPR_Epithelial cells"    "Pre_MPR_Epithelial cells"     "Pre_NMPR_Epithelial cells"   
# "Post_MPR_Fibroblasts"         "Pre_MPR_Fibroblasts"          "Pre_NMPR_Fibroblasts"        
# "Post_MPR_Mast cells"          "Pre_MPR_Mast cells"           "Pre_NMPR_Mast cells"         
# "Post_MPR_Myeloid cells"       "Pre_MPR_Myeloid cells"        "Pre_NMPR_Myeloid cells"      
# "Post_MPR_Pericytes"           "Pre_MPR_Pericytes"            "Pre_NMPR_Pericytes"          
# "Post_MPR_Plasma cells"        "Pre_MPR_Plasma cells"         "Pre_NMPR_Plasma cells"       
# "Post_MPR_Smooth muscle cells" "Pre_MPR_Smooth muscle cells"  "Pre_NMPR_Smooth muscle cells"
# "Post_MPR_T cells"             "Pre_MPR_T cells"              "Pre_NMPR_T cells"   

write.csv(alpha,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/INVADE-seq/aMiAD_score/不同细胞类型alpha11.csv")

alpha1<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/INVADE-seq/aMiAD_score/不同细胞类型alpha11.csv",header = T)
alpha1<- sapply(alpha, as.numeric)
View(alpha1)
type(alpha1)

aMiAD_score <- aMiAD(alpha1, y, cov=NULL, model=c("gaussian","binomial"), n.perm=5000)


####alpha箱装图----
sample_data$orig.ident
library(reshape2)
alpha$MPR_response<-sample_data$MPR_Response2

alpha$cell_type <-sub("^[^_]*_[^_]*_[^_]*_", "", rownames(alpha))

nrow(alpha)
alpha$MPR_Response2 <-c(rep("Pre_NMPR",10),
                        rep("Pre_NMPR",9),
                        rep("Post_MPR",9),
                        rep("Pre_MPR",10),
                        rep("Pre_MPR",9),
                        rep("Pre_MPR",9),
                        rep("Post_MPR",10),
                        rep("Pre_MPR",9),
                        rep("Post_MPR",9),
                        rep("Pre_MPR",8),
                        rep("Post_MPR",9),
                        rep("Pre_MPR",4),
                        rep("Pre_NMPR",9),
                        rep("Pre_NMPR",4))
names(alpha)
alpha$merged_group<-rownames(alpha)
dat <- melt(alpha,id.vars = c(7,8,9),variable.name = "Alpha")
head(dat)
library(rstatix)
library(ggpubr)
dat$value<- as.numeric(dat$value)
df_p_val1 <- dat %>%group_by(Alpha)%>%
  wilcox_test(value~ MPR_Response2) %>%
  adjust_pvalue(p.col = "p", method = "bonferroni") %>%
  add_significance(p.col = "p.adj") %>%
  add_xy_position(x = "MPR_response", dodge = 0.8)

response_color <- c("Pre_NMPR"="#4974a4",
                    "Pre_MPR"="#4dae47",
                    "Post_MPR"="#f29600")

levels(dat$MPR_response)
dat$MPR_Response2 <- factor(dat$MPR_Response2,levels = c("Pre_NMPR","Pre_MPR","Post_MPR"))

PMCMR_compare1 <- function(data,group,compare,value){
  library(multcompView)
  library(multcomp)
  library(PMCMRplus)
  #library(PMCMR)##多组两两比较函数用到的R包
  a <- data.frame(stringsAsFactors = F)
  type <- unique(data[,group])
  for (i in type)
  {
    g1=compare
    sub_dat <- data[data[,group]==i,]
    names(sub_dat)[names(sub_dat)==compare] <- 'g1'
    names(sub_dat)[names(sub_dat)==value] <- 'value'
    sub_dat$g1 <- factor(sub_dat$g1)
    options(warn = -1)
    
    k <- PMCMRplus::kwAllPairsNemenyiTest(value ~ g1,data=sub_dat)
    n <- as.data.frame(k$p.value)
    h <- n %>%
      mutate(compare=rownames(n)) %>%
      gather(group,p,-compare,na.rm = TRUE) %>%
      unite(compare,group,col="G",sep="-")
    dif <- h$p
    names(dif) <- h$G
    difL <- multcompLetters(dif)
    K.labels <- data.frame(difL['Letters'], stringsAsFactors = FALSE)
    K.labels$compare = rownames(K.labels)
    K.labels$type <- i
    
    mean_sd <- merge(aggregate(sub_dat[['value']],by=list(sub_dat[,'g1']),FUN=sd),
                     aggregate(sub_dat[['value']],by=list(sub_dat[,'g1']),FUN=mean),by="Group.1"
    )
    names(mean_sd) <- c('compare','std','mean')
    a <- rbind(a,merge(mean_sd,K.labels,by='compare'))
  }
  names(a) <- c(compare,'std','mean','Letters',group)
  return(a)
}



names(dat)
unique(dat$Alpha)
dat1<- dat%>% filter(Alpha %in% c("Simpson"))
dat1$cell_type<-factor(dat1$cell_type,levels= c('T cells', 'B cells','Plasma cells', 'Myeloid cells',  
                                                'Mast cells',"Endothelial cells","Fibroblasts","Smooth muscle cells",
                                                "Pericytes","Epithelial cells"))

df1 <- PMCMR_compare1(dat1,'Alpha','cell_type','value')
dat1$MPR_Response2 <- factor(dat1$MPR_Response2,levels = c("Pre_NMPR","Pre_MPR","Post_MPR"))

alpha_plot<- ggplot(dat1, aes(x = cell_type, y = value, fill = cell_type)) + 
  stat_boxplot(geom = "errorbar", linewidth=0.5)+
  geom_boxplot(linewidth=0.5)+
  scale_fill_manual(values = cell_type_color) +
  geom_text(data=df1,aes(x=cell_type,y=mean+1.5*std,label=Letters))+
  labs(title = "",
       x = "",
       y = "Simpson") +
  facet_wrap(.~MPR_Response2, scales = "free_y",nrow =3)+
  theme_minimal()+
  theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),units=,"cm"),
        plot.title = element_text(hjust = 0.5, size = 13),
        panel.background = element_blank(),
        plot.background = element_blank(),
        strip.text = element_text(size=12),
        axis.title.y = element_text(size = 14),
        axis.line = element_line(color = "black",size = 0.4),
        axis.text.y = element_text(color="black",size=12),
        axis.text.x = element_text(size=12,colour ="black",angle =45,vjust=0.7,hjust=0.5),
        legend.position = "none",
        panel.spacing = unit(0,"lines"),
        axis.ticks.x=element_line(color="black",size=0.5,lineend = 1),
        axis.line.y=element_line(linetype=1,color="black",size=0.7),
        axis.ticks.y=element_line(color="black",size=0.5,lineend = 1),
        axis.line.x=element_line(linetype=1,color="black",size=0.7),
        panel.grid.major = element_blank(),  # 去除主网格线
        panel.grid.minor = element_blank())+
  coord_cartesian()
alpha_plot
ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/INVADE-seq/aMiAD_score/每个细胞_alpha_diversity_Simpson.pdf",alpha_plot,width=6,height=8)


####单细胞代谢通路分析----

library(stringr)
library(reshape2)
library(scales)
library(scater)
library(pheatmap)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(RColorBrewer)
library(fgsea)

MPRNMPR_object_miMPRobe_remove<- readRDS("../MPRNMPR_object_miMPRobe_remove.rds")
MPRNMPR_object_miMPRobe_remove = UpdateSeuratObject(MPRNMPR_object_miMPRobe_remove)
MPRNMPR_object_miMPRobe_remove

pathway_file <- "D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/INVADE-seq/aMiAD_score/代谢通路分析/KEGG_metabolism_nc.gmt"
pathway_file <-"./KEGG_metabolism_nc.gmt"
pathways <- gmtPathways(pathway_file)
pathway_names <- names(pathways)


{Pathway_act_Score <- function(data,
                               pathways,
                               assay,
                               filterGene=F,
                               Mean_cut=NULL,
                               percent_cut=NULL,
                               all_cell_types,
                               cell_types
){
  
  
  DefaultAssay(data) <- assay
  norm_tpm <- GetAssayData(data, layer = "data")
  norm_tpm <- as.matrix(norm_tpm)
  norm_tpm <- as.data.frame(norm_tpm)
  
  ###filter genes
  if(filterGene==F){
    
    norm_tpm = norm_tpm
    
  }else{
    
    pctExpr <- function(mat){
      nexpr <- apply(mat, 1, function(x) sum(x > 0))
      nexpr / ncol(mat)
    }
    
    norm_tpm$percent_exp = pctExpr(norm_tpm)
    norm_tpm$mean_exp = rowMeans(norm_tpm[,-ncol(norm_tpm)])
    
    norm_tpm <- norm_tpm[which(norm_tpm$mean_exp >Mean_cut & norm_tpm$percent_exp >percent_cut),]
    norm_tpm <- norm_tpm[,-ncol(norm_tpm)]
    norm_tpm <- norm_tpm[,-ncol(norm_tpm)]
  }
  
  
  ##calculate how many pathways of one gene involved.
  num_of_pathways <- function (pathway_file,overlapgenes){
    pathway_names <- names(pathway_file)
    filter_pathways <- list()
    for (p in pathway_names){
      genes <- pathway_file[[p]]
      common_genes <- intersect(genes,overlapgenes)
      if(length(common_genes>=5)){
        filter_pathways[[p]] <- common_genes
      }
    }
    
    all_genes <- unique(as.vector(unlist(filter_pathways)))
    gene_times <- data.frame(num =rep(0,length(all_genes)),row.names = all_genes)
    for(p in pathway_names){
      for(g in filter_pathways[[p]]){
        gene_times[g,"num"] = gene_times[g,"num"]+1
      }
    }
    gene_times
  } 
  
  
  #some genes occur in multiple pathways.
  gene_pathway_number <- num_of_pathways(pathways,rownames(norm_tpm))
  
  ##Calculate the pathway activities
  #mean ratio of genes in each pathway for each cell type
  pathway_names <- names(pathways)
  mean_expression_shuffle <- matrix(NA,nrow=length(pathway_names),ncol=length(cell_types),dimnames = list(pathway_names,cell_types))
  mean_expression_noshuffle <- matrix(NA,nrow=length(pathway_names),ncol=length(cell_types),dimnames = list(pathway_names,cell_types))
  ###calculate the pvalues using shuffle method
  pvalues_mat <- matrix(NA,nrow=length(pathway_names),ncol=length(cell_types),dimnames = (list(pathway_names, cell_types)))
  
  
  for(p in pathway_names){
    genes <- pathways[[p]]
    genes_comm <- intersect(genes, rownames(norm_tpm))
    if(length(genes_comm) < 5) next
    
    pathway_metabolic_tpm <- norm_tpm[genes_comm, ]
    pathway_metabolic_tpm <- pathway_metabolic_tpm[rowSums(pathway_metabolic_tpm)>0,]
    
    mean_exp_eachCellType <- apply(pathway_metabolic_tpm, 1, function(x)by(x, all_cell_types, mean))
    
    #remove genes which are zeros in any celltype to avoid extreme ratio value
    keep <- colnames(mean_exp_eachCellType)[colAlls(mean_exp_eachCellType>0.001)]
    
    if(length(keep)<3) next
    
    #using the loweset value to replace zeros for avoiding extreme ratio value
    pathway_metabolic_tpm <- pathway_metabolic_tpm[keep,]
    pathway_metabolic_tpm <- t( apply(pathway_metabolic_tpm,1,function(x) {x[x<=0] <- min(x[x>0]);x} ))
    
    
    pathway_number_weight = 1 / gene_pathway_number[keep,]
    #
    mean_exp_eachCellType <- apply(pathway_metabolic_tpm, 1, function(x)by(x, all_cell_types, mean))
    ratio_exp_eachCellType <- t(mean_exp_eachCellType) / colMeans(mean_exp_eachCellType)
    #exclude the extreme ratios
    col_quantile <- apply(ratio_exp_eachCellType,2,function(x) quantile(x,na.rm=T))
    col_q1 <- col_quantile["25%",]
    col_q3 <- col_quantile["75%",]
    col_upper <- col_q3 * 3
    col_lower <- col_q1 / 3
    outliers <- apply(ratio_exp_eachCellType,1,function(x) {any( (x>col_upper)|(x<col_lower) )} )
    
    if(sum(!outliers) < 3) next
    
    keep <- names(outliers)[!outliers]
    pathway_metabolic_tpm <- pathway_metabolic_tpm[keep,]
    pathway_number_weight = 1 / gene_pathway_number[keep,]
    mean_exp_eachCellType <- apply(pathway_metabolic_tpm, 1, function(x)by(x, all_cell_types, mean))
    ratio_exp_eachCellType <- t(mean_exp_eachCellType) / colMeans(mean_exp_eachCellType)
    mean_exp_pathway <- apply(ratio_exp_eachCellType,2, function(x) weighted.mean(x, pathway_number_weight/sum(pathway_number_weight)))
    mean_expression_shuffle[p, ] <-  mean_exp_pathway[cell_types]
    mean_expression_noshuffle[p, ] <-  mean_exp_pathway[cell_types]
    
    ##shuffle 5000 times:  
    ##define the functions 
    group_mean <- function(x){
      sapply(cell_types,function(y) rowMeans(pathway_metabolic_tpm[,shuffle_cell_types_list[[x]]==y,drop=F]))
    }
    column_weigth_mean <- function(x){
      apply(ratio_exp_eachCellType_list[[x]],2, function(y) weighted.mean(y, weight_values))
    }
    #####  
    times <- 1:5000
    weight_values <- pathway_number_weight/sum(pathway_number_weight)
    shuffle_cell_types_list <- lapply(times,function(x) sample(all_cell_types)) 
    names(shuffle_cell_types_list) <- times
    mean_exp_eachCellType_list <- lapply(times,function(x) group_mean(x))
    ratio_exp_eachCellType_list <- lapply(times,function(x) mean_exp_eachCellType_list[[x]] / rowMeans(mean_exp_eachCellType_list[[x]]))
    mean_exp_pathway_list <- lapply(times,function(x) column_weigth_mean(x))
    
    shuffle_results <- matrix(unlist(mean_exp_pathway_list),ncol=length(cell_types),byrow = T) 
    rownames(shuffle_results) <- times
    colnames(shuffle_results) <- cell_types
    for(c in cell_types){
      if(is.na(mean_expression_shuffle[p,c])) next
      if(mean_expression_shuffle[p,c]>1){
        pval <- sum(shuffle_results[,c] > mean_expression_shuffle[p,c]) / 5000 
      }else if(mean_expression_shuffle[p,c]<1){
        pval <- sum(shuffle_results[,c] < mean_expression_shuffle[p,c]) / 5000
      }
      if(pval>0.01) mean_expression_shuffle[p, c] <- NA  ### NA is  blank in heatmap
      pvalues_mat[p,c] <- pval
    }
  }
  
  
  analysis_res <- list(mean_expression_shuffle,
                       mean_expression_noshuffle,
                       pvalues_mat)
  
  names(analysis_res) <- c("mean_expression_shuffle",
                           "mean_expression_noshuffle",
                           "pvalues_mat")
  return(analysis_res)
  
}
}

MPRNMPR_object_miMPRobe_remove$newcelltype <- paste0(MPRNMPR_object_miMPRobe_remove$MPR_Response2,"_", MPRNMPR_object_miMPRobe_remove$cell_type_new)

all_cell_types <- as.vector(MPRNMPR_object_miMPRobe_remove$newcelltype)
cell_types <- unique(all_cell_types)

metabolism_activaty <- Pathway_act_Score(MPRNMPR_object_miMPRobe_remove,
                                         pathways=pathways,
                                         assay = "RNA",
                                         filterGene=T,
                                         Mean_cut = 0.001,
                                         percent_cut =0.1,
                                         all_cell_types = all_cell_types,
                                         cell_types = cell_types)


metabolism_activaty <- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/INVADE-seq/aMiAD_score/代谢通路分析/metabolism_activaty.csv",row.names=1)
all_NA <- rowAlls(is.na(metabolism_activaty ))
metabolism_activaty  <- metabolism_activaty [!all_NA,]
View(metabolism_activaty)

dat <- metabolism_activaty

KEGG_anno <- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/INVADE-seq/aMiAD_score/代谢通路分析/kegg_annotation.csv")
rownames(KEGG_anno) <- KEGG_anno$directory3
Anno_pathway <- KEGG_anno[rownames(dat),]
dat$pathway_anno <- Anno_pathway$directory2


#宽数据转化为长数据
dat$pathway <- rownames(dat)
melt_dat = melt(dat,
                id.vars = c("pathway","pathway_anno"),
                measure.vars = 1:(length(dat)-2),
                variable.name = c('celltype'),
                value.name = 'score')

head(melt_dat)


library(dplyr)
library(tidyr)

# melt_dat$MPR_Response2<-sub("_[^_]+$", "", melt_dat$celltype) 
# melt_dat$cell_type_new<-sub("^[^_]*_[^_]*_", "", melt_dat$celltype)


#按照分组排个序
melt_dat <- melt_dat %>%arrange(pathway_anno, desc(score))  

unique(melt_dat$celltype)
melt_dat$celltype <- factor(melt_dat$celltype, levels = c("Pre_NMPR_T.cells","Pre_NMPR_B.cells","Pre_NMPR_Plasma.cells","Pre_NMPR_Myeloid.cells","Pre_NMPR_Mast.cells","Pre_NMPR_Endothelial.cells","Pre_NMPR_Fibroblasts","Pre_NMPR_Smooth.muscle.cells","Pre_NMPR_Pericytes","Pre_NMPR_Epithelial.cells",
                                                          
                                                          "Pre_MPR_T.cells","Pre_MPR_B.cells","Pre_MPR_Plasma.cells","Pre_MPR_Myeloid.cells","Pre_MPR_Mast.cells","Pre_MPR_Endothelial.cells","Pre_MPR_Fibroblasts","Pre_MPR_Smooth.muscle.cells", "Pre_MPR_Pericytes","Pre_MPR_Epithelial.cells",
                                                            
                                                          "Post_MPR_T.cells","Post_MPR_B.cells","Post_MPR_Plasma.cells","Post_MPR_Myeloid.cells","Post_MPR_Mast.cells","Post_MPR_Endothelial.cells","Post_MPR_Fibroblasts", "Post_MPR_Smooth.muscle.cells","Post_MPR_Pericytes","Post_MPR_Epithelial.cells"
                                                          ))


melt_dat$pathway <- factor(melt_dat$pathway, levels = unique(melt_dat$pathway))


p = ggplot(melt_dat,aes(x=celltype,y=pathway, fill=score))+
  geom_tile(color='black')+
  theme(panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 90,hjust=1,vjust=0.5,color = "black", size=8),
        axis.text.y = element_text(size=8,vjust=1,color = "black"),
        legend.margin = margin(-0.2,-0.2,0,0,'cm'))+
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.5, linetype="solid"))+
  scale_y_discrete(position = 'right',expand = c(0,0)) + 
  scale_x_discrete(expand = c(0,0))+
  scale_fill_gradientn(colors=colorRampPalette(c("#1AA3FF","white","#FF6B67"))(100),
                       guide = guide_colorbar(ticks.colour = "black",
                                              frame.colour = "black"),
                       name = "Pathway activity score",
                       values = c(0,0.15,0.3,0.5,1))

p

Anno_row <- rownames(dat) %>% as.data.frame() %>%
  mutate(group=dat$pathway_anno) %>%
  mutate(p="")
colnames(Anno_row)[1] <- 'pathway'
Anno_row$pathway <- factor(Anno_row$pathway, levels = Anno_row$pathway)
Anno_rows <- ggplot(Anno_row, aes(p,pathway,fill=group))+
  geom_tile() + 
  scale_y_discrete(position="right") +
  theme_minimal()+xlab(NULL) + ylab(NULL) +
  theme(axis.text = element_blank())+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))+
  scale_fill_manual(values = c('#855C75', '#D9AF6B', '#AF6458', '#736F4C', '#526A83', '#625377', '#68855C','#9C9C5E', '#8C785D', '#467378'))+
  scale_y_discrete(expand = c(0,0)) + 
  scale_x_discrete(expand = c(0,0))+
  guides(fill=guide_legend(title="Metabolic pathway"))



library(patchwork)
plot <- Anno_rows+p+plot_layout(ncol  = 2, widths  = c(0.2,5),guides = 'collect')
plot & theme(plot.margin = margin(0.5,0.5,0.5,0.5))

##-------------------------------------------------------------------------------------------------
##--------------------------------------------------------------------------------------------------
library(reshape2)
library(pheatmap)
library(RColorBrewer)
library(ggtree)
library(aplot)
library(dplyr)
library(ggplot2)
library(ggheatmap)
BiocManager::install("dittoSeq")
library(dittoSeq)
#####绘图函数---
{ks_gg_heatmap <- function(mat,#矩阵，类似于行是基因，列是分组样本
                           scaleData,#是否需要scale data，不需要“none", 按照行"row", 按照列“column"
                           cluster_method,#聚类方法，参考hclust
                           cluster_row=T,#是否聚类，默认T
                           cluster_col=T,#是否聚类，默认T
                           cols,#热图颜色
                           show_cluster_rows=T,#是否展示聚类树，默认T
                           show_cluster_cols=T,#是否展示聚类树，默认T
                           text_size,#x,y文字大小
                           legend_title,#热图legend标题
                           Discrete_data=F,#输入数据类型，是连续性还是离散型
                           rowname_pos="right",#rowname的位置，默认右侧
                           shape="S",#热图形状，默认为S，S表示矩形，C表示⚪
                           myorder=F,#是否自行设定行列顺序，默认F
                           level_row=NULL,#myorder=T时，行的排序
                           level_col=NULL,#myorder=T时，列的排序
                           anno_col=NULL,#列注释
                           anno_row=NULL,#行注释
                           anno_col_color=NULL,#列注释颜色
                           anno_row_color=NULL,#行注释颜色
                           anno_col_title=NULL,#列注释标题
                           anno_row_title=NULL,#行注释标题
                           label_color=F,#是否标记行名特殊颜色，默认F
                           label1=NULL,#label_color=T时，标记的第一组行名
                           label2=NULL,#label_color=T时，标记的第2组行名
                           label3=NULL,#label_color=T时，标记的第3组行名
                           values
                           
){
  
  scale_rows = function(x){
    m = apply(x, 1, mean, na.rm = T)
    s = apply(x, 1, sd, na.rm = T)
    return((x - m) / s)}
  
  df = switch(scaleData, 
              none = mat, 
              row = scale_rows(mat), 
              column = t(scale_rows(t(mat))))  
  
  
  if(cluster_row==T){
    
    row_clust = hclust(dist(df), method = cluster_method)
    roword = row_clust$order
    
    if(show_cluster_rows==T){
      tree <- ape::as.phylo(row_clust)
      row_tree <- ggtree(tree)+
        theme(legend.position = "none")
    }else{
      
      row_tree=NULL
      
    }
    
    
  }else{
    
    roword=seq(1:nrow(df))
    
  }
  
  if(cluster_col==T){
    
    row_clust = hclust(dist(t(df)), method = cluster_method)
    colord = row_clust$order
    
    if(show_cluster_cols==T){
      tree <- ape::as.phylo(row_clust)
      col_tree <- ggtree(tree,layout = "dendrogram")+
        theme(legend.position = "none")
    }else{
      
      col_tree=NULL
      
    }
    
    
  }else{
    
    colord = seq(1:length(df))
    
  }
  
  
  
  dat=df[roword,colord]
  
  dat <- cbind(rownames(dat),dat)
  colnames(dat)[1] <- 'gene'
  
  melt_dat = melt(dat,
                  id.vars = c("gene"),
                  measure.vars = 2:length(dat),
                  variable.name = c('group'),
                  value.name = 'Exp')
  
  
  if(myorder==F){
    
    melt_dat$gene <- factor(melt_dat$gene, levels = rownames(dat))
    melt_dat$group <- factor(melt_dat$group, levels = colnames(dat))
    
  }else{
    
    melt_dat$gene <- factor(melt_dat$gene, levels = level_row)
    melt_dat$group <- factor(melt_dat$group, levels = level_col)
    
  }
  
  
  if(label_color == T){
    
    melt_dat$color <-   ifelse(melt_dat$gene %in% label1, "#D52126",
                               ifelse(melt_dat$gene %in% label2, "#117733",
                                      ifelse(melt_dat$gene %in% label3,"#E68316","black")))
    
    
  }else{
    
    melt_dat$color <- 'black'
  }
  
  
  
  
  if(Discrete_data==F){
    if(shape=="S"){
      p = ggplot(melt_dat,aes(x=group,y=gene, fill=Exp))+
        geom_tile(color='black')+
        theme(panel.background = element_blank(),
              panel.border = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.title = element_blank(),
              axis.line = element_blank(),
              axis.ticks = element_blank(),
              axis.text.x = element_text(angle = 90,hjust=1,vjust=0.5,color = "black", size=text_size),
              axis.text.y = element_text(size=text_size,vjust=1,
                                         colour =melt_dat$color),
              legend.margin = margin(-0.2,-0.2,0,0,'cm'),
              legend.key.height = unit(0.5,'cm'),
              legend.key.width = unit(0.3,'cm'))+
        scale_y_discrete(position = rowname_pos,expand = c(0,0)) + 
        scale_x_discrete(expand = c(0,0))+
        scale_fill_gradientn(colors=cols,
                             guide = guide_colorbar(ticks.colour = "black",
                                                    frame.colour = "black"),
                             name = legend_title,
                             values = values)
      
    }
    
    if(shape=="C"){
      p = ggplot(melt_dat,aes(x=group,y=gene))+
        geom_point(aes(size=abs(Exp), color=Exp))+
        theme(panel.background = element_blank(),
              panel.border = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.title = element_blank(),
              axis.line = element_blank(),
              axis.ticks = element_blank(),
              axis.text.x = element_text(angle = 90,hjust=1,vjust=0.5,color = "black", size=text_size),
              axis.text.y = element_text(size=text_size,hjust=1,
                                         colour =melt_dat$color),
              legend.margin = margin(-0.2,-0.2,0,0,'cm'),
              legend.key.height = unit(0.5,'cm'),
              legend.key.width = unit(0.3,'cm'))+
        scale_y_discrete(position = rowname_pos, expand = c(0,0)) + 
        scale_x_discrete(expand = c(0,0))+
        scale_color_gradientn(colors=cols,
                              guide = guide_colorbar(ticks.colour = "black",
                                                     frame.colour = "black"),
                              name = legend_title,
                              values = values)
    }
    
    
    
  }else{
    
    melt_dat$Exp <- as.factor(melt_dat$Exp)
    if(shape=="S"){
      p = ggplot(melt_dat,aes(x=group,y=gene,fill=Exp))+
        geom_tile(color='black')+
        theme(panel.background = element_blank(),
              panel.border = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.title = element_blank(),
              axis.line = element_blank(),
              axis.ticks = element_blank(),
              axis.text.x = element_text(angle = 90,hjust=1,vjust=0.5,color = "black", size=text_size),
              axis.text.y = element_text(size=text_size,hjust=1,
                                         colour =melt_dat$color),
              legend.margin = margin(-0.2,-0.2,0,0,'cm'),
              legend.key.height = unit(0.5,'cm'),
              legend.key.width = unit(0.3,'cm'))+
        scale_y_discrete(position = rowname_pos,expand = c(0,0)) + 
        scale_x_discrete(expand = c(0,0))+
        scale_fill_manual(values =cols,
                          name = legend_title)
      
    }
    
    if(shape=="C"){
      p = ggplot(melt_dat,aes(x=group,y=gene,color=Exp))+
        geom_point(aes(size=5))+
        theme(panel.background = element_blank(),
              panel.border = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.title = element_blank(),
              axis.line = element_blank(),
              axis.ticks = element_blank(),
              axis.text.x = element_text(angle = 90,hjust=1,vjust=0.5,color = "black", size=text_size),
              axis.text.y = element_text(size=text_size,vjust=1,
                                         colour =melt_dat$color),
              legend.margin = margin(-0.2,-0.2,0,0,'cm'),
              legend.key.height = unit(0.5,'cm'),
              legend.key.width = unit(0.3,'cm'))+
        scale_y_discrete(position = rowname_pos,expand = c(0,0)) + 
        scale_x_discrete(expand = c(0,0))+
        scale_color_manual(values =cols,
                           name = legend_title)
      
    }
    
    
  }
  # 
  # # plot =p%>%insert_top(col_tree,height = 0.1)
  # # plot = ggplotify::as.ggplot(plot)
  # plot =row_tree+p+plot_layout(ncol = 2, widths  = c(0.5,4),guides = 'collect')
  # plot = plot & theme(plot.margin = margin(0.8,0.8,0.8,0.8))
  
  if(!is.null(anno_col)){
    
    Anno_col <- colnames(mat) %>% as.data.frame() %>% 
      mutate(group=anno_col) %>%
      mutate(p="") 
    Anno_cols <- ggplot(Anno_col, aes(.,y=p,fill=group))+
      geom_tile() + 
      scale_y_discrete(position="right") +
      theme_minimal()+xlab(NULL) + ylab(NULL) +
      theme(axis.text = element_blank())+
      theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))+
      scale_fill_manual(values = anno_col_color)+
      scale_y_discrete(expand = c(0,0)) + 
      scale_x_discrete(expand = c(0,0))+
      guides(fill=guide_legend(title=anno_col_title))
    
  }
  
  
  
  
  if(!is.null(anno_row)){
    Anno_row <- rownames(mat) %>% as.data.frame() %>%
      mutate(group=anno_row) %>%
      mutate(p="")
    Anno_rows <- ggplot(Anno_row, aes(p,.,fill=group))+
      geom_tile() + 
      scale_y_discrete(position="right") +
      theme_minimal()+xlab(NULL) + ylab(NULL) +
      theme(axis.text = element_blank())+
      theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))+
      scale_fill_manual(values = anno_row_color)+
      scale_y_discrete(expand = c(0,0)) + 
      scale_x_discrete(expand = c(0,0))+
      guides(fill=guide_legend(title=anno_row_title))
    
  }
  
  
  # 
  # 
  # if(show_cluster_cols==F & show_cluster_cols==F){
  #   
  #   p <- Anno_cols+p+plot_layout(nrow = 2, heights = c(0.15,5),guides = 'collect')
  #   p <- Anno_rows+p+plot_layout(ncol  = 2, widths  = c(0.15,5),guides = 'collect')
  #   p = p & theme(plot.margin = margin(0.5,0.5,0.5,0.5))
  #   
  # }else{
  
  #layout
  if(show_cluster_rows==F&show_cluster_cols==T){
    if(!is.null(anno_col) &!is.null(anno_row)){
      
      p <-  p%>%insert_left(Anno_rows,width = 0.05)%>%
        insert_top(Anno_cols,height = 0.02)%>%
        insert_top(col_tree,height = 0.1)
    }
    
    if(!is.null(anno_col) & is.null(anno_row)){
      
      p <-  p%>%insert_top(Anno_cols,height = 0.02)%>%
        insert_top(col_tree,height = 0.1)
    }
    
    
    if(is.null(anno_col) &!is.null(anno_row)){
      
      p <-  p%>%insert_left(Anno_rows,width = 0.05)%>%
        insert_top(col_tree,height = 0.1)
    }
    
    
    if(is.null(anno_col) &is.null(anno_row)){
      
      p <-  p%>%insert_top(col_tree,height = 0.1)
    }
    
    
    
  }
  
  
  
  if(show_cluster_rows==T&show_cluster_cols==F){
    if(!is.null(anno_col) &!is.null(anno_row)){
      p <-  p%>%insert_left(Anno_rows,width = 0.05)%>%
        insert_left(row_tree,width = 0.1)%>%
        insert_top(Anno_cols,height = 0.02)
    }
    
    if(!is.null(anno_col) & is.null(anno_row)){
      p <-  p%>%insert_left(row_tree,width = 0.1)%>%
        insert_top(Anno_cols,height = 0.02)
    }
    
    
    if(is.null(anno_col) &!is.null(anno_row)){
      
      p <-  p%>%insert_left(Anno_rows,width = 0.05)%>%
        insert_left(row_tree,width = 0.1)
    }
    
    
    if(is.null(anno_col) &is.null(anno_row)){
      
      p <-  p%>%insert_left(row_tree,width = 0.1)
    }
    
    
    
  }
  
  
  if(show_cluster_rows==F&show_cluster_cols==F){
    
    if(!is.null(anno_col) &!is.null(anno_row)){
      p <-  p%>%insert_left(Anno_rows,width = 0.05)%>%
        insert_top(Anno_cols,height = 0.02)
    }
    
    if(!is.null(anno_col) & is.null(anno_row)){
      p <-  p%>%insert_top(Anno_cols,height = 0.02)
    }
    
    
    if(is.null(anno_col) &!is.null(anno_row)){
      
      p <-  p%>%insert_left(Anno_rows,width = 0.05)
    }
    
    
    if(is.null(anno_col) &is.null(anno_row)){
      
      p <-  p
    }
    
    
  }
  
  
  
  if(show_cluster_rows==T&show_cluster_cols==T){
    if(!is.null(anno_col) &!is.null(anno_row)){
      
      p <-  p%>%insert_left(Anno_rows,width = 0.05)%>%
        insert_left(row_tree,width = 0.1)%>%
        insert_top(Anno_cols,height = 0.02)%>%
        insert_top(col_tree,height = 0.1)
    }
    
    if(!is.null(anno_col) & is.null(anno_row)){
      
      p <-  p%>%insert_left(row_tree,width = 0.1)%>%
        insert_top(Anno_cols,height = 0.02)%>%
        insert_top(col_tree,height = 0.1)
    }
    
    
    if(is.null(anno_col) &!is.null(anno_row)){
      p <-  p%>%insert_left(Anno_rows,width = 0.05)%>%
        insert_left(row_tree,width = 0.1)%>%
        insert_top(col_tree,height = 0.1)
    }
    
    
    if(is.null(anno_col) &is.null(anno_row)){
      
      p <-  p%>%insert_left(row_tree,width = 0.1)%>%
        insert_top(col_tree,height = 0.1)
    }
    
    
  }
  
  
  
  return(p)
  
}
}
#####

anno_col1=c("Pre_NMPR","Pre_NMPR","Pre_NMPR","Pre_NMPR","Pre_NMPR","Pre_NMPR","Pre_NMPR","Pre_NMPR","Pre_NMPR","Pre_NMPR",
            
            "Pre_MPR","Pre_MPR","Pre_MPR","Pre_MPR","Pre_MPR","Pre_MPR","Pre_MPR","Pre_MPR", "Pre_MPR","Pre_MPR",
            
            "Post_MPR","Post_MPR","Post_MPR","Post_MPR","Post_MPR","Post_MPR","Post_MPR", "Post_MPR","Post_MPR","Post_MPR"
)
level_col1= c("Pre_NMPR_T.cells","Pre_NMPR_B.cells","Pre_NMPR_Plasma.cells","Pre_NMPR_Myeloid.cells","Pre_NMPR_Mast.cells","Pre_NMPR_Endothelial.cells","Pre_NMPR_Fibroblasts","Pre_NMPR_Smooth.muscle.cells","Pre_NMPR_Pericytes","Pre_NMPR_Epithelial.cells",
              
              "Pre_MPR_T.cells","Pre_MPR_B.cells","Pre_MPR_Plasma.cells","Pre_MPR_Myeloid.cells","Pre_MPR_Mast.cells","Pre_MPR_Endothelial.cells","Pre_MPR_Fibroblasts","Pre_MPR_Smooth.muscle.cells", "Pre_MPR_Pericytes","Pre_MPR_Epithelial.cells",
              
              "Post_MPR_T.cells","Post_MPR_B.cells","Post_MPR_Plasma.cells","Post_MPR_Myeloid.cells","Post_MPR_Mast.cells","Post_MPR_Endothelial.cells","Post_MPR_Fibroblasts", "Post_MPR_Smooth.muscle.cells","Post_MPR_Pericytes","Post_MPR_Epithelial.cells"
)

library(ggplot2)
library(dplyr)
library(tidyr)

dat <- dat[,c(1:30)]  %>% mutate(across(everything(), ~ replace_na(., 0))) 
names(dat)
library(ComplexHeatmap)
library(circlize)

anno_col1=c("Pre_NMPR","Pre_NMPR","Pre_NMPR","Pre_NMPR","Pre_NMPR","Pre_NMPR","Pre_NMPR","Pre_NMPR","Pre_NMPR","Pre_NMPR",
            
            "Pre_MPR","Pre_MPR","Pre_MPR","Pre_MPR","Pre_MPR","Pre_MPR","Pre_MPR","Pre_MPR", "Pre_MPR","Pre_MPR",
            
            "Post_MPR","Post_MPR","Post_MPR","Post_MPR","Post_MPR","Post_MPR","Post_MPR", "Post_MPR","Post_MPR","Post_MPR"
)
level_col1= c("Pre_NMPR_T.cells","Pre_NMPR_B.cells","Pre_NMPR_Plasma.cells","Pre_NMPR_Myeloid.cells","Pre_NMPR_Mast.cells","Pre_NMPR_Endothelial.cells","Pre_NMPR_Fibroblasts","Pre_NMPR_Smooth.muscle.cells","Pre_NMPR_Pericytes","Pre_NMPR_Epithelial.cells",
              
              "Pre_MPR_T.cells","Pre_MPR_B.cells","Pre_MPR_Plasma.cells","Pre_MPR_Myeloid.cells","Pre_MPR_Mast.cells","Pre_MPR_Endothelial.cells","Pre_MPR_Fibroblasts","Pre_MPR_Smooth.muscle.cells", "Pre_MPR_Pericytes","Pre_MPR_Epithelial.cells",
              
              "Post_MPR_T.cells","Post_MPR_B.cells","Post_MPR_Plasma.cells","Post_MPR_Myeloid.cells","Post_MPR_Mast.cells","Post_MPR_Endothelial.cells","Post_MPR_Fibroblasts", "Post_MPR_Smooth.muscle.cells","Post_MPR_Pericytes","Post_MPR_Epithelial.cells"
)

ks_gg_heatmap_p<-ks_gg_heatmap(dat, scaleData="none",
              # cluster_method='ward.D2',
              cluster_row = F, 
              cluster_col = F,
              show_cluster_rows = F,
              show_cluster_cols =F,
              anno_col = anno_col1,
              text_size = 8,
              legend_title = "Pathway activity score",
              cols = colorRampPalette(c("#1AA3FF","#C1FFC1","white","#EED5D2","brown"))(100),
              #cols = rev(colorRampPalette(brewer.pal(11, "PiYG"))(25)),
              #cols =colorRamp2(c(0, 1, 2), c("white", "yellow", "red")),
              rowname_pos = "right",
              shape="S",
              anno_col_color  = c("#4974a4","#4dae47","#f29600"),
              anno_row = Anno_pathway$directory2,
              anno_row_color = dittoColors(),
              anno_col_title = "group",
              anno_row_title = "Metabolic pathway",
              label_color = F,
              myorder = T,
              level_col = level_col1,
              level_row = rownames(dat),
              values = c(0,0.15,0.3,0.5,1))

ks_gg_heatmap_p

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/INVADE-seq/aMiAD_score/代谢通路分析/ks_gg_heatmap.pdf",ks_gg_heatmap_p,width=12,height=9)

metabolism_activaty_2<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/INVADE-seq/aMiAD_score/代谢通路分析/metabolism_activaty_2.csv",row.names = 1)
metabolism_activaty<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/INVADE-seq/aMiAD_score/代谢通路分析/metabolism_activaty.csv",row.names = 1)
scRNA_dat <- as.data.frame(metabolism_activaty)
scRNA_dat$X <- NULL


scRNA_df <- melt(scRNA_dat)
scRNA_df <- scRNA_df[!is.na(scRNA_df$value),]

scRNA_df$MPR_Response2<-sub("_[^_]+$", "", scRNA_df$variable) 
scRNA_df$cell_type_new<-sub("^[^_]*_[^_]*_", "", scRNA_df$variable)
cell_type_color <- c("#FF6347","#ffc556","#1E90FF","#0000FF","#9370DB","#00FFFF","#32cd32","#006400","#FF00ff","#20B2AA",
                     "#FF6347","#ffc556","#1E90FF","#0000FF","#9370DB","#00FFFF","#32cd32","#006400","#FF00ff","#20B2AA",
                     "#FF6347","#ffc556","#1E90FF","#0000FF","#9370DB","#00FFFF","#32cd32","#006400","#FF00ff","#20B2AA")

scRNA_df$variable<- factor(scRNA_df$variable,levels=c("Pre_NMPR_T.cells","Pre_NMPR_B.cells","Pre_NMPR_Plasma.cells","Pre_NMPR_Myeloid.cells","Pre_NMPR_Mast.cells","Pre_NMPR_Endothelial.cells","Pre_NMPR_Fibroblasts","Pre_NMPR_Smooth.muscle.cells","Pre_NMPR_Pericytes","Pre_NMPR_Epithelial.cells",
                                                      
                                                      "Pre_MPR_T.cells","Pre_MPR_B.cells","Pre_MPR_Plasma.cells","Pre_MPR_Myeloid.cells","Pre_MPR_Mast.cells","Pre_MPR_Endothelial.cells","Pre_MPR_Fibroblasts","Pre_MPR_Smooth.muscle.cells", "Pre_MPR_Pericytes","Pre_MPR_Epithelial.cells",
                                                      
                                                      "Post_MPR_T.cells","Post_MPR_B.cells","Post_MPR_Plasma.cells","Post_MPR_Myeloid.cells","Post_MPR_Mast.cells","Post_MPR_Endothelial.cells","Post_MPR_Fibroblasts", "Post_MPR_Smooth.muscle.cells","Post_MPR_Pericytes","Post_MPR_Epithelial.cells"))

scRNA_df$cell_type_new<- gsub("\\.", " ", scRNA_df$cell_type_new)

scRNA_df$cell_type_new <- factor(scRNA_df$cell_type_new ,levels=c('T cells', 'B cells','Plasma cells', 'Myeloid cells',  
  'Mast cells',"Endothelial cells","Fibroblasts","Smooth muscle cells",
  "Pericytes","Epithelial cells"))


scRNA_df$MPR_Response2<- factor(scRNA_df$MPR_Response2,levels=c("Pre_NMPR","Pre_MPR","Post_MPR"))
g <- ggplot(scRNA_df,aes(x=cell_type_new,y=value,fill=cell_type_new)) +
  # scale_y_continuous(limits=c(0,2),breaks=0:2,labels=0:2)+
  geom_violin(trim=F,size=0.2,show.legend = F,width=1.0) + labs(y=NULL,x=NULL) + 
  stat_summary(fun.y = median,geom="point",size=1,color="blue")+
  scale_fill_manual(values=cell_type_color) +
  facet_wrap(MPR_Response2~.,nrow=3)+
  theme_classic() + 
  theme(legend.position="none",
        axis.text.x=element_text(colour="black", size = 10,angle=45,hjust=1,vjust=1),
        axis.text.y=element_text(colour="black", size = 10),
        axis.line=element_line(size=0.2,color="black"),
        axis.ticks = element_line(colour = "black",size=0.2),
        panel.border = element_blank(), panel.background = element_blank(),
        axis.ticks.length= unit(.5, "mm"))+
  geom_hline(yintercept = 1)+
  labs(y="Pathway activity score")

g
ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/INVADE-seq/aMiAD_score/代谢通路分析/Pathway activity score_小提琴.pdf",g,width=6,height=7)
head(scRNA_df)

###boxplot----------------------------

library(ggpubr)

scRNA_df_T_cells<-scRNA_df %>% filter(cell_type_new %in% c("T cells")) 
names(scRNA_df_T_cells)

gg1 <- ggplot(data=scRNA_df_T_cells, aes(x = variable, y = value,
                            fill=variable)) +
  geom_boxplot(alpha =0.5,size=0.5,outlier.shape = NA)+
  scale_fill_manual(values=c("#4974a4","#4dae47","#f29600"))+
  stat_compare_means(method = "t.test",paired = F, 
                     comparisons=list(c("Pre_NMPR_T.cells", "Pre_MPR_T.cells"),
                                      c("Pre_NMPR_T.cells", "Post_MPR_T.cells"),
                                      c("Pre_MPR_T.cells", "Post_MPR_T.cells")))+
  geom_jitter(alpha = 0.3,size=3, shape=21)+
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))+
  theme_bw() + 
  theme(panel.grid =element_blank(),
        axis.text = element_text(size = 10,colour = "black"),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = 'none')+
  labs(y='Pathway activity score')


gg2 <- ggplot(data=scRNA_df_T_cells, aes(x = variable, y = value,
                                         fill=variable)) +
  geom_violin(alpha =0.5,size=0.5,outlier.shape = NA)+
  geom_boxplot(alpha =0.5,size=0.5,outlier.shape = NA, width=0.2, fill="white",color="black")+
  scale_fill_manual(values=c("#4974a4","#4dae47","#f29600"))+
  stat_compare_means(method = "t.test",paired = F, 
                     comparisons=list(c("Pre_NMPR_T.cells", "Pre_MPR_T.cells"),
                                      c("Pre_NMPR_T.cells", "Post_MPR_T.cells"),
                                      c("Pre_MPR_T.cells", "Post_MPR_T.cells")))+
  geom_jitter(alpha = 0.3,size=3, shape=21)+
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))+
  theme_bw() + 
  labs(y='Pathway activity score',title = "T cells")+
  theme(panel.grid =element_blank(),
        axis.text = element_text(size = 12,colour = "black"),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(colour = "black", size = 12, angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(colour = "black", size = 12),
        axis.line = element_line(size = 0.2, color = "black"),
        legend.position = 'none',
        plot.title = element_text(
          hjust = 0.5,
          size = 14,  
          face = "bold"
        ))
  
gg2


#####循环画图
library(ggplot2)
library(cowplot)
library(dplyr)

# 细胞类型组和颜色
cell_type_groups <- c('T cells', 'B cells', 'Plasma cells', 'Myeloid cells',  
                      'Mast cells', 'Endothelial cells', 'Fibroblasts', 'Smooth muscle cells',
                      'Pericytes', 'Epithelial cells')

treatment_color <- c("#4974a4", "#4dae47", "#f29600")

comparison_list <- list(
  'T cells' = list(c("Pre_NMPR_T.cells", "Pre_MPR_T.cells"), 
                   c("Pre_NMPR_T.cells", "Post_MPR_T.cells"),
                   c("Pre_MPR_T.cells", "Post_MPR_T.cells")),
  'B cells' = list(c("Pre_NMPR_B.cells", "Pre_MPR_B.cells"), 
                   c("Pre_NMPR_B.cells", "Post_MPR_B.cells"),
                   c("Pre_MPR_B.cells", "Post_MPR_B.cells")),
  'Plasma cells' = list(c("Pre_NMPR_Plasma.cells", "Pre_MPR_Plasma.cells"), 
                        c("Pre_NMPR_Plasma.cells", "Post_MPR_Plasma.cells"),
                        c("Pre_MPR_Plasma.cells", "Post_MPR_Plasma.cells")),
  'Myeloid cells' =list(c("Pre_NMPR_Myeloid cells", "Pre_MPR_Myeloid cells"), 
                        c("Pre_NMPR_Myeloid cells", "Post_MPR_Myeloid cells"),
                        c("Pre_MPR_Myeloid cells", "Post_MPR_Myeloid cells")))

# 用于存储图形的列表
pplist <- list()

# 循环绘制每种细胞类型的图表
for (cell_type in cell_type_groups) {
  # 过滤数据
  scRNA_df_cell_type <- scRNA_df %>% filter(cell_type_new == cell_type)
  
  # 创建图表
  gg2 <- ggplot(data = scRNA_df_cell_type, aes(x = MPR_Response2, y = value, fill = variable)) +
    geom_violin(alpha = 0.5, size = 0.5, outlier.shape = NA) +
    geom_boxplot(alpha = 0.5, size = 0.5, outlier.shape = NA, width = 0.2, fill = "white", color = "black") +
    scale_fill_manual(values = treatment_color) +
    stat_compare_means(method = "t.test", paired = FALSE, 
                       comparisons =list(c("Pre_NMPR", "Pre_MPR"),
                                                       c("Pre_NMPR", "Post_MPR"),
                                                       c("Pre_MPR", "Post_MPR"))
                       ) +
    geom_jitter(alpha = 0.3, size = 2, shape = 21) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
    theme_bw() + 
    labs(y = 'Pathway activity score', title = cell_type) +
    theme(
      panel.grid = element_blank(),
      axis.text = element_text(size = 12, colour = "black"),
      axis.title.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.x=element_blank(),
      #axis.text.x = element_text(colour = "black", size = 12, angle = 45, hjust = 1, vjust = 1),
      axis.text.y = element_text(colour = "black", size = 12),
      axis.line = element_line(size = 0.2, color = "black"),
      legend.position = 'none',
      plot.title = element_text(hjust = 0.5, size = 12, face = "bold")
    )
  
  # 将图表存储到列表中
  pplist[[cell_type]] <- gg2
}

# 使用 plot_grid 将所有图表组合在一起
cells_plot <- plot_grid(
  plotlist = pplist,
  align = "h",
  nrow = 2
)

# 显示合并后的图表
print(cells_plot)






######scMetabolism----

#修改函数KS.sc.metabolism.score
{KS.sc.metabolism.score <- function(obj, 
                                    method = "VISION", #AUCell,ssGSEA,GSVA
                                    imputation = F, 
                                    ncores = 2, 
                                    metabolism.type = "KEGG") {
  
  require(Seurat)
  require(GSVA)
  require(VISION)
  require(AUCell)
  require(scMetabolism)
  
  data_type <- class(obj)
  
  if(data_type == "Seurat"){
    
    countexp<-GetAssayData(obj, layer='counts')
    countexp<-as.matrix(countexp)
    
  }else{
    
    countexp <- obj
    
  }
  
  
  #imputation
  if (imputation == F) {
    countexp2<-countexp
  }
  if (imputation == T) {
    
    result.completed <- alra(as.matrix(countexp))
    countexp2 <- result.completed[[3]]; row.names(countexp2) <- row.names(countexp)
  }
  
  
  signatures_KEGG_metab <- system.file("data", "KEGG_metabolism_nc.gmt", package = "scMetabolism")
  signatures_REACTOME_metab <- system.file("data", "REACTOME_metabolism.gmt", package = "scMetabolism")
  
  
  if (metabolism.type == "KEGG")  {gmtFile<-signatures_KEGG_metab; cat("Your choice is: KEGG\n")}
  if (metabolism.type == "REACTOME")  {gmtFile<-signatures_REACTOME_metab; cat("Your choice is: REACTOME\n")}
  
  
  
  
  
  #VISION
  if (method == "VISION") {
    library(VISION)
    n.umi <- colSums(countexp2)
    scaled_counts <- t(t(countexp2) / n.umi) * median(n.umi)
    vis <- Vision(scaled_counts, signatures = gmtFile)
    
    options(mc.cores = ncores)
    
    vis <- analyze(vis)
    
    signature_exp<-data.frame(t(vis@SigScores))
  }
  
  
  
  #AUCell
  if (method == "AUCell") {
    library(AUCell)
    library(GSEABase)
    cells_rankings <- AUCell_buildRankings(as.matrix(countexp2), nCores=ncores, plotStats=F) #rank
    geneSets <- getGmt(gmtFile) #signature read
    cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings) #calc
    signature_exp <- data.frame(getAUC(cells_AUC))
  }
  
  #ssGSEA
  if (method == "ssGSEA") {
    library(GSVA)
    library(GSEABase)
    geneSets <- getGmt(gmtFile) #signature read
    gsva_es <- gsva(as.matrix(countexp2), geneSets, method=c("ssgsea"), kcdf=c("Poisson"), parallel.sz=ncores) #
    signature_exp<-data.frame(gsva_es)
  }
  
  #GSVA
  if (method == "ssGSEA") {
    library(GSVA)
    library(GSEABase)
    geneSets <- getGmt(gmtFile) #signature read
    gsva_es <- gsva(as.matrix(countexp2), geneSets, method=c("gsva"), kcdf=c("Poisson"), parallel.sz=ncores) #
    signature_exp<-data.frame(gsva_es)
  }
  
  
  
  
  
  if(data_type == "Seurat"){
    
    obj@assays$METABOLISM$score<-signature_exp
    return(obj)
    
    
  }else{
    
    return(signature_exp)
    
  }
  
  
}
}

library(scMetabolism)
library(Seurat)
library(ggplot2)
library(rsvd)
library(pheatmap)

#V4

MPRNMPR_object_miMPRobe_remove<- readRDS("../MPRNMPR_object_miMPRobe_remove.rds")
MPRNMPR_object_miMPRobe_remove = UpdateSeuratObject(MPRNMPR_object_miMPRobe_remove)

Macrophage_metabolism <- KS.sc.metabolism.score(obj = Macrophage, 
                                                method = 'VISION')

Macrophage_metabolism <- KS.sc.metabolism.score(obj = Macrophage, 
                                                method = 'VISION')

library(devtools)
install_github("YosefLab/VISION")

devtools::install_github("wu-yc/scMetabolism")
#v5
DefaultAssay(MPRNMPR_object_miMPRobe_remove)<- 'RNA'
sceV5_T_cells <- subset(MPRNMPR_object_miMPRobe_remove, cell_type_new=="T cells") 

T_cells_metabolism <- KS.sc.metabolism.score(obj = sceV5_T_cells, method = 'AUCell')
View(T_cells_metabolism)
colnames(T_cells_metabolism @assays[["METABOLISM"]][["score"]])

####每个样本----

DefaultAssay(MPRNMPR_object_miMPRobe_remove)<- 'RNA'

MPRNMPR_object_miMPRobe_remove_metabolism <- KS.sc.metabolism.score(obj = MPRNMPR_object_miMPRobe_remove, method = 'AUCell')

metabolism_score_df <- t(MPRNMPR_object_miMPRobe_remove_metabolism @assays[["METABOLISM"]][["score"]]) %>% 
  cbind("orig.ident"=MPRNMPR_object_miMPRobe_remove$orig.ident,"cell_type_new"=MPRNMPR_object_miMPRobe_remove$cell_type_new,"MPR_Response2"=MPRNMPR_object_miMPRobe_remove$MPR_Response2)%>% as.data.frame()

MPRNMPR_object_miMPRobe_remove$merge_cell_type<-paste0(MPRNMPR_object_miMPRobe_remove$MPR_Response2,"_",MPRNMPR_object_miMPRobe_remove$cell_type_new)

MPRNMPR_object_miMPRobe_remove_metabolism$merge_cell_type<-paste0(MPRNMPR_object_miMPRobe_remove_metabolism$MPR_Response2,"_",MPRNMPR_object_miMPRobe_remove_metabolism$cell_type_new)

MPRNMPR_object_miMPRobe_remove_metabolism$merge_cell_type<- factor(MPRNMPR_object_miMPRobe_remove_metabolism$merge_cell_type,levels=c("Pre_NMPR_T cells","Pre_NMPR_B cells","Pre_NMPR_Plasma cells","Pre_NMPR_Myeloid cells","Pre_NMPR_Mast cells","Pre_NMPR_Endothelial cells","Pre_NMPR_Fibroblasts","Pre_NMPR_Smooth muscle cells","Pre_NMPR_Pericytes","Pre_NMPR_Epithelial cells",
                                                      
                                                      "Pre_MPR_T cells","Pre_MPR_B cells","Pre_MPR_Plasma cells","Pre_MPR_Myeloid cells","Pre_MPR_Mast cells","Pre_MPR_Endothelial cells","Pre_MPR_Fibroblasts","Pre_MPR_Smooth muscle cells", "Pre_MPR_Pericytes","Pre_MPR_Epithelial cells",
                                                      
                                                      "Post_MPR_T cells","Post_MPR_B cells","Post_MPR_Plasma cells","Post_MPR_Myeloid cells","Post_MPR_Mast cells","Post_MPR_Endothelial cells","Post_MPR_Fibroblasts", "Post_MPR_Smooth muscle cells","Post_MPR_Pericytes","Post_MPR_Epithelial cells"))

saveRDS(MPRNMPR_object_miMPRobe_remove_metabolism,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/INVADE-seq/aMiAD_score/代谢通路分析/MPRNMPR_object_miMPRobe_remove_metabolism.rds")
MPRNMPR_object_miMPRobe_remove_metabolism<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/INVADE-seq/aMiAD_score/代谢通路分析/MPRNMPR_object_miMPRobe_remove_metabolism.rds")

MPRNMPR_object_miMPRobe_remove_metabolism$MPR_Response2<- factor(MPRNMPR_object_miMPRobe_remove_metabolism$MPR_Response2,levels=c("Pre_NMPR","Pre_MPR","Post_MPR"))
label<-c("Pre_NMPR","Pre_MPR","Post_MPR")

rownames(MPRNMPR_object_miMPRobe_remove_metabolism @assays[["METABOLISM"]][["score"]])


DotPlot <- DotPlot.metabolism(obj = MPRNMPR_object_miMPRobe_remove_metabolism, 
                   pathway = pathway_name, 
                   phenotype = "MPR_Response2", norm = "y")


ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/INVADE-seq/aMiAD_score/代谢通路分析/DotPlot.metabolism_MPR.pdf",DotPlot,width=10,height=15)



saveRDS(MPRNMPR_object_miMPRobe_remove, file = "D:/wangys_uestc_data2/OSCC_data/oscc_process/rds数据/MPRNMPR_object_miMPRobe_remove.rds")


metabolism_score_df <- t(MPRNMPR_object_miMPRobe_remove_metabolism @assays[["METABOLISM"]][["score"]]) %>% 
  cbind("orig.ident"=MPRNMPR_object_miMPRobe_remove$orig.ident,
        "cell_type_new"=MPRNMPR_object_miMPRobe_remove$cell_type_new,
        "MPR_Response2"=MPRNMPR_object_miMPRobe_remove$MPR_Response2,
        "merge_cell_type"=MPRNMPR_object_miMPRobe_remove$merge_cell_type)%>% as.data.frame()
names(metabolism_score_df)


metabolism_score_df[1:85] <-  metabolism_score_df[1:85] %>% 
  mutate(across(everything(), as.numeric))

metabolism_score_avg  <- metabolism_score_df %>%
  group_by(orig.ident) %>%
  summarise(across(1:85, mean, na.rm = TRUE))
dim(metabolism_score_avg)
rownames(metabolism_score_avg)
metabolism_score_avg$orig.ident

metabolism_score_avg  <- metabolism_score_df %>%
  group_by(MPR_Response2) %>%
  summarise(across(1:85, mean, na.rm = TRUE))
rownames(metabolism_score_avg)<-metabolism_score_avg$MPR_Response2
metabolism_score_avg<-metabolism_score_avg[,-1] %>% as.data.frame()

anno_col1=c("Pre_NMPR","Pre_MPR","Post_MPR")
metabolism_score_avg_t<-t(metabolism_score_avg) %>% as.data.frame()
colnames(metabolism_score_avg_t)
metabolism_score_avg_t<-metabolism_score_avg_t[-1,] %>% as.data.frame()
dat
metabolism_score_avg_t <- metabolism_score_avg_t[1:3]  %>% mutate(across(everything(), ~ replace_na(., 0)))

metabolism_score_avg_t$pathway <- rownames(metabolism_score_avg_t)


library(aMiAD)
library(phyloseq)
library(picante)
library(entropart)
library(vegan)

aMiAD_scores <- list()

columns_to_analyze <- rownames(MPRNMPR_object_miMPRobe_remove_metabolism @assays[["METABOLISM"]][["score"]])[1:85]

dim(metabolism_score_avg)

for (col_name in columns_to_analyze) {
  y <- metabolism_score_avg[[col_name]]  # 提取当前列的数据
  # 运行 aMiAD 函数
  aMiAD_scores[[col_name]] <- aMiAD(alpha1, y, cov = NULL, model = c("gaussian", "binomial"), n.perm = 5000)
  # 提取结果
  print(col_name)
}


all_item_by_item <- data.frame()
all_a_miad <- data.frame()

for (col_name in columns_to_analyze) {
  # 提取当前列的结果
  item_by_item_out <- aMiAD_scores[[col_name]]$ItembyItem.out
  a_miad_out <- aMiAD_scores[[col_name]]$aMiAD.out
  
  # 将结果转换为数据框
  df_item_by_item <- as.data.frame(item_by_item_out)
  df_a_miad <- as.data.frame(a_miad_out)
  
  # 添加列名标识
  df_item_by_item$col_name <- col_name
  df_a_miad$col_name <- col_name
  
  # 将当前迭代的结果添加到总的合并数据框中
  all_item_by_item <- bind_rows(all_item_by_item, df_item_by_item)
  all_a_miad <- bind_rows(all_a_miad, df_a_miad)
}


# 打印结果
print(all_item_by_item)
print(all_a_miad)

write.csv(all_item_by_item,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/INVADE-seq/aMiAD_score/aMiAD_metabolism_score_all_item_by_item2样本.csv")

write.csv(all_a_miad,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/INVADE-seq/aMiAD_score/aMiAD_metabolism_score_all_a_miad2样本.csv")

all_a_miad$aMiDivES <- rownames(all_a_miad)
all_a_miad_filter <- all_a_miad %>% filter(str_detect(aMiDivES, "aMiDivES"))
all_a_miad_filter$a_miad_out <- as.numeric(all_a_miad_filter$a_miad_out)
View(all_a_miad_filter)
all_a_miad_filter$col_name<- factor(all_a_miad_filter$col_name,levels=rev(order))
all_a_miad_filter$col_name
all_a_miad_filter <- na.omit(all_a_miad_filter)

{order<-c("Vitamin B6 metabolism",
"Alanine, aspartate and glutamate metabolism",
"alpha−Linolenic acid metabolism",
"Amino sugar and nucleotide sugar metabolism",
"Arachidonic acid metabolism",
"Arginine and proline metabolism",
"Arginine biosynthesis",
"Ascorbate and aldarate metabolism",
"beta−Alanine metabolism",
"Biosynthesis of unsaturated fatty acids",
"Biotin metabolism",
"Butanoate metabolism",
"Caffeine metabolism",
"Citrate cycle (TCA cycle)",
"Cysteine and methionine metabolism",
"D−Arginine and D−ornithine metabolism",
"D−Glutamine and D−glutamate metabolism",
"Drug metabolism − cytochrome P450",
"Drug metabolism − other enzymes",
"Ether lipid metabolism",
"Fatty acid biosynthesis",
"Fatty acid degradation",
"Fatty acid elongation",
"Folate biosynthesis",
"Fructose and mannose metabolism",
"Galactose metabolism",
"Glutathione metabolism",
"Glycerolipid metabolism",
"Glycerophospholipid metabolism",
"Glycine, serine and threonine metabolism",
"Glycolysis / Gluconeogenesis",
"Glycosaminoglycan biosynthesis − chondroitin sulfate / dermatan sulfate",
"Glycosaminoglycan biosynthesis − heparan sulfate / heparin",
"Glycosaminoglycan biosynthesis − keratan sulfate",
"Glycosaminoglycan degradation",
"Glycosphingolipid biosynthesis − ganglio series",
"Glycosphingolipid biosynthesis − globo and isoglobo series",
"Glycosphingolipid biosynthesis − lacto and neolacto series",
"Glycosylphosphatidylinositol (GPI)−anchor biosynthesis",
"Glyoxylate and dicarboxylate metabolism",
"Histidine metabolism",
"Inositol phosphate metabolism",
"Linoleic acid metabolism",
"Lipoic acid metabolism",
"Lysine degradation",
"Mannose type O−glycan biosynthesis",
"Metabolism of xenobiotics by cytochrome P450",
"Mucin type O−glycan biosynthesis",
"N−Glycan biosynthesis",
"Neomycin, kanamycin and gentamicin biosynthesis",
"Nicotinate and nicotinamide metabolism",
"Nitrogen metabolism",
"One carbon pool by folate",
"Other glycan degradation",
"Other types of O−glycan biosynthesis",
"Oxidative phosphorylation",
"Pantothenate and CoA biosynthesis",
"Pentose and glucuronate interconversions",
"Pentose phosphate pathway",
"Phenylalanine metabolism",
"Phenylalanine, tyrosine and tryptophan biosynthesis",
"Phosphonate and phosphinate metabolism",
"Porphyrin and chlorophyll metabolism",
"Primary bile acid biosynthesis",
"Propanoate metabolism",
"Purine metabolism",
"Pyrimidine metabolism",
"Pyruvate metabolism",
"Retinol metabolism",
"Riboflavin metabolism",
"Selenocompound metabolism",
"Sphingolipid metabolism",
"Starch and sucrose metabolism",
"Steroid biosynthesis",
"Steroid hormone biosynthesis",
"Sulfur metabolism",
"Synthesis and degradation of ketone bodies",
"Taurine and hypotaurine metabolism",
"Terpenoid backbone biosynthesis",
"Thiamine metabolism",
"Tryptophan metabolism",
"Tyrosine metabolism",
"Ubiquinone and other terpenoid−quinone biosynthesis",
"Valine, leucine and isoleucine biosynthesis",
"Valine, leucine and isoleucine degradation")
}

p<- ggplot(all_a_miad_filter, aes(x =a_miad_out, y = col_name, fill = a_miad_out)) +
  geom_bar(stat = "identity") +
  #ylim(-1,1) +
  scale_fill_gradient2(low = "blue", mid = "grey", high = "red", midpoint = 0) +
  labs(title = "Metabolism_pathways",
       x = "aMiDivES",
       y = "",
       fill = "aMiDivES") +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = NA, color = "black"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        strip.text.y.left = element_text(angle = 0),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(angle = 90, size=8, hjust = 1, vjust = 0.5, color = "black"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7))
 # coord_polar(start = 0)

p  
#p+geom_text(data=label_data, aes(x=id, y=value+1, label=individual, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=3, angle= label_data$angle, inherit.aes = FALSE ) 

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/INVADE-seq/aMiAD_score/代谢通路分析/metabolism_scorey_aMiDivES.pdf",p,width =7,height=9)


data <- data.frame(individual=rownames(MPRNMPR_object_miMPRobe_remove_metabolism @assays[["METABOLISM"]][["score"]])[1:85],
                   group=c( rep('A', 10), rep('B', 30), rep('C', 14), rep('D', 6),rep("E",25)),
                   value=sample(-1, 85, replace=T))

data <- data.frame(individual=paste( "Mister ", seq(1,60), sep=""),
                   group=c( rep('A', 10), rep('B', 30), rep('C', 14), rep('D', 6)) ,
                   value=sample( seq(10,100), 60, replace=T))

# 设置空值作为组别间隙
empty_bar <- 4
to_add <- data.frame( matrix(NA, empty_bar*nlevels(data$group), ncol(data)) )
colnames(to_add) <- colnames(data)
to_add$group <- rep(levels(data$group), each=empty_bar)
data <- rbind(data, to_add)
data <- data %>% arrange(group)
data$id <- seq(1, nrow(data))

# 设置标签
label_data <- data
number_of_bar <- nrow(label_data)
angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
label_data$hjust <- ifelse( angle < -90, 1, 0)
label_data$angle <- ifelse(angle < -90, angle+180, angle)

# Fib_Seurat<-scMetabolism::sc.metabolism.Seurat(obj = sceV5_Fibs,
#                                             method = "AUCell", 
#                                             imputation =F, 
#                                             ncores = 2, 
#                                             metabolism.type = "KEGG")


DimPlot.metabolismV5 <- function (obj, pathway, dimention.reduction.type = "umap.harmony", dimention.reduction.run = T, 
                                  size = 1) 
{
  cat("\nPlease Cite: \nYingcheng Wu, Qiang Gao, et al. Cancer Discovery. 2021. \nhttps://pubmed.ncbi.nlm.nih.gov/34417225/   \n\n")
  if (dimention.reduction.type == "umap.harmony") {
    if (dimention.reduction.run == T) 
    umap.loc <- obj@reductions$umap.harmony@cell.embeddings
    row.names(umap.loc) <- colnames(obj)
    signature_exp <- obj@assays$METABOLISM$score
    input.pathway <- pathway
    signature_ggplot <- data.frame(umap.loc,  t(signature_exp[input.pathway, 
    ]),obj$cell_type_new)
    library(wesanderson)
    pal <- wes_palette("Zissou1", 100, type = "continuous")
    library(ggplot2)
    plot1 <- ggplot(data = signature_ggplot, aes(x = umapharmony_1, 
                                                y = umapharmony_2, color = signature_ggplot[, 3])) + geom_point(size = size) + 
      scale_fill_gradientn(colours = pal) + scale_color_gradientn(colours = pal) + 
      labs(color = input.pathway) + xlab("UMAP_1") + ylab("UMAP_2") + 
      theme(aspect.ratio = 1) + theme(panel.grid.major = element_blank(), 
                                      panel.grid.minor = element_blank(), panel.background = element_blank(), 
                                      axis.line = element_line(colour = "black"))
    
    signature_ggplot_med <- signature_ggplot %>%
      group_by(obj.cell_type_new) %>%
      summarise(
        umapharmony_1 = median(umapharmony_1),
        umapharmony_2 = median(umapharmony_2)
      )
    
    library(ggrepel)
  
    plot<-plot1 +geom_label_repel(aes(label=obj.cell_type_new),size=4,color="black",fontface="bold",data = signature_ggplot_med,
                                                            point.padding = NA,label.size = NA, fill = alpha(c("white"),0.5),
                                                            segment.size=0.5,force = 1,nudge_x=0.5, nudge_y = 0,direction="y",max.overlaps=50)+
      theme(legend.position = "none")
    
  }
  if (dimention.reduction.type == "tsne") {
    if (dimention.reduction.run == T) 
      obj <- Seurat::RunTSNE(obj, reduction = "pca", dims = 1:40)
    tsne.loc <- obj@reductions$tsne@cell.embeddings
    row.names(tsne.loc) <- colnames(obj)
    signature_exp <- obj@assays$METABOLISM$score
    input.pathway <- pathway
    signature_ggplot <- data.frame(tsne.loc, t(signature_exp[input.pathway, 
    ]))
    pal <- wes_palette("Zissou1", 100, type = "continuous")
    library(ggplot2)
    plot <- ggplot(data = signature_ggplot, aes(x = tSNE_1, 
                                                y = tSNE_2, color = signature_ggplot[, 3])) + geom_point(size = size) + 
      scale_fill_gradientn(colours = pal) + scale_color_gradientn(colours = pal) + 
      labs(color = input.pathway) + xlab("tSNE 1") + ylab("tSNE 2") 
      theme(aspect.ratio = 1) + theme(panel.grid.major = element_blank(), 
                                      panel.grid.minor = element_blank(), panel.background = element_blank(), 
                                      axis.line = element_line(colour = "black"))
  }
  plot
}

DotPlot <- DotPlot.metabolism(obj = MPRNMPR_object_miMPRobe_remove_metabolism, 
                              pathway = pathway_name, 
                              phenotype = "MPR_Response2", norm = "y")


ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/INVADE-seq/aMiAD_score/代谢通路分析/DotPlot.metabolism_MPR.pdf",DotPlot,width=10,height=12)
MPRNMPR_object_miMPRobe_remove_metabolism@reductions$umap.harmony@cell.embeddings
MPRNMPR_object_miMPRobe_remove_metabolism@reductions
MPRNMPR_object_miMPRobe_remove_metabolism$cell_type_new

p1<-DimPlot.metabolismV5(obj = MPRNMPR_object_miMPRobe_remove_metabolism, 
                     pathway = c("Nicotinate and nicotinamide metabolism"), 
                     dimention.reduction.type = "umap.harmony",  
                     dimention.reduction.run =T, size = 1)

p2<-DimPlot.metabolismV5(obj = MPRNMPR_object_miMPRobe_remove_metabolism, 
                     pathway = c("Riboflavin metabolism"), 
                     dimention.reduction.type = "umap.harmony",  
                     dimention.reduction.run =T, size = 1)

p3<-DimPlot.metabolismV5(obj = MPRNMPR_object_miMPRobe_remove_metabolism, 
                     pathway = c("Porphyrin and chlorophyll metabolism"), 
                     dimention.reduction.type = "umap.harmony",  
                     dimention.reduction.run =T, size = 1)

p<-p1+p2+p3

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/INVADE-seq/aMiAD_score/代谢通路分析/Nicotinate代谢通路umap.pdf",p1,width=6,height=6)

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/INVADE-seq/aMiAD_score/代谢通路分析/Riboflavin代谢通路umap.pdf",p2,width=6,height=6)

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/INVADE-seq/aMiAD_score/代谢通路分析/Porphyrin代谢通路umap.pdf",p3,width=6,height=6)




MPRNMPR_object_miMPRobe_remove_metabolism@reductions

input.pathway<-c("Nicotinate and nicotinamide metabolism",
                 "Porphyrin and chlorophyll metabolism",
                 "Riboflavin metabolism")
MPRNMPR_object_miMPRobe_remove_metabolism$cell_type_new<- factor(MPRNMPR_object_miMPRobe_remove_metabolism$cell_type_new,levels=c('T cells', 'B cells','Plasma cells', 'Myeloid cells',  
                                                                                                                                  'Mast cells',"Endothelial cells","Fibroblasts","Smooth muscle cells",
                                                                                                                                  "Pericytes","Epithelial cells"))
DotPlot2 <- DotPlot.metabolism(obj = MPRNMPR_object_miMPRobe_remove_metabolism,
                   pathway =input.pathway,
                   phenotype = "cell_type_new", #这个参数需按需修改
                   norm = "y")

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/INVADE-seq/aMiAD_score/代谢通路分析/三条代谢通路不同细胞类型dotplot.pdf",DotPlot2,width=9,height=4)


countexp.Seurat <- MPRNMPR_object_miMPRobe_remove_metabolism
BoxPlot.metabolism(obj = countexp.Seurat,
                   pathway = input.pathway, 
                   phenotype = "cell_type_new", #这个参数需按需修改
                   ncol = 1)

###########################################

b_top30 <-c("Ralstonia",
            "Basfia",
            "Rothia",
            "Moraxella",
            "Aggregatibacter",
            "Pelomonas",
            "Bacteroides",
            "Escherichia",
            "Mycobacterium",
            "Pseudomonas",
            "Alloprevotella",
            "Haemophilus",
            "Capnocytophaga",
            "Corynebacterium",
            "Nesterenkonia",
            "Bacillus_A",
            "Halomonas",
            "Pauljensenia",
            "Veillonella",
            "Leptotrichia",
            "Haemophilus",
            "Acinetobacter",
            "Geobacillus",
            "Klebsiella",
            "Fusobacterium",
            "Porphyromonas",
            "Cutibacterium",
            "Neisseria",
            "Prevotella",
            "Streptococcus")

s16_top30<- c("Streptococcus",
              "Muribaculaceae",
              "Rothia", 
              "Bacteroides",
              "Neisseria" ,
              "Others",
              "Fusobacterium",
              "Haemophilus",
              "Prevotella", 
              "Lachnospiraceae_NK4A136_group",
              "Ralstonia",
              "Actinomyces", 
              "Alistipes",
              "Porphyromonas",
              "Lactobacillus",
              "Parabacteroides",
              "Enterococcus",
              "Granulicatella",
              "Alloprevotella",
              "Leptotrichia" ,
              "Capnocytophaga",
              "[Eubacterium]_coprostanoligenes_group",
              "Helicobacter",
              "Oscillospiraceae UCG-002", 
              "Geobacillus",
              "Lentimicrobium",
              "Clostridia_UCG-014",
              "Gemella",
              "Odoribacter",
              "Bifidobacterium" )


invade_top30<- c("Massilia",
                 "Herbaspirillum",
                 "Corynebacterium",
                 "Sphingomonas",
                 "Treponema",
                 "Acinetobacter",
                 "Achromobacter",
                 "Fusobacterium",
                 "Arthrobacter",
                 "Rothia",
                 "Lawsonella",
                 "Rhodococcus",
                 "Blattabacterium",
                 "Methylobacterium",
                 "Duganella",
                 "Dietzia",
                 "unclassified.Prevotellaceae",
                 "Streptococcus",
                 "Prevotella",
                 "Aerococcus",
                 "Planococcus",
                 "Staphylococcus",
                 "Campylobacter",
                 "Capnocytophaga",
                "Bradyrhizobium",
                "Propionibacterium",
                "Lactobacillus",
                "Actinomyces",
                "Blastococcus",
                "Belnapia")

st_bacteria<-c(
  "Pseudomonas",
  "Thiomicrorhabdus",
  "Stenotrophomonas",
  "Aeromonas",
  "Lysinibacillus",
  "Vibrio",
  "Bacillus",
  "Streptococcus",
  "Fusobacterium",
  "Proteus",
  "Desulfohalotomaculum",
  "Neisseria",
  "Weissella",
  "Klebsiella",
  "Myroides",
  "Staphylococcus",
  "Mycoplasma",
  "Capnocytophaga",
  "Endomicrobium",
  "Arsenophonus",
  "Serratia",
  "Epibacterium",
  "Leptotrichia",
  "Paraburkholderia",
  "Ureaplasma",
  "NS2b.marine.group",
  "Veillonella",
  "Haemophilus",
  "Lachnoanaerobaculum",
  "Desulfobulbus"
)

intersect(invade_top30, st_bacteria)


common_elements <- Reduce(intersect, list(b_top30, s16_top30, invade_top30, st_bacteria))


######alpha 相关性图----


alpha_16s<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/16S_结题报告/MPR_group/alpha_diversity/alpha_estimator_summary.xls",header = T, row.names=1,sep="\t")

alpha_2b<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/2bRAD_分析/MPR_NMPR_2b/alpha_diveristy/alpha_estimator_summary.xls",header = T,row.names = 1,sep="\t")

alpha_invade <-alpha 

rownames(alpha_2b)

sc_rows <- rownames(alpha_invade)[grepl("SC", rownames(alpha_invade))]
alpha_invade_sc <- alpha_invade[sc_rows, ]
alpha_2b_sorted <- alpha_2b[sub("SC", "MC", rownames(alpha_invade_sc)), ]
alpha_2b_sorted<- na.omit(alpha_2b_sorted)

alpha_16s_sorted <- alpha_16s[sub("SC", "MC", rownames(alpha_invade_sc)), ]
alpha_16s_sorted<- na.omit(alpha_16s_sorted)

alpha_invade_sc_sorted <- alpha_invade_sc[sub("MC", "SC", rownames(alpha_2b_sorted)), ]


df12 <- merge(alpha_2b_sorted, alpha_16s_sorted, by = "row.names", all = TRUE)
rownames(df12)<-rownames(alpha_invade_sc_sorted)
df_all <- merge(df12, alpha_invade_sc_sorted, by = "row.names", all = TRUE)

df_all <- df_all[, -c(1,2)]

df_all <- df_all %>%
  mutate(across(everything(), as.character))

write.csv(df_all,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/INVADE-seq/alpha线性回归/df_alpha.csv")


df_long <- df_all %>%
  pivot_longer(cols =everything(),  # 需要转换为长格式的列
               names_to = "variable",         # 新列，用来存放原列的名字
               values_to = "value")  

write.csv(df_long,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/INVADE-seq/alpha线性回归/df_long_alpha.csv")

data("iris")
plot(alpha_2b_sorted$chao1,alpha_invade_sc_sorted$chao1,lwd = 1,
     xlab = "alpha_16s_sorted_shannon",ylab = "alpha_invade_sc_sorted_Shannon")
fit <- lm(alpha_16s_sorted$observed_species~alpha_invade_sc_sorted$observed_species
          )#线性回归拟合
summary(fit)#查看计算结果


df_alpha<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/INVADE-seq/alpha线性回归/df_alpha.csv")

names(df_alpha)

A <- read.csv("F:/生信分析代码/KS-10.5合集/KS-account-codes/6-KS个性化作图（不止单细胞）/9-回归及置信区间/回归置信.csv", header = T,row.names = 1)
library(ggpubr)
library(ggplot2)
install.packages("Hmisc")
library(Hmisc)
#计算相关系数和p
cor <- rcorr(df_alpha$invade,df_alpha$X2b, type = "pearson")
pvalue <- round(cor$P[1,2],3)
corvalue <- round(cor$r[1,2],2)
names(df_alpha)
chao_df<-df_alpha %>% filter(alpha %in% "chao1")

library(reshape2)

data_long_A<-melt(A, id.vars = c("KRAS_VAF"), #需保留的不参与聚合的变量列名
                  measure.vars = 2:4,#选择需要转化的列
                  variable.name = c('sample'),#聚合变量的新列名
                  value.name = 'value')#聚合值的新列名


data_long_A<-melt(chao_df, id.vars = c("invade"), #需保留的不参与聚合的变量列名
                  measure.vars = 2:3,#选择需要转化的列
                  variable.name = c('sample'),#聚合变量的新列名
                  value.name = 'value')#聚合值的新列名

names(data_long_A)

#计算相关系数和P
a <- vector()
b <- vector()
c <- vector()


for (i in 2:ncol(chao_df)-1) {
  cor <- rcorr(chao_df$invade, chao_df[,i], type = "spearman")
  pvalue <- round(cor$P[1,2],5)
  corvalue <- round(cor$r[1,2],2)
  
  a[[i]] <- pvalue
  b[[i]] <- corvalue
  c[[i]] <- paste0('R=', b[i], ', P=', a[i])
}

cor_data <- data.frame(a,b,c)
cor_data <- cor_data[-1,]
rownames(cor_data) <- colnames(chao_df[, 2:3])
colnames(cor_data) <- c("pvalue", "corvalue", 'labe')

cor_data$y <- c(250,750)#y坐标


#作图
ggplot(data_long_A,aes(invade,value))+
  geom_point(size = 3,aes(color = sample,fill = sample))+
  geom_smooth(aes(color = sample,fill = sample),method = "lm",level = 0.95,formula = y~x,linetype = 2,alpha = 0.2)+
  scale_color_manual(values = c('#e8351e',
                                '#0f8096'))+
  geom_text(data = cor_data,
            aes(100,y,label = labe),
            size = 4,
            color = c('#e8351e',
                      '#0f8096'))+
  theme_bw()+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black",size = 1),
        axis.text.x=element_text(vjust=1, hjust = 1, size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 15, colour = 'black'))



library(ggplot2)
ggplot(chao_df, aes(x=invade, y=X2b)) + 
  geom_point(color='#852f88', size=3)+
  geom_smooth(method = "lm",level = 0.95,formula = y~x, color='#852f88') +
  ggpubr::stat_cor(method = "spearman", 
                   aes(label = paste(..rr.label.., ..p.label.., sep = "~")),
                   color = "black", geom = "label", label.x = 0.2, label.y = 0) + #这里可以直接添加R、P，可以应用到上边的ggplot线性归回
  theme_bw()+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black",size = 1),
        axis.text.x=element_text(vjust=1, hjust = 1, size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 15, colour = 'black'))

#=======================================================================
#=================================================================
  

unique(df_alpha$alpha)
chao_df<-df_alpha %>% filter(alpha %in% "observed_species")

library(reshape2)


data_long_A<-melt(chao_df, id.vars = c("invade"), #需保留的不参与聚合的变量列名
                  measure.vars = 2:3,#选择需要转化的列
                  variable.name = c('sample'),#聚合变量的新列名
                  value.name = 'value')#聚合值的新列名

names(data_long_A)

#计算相关系数和P
a <- vector()
b <- vector()
c <- vector()


for (i in 2:ncol(chao_df)-1) {
  cor <- rcorr(chao_df$invade, chao_df[,i], type = "spearman")
  pvalue <- round(cor$P[1,2],5)
  corvalue <- round(cor$r[1,2],2)
  
  a[[i]] <- pvalue
  b[[i]] <- corvalue
  c[[i]] <- paste0('R=', b[i], ', P=', a[i])
}

cor_data <- data.frame(a,b,c)
cor_data <- cor_data[-1,]
rownames(cor_data) <- colnames(chao_df[, 2:3])
colnames(cor_data) <- c("pvalue", "corvalue", 'labe')

cor_data$y <- c(300,600)#y坐标


#作图
ggplot(data_long_A,aes(invade,value))+
  geom_point(size = 3,aes(color = sample,fill = sample))+
  geom_smooth(aes(color = sample,fill = sample),method = "lm",level = 0.95,formula = y~x,linetype = 2,alpha = 0.2)+
  scale_color_manual(values = c('#e8351e',
                                '#0f8096'))+
  geom_text(data = cor_data,
            aes(3,y,label = labe),
            size = 4,
            color = c('#e8351e',
                      '#0f8096'))+
  theme_bw()+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black",size = 1),
        axis.text.x=element_text(vjust=1, hjust = 1, size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 15, colour = 'black'))


