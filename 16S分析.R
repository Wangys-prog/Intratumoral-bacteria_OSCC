
library(microeco)
library(ggplot2)
library(dplyr)

abundance<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/16S_结题报告/MPR_group/feature_tax_even.xls",sep="\t",header = T,row.names = 1,check.names = FALSE)


###处理一下species这一列
abundance$Kingdom <- gsub("k__", "", abundance$Kingdom)
abundance$Phylum<- gsub("p__", "", abundance$Phylum)
abundance$Class<- gsub("c__", "", abundance$Class)
abundance$Order<- gsub("o__", "", abundance$Order)
abundance$Family<- gsub("f__", "", abundance$Family)
abundance$Genus<- gsub("g__", "", abundance$Genus)
abundance$Species <- gsub("_g_.*", "", abundance$Species)
abundance$Species <- gsub("s__", "", abundance$Species)

names(abundance)

genus0 <- abundance[,c(1:38,44)]
names(genus0)
genus1 <- genus0 %>%
  group_by(!!sym("Genus")) %>%  # 将 V53 替换为你的实际列名
  summarise(across(1:38, sum, na.rm = TRUE))  %>% as.data.frame()
names(genus1)
rownames(genus1)[1]<-"Others"
View(genus1)
rownames(genus1)<- genus1$Genus 
genus2<- genus1[,-1] %>% as.data.frame
View(genus2)
rownames(genus2)[1]<-"Others"
genus2[is.na(genus2)] <- 0
genus2 <- genus2[rowSums(genus2) != 0, ]
genus2 <- genus2[, colSums(genus2) != 0]

genus2[] <- lapply(genus2, function(x) as.numeric(as.character(x)))
write.csv(genus2,"D:/wangys_uestc_data2/OSCC_data/oscc_process/16S_结题报告/MPR_group/genus2.csv")

###分组信息
sample <- read.table("D:/wangys_uestc_data2/OSCC_data/oscc_process/16S_结题报告/MPR_group/group_mpr.csv", header=T, sep=",",row.names = 1, comment.char="",stringsAsFactors = TRUE,quote = "")
View(sample)
sample$group_id <- rownames(sample)


#物种分类表
tax<- rownames(genus2) %>% as.data.frame()
View(tax)
colnames(tax)<- "genus"
rownames(tax) <- tax$genus

df<-microtable$new(sample_table=sample,
                   otu_table=genus2,
                   tax_table=tax,
                   auto_tidy=F)
df

rownames(sample)
df$sample_table$MPR_response2
df$sample_table$MPR_response2
write.csv(genus_df$data_abund,"D:/wangys_uestc_data2/OSCC_data/oscc_process/16S_结题报告/MPR_group/物种堆砌图.csv")


genus_df<-trans_abund$new(dataset=df,taxrank="genus",ntaxa=30,groupmean = "MPR_response2")

genus_df<-trans_abund$new(dataset=df,taxrank="genus",ntaxa=30)
genus_df$data_abund$Sample<- factor(genus_df$data_abund$Sample,levels=c("ME","Pre_NMPR","Pre_MPR","Post_MPR"))
p1<-genus_df$plot_bar(others_color="grey70",#剩余分类的填充色
                      #facet="MPR_Response2",#根据组进行分面
                      xtext_keep=F,#是否显示样本名称
                      legend_text_italic=F)+
  theme(
    panel.background = element_blank(),  # 去除面板背景
    plot.background = element_blank(),    # 去除整个图形的背景
    panel.grid.major = element_blank(),   # 去除主要网格线
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 14,angle=30,vjust = 0.5,hjust=0.5),
    axis.text.y = element_text(size = 12), 
    axis.title.y = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 13, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA,linewidth = 1)
  )

p1

library(ggsci)
col=pal_d3("category20")(20)
col2 = pal_d3("category20",alpha = 0.5)(20)
mypal=c(col,col2[-8])
df1<- genus_df$data_abund 
taxa <- genus_df$data_taxanames
df2<- df1 %>% filter(df1$Taxonomy%in% taxa) 
df2$Sample<- factor(df2$Sample,levels=c("ME","Pre_NMPR","Pre_MPR","Post_MPR"))
ppp<-ggplot()+  geom_bar(data=df2,          
                    aes(x=Sample,              
                        weight=Abundance,              
                        fill=reorder(Taxonomy ,-Abundance)),          
                    position = "fill",# 百分比堆叠图position = "fill"          
                    width=0.5)+  
  scale_fill_manual(values = mypal[-18])+  
  theme_bw()+  
  guides(fill=guide_legend(title = "Phylum",ncol = 2))+  
  labs(y = "Absolute abundance",x=NULL)+  
  theme(legend.position="right",        
                                                
        axis.title = element_text(face = "bold",                                   
                                  size = 12,colour = "black"))+  
  theme(axis.text = element_text(face = "bold", 
                                 size = 10,color="black"), 
        strip.text.x = element_text(face = "bold", 
                                    size =12,color="black"))+  
  theme(panel.grid=element_blank())+  
  theme(legend.title = element_text(face = "bold",  
                                    size =12,color="black"), 
        legend.text = element_text(face = "bold", 
                                   size =10,color="black"))

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/16S_结题报告/MPR_group/16s_MPR_Response2丰度堆砌图.pdf",ppp,width = 11,height=8)

######箱状图----

genus_df<-trans_abund$new(dataset=df,taxrank="genus",input_taxaname=c("Streptococcus","Fusobacterium"))
genus_df$data_abund$MPR_response2<- factor(genus_df$data_abund$MPR_response2,levels=c("ME","Pre_NMPR","Pre_MPR","Post_MPR"))
p3<- genus_df$plot_box(color_values=c("black","#4974a4","#4dae47","#f29600"),group="MPR_response2", 
                       position = position_dodge(0.9),
                       xtext_angle=30,
                       show_point=FALSE,
                       point_size = 6)+
  #scale_y_continuous(limits = c(0, 3))+
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
library(ggpubr)
names(genus_df$data_abund)
genus_df1<- genus_df$data_abund %>% filter(Taxonomy %in% c("Streptococcus", "Fusobacterium"))
genus_df1$Taxonomy <- factor(genus_df1$Taxonomy,levels = c("Streptococcus", "Fusobacterium"))
genus_df1$Abundance
View(genus_df1)
head(genus_df1)
library(ggplot2)
library(ggpubr)
library(rstatix)
genus_df1$Abundance <- as.numeric(genus_df1$Abundance)

df_p_val1 <- genus_df1 %>%group_by(Taxonomy) %>%
  wilcox_test(Abundance~MPR_response2) %>%
  adjust_pvalue(p.col="p",method="bonferroni") %>%
  add_significance(p.col="p.adj") %>%
  add_xy_position(x="MPR_response2",dodge=0.8)


p3<- genus_df1%>% ggplot(aes(Taxonomy,Abundance))+
  #ggplot(genus_df1, aes(x =Taxonomy, y = Abundance,fill=MPR_response2)) +
  geom_boxplot(aes(fill=MPR_response2),outlier.shape = NA) +  
  labs(title = "", x = "", y = "Relative abundance(%)") +  # 添加标题和轴标签
  stat_pvalue_manual(df_p_val1,label="p.adj.signif",hide.ns=T,
                     tip.length = 0,label.size = 5,color="black")+
  theme_minimal() +  # 使用简单主题
  scale_fill_manual(values = c("black","#4974a4","#4dae47","#f29600"))+
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

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/16S_结题报告/MPR_group/16s_物种箱状图.pdf",p3,width = 6,height=5)

####heatmap---
genus_df<-trans_abund$new(dataset=df,taxrank="genus",ntaxa=30,groupmean = "MPR_response2")
genus_df$data_abund$Sample<- factor(genus_df$data_abund$Sample,levels=c("ME","Pre_NMPR","Pre_MPR","Post_MPR"))

genus_df<-trans_abund$new(dataset=df,taxrank="genus",ntaxa=30)
genus_df$data_abund$MPR_response2<- factor(genus_df$data_abund$MPR_response2,levels = c("ME","Pre_NMPR","Pre_MPR","Post_MPR"))

p4<- genus_df$plot_heatmap(
  color_values = rev(RColorBrewer::brewer.pal(n = 11, name = "RdYlBu")),
  #facet ="MPR_response2",
  x_axis_name = NULL,
  order_x = NULL,
  xtext_size = 10,
  ytext_size = 14,
  xtext_angle = 30)
p4
ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/16S_结题报告/MPR_group/microeco/16s_heatmap_MPR.pdf",p4,width = 8,height = 7)

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/16S_结题报告/MPR_group/microeco/16s_heatmap_patient.pdf",p4,width = 10,height = 7)


###物种丰度###----
setwd("D:/wangys_uestc_data2/OSCC_data/oscc_process/16S_结题报告/")
dir()

abundance<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/16S_结题报告/MPR_group/feature_tax_even.xls",sep="\t",header = T,row.names = 1)
View(abundance)
colnames(abundance)
abundance$Kingdom
colnames(abundance)<-c("WRY_MP_post","XL_MP_post","XYH_MP_post","LXD_MP_post","QHH_MP_post",
                       "YMS_MP_post","WRY_MC_post","XL_MC_post","XYH_MC_post","LXD_MC_post","QHH_MC_post",
                       "YMS_MC_post","LXD_MP_pre","XYC_MP_pre","XL_MP_pre","XYH_MP_pre","YMS_MP_pre", 
                       "WRY_MP_pre","LM_MP_pre","HTX_MP_pre","LXD_MC_pre","XYC_MC_pre","XL_MC_pre",  
                       "XYH_MC_pre","YMS_MC_pre","WRY_MC_pre","LM_MC_pre","HTX_MC_pre","WRY_ME_post",
                       "WRY_ME_pre","YMS_ME_post","YMS_ME_pre","LXD_ME_post","XL_ME_post","XL_ME_pre",  
                       "XYH_ME_post","XYH_ME_pre","YXC_ME_pre","Kingdom","Phylum","Class","Order",      
                       "Family","Genus","Species")
###处理一下species这一列
abundance$Kingdom <- gsub("k__", "", abundance$Kingdom)
abundance$Phylum<- gsub("p__", "", abundance$Phylum)
abundance$Class<- gsub("c__", "", abundance$Class)
abundance$Order<- gsub("o__", "", abundance$Order)
abundance$Family<- gsub("f__", "", abundance$Family)
abundance$Genus<- gsub("g__", "", abundance$Genus)
abundance$Species <- gsub("_g_.*", "", abundance$Species)
abundance$Species <- gsub("s__", "", abundance$Species)

names(abundance)

abundance_mat <- abundance[1:38]
tax_name <- abundance[39:45]
names(tax_name)
head(abundance_mat)
head(tax_name)

## 样本分组数据
group <- read.table("D:/wangys_uestc_data2/OSCC_data/oscc_process/16S_结题报告/MPR_group/group_mpr.csv", header=T, sep=",",row.names = 1, comment.char="",stringsAsFactors = TRUE,quote = "")
names(group)

#####species相对丰度堆砌图----
species<- tax_name$Species
####筛选MPR 分组所需要的样本
sample_mpr<- intersect(rownames(group), colnames(abundance_mat))
abundance_mat_mpr<- abundance_mat[,sample_mpr]
####将reads数目转化成相对丰度
head(abundance_mat_mpr)
# 计算每个样本的总计数
total_counts <- colSums(abundance_mat_mpr)
# 计算相对丰度
relative_abundance <- sweep(abundance_mat_mpr, 2, total_counts, FUN = "/")
# 乘以100可以得到百分比丰度
relative_abundance <- relative_abundance * 100

# 查看结果
print(relative_abundance)

####合并species 物种丰度
relative_abundance$species<- tax_name$Species
names(relative_abundance)
relative_abundance$species <- replace(relative_abundance$species, relative_abundance$species == "", "Others")
relative_abundance$species <- replace(relative_abundance$species, relative_abundance$species == "other", "Others")
relative_abundance$species <- replace(relative_abundance$species, relative_abundance$species == "uncultured_bacterium", "Others")
relative_abundance$species <- replace(relative_abundance$species, relative_abundance$species == "uncultured_organism", "Others")
relative_abundance$species <- replace(relative_abundance$species, relative_abundance$species == "unidentified", "Others")
relative_abundance$species <- replace(relative_abundance$species, relative_abundance$species == "gut_metagenome", "Others")
relative_abundance$species <- replace(relative_abundance$species, relative_abundance$species == "metagenome", "Others")

relative_abundance_aggregated <- aggregate(. ~ species, data = relative_abundance[, -1], FUN = sum)
unique(relative_abundance_aggregated$species)

rownames(relative_abundance_aggregated) <- relative_abundance_aggregated$species

names(relative_abundance_aggregated)
relative_abundance_aggregated<-relative_abundance_aggregated[,-1]

relative_abundance_aggregated$sum <- rowSums(relative_abundance_aggregated)
names(relative_abundance_aggregated)

####降序排序
relative_abundance_aggregated<- relative_abundance_aggregated[order(relative_abundance_aggregated$sum,decreasing = TRUE),]

write.csv(relative_abundance_aggregated,"D:/wangys_uestc_data2/OSCC_data/oscc_process/16S_结题报告/MPR_group/16_物种堆砌图/species_abundance.csv")

View(relative_abundance_aggregated)

#relative_abundance_aggregated<- relative_abundance_aggregated[rownames(relative_abundance_aggregated) != "Others", ]

species_top50<- relative_abundance_aggregated[1:50,-ncol(relative_abundance_aggregated)]

View(species_top50)
species_top50["Others",] <- 100-colSums(species_top50[-1,])

library(reshape2)
species_top50$species <- factor(rownames(species_top50),levels = rev(rownames(species_top50)))
species_top50<- melt(species_top50,id="species")
names(species_top50)[2]<- "samples"
View(species_top50)

##补全分组信息###
head(group)
species_top50$group1<-group[match(species_top50$samples, rownames(group)),3]
species_top50$treatments<-group[match(species_top50$samples, rownames(group)),4]
species_top50$MPR_response1<-group[match(species_top50$samples, rownames(group)),5]
species_top50$MPR_response2<-group[match(species_top50$samples, rownames(group)),6]
species_top50$CR_response1<-group[match(species_top50$samples, rownames(group)),7]
species_top50$CR_response2<-group[match(species_top50$samples, rownames(group)),8]
View(species_top50)

new_data = species_top50[species_top50$group1=="ME"&species_top50$value>0, ]
me_species <- new_data$species

species_top50_filter<-  species_top50[-which(species_top50$species %in% me_species),]

species_top50_filter<-  species_top50_filter[-which(species_top50_filter$group1 == "ME"),]
View(species_top50_filter)

###调整坐标轴位置
species_top50_filter$MPR_response2 <- factor(species_top50_filter$MPR_response2,levels = c('Pre_NMPR','Pre_MPR','Post_MPR'))

species_top50$MPR_response2 <- factor(species_top50$MPR_response2,levels = c('Pre_NMPR','Pre_MPR','Post_MPR'))

color2 <- c("#FF6347","#ffc556","#20B2AA")
# 2.1.2 颜色
library(ggsci)
col=pal_d3("category20")(20)
col2 = pal_d3("category20",alpha = 0.5)(20)
mypal=c(col,col2[-8])
###随机生成颜色
library(randomcoloR)
palette <- randomColor(count = 60)  #随机生成60种颜色，其实里面有重复的
palette <- distinctColorPalette(51) #差异明显的60种

# 2.1.3  物种组成堆叠柱形图-绝对丰度
library(ggplot2)

species_filter<-ggplot()+geom_bar(data=species_top50_filter,
                                  aes(x=MPR_response2,
                                      weight=value,
                                      fill=reorder(species,-value)),
                                  position = "fill", # ggplot2会自行计算相对丰度，无需提前计算。
                                  width=0.5)+
  #facet_grid(.~depth)+
  scale_fill_manual(values = mypal[-18])+ # 颜色与绝对丰度堆叠柱形图保持一致
  #scale_fill_manual(values = palette)+
  scale_y_continuous(#expand=c(0,0), #设置横坐标轴紧挨柱形图
    name = "Relative abundance (%)",
    limits = c(0,1),
    breaks = seq(0,1,0.25),
    labels = paste(seq(0,100,25),"%")
  )+
  theme_bw()+
  guides(fill=guide_legend(title = "Species",ncol = 1))+
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
species_filter

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/16S_结题报告/MPR_group/16_物种堆砌图/species_filter_abundance_top100.jpg",species_filter,width =9,height =7)
ggsave("species_filter_abundance_top100.pdf",species_filter,width =9,height =7)


#####################
#######################
library(ggplot2)
species_top50$MPR_response2 <- factor(species_top50$MPR_response2,levels = c("ME",'Pre_NMPR','Pre_MPR','Post_MPR'))

######
unique(species_top50$MPR_response2)
species_abun<-ggplot()+geom_bar(data=species_top50,
                                aes(x=MPR_response2,
                                    weight=value,
                                    fill=reorder(species,-value)),
                                position = "fill", # ggplot2会自行计算相对丰度，无需提前计算。
                                width=0.5)+
  #facet_grid(.~depth)+
  scale_fill_manual(values = palette)+ # 颜色与绝对丰度堆叠柱形图保持一致
  scale_y_continuous(#expand=c(0,0), #设置横坐标轴紧挨柱形图
    name = "Relative abundance (%)",
    limits = c(0,1),
    breaks = seq(0,1,0.25),
    labels = paste(seq(0,100,25),"%")
  )+
  theme_bw()+
  guides(fill=guide_legend(title = "Species",ncol = 2))+
  labs(x=NULL)+
  theme(panel.border = element_rect(color = "black", size = 1, fill = NA),
        legend.position="right",
        axis.title = element_text(face = "bold", 
                                  size = 20,colour = "black"))+
  theme(axis.text = element_text(face = "bold", 
                                 size = 20,color="black"),
        strip.text.x = element_text(face = "bold", 
                                    size =20,color="black"))+
  theme(panel.grid=element_blank())+
  theme(legend.title = element_text(face = "bold", 
                                    size =20,color="black"),
        legend.text = element_text(face = "bold", 
                                   size =16,color="black"))
species_abun

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/16S_结题报告/MPR_group/16_物种堆砌图/species_abundance_top50.jpg",species_abun,width =15,height =9)
ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/16S_结题报告/MPR_group/16_物种堆砌图/species_abundance_top50.pdf",species_abun,width =15,height =9)

#####species###----


#####genus相对丰度堆砌图----
names(relative_abundance)

relative_abundance <- sweep(abundance_mat_mpr, 2, total_counts, FUN = "/")
# 乘以100可以得到百分比丰度
relative_abundance <- relative_abundance * 100
# 查看结果
print(relative_abundance)

relative_abundance <- relative_abundance[,-24]

relative_abundance$genus<- tax_name$Genus

relative_abundance$genus <- replace(relative_abundance$genus, relative_abundance$genus == "", "Others")


relative_abundance_aggregated <- aggregate(. ~ genus, data = relative_abundance[, -1], FUN = sum)
unique(relative_abundance_aggregated$genus)

rownames(relative_abundance_aggregated) <- relative_abundance_aggregated$genus

names(relative_abundance_aggregated)

relative_abundance_aggregated<-relative_abundance_aggregated[,-1]

relative_abundance_aggregated$sum <- rowSums(relative_abundance_aggregated)
names(relative_abundance_aggregated)
head
####降序排序
relative_abundance_aggregated<- relative_abundance_aggregated[order(relative_abundance_aggregated$sum,decreasing = TRUE),]

write.csv(relative_abundance_aggregated,"D:/wangys_uestc_data2/OSCC_data/oscc_process/16S_结题报告/MPR_group/16_物种堆砌图/genus_abundance.csv")

View(relative_abundance_aggregated)

#relative_abundance_aggregated<- relative_abundance_aggregated[rownames(relative_abundance_aggregated) != "Others", ]

genus_top50<- relative_abundance_aggregated[1:50,-ncol(relative_abundance_aggregated)]

View(genus_top50)
genus_top50["Others",] <- 100-colSums(genus_top50)


library(reshape2)
genus_top50$genus<- factor(rownames(genus_top50),levels = rev(rownames(genus_top50)))
genus_top50<- melt(genus_top50,id="genus")
names(genus_top50)[2]<- "samples"
View(genus_top50)

##补全分组信息###
head(group)
genus_top50$group1<-group[match(genus_top50$samples, rownames(group)),3]
genus_top50$treatments<-group[match(genus_top50$samples, rownames(group)),4]
genus_top50$MPR_response1<-group[match(genus_top50$samples, rownames(group)),5]
genus_top50$MPR_response2<-group[match(genus_top50$samples, rownames(group)),6]
genus_top50$CR_response1<-group[match(genus_top50$samples, rownames(group)),7]
genus_top50$CR_response2<-group[match(genus_top50$samples, rownames(group)),8]


new_data = genus_top100[genus_top100$group1=="ME"&genus_top100$value>0, ]
me_genus <- new_data$genus

genus_top100_filter<- genus_top100[-which(genus_top100$genus %in% me_genus),]

genus_top100_filter<-  genus_top100_filter[-which(genus_top100_filter$group1 == "ME"),]
View(genus_top100_filter)

###调整坐标轴位置
genus_top100_filter$MPR_response2 <- factor(genus_top100_filter$MPR_response2,levels = c('Pre_NMPR','Pre_MPR','Post_MPR'))

# 2.1.2 颜色
library(ggsci)
col=pal_d3("category20")(20)
col2 = pal_d3("category20",alpha = 0.5)(20)
mypal=c(col,col2[-8])
# 2.1.3  物种组成堆叠柱形图-绝对丰度
library(ggplot2)

genus_filter<-ggplot()+geom_bar(data=genus_top100_filter,
                                aes(x=MPR_response2,
                                    weight=value,
                                    fill=reorder(genus,-value)),
                                position = "fill", # ggplot2会自行计算相对丰度，无需提前计算。
                                width=0.5)+
  #facet_grid(.~depth)+
  scale_fill_manual(values = mypal[-18])+ # 颜色与绝对丰度堆叠柱形图保持一致
  scale_y_continuous(#expand=c(0,0), #设置横坐标轴紧挨柱形图
    name = "Relative abundance (%)",
    limits = c(0,1),
    breaks = seq(0,1,0.25),
    labels = paste(seq(0,100,25),"%")
  )+
  theme_bw()+
  guides(fill=guide_legend(title = "Genus",ncol = 1))+
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
genus_filter

ggsave("genus_filter_abundance_top200.jpg",genus_filter,width =9,height =7)
ggsave("genus_filter_abundance_top200.pdf",genus_filter,width =9,height =7)

library(randomcoloR)
palette <- distinctColorPalette(51)

genus_top50$MPR_response2 <- factor(genus_top50$MPR_response2,levels = c("ME",'Pre_NMPR','Pre_MPR','Post_MPR'))

genus_abun<-ggplot()+geom_bar(data=genus_top50,
                              aes(x=MPR_response2,
                                  weight=value,
                                  fill=reorder(genus,-value)),
                              position = "fill", # ggplot2会自行计算相对丰度，无需提前计算。
                              width=0.5)+
  #facet_grid(.~depth)+
  scale_fill_manual(values = palette)+ # 颜色与绝对丰度堆叠柱形图保持一致
  scale_y_continuous(#expand=c(0,0), #设置横坐标轴紧挨柱形图
    name = "Relative abundance (%)",
    limits = c(0,1),
    breaks = seq(0,1,0.25),
    labels = paste(seq(0,100,25),"%")
  )+
  theme_bw()+
  guides(fill=guide_legend(title = "Genus",ncol = 2))+
  labs(x=NULL)+
  theme(panel.border = element_rect(color = "black", size = 1, fill = NA),
        legend.position="right",
        axis.title = element_text(face = "bold", 
                                  size = 18,colour = "black"))+
  theme(axis.text = element_text(face = "bold", 
                                 size = 18,color="black"),
        strip.text.x = element_text(face = "bold", 
                                    size =18,color="black"))+
  theme(panel.grid=element_blank())+
  theme(legend.title = element_text(face = "bold", 
                                    size =15,color="black"),
        legend.text = element_text(face = "bold", 
                                   size =15,color="black"))
genus_abun

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/16S_结题报告/MPR_group/16_物种堆砌图/genus_abundance_top50.jpg",genus_abun,width =15,height =9)
ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/16S_结题报告/MPR_group/16_物种堆砌图/genus_abundance_top50.pdf",genus_abun,width =15,height =9)

####

####tax_circlize()函数绘制弦图，方法和参数同tax_stackplot()。----

#####使用R的circlize包绘制样本-物种丰度关联弦状图----
###https://www.jianshu.com/p/11c997dee3f7
###----
relative_abundance <- sweep(abundance_mat_mpr, 2, total_counts, FUN = "/")
# 乘以100可以得到百分比丰度
relative_abundance <- relative_abundance * 100
# 查看结果
print(relative_abundance)

####合并species 物种丰度
relative_abundance$species<- tax_name$Species
library(circlize)
library(tidyverse)
?circlize
unique(species_top50$species)

###species 弦图----
library(dplyr)
names(species_top50)
species_top50<-species_top50 %>% select("species","MPR_response2","value")
names(species_top50)
species <- as.data.frame(unique(species_top50$species))

colnames(species) <- "species_names"
rownames(species) <- species$species_names

View(species_top50)
species_top50 <- species_top50 %>% filter(!species_top50$species =="Others")

# 根据 MPR_response2 分组，然后选择每个分组中 value 值最大的前 10 个物种
top_species_per_group <- species_top50 %>%
  group_by(MPR_response2) %>%
  top_n(10, wt = value) %>%
  ungroup()

top_species_per_group$MPR_response2<- factor(top_species_per_group$MPR_response2,levels=c("ME","Pre_NMPR","Pre_MPR","Post_MPR"))

color=NULL
color[rownames(top_species_per_group)] =distinctColorPalette(4)
color[c("Pre_NMPR","Pre_MPR","Post_MPR")]<-c("#4974a4","#4dae47","#f29600")


circos.clear()

pdf(file="D:/wangys_uestc_data2/OSCC_data/oscc_process/16S_结题报告/MPR_group/物种丰度关联弦状图/species_5.pdf", width=14, height=8, pointsize=10)
library(circlize)
chordDiagram(top_species_per_group,
             grid.col =color,#颜色
             annotationTrack = "grid",
             transparency = 0.2,#透明度
             link.lwd = 0.00001,#线条宽度
             link.lty = 1,    # 线路类型
             link.border = 0,#边框颜色
             directional = -1,#表示线条的方向，0代表没有方向，1代表正向，-1代表反向，2代表双向
             diffHeight = mm_h(3),#外圈和中间连线的间隔
             direction.type = c("diffHeight","arrows"), #线条是否带有箭头
             link.arr.type = "big.arrow",#箭头类型
             annotationTrackHeight = c(0.04, 0.1))#网格高度


circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
              facing = "clockwise", niceFacing = T, adj = c(-0.5, 0.5), cex = 0.8)
  circos.axis(h = "top", labels.cex = 0.4,labels.niceFacing = F, labels.pos.adjust =F)
}, bg.border = NA)


legend("right",pch=20,legend=top_species_per_group$species,
       col=color[top_species_per_group$species],bty="n",
       cex=1,pt.cex=3,border="black")

dev.off()


####弦图另一种画法----
library(tidyverse)
library(circlize)
library(cowplot)
library(grid)
library(ggplotify)
library(MetBrewer)
install.packages("cropcircles")
library(cropcircles) # Version:0.2.4
library(ggimage)
library(magick)

head(species_top50)
library(dplyr)
# 根据 MPR_response2 分组，然后选择每个分组中 value 值最大的前 10 个物种
top_species_per_group <- species_top50 %>%
  group_by(MPR_response2) %>%
  top_n(10, wt = value) %>%
  ungroup()

top_species_per_group$MPR_response2<- factor(top_species_per_group$MPR_response2,levels=c("ME","Pre_NMPR","Pre_MPR","Post_MPR"))

pal <- colorRampPalette(pal)(nrow(top_species_per_group))

grid_colors <- setNames(pal, rownames(top_species_per_group))
# 绘制 chord diagram

chordDiagram(top_species_per_group, grid.col = grid_colors, annotationTrack = "grid")

# 使用 circos.trackPlotRegion 来手动添加扇区标签
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  sector.index = get.cell.meta.data("sector.index")
  circos.text(
    x = mean(get.cell.meta.data("xlim")),
    y = get.cell.meta.data("ylim")[1] + 2,  # 调整 y 位置，避免重叠
    labels = sector.index,
    facing = "bending.inside",  # 文字向内弯曲
    niceFacing = TRUE,          # 自动调整文字角度
    cex = 0.4,                  # 字体大小
    adj = c(0.5, 0.5)           # 文本对齐方式，居中
  )
}, bg.border = NA)  # 不显示背景边框


######################3
#######################

genus_top100_filter<-genus_top100_filter %>% select("genus","MPR_response2","value")

genus <- as.data.frame(unique(genus_top100_filter$genus))

colnames(genus) <- "genus_names"
rownames(genus) <- genus$genus_names

color=NULL
color[rownames(genus)] =distinctColorPalette(7)
color[c("Pre_NMPR","Pre_MPR","Post_MPR")]<-c("blue","red","green")


circos.clear()

pdf(file="genus.pdf", width=14, height=8, pointsize=10)

chordDiagram(genus_top100_filter,
             grid.col =color,#颜色
             annotationTrack = "grid",
             transparency = 0.2,#透明度
             link.lwd = 0.00001,#线条宽度
             link.lty = 1,    # 线路类型
             link.border = 0,#边框颜色
             directional = -1,#表示线条的方向，0代表没有方向，1代表正向，-1代表反向，2代表双向
             diffHeight = mm_h(3),#外圈和中间连线的间隔
             direction.type = c("diffHeight","arrows"), #线条是否带有箭头
             link.arr.type = "big.arrow",#箭头类型
             annotationTrackHeight = c(0.04, 0.1))#网格高度


circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
              facing = "clockwise", niceFacing = T, adj = c(-0.5, 0.5), cex = 0.8)
  circos.axis(h = "top", labels.cex = 0.4,labels.niceFacing = F, labels.pos.adjust =F)
}, bg.border = NA)


legend("right",pch=20,legend=rownames(genus),
       col=color[rownames(genus)],bty="n",
       cex=1,pt.cex=3,border="black")

dev.off()

#######alpha_diversity###----
setwd("../alpha_diveristy")
alpha<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/2bRAD_分析/MPR_NMPR_2b/alpha_diveristy/alpha_estimator_summary.xls",header = T,row.names = 1)
## 样本分组数据
group <- read.table("D:/wangys_uestc_data2/OSCC_data/oscc_process/2bRAD_分析/MPR_NMPR_2b/alpha_diveristy/group_mpr.csv", header=T, sep=",",row.names = 1, comment.char="",stringsAsFactors = TRUE,quote = "")

alpha_MPR1<- alpha_MPR[match(rownames(group), rownames(alpha_MPR)),]

View(alpha_MPR1)
View(group)
alpha_MPR_group <- cbind(alpha_MPR1,group)
View(alpha_MPR_group)

alpha_MPR_group_long$MPR_response2<-factor(alpha_MPR_group_long$MPR_response2,levels = c("ME",'Pre_NMPR','Pre_MPR','Post_MPR'))

color<- c("#4974a4","#4dae47","#f29600")

library(tidyr)


alpha_MPR_group_long <- tidyr::gather(alpha_MPR_group,alpha,alpha_value,Shannon,Simpson,Chao1)

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

df1 <- PMCMR_compare1(alpha_MPR_group_long,'alpha','MPR_response2','alpha_value')


View(alpha_MPR_group_long)
library(ggplot2)
library(ggsignif)
library(ggpubr)
total<-ggplot(alpha_MPR_group_long, aes(x = MPR_response2, y = alpha_value,color=MPR_response2)) +
  stat_boxplot(geom = "errorbar", linewidth=1)+
  geom_boxplot(linewidth=1)+
  scale_color_manual(values =c("black","#4974a4","#4dae47","#f29600"))+
  #geom_jitter(width = 0.1,alpha =0.5,size=2)+
  geom_text(data=df1,aes(x=MPR_response2,y=mean+1.5*std,label=Letters))+
  facet_grid(alpha~., scales = "free_y") + #按variable纵向分面
  labs(title = "", x = "", y = "")+
  theme_bw()+
  theme(panel.background = element_blank(),
        panel.grid = element_blank(), 
        strip.text.y = element_text(size = 12, colour = "black",face = "bold"), #分面标题文字设置
        plot.title = element_text(hjust = 0.5,size=13,face = "bold"),
        axis.text.x = element_text(size=14,colour ="black",face = "bold",angle = 10,vjust=-0.001),
        axis.text.y = element_text(size=15,colour ="black"),
        legend.position = "none")

total
ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/2bRAD_分析/MPR_NMPR_2b/alpha_diveristy/alpha_diversity_MPR.jpg",total, height = 10,width = 5)

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/2bRAD_分析/MPR_NMPR_2b/alpha_diveristy/alpha_diversity_MPR.pdf",total, height = 10,width = 5)


#######显著性差异----

# library(tidyverse)
# library(ggprism)
# dat <- read.delim('./vegan.txt')#读入α多样性数据
# head(dat, n = 3)
# design <- read.delim('./metadata.txt')#读入试验设计文件
# head(design, n = 3)
# dat <- design %>%
#   dplyr::select(SampleID,Group) %>%
#   inner_join(dat,by='SampleID')#按照SampleID内连接文件
# head(dat, n = 3)
# dat <- gather(dat,alpha,v,-(SampleID:Group))#宽数据变长数据
# head(dat, n = 3)
# 
# library(pgirmess)
# # 1 -----------------------------------------------------------------------
# kruskalmc_compare1 <- function(data,group,compare,value){
#   library(multcompView)
#   library(pgirmess)##多组两两比较函数用到的R包
#   library(multcomp)
#   a <- data.frame(stringsAsFactors = F)
#   type <- unique(data[,group])
#   for (i in type)
#   {
#     g1=compare
#     sub_dat <- data[data[,group]==i,]
#     names(sub_dat)[names(sub_dat)==compare] <- 'g1'
#     names(sub_dat)[names(sub_dat)==value] <- 'value'
#     sub_dat$g1 <- factor(sub_dat$g1)
#     options(warn = -1)
#     
#     k <- kruskalmc(value ~ g1, data=sub_dat, probs=0.05)
#     dif <- k$dif.com[['difference']]
#     names(dif) <- rownames(k$dif.com)
#     difL <- multcompLetters(dif)
#     label <- data.frame(difL['Letters'], stringsAsFactors = FALSE)
#     label$compare = rownames(label)
#     label$type <- i
#     
#     mean_sd <- merge(aggregate(sub_dat[['value']],by=list(sub_dat[,'g1']),FUN=sd),
#                      aggregate(sub_dat[['value']],by=list(sub_dat[,'g1']),FUN=mean),by="Group.1"
#     )
#     names(mean_sd) <- c('compare','std','mean')
#     a <- rbind(a,merge(mean_sd,label,by='compare'))
#   }
#   names(a) <- c(compare,'std','mean','Letters',group)
#   return(a)
# }
# # 2 -----------------------------------------------------------------------
# PMCMR_compare1 <- function(data,group,compare,value){
#   library(multcompView)
#   library(multcomp)
#   library(PMCMRplus)
#   #library(PMCMR)##多组两两比较函数用到的R包
#   a <- data.frame(stringsAsFactors = F)
#   type <- unique(data[,group])
#   for (i in type)
#   {
#     g1=compare
#     sub_dat <- data[data[,group]==i,]
#     names(sub_dat)[names(sub_dat)==compare] <- 'g1'
#     names(sub_dat)[names(sub_dat)==value] <- 'value'
#     sub_dat$g1 <- factor(sub_dat$g1)
#     options(warn = -1)
#     
#     k <- PMCMRplus::kwAllPairsNemenyiTest(value ~ g1,data=sub_dat)
#     n <- as.data.frame(k$p.value)
#     h <- n %>%
#       mutate(compare=rownames(n)) %>%
#       gather(group,p,-compare,na.rm = TRUE) %>%
#       unite(compare,group,col="G",sep="-")
#     dif <- h$p
#     names(dif) <- h$G
#     difL <- multcompLetters(dif)
#     K.labels <- data.frame(difL['Letters'], stringsAsFactors = FALSE)
#     K.labels$compare = rownames(K.labels)
#     K.labels$type <- i
#     
#     mean_sd <- merge(aggregate(sub_dat[['value']],by=list(sub_dat[,'g1']),FUN=sd),
#                      aggregate(sub_dat[['value']],by=list(sub_dat[,'g1']),FUN=mean),by="Group.1"
#     )
#     names(mean_sd) <- c('compare','std','mean')
#     a <- rbind(a,merge(mean_sd,K.labels,by='compare'))
#   }
#   names(a) <- c(compare,'std','mean','Letters',group)
#   return(a)
# }
# 
# # 3 -----------------------------------------------------------------------
# nparcomp_compare1 <- function(data,group,compare,value){
#   library(nparcomp)##多组两两比较函数用到的R包
#   library(multcompView)
#   library(do)
#   a <- data.frame(stringsAsFactors = F)
#   type <- unique(data[,group])
#   for (i in type)
#   {
#     g1=compare
#     sub_dat <- data[data[,group]==i,]
#     names(sub_dat)[names(sub_dat)==compare] <- 'g1' 
#     names(sub_dat)[names(sub_dat)==value] <- 'value' 
#     sub_dat$g1 <- factor(sub_dat$g1)
#     
#     k <- nparcomp::nparcomp(value ~ g1,data=sub_dat, alternative = "two.sided")
#     k$Analysis
#     b <- k$Analysis
#     dif <- b$p.Value
#     names(dif) <- b$Comparison
#     difname <- names(dif) %>% 
#       Replace(from = ' , ',to='-') %>% 
#       Replace(from=c('p\\(','\\)'),to='')#正则表达式替换
#     temp_name <- data.frame(tn=difname) %>% 
#       separate(col = 'tn',sep = '-',into = c('n1','n2'))
#     difname <- paste(temp_name$n2,temp_name$n1,sep = '-') %>% 
#       Replace(from = ' - ',to='-')#正则表达式替换
#     names(dif) <- difname
#     difL <- multcompLetters(dif)
#     labels <- data.frame(difL['Letters'], stringsAsFactors = FALSE)
#     labels$compare = rownames(labels)
#     labels$type <- i
#     
#     mean_sd <- merge(aggregate(sub_dat[['value']],by=list(sub_dat[,'g1']),FUN=sd),
#                      aggregate(sub_dat[['value']],by=list(sub_dat[,'g1']),FUN=mean),by="Group.1"
#     )
#     names(mean_sd) <- c('compare','std','mean')
#     
#     a <- rbind(a,merge(mean_sd,labels,by='compare'))
#   }
#   names(a) <- c(compare,'std','mean','Letters',group)
#   return(a)
# }
# 
# df1 <- kruskalmc_compare1(dat,'alpha','Group','v')
# 
# df1 <- PMCMR_compare1(dat,'alpha','Group','v')
# df1 <- nparcomp_compare1(dat,'alpha','Group','v')
# 
# p1 = ggplot(dat)+geom_boxplot(aes(x=Group,y=v,fill=Group))+
#   geom_text(data=df1,aes(x=Group,y=mean+1.3*std,label=Letters))+
#   facet_wrap(~alpha,scales = "free_y")+ labs(x='',y='Alpha-diversity')+
#   ggprism::theme_prism()+theme(axis.text.x = element_text(angle = 45))
# p1
# 
# 
# alpha_MPR_group_long <- tidyr::gather(alpha_MPR_group,alpha,alpha_value,shannon,simpson,chao1)
# 
# df1 <- PMCMR_compare1(alpha_MPR_group_long,'alpha','MPR_response2','alpha_value')
# 
# p1 = ggplot(alpha_MPR_group_long)+geom_boxplot(aes(x=MPR_response2,y=alpha_value,fill=MPR_response2))+
#   geom_text(data=df1,aes(x=MPR_response2,y=mean+1.3*std,label=Letters))+
#   scale_color_manual(values =c('black',"blue","red","green"))+
#   facet_wrap(~alpha,scales = "free_y")+ labs(x='',y='Alpha-diversity')+
#   ggprism::theme_prism()+theme(axis.text.x = element_text(angle = 45))
# p1
# 
#######显著性差异----

####过滤掉环境因素计算alpha 多样性####
library(vegan)
library(ggplot2)

abundance_table<-  read.csv("Abundance.filtered.anno.csv", header=T, sep=",", row.names=1, stringsAsFactors = FALSE)

Transp_abundance_table <- t(abundance_table)

group <- read.table("group_mpr.csv", header=T, sep=",",row.names = 1, comment.char="",stringsAsFactors = TRUE,quote = "")
names(group)
Transp_abundance_table<- Transp_abundance_table[match(rownames(group), rownames(Transp_abundance_table)),]

Transp_abundance_table<- as.data.frame(Transp_abundance_table)

View(Transp_abundance_table)

Transp_abundance_table$group1<-group[match(rownames(Transp_abundance_table), rownames(group)),3]



genus_top100_filter<- genus_top100[-which(genus_top100$genus %in% me_genus),]


Alpha_diversity_index <- function(x, tree = NULL, base = exp(1)) {
  est <- estimateR(x)
  Obs <-  est[1, ]
  Shannon <- diversity(x, index = 'shannon', base = base)
  Simpson <- diversity(x, index = 'simpson')
  Pielou <- Shannon / log(Obs, base)
  goods_coverage <- 1 - rowSums(x == 1) / rowSums(x)
  
  result <- rbind(est, Shannon, Simpson,
                  Pielou, goods_coverage)
  if (!is.null(tree)) {
    Pd <- pd(x, tree, include.root = FALSE)[1]
    Pd <- t(Pd)
    result <- rbind(result, Pd)
  }
  result <- as.data.frame(t(result))
  return(result)
}

alpha_diversity <- Alpha_diversity_index(Transp_otu)
alpha_diversity[c(1:3), ]





####beta-diversity----

library(vegan)

library(ggplot2)

abundance<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/2bRAD_分析/MPR_NMPR_2b/beta-diversity/Abundance.filtered.anno.csv",header = T)

group <- read.table("D:/wangys_uestc_data2/OSCC_data/oscc_process/2bRAD_分析/MPR_NMPR_2b/beta-diversity/group_mpr.csv", header=T, sep=",",row.names = 1, comment.char="",stringsAsFactors = TRUE,quote = "")

abundance1<- abundance[,match(rownames(group), colnames(abundance))]

View(abundance1)

library(vegan)
library(ape)
abundnce_dist<- vegdist(t(abundance1), method="bray")
pcoa <- cmdscale(abundnce_dist, k=3, eig=T)#高维度数据降低到低维度数据的操作
pcoa_points <- as.data.frame(pcoa$points)#提取三个维度的PCoA值
sum_eig <- sum(pcoa$eig)#提取特征值作为
eig_percent <- round(pcoa$eig/sum_eig*100,1)
colnames(pcoa_points) <- paste0("PCoA", 1:2)
pcoa_result <- cbind(pcoa_points,t(abundance1))#PCoA值和OTU合并
View(pcoa_result)
abundance1_t <- as.data.frame(t(abundance1))

MPR_response2 <- group$MPR_response2#提取Moisture分组信息

div <- adonis2(abundance1_t ~ MPR_response2, data = abundance1_t, permutations = 999, method="bray")#adonis
adonis <- paste0("adonis R2: ",round(div$R2,2), "; P-value: ", div$`Pr(>F)`)#提取p值，组成一个字符串

group$MPR_response2<-factor(group$MPR_response2,levels = c("ME",'Pre_NMPR','Pre_MPR','Post_MPR'))

# matched_cols <- match(rownames(group), rownames(pcoa_result))
# matched_cols <- matched_cols[!is.na(matched_cols)]
# pcoa_result$group <- group[matched_cols,][,6]

#pcoa_result$group <- group[match(rownames(group),rownames(pcoa_result)),6]

p<-ggplot(pcoa_result, aes(x=PCoA1, y=PCoA2, color=group$MPR_response2)) +
  scale_color_manual(values =c("black","#4974a4","#4dae47","#f29600"))+
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

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/2bRAD_分析/MPR_NMPR_2b/beta-diversity/MPR_beta_diversity2.jpg",p,height = 7,width = 9)
ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/2bRAD_分析/MPR_NMPR_2b/beta-diversity/MPR_beta_diversity2.pdf",p,height = 7,width = 9)


####adonis分析####
pairwise.adonis <-function(x,factors, sim.method, p.adjust.m)
{library(vegan)
  co = as.matrix(combn(unique(factors),2))
  pairs = c()
  F.Model =c()
  R2 = c()
  p.value = c()
  for(elem in 1:ncol(co)){
    ad = adonis(x[factors %in%c(as.character(co[1,elem]),as.character(co[2,elem])),] ~
                  factors[factors %in%c(as.character(co[1,elem]),as.character(co[2,elem]))] , method =sim.method);
    pairs =c(pairs,paste(co[1,elem],'vs',co[2,elem]));
    F.Model =c(F.Model,ad$aov.tab[1,4]);
    R2 = c(R2,ad$aov.tab[1,5]);
    p.value = c(p.value,ad$aov.tab[1,6])}
  p.adjusted =p.adjust(p.value,method=p.adjust.m)
  pairw.res = data.frame(pairs,F.Model,R2,p.value,p.adjusted)
  return(pairw.res)
}

pairwise.adonis(abundance1_t, group$MPR_response2, sim.method="bray", p.adjust.m= "BH")


####共有特有微生物----


#####联合INVADE-seq分析----
setwd("D:/wangys_uestc_data2/OSCC_data/oscc_process/2bRAD_分析/MPR_NMPR_2b/INVADEseq联合top30弦图")
dir()
invade_top <-  read.csv("./umi_16s_MPR_heatmap_top30_scaled.csv")
View(invade_top)
genus_top30<-colnames(invade_top)[c(-1,-2,-3)] %>% as.data.frame() %>% rename(genus_id=".")
View(genus_top30)

genus_2b<- read.csv("Abundance.filtered.anno.xls",sep="\t",header = T)
colnames(genus_2b)
genus_2b<-genus_2b[,-c(1:5,7:9)]
View(genus_2b)

## 样本分组数据
group <- read.table("group_mpr.csv", header=T, sep=",",row.names = 1, comment.char="",stringsAsFactors = TRUE,quote = "")
names(group)

genus<- genus_2b$Genus
genus_2b<- genus_2b[,match(rownames(group), colnames(genus_2b))]
View(genus_2b)
genus_2b$Genus <- genus

genus_2b_top30_2 <- genus_2b%>% 
  filter(Genus %in% genus_top30$genus_id) %>% group_by(Genus) %>% 
  summarise_at(vars(matches(rownames(group))), funs(sum))%>% 
  as.data.frame()

View(genus_2b_top30_2)

genus_2b_top30_2<- genus_2b_top30_2[match(genus_top30$genus_id,genus_2b_top30_2$Genus),]

write.csv(genus_2b_top30_2,"根据inade-seq提取的top30丰度表.csv")



genus_2b_top30_2<- read.csv("根据inade-seq提取的top30丰度表.csv")%>% dplyr::select(-X) %>% drop_na()
View(genus_2b_top30_2)

library(reshape2)


genus_2b_top30_long<- melt(genus_2b_top30_2,id="Genus")
colnames(genus_2b_top30_long) <- c("Genus","samples","abundance")

genus_2b_top30_long$group1<-group[match(genus_2b_top30_long$samples, rownames(group)),3]
genus_2b_top30_long$treatments<-group[match(genus_2b_top30_long$samples, rownames(group)),4]
genus_2b_top30_long$MPR_response1<-group[match(genus_2b_top30_long$samples, rownames(group)),5]
genus_2b_top30_long$MPR_response2<-group[match(genus_2b_top30_long$samples, rownames(group)),6]
genus_2b_top30_long$CR_response1<-group[match(genus_2b_top30_long$samples, rownames(group)),7]
genus_2b_top30_long$CR_response2<-group[match(genus_2b_top30_long$samples, rownames(group)),8]
View(genus_2b_top30_long)

# 
# colnames(list1)
# 
# list1<- list(unique(genus_2b_top30_2$Genus)) %>% as.data.frame()
# colnames(list1) <-"genus_id"
# list2<- list(genus_top30)
# setdiff(list2,list1)


genus <- as.data.frame(unique(genus_2b_top30_long$Genus))

colnames(genus) <- "genus_names"
rownames(genus) <- genus$genus_names

color=NULL
color[rownames(genus)] =distinctColorPalette(26)
color[c("Pre_NMPR","Pre_MPR","Post_MPR")]<-c("blue","red","green")

genus_2b_top30_long<-genus_2b_top30_long[-which(genus_2b_top30_long$group1 == "ME"),]
#genus_2b_top30_2$Genus <- factor(genus_2b_top30_2$Genus,levels = rev(genus_top30$genus_id))
genus_2b_top30_long_filter<-genus_2b_top30_long %>% dplyr::select("Genus","MPR_response2","abundance")
pattern <- c("Massilia",
             "Herbaspirillum",
             "Sphingomonas",
             "Lawsonella",
             "Rhodococcus",
             "Methylobacterium")
genus_2b_top30_long_filter_top<- genus_2b_top30_long_filter %>% filter(!grepl(pattern ,Genus))



write.csv(genus_2b_top30_long_filter_top,"genus_2b_top30_long_filter_top5.csv") 
View(genus_2b_top30_long_filter_top)
colnames(genus_2b_top30_long_filter)
View(genus_2b_top30_long_filter)
circos.clear()

pdf(file="genus_2b_invade_top.pdf", width=14, height=8, pointsize=10)

chordDiagram(genus_2b_top30_long_filter_top,
             grid.col =color,#颜色
             annotationTrack = "grid",
             transparency = 0.2,#透明度
             link.lwd = 0.00001,#线条宽度
             link.lty = 1,    # 线路类型
             link.border = 0,#边框颜色
             directional = -1,#表示线条的方向，0代表没有方向，1代表正向，-1代表反向，2代表双向
             diffHeight = mm_h(3),#外圈和中间连线的间隔
             direction.type = c("diffHeight","arrows"), #线条是否带有箭头
             link.arr.type = "big.arrow",#箭头类型
             annotationTrackHeight = c(0.04, 0.1))#网格高度


circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
              facing = "clockwise", niceFacing = T, adj = c(-0.5, 0.5), cex = 0.8)
  circos.axis(h = "top", labels.cex = 0.4,labels.niceFacing = F, labels.pos.adjust =F)
}, bg.border = NA)


legend("right",pch=20,legend=rownames(genus),
       col=color[rownames(genus)],bty="n",
       cex=1,pt.cex=3,border="black")

dev.off()

####alpha_diversity----

library(microeco)
library(ggplot2)
library(dplyr)

abundance<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/16S_结题报告/MPR_group/feature_tax_even.xls",sep="\t",header = T,check.names = FALSE)


###处理一下species这一列
abundance$Kingdom <- gsub("k__", "", abundance$Kingdom)
abundance$Phylum<- gsub("p__", "", abundance$Phylum)
abundance$Class<- gsub("c__", "", abundance$Class)
abundance$Order<- gsub("o__", "", abundance$Order)
abundance$Family<- gsub("f__", "", abundance$Family)
abundance$Genus<- gsub("g__", "", abundance$Genus)
abundance$Species <- gsub("_g_.*", "", abundance$Species)
abundance$Species <- gsub("s__", "", abundance$Species)

names(abundance)

genus0 <- abundance[,c(1:39)]
names(genus0)

rownames(genus0)<- genus0$ASV_ID 
genus1<- genus0[,-1] %>% as.data.frame
genus1[is.na(genus1)] <- 0
genus1 <- genus1[rowSums(genus1) != 0, ]
genus1 <- genus1[, colSums(genus1) != 0]

genus1[] <- lapply(genus1, function(x) as.numeric(as.character(x)))


library(vegan)
library(picante)

otu <- t(genus1)


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

write.csv(alpha,"D:/wangys_uestc_data2/OSCC_data/oscc_process/16S_结题报告/MPR_group/alpha_diversity/alpha_16s.csv")

sample_data<-read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/16S_结题报告/MPR_group/group_mpr.csv")

sample_data$orig.ident

library(reshape2)

matched_indices <- match(rownames(alpha), sample_data$samples)

sample_data_aligned <- sample_data[matched_indices, ]
sample_data_aligned <- na.omit(sample_data_aligned)

all(rownames(alpha) == sample_data_aligned$samples) # 结果应为TRUE

sample_data_aligned <- sample_data_aligned[sample_data_aligned$samples %in% rownames(alpha), ]

# 筛选出 alpha 中行名与 sample_data_aligned 的 samples 列匹配的行
alpha <- alpha[rownames(alpha) %in% sample_data_aligned$samples, ]

# 重新排序，使 alpha 和 sample_data_aligned 顺序一致
sample_data_aligned <- sample_data_aligned[match(rownames(alpha), sample_data_aligned$samples), ]

# 检查对齐结果
all(rownames(alpha) == sample_data_aligned$samples) # 结果应为 TRUE
View(alpha)

library(reshape2)

alpha$MPR_response2<-sample_data_aligned$MPR_response2
dat <- melt(alpha,id.vars = c(5),variable.name = "Alpha")


library(rstatix)
library(ggpubr)
dat$value<- as.numeric(dat$value)
df_p_val1 <- dat %>%group_by(Alpha)%>%
  wilcox_test(value~ MPR_response2) %>%
  adjust_pvalue(p.col = "p", method = "bonferroni") %>%
  add_significance(p.col = "p.adj") %>%
  add_xy_position(x = "MMPR_response2", dodge = 0.8)

response_color <- c("ME"="black",
                    "Pre_NMPR"="#4974a4",
                    "Pre_MPR"="#4dae47",
                    "Post_MPR"="#f29600")

levels(dat$MPR_response2)
dat$MPR_response2 <- factor(dat$MPR_response2,levels = c("ME","Pre_NMPR","Pre_MPR","Post_MPR"))

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

df1 <- PMCMR_compare1(dat,'Alpha','MPR_response2','value')

names(dat)


alpha_plot<- ggplot(dat, aes(x = MPR_response2, y = value, fill = MPR_response2)) + 
  stat_boxplot(geom = "errorbar", linewidth=1)+
  geom_boxplot(linewidth=1,outlier.shape = NA)+
  scale_fill_manual(values = response_color) +
  geom_text(data=df1,aes(x=MPR_response2,y=mean+1.5*std,label=Letters))+
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
alpha_plot
ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/16S_结题报告/MPR_group/alpha_diversity/16s_seq_alpha_diversity.pdf",alpha_plot,width=12,height=4)

####betadiveristy----

library(vegan)
library(ggplot2)

abundance<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/16S_结题报告/MPR_group/feature_tax_even2.xls",sep="\t",header = T,row.names = 1,check.names = FALSE)

###处理一下species这一列
abundance$Kingdom <- gsub("k__", "", abundance$Kingdom)
abundance$Phylum<- gsub("p__", "", abundance$Phylum)
abundance$Class<- gsub("c__", "", abundance$Class)
abundance$Order<- gsub("o__", "", abundance$Order)
abundance$Family<- gsub("f__", "", abundance$Family)
abundance$Genus<- gsub("g__", "", abundance$Genus)
abundance$Species <- gsub("_g_.*", "", abundance$Species)
abundance$Species <- gsub("s__", "", abundance$Species)

names(abundance)

abundance <- abundance[,c(1:24)]
ncol(abundance)
head(abundance)
View(abundance)

group <- read.table("D:/wangys_uestc_data2/OSCC_data/oscc_process/16S_结题报告/MPR_group/group_mpr.csv", header=T, sep=",",row.names = 1, comment.char="",stringsAsFactors = TRUE,quote = "")
nrow(group)

abundance1<- abundance[,match(rownames(group), colnames(abundance))]



View(abundance1)

library(vegan)
library(ape)
abundnce_dist<- vegdist(t(abundance1), method="bray")
pcoa <- cmdscale(abundnce_dist, k=3, eig=T)#高维度数据降低到低维度数据的操作
pcoa_points <- as.data.frame(pcoa$points)#提取三个维度的PCoA值
sum_eig <- sum(pcoa$eig)#提取特征值作为
eig_percent <- round(pcoa$eig/sum_eig*100,1)
colnames(pcoa_points) <- paste0("PCoA", 1:2)
pcoa_result <- cbind(pcoa_points,t(abundance1))#PCoA值和OTU合并
View(pcoa_result)
abundance1_t <- as.data.frame(t(abundance1))
View(abundance1_t)

nrow(abundance1_t)
nrow(group)

MPR_response2 <- group$MPR_response2#提取Moisture分组信息

div <- adonis2(abundance1_t ~ MPR_response2, data = abundance1_t, permutations = 999, method="bray")#adonis
adonis <- paste0("adonis R2: ",round(div$R2,2), "; P-value: ", div$`Pr(>F)`)#提取p值，组成一个字符串

group$MPR_response2<-factor(group$MPR_response2,levels = c("ME",'Pre_NMPR','Pre_MPR','Post_MPR'))

# matched_cols <- match(rownames(group), rownames(pcoa_result))
# matched_cols <- matched_cols[!is.na(matched_cols)]
# pcoa_result$group <- group[matched_cols,][,6]

#pcoa_result$group <- group[match(rownames(group),rownames(pcoa_result)),6]

p<-ggplot(pcoa_result, aes(x=PCoA1, y=PCoA2, color=group$MPR_response2)) +
  scale_color_manual(values =c("black","#4974a4","#4dae47","#f29600"))+
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

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/16S_结题报告/MPR_group/alpha_diversity/16s_MPR_beta_diversity2.pdf",p,height = 7,width = 9)


####adonis分析####
pairwise.adonis <-function(x,factors, sim.method, p.adjust.m)
{library(vegan)
  co = as.matrix(combn(unique(factors),2))
  pairs = c()
  F.Model =c()
  R2 = c()
  p.value = c()
  for(elem in 1:ncol(co)){
    ad = adonis(x[factors %in%c(as.character(co[1,elem]),as.character(co[2,elem])),] ~
                  factors[factors %in%c(as.character(co[1,elem]),as.character(co[2,elem]))] , method =sim.method);
    pairs =c(pairs,paste(co[1,elem],'vs',co[2,elem]));
    F.Model =c(F.Model,ad$aov.tab[1,4]);
    R2 = c(R2,ad$aov.tab[1,5]);
    p.value = c(p.value,ad$aov.tab[1,6])}
  p.adjusted =p.adjust(p.value,method=p.adjust.m)
  pairw.res = data.frame(pairs,F.Model,R2,p.value,p.adjusted)
  return(pairw.res)
}

pairwise.adonis(abundance1_t, group$MPR_response2, sim.method="bray", p.adjust.m= "BH")




