###预后分析数据收集

####https://github.com/twbattaglia/MicrobeDS
####This package contains datasets provided by large-scale microbiome studies. Each dataset is formatted for use with phyloseq.
####(https://joey711.github.io/phyloseq/).

remotes::install_github("twbattaglia/MicrobeDS")

TwinsUK
?MicrobeDS
# Load Library
library(MicrobeDS)
library(sample_data)
# Load selected datasets
data('HMPv35')

# Check number of samples
nsamples(HMPv35)

# Check sample metadata
sample_data(HMPv35)
View(HMPv35)
sample_table <- HMPv35@sam_data %>% as.data.frame()


####qa10394----

data(qa10394)
View(qa10394)
sample_table <- qa10394@sam_data %>% as.data.frame()


####TwinsUK----

data(TwinsUK)
View(TwinsUK)
sample_table <- TwinsUK@sam_data %>% as.data.frame()

####RISK_CCFA----

data(RISK_CCFA)
View(RISK_CCFA)
sample_table <- RISK_CCFA@sam_data %>% as.data.frame()

####MovingPictures----

data(MovingPictures)
View(MovingPictures)
sample_table <- MovingPictures@sam_data %>% as.data.frame()


####Hartwig----

#可参考 ICBatlas 或 Cancer Immunology Research 的相关研究，

####Hartwig: Description: Hartwig Medical Foundation (HMF)

library(phyloseq)

data(Hartwig)

# Check number of samples
nsamples(Hartwig)

# Check sample metadata
ntaxa(Hartwig)
View(Hartwig)

####Hartwig, Capnocytophaga,Fusobacterium and Streptococcus----

otu_table <- as.data.frame(phyloseq::otu_table(Hartwig))
otu_table_data <- otu_table(Hartwig)
tax_table_data <- phyloseq::tax_table(Hartwig)

####Hartwig Fusobacterium hight low 分组密度图绘制----

# 提取 Fusobacterium 的丰度

fusobacterium_OTUs <- rownames(tax_table_data)[tax_table_data[, "Genus"] == "Fusobacterium"]

unique(tax_table_data[, "Genus"])

fusobacterium_abundance <- otu_table_data[fusobacterium_OTUs,drop = FALSE]

head(fusobacterium_abundance)
View(fusobacterium_abundance)
length(fusobacterium_abundance)
# 计算每个样本中Fusobacterium的百分比
total_abundance <- colSums(otu_table_data)  # 每个样本的总丰度
length(total_abundance)

fusobacterium_percentage <- fusobacterium_abundance / total_abundance * 100

fusobacterium_percentage<- as.data.frame(fusobacterium_percentage)

library(ggplot2)
library(dplyr)
set.seed(123)

# 计算 75% 分位点
threshold <- quantile(data$Percentage, 0.75,na.rm = TRUE)

# 添加分组信息
data$fusobacterium_group <- ifelse(data$Percentage >= threshold, "High", "Low")

unique(data$fusobacterium_group)

pp1<- ggplot(data, aes(x = Percentage, fill = fusobacterium_group)) +
  geom_density(alpha = 1) +  # 密度图，透明度为0.6
  geom_vline(xintercept = threshold, linetype = "dashed") +  # 添加分位点虚线
  scale_fill_manual(values = c("Low" = "skyblue", "High" = "red","NaN"="skyblue")) +  # 自定义填充颜色
  labs(
    x = "Fusobacterium (%)",
    y = "Density",
    fill = "Fusobacterium"
  ) +
  theme_minimal() +  # 设置主题
  theme(
    legend.position = "top",  # 调整图例位置
    axis.ticks.x = element_line(linewidth = 0.8), # 加粗 X 轴刻度线
    axis.ticks.y = element_line(linewidth = 0.8), # 加粗 Y 轴刻度线
    axis.title = element_text(size = 14),   # 调整坐标轴标题的大小
    axis.text = element_text(size = 12),     # 调整坐标轴文字的大小
    axis.line.x = element_line(linewidth = 0.8),  # 加粗 X 轴线条
    axis.line.y = element_line(linewidth = 0.8)  # 加粗 Y 轴线条
  )+
  # scale_x_continuous(expand = expansion(mult = c(0.1, 0.1))) + 
  # scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) + 
  coord_cartesian(xlim = c(-1, 80), ylim = c(0, 0.3)) 
pp1

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Apan-cancer_analysis_of_the_microbiome/Hartwig_Fusobacterium_density.pdf",pp1,width = 5,height = 5)



taxa_table <- as(sample_data(taxa_table), "data.frame")

names(taxa_table)
View(Hartwig)

sample_table <- Hartwig@sam_data %>% as.data.frame()


dim(sample_table)
View(sample_table)
df <- sample_table
df$survival_time <- ifelse(
  df$vital_status == "Dead", 
  df$days_to_death, 
  df$days_to_last_followup
)

df$disease_code
df$event_status <- ifelse(df$vital_status == "Dead", 1, 0)
names(df)
head(df)
df <- as(sample_data(df), "data.frame")

df <- df %>%
  tibble::rownames_to_column(var = "SampleID")

View(df)

survival_data <- df %>%
  dplyr::select(SampleID, survival_time, event_status,gender,
                age_at_initial_pathologic_diagnosis,
                height,weight,
                tobacco_smoking_history,venous_invasion,lymphatic_invasion,residual_tumor,
                tumor_tissue_site,
                history_of_neoadjuvant_treatment,histological_type)
View(survival_data)
unique(survival_data$tobacco_smoking_history)


survival_data_numeric <- survival_data %>%
  mutate(
    gender=ifelse(gender == "MALE", 1, 0),
    tobacco_smoking_history <- case_when(
      tobacco_smoking_history == "Lifelong Non-smoker" ~ 0,
      tobacco_smoking_history == "Current reformed smoker for <= 15 years" ~ 1,
      tobacco_smoking_history == "Current reformed smoker for > 15 years" ~ 2,
      tobacco_smoking_history == "Current smoker" ~ 3
    ),
    venous_invasion = ifelse(venous_invasion == "True", 1, 0),
    lymphatic_invasion = ifelse(lymphatic_invasion == "True", 1, 0),
    residual_tumor = case_when(
      residual_tumor == "R0" ~ 0,
      residual_tumor == "R1" ~ 1,
      residual_tumor == "R2" ~ 2,
      residual_tumor == "RX" ~ 3
    )
  )
survival_data_numeric <- survival_data_numeric[, -ncol(survival_data_numeric)]
View(survival_data_numeric )
names(data)  
data$Sample
names(survival_data_numeric)
names(df)
unique(df$histological_type)

merged_data <- merge(data, df, 
                     by.x = "Sample", by.y = "SampleID", 
                     all = FALSE)

names(merged_data)

merged_data2<- merged_data %>% filter(tumor_tissue_site %in% "Head and Neck")
View(merged_data2)
names(merged_data2)
names(merged_data2)

survival_data_numeric
library(survival)


merged_data2$fusobacterium_group <- na_if(merged_data2$fusobacterium_group, "NaN")
merged_data2$fusobacterium_group <- factor(
  merged_data2$fusobacterium_group,
  levels = c("Low", "High", "NaN") # 包含所有可能值
)
merged_data2$fusobacterium_group <- relevel(merged_data2$fusobacterium_group, ref = "Low")


merged_data$fusobacterium_cut <- na_if(merged_data$fusobacterium_cut, "NaN")
merged_data$fusobacterium_cut <- factor(
  merged_data$fusobacterium_cut,
  levels = c("Low", "High", "NaN") # 包含所有可能值
)
merged_data$fusobacterium_cut <- relevel(merged_data$fusobacterium_cut, ref = "Low")

cox_model <- coxph(Surv(survival_time, event_status) ~ fusobacterium_cut, data = merged_data)

summary(cox_model)
cox.zph(cox_model)

library(survival)
library(survminer)
sur.cut <- surv_cutpoint(merged_data, time = "survival_time", 
                         event = "event_status",
                         variables = "Percentage")
summary(sur.cut)
#对截断值两侧的数据分布可视化。
plot(sur.cut, "Percentage", palette = "npg")

#使用surv_categorize函数，根据截断值转为二分类变量。
sur.cat <- surv_categorize(sur.cut, labels = c(0, 1))


data$Percentage <- as.factor(data$Percentage)
#0为较低得分，1为较高得分，指定较高得分为参考水平。参考水平需要加引号。
data$Percentage <- relevel(data$Percentage, ref = "1")


merged_data$mayoscore5_cut <- as.factor(merged_data$mayoscore5_cut)
#0为较低得分，1为较高得分，指定较高得分为参考水平。参考水平需要加引号。
data$mayoscore5_cut <- relevel(data$mayoscore5_cut, ref = "1")


merged_data$fusobacterium_cut <-
  ifelse(merged_data$Percentage >= 0.5347594, "High", "Low")

cox_model_1 <- coxph(Surv(survival_time, event_status) ~ fusobacterium_cut, data = merged_data)
summary(cox_model_1)

merged_data$fusobacterium_cut <- relevel(merged_data$fusobacterium_cut, ref = "0")

cox_model_2 <- coxph(Surv(time, censor) ~  fusobacterium_cut, data = data)
summary(cox_model_2)



####TCGA----
####Description: The Cancer Genome Atlas (TCGA)
##### Taxonomy contains only Kingdom & Genus level 
# due to missing intermediate taxonomy annotations
library(MicrobeDS)
library(sample_data)
library(ggplot2)
library(dplyr)
data(TCGA)
View(TCGA)
sample_table <- TCGA@sam_data %>% as.data.frame()
otu_table <- as.data.frame(phyloseq::otu_table(TCGA))

otu_table_data <- phyloseq::otu_table(TCGA)
tax_table_data <- phyloseq::tax_table(TCGA)

####TCGA Fusobacterium hight low 分组密度图绘制----

# 提取 Fusobacterium 的丰度

fusobacterium_OTUs <- rownames(tax_table_data)[tax_table_data[, "Genus"] == "Fusobacterium"]

unique(tax_table_data[, "Genus"])

fusobacterium_abundance <- otu_table_data[fusobacterium_OTUs,drop = FALSE]
head(fusobacterium_abundance)
View(fusobacterium_abundance)
length(fusobacterium_abundance)
# 计算每个样本中Fusobacterium的百分比
total_abundance <- colSums(otu_table_data)  # 每个样本的总丰度
length(total_abundance)

fusobacterium_percentage <- fusobacterium_abundance / total_abundance * 100
length(fusobacterium_percentage )
fusobacterium_percentage<- as.data.frame(fusobacterium_percentage)
View(fusobacterium_percentage)
# 提取 Capnocytophaga 的丰度

Capnocytophaga_OTUs <- rownames(tax_table_data)[tax_table_data[, "Genus"] == "Capnocytophaga"]

unique(tax_table_data[, "Genus"])

Capnocytophaga_abundance <- otu_table_data[Capnocytophaga_OTUs,drop = FALSE]
head(Capnocytophaga_abundance)
View(Capnocytophaga_abundance)
length(Capnocytophaga_abundance)
# 计算每个样本中Fusobacterium的百分比
total_abundance <- colSums(otu_table_data)  # 每个样本的总丰度
length(total_abundance)

Capnocytophaga_percentage <- Capnocytophaga_abundance / total_abundance * 100

Capnocytophaga_percentage<- as.data.frame(Capnocytophaga_percentage)


# 提取 Streptococcus 的丰度

Streptococcus_OTUs <- rownames(tax_table_data)[tax_table_data[, "Genus"] == "Streptococcus"]

unique(tax_table_data[, "Genus"])

Streptococcus_abundance <- otu_table_data[Streptococcus_OTUs,drop = FALSE]
head(Streptococcus_abundance)
View(Streptococcus_abundance)
length(Streptococcus_abundance)
# 计算每个样本中Streptococcus的百分比
total_abundance <- colSums(otu_table_data)  # 每个样本的总丰度
length(total_abundance)

Streptococcus_percentage <- Streptococcus_abundance / total_abundance * 100

Streptococcus_percentage<- as.data.frame(Streptococcus_percentage)


total_abundance<- as.data.frame(total_abundance)
data <- data.frame(
  Sample = rownames(total_abundance),
  Fusobacterium_percentage = unlist(fusobacterium_percentage),
  Capnocytophaga_percentage=unlist(Capnocytophaga_percentage),
  Streptococcus_percentage =unlist(Streptococcus_percentage)
  
)



library(ggplot2)
library(dplyr)
set.seed(123)

# 计算 75% 分位点
threshold1 <- quantile(data$Fusobacterium_percentage, 0.75,na.rm = TRUE)

# 添加分组信息
data$fusobacterium_group <- ifelse(data$Fusobacterium_percentage >= threshold, "High", "Low")

unique(data$fusobacterium_group)

View(data)
pp1<- ggplot(data, aes(x = Fusobacterium_percentage, fill = fusobacterium_group)) +
  geom_density(alpha = 1) +  # 密度图，透明度为0.6
  geom_vline(xintercept = threshold, linetype = "dashed") +  # 添加分位点虚线
  scale_fill_manual(values = c("Low" = "skyblue", "High" = "red","NaN"="skyblue")) +  # 自定义填充颜色
  labs(
    x = "Fusobacterium (%)",
    y = "Density",
    fill = "Fusobacterium"
  ) +
  theme_minimal() +  # 设置主题
  theme(
    legend.position = "top",  # 调整图例位置
    axis.ticks.x = element_line(linewidth = 0.8), # 加粗 X 轴刻度线
    axis.ticks.y = element_line(linewidth = 0.8), # 加粗 Y 轴刻度线
    axis.title = element_text(size = 14),   # 调整坐标轴标题的大小
    axis.text = element_text(size = 12),     # 调整坐标轴文字的大小
    axis.line.x = element_line(linewidth = 0.8),  # 加粗 X 轴线条
    axis.line.y = element_line(linewidth = 0.8)  # 加粗 Y 轴线条
  )+
  # scale_x_continuous(expand = expansion(mult = c(0.1, 0.1))) + 
  # scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) + 
  coord_cartesian(xlim = c(-0.1,3), ylim = c(0, 0.75))
pp1

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Apan-cancer_analysis_of_the_microbiome/TCGA_Fusobacterium_density.pdf",pp1,width = 5,height = 5)


library(ggplot2)
library(dplyr)
set.seed(123)

# 计算 75% 分位点
threshold <- quantile(data$Capnocytophaga_percentage, 0.75,na.rm = TRUE)

# 添加分组信息
data$Capnocytophaga_group <- ifelse(data$Capnocytophaga_percentage >= threshold, "High", "Low")

unique(data$Capnocytophaga_group)

View(data)
pp1<- ggplot(data, aes(x = Capnocytophaga_percentage, fill = Capnocytophaga_group)) +
  geom_density(alpha = 1) +  # 密度图，透明度为0.6
  geom_vline(xintercept = threshold, linetype = "dashed") +  # 添加分位点虚线
  scale_fill_manual(values = c("Low" = "skyblue", "High" = "red")) +  # 自定义填充颜色
  labs(
    x = "Capnocytophaga (%)",
    y = "Density",
    fill = "Capnocytophaga"
  ) +
  theme_minimal() +  # 设置主题
  theme(
    legend.position = "top",  # 调整图例位置
    axis.ticks.x = element_line(linewidth = 0.8), # 加粗 X 轴刻度线
    axis.ticks.y = element_line(linewidth = 0.8), # 加粗 Y 轴刻度线
    axis.title = element_text(size = 14),   # 调整坐标轴标题的大小
    axis.text = element_text(size = 12),     # 调整坐标轴文字的大小
    axis.line.x = element_line(linewidth = 0.8),  # 加粗 X 轴线条
    axis.line.y = element_line(linewidth = 0.8)  # 加粗 Y 轴线条
  )+
  # scale_x_continuous(expand = expansion(mult = c(0.1, 0.1))) + 
  # scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) + 
  coord_cartesian(xlim = c(-0.1,3), ylim = c(0, 0.75))
pp1

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Apan-cancer_analysis_of_the_microbiome/TCGA_Capnocytophaga_density.pdf",pp1,width = 5,height = 5)

library(ggplot2)
library(dplyr)
set.seed(123)

# 计算 75% 分位点
threshold <- quantile(data$Streptococcus_percentage, 0.75,na.rm = TRUE)

# 添加分组信息
data$Streptococcus_group <- ifelse(data$Streptococcus_percentage >= threshold, "High", "Low")

unique(data$Streptococcus_group)

View(data)
pp1<- ggplot(data, aes(x = Streptococcus_percentage, fill = Streptococcus_group)) +
  geom_density(alpha = 1) +  # 密度图，透明度为0.6
  geom_vline(xintercept = threshold, linetype = "dashed") +  # 添加分位点虚线
  scale_fill_manual(values = c("Low" = "skyblue", "High" = "red")) +  # 自定义填充颜色
  labs(
    x = "Streptococcus (%)",
    y = "Density",
    fill = "Streptococcus"
  ) +
  theme_minimal() +  # 设置主题
  theme(
    legend.position = "top",  # 调整图例位置
    axis.ticks.x = element_line(linewidth = 0.8), # 加粗 X 轴刻度线
    axis.ticks.y = element_line(linewidth = 0.8), # 加粗 Y 轴刻度线
    axis.title = element_text(size = 14),   # 调整坐标轴标题的大小
    axis.text = element_text(size = 12),     # 调整坐标轴文字的大小
    axis.line.x = element_line(linewidth = 0.8),  # 加粗 X 轴线条
    axis.line.y = element_line(linewidth = 0.8)  # 加粗 Y 轴线条
  )+
  # scale_x_continuous(expand = expansion(mult = c(0.1, 0.1))) + 
  # scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) + 
  coord_cartesian(xlim = c(-0.1,60), ylim = c(0, 0.75))
pp1

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Apan-cancer_analysis_of_the_microbiome/TCGA_Streptococcus_density.pdf",pp1,width = 5,height = 5)



####TCGA预后分析----

taxa_table <-TCGA@tax_table %>% as.data.frame()
taxa_table <- as(sample_data(taxa_table), "data.frame")
names(taxa_table)

sample_table <- TCGA@sam_data %>% as.data.frame()
sample_table$
sample_table <- as(phyloseq::sample_data(sample_table), "data.frame")

write.csv(sample_table,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Apan-cancer_analysis_of_the_microbiome/TCGA_clin_data.csv")

sample_table<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Apan-cancer_analysis_of_the_microbiome/TCGA_clin_data.csv")
names(sample_table)


library(dplyr)
library(stringr)
library(TCGAbiolinks)
?GDCquery_clinic
clin_data <- GDCquery_clinic(project = "TCGA-COAD", type = "clinical")
print(dim(clin_data))

merged_df <- data.frame()

df1 <- sample_table %>%
  mutate(submitter_id = str_extract(filename, "TCGA-\\w+-\\w+"))

clin_data <- GDCquery_clinic(project = "TCGA-HNSC", type = "clinical")
clin_data$days_to_last_follow_up

library(TCGAbiolinks)
library(dplyr)

# 定义项目列表
project_ids <- c("TCGA-ACC",
  "TCGA-BLCA", "TCGA-BRCA", "TCGA-CESC", "TCGA-CHOL", "TCGA-COAD", "TCGA-DLBC", 
  "TCGA-ESCA", "TCGA-GBM", "TCGA-HNSC", "TCGA-KICH", "TCGA-KIRC", "TCGA-KIRP", 
  "TCGA-LAML", "TCGA-LGG", "TCGA-LIHC", "TCGA-LUAD", "TCGA-LUSC", "TCGA-MESO", 
  "TCGA-OV", "TCGA-PAAD", "TCGA-PCPG", "TCGA-PRAD", "TCGA-READ", "TCGA-SARC", 
  "TCGA-SKCM", "TCGA-STAD", "TCGA-TGCT", "TCGA-THCA", "TCGA-THYM", "TCGA-UCEC", 
  "TCGA-UCS", "TCGA-UVM"
)

# 初始化合并数据框
final_clin_data <- data.frame()


for (project_id in project_ids) {
  cat("Processing project:", project_id, "\n")
  
  # 捕获下载数据时的错误
  clin_data <- tryCatch(
    {
      GDCquery_clinic(project = project_id, type = "clinical")
    },
    error = function(e) {
      message("Error with project:", project_id, "-", e$message)
      return(NULL)
    }
  )
  
  # 如果 clin_data 为 NULL 或为空，则跳过当前项目
  if (is.null(clin_data) || nrow(clin_data) == 0) {
    cat("No data or error for project:", project_id, "\n")
    next
  }
  
  # 合并数据
  tryCatch({
    if (nrow(final_clin_data) == 0) {
      # 如果 final_clin_data 是空的，直接赋值
      final_clin_data <- clin_data
    } else {
      # 找到公共列
      common_cols <- intersect(colnames(final_clin_data), colnames(clin_data))
      
      if (length(common_cols) == 0) {
        # 如果没有公共列，跳过该项目
        cat("No common columns for project:", project_id, "\n")
        next
      }
      
      # 按公共列合并
      final_clin_data <- bind_rows(
        final_clin_data[, common_cols, drop = FALSE],
        clin_data[, common_cols, drop = FALSE]
      )
    }
    
    cat("Successfully processed project:", project_id, "\n")
    
  }, error = function(e) {
    # 捕获合并错误，继续处理下一个项目
    cat("Error during merge for project:", project_id, "-", e$message, "\n")
  })
}

# 查看最终合并数据的维度
dim(final_clin_data)

# 检查合并数据的前几行
head(final_clin_data)

write.csv(final_clin_data,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Apan-cancer_analysis_of_the_microbiome/TCGA_clin_data_合并.csv")
final_clin_data$days_to_last_follow_up

merged_df <- df1 %>%
  left_join(final_clin_data, by = "submitter_id") %>%
  select(-ends_with(".y"))  # 移除重复列
merged_df$days_to_last_follow_up
View(merged_df)

df <- merged_df

df$survival_time <- ifelse(
  df$vital_status_label == "Dead", 
  df$days_to_death, 
  df$days_to_last_follow_up
)

df$event_status <- ifelse(df$vital_status_label == "Dead", 1, 0)

names(df)[1] <- "sampleid"
View(df)
View(data)
View(df)
data$Sampleid<- rownames(data)

merged_data <- merge(data, df, 
                     by.x = "Sampleid", by.y = "sample", 
                     all = FALSE)

names(merged_data)

write.csv(merged_data,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Apan-cancer_analysis_of_the_microbiome/TCGA_merged_data_预后分析数据.csv")

merged_data<-read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Apan-cancer_analysis_of_the_microbiome/TCGA_merged_data_预后分析数据.csv")
  
# library(survminer)
# sur.cut <- surv_cutpoint(merged_data, time = "survival_time", 
#                          event = "event_status",
#                          variables = "Percentage")
# summary(sur.cut)
# #对截断值两侧的数据分布可视化。
# plot(sur.cut, "Percentage", palette = "npg")
# 
# #使用surv_categorize函数，根据截断值转为二分类变量。
# sur.cat <- surv_categorize(sur.cut, labels = c("High", "Low"))
# data

# merged_data$fusobacterium_cut <-
#   ifelse(merged_data$Percentage >= 0.001779528, "High", "Low")

library(survival)

merged_data$fusobacterium_group <- factor(
  merged_data$fusobacterium_group,
  levels = c("Low", "High") # 包含所有可能值
)
merged_data$fusobacterium_group <- relevel(merged_data$fusobacterium_group, ref = "Low")

merged_data$fusobacterium_cut <- factor(
  merged_data$fusobacterium_cut,
  levels = c("Low", "High") # 包含所有可能值
)
merged_data$fusobacterium_cut <- relevel(merged_data$fusobacterium_cut, ref = "Low")

names(merged_data)

cox_model <- coxph(Surv(survival_time, event_status) ~ Percentage, data = merged_data)
summary(cox_model)
cox.zph(cox_model)

View(merged_data)
merged_data$age_at_index
merged_data$project


merged_data$fusobacterium_group <- factor(
  merged_data$fusobacterium_group,
  levels = c("Low", "Medium","High") # 包含所有可能值
)
merged_data$fusobacterium_group <- relevel(merged_data$fusobacterium_group, ref = "Low")

cox_model <- coxph(Surv(survival_time, event_status) ~ fusobacterium_group+age_at_index+
                     alcohol_history+
                     portion_weight+gender.x, data = merged_data)



summary(cox_model)
cox.zph(cox_model)

#####fusobacterium HNSC 头颈鳞癌----
merged_data2<- merged_data %>% filter(project %in% "TCGA-HNSC")

# 计算 75% 分位点
threshold <- quantile(merged_data2$Fusobacterium_percentage, 0.75,na.rm = TRUE)

merged_data2$Fusobacterium_group <- ifelse(merged_data2$Fusobacterium_percentage <= threshold, "Low", "High")


pp1<- ggplot(merged_data2, aes(x = Fusobacterium_percentage, fill = Fusobacterium_group)) +
  geom_density(alpha = 1) +  # 密度图，透明度为0.6
  geom_vline(xintercept = threshold, linetype = "dashed") +  # 添加分位点虚线
  scale_fill_manual(values = c("Low" = "skyblue", "High" = "red")) +  # 自定义填充颜色
  labs(
    x = "Fusobacterium (%)",
    y = "Density",
    fill = "Fusobacterium"
  ) +
  theme_minimal() +  # 设置主题
  theme(
    legend.position = "top",  # 调整图例位置
    axis.ticks.x = element_line(linewidth = 0.8), # 加粗 X 轴刻度线
    axis.ticks.y = element_line(linewidth = 0.8), # 加粗 Y 轴刻度线
    axis.title = element_text(size = 14),   # 调整坐标轴标题的大小
    axis.text = element_text(size = 12),     # 调整坐标轴文字的大小
    axis.line.x = element_line(linewidth = 0.8),  # 加粗 X 轴线条
    axis.line.y = element_line(linewidth = 0.8)  # 加粗 Y 轴线条
  )+
  # scale_x_continuous(expand = expansion(mult = c(0.1, 0.1))) + 
  # scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) + 
  coord_cartesian(xlim = c(-0.1,10), ylim = c(0, 0.75))
pp1

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Apan-cancer_analysis_of_the_microbiome/TCGA-HNSC头颈/TCGA_Fusobacterium_density.pdf",pp1,width = 5,height = 5)


merged_data2$Fusobacterium_group <- factor(
  merged_data2$Fusobacterium_group,
  levels = c("Low", "High") # 包含所有可能值
)
merged_data2$Fusobacterium_group <- relevel(merged_data2$Fusobacterium_group, ref = "Low")

merged_data2$gender.x <- factor(merged_data2$gender.x)
sum(is.na(merged_data2$gender.x))

names(merged_data2)
unique(merged_data2$years_smoked)
merged_data2$years_smoked_history<-ifelse(merged_data2$years_smoked <= 20, "Low", "High")
merged_data2$portion_weight2 <- log(merged_data2$portion_weight + 1)
merged_data2$years_smoked2 <- log(merged_data2$years_smoked + 1)
merged_data2$years_smoked2[is.na(merged_data2$years_smoked2)] <- mean(merged_data2$years_smoked2, na.rm = TRUE)



cox_model <- coxph(Surv(survival_time, event_status) ~ Fusobacterium_cut+age_at_index+
                     alcohol_history+prior_malignancy+
                     portion_weight2+gender.x, data = merged_data2)

summary(cox_model)
cox.zph(cox_model)

cox_model <- coxph(Surv(survival_time, event_status) ~ Fusobacterium_percentage, data = merged_data2)

merged_data2$Fusobacterium_cut <- factor(
  merged_data2$Fusobacterium_cut,
  levels = c("Low", "High") # 包含所有可能值
)
merged_data2$Fusobacterium_cut <- relevel(merged_data2$Fusobacterium_cut, ref = "Low")
cox_model <- coxph(Surv(survival_time, event_status) ~ Fusobacterium_cut, data = merged_data2)

summary(cox_model)
f1_Fusobacterium<-ggforest(cox_model,
                  data=merged_data2,
                  main = "Hazard ratio",        # 标题
                  cpositions = c(0.01, 0.1, 0.3), #左边三列的相对距离
                  fontsize = 1, # 字体大小
                  refLabel = "reference", #显示因子的参考水平
                  noDigits = 3)  #小数点位数

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Apan-cancer_analysis_of_the_microbiome/TCGA-HNSC头颈/fusobacterium_HNSC_ggforest_cut_多因素.pdf",f1_Fusobacterium,width =13,height=6)

library(survminer)
sur.cut <- surv_cutpoint(merged_data2, time = "survival_time",
                         event = "event_status",
                         variables = "Fusobacterium_percentage")
summary(sur.cut)
#对截断值两侧的数据分布可视化。
plot(sur.cut, "Fusobacterium_percentage", palette = "npg")

#使用surv_categorize函数，根据截断值转为二分类变量。
sur.cat <- surv_categorize(sur.cut, labels = c("High", "Low"))


merged_data2$Fusobacterium_cut <-
  ifelse(merged_data2$Fusobacterium_percentage >= 0.003430682, "High", "Low")


fit <- survfit(Surv(survival_time, event_status) ~Fusobacterium_cut,  # 创建生存对象 
               data = merged_data2)

pdf("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Apan-cancer_analysis_of_the_microbiome/TCGA-HNSC头颈/Fusobacterium_HNSC_curvep_cut.pdf", width = 6, height =7)
ggsurvplot(fit, # 创建的拟合对象
           data = merged_data2,  # 指定变量数据来源
           conf.int = TRUE, # 显示置信区间
           pval = TRUE, # 添加P值
           surv.median.line = "hv",  # 添加中位生存时间线
           risk.table = TRUE, # 添加风险表
           ylab="Overall Survival(%)",
           xlab = "Follow up time(d)", # 指定x轴标签
           legend = c(0.8,0.75), # 指定图例位置
           legend.title = "", # 设置图例标题，这里设置不显示标题，用空格替代
           palette = c("skyblue", "red"),
           legend.labs = c("Low", "High"), # 指定图例分组标签
           break.x.by = 1000,
           title = "Fusobacterium HNSC" )  # 设置x轴刻度间距
dev.off()


#####Capnocytophaga HNSC----
merged_data2<- merged_data %>% filter(project %in% "TCGA-HNSC")

# 计算 75% 分位点
threshold <- quantile(merged_data2$Capnocytophaga_percentage, 0.75,na.rm = TRUE)

merged_data2$Capnocytophaga_group <- ifelse(merged_data2$Capnocytophaga_percentage <= threshold, "Low", "High")

pp1<- ggplot(merged_data2, aes(x = Capnocytophaga_percentage, fill = Capnocytophaga_group)) +
  geom_density(alpha = 1) +  # 密度图，透明度为0.6
  geom_vline(xintercept = threshold, linetype = "dashed") +  # 添加分位点虚线
  scale_fill_manual(values = c("Low" = "skyblue", "High" = "red")) +  # 自定义填充颜色
  labs(
    x = "Capnocytophaga (%)",
    y = "Density",
    fill = "Capnocytophaga"
  ) +
  theme_minimal() +  # 设置主题
  theme(
    legend.position = "top",  # 调整图例位置
    axis.ticks.x = element_line(linewidth = 0.8), # 加粗 X 轴刻度线
    axis.ticks.y = element_line(linewidth = 0.8), # 加粗 Y 轴刻度线
    axis.title = element_text(size = 14),   # 调整坐标轴标题的大小
    axis.text = element_text(size = 12),     # 调整坐标轴文字的大小
    axis.line.x = element_line(linewidth = 0.8),  # 加粗 X 轴线条
    axis.line.y = element_line(linewidth = 0.8)  # 加粗 Y 轴线条
  )+
  # scale_x_continuous(expand = expansion(mult = c(0.1, 0.1))) + 
  # scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) + 
  coord_cartesian(xlim = c(-0.1,10), ylim = c(0, 0.75))
pp1

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Apan-cancer_analysis_of_the_microbiome/TCGA-HNSC头颈/TCGA_Capnocytophaga_density.pdf",pp1,width = 5,height = 5)


merged_data2$Capnocytophaga_group <- factor(
  merged_data2$Capnocytophaga_group,
  levels = c("Low", "High") # 包含所有可能值
)
merged_data2$Capnocytophaga_group <- relevel(merged_data2$Capnocytophaga_group, ref = "Low")

merged_data2$gender.x <- factor(merged_data2$gender.x)
sum(is.na(merged_data2$gender.x))

names(merged_data2)
unique(merged_data2$years_smoked)
merged_data2$years_smoked_history<-ifelse(merged_data2$years_smoked <= 20, "Low", "High")
merged_data2$portion_weight2 <- log(merged_data2$portion_weight + 1)
merged_data2$years_smoked2 <- log(merged_data2$years_smoked + 1)
merged_data2$years_smoked2[is.na(merged_data2$years_smoked2)] <- mean(merged_data2$years_smoked2, na.rm = TRUE)

cox_model <- coxph(Surv(survival_time, event_status) ~ Capnocytophaga_cut+age_at_index+
                     alcohol_history+prior_malignancy+
                     portion_weight2+gender.x, data = merged_data2)

cox_model <- coxph(Surv(survival_time, event_status) ~ Capnocytophaga_cut, data = merged_data2)
summary(cox_model)
cox.zph(cox_model)


f1_Capnocytophaga<-ggforest(cox_model,
                           data=merged_data2,
                           main = "Hazard ratio",        # 标题
                           cpositions = c(0.01, 0.1, 0.3), #左边三列的相对距离
                           fontsize = 1, # 字体大小
                           refLabel = "reference", #显示因子的参考水平
                           noDigits = 3)  #小数点位数

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Apan-cancer_analysis_of_the_microbiome/TCGA-HNSC头颈/Capnocytophaga_HNSC_ggforest_cut_多因素.pdf",f1_Capnocytophaga,width =13,height=6)

library(survminer)
sur.cut <- surv_cutpoint(merged_data2, time = "survival_time",
                         event = "event_status",
                         variables = "Capnocytophaga_percentage")
summary(sur.cut)
#对截断值两侧的数据分布可视化。
plot(sur.cut, "Capnocytophaga_percentage", palette = "npg")

#使用surv_categorize函数，根据截断值转为二分类变量。
sur.cat <- surv_categorize(sur.cut, labels = c("High", "Low"))


merged_data2$Capnocytophaga_cut <-
  ifelse(merged_data2$Capnocytophaga_percentage >= 0.2701088, "High", "Low")



fit <- survfit(Surv(survival_time, event_status) ~Capnocytophagas_group,  # 创建生存对象 
               data = merged_data2)


fit <- survfit(Surv(survival_time, event_status) ~Capnocytophaga_cut,  # 创建生存对象 
               data = merged_data2)

pdf("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Apan-cancer_analysis_of_the_microbiome/TCGA-HNSC头颈/Capnocytophagas_HNSC_curvep_cut.pdf", width = 6, height =7)
ggsurvplot(fit, # 创建的拟合对象
           data = merged_data2,  # 指定变量数据来源
           conf.int = TRUE, # 显示置信区间
           pval = TRUE, # 添加P值
           surv.median.line = "hv",  # 添加中位生存时间线
           risk.table = TRUE, # 添加风险表
           ylab="Overall Survival(%)",
           xlab = "Follow up time(d)", # 指定x轴标签
           legend = c(0.8,0.75), # 指定图例位置
           legend.title = "", # 设置图例标题，这里设置不显示标题，用空格替代
           palette = c("skyblue", "red"),
           legend.labs = c("Low", "High"), # 指定图例分组标签
           break.x.by = 1000,
           title = "Capnocytophagas HNSC" )  # 设置x轴刻度间距
dev.off()


#####Streptococcus HNSC----
merged_data2<- merged_data %>% filter(project %in% "TCGA-HNSC")

# 计算 75% 分位点
threshold <- quantile(merged_data2$Streptococcus_percentage, 0.75,na.rm = TRUE)

merged_data2$Streptococcus_group <- ifelse(merged_data2$Streptococcus_percentage <= threshold, "Low", "High")

pp1<- ggplot(merged_data2, aes(x = Streptococcus_percentage, fill = Streptococcus_group)) +
  geom_density(alpha = 1) +  # 密度图，透明度为0.6
  geom_vline(xintercept = threshold, linetype = "dashed") +  # 添加分位点虚线
  scale_fill_manual(values = c("Low" = "skyblue", "High" = "red")) +  # 自定义填充颜色
  labs(
    x = "Streptococcus (%)",
    y = "Density",
    fill = "Streptococcus"
  ) +
  theme_minimal() +  # 设置主题
  theme(
    legend.position = "top",  # 调整图例位置
    axis.ticks.x = element_line(linewidth = 0.8), # 加粗 X 轴刻度线
    axis.ticks.y = element_line(linewidth = 0.8), # 加粗 Y 轴刻度线
    axis.title = element_text(size = 14),   # 调整坐标轴标题的大小
    axis.text = element_text(size = 12),     # 调整坐标轴文字的大小
    axis.line.x = element_line(linewidth = 0.8),  # 加粗 X 轴线条
    axis.line.y = element_line(linewidth = 0.8)  # 加粗 Y 轴线条
  )+
  # scale_x_continuous(expand = expansion(mult = c(0.1, 0.1))) + 
  # scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) + 
  coord_cartesian(xlim = c(0,60), ylim = c(0, 0.3))
pp1

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Apan-cancer_analysis_of_the_microbiome/TCGA-HNSC头颈/TCGA_Streptococcus_density.pdf",pp1,width = 5,height = 5)


merged_data2$Streptococcus_group <- factor(
  merged_data2$Streptococcus_group,
  levels = c("Low", "High") # 包含所有可能值
)
merged_data2$Streptococcus_group <- relevel(merged_data2$Streptococcus_group, ref = "Low")

merged_data2$gender.x <- factor(merged_data2$gender.x)
sum(is.na(merged_data2$gender.x))

names(merged_data2)
unique(merged_data2$years_smoked)
merged_data2$years_smoked_history<-ifelse(merged_data2$years_smoked <= 20, "Low", "High")
merged_data2$portion_weight2 <- log(merged_data2$portion_weight + 1)
merged_data2$years_smoked2 <- log(merged_data2$years_smoked + 1)
merged_data2$years_smoked2[is.na(merged_data2$years_smoked2)] <- mean(merged_data2$years_smoked2, na.rm = TRUE)

cox_model <- coxph(Surv(survival_time, event_status) ~ Streptococcus_cut+age_at_index+
                     alcohol_history+prior_malignancy+
                     portion_weight2+gender.x, data = merged_data2)
cox_model <- coxph(Surv(survival_time, event_status) ~ Streptococcus_cut, data = merged_data2)
summary(cox_model)
cox.zph(cox_model)


f1_Streptococcus<-ggforest(cox_model,
                            data=merged_data2,
                            main = "Hazard ratio",        # 标题
                            cpositions = c(0.01, 0.1, 0.3), #左边三列的相对距离
                            fontsize = 1, # 字体大小
                            refLabel = "reference", #显示因子的参考水平
                            noDigits = 3)  #小数点位数

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Apan-cancer_analysis_of_the_microbiome/TCGA-HNSC头颈/Streptococcus_HNSC_ggforest_cut_多因素.pdf",f1_Streptococcus,width =13,height=6)

library(survminer)
sur.cut <- surv_cutpoint(merged_data2, time = "survival_time",
                         event = "event_status",
                         variables = "Streptococcus_percentage")
summary(sur.cut)
#对截断值两侧的数据分布可视化。
plot(sur.cut, "Streptococcus_percentage", palette = "npg")

#使用surv_categorize函数，根据截断值转为二分类变量。
sur.cat <- surv_categorize(sur.cut, labels = c("High", "Low"))


merged_data2$Streptococcus_cut <-
  ifelse(merged_data2$Streptococcus_percentage >= 35.52593, "High", "Low")


fit <- survfit(Surv(survival_time, event_status) ~Streptococcus_group,  # 创建生存对象 
               data = merged_data2)

fit <- survfit(Surv(survival_time, event_status) ~Streptococcus_cut,  # 创建生存对象 
               data = merged_data2)

summary(fit)

pdf("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Apan-cancer_analysis_of_the_microbiome/TCGA-HNSC头颈/Streptococcus_HNSC_curvep_cut.pdf", width = 6, height =7)
ggsurvplot(fit, # 创建的拟合对象
           data = merged_data2,  # 指定变量数据来源
           conf.int = TRUE, # 显示置信区间
           pval = TRUE, # 添加P值
           surv.median.line = "hv",  # 添加中位生存时间线
           risk.table = TRUE, # 添加风险表
           ylab="Overall Survival(%)",
           xlab = "Follow up time(d)", # 指定x轴标签
           legend = c(0.8,0.75), # 指定图例位置
           legend.title = "", # 设置图例标题，这里设置不显示标题，用空格替代
           palette = c("skyblue", "red"),
           legend.labs = c("Low", "High"), # 指定图例分组标签
           break.x.by = 1000,
           title = "Streptococcus HNSC" )  # 设置x轴刻度间距
dev.off()
dev.new()

# 
# ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Apan-cancer_analysis_of_the_microbiome/TCGA-HNSC头颈/Streptococcus_HNSC_curvep.pdf",plot = curvep$plot,width =6,height=6)



####TCGA-ESCA食管癌----

#####fusobacterium ESCA----
merged_data2<- merged_data %>% filter(project %in% "TCGA-ESCA")


# 计算 75% 分位点
threshold <- quantile(merged_data2$Fusobacterium_percentage, 0.75,na.rm = TRUE)

merged_data2$Fusobacterium_group <- ifelse(merged_data2$Fusobacterium_percentage <= threshold, "Low", "High")

merged_data3<-  merged_data2 %>%
  filter(survival_time <= 800)

pp1<- ggplot(merged_data2, aes(x = Fusobacterium_percentage, fill = Fusobacterium_group)) +
  geom_density(alpha = 1) +  # 密度图，透明度为0.6
  geom_vline(xintercept = threshold, linetype = "dashed") +  # 添加分位点虚线
  scale_fill_manual(values = c("Low" = "skyblue", "High" = "red")) +  # 自定义填充颜色
  labs(
    x = "Fusobacterium (%)",
    y = "Density",
    fill = "Fusobacterium"
  ) +
  theme_minimal() +  # 设置主题
  theme(
    legend.position = "top",  # 调整图例位置
    axis.ticks.x = element_line(linewidth = 0.8), # 加粗 X 轴刻度线
    axis.ticks.y = element_line(linewidth = 0.8), # 加粗 Y 轴刻度线
    axis.title = element_text(size = 14),   # 调整坐标轴标题的大小
    axis.text = element_text(size = 12),     # 调整坐标轴文字的大小
    axis.line.x = element_line(linewidth = 0.8),  # 加粗 X 轴线条
    axis.line.y = element_line(linewidth = 0.8)  # 加粗 Y 轴线条
  )+
  # scale_x_continuous(expand = expansion(mult = c(0.1, 0.1))) + 
  # scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) + 
  coord_cartesian(xlim = c(-0,0.2), ylim = c(0, 30))
pp1

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Apan-cancer_analysis_of_the_microbiome/TCGA-ESCA食管/TCGA_Fusobacterium_density.pdf",pp1,width = 5,height = 5)


merged_data2$Fusobacterium_group <- factor(
  merged_data2$Fusobacterium_group,
  levels = c("Low", "High") # 包含所有可能值
)
merged_data2$Fusobacterium_group <- relevel(merged_data2$Fusobacterium_group, ref = "Low")

merged_data2$gender.x <- factor(merged_data2$gender.x)
sum(is.na(merged_data2$gender.x))

names(merged_data2)
unique(merged_data2$years_smoked)
merged_data2$years_smoked_history<-ifelse(merged_data2$years_smoked <= 20, "Low", "High")
merged_data2$portion_weight2 <- log(merged_data2$portion_weight + 1)
merged_data2$years_smoked2 <- log(merged_data2$years_smoked + 1)
merged_data2$years_smoked2[is.na(merged_data2$years_smoked2)] <- mean(merged_data2$years_smoked2, na.rm = TRUE)

cox_model <- coxph(Surv(survival_time, event_status) ~ Fusobacterium_cut+age_at_index+
                     prior_malignancy+
                     portion_weight2+gender.x, data = merged_data2)

cox_model <- coxph(Surv(survival_time, event_status) ~ Fusobacterium_percentage+age_at_index+
                     prior_malignancy+
                     portion_weight2+gender.x, data = merged_data2)

cox_model <- coxph(Surv(survival_time, event_status) ~ Fusobacterium_cut, data = merged_data2)


summary(cox_model)
cox.zph(cox_model)

library(survminer)

f1_Fusobacterium<-ggforest(cox_model,
                           data=merged_data2,
                           main = "Hazard ratio",        # 标题
                           cpositions = c(0.01, 0.1, 0.3), #左边三列的相对距离
                           fontsize = 1, # 字体大小
                           refLabel = "reference", #显示因子的参考水平
                           noDigits = 3)  #小数点位数

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Apan-cancer_analysis_of_the_microbiome/TCGA-ESCA食管/fusobacterium_ESCA_ggforest_cut_多因素.pdf",f1_Fusobacterium,width =13,height=6)


library(survminer)
sur.cut <- surv_cutpoint(merged_data2, time = "survival_time",
                         event = "event_status",
                         variables = "Fusobacterium_percentage")
summary(sur.cut)
#对截断值两侧的数据分布可视化。
plot(sur.cut, "Fusobacterium_percentage", palette = "npg")

#使用surv_categorize函数，根据截断值转为二分类变量。
sur.cat <- surv_categorize(sur.cut, labels = c("High", "Low"))


merged_data2$Fusobacterium_cut <-
  ifelse(merged_data2$Fusobacterium_percentage >= 0.001484865, "High", "Low")


fit <- survfit(Surv(survival_time, event_status) ~Fusobacterium_group,  # 创建生存对象 
               data = merged_data2)

fit <- survfit(Surv(survival_time, event_status) ~Fusobacterium_cut,  # 创建生存对象 
               data = merged_data2)

merged_data3<-  merged_data2 %>%
  filter(survival_time <= 1000)

fit <- survfit(Surv(survival_time, event_status) ~Fusobacterium_group,  # 创建生存对象 
               data = merged_data3)


fit <- survfit(Surv(survival_time, event_status) ~Fusobacterium_cut,  # 创建生存对象 
               data = merged_data3)

pdf("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Apan-cancer_analysis_of_the_microbiome/TCGA-ESCA食管/Fusobacterium_ESCA_curvep_cut2.pdf", width = 6, height =7)
ggsurvplot(fit, # 创建的拟合对象
           data = merged_data2,  # 指定变量数据来源
           conf.int = TRUE, # 显示置信区间
           pval = TRUE, # 添加P值
           surv.median.line = "hv",  # 添加中位生存时间线
           risk.table = TRUE, # 添加风险表
           ylab="Overall Survival(%)",
           xlab = "Follow up time(d)", # 指定x轴标签
           legend = c(0.8,0.75), # 指定图例位置
           legend.title = "", # 设置图例标题，这里设置不显示标题，用空格替代
           palette = c("skyblue", "red"),
           legend.labs = c("Low", "High"), # 指定图例分组标签
           break.x.by = 1000,
           title = "Fusobacterium ESCA" )  # 设置x轴刻度间距
dev.off()


#####Capnocytophaga ESCA----
merged_data2<- merged_data %>% filter(project %in% "TCGA-ESCA")

# 计算 75% 分位点
threshold <- quantile(merged_data2$Capnocytophaga_percentage, 0.75,na.rm = TRUE)

merged_data2$Capnocytophaga_group <- ifelse(merged_data2$Capnocytophaga_percentage <= threshold, "Low", "High")

pp1<- ggplot(merged_data2, aes(x = Capnocytophaga_percentage, fill = Capnocytophaga_group)) +
  geom_density(alpha = 1) +  # 密度图，透明度为0.6
  geom_vline(xintercept = threshold, linetype = "dashed") +  # 添加分位点虚线
  scale_fill_manual(values = c("Low" = "skyblue", "High" = "red")) +  # 自定义填充颜色
  labs(
    x = "Capnocytophaga (%)",
    y = "Density",
    fill = "Capnocytophaga"
  ) +
  theme_minimal() +  # 设置主题
  theme(
    legend.position = "top",  # 调整图例位置
    axis.ticks.x = element_line(linewidth = 0.8), # 加粗 X 轴刻度线
    axis.ticks.y = element_line(linewidth = 0.8), # 加粗 Y 轴刻度线
    axis.title = element_text(size = 14),   # 调整坐标轴标题的大小
    axis.text = element_text(size = 12),     # 调整坐标轴文字的大小
    axis.line.x = element_line(linewidth = 0.8),  # 加粗 X 轴线条
    axis.line.y = element_line(linewidth = 0.8)  # 加粗 Y 轴线条
  )+
  # scale_x_continuous(expand = expansion(mult = c(0.1, 0.1))) + 
  # scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) + 
  coord_cartesian(xlim = c(-0.1,0.02), ylim = c(0, 0.75))
pp1

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Apan-cancer_analysis_of_the_microbiome/TCGA-ESCA食管/TCGA_Capnocytophaga_density.pdf",pp1,width = 5,height = 5)


merged_data2$Capnocytophaga_group <- factor(
  merged_data2$Capnocytophaga_group,
  levels = c("Low", "High") # 包含所有可能值
)
merged_data2$Capnocytophaga_group <- relevel(merged_data2$Capnocytophaga_group, ref = "Low")

merged_data2$gender.x <- factor(merged_data2$gender.x)
sum(is.na(merged_data2$gender.x))

names(merged_data2)
unique(merged_data2$years_smoked)
merged_data2$years_smoked_history<-ifelse(merged_data2$years_smoked <= 20, "Low", "High")
merged_data2$portion_weight2 <- log(merged_data2$portion_weight + 1)
merged_data2$years_smoked2 <- log(merged_data2$years_smoked + 1)
merged_data2$years_smoked2[is.na(merged_data2$years_smoked2)] <- mean(merged_data2$years_smoked2, na.rm = TRUE)

cox_model <- coxph(Surv(survival_time, event_status) ~ Capnocytophaga_cut+age_at_index+
                     prior_malignancy+
                     portion_weight2+gender.x, data = merged_data2)

cox_model <- coxph(Surv(survival_time, event_status) ~ Capnocytophaga_group+age_at_index+
                     prior_malignancy+
                     portion_weight2+gender.x, data = merged_data2)

cox_model <- coxph(Surv(survival_time, event_status) ~ Capnocytophaga_cut, data = merged_data2)

summary(cox_model)
cox.zph(cox_model)

dev.new()
f1_Capnocytophaga<-ggforest(cox_model,
                            data=merged_data2,
                            main = "Hazard ratio",        # 标题
                            cpositions = c(0.01, 0.1, 0.3), #左边三列的相对距离
                            fontsize = 1, # 字体大小
                            refLabel = "reference", #显示因子的参考水平
                            noDigits = 3)  #小数点位数

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Apan-cancer_analysis_of_the_microbiome/TCGA-ESCA食管/Capnocytophaga_ESCA_ggforest_cut_多因素.pdf",f1_Capnocytophaga,width =13,height=6)

library(survminer)
sur.cut <- surv_cutpoint(merged_data2, time = "survival_time",
                         event = "event_status",
                         variables = "Capnocytophaga_percentage")
summary(sur.cut)
#对截断值两侧的数据分布可视化。
plot(sur.cut, "Capnocytophaga_percentage", palette = "npg")

#使用surv_categorize函数，根据截断值转为二分类变量。
sur.cat <- surv_categorize(sur.cut, labels = c("High", "Low"))


merged_data2$Capnocytophaga_cut <-
  ifelse(merged_data2$Capnocytophaga_percentage >= 0.001608721, "High", "Low")



fit <- survfit(Surv(survival_time, event_status) ~Capnocytophaga_group,  # 创建生存对象 
               data = merged_data2)


fit <- survfit(Surv(survival_time, event_status) ~Capnocytophaga_cut,  # 创建生存对象 
               data = merged_data2)

merged_data3<-  merged_data2 %>%
  filter(survival_time <= 800)

fit <- survfit(Surv(survival_time, event_status) ~Capnocytophaga_group,  # 创建生存对象 
               data = merged_data3)

fit <- survfit(Surv(survival_time, event_status) ~Capnocytophaga_cut,  # 创建生存对象 
               data = merged_data3)


pdf("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Apan-cancer_analysis_of_the_microbiome/TCGA-ESCA食管/Capnocytophagas_ESCA_curvep_cut1000.pdf", width = 6, height =7)
ggsurvplot(fit, # 创建的拟合对象
           data = merged_data3,  # 指定变量数据来源
           conf.int = TRUE, # 显示置信区间
           pval = TRUE, # 添加P值
           surv.median.line = "hv",  # 添加中位生存时间线
           risk.table = TRUE, # 添加风险表
           ylab="Overall Survival(%)",
           xlab = "Follow up time(d)", # 指定x轴标签
           legend = c(0.8,0.75), # 指定图例位置
           legend.title = "", # 设置图例标题，这里设置不显示标题，用空格替代
           palette = c("skyblue", "red"),
           legend.labs = c("Low", "High"), # 指定图例分组标签
           break.x.by = 200,
           title = "Capnocytophagas ESCA" )  # 设置x轴刻度间距
dev.off()


#####Streptococcus ESCA----
merged_data2<- merged_data %>% filter(project %in% "TCGA-ESCA")

# 计算 75% 分位点
threshold <- quantile(merged_data2$Streptococcus_percentage, 0.75,na.rm = TRUE)

merged_data2$Streptococcus_group <- ifelse(merged_data2$Streptococcus_percentage <= threshold, "Low", "High")

pp1<- ggplot(merged_data2, aes(x = Streptococcus_percentage, fill = Streptococcus_group)) +
  geom_density(alpha = 0.6) +  # 密度图，透明度为0.6
  geom_vline(xintercept = threshold, linetype = "dashed") +  # 添加分位点虚线
  scale_fill_manual(values = c("Low" = "skyblue", "High" = "red")) +  # 自定义填充颜色
  labs(
    x = "Streptococcus (%)",
    y = "Density",
    fill = "Streptococcus"
  ) +
  theme_minimal() +  # 设置主题
  theme(
    legend.position = "top",  # 调整图例位置
    axis.ticks.x = element_line(linewidth = 0.8), # 加粗 X 轴刻度线
    axis.ticks.y = element_line(linewidth = 0.8), # 加粗 Y 轴刻度线
    axis.title = element_text(size = 14),   # 调整坐标轴标题的大小
    axis.text = element_text(size = 12),     # 调整坐标轴文字的大小
    axis.line.x = element_line(linewidth = 0.8),  # 加粗 X 轴线条
    axis.line.y = element_line(linewidth = 0.8)  # 加粗 Y 轴线条
  )+
  # scale_x_continuous(expand = expansion(mult = c(0.1, 0.1))) + 
  # scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) + 
  coord_cartesian(xlim = c(0,60), ylim = c(0, 0.3))
pp1

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Apan-cancer_analysis_of_the_microbiome/TCGA-ESCA食管/TCGA_Streptococcus_density.pdf",pp1,width = 5,height = 5)


merged_data2$Streptococcus_group <- factor(
  merged_data2$Streptococcus_group,
  levels = c("Low", "High") # 包含所有可能值
)
merged_data2$Streptococcus_group <- relevel(merged_data2$Streptococcus_group, ref = "Low")

merged_data2$gender.x <- factor(merged_data2$gender.x)
sum(is.na(merged_data2$gender.x))

names(merged_data2)
unique(merged_data2$years_smoked)
merged_data2$years_smoked_history<-ifelse(merged_data2$years_smoked <= 20, "Low", "High")
merged_data2$portion_weight2 <- log(merged_data2$portion_weight + 1)
merged_data2$years_smoked2 <- log(merged_data2$years_smoked + 1)
merged_data2$years_smoked2[is.na(merged_data2$years_smoked2)] <- mean(merged_data2$years_smoked2, na.rm = TRUE)

cox_model <- coxph(Surv(survival_time, event_status) ~ Streptococcus_cut+age_at_index+
                     prior_malignancy+
                     portion_weight2+gender.x, data = merged_data2)
cox_model <- coxph(Surv(survival_time, event_status) ~ Streptococcus_cut, data = merged_data2)
summary(cox_model)
cox.zph(cox_model)


f1_Streptococcus<-ggforest(cox_model,
                           data=merged_data2,
                           main = "Hazard ratio",        # 标题
                           cpositions = c(0.01, 0.1, 0.3), #左边三列的相对距离
                           fontsize = 1, # 字体大小
                           refLabel = "reference", #显示因子的参考水平
                           noDigits = 3)  #小数点位数

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Apan-cancer_analysis_of_the_microbiome/TCGA-ESCA食管/Streptococcus_ESCA_ggforest_cut_多因素.pdf",f1_Streptococcus,width =13,height=6)

library(survminer)
sur.cut <- surv_cutpoint(merged_data2, time = "survival_time",
                         event = "event_status",
                         variables = "Streptococcus_percentage")
summary(sur.cut)
#对截断值两侧的数据分布可视化。
plot(sur.cut, "Streptococcus_percentage", palette = "npg")

#使用surv_categorize函数，根据截断值转为二分类变量。
sur.cat <- surv_categorize(sur.cut, labels = c("High", "Low"))


merged_data2$Streptococcus_cut <-
  ifelse(merged_data2$Streptococcus_percentage >= 23.19544, "High", "Low")



fit <- survfit(Surv(survival_time, event_status) ~Streptococcus_group,  # 创建生存对象 
               data = merged_data2)

fit <- survfit(Surv(survival_time, event_status) ~Streptococcus_cut,  # 创建生存对象 
               data = merged_data2)

summary(fit)

merged_data3<-  merged_data2 %>%
  filter(survival_time <= 1000)

fit <- survfit(Surv(survival_time, event_status) ~Streptococcus_group,  # 创建生存对象 
               data = merged_data3)

fit <- survfit(Surv(survival_time, event_status) ~Streptococcus_cut,  # 创建生存对象 
               data = merged_data3)

pdf("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Apan-cancer_analysis_of_the_microbiome/TCGA-ESCA食管/Streptococcus_ESCA_curvep_cut.pdf", width = 6, height =7)
ggsurvplot(fit, # 创建的拟合对象
           data = merged_data2,  # 指定变量数据来源
           conf.int = TRUE, # 显示置信区间
           pval = TRUE, # 添加P值
           surv.median.line = "hv",  # 添加中位生存时间线
           risk.table = TRUE, # 添加风险表
           ylab="Overall Survival(%)",
           xlab = "Follow up time(d)", # 指定x轴标签
           legend = c(0.8,0.75), # 指定图例位置
           legend.title = "", # 设置图例标题，这里设置不显示标题，用空格替代
           palette = c("skyblue", "red"),
           legend.labs = c("Low", "High"), # 指定图例分组标签
           break.x.by = 1000,
           title = "Streptococcus ESCA" )  # 设置x轴刻度间距
dev.off()
dev.new()

# 
# ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Apan-cancer_analysis_of_the_microbiome/TCGA-ESCA头颈/Streptococcus_ESCA_curvep.pdf",plot = curvep$plot,width =6,height=6)


####TCGA-LUSC – 非小细胞肺癌（鳞癌） (Lung squamous cell carcinoma)----

#####fusobacterium LUSC----
merged_data2<- merged_data %>% filter(project %in% "TCGA-LUSC")

# 计算 75% 分位点
threshold <- quantile(merged_data2$Fusobacterium_percentage, 0.75,na.rm = TRUE)

merged_data2$Fusobacterium_group <- ifelse(merged_data2$Fusobacterium_percentage <= threshold, "Low", "High")

pp1<- ggplot(merged_data2, aes(x = Fusobacterium_percentage, fill = Fusobacterium_group)) +
  geom_density(alpha = 1) +  # 密度图，透明度为0.6
  geom_vline(xintercept = threshold, linetype = "dashed") +  # 添加分位点虚线
  scale_fill_manual(values = c("Low" = "skyblue", "High" = "red")) +  # 自定义填充颜色
  labs(
    x = "Fusobacterium (%)",
    y = "Density",
    fill = "Fusobacterium"
  ) +
  theme_minimal() +  # 设置主题
  theme(
    legend.position = "top",  # 调整图例位置
    axis.ticks.x = element_line(linewidth = 0.8), # 加粗 X 轴刻度线
    axis.ticks.y = element_line(linewidth = 0.8), # 加粗 Y 轴刻度线
    axis.title = element_text(size = 14),   # 调整坐标轴标题的大小
    axis.text = element_text(size = 12),     # 调整坐标轴文字的大小
    axis.line.x = element_line(linewidth = 0.8),  # 加粗 X 轴线条
    axis.line.y = element_line(linewidth = 0.8)  # 加粗 Y 轴线条
  )+
  # scale_x_continuous(expand = expansion(mult = c(0.1, 0.1))) + 
  # scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) + 
  coord_cartesian(xlim = c(-0,0.006), ylim = c(0, 30))
pp1

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Apan-cancer_analysis_of_the_microbiome/TCGA-LUSC肺麟癌/TCGA_Fusobacterium_density.pdf",pp1,width = 5,height = 5)


merged_data2$Fusobacterium_group <- factor(
  merged_data2$Fusobacterium_group,
  levels = c("Low", "High") # 包含所有可能值
)
merged_data2$Fusobacterium_group <- relevel(merged_data2$Fusobacterium_group, ref = "Low")

merged_data2$gender.x <- factor(merged_data2$gender.x)
sum(is.na(merged_data2$gender.x))

names(merged_data2)
unique(merged_data2$years_smoked)
merged_data2$years_smoked_history<-ifelse(merged_data2$years_smoked <= 20, "Low", "High")
merged_data2$portion_weight2 <- log(merged_data2$portion_weight + 1)
merged_data2$years_smoked2 <- log(merged_data2$years_smoked + 1)
merged_data2$years_smoked2[is.na(merged_data2$years_smoked2)] <- mean(merged_data2$years_smoked2, na.rm = TRUE)
library(survminer)

cox_model <- coxph(Surv(survival_time, event_status) ~ Fusobacterium_group+age_at_index+
                     prior_malignancy+
                     portion_weight2+gender.x, data = merged_data2)

cox_model <- coxph(Surv(survival_time, event_status) ~ Fusobacterium_cut+age_at_index+
                     prior_malignancy+
                     portion_weight2+gender.x, data = merged_data2)
cox_model <- coxph(Surv(survival_time, event_status) ~ Fusobacterium_cut, data = merged_data2)

summary(cox_model)
cox.zph(cox_model)


f1_Fusobacterium<-ggforest(cox_model,
                           data=merged_data2,
                           main = "Hazard ratio",        # 标题
                           cpositions = c(0.01, 0.1, 0.3), #左边三列的相对距离
                           fontsize = 1, # 字体大小
                           refLabel = "reference", #显示因子的参考水平
                           noDigits = 3)  #小数点位数

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Apan-cancer_analysis_of_the_microbiome/TCGA-LUSC肺鳞癌/fusobacterium_LUSC_ggforest_cut_多因素.pdf",f1_Fusobacterium,width =13,height=6)

library(survminer)
sur.cut <- surv_cutpoint(merged_data2, time = "survival_time",
                         event = "event_status",
                         variables = "Fusobacterium_percentage")
summary(sur.cut)
#对截断值两侧的数据分布可视化。
plot(sur.cut, "Fusobacterium_percentage", palette = "npg")

#使用surv_categorize函数，根据截断值转为二分类变量。
sur.cat <- surv_categorize(sur.cut, labels = c("High", "Low"))


merged_data2$Fusobacterium_cut <-
  ifelse(merged_data2$Fusobacterium_percentage >= 0.005432302, "High", "Low")



fit <- survfit(Surv(survival_time, event_status) ~Fusobacterium_group,  # 创建生存对象 
               data = merged_data2)

fit <- survfit(Surv(survival_time, event_status) ~Fusobacterium_cut,  # 创建生存对象 
               data = merged_data2)
dev.new()
pdf("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Apan-cancer_analysis_of_the_microbiome/TCGA-LUSC肺鳞癌/Fusobacterium_LUSC_curvep_cut.pdf", width = 6, height =7)
ggsurvplot(fit, # 创建的拟合对象
           data = merged_data2,  # 指定变量数据来源
           conf.int = TRUE, # 显示置信区间
           pval = TRUE, # 添加P值
           surv.median.line = "hv",  # 添加中位生存时间线
           risk.table = TRUE, # 添加风险表
           ylab="Overall Survival(%)",
           xlab = "Follow up time(d)", # 指定x轴标签
           legend = c(0.8,0.75), # 指定图例位置
           legend.title = "", # 设置图例标题，这里设置不显示标题，用空格替代
           palette = c("skyblue", "red"),
           legend.labs = c("Low", "High"), # 指定图例分组标签
           break.x.by = 1000,
           title = "Fusobacterium LUSC" )  # 设置x轴刻度间距
dev.off()


#####Capnocytophaga LUSC----
merged_data2<- merged_data %>% filter(project %in% "TCGA-LUSC")

# 计算 75% 分位点
threshold <- quantile(merged_data2$Capnocytophaga_percentage, 0.75,na.rm = TRUE)

merged_data2$Capnocytophaga_group <- ifelse(merged_data2$Capnocytophaga_percentage <= threshold, "Low", "High")

pp1<- ggplot(merged_data2, aes(x = Capnocytophaga_percentage, fill = Capnocytophaga_group)) +
  geom_density(alpha = 1) +  # 密度图，透明度为0.6
  geom_vline(xintercept = threshold, linetype = "dashed") +  # 添加分位点虚线
  scale_fill_manual(values = c("Low" = "skyblue", "High" = "red")) +  # 自定义填充颜色
  labs(
    x = "Capnocytophaga (%)",
    y = "Density",
    fill = "Capnocytophaga"
  ) +
  theme_minimal() +  # 设置主题
  theme(
    legend.position = "top",  # 调整图例位置
    axis.ticks.x = element_line(linewidth = 0.8), # 加粗 X 轴刻度线
    axis.ticks.y = element_line(linewidth = 0.8), # 加粗 Y 轴刻度线
    axis.title = element_text(size = 14),   # 调整坐标轴标题的大小
    axis.text = element_text(size = 12),     # 调整坐标轴文字的大小
    axis.line.x = element_line(linewidth = 0.8),  # 加粗 X 轴线条
    axis.line.y = element_line(linewidth = 0.8)  # 加粗 Y 轴线条
  )+
  # scale_x_continuous(expand = expansion(mult = c(0.1, 0.1))) + 
  # scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) + 
  coord_cartesian(xlim = c(-0.1,10), ylim = c(0, 0.75))
pp1

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Apan-cancer_analysis_of_the_microbiome/TCGA-LUSC肺鳞癌/TCGA_Capnocytophaga_density.pdf",pp1,width = 5,height = 5)


merged_data2$Capnocytophaga_group <- factor(
  merged_data2$Capnocytophaga_group,
  levels = c("Low", "High") # 包含所有可能值
)
merged_data2$Capnocytophaga_group <- relevel(merged_data2$Capnocytophaga_group, ref = "Low")

merged_data2$gender.x <- factor(merged_data2$gender.x)
sum(is.na(merged_data2$gender.x))

names(merged_data2)
unique(merged_data2$years_smoked)
merged_data2$years_smoked_history<-ifelse(merged_data2$years_smoked <= 20, "Low", "High")
merged_data2$portion_weight2 <- log(merged_data2$portion_weight + 1)
merged_data2$years_smoked2 <- log(merged_data2$years_smoked + 1)
merged_data2$years_smoked2[is.na(merged_data2$years_smoked2)] <- mean(merged_data2$years_smoked2, na.rm = TRUE)

cox_model <- coxph(Surv(survival_time, event_status) ~ Capnocytophaga_percentage+age_at_index+
                     prior_malignancy+
                     portion_weight2+gender.x, data = merged_data2)

cox_model <- coxph(Surv(survival_time, event_status) ~ Capnocytophaga_cut+age_at_index+
                     prior_malignancy+
                     portion_weight2+gender.x, data = merged_data2)
cox_model <- coxph(Surv(survival_time, event_status) ~ Capnocytophaga_cut, data = merged_data2)
summary(cox_model)
cox.zph(cox_model)


f1_Capnocytophaga<-ggforest(cox_model,
                            data=merged_data2,
                            main = "Hazard ratio",        # 标题
                            cpositions = c(0.01, 0.1, 0.3), #左边三列的相对距离
                            fontsize = 1, # 字体大小
                            refLabel = "reference", #显示因子的参考水平
                            noDigits = 3)  #小数点位数

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Apan-cancer_analysis_of_the_microbiome/TCGA-LUSC肺鳞癌/Capnocytophaga_LUSC_ggforest_cut_多因素.pdf",f1_Capnocytophaga,width =13,height=6)


library(survminer)
sur.cut <- surv_cutpoint(merged_data2, time = "survival_time",
                         event = "event_status",
                         variables = "Capnocytophaga_percentage")
summary(sur.cut)
#对截断值两侧的数据分布可视化。
plot(sur.cut, "Capnocytophaga_percentage", palette = "npg")

#使用surv_categorize函数，根据截断值转为二分类变量。
sur.cat <- surv_categorize(sur.cut, labels = c("High", "Low"))


merged_data2$Capnocytophaga_cut <-
  ifelse(merged_data2$Capnocytophaga_percentage >= 0.002396669, "High", "Low")


fit <- survfit(Surv(survival_time, event_status) ~Capnocytophaga_group,  # 创建生存对象 
               data = merged_data2)

fit <- survfit(Surv(survival_time, event_status) ~Capnocytophaga_cut,  # 创建生存对象 
               data = merged_data2)

pdf("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Apan-cancer_analysis_of_the_microbiome/TCGA-LUSC肺鳞癌/Capnocytophagas_LUSC_curvep_cut.pdf", width = 6, height =7)
ggsurvplot(fit, # 创建的拟合对象
           data = merged_data2,  # 指定变量数据来源
           conf.int = TRUE, # 显示置信区间
           pval = TRUE, # 添加P值
           surv.median.line = "hv",  # 添加中位生存时间线
           risk.table = TRUE, # 添加风险表
           ylab="Overall Survival(%)",
           xlab = "Follow up time(d)", # 指定x轴标签
           legend = c(0.8,0.75), # 指定图例位置
           legend.title = "", # 设置图例标题，这里设置不显示标题，用空格替代
           palette = c("skyblue", "red"),
           legend.labs = c("Low", "High"), # 指定图例分组标签
           break.x.by = 1000,
           title = "Capnocytophagas LUSC" )  # 设置x轴刻度间距
dev.off()


#####Streptococcus LUSC----
merged_data2<- merged_data %>% filter(project %in% "TCGA-LUSC")

# 计算 75% 分位点
threshold <- quantile(merged_data2$Streptococcus_percentage, 0.75,na.rm = TRUE)

merged_data2$Streptococcus_group <- ifelse(merged_data2$Streptococcus_percentage <= threshold, "Low", "High")

pp1<- ggplot(merged_data2, aes(x = Streptococcus_percentage, fill = Streptococcus_group)) +
  geom_density(alpha = 1) +  # 密度图，透明度为0.6
  geom_vline(xintercept = threshold, linetype = "dashed") +  # 添加分位点虚线
  scale_fill_manual(values = c("Low" = "skyblue", "High" = "red")) +  # 自定义填充颜色
  labs(
    x = "Streptococcus (%)",
    y = "Density",
    fill = "Streptococcus"
  ) +
  theme_minimal() +  # 设置主题
  theme(
    legend.position = "top",  # 调整图例位置
    axis.ticks.x = element_line(linewidth = 0.8), # 加粗 X 轴刻度线
    axis.ticks.y = element_line(linewidth = 0.8), # 加粗 Y 轴刻度线
    axis.title = element_text(size = 14),   # 调整坐标轴标题的大小
    axis.text = element_text(size = 12),     # 调整坐标轴文字的大小
    axis.line.x = element_line(linewidth = 0.8),  # 加粗 X 轴线条
    axis.line.y = element_line(linewidth = 0.8)  # 加粗 Y 轴线条
  )+
  # scale_x_continuous(expand = expansion(mult = c(0.1, 0.1))) + 
  # scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) + 
  coord_cartesian(xlim = c(0,60), ylim = c(0, 0.3))
pp1

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Apan-cancer_analysis_of_the_microbiome/TCGA-LUSC肺鳞癌/TCGA_Streptococcus_density.pdf",pp1,width = 5,height = 5)


merged_data2$Streptococcus_group <- factor(
  merged_data2$Streptococcus_group,
  levels = c("Low", "High") # 包含所有可能值
)
merged_data2$Streptococcus_group <- relevel(merged_data2$Streptococcus_group, ref = "Low")

merged_data2$gender.x <- factor(merged_data2$gender.x)
sum(is.na(merged_data2$gender.x))

names(merged_data2)
unique(merged_data2$years_smoked)
merged_data2$years_smoked_history<-ifelse(merged_data2$years_smoked <= 20, "Low", "High")
merged_data2$portion_weight2 <- log(merged_data2$portion_weight + 1)
merged_data2$years_smoked2 <- log(merged_data2$years_smoked + 1)
merged_data2$years_smoked2[is.na(merged_data2$years_smoked2)] <- mean(merged_data2$years_smoked2, na.rm = TRUE)

cox_model <- coxph(Surv(survival_time, event_status) ~ Streptococcus_cut+age_at_index+
                     prior_malignancy+
                     portion_weight2+gender.x, data = merged_data2)

cox_model <- coxph(Surv(survival_time, event_status) ~ Streptococcus_cut, data = merged_data2)
summary(cox_model)
cox.zph(cox_model)


f1_Streptococcus<-ggforest(cox_model,
                           data=merged_data2,
                           main = "Hazard ratio",        # 标题
                           cpositions = c(0.01, 0.1, 0.3), #左边三列的相对距离
                           fontsize = 1, # 字体大小
                           refLabel = "reference", #显示因子的参考水平
                           noDigits = 3)  #小数点位数

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Apan-cancer_analysis_of_the_microbiome/TCGA-LUSC肺鳞癌/Streptococcus_LUSC_ggforest_cut_多因素.pdf",f1_Streptococcus,width =13,height=6)

library(survminer)
sur.cut <- surv_cutpoint(merged_data2, time = "survival_time",
                         event = "event_status",
                         variables = "Streptococcus_percentage")
summary(sur.cut)
#对截断值两侧的数据分布可视化。
plot(sur.cut, "Streptococcus_percentage", palette = "npg")

#使用surv_categorize函数，根据截断值转为二分类变量。
sur.cat <- surv_categorize(sur.cut, labels = c("High", "Low"))


merged_data2$Streptococcus_cut <-
  ifelse(merged_data2$Streptococcus_percentage >= 22.83706, "High", "Low")

fit <- survfit(Surv(survival_time, event_status) ~Streptococcus_group,  # 创建生存对象 
               data = merged_data2)

fit <- survfit(Surv(survival_time, event_status) ~Streptococcus_cut,  # 创建生存对象 
               data = merged_data2)

summary(fit)

pdf("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Apan-cancer_analysis_of_the_microbiome/TCGA-LUSC肺鳞癌/Streptococcus_LUSC_curvep_cut.pdf", width = 6, height =7)
ggsurvplot(fit, # 创建的拟合对象
           data = merged_data2,  # 指定变量数据来源
           conf.int = TRUE, # 显示置信区间
           pval = TRUE, # 添加P值
           surv.median.line = "hv",  # 添加中位生存时间线
           risk.table = TRUE, # 添加风险表
           ylab="Overall Survival(%)",
           xlab = "Follow up time(d)", # 指定x轴标签
           legend = c(0.8,0.75), # 指定图例位置
           legend.title = "", # 设置图例标题，这里设置不显示标题，用空格替代
           palette = c("skyblue", "red"),
           legend.labs = c("Low", "High"), # 指定图例分组标签
           break.x.by = 1000,
           title = "Streptococcus LUSC" )  # 设置x轴刻度间距
dev.off()
dev.new()







#############################################################
#############################################################

##git clone https://github.com/twbattaglia/tumor-microbiome

install.packages("renv")
renv::restore()


#########################################################
########################################################
data("Hartwig")
nsamples(Hartwig)
sample_data(Hartwig) %>% 
  head()

aerophilicity = read_csv("resources/aerophilicity.csv") %>% 
  filter(Score > 0.5) %>% 
  mutate(NCBI_ID = as.character(NCBI_ID)) %>% 
  mutate(Attribute2 = case_when(
    Attribute %in% c(":aerobic", ":obligately_aerobic") ~ "Aerobic",
    Attribute %in% c(":anaerobic", ":obligately_anaerobic") ~ "Anaerobic",
    Attribute %in% c(":facultatively_anaerobic") ~ "Facultatively anaerobic",
    Attribute %in% c(":microaerophilic") ~ "microaerophilic",
    Attribute %in% c("missing") ~ "missing"
  ))

gram_status = read_csv("resources/gram_status.csv") %>% 
  filter(Score > 0.5) %>% 
  mutate(NCBI_ID = as.character(NCBI_ID))

pseq.fil.taxtable = tax_table(Hartwig) %>% 
  as.data.frame() %>% 
  rownames_to_column("taxId")



###################################################
################################################
########Data from: The cancer microbiome atlas (TCMA):
####A pan-cancer comparative analysis to distinguish organ-associated microbiota from contaminants.----

####读取数据WGS----
library(dplyr)
library(phyloseq)
physeq.WGS.solid.case.relabund<- readRDS("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/The_cancer_microbiome/tb09j6496/physeq.WXS.solid.case.relabund.rds")
is.null(tax_table(physeq.WGS.solid.case.relabund))
View(physeq.WGS.solid.case.relabund)
otu_table <- as.data.frame(phyloseq::otu_table(physeq.WGS.solid.case.relabund))

otu_table_data <- otu_table(physeq.WGS.solid.case.relabund)
tax_table_data <- phyloseq::tax_table(physeq.WGS.solid.case.relabund)

####Fusobacterium hight low 分组密度图绘制----

# 提取 Fusobacterium 的丰度

fusobacterium_OTUs <- rownames(tax_table_data)[tax_table_data[, "genus"] == "Capnocytophaga"]

unique(tax_table_data[, "genus"])
# 提取OTU_table中与Fusobacterium属相关的丰度数据
fusobacterium_abundance <- otu_table_data[fusobacterium_OTUs, , drop = FALSE]
View(fusobacterium_abundance)
dim(fusobacterium_abundance)
dim(otu_table_data)

# 计算每个样本的Fusobacterium丰度总和
fusobacterium_total <- colSums(fusobacterium_abundance)
length(fusobacterium_total)
# 计算每个样本中Fusobacterium的百分比
total_abundance <- colSums(otu_table_data)  # 每个样本的总丰度
length(total_abundance)
fusobacterium_percentage <- fusobacterium_total / total_abundance * 100

fusobacterium_percentage<- as.data.frame(fusobacterium_percentage)
data <- data.frame(
  Sample = rownames(fusobacterium_percentage),
  Percentage = unlist(fusobacterium_percentage)
)


library(ggplot2)
library(dplyr)
set.seed(123)

# 计算 75% 分位点
threshold <- quantile(data$Percentage, 0.75,na.rm = TRUE)

# 添加分组信息
data$fusobacterium_group <- ifelse(is.na(data$Percentage), "NaN",
                     ifelse(data$Percentage >= threshold, "High", "Low"))

unique(data$fusobacterium_group)

View(data)
pp1<- ggplot(data, aes(x = Percentage, fill = fusobacterium_group)) +
  geom_density(alpha = 1) +  # 密度图，透明度为0.6
  geom_vline(xintercept = threshold, linetype = "dashed") +  # 添加分位点虚线
  scale_fill_manual(values = c("Low" = "skyblue", "High" = "red","NaN"="skyblue")) +  # 自定义填充颜色
  labs(
    x = "Fusobacterium (%)",
    y = "Density",
    fill = "Fusobacterium"
  ) +
  theme_minimal() +  # 设置主题
  theme(
    legend.position = "top",  # 调整图例位置
    axis.ticks.x = element_line(linewidth = 0.8), # 加粗 X 轴刻度线
    axis.ticks.y = element_line(linewidth = 0.8), # 加粗 Y 轴刻度线
    axis.title = element_text(size = 14),   # 调整坐标轴标题的大小
    axis.text = element_text(size = 12),     # 调整坐标轴文字的大小
    axis.line.x = element_line(linewidth = 0.8),  # 加粗 X 轴线条
    axis.line.y = element_line(linewidth = 0.8)  # 加粗 Y 轴线条
  )+
  # scale_x_continuous(expand = expansion(mult = c(0.1, 0.1))) + 
  # scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) + 
  coord_cartesian(xlim = c(-1, 80), ylim = c(0, 0.3)) 
pp1

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/The_cancer_microbiome/Fusobacterium_density_physeq.WGS.solid.case.relabund.pdf",pp1,width = 5,height = 5)

###########################################
####预后分析----


taxa_table <-physeq.WGS.solid.case.relabund@tax_table %>% as.data.frame()

taxa_table <- as(sample_data(taxa_table), "data.frame")


names(taxa_table)

sample_table <- physeq.WGS.solid.case.relabund@sam_data %>% as.data.frame()
dim(sample_table)
View(sample_table)
df <- sample_table
df$survival_time <- ifelse(
  df$vital_status == "Dead", 
  df$days_to_death, 
  df$days_to_last_followup
)

df$disease_code
df$event_status <- ifelse(df$vital_status == "Dead", 1, 0)
names(df)
head(df)
df <- as(sample_data(df), "data.frame")

df <- df %>%
  tibble::rownames_to_column(var = "SampleID")

View(df)

survival_data <- df %>%
  dplyr::select(SampleID, survival_time, event_status,gender,
                age_at_initial_pathologic_diagnosis,
                height,weight,
                tobacco_smoking_history,venous_invasion,lymphatic_invasion,residual_tumor,
                tumor_tissue_site,
                history_of_neoadjuvant_treatment,histological_type)
View(survival_data)
unique(survival_data$tobacco_smoking_history)


survival_data_numeric <- survival_data %>%
  mutate(
    gender=ifelse(gender == "MALE", 1, 0),
    tobacco_smoking_history <- case_when(
      tobacco_smoking_history == "Lifelong Non-smoker" ~ 0,
      tobacco_smoking_history == "Current reformed smoker for <= 15 years" ~ 1,
      tobacco_smoking_history == "Current reformed smoker for > 15 years" ~ 2,
      tobacco_smoking_history == "Current smoker" ~ 3
    ),
    venous_invasion = ifelse(venous_invasion == "True", 1, 0),
    lymphatic_invasion = ifelse(lymphatic_invasion == "True", 1, 0),
    residual_tumor = case_when(
      residual_tumor == "R0" ~ 0,
      residual_tumor == "R1" ~ 1,
      residual_tumor == "R2" ~ 2,
      residual_tumor == "RX" ~ 3
    )
  )
survival_data_numeric <- survival_data_numeric[, -ncol(survival_data_numeric)]
View(survival_data_numeric )
names(data)  
data$Sample
names(survival_data_numeric)
names(df)
unique(df$histological_type)

merged_data <- merge(data, df, 
                     by.x = "Sample", by.y = "SampleID", 
                     all = FALSE)

names(merged_data)

merged_data2<- merged_data %>% filter(tumor_tissue_site %in% "Head and Neck")
View(merged_data2)
names(merged_data2)
names(merged_data2)

survival_data_numeric
library(survival)


merged_data2$fusobacterium_group <- na_if(merged_data2$fusobacterium_group, "NaN")
merged_data2$fusobacterium_group <- factor(
  merged_data2$fusobacterium_group,
  levels = c("Low", "High", "NaN") # 包含所有可能值
)
merged_data2$fusobacterium_group <- relevel(merged_data2$fusobacterium_group, ref = "Low")


merged_data$fusobacterium_cut <- na_if(merged_data$fusobacterium_cut, "NaN")
merged_data$fusobacterium_cut <- factor(
  merged_data$fusobacterium_cut,
  levels = c("Low", "High", "NaN") # 包含所有可能值
)
merged_data$fusobacterium_cut <- relevel(merged_data$fusobacterium_cut, ref = "Low")

cox_model <- coxph(Surv(survival_time, event_status) ~ fusobacterium_cut, data = merged_data)

summary(cox_model)
cox.zph(cox_model)


library(survminer)
sur.cut <- surv_cutpoint(merged_data, time = "survival_time", 
                         event = "event_status",
                         variables = "Percentage")
summary(sur.cut)
#对截断值两侧的数据分布可视化。
plot(sur.cut, "Percentage", palette = "npg")

#使用surv_categorize函数，根据截断值转为二分类变量。
sur.cat <- surv_categorize(sur.cut, labels = c(0, 1))


data$Percentage <- as.factor(data$Percentage)
#0为较低得分，1为较高得分，指定较高得分为参考水平。参考水平需要加引号。
data$Percentage <- relevel(data$Percentage, ref = "1")


merged_data$mayoscore5_cut <- as.factor(merged_data$mayoscore5_cut)
#0为较低得分，1为较高得分，指定较高得分为参考水平。参考水平需要加引号。
data$mayoscore5_cut <- relevel(data$mayoscore5_cut, ref = "1")


merged_data$fusobacterium_cut <-
              ifelse(merged_data$Percentage >= 0.5347594, "High", "Low")

cox_model_1 <- coxph(Surv(survival_time, event_status) ~ fusobacterium_cut, data = merged_data)
summary(cox_model_1)

merged_data$fusobacterium_cut <- relevel(merged_data$fusobacterium_cut, ref = "0")

cox_model_2 <- coxph(Surv(time, censor) ~  fusobacterium_cut, data = data)
summary(cox_model_2)



library(survival)
data(lung)


library(survival)
expFile="surSigExp.txt"
View(expFile)
clinicalFile="clinicalNum.txt"

exp=read.table(expFile,sep="\t",header=T,check.names=F,row.names=1) #读取表达数据文件
cli=read.table(clinicalFile,sep="\t",header=T,check.names=F,row.names=1) #读取临床数据文件
samSample=intersect(row.names(exp),row.names(cli))
exp=exp[samSample,]
cli=cli[samSample,]
pFilter=0.05 #p值过滤条件

bioForest=function(coxFile=null,forestFile=null,forestCol=null){
  # 定义自定义函数bioForest，用于从Cox比例风险模型结果创建森林图。
  # 此函数接受参数，如coxFile（Cox模型结果文件）、forestFile（输出PDF文件）和forestCol（森林图框的颜色）。
  #读取输入文件
  rt <- read.table(coxFile,header=T,sep="\t",row.names=1,check.names=F)
  gene <- rownames(rt)
  #循环遍历表达数据中的每个基因。
  hr <- sprintf("%.3f",rt$"HR")
  hrLow <- sprintf("%.3f",rt$"HR.95L")
  hrHigh <- sprintf("%.3f",rt$"HR.95H")
  Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
  pVal <- ifelse(rt$pvalue<0.001, "<0.001", sprintf("%.3f", rt$pvalue))
  
  #输出图形
  pdf(file=forestFile, width = 7,height = 4.5)
  n <- nrow(rt)
  nRow <- n+1
  ylim <- c(1,nRow)
  layout(matrix(c(1,2),nc=2),width=c(3,2.5))
  
  #绘制森林图左边的临床信息
  xlim = c(0,3)
  par(mar=c(4,2.5,2,1))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
  text.cex=0.8
  text(0,n:1,gene,adj=0,cex=text.cex)
  text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
  text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1,)
  
  #绘制森林图
  par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
  xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
  arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=2.5)
  abline(v=1,col="black",lty=2,lwd=2)
  boxcolor = ifelse(as.numeric(hr) > 1, forestCol, forestCol)
  points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.3)
  axis(1)
  dev.off()
}
############绘制森林图函数############
#独立预后分析，输出表格
allOutTab=data.frame()
#初始化一个空的数据框 allOutTab，用于存储独立预后分析的结果
sigGenes=c("futime","fustat")
#一个字符向量 sigGenes，其中包含初始的两个列名 "futime" 和 "fustat"，这些列用于构建 Cox 比例风险模型。
for(i in colnames(exp[,3:ncol(exp)])){
  #循环遍历表达数据的列，其中包括具体的基因表达信息。
  #针对每个基因，构建 Cox 比例风险模型，并对模型的结果进行总结（summary）。
  rt=cbind(exp[,1:2],cli,gene=exp[,i])
  colnames(rt)[ncol(rt)]=i
  cox <- coxph(Surv(futime, fustat) ~ rt[,ncol(rt)], data = rt)
  coxSummary = summary(cox)
  outTab=data.frame()
  # p 值小于预设的过滤条件（pFilter），则继续在内部循环中，分别构建只包含一个基因表达变量的 Cox 比例风险模型。
  if(coxSummary$coefficients[,"Pr(>|z|)"]){
     for(j in colnames(rt[,3:ncol(rt)])){
       coxJ <- coxph(Surv(futime, fustat) ~ rt[,j], data = rt)
       coxSummaryJ = summary(coxJ)
       coxP=coxSummaryJ$coefficients[,"Pr(>|z|)"]
       if(coxSummaryJ$conf.int[,"upper .95"] != Inf){
         outTab=rbind(outTab,
                      cbind(id=j,
                            HR=coxSummaryJ$conf.int[,"exp(coef)"],
                            HR.95L=coxSummaryJ$conf.int[,"lower .95"],
                            HR.95H=coxSummaryJ$conf.int[,"upper .95"],
                            pvalue=coxSummaryJ$coefficients[,"Pr(>|z|)"])
         )
       }
     }
     write.table(outTab,file="uniCox.txt",sep="\t",row.names=F,quote=F)
     bioForest(coxFile="uniCox.txt",forestFile=paste0("uni.",i,".pdf"),forestCol="green") unlink("uniCox.txt")
     
     sigGenes=c(sigGenes,i)
     allOutTab=rbind(allOutTab,
                     cbind(id=i,
                           HR=coxSummary$conf.int[,"exp(coef)"],
                           HR.95L=coxSummary$conf.int[,"lower .95"],
                           HR.95H=coxSummary$conf.int[,"upper .95"],
                           pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
     )
}
}
#输出单因素独立预后分析结果
write.table(allOutTab,file="uniIndep.xls",sep="\t",row.names=F,quote=F)
#输出单因素独立预后基因表达量
indepSigExp=exp[,sigGenes]
indepSigExp=cbind(id=row.names(indepSigExp),indepSigExp)
write.table(indepSigExp,file="uniIndepSigExp.txt",sep="\t",row.names=F,quote=F)





####分组

# 计算 Fuso_abundance 的上四分位数 (upper quartile)
upper_quartile <- quantile(data$Fuso_abundance, 0.75)

# 基于上四分位数将样本分为 Fuso-high 和 Fuso-low
data$Fuso_group <- ifelse(data$Fuso_abundance >= upper_quartile, "Fuso-high", "Fuso-low")

# 查看结果
head(data)


####不同分组差异物种----
####数据准备

invade_seq<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/INVADE-seq/不同分组差异物种genus/umi_16s_MPR_group.csv",row.names=1) 

names(invade_seq)
otu_table <- invade_seq %>%
  group_by(samples) %>%
  summarise(across(1:513, ~ sum(.x, na.rm = TRUE))) 

otu_table<- t(otu_table)%>% as.data.frame()  
View(otu_table)
colnames(otu_table)<-otu_table[1,]
otu_table<-otu_table[-1,]
colnames(otu_table)
sample_table<-invade_seq[514:518] 
sample_data <- sample_table %>%
  distinct(samples, .keep_all = TRUE)
sample_data$samples

sample_data$samples
sample_order <- sample_data$samples

rownames(sample_data)<-sample_data$SampleID

colnames(sample_data)<-c("SampleID","cell_type_new","group_microbe2","Group","CRNCR_Response2")


# 重新排列 otu_table 的列顺序，使其与 sample_order 一致
otu_table <- otu_table[, match(sample_order, colnames(otu_table))]
colnames(otu_table)
rownames(otu_table)

Genus_taxonomy <-rownames(otu_table) %>% data.frame()
rownames(Genus_taxonomy)<-Genus_taxonomy$genus

View(Genus_taxonomy)
colnames(Genus_taxonomy)<- "genus"

otu_table[is.na(otu_table)] <- 0
View(otu_table)
otu_table[] <- lapply(otu_table, function(x) as.numeric(as.character(x)))

column_sums <- colSums(otu_table)

otu_percentage <- sweep(otu_table, 2, column_sums, FUN = "/") * 100
head(otu_percentage)
rownames(otu_percentage)


otu_table$taxonomy<- rownames(otu_table)
names(otu_table)
otu_table<- otu_table[,-15]

###########################
########################
library(readxl)
data <- read_excel("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/nature medicine 微生物生存分析/41591_2022_1698_MOESM3_ESM (1).xlsx")

metadata <- read_excel("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/nature medicine 微生物生存分析/41591_2022_1698_MOESM3_ESM (1).xlsx", sheet = "Metadata")
head(metadata)

LKT_PPM<-read_excel("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/nature medicine 微生物生存分析/41591_2022_1698_MOESM3_ESM (1).xlsx", sheet = "LKT_PPM")
head(LKT_PPM)
View(LKT_PPM)

LKT_featuretable<- read_excel("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/nature medicine 微生物生存分析/41591_2022_1698_MOESM3_ESM (1).xlsx", sheet = "LKT_featuretable")

head(LKT_featuretable)


LKT_featuretable <- LKT_featuretable[rownames(LKT_PPM), , drop = FALSE]

View(LKT_featuretable)

LKT_data <- cbind(LKT_featuretable, LKT_PPM)

LKT_data <- merge(LKT_featuretable, LKT_PPM, by = "row.names", all = TRUE)

LKT_featuretable$Row.names <- rownames(LKT_featuretable)
LKT_PPM$Row.names <- rownames(LKT_PPM)

LKT_data <- merge(LKT_featuretable, LKT_PPM, by = "Row.names", all = TRUE)

View(LKT_data)

Streptococcus_OTUs <- rownames(LKT_featuretable)[LKT_featuretable[, "Genus"] == "g__Streptococcus"]

unique(LKT_featuretable[, "Genus"])

Streptococcus_abundance <- LKT_data[LKT_data[, "Genus"] == "g__Streptococcus", ]

View(Streptococcus_abundance)

numeric_data <- Streptococcus_abundance[, sapply(Streptococcus_abundance, is.numeric)]

# 按行求和
col_sum <- colSums(numeric_data)

numeric_data2 <- LKT_data[, sapply(LKT_data, is.numeric)]

total_abundance <- colSums(numeric_data2)  # 每个样本的总丰度

length(total_abundance)

Streptococcus_percentage <- col_sum / total_abundance * 100

Streptococcus_percentage<- as.data.frame(Streptococcus_percentage)



Capnocytophaga_OTUs <- rownames(LKT_featuretable)[LKT_featuretable[, "Genus"] == "g__Capnocytophaga"]

unique(LKT_featuretable[, "Genus"])

Capnocytophaga_abundance <- LKT_data[LKT_data[, "Genus"] == "g__Capnocytophaga", ]

View(Capnocytophaga_abundance)

numeric_data <- Capnocytophaga_abundance[, sapply(Capnocytophaga_abundance, is.numeric)]

# 按行求和
col_sum <- colSums(numeric_data)

numeric_data2 <- LKT_data[, sapply(LKT_data, is.numeric)]

total_abundance <- colSums(numeric_data2)  # 每个样本的总丰度

length(total_abundance)

Capnocytophaga_percentage <- col_sum / total_abundance * 100

Capnocytophaga_percentage<- as.data.frame(Capnocytophaga_percentage)





library(ggplot2)
library(dplyr)
set.seed(123)

row_order <- rownames(Streptococcus_percentage)
metadata_sorted <- metadata[match(row_order, metadata$Sample ), ]
metadata_sorted$Sample  <- factor(metadata_sorted$Sample , levels = row_order)


data<- cbind(metadata_sorted,Streptococcus_percentage)

data<- cbind(metadata_sorted,Streptococcus_percentage,Capnocytophaga_percentage)
# Progression_Event
#PFS_days
# 计算 75% 分位点
threshold <- quantile(data$Streptococcus_percentage, 0.75,na.rm = TRUE)

# 添加分组信息
data$Streptococcus_group <- ifelse(data$Streptococcus_percentage >= threshold, "High", "Low")

unique(data$fusobacterium_group)

pp1<- ggplot(data, aes(x =Streptococcus_percentage, fill = Streptococcus_group)) +
  geom_density(alpha = 1) +  # 密度图，透明度为0.6
  geom_vline(xintercept = threshold, linetype = "dashed") +  # 添加分位点虚线
  scale_fill_manual(values = c("Low" = "skyblue", "High" = "red","NaN"="skyblue")) +  # 自定义填充颜色
  labs(
    x = "Streptococcus(%)",
    y = "Density",
    fill = "Streptococcus"
  ) +
  theme_minimal() +  # 设置主题
  theme(
    legend.position = "top",  # 调整图例位置
    axis.ticks.x = element_line(linewidth = 0.8), # 加粗 X 轴刻度线
    axis.ticks.y = element_line(linewidth = 0.8), # 加粗 Y 轴刻度线
    axis.title = element_text(size = 14),   # 调整坐标轴标题的大小
    axis.text = element_text(size = 12),     # 调整坐标轴文字的大小
    axis.line.x = element_line(linewidth = 0.8),  # 加粗 X 轴线条
    axis.line.y = element_line(linewidth = 0.8)  # 加粗 Y 轴线条
  )+
  # scale_x_continuous(expand = expansion(mult = c(0.1, 0.1))) + 
  # scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) + 
  coord_cartesian(xlim = c(-1, 20), ylim = c(0, 0.3)) 
pp1

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/nature medicine 微生物生存分析/Streptococcus_density.pdf",pp1,width = 5,height = 5)

library(survival)
library(survminer)
data$PFS_days

data$PFS_days <- as.numeric(as.character(data$PFS_days))
data$Progression_Event <- as.numeric(as.character(data$Progression_Event))
data$Age <- as.numeric(as.character(data$Age))
data$BMI <- as.numeric(as.character(data$BMI))

data_clean <- na.omit(data)

cox_model <- coxph(Surv(PFS_days, Progression_Event) ~ Streptococcus_percentage+Age+
                     Immunotherapy_Drug+Sex+
                     BMI+Pre_Treatment_Disease_Stage, data =data_clean)

cox_model <- coxph(Surv(PFS_days, Progression_Event) ~ Streptococcus_group+Age+
                     Sex+
                     BMI, data =data)

summary(cox_model)
cox.zph(cox_model)


f1_Streptococcus<-ggforest(cox_model,
                           data=data,
                           main = "Hazard ratio",        # 标题
                           cpositions = c(0.01, 0.1, 0.3), #左边三列的相对距离
                           fontsize = 1, # 字体大小
                           refLabel = "reference", #显示因子的参考水平
                           noDigits = 3)  #小数点位数

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/nature medicine 微生物生存分析/Streptococcus_ggforest.pdf",f1_Fusobacterium,width =13,height=6)

sur.cut <- surv_cutpoint(data, time = "PFS_days", 
                         event = "Progression_Event",
                         variables = "Streptococcus_percentage")
summary(sur.cut)
#对截断值两侧的数据分布可视化。
plot(sur.cut, "Streptococcus_percentage", palette = "npg")

#使用surv_categorize函数，根据截断值转为二分类变量。
sur.cat <- surv_categorize(sur.cut, labels = c(0, 1))


data$Streptococcus_cut <-ifelse(data$Streptococcus_percentage >=0.1390079, "High", "Low")
data$Streptococcus_percentage <- as.numeric(as.character(data$Streptococcus_percentage))

data$Streptococcus_cut <- ifelse(data$Streptococcus_percentage >= 0.1390079, "High", "Low")


fit <- survfit(Surv(PFS_days, Progression_Event) ~Streptococcus_cut,  # 创建生存对象 
               data = data)

pdf("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/nature medicine 微生物生存分析/Streptococcus_curvep.pdf", width = 6, height =7)
ggsurvplot(fit, # 创建的拟合对象
           data = data,  # 指定变量数据来源
           conf.int = TRUE, # 显示置信区间
           pval = TRUE, # 添加P值
           surv.median.line = "hv",  # 添加中位生存时间线
           risk.table = TRUE, # 添加风险表
           ylab="PFS probability(%)",
           xlab = "PFS_days(d)", # 指定x轴标签
           legend = c(0.8,0.75), # 指定图例位置
           legend.title = "", # 设置图例标题，这里设置不显示标题，用空格替代
           palette = c("skyblue", "red"),
           legend.labs = c("Low", "High"), # 指定图例分组标签
           break.x.by = 500,
           title = "Streptococcus" )  # 设置x轴刻度间距
dev.off()



###################################################

threshold <- quantile(data$Capnocytophaga_percentage, 0.75,na.rm = TRUE)

# 添加分组信息
data$Capnocytophaga_group <- ifelse(data$Capnocytophaga_percentage >= threshold, "High", "Low")

unique(data$fusobacterium_group)

pp1<- ggplot(data, aes(x =Capnocytophaga_percentage, fill = Capnocytophaga_group)) +
  geom_density(alpha = 1) +  # 密度图，透明度为0.6
  geom_vline(xintercept = threshold, linetype = "dashed") +  # 添加分位点虚线
  scale_fill_manual(values = c("Low" = "skyblue", "High" = "red","NaN"="skyblue")) +  # 自定义填充颜色
  labs(
    x = "Capnocytophaga(%)",
    y = "Density",
    fill = "Capnocytophaga"
  ) +
  theme_minimal() +  # 设置主题
  theme(
    legend.position = "top",  # 调整图例位置
    axis.ticks.x = element_line(linewidth = 0.8), # 加粗 X 轴刻度线
    axis.ticks.y = element_line(linewidth = 0.8), # 加粗 Y 轴刻度线
    axis.title = element_text(size = 14),   # 调整坐标轴标题的大小
    axis.text = element_text(size = 12),     # 调整坐标轴文字的大小
    axis.line.x = element_line(linewidth = 0.8),  # 加粗 X 轴线条
    axis.line.y = element_line(linewidth = 0.8)  # 加粗 Y 轴线条
  )+
  # scale_x_continuous(expand = expansion(mult = c(0.1, 0.1))) + 
  # scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) + 
  coord_cartesian(xlim = c(-1, 20), ylim = c(0, 0.3)) 
pp1

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/nature medicine 微生物生存分析/Capnocytophaga_density.pdf",pp1,width = 5,height = 5)

library(survival)
library(survminer)
data$PFS_days

data$PFS_days <- as.numeric(as.character(data$PFS_days))
data$Progression_Event <- as.numeric(as.character(data$Progression_Event))
data$Age <- as.numeric(as.character(data$Age))
data$BMI <- as.numeric(as.character(data$BMI))

data_clean <- na.omit(data)

cox_model <- coxph(Surv(PFS_days, Progression_Event) ~ Capnocytophaga_percentage+Age+
                     Immunotherapy_Drug+Sex+
                     BMI+Pre_Treatment_Disease_Stage, data =data_clean)

cox_model <- coxph(Surv(PFS_days, Progression_Event) ~ Capnocytophaga_group+Age+
                     Sex+
                     BMI, data =data)

summary(cox_model)
cox.zph(cox_model)


f1_Capnocytophaga<-ggforest(cox_model,
                           data=data,
                           main = "Hazard ratio",        # 标题
                           cpositions = c(0.01, 0.1, 0.3), #左边三列的相对距离
                           fontsize = 1, # 字体大小
                           refLabel = "reference", #显示因子的参考水平
                           noDigits = 3)  #小数点位数

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/nature medicine 微生物生存分析/Capnocytophaga_ggforest.pdf",f1_Fusobacterium,width =13,height=6)

sur.cut <- surv_cutpoint(data, time = "PFS_days", 
                         event = "Progression_Event",
                         variables = "Capnocytophaga_percentage")
summary(sur.cut)
#对截断值两侧的数据分布可视化。
plot(sur.cut, "Capnocytophaga_percentage", palette = "npg")

#使用surv_categorize函数，根据截断值转为二分类变量。
sur.cat <- surv_categorize(sur.cut, labels = c(0, 1))
data$Capnocytophaga_group

data$Capnocytophaga_cut <-ifelse(data$Capnocytophaga_percentage >=0.1390079, "High", "Low")
data$Capnocytophaga_percentage <- as.numeric(as.character(data$Capnocytophaga_percentage))

data$Capnocytophaga_cut <- ifelse(data$Capnocytophaga_percentage >= 0.1390079, "High", "Low")


fit <- survfit(Surv(PFS_days, Progression_Event) ~Capnocytophaga_group,  # 创建生存对象 
               data = data)

pdf("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/nature medicine 微生物生存分析/Capnocytophaga_curvep.pdf", width = 6, height =7)
ggsurvplot(fit, # 创建的拟合对象
           data = data,  # 指定变量数据来源
           conf.int = TRUE, # 显示置信区间
           pval = TRUE, # 添加P值
           surv.median.line = "hv",  # 添加中位生存时间线
           risk.table = TRUE, # 添加风险表
           ylab="PFS probability(%)",
           xlab = "PFS_days(d)", # 指定x轴标签
           legend = c(0.8,0.75), # 指定图例位置
           legend.title = "", # 设置图例标题，这里设置不显示标题，用空格替代
           palette = c("skyblue", "red"),
           legend.labs = c("Low", "High"), # 指定图例分组标签
           break.x.by = 500,
           title = "Capnocytophaga" )  # 设置x轴刻度间距
dev.off()



