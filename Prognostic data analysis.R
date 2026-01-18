###预后分析数据收集

####https://github.com/twbattaglia/MicrobeDS
####This package contains datasets provided by large-scale microbiome studies. Each dataset is formatted for use with phyloseq.
####(https://joey711.github.io/phyloseq/).

remotes::install_github("twbattaglia/MicrobeDS")
library(phyloseq)
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

otu_table_data <- phyloseq::otu_table(Hartwig)
otu_table_data <- phyloseq::otu_table(HMPv13)


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
# 查看 TCGA 中的癌症类型（可用的项目）
tcga_projects <- TCGAbiolinks::getGDCprojects()
tcga_cancer_types <- tcga_projects[grep("TCGA", tcga_projects$project_id), c("project_id", "name")]
print(tcga_cancer_types)
print(tcga_cancer_types)
#project_id  name

#2   TCGA-KIRC                                Kidney Renal Clear Cell Carcinoma
#4   TCGA-BRCA                                        Breast Invasive Carcinoma
#55  TCGA-DLBC                  Lymphoid Neoplasm Diffuse Large B-cell Lymphoma
#56  TCGA-PCPG                               Pheochromocytoma and Paraganglioma
#57  TCGA-READ                                            Rectum Adenocarcinoma
#58  TCGA-CHOL                                               Cholangiocarcinoma
#59  TCGA-ESCA                                             Esophageal Carcinoma
#60  TCGA-THYM                                                          Thymoma
#61  TCGA-MESO                                                     Mesothelioma
#62   TCGA-UVM                                                   Uveal Melanoma
#63  TCGA-KICH                                               Kidney Chromophobe
#64  TCGA-SARC                                                          Sarcoma
#65  TCGA-PAAD                                        Pancreatic Adenocarcinoma
#66   TCGA-UCS                                           Uterine Carcinosarcoma
#67  TCGA-LIHC                                   Liver Hepatocellular Carcinoma
#68  TCGA-KIRP                            Kidney Renal Papillary Cell Carcinoma
#69  TCGA-LAML                                           Acute Myeloid Leukemia
#70   TCGA-ACC                                         Adrenocortical Carcinoma
#71  TCGA-TGCT                                      Testicular Germ Cell Tumors
#72  TCGA-SKCM                                          Skin Cutaneous Melanoma
#73  TCGA-CESC Cervical Squamous Cell Carcinoma and Endocervical Adenocarcinoma
#74  TCGA-COAD                                             Colon Adenocarcinoma
#75  TCGA-BLCA                                     Bladder Urothelial Carcinoma
#76  TCGA-UCEC                             Uterine Corpus Endometrial Carcinoma
#77  TCGA-PRAD                                          Prostate Adenocarcinoma
#78  TCGA-STAD                                           Stomach Adenocarcinoma
#80   TCGA-LGG                                         Brain Lower Grade Glioma
#81  TCGA-HNSC                            Head and Neck Squamous Cell Carcinoma
#82  TCGA-LUAD                                              Lung Adenocarcinoma
#83  TCGA-THCA                                                Thyroid Carcinoma
#84   TCGA-GBM                                          Glioblastoma Multiforme
#85  TCGA-LUSC                                     Lung Squamous Cell Carcinoma
#86    TCGA-OV                                Ovarian Serous Cystadenocarcinoma


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

###########################################

project_ids <- c("TCGA-ACC",
                 "TCGA-BLCA", "TCGA-BRCA", "TCGA-CESC", "TCGA-CHOL", "TCGA-COAD", "TCGA-DLBC", 
                 "TCGA-ESCA", "TCGA-GBM", "TCGA-HNSC", "TCGA-KICH", "TCGA-KIRC", "TCGA-KIRP", 
                 "TCGA-LAML", "TCGA-LGG", "TCGA-LIHC", "TCGA-LUAD", "TCGA-LUSC", "TCGA-MESO", 
                 "TCGA-OV", "TCGA-PAAD", "TCGA-PCPG", "TCGA-PRAD", "TCGA-READ", "TCGA-SARC", 
                 "TCGA-SKCM", "TCGA-STAD", "TCGA-TGCT", "TCGA-THCA", "TCGA-THYM", "TCGA-UCEC", 
                 "TCGA-UCS", "TCGA-UVM"
)

####fusobacterium TCGA-BLCA ----

merged_data2<- merged_data %>% filter(project %in% "TCGA-BLCA")
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
  ifelse(merged_data2$Fusobacterium_percentage >= 0.004534258, "High", "Low")


fit <- survfit(Surv(survival_time, event_status) ~Fusobacterium_cut,  # 创建生存对象 
               data = merged_data2)

pdf("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Apan-cancer_analysis_of_the_microbiome/TCGA-BLCA/Fusobacterium_BLCA_curvep_cut.pdf", width = 6, height =7)
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
           title = "Fusobacterium BLCA" )  # 设置x轴刻度间距
dev.off()


merged_data2$Fusobacterium_cut <- factor(
  merged_data2$Fusobacterium_cut,
  levels = c("Low", "High") # 包含所有可能值
)

merged_data2$gender.x <- factor(merged_data2$gender.x)
sum(is.na(merged_data2$gender.x))

names(merged_data2)
unique(merged_data2$years_smoked)
merged_data2$years_smoked_history<-ifelse(merged_data2$years_smoked <= 20, "Low", "High")
merged_data2$portion_weight2 <- log(merged_data2$portion_weight + 1)
merged_data2$years_smoked2 <- log(merged_data2$years_smoked + 1)
merged_data2$years_smoked2[is.na(merged_data2$years_smoked2)] <- mean(merged_data2$years_smoked2, na.rm = TRUE)
merged_data2$Fusobacterium_cut <- relevel(merged_data2$Fusobacterium_cut, ref = "Low")

merged_data2$age_at_index
merged_data2$alcohol_history


cox_model <- coxph(Surv(survival_time, event_status) ~ Fusobacterium_cut+age_at_index+
                   prior_malignancy+
                     portion_weight2+gender.x, data = merged_data2)

summary(cox_model)
f1_Fusobacterium<-ggforest(cox_model,
                           data=merged_data2,
                           main = "Hazard ratio",        # 标题
                           cpositions = c(0.01, 0.1, 0.3), #左边三列的相对距离
                           fontsize = 1, # 字体大小
                           refLabel = "reference", #显示因子的参考水平
                           noDigits = 3)  #小数点位数

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Apan-cancer_analysis_of_the_microbiome/TCGA-BLCA/fusobacterium_BLCA_ggforest_cut_多因素.pdf",f1_Fusobacterium,width =13,height=6)



#####Capnocytophaga TCGA-BLCA----

merged_data2<- merged_data %>% filter(project %in% "TCGA-BLCA")
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
  ifelse(merged_data2$Capnocytophaga_percentage >= 0.001787426, "High", "Low")


fit <- survfit(Surv(survival_time, event_status) ~Capnocytophaga_cut,  # 创建生存对象 
               data = merged_data2)

pdf("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Apan-cancer_analysis_of_the_microbiome/TCGA-BLCA/Capnocytophaga_BLCA_curvep_cut.pdf", width = 6, height =7)
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
           break.x.by = 500,
           xlim = c(0, 2000),
           title = "Capnocytophaga BLCA" )  # 设置x轴刻度间距
dev.off()


merged_data2$Capnocytophaga_cut <- factor(
  merged_data2$Capnocytophaga_cut,
  levels = c("Low", "High") # 包含所有可能值
)

merged_data2$gender.x <- factor(merged_data2$gender.x)
sum(is.na(merged_data2$gender.x))

names(merged_data2)
unique(merged_data2$years_smoked)
merged_data2$years_smoked_history<-ifelse(merged_data2$years_smoked <= 20, "Low", "High")
merged_data2$portion_weight2 <- log(merged_data2$portion_weight + 1)
merged_data2$years_smoked2 <- log(merged_data2$years_smoked + 1)
merged_data2$years_smoked2[is.na(merged_data2$years_smoked2)] <- mean(merged_data2$years_smoked2, na.rm = TRUE)
merged_data2$Capnocytophaga_cut <- relevel(merged_data2$Capnocytophaga_cut, ref = "Low")

merged_data2$age_at_index
merged_data2$alcohol_history


cox_model <- coxph(Surv(survival_time, event_status) ~ Capnocytophaga_cut+age_at_index+
                     prior_malignancy+
                     portion_weight2+gender.x, data = merged_data2)

summary(cox_model)
f1_Capnocytophaga<-ggforest(cox_model,
                           data=merged_data2,
                           main = "Hazard ratio",        # 标题
                           cpositions = c(0.01, 0.1, 0.3), #左边三列的相对距离
                           fontsize = 1, # 字体大小
                           refLabel = "reference", #显示因子的参考水平
                           noDigits = 3)  #小数点位数

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Apan-cancer_analysis_of_the_microbiome/TCGA-BLCA/Capnocytophaga_BLCA_ggforest_cut_多因素.pdf",f1_Capnocytophaga,width =13,height=6)


#####Streptococcus TCGA-BLCA----

merged_data2<- merged_data %>% filter(project %in% "TCGA-BLCA")
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
  ifelse(merged_data2$Streptococcus_percentage >= 24.31593, "High", "Low")


fit <- survfit(Surv(survival_time, event_status) ~Streptococcus_cut,  # 创建生存对象 
               data = merged_data2)

pdf("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Apan-cancer_analysis_of_the_microbiome/TCGA-BLCA/Streptococcus_BLCA_curvep_cut.pdf", width = 6, height =7)
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
           break.x.by = 500,
           xlim = c(0, 2000),
           title = "Streptococcus BLCA" )  # 设置x轴刻度间距
dev.off()


merged_data2$Streptococcus_cut <- factor(
  merged_data2$Streptococcus_cut,
  levels = c("Low", "High") # 包含所有可能值
)

merged_data2$gender.x <- factor(merged_data2$gender.x)
sum(is.na(merged_data2$gender.x))

names(merged_data2)
unique(merged_data2$years_smoked)
merged_data2$years_smoked_history<-ifelse(merged_data2$years_smoked <= 20, "Low", "High")
merged_data2$portion_weight2 <- log(merged_data2$portion_weight + 1)
merged_data2$years_smoked2 <- log(merged_data2$years_smoked + 1)
merged_data2$years_smoked2[is.na(merged_data2$years_smoked2)] <- mean(merged_data2$years_smoked2, na.rm = TRUE)
merged_data2$Streptococcus_cut <- relevel(merged_data2$Streptococcus_cut, ref = "Low")

merged_data2$age_at_index
merged_data2$alcohol_history


cox_model <- coxph(Surv(survival_time, event_status) ~ Streptococcus_cut+age_at_index+
                     prior_malignancy+
                     portion_weight2+gender.x, data = merged_data2)

summary(cox_model)
f1_Streptococcus<-ggforest(cox_model,
                            data=merged_data2,
                            main = "Hazard ratio",        # 标题
                            cpositions = c(0.01, 0.1, 0.3), #左边三列的相对距离
                            fontsize = 1, # 字体大小
                            refLabel = "reference", #显示因子的参考水平
                            noDigits = 3)  #小数点位数

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Apan-cancer_analysis_of_the_microbiome/TCGA-BLCA/Streptococcus_BLCA_ggforest_cut_多因素.pdf",f1_Streptococcus,width =13,height=6)


project_ids <- c("TCGA-ACC",
                 "TCGA-BLCA", "TCGA-BRCA", "TCGA-CESC", "TCGA-CHOL", "TCGA-COAD", "TCGA-DLBC", 
                 "TCGA-ESCA", "TCGA-GBM", "TCGA-HNSC", "TCGA-KICH", "TCGA-KIRC", "TCGA-KIRP", 
                 "TCGA-LAML", "TCGA-LGG", "TCGA-LIHC", "TCGA-LUAD", "TCGA-LUSC", "TCGA-MESO", 
                 "TCGA-OV", "TCGA-PAAD", "TCGA-PCPG", "TCGA-PRAD", "TCGA-READ", "TCGA-SARC", 
                 "TCGA-SKCM", "TCGA-STAD", "TCGA-TGCT", "TCGA-THCA", "TCGA-THYM", "TCGA-UCEC", 
                 "TCGA-UCS", "TCGA-UVM"
)

####fusobacterium TCGA-LIHC 肝细胞癌----

merged_data2<- merged_data %>% filter(project %in% "TCGA-LIHC")
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
  ifelse(merged_data2$Fusobacterium_percentage >= 0.002787488, "High", "Low")


fit <- survfit(Surv(survival_time, event_status) ~Fusobacterium_cut,  # 创建生存对象 
               data = merged_data2)

pdf("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Apan-cancer_analysis_of_the_microbiome/TCGA-LIHC肝细胞癌/Fusobacterium_TCGA-LIHC_curvep_cut.pdf", width = 6, height =7)
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
           break.x.by = 200,
           xlim = c(0, 1000),
           title = "Fusobacterium LIHC" )  # 设置x轴刻度间距
dev.off()

merged_data2$Fusobacterium_cut <- factor(
  merged_data2$Fusobacterium_cut,
  levels = c("Low", "High") # 包含所有可能值
)

merged_data2$Fusobacterium_cut <- factor(
  merged_data2$Fusobacterium_cut,
  levels = c("High", "Low") # 包含所有可能值
)


merged_data2$Fusobacterium_cut <- relevel(merged_data2$Fusobacterium_cut, ref = "High")

library(survival)

# 创建Cox模型
cox_model <- coxph(Surv(survival_time, event_status) ~ Fusobacterium_cut, data = merged_data2)

# 查看模型摘要
summary(cox_model)



merged_data2$gender.x <- factor(merged_data2$gender.x)
sum(is.na(merged_data2$gender.x))

names(merged_data2)
unique(merged_data2$years_smoked)
merged_data2$years_smoked_history<-ifelse(merged_data2$years_smoked <= 20, "Low", "High")
merged_data2$portion_weight2 <- log(merged_data2$portion_weight + 1)
merged_data2$years_smoked2 <- log(merged_data2$years_smoked + 1)
merged_data2$years_smoked2[is.na(merged_data2$years_smoked2)] <- mean(merged_data2$years_smoked2, na.rm = TRUE)
merged_data2$Fusobacterium_cut <- relevel(merged_data2$Fusobacterium_cut, ref = "High")

merged_data2$age_at_index
merged_data2$alcohol_history


cox_model <- coxph(Surv(survival_time, event_status) ~ Fusobacterium_cut+age_at_index+
                     prior_malignancy+
                     portion_weight2+gender.x, data = merged_data2)

summary(cox_model)
f1_Fusobacterium<-ggforest(cox_model,
                           data=merged_data2,
                           main = "Hazard ratio",        # 标题
                           cpositions = c(0.01, 0.1, 0.3), #左边三列的相对距离
                           fontsize = 1, # 字体大小
                           refLabel = "reference", #显示因子的参考水平
                           noDigits = 3)  #小数点位数

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Apan-cancer_analysis_of_the_microbiome/TCGA-LIHC肝细胞癌/fusobacterium_LIHC_ggforest_cut_多因素.pdf",f1_Fusobacterium,width =13,height=6)



#####Capnocytophaga TCGA-LIHC----

merged_data2<- merged_data %>% filter(project %in% "TCGA-LIHC")
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
  ifelse(merged_data2$Capnocytophaga_percentage >= 0.002787488, "High", "Low")


fit <- survfit(Surv(survival_time, event_status) ~Capnocytophaga_cut,  # 创建生存对象 
               data = merged_data2)

pdf("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Apan-cancer_analysis_of_the_microbiome/TCGA-LIHC肝细胞癌/Capnocytophaga_TCGA-LIHC_curvep_cut.pdf", width = 6, height =7)
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
           break.x.by = 200,
           xlim = c(0, 1000),
           title = "Capnocytophaga LIHC" )  # 设置x轴刻度间距
dev.off()

merged_data2$Capnocytophaga_cut <- factor(
  merged_data2$Capnocytophaga_cut,
  levels = c("Low", "High") # 包含所有可能值
)

merged_data2$Capnocytophaga_cut <- factor(
  merged_data2$Capnocytophaga_cut,
  levels = c("High", "Low") # 包含所有可能值
)


merged_data2$Capnocytophaga_cut <- relevel(merged_data2$Capnocytophaga_cut, ref = "High")

library(survival)

# 创建Cox模型
cox_model <- coxph(Surv(survival_time, event_status) ~ Capnocytophaga_cut, data = merged_data2)

# 查看模型摘要
summary(cox_model)



merged_data2$gender.x <- factor(merged_data2$gender.x)
sum(is.na(merged_data2$gender.x))

names(merged_data2)
unique(merged_data2$years_smoked)
merged_data2$years_smoked_history<-ifelse(merged_data2$years_smoked <= 20, "Low", "High")
merged_data2$portion_weight2 <- log(merged_data2$portion_weight + 1)
merged_data2$years_smoked2 <- log(merged_data2$years_smoked + 1)
merged_data2$years_smoked2[is.na(merged_data2$years_smoked2)] <- mean(merged_data2$years_smoked2, na.rm = TRUE)
merged_data2$Capnocytophaga_cut <- relevel(merged_data2$Capnocytophaga_cut, ref = "High")

merged_data2$age_at_index
merged_data2$alcohol_history


cox_model <- coxph(Surv(survival_time, event_status) ~ Capnocytophaga_cut+age_at_index+
                     prior_malignancy+
                     portion_weight2+gender.x, data = merged_data2)

summary(cox_model)
f1_Capnocytophaga<-ggforest(cox_model,
                           data=merged_data2,
                           main = "Hazard ratio",        # 标题
                           cpositions = c(0.01, 0.1, 0.3), #左边三列的相对距离
                           fontsize = 1, # 字体大小
                           refLabel = "reference", #显示因子的参考水平
                           noDigits = 3)  #小数点位数

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Apan-cancer_analysis_of_the_microbiome/TCGA-LIHC肝细胞癌/Capnocytophaga_LIHC_ggforest_cut_多因素.pdf",f1_Capnocytophaga,width =13,height=6)



#####Streptococcus TCGA-LIHC----

merged_data2<- merged_data %>% filter(project %in% "TCGA-LIHC")
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
  ifelse(merged_data2$Streptococcus_percentage >= 24.89075, "High", "Low")


fit <- survfit(Surv(survival_time, event_status) ~Streptococcus_cut,  # 创建生存对象 
               data = merged_data2)

pdf("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Apan-cancer_analysis_of_the_microbiome/TCGA-LIHC肝细胞癌/Streptococcus_LIHC_curvep_cut.pdf", width = 6, height =7)
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
           break.x.by = 200,
           xlim = c(0, 1000),
           title = "Streptococcus LIHC" )  # 设置x轴刻度间距
dev.off()


merged_data2$Streptococcus_cut <- factor(
  merged_data2$Streptococcus_cut,
  levels = c("High", "Low") # 包含所有可能值
)

library(survival)

# 创建Cox模型
cox_model <- coxph(Surv(survival_time, event_status) ~ Streptococcus_cut, data = merged_data2)

# 查看模型摘要
summary(cox_model)


merged_data2$gender.x <- factor(merged_data2$gender.x)
sum(is.na(merged_data2$gender.x))

names(merged_data2)
unique(merged_data2$years_smoked)
merged_data2$years_smoked_history<-ifelse(merged_data2$years_smoked <= 20, "Low", "High")
merged_data2$portion_weight2 <- log(merged_data2$portion_weight + 1)
merged_data2$years_smoked2 <- log(merged_data2$years_smoked + 1)
merged_data2$years_smoked2[is.na(merged_data2$years_smoked2)] <- mean(merged_data2$years_smoked2, na.rm = TRUE)
merged_data2$Streptococcus_cut <- relevel(merged_data2$Streptococcus_cut, ref = "High")

merged_data2$age_at_index
merged_data2$alcohol_history


cox_model <- coxph(Surv(survival_time, event_status) ~ Streptococcus_cut+age_at_index+
                     prior_malignancy+
                     portion_weight2+gender.x, data = merged_data2)

summary(cox_model)
f1_Streptococcus<-ggforest(cox_model,
                           data=merged_data2,
                           main = "Hazard ratio",        # 标题
                           cpositions = c(0.01, 0.1, 0.3), #左边三列的相对距离
                           fontsize = 1, # 字体大小
                           refLabel = "reference", #显示因子的参考水平
                           noDigits = 3)  #小数点位数

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Apan-cancer_analysis_of_the_microbiome/TCGA-LIHC/Streptococcus_LIHC_ggforest_cut_多因素.pdf",f1_Streptococcus,width =13,height=6)


project_ids <- c("TCGA-ACC",
                 "TCGA-BLCA", "TCGA-BRCA", "TCGA-CESC", "TCGA-CHOL", "TCGA-COAD", "TCGA-DLBC", 
                 "TCGA-ESCA", "TCGA-GBM", "TCGA-HNSC", "TCGA-KICH", "TCGA-KIRC", "TCGA-KIRP", 
                 "TCGA-LAML", "TCGA-LGG", "TCGA-LIHC", "TCGA-LUAD", "TCGA-LUSC", "TCGA-MESO", 
                 "TCGA-OV", "TCGA-PAAD", "TCGA-PCPG", "TCGA-PRAD", "TCGA-READ", "TCGA-SARC", 
                 "TCGA-SKCM", "TCGA-STAD", "TCGA-TGCT", "TCGA-THCA", "TCGA-THYM", "TCGA-UCEC", 
                 "TCGA-UCS", "TCGA-UVM"
)

####fusobacterium TCGA-LUAD ----

merged_data2<- merged_data %>% filter(project %in% "TCGA-LUAD")
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
  ifelse(merged_data2$Fusobacterium_percentage >= 0.003573513, "High", "Low")


fit <- survfit(Surv(survival_time, event_status) ~Fusobacterium_cut,  # 创建生存对象 
               data = merged_data2)

pdf("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Apan-cancer_analysis_of_the_microbiome/TCGA-LUAD/Fusobacterium_TCGA-LUAD_curvep_cut.pdf", width = 6, height =7)
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
           xlim = c(0, 3500),
           title = "Fusobacterium LUAD" )  # 设置x轴刻度间距
dev.off()


merged_data2$Fusobacterium_cut <- factor(
  merged_data2$Fusobacterium_cut,
  levels = c("Low", "High") # 包含所有可能值
)

merged_data2$gender.x <- factor(merged_data2$gender.x)
sum(is.na(merged_data2$gender.x))

names(merged_data2)
unique(merged_data2$years_smoked)
merged_data2$years_smoked_history<-ifelse(merged_data2$years_smoked <= 20, "Low", "High")
merged_data2$portion_weight2 <- log(merged_data2$portion_weight + 1)
merged_data2$years_smoked2 <- log(merged_data2$years_smoked + 1)
merged_data2$years_smoked2[is.na(merged_data2$years_smoked2)] <- mean(merged_data2$years_smoked2, na.rm = TRUE)
merged_data2$Fusobacterium_cut <- relevel(merged_data2$Fusobacterium_cut, ref = "Low")

merged_data2$age_at_index
merged_data2$alcohol_history


cox_model <- coxph(Surv(survival_time, event_status) ~ Fusobacterium_cut+age_at_index+
                     prior_malignancy+
                     portion_weight2+gender.x, data = merged_data2)

summary(cox_model)
f1_Fusobacterium<-ggforest(cox_model,
                           data=merged_data2,
                           main = "Hazard ratio",        # 标题
                           cpositions = c(0.01, 0.1, 0.3), #左边三列的相对距离
                           fontsize = 1, # 字体大小
                           refLabel = "reference", #显示因子的参考水平
                           noDigits = 3)  #小数点位数

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Apan-cancer_analysis_of_the_microbiome/TCGA-LUAD/fusobacterium_LUAD_ggforest_cut_多因素.pdf",f1_Fusobacterium,width =13,height=6)



#####Capnocytophaga TCGA-LUAD----

merged_data2<- merged_data %>% filter(project %in% "TCGA-LUAD")
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
  ifelse(merged_data2$Capnocytophaga_percentage >= 0.003779839, "High", "Low")


fit <- survfit(Surv(survival_time, event_status) ~Capnocytophaga_cut,  # 创建生存对象 
               data = merged_data2)

pdf("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Apan-cancer_analysis_of_the_microbiome/TCGA-LUAD/Capnocytophaga_LUAD_curvep_cut.pdf", width = 6, height =7)
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
           break.x.by = 500,
           xlim = c(0, 2000),
           title = "Capnocytophaga LUAD" )  # 设置x轴刻度间距
dev.off()


merged_data2$Capnocytophaga_cut <- factor(
  merged_data2$Capnocytophaga_cut,
  levels = c("Low", "High") # 包含所有可能值
)

merged_data2$gender.x <- factor(merged_data2$gender.x)
sum(is.na(merged_data2$gender.x))

names(merged_data2)
unique(merged_data2$years_smoked)
merged_data2$years_smoked_history<-ifelse(merged_data2$years_smoked <= 20, "Low", "High")
merged_data2$portion_weight2 <- log(merged_data2$portion_weight + 1)
merged_data2$years_smoked2 <- log(merged_data2$years_smoked + 1)
merged_data2$years_smoked2[is.na(merged_data2$years_smoked2)] <- mean(merged_data2$years_smoked2, na.rm = TRUE)
merged_data2$Capnocytophaga_cut <- relevel(merged_data2$Capnocytophaga_cut, ref = "Low")

merged_data2$age_at_index
merged_data2$alcohol_history


cox_model <- coxph(Surv(survival_time, event_status) ~ Capnocytophaga_cut+age_at_index+
                     prior_malignancy+
                     portion_weight2+gender.x, data = merged_data2)

summary(cox_model)
f1_Capnocytophaga<-ggforest(cox_model,
                            data=merged_data2,
                            main = "Hazard ratio",        # 标题
                            cpositions = c(0.01, 0.1, 0.3), #左边三列的相对距离
                            fontsize = 1, # 字体大小
                            refLabel = "reference", #显示因子的参考水平
                            noDigits = 3)  #小数点位数

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Apan-cancer_analysis_of_the_microbiome/TCGA-LUAD/Capnocytophaga_LUAD_ggforest_cut_多因素.pdf",f1_Capnocytophaga,width =13,height=6)


#####Streptococcus TCGA-LUAD----

merged_data2<- merged_data %>% filter(project %in% "TCGA-LUAD")
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
  ifelse(merged_data2$Streptococcus_percentage >= 17.31877, "High", "Low")


fit <- survfit(Surv(survival_time, event_status) ~Streptococcus_cut,  # 创建生存对象 
               data = merged_data2)

pdf("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Apan-cancer_analysis_of_the_microbiome/TCGA-LUAD/Streptococcus_LUAD_curvep_cut.pdf", width = 6, height =7)
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
           break.x.by = 200,
           xlim = c(0, 1000),
           title = "Streptococcus LUAD" )  # 设置x轴刻度间距
dev.off()


merged_data2$Streptococcus_cut <- factor(
  merged_data2$Streptococcus_cut,
  levels = c("Low", "High") # 包含所有可能值
)

merged_data2$gender.x <- factor(merged_data2$gender.x)
sum(is.na(merged_data2$gender.x))

names(merged_data2)
unique(merged_data2$years_smoked)
merged_data2$years_smoked_history<-ifelse(merged_data2$years_smoked <= 20, "Low", "High")
merged_data2$portion_weight2 <- log(merged_data2$portion_weight + 1)
merged_data2$years_smoked2 <- log(merged_data2$years_smoked + 1)
merged_data2$years_smoked2[is.na(merged_data2$years_smoked2)] <- mean(merged_data2$years_smoked2, na.rm = TRUE)
merged_data2$Streptococcus_cut <- relevel(merged_data2$Streptococcus_cut, ref = "Low")

merged_data2$age_at_index
merged_data2$alcohol_history


cox_model <- coxph(Surv(survival_time, event_status) ~ Streptococcus_cut+age_at_index+
                     prior_malignancy+
                     portion_weight2+gender.x, data = merged_data2)

summary(cox_model)
f1_Streptococcus<-ggforest(cox_model,
                           data=merged_data2,
                           main = "Hazard ratio",        # 标题
                           cpositions = c(0.01, 0.1, 0.3), #左边三列的相对距离
                           fontsize = 1, # 字体大小
                           refLabel = "reference", #显示因子的参考水平
                           noDigits = 3)  #小数点位数

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Apan-cancer_analysis_of_the_microbiome/TCGA-LUAD/Streptococcus_LUAD_ggforest_cut_多因素.pdf",f1_Streptococcus,width =13,height=6)


project_ids <- c("TCGA-ACC",
                 "TCGA-BLCA", "TCGA-BRCA", "TCGA-CESC", "TCGA-CHOL", "TCGA-COAD", "TCGA-DLBC", 
                 "TCGA-ESCA", "TCGA-GBM", "TCGA-HNSC", "TCGA-KICH", "TCGA-KIRC", "TCGA-KIRP", 
                 "TCGA-LAML", "TCGA-LGG", "TCGA-LIHC", "TCGA-LUAD", "TCGA-LUSC", "TCGA-MESO", 
                 "TCGA-OV", "TCGA-PAAD", "TCGA-PCPG", "TCGA-PRAD", "TCGA-READ", "TCGA-SARC", 
                 "TCGA-SKCM", "TCGA-STAD", "TCGA-TGCT", "TCGA-THCA", "TCGA-THYM", "TCGA-UCEC", 
                 "TCGA-UCS", "TCGA-UVM"
)

####fusobacterium TCGA-STAD 胃腺癌----

merged_data2<- merged_data %>% filter(project %in% "TCGA-STAD")
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
  ifelse(merged_data2$Fusobacterium_percentage >= 0.001874307, "High", "Low")


fit <- survfit(Surv(survival_time, event_status) ~Fusobacterium_cut,  # 创建生存对象 
               data = merged_data2)

pdf("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Apan-cancer_analysis_of_the_microbiome/TCGA-STAD胃腺癌/Fusobacterium_TCGA-STAD_curvep_cut.pdf", width = 6, height =7)
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
           break.x.by = 500,
           xlim = c(0, 3500),
           title = "Fusobacterium STAD" )  # 设置x轴刻度间距
dev.off()


merged_data2$Fusobacterium_cut <- factor(
  merged_data2$Fusobacterium_cut,
  levels = c("High", "Low") # 包含所有可能值
)

library(survival)

# 创建Cox模型
cox_model <- coxph(Surv(survival_time, event_status) ~ Fusobacterium_cut, data = merged_data2)

# 查看模型摘要
summary(cox_model)


merged_data2$gender.x <- factor(merged_data2$gender.x)
sum(is.na(merged_data2$gender.x))

names(merged_data2)
unique(merged_data2$years_smoked)
merged_data2$years_smoked_history<-ifelse(merged_data2$years_smoked <= 20, "Low", "High")
merged_data2$portion_weight2 <- log(merged_data2$portion_weight + 1)
merged_data2$years_smoked2 <- log(merged_data2$years_smoked + 1)
merged_data2$years_smoked2[is.na(merged_data2$years_smoked2)] <- mean(merged_data2$years_smoked2, na.rm = TRUE)
merged_data2$Fusobacterium_cut <- relevel(merged_data2$Fusobacterium_cut, ref = "High")

merged_data2$age_at_index
merged_data2$alcohol_history


cox_model <- coxph(Surv(survival_time, event_status) ~ Fusobacterium_cut+age_at_index+
                     prior_malignancy+
                     portion_weight2+gender.x, data = merged_data2)

summary(cox_model)
f1_Fusobacterium<-ggforest(cox_model,
                           data=merged_data2,
                           main = "Hazard ratio",        # 标题
                           cpositions = c(0.01, 0.1, 0.3), #左边三列的相对距离
                           fontsize = 1, # 字体大小
                           refLabel = "reference", #显示因子的参考水平
                           noDigits = 3)  #小数点位数

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Apan-cancer_analysis_of_the_microbiome/TCGA-STAD胃腺癌/fusobacterium_STAD_ggforest_cut_多因素.pdf",f1_Fusobacterium,width =13,height=6)



#####Capnocytophaga TCGA-STAD----

merged_data2<- merged_data %>% filter(project %in% "TCGA-STAD")
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
  ifelse(merged_data2$Capnocytophaga_percentage >= 0.01048088, "High", "Low")


fit <- survfit(Surv(survival_time, event_status) ~Capnocytophaga_cut,  # 创建生存对象 
               data = merged_data2)

pdf("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Apan-cancer_analysis_of_the_microbiome/TCGA-STAD胃腺癌/Capnocytophaga_TCGA-STAD_curvep_cut.pdf", width = 6, height =7)
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
           break.x.by = 500,
           xlim = c(0, 3500),
           title = "Capnocytophaga STAD" )  # 设置x轴刻度间距
dev.off()


merged_data2$Capnocytophaga_cut <- factor(
  merged_data2$Capnocytophaga_cut,
  levels = c("High", "Low") # 包含所有可能值
)

library(survival)

# 创建Cox模型
cox_model <- coxph(Surv(survival_time, event_status) ~ Capnocytophaga_cut, data = merged_data2)

# 查看模型摘要
summary(cox_model)


merged_data2$gender.x <- factor(merged_data2$gender.x)
sum(is.na(merged_data2$gender.x))

names(merged_data2)
unique(merged_data2$years_smoked)
merged_data2$years_smoked_history<-ifelse(merged_data2$years_smoked <= 20, "Low", "High")
merged_data2$portion_weight2 <- log(merged_data2$portion_weight + 1)
merged_data2$years_smoked2 <- log(merged_data2$years_smoked + 1)
merged_data2$years_smoked2[is.na(merged_data2$years_smoked2)] <- mean(merged_data2$years_smoked2, na.rm = TRUE)
merged_data2$Fusobacterium_cut <- relevel(merged_data2$Fusobacterium_cut, ref = "High")

merged_data2$age_at_index
merged_data2$alcohol_history


cox_model <- coxph(Surv(survival_time, event_status) ~ Capnocytophaga_cut+age_at_index+
                     prior_malignancy+
                     portion_weight2+gender.x, data = merged_data2)

summary(cox_model)
f1_Capnocytophaga<-ggforest(cox_model,
                           data=merged_data2,
                           main = "Hazard ratio",        # 标题
                           cpositions = c(0.01, 0.1, 0.3), #左边三列的相对距离
                           fontsize = 1, # 字体大小
                           refLabel = "reference", #显示因子的参考水平
                           noDigits = 3)  #小数点位数

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Apan-cancer_analysis_of_the_microbiome/TCGA-STAD胃腺癌/Capnocytophaga_STAD_ggforest_cut_多因素.pdf",f1_Capnocytophaga,width =13,height=6)


#####Streptococcus TCGA-STAD----

merged_data2<- merged_data %>% filter(project %in% "TCGA-STAD")
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
  ifelse(merged_data2$Streptococcus_percentage >=6.525977, "High", "Low")


fit <- survfit(Surv(survival_time, event_status) ~Streptococcus_cut,  # 创建生存对象 
               data = merged_data2)

pdf("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Apan-cancer_analysis_of_the_microbiome/TCGA-STAD胃腺癌/Streptococcus_STAD_curvep_cut.pdf", width = 6, height =7)
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
           break.x.by = 500,
           xlim = c(0, 2000),
           title = "Streptococcus STAD" )  # 设置x轴刻度间距
dev.off()


merged_data2$Streptococcus_cut <- factor(
  merged_data2$Streptococcus_cut,
  levels = c("High", "Low") # 包含所有可能值
)

library(survival)

# 创建Cox模型
cox_model <- coxph(Surv(survival_time, event_status) ~ Streptococcus_cut, data = merged_data2)

# 查看模型摘要
summary(cox_model)


merged_data2$gender.x <- factor(merged_data2$gender.x)
sum(is.na(merged_data2$gender.x))

names(merged_data2)
unique(merged_data2$years_smoked)
merged_data2$years_smoked_history<-ifelse(merged_data2$years_smoked <= 20, "Low", "High")
merged_data2$portion_weight2 <- log(merged_data2$portion_weight + 1)
merged_data2$years_smoked2 <- log(merged_data2$years_smoked + 1)
merged_data2$years_smoked2[is.na(merged_data2$years_smoked2)] <- mean(merged_data2$years_smoked2, na.rm = TRUE)
merged_data2$Streptococcus_cut <- relevel(merged_data2$Streptococcus_cut, ref = "High")

merged_data2$age_at_index
merged_data2$alcohol_history


cox_model <- coxph(Surv(survival_time, event_status) ~ Streptococcus_cut+age_at_index+
                     prior_malignancy+
                     portion_weight2+gender.x, data = merged_data2)

summary(cox_model)
f1_Streptococcus<-ggforest(cox_model,
                           data=merged_data2,
                           main = "Hazard ratio",        # 标题
                           cpositions = c(0.01, 0.1, 0.3), #左边三列的相对距离
                           fontsize = 1, # 字体大小
                           refLabel = "reference", #显示因子的参考水平
                           noDigits = 3)  #小数点位数

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Apan-cancer_analysis_of_the_microbiome/TCGA-STAD胃腺癌/Streptococcus_STAD_ggforest_cut_多因素.pdf",f1_Streptococcus,width =13,height=6)


project_ids <- c("TCGA-ACC",
                 "TCGA-BLCA", "TCGA-BRCA", "TCGA-CESC", "TCGA-CHOL", "TCGA-COAD", "TCGA-DLBC", 
                 "TCGA-ESCA", "TCGA-GBM", "TCGA-HNSC", "TCGA-KICH", "TCGA-KIRC", "TCGA-KIRP", 
                 "TCGA-LAML", "TCGA-LGG", "TCGA-LIHC", "TCGA-LUAD", "TCGA-LUSC", "TCGA-MESO", 
                 "TCGA-OV", "TCGA-PAAD", "TCGA-PCPG", "TCGA-PRAD", "TCGA-READ", "TCGA-SARC", 
                 "TCGA-SKCM", "TCGA-STAD", "TCGA-TGCT", "TCGA-THCA", "TCGA-THYM", "TCGA-UCEC", 
                 "TCGA-UCS", "TCGA-UVM"
)


####fusobacterium HNSC 头颈鳞癌----
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
                     prior_malignancy+
                     portion_weight2+gender.x, data = merged_data2)

summary(cox_model)
cox.zph(cox_model)

cox_model <- coxph(Surv(survival_time, event_status) ~ Fusobacterium_percentage, data = merged_data2)

merged_data2$Fusobacterium_cut <- factor(
  merged_data2$Fusobacterium_cut,
  levels = c("Low", "High") # 包含所有可能值
)
merged_data2$Fusobacterium_cut <- factor(
  merged_data2$Fusobacterium_cut,
  levels = c("High", "Low") # 包含所有可能值
)

#merged_data2$Fusobacterium_cut <- relevel(merged_data2$Fusobacterium_cut, ref = "Low")
#merged_data2$Fusobacterium_cut <- relevel(merged_data2$Fusobacterium_cut, ref = "High")
cox_model <- coxph(Surv(survival_time, event_status) ~ Fusobacterium_cut, data = merged_data2)

summary(cox_model)

#coef exp(coef) se(coef)     z Pr(>|z|)   
#Fusobacterium_cutLow 0.6944    2.0026   0.2181 3.184  0.00145 **


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


####TCGAbiolinks + TCGAbiolinksGUI----

library(TCGAbiolinks)
library(survival)
library(survminer)

# 1. 下载临床数据
clin_query <- GDCquery(project = "TCGA-COAD", data.category = "Clinical", file.type = "xml")
GDCdownload(clin_query)
clin <- GDCprepare_clinic(clin_query, clinical.info = "patient")

# 2. 加载微生物丰度（假设已有 CSV）
mb <- read.csv("TCGA_COAD_microbiome_kraken.csv", row.names = 1)

# 3. 合并数据
df <- merge(clin, mb, by.x="bcr_patient_barcode", by.y="SampleID")
df$OS <- as.numeric(df$days_to_death)
df$status <- ifelse(df$vital_status == "Dead", 1, 0)

# 4. 以某 genus 举例做 Cox 模型
cox_res <- coxph(Surv(OS, status) ~ mb$Fusobacterium + age + gender, data=df)
summary(cox_res)


####TCGA重新换参考文献分析----
####Description: The Cancer Genome Atlas (TCGA)
##### Taxonomy contains only Kingdom & Genus level 
# due to missing intermediate taxonomy annotations
## https://doi.org/10.1016/j.canlet.2024.217077

library(sample_data)
library(ggplot2)
library(dplyr)
data(TCGA)
View(TCGA)

####ACC----

otu_table<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Pan-cancer atlas of tumor-resident microbiome, immunity and prognosis/ACC_otutable_bacteria_genus.csv",row.names=1,header=T)

View(otu_table)

otu_rel_abund <- otu_table / rowSums(otu_table)  

otu_rel_abund[is.na(otu_rel_abund)] <- 0

write.csv(otu_rel_abund,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Pan-cancer atlas of tumor-resident microbiome, immunity and prognosis/ACC_otu_rel_abund_genus.csv") 

View(otu_rel_abund)
# 提取 Fusobacterium 的丰度
library(dplyr)
fusobacterium_percentage <- otu_rel_abund %>% select(g__Fusobacterium)
# 提取 Capnocytophaga 的丰度
Capnocytophaga_percentage <- otu_rel_abund %>% select(g__Capnocytophaga)
# 提取 Streptococcus 的丰度
Streptococcus_percentage <- otu_rel_abund %>% select(g__Streptococcus)
data <- data.frame(
  Sample = rownames(otu_rel_abund),
  Fusobacterium_percentage = unlist(fusobacterium_percentage),
  Capnocytophaga_percentage=unlist(Capnocytophaga_percentage),
  Streptococcus_percentage =unlist(Streptococcus_percentage)
  
)

rownames(data)<- data$Sample

ACC_data<- data[,-1]



####BLCA----

otu_table<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Pan-cancer atlas of tumor-resident microbiome, immunity and prognosis/BLCA_otutable_bacteria_genus.csv",row.names=1,header=T)

View(otu_table)

otu_table <- otu_table[, !(colnames(otu_table) %in% "cancertype")]

otu_rel_abund <- otu_table / rowSums(otu_table)  

otu_rel_abund[is.na(otu_rel_abund)] <- 0

# 提取 Fusobacterium 的丰度
library(dplyr)
fusobacterium_percentage <- otu_rel_abund %>% select(g__Fusobacterium)
# 提取 Capnocytophaga 的丰度
Capnocytophaga_percentage <- otu_rel_abund %>% select(g__Capnocytophaga)
# 提取 Streptococcus 的丰度
Streptococcus_percentage <- otu_rel_abund %>% select(g__Streptococcus)
data <- data.frame(
  Sample = rownames(otu_rel_abund),
  Fusobacterium_percentage = unlist(fusobacterium_percentage),
  Capnocytophaga_percentage=unlist(Capnocytophaga_percentage),
  Streptococcus_percentage =unlist(Streptococcus_percentage)
  
)

rownames(data)<- data$Sample

BLCA_data<- data[,-1]


####BRCA----

otu_table<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Pan-cancer atlas of tumor-resident microbiome, immunity and prognosis/BRCA_otutable_bacteria_genus.csv",row.names=1,header=T)

View(otu_table)

otu_table <- otu_table[, !(colnames(otu_table) %in% "cancertype")]

otu_rel_abund <- otu_table / rowSums(otu_table)  

otu_rel_abund[is.na(otu_rel_abund)] <- 0

# 提取 Fusobacterium 的丰度
library(dplyr)
fusobacterium_percentage <- otu_rel_abund %>% select(g__Fusobacterium)
# 提取 Capnocytophaga 的丰度
Capnocytophaga_percentage <- otu_rel_abund %>% select(g__Capnocytophaga)
# 提取 Streptococcus 的丰度
Streptococcus_percentage <- otu_rel_abund %>% select(g__Streptococcus)
data <- data.frame(
  Sample = rownames(otu_rel_abund),
  Fusobacterium_percentage = unlist(fusobacterium_percentage),
  Capnocytophaga_percentage=unlist(Capnocytophaga_percentage),
  Streptococcus_percentage =unlist(Streptococcus_percentage)
  
)

rownames(data)<- data$Sample

BRCA_data<- data[,-1]



####CESC----

otu_table<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Pan-cancer atlas of tumor-resident microbiome, immunity and prognosis/CESC_otutable_bacteria_genus.csv",row.names=1,header=T)

View(otu_table)

otu_table <- otu_table[, !(colnames(otu_table) %in% "cancertype")]

otu_rel_abund <- otu_table / rowSums(otu_table)  

otu_rel_abund[is.na(otu_rel_abund)] <- 0

# 提取 Fusobacterium 的丰度
library(dplyr)
fusobacterium_percentage <- otu_rel_abund %>% select(g__Fusobacterium)
# 提取 Capnocytophaga 的丰度
Capnocytophaga_percentage <- otu_rel_abund %>% select(g__Capnocytophaga)
# 提取 Streptococcus 的丰度
Streptococcus_percentage <- otu_rel_abund %>% select(g__Streptococcus)
data <- data.frame(
  Sample = rownames(otu_rel_abund),
  Fusobacterium_percentage = unlist(fusobacterium_percentage),
  Capnocytophaga_percentage=unlist(Capnocytophaga_percentage),
  Streptococcus_percentage =unlist(Streptococcus_percentage)
  
)

rownames(data)<- data$Sample

CESC_data<- data[,-1]





####CHOL----

otu_table<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Pan-cancer atlas of tumor-resident microbiome, immunity and prognosis/CHOL_otutable_bacteria_genus.csv",row.names=1,header=T)

View(otu_table)

otu_table <- otu_table[, !(colnames(otu_table) %in% "cancertype")]

otu_rel_abund <- otu_table / rowSums(otu_table)  

otu_rel_abund[is.na(otu_rel_abund)] <- 0

# 提取 Fusobacterium 的丰度
library(dplyr)
fusobacterium_percentage <- otu_rel_abund %>% select(g__Fusobacterium)
# 提取 Capnocytophaga 的丰度
Capnocytophaga_percentage <- otu_rel_abund %>% select(g__Capnocytophaga)
# 提取 Streptococcus 的丰度
Streptococcus_percentage <- otu_rel_abund %>% select(g__Streptococcus)
data <- data.frame(
  Sample = rownames(otu_rel_abund),
  Fusobacterium_percentage = unlist(fusobacterium_percentage),
  Capnocytophaga_percentage=unlist(Capnocytophaga_percentage),
  Streptococcus_percentage =unlist(Streptococcus_percentage)
  
)

rownames(data)<- data$Sample

CHOL_data<- data[,-1]



####COAD----

otu_table<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Pan-cancer atlas of tumor-resident microbiome, immunity and prognosis/COAD_otutable_bacteria_genus.csv",row.names=1,header=T)

View(otu_table)

otu_table <- otu_table[, !(colnames(otu_table) %in% "cancertype")]

otu_rel_abund <- otu_table / rowSums(otu_table)  

otu_rel_abund[is.na(otu_rel_abund)] <- 0

# 提取 Fusobacterium 的丰度
library(dplyr)
fusobacterium_percentage <- otu_rel_abund %>% select(g__Fusobacterium)
# 提取 Capnocytophaga 的丰度
Capnocytophaga_percentage <- otu_rel_abund %>% select(g__Capnocytophaga)
# 提取 Streptococcus 的丰度
Streptococcus_percentage <- otu_rel_abund %>% select(g__Streptococcus)
data <- data.frame(
  Sample = rownames(otu_rel_abund),
  Fusobacterium_percentage = unlist(fusobacterium_percentage),
  Capnocytophaga_percentage=unlist(Capnocytophaga_percentage),
  Streptococcus_percentage =unlist(Streptococcus_percentage)
  
)

rownames(data)<- data$Sample

COAD_data<- data[,-1]



####DLBC----

otu_table<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Pan-cancer atlas of tumor-resident microbiome, immunity and prognosis/DLBC_otutable_bacteria_genus.csv",row.names=1,header=T)

View(otu_table)

otu_table <- otu_table[, !(colnames(otu_table) %in% "cancertype")]

otu_rel_abund <- otu_table / rowSums(otu_table)  

otu_rel_abund[is.na(otu_rel_abund)] <- 0

# 提取 Fusobacterium 的丰度
library(dplyr)
fusobacterium_percentage <- otu_rel_abund %>% select(g__Fusobacterium)
# 提取 Capnocytophaga 的丰度
Capnocytophaga_percentage <- otu_rel_abund %>% select(g__Capnocytophaga)
# 提取 Streptococcus 的丰度
Streptococcus_percentage <- otu_rel_abund %>% select(g__Streptococcus)
data <- data.frame(
  Sample = rownames(otu_rel_abund),
  Fusobacterium_percentage = unlist(fusobacterium_percentage),
  Capnocytophaga_percentage=unlist(Capnocytophaga_percentage),
  Streptococcus_percentage =unlist(Streptococcus_percentage)
  
)

rownames(data)<- data$Sample

DLBC_data<- data[,-1]


####ESCA----

otu_table<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Pan-cancer atlas of tumor-resident microbiome, immunity and prognosis/ESCA_otutable_bacteria_genus.csv",row.names=1,header=T)

View(otu_table)

otu_table <- otu_table[, !(colnames(otu_table) %in% "cancertype")]

otu_rel_abund <- otu_table / rowSums(otu_table)  

otu_rel_abund[is.na(otu_rel_abund)] <- 0

# 提取 Fusobacterium 的丰度
library(dplyr)
fusobacterium_percentage <- otu_rel_abund %>% select(g__Fusobacterium)
# 提取 Capnocytophaga 的丰度
Capnocytophaga_percentage <- otu_rel_abund %>% select(g__Capnocytophaga)
# 提取 Streptococcus 的丰度
Streptococcus_percentage <- otu_rel_abund %>% select(g__Streptococcus)
data <- data.frame(
  Sample = rownames(otu_rel_abund),
  Fusobacterium_percentage = unlist(fusobacterium_percentage),
  Capnocytophaga_percentage=unlist(Capnocytophaga_percentage),
  Streptococcus_percentage =unlist(Streptococcus_percentage)
  
)

rownames(data)<- data$Sample
ESCA_data<- data[,-1]


####GBM----

otu_table<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Pan-cancer atlas of tumor-resident microbiome, immunity and prognosis/GBM_otutable_bacteria_genus.csv",row.names=1,header=T)

View(otu_table)

otu_table <- otu_table[, !(colnames(otu_table) %in% "cancertype")]

otu_rel_abund <- otu_table / rowSums(otu_table)  

otu_rel_abund[is.na(otu_rel_abund)] <- 0

# 提取 Fusobacterium 的丰度
library(dplyr)
fusobacterium_percentage <- otu_rel_abund %>% select(g__Fusobacterium)
# 提取 Capnocytophaga 的丰度
Capnocytophaga_percentage <- otu_rel_abund %>% select(g__Capnocytophaga)
# 提取 Streptococcus 的丰度
Streptococcus_percentage <- otu_rel_abund %>% select(g__Streptococcus)
data <- data.frame(
  Sample = rownames(otu_rel_abund),
  Fusobacterium_percentage = unlist(fusobacterium_percentage),
  Capnocytophaga_percentage=unlist(Capnocytophaga_percentage),
  Streptococcus_percentage =unlist(Streptococcus_percentage)
  
)

rownames(data)<- data$Sample
GBM_data<- data[,-1]


####HNSC----

otu_table<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Pan-cancer atlas of tumor-resident microbiome, immunity and prognosis/HNSC_otutable_bacteria_genus.csv",row.names=1,header=T)

View(otu_table)

otu_table <- otu_table[, !(colnames(otu_table) %in% "cancertype")]

otu_rel_abund <- otu_table / rowSums(otu_table)  

otu_rel_abund[is.na(otu_rel_abund)] <- 0

# 提取 Fusobacterium 的丰度
library(dplyr)
fusobacterium_percentage <- otu_rel_abund %>% select(g__Fusobacterium)
# 提取 Capnocytophaga 的丰度
Capnocytophaga_percentage <- otu_rel_abund %>% select(g__Capnocytophaga)
# 提取 Streptococcus 的丰度
Streptococcus_percentage <- otu_rel_abund %>% select(g__Streptococcus)
data <- data.frame(
  Sample = rownames(otu_rel_abund),
  Fusobacterium_percentage = unlist(fusobacterium_percentage),
  Capnocytophaga_percentage=unlist(Capnocytophaga_percentage),
  Streptococcus_percentage =unlist(Streptococcus_percentage)
  
)

rownames(data)<- data$Sample
HNSC_data<- data[,-1]


####KICH----

otu_table<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Pan-cancer atlas of tumor-resident microbiome, immunity and prognosis/KICH_otutable_bacteria_genus.csv",row.names=1,header=T)

View(otu_table)

otu_table <- otu_table[, !(colnames(otu_table) %in% "cancertype")]

otu_rel_abund <- otu_table / rowSums(otu_table)  

otu_rel_abund[is.na(otu_rel_abund)] <- 0

# 提取 Fusobacterium 的丰度
library(dplyr)
fusobacterium_percentage <- otu_rel_abund %>% select(g__Fusobacterium)
# 提取 Capnocytophaga 的丰度
Capnocytophaga_percentage <- otu_rel_abund %>% select(g__Capnocytophaga)
# 提取 Streptococcus 的丰度
Streptococcus_percentage <- otu_rel_abund %>% select(g__Streptococcus)
data <- data.frame(
  Sample = rownames(otu_rel_abund),
  Fusobacterium_percentage = unlist(fusobacterium_percentage),
  Capnocytophaga_percentage=unlist(Capnocytophaga_percentage),
  Streptococcus_percentage =unlist(Streptococcus_percentage)
  
)

rownames(data)<- data$Sample
KICH_data<- data[,-1]


####KIRC----

otu_table<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Pan-cancer atlas of tumor-resident microbiome, immunity and prognosis/KIRC_otutable_bacteria_genus.csv",row.names=1,header=T)

View(otu_table)

otu_table <- otu_table[, !(colnames(otu_table) %in% "cancertype")]

otu_rel_abund <- otu_table / rowSums(otu_table)  

otu_rel_abund[is.na(otu_rel_abund)] <- 0

# 提取 Fusobacterium 的丰度
library(dplyr)
fusobacterium_percentage <- otu_rel_abund %>% select(g__Fusobacterium)
# 提取 Capnocytophaga 的丰度
Capnocytophaga_percentage <- otu_rel_abund %>% select(g__Capnocytophaga)
# 提取 Streptococcus 的丰度
Streptococcus_percentage <- otu_rel_abund %>% select(g__Streptococcus)
data <- data.frame(
  Sample = rownames(otu_rel_abund),
  Fusobacterium_percentage = unlist(fusobacterium_percentage),
  Capnocytophaga_percentage=unlist(Capnocytophaga_percentage),
  Streptococcus_percentage =unlist(Streptococcus_percentage)
  
)

rownames(data)<- data$Sample
KIRC_data<- data[,-1]


####KIRP----

otu_table<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Pan-cancer atlas of tumor-resident microbiome, immunity and prognosis/KIRP_otutable_bacteria_genus.csv",row.names=1,header=T)

View(otu_table)

otu_table <- otu_table[, !(colnames(otu_table) %in% "cancertype")]

otu_rel_abund <- otu_table / rowSums(otu_table)  

otu_rel_abund[is.na(otu_rel_abund)] <- 0

# 提取 Fusobacterium 的丰度
library(dplyr)
fusobacterium_percentage <- otu_rel_abund %>% select(g__Fusobacterium)
# 提取 Capnocytophaga 的丰度
Capnocytophaga_percentage <- otu_rel_abund %>% select(g__Capnocytophaga)
# 提取 Streptococcus 的丰度
Streptococcus_percentage <- otu_rel_abund %>% select(g__Streptococcus)
data <- data.frame(
  Sample = rownames(otu_rel_abund),
  Fusobacterium_percentage = unlist(fusobacterium_percentage),
  Capnocytophaga_percentage=unlist(Capnocytophaga_percentage),
  Streptococcus_percentage =unlist(Streptococcus_percentage)
  
)

rownames(data)<- data$Sample
KIRP_data<- data[,-1]




####LAML----

otu_table<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Pan-cancer atlas of tumor-resident microbiome, immunity and prognosis/LAML_otutable_bacteria_genus.csv",row.names=1,header=T)

View(otu_table)

otu_table <- otu_table[, !(colnames(otu_table) %in% "cancertype")]

otu_rel_abund <- otu_table / rowSums(otu_table)  

otu_rel_abund[is.na(otu_rel_abund)] <- 0

# 提取 Fusobacterium 的丰度
library(dplyr)
fusobacterium_percentage <- otu_rel_abund %>% select(g__Fusobacterium)
# 提取 Capnocytophaga 的丰度
Capnocytophaga_percentage <- otu_rel_abund %>% select(g__Capnocytophaga)
# 提取 Streptococcus 的丰度
Streptococcus_percentage <- otu_rel_abund %>% select(g__Streptococcus)
data <- data.frame(
  Sample = rownames(otu_rel_abund),
  Fusobacterium_percentage = unlist(fusobacterium_percentage),
  Capnocytophaga_percentage=unlist(Capnocytophaga_percentage),
  Streptococcus_percentage =unlist(Streptococcus_percentage)
  
)

rownames(data)<- data$Sample
LAML_data<- data[,-1]



####LGG----

otu_table<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Pan-cancer atlas of tumor-resident microbiome, immunity and prognosis/LGG_otutable_bacteria_genus.csv",row.names=1,header=T)

View(otu_table)

otu_table <- otu_table[, !(colnames(otu_table) %in% "cancertype")]

otu_rel_abund <- otu_table / rowSums(otu_table)  

otu_rel_abund[is.na(otu_rel_abund)] <- 0

# 提取 Fusobacterium 的丰度
library(dplyr)
fusobacterium_percentage <- otu_rel_abund %>% select(g__Fusobacterium)
# 提取 Capnocytophaga 的丰度
Capnocytophaga_percentage <- otu_rel_abund %>% select(g__Capnocytophaga)
# 提取 Streptococcus 的丰度
Streptococcus_percentage <- otu_rel_abund %>% select(g__Streptococcus)
data <- data.frame(
  Sample = rownames(otu_rel_abund),
  Fusobacterium_percentage = unlist(fusobacterium_percentage),
  Capnocytophaga_percentage=unlist(Capnocytophaga_percentage),
  Streptococcus_percentage =unlist(Streptococcus_percentage)
  
)

rownames(data)<- data$Sample
LGG_data<- data[,-1]


####LIHC----

otu_table<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Pan-cancer atlas of tumor-resident microbiome, immunity and prognosis/LIHC_otutable_bacteria_genus.csv",row.names=1,header=T)

View(otu_table)

otu_table <- otu_table[, !(colnames(otu_table) %in% "cancertype")]

otu_rel_abund <- otu_table / rowSums(otu_table)  

otu_rel_abund[is.na(otu_rel_abund)] <- 0

# 提取 Fusobacterium 的丰度
library(dplyr)
fusobacterium_percentage <- otu_rel_abund %>% select(g__Fusobacterium)
# 提取 Capnocytophaga 的丰度
Capnocytophaga_percentage <- otu_rel_abund %>% select(g__Capnocytophaga)
# 提取 Streptococcus 的丰度
Streptococcus_percentage <- otu_rel_abund %>% select(g__Streptococcus)
data <- data.frame(
  Sample = rownames(otu_rel_abund),
  Fusobacterium_percentage = unlist(fusobacterium_percentage),
  Capnocytophaga_percentage=unlist(Capnocytophaga_percentage),
  Streptococcus_percentage =unlist(Streptococcus_percentage)
  
)

rownames(data)<- data$Sample
LIHC_data<- data[,-1]



####LUAD----

otu_table<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Pan-cancer atlas of tumor-resident microbiome, immunity and prognosis/LUAD_otutable_bacteria_genus.csv",row.names=1,header=T)

View(otu_table)

otu_table <- otu_table[, !(colnames(otu_table) %in% "cancertype")]

otu_rel_abund <- otu_table / rowSums(otu_table)  

otu_rel_abund[is.na(otu_rel_abund)] <- 0

# 提取 Fusobacterium 的丰度
library(dplyr)
fusobacterium_percentage <- otu_rel_abund %>% select(g__Fusobacterium)
# 提取 Capnocytophaga 的丰度
Capnocytophaga_percentage <- otu_rel_abund %>% select(g__Capnocytophaga)
# 提取 Streptococcus 的丰度
Streptococcus_percentage <- otu_rel_abund %>% select(g__Streptococcus)
data <- data.frame(
  Sample = rownames(otu_rel_abund),
  Fusobacterium_percentage = unlist(fusobacterium_percentage),
  Capnocytophaga_percentage=unlist(Capnocytophaga_percentage),
  Streptococcus_percentage =unlist(Streptococcus_percentage)
  
)

rownames(data)<- data$Sample
LUAD_data<- data[,-1]



####LUSC----

otu_table<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Pan-cancer atlas of tumor-resident microbiome, immunity and prognosis/LUSC_otutable_bacteria_genus.csv",row.names=1,header=T)

View(otu_table)

otu_table <- otu_table[, !(colnames(otu_table) %in% "cancertype")]

otu_rel_abund <- otu_table / rowSums(otu_table)  

otu_rel_abund[is.na(otu_rel_abund)] <- 0

# 提取 Fusobacterium 的丰度
library(dplyr)
fusobacterium_percentage <- otu_rel_abund %>% select(g__Fusobacterium)
# 提取 Capnocytophaga 的丰度
Capnocytophaga_percentage <- otu_rel_abund %>% select(g__Capnocytophaga)
# 提取 Streptococcus 的丰度
Streptococcus_percentage <- otu_rel_abund %>% select(g__Streptococcus)
data <- data.frame(
  Sample = rownames(otu_rel_abund),
  Fusobacterium_percentage = unlist(fusobacterium_percentage),
  Capnocytophaga_percentage=unlist(Capnocytophaga_percentage),
  Streptococcus_percentage =unlist(Streptococcus_percentage)
  
)

rownames(data)<- data$Sample
LUSC_data<- data[,-1]



####MESO----

otu_table<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Pan-cancer atlas of tumor-resident microbiome, immunity and prognosis/MESO_otutable_bacteria_genus.csv",row.names=1,header=T)

View(otu_table)

otu_table <- otu_table[, !(colnames(otu_table) %in% "cancertype")]

otu_rel_abund <- otu_table / rowSums(otu_table)  

otu_rel_abund[is.na(otu_rel_abund)] <- 0

# 提取 Fusobacterium 的丰度
library(dplyr)
fusobacterium_percentage <- otu_rel_abund %>% select(g__Fusobacterium)
# 提取 Capnocytophaga 的丰度
Capnocytophaga_percentage <- otu_rel_abund %>% select(g__Capnocytophaga)
# 提取 Streptococcus 的丰度
Streptococcus_percentage <- otu_rel_abund %>% select(g__Streptococcus)
data <- data.frame(
  Sample = rownames(otu_rel_abund),
  Fusobacterium_percentage = unlist(fusobacterium_percentage),
  Capnocytophaga_percentage=unlist(Capnocytophaga_percentage),
  Streptococcus_percentage =unlist(Streptococcus_percentage)
  
)

rownames(data)<- data$Sample
MESO_data<- data[,-1]


####OV----

otu_table<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Pan-cancer atlas of tumor-resident microbiome, immunity and prognosis/OV_otutable_bacteria_genus.csv",row.names=1,header=T)

View(otu_table)

otu_table <- otu_table[, !(colnames(otu_table) %in% "cancertype")]

otu_rel_abund <- otu_table / rowSums(otu_table)  

otu_rel_abund[is.na(otu_rel_abund)] <- 0

# 提取 Fusobacterium 的丰度
library(dplyr)
fusobacterium_percentage <- otu_rel_abund %>% select(g__Fusobacterium)
# 提取 Capnocytophaga 的丰度
Capnocytophaga_percentage <- otu_rel_abund %>% select(g__Capnocytophaga)
# 提取 Streptococcus 的丰度
Streptococcus_percentage <- otu_rel_abund %>% select(g__Streptococcus)
data <- data.frame(
  Sample = rownames(otu_rel_abund),
  Fusobacterium_percentage = unlist(fusobacterium_percentage),
  Capnocytophaga_percentage=unlist(Capnocytophaga_percentage),
  Streptococcus_percentage =unlist(Streptococcus_percentage)
  
)

rownames(data)<- data$Sample
OV_data<- data[,-1]


####PAAD----

otu_table<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Pan-cancer atlas of tumor-resident microbiome, immunity and prognosis/PAAD_otutable_bacteria_genus.csv",row.names=1,header=T)

View(otu_table)

otu_table <- otu_table[, !(colnames(otu_table) %in% "cancertype")]

otu_rel_abund <- otu_table / rowSums(otu_table)  

otu_rel_abund[is.na(otu_rel_abund)] <- 0

# 提取 Fusobacterium 的丰度
library(dplyr)
fusobacterium_percentage <- otu_rel_abund %>% select(g__Fusobacterium)
# 提取 Capnocytophaga 的丰度
Capnocytophaga_percentage <- otu_rel_abund %>% select(g__Capnocytophaga)
# 提取 Streptococcus 的丰度
Streptococcus_percentage <- otu_rel_abund %>% select(g__Streptococcus)
data <- data.frame(
  Sample = rownames(otu_rel_abund),
  Fusobacterium_percentage = unlist(fusobacterium_percentage),
  Capnocytophaga_percentage=unlist(Capnocytophaga_percentage),
  Streptococcus_percentage =unlist(Streptococcus_percentage)
  
)

rownames(data)<- data$Sample
PAAD_data<- data[,-1]



####PCPG----

otu_table<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Pan-cancer atlas of tumor-resident microbiome, immunity and prognosis/PCPG_otutable_bacteria_genus.csv",row.names=1,header=T)

View(otu_table)

otu_table <- otu_table[, !(colnames(otu_table) %in% "cancertype")]

otu_rel_abund <- otu_table / rowSums(otu_table)  

otu_rel_abund[is.na(otu_rel_abund)] <- 0

# 提取 Fusobacterium 的丰度
library(dplyr)
fusobacterium_percentage <- otu_rel_abund %>% select(g__Fusobacterium)
# 提取 Capnocytophaga 的丰度
Capnocytophaga_percentage <- otu_rel_abund %>% select(g__Capnocytophaga)
# 提取 Streptococcus 的丰度
Streptococcus_percentage <- otu_rel_abund %>% select(g__Streptococcus)
data <- data.frame(
  Sample = rownames(otu_rel_abund),
  Fusobacterium_percentage = unlist(fusobacterium_percentage),
  Capnocytophaga_percentage=unlist(Capnocytophaga_percentage),
  Streptococcus_percentage =unlist(Streptococcus_percentage)
  
)

rownames(data)<- data$Sample
PCPG_data<- data[,-1]


####PRAD----

otu_table<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Pan-cancer atlas of tumor-resident microbiome, immunity and prognosis/PRAD_otutable_bacteria_genus.csv",row.names=1,header=T)

View(otu_table)

otu_table <- otu_table[, !(colnames(otu_table) %in% "cancertype")]

otu_rel_abund <- otu_table / rowSums(otu_table)  

otu_rel_abund[is.na(otu_rel_abund)] <- 0

# 提取 Fusobacterium 的丰度
library(dplyr)
fusobacterium_percentage <- otu_rel_abund %>% select(g__Fusobacterium)
# 提取 Capnocytophaga 的丰度
Capnocytophaga_percentage <- otu_rel_abund %>% select(g__Capnocytophaga)
# 提取 Streptococcus 的丰度
Streptococcus_percentage <- otu_rel_abund %>% select(g__Streptococcus)
data <- data.frame(
  Sample = rownames(otu_rel_abund),
  Fusobacterium_percentage = unlist(fusobacterium_percentage),
  Capnocytophaga_percentage=unlist(Capnocytophaga_percentage),
  Streptococcus_percentage =unlist(Streptococcus_percentage)
  
)

rownames(data)<- data$Sample
PRAD_data<- data[,-1]



####READ----

otu_table<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Pan-cancer atlas of tumor-resident microbiome, immunity and prognosis/READ_otutable_bacteria_genus.csv",row.names=1,header=T)

View(otu_table)

otu_table <- otu_table[, !(colnames(otu_table) %in% "cancertype")]

otu_rel_abund <- otu_table / rowSums(otu_table)  

otu_rel_abund[is.na(otu_rel_abund)] <- 0

# 提取 Fusobacterium 的丰度
library(dplyr)
fusobacterium_percentage <- otu_rel_abund %>% select(g__Fusobacterium)
# 提取 Capnocytophaga 的丰度
Capnocytophaga_percentage <- otu_rel_abund %>% select(g__Capnocytophaga)
# 提取 Streptococcus 的丰度
Streptococcus_percentage <- otu_rel_abund %>% select(g__Streptococcus)
data <- data.frame(
  Sample = rownames(otu_rel_abund),
  Fusobacterium_percentage = unlist(fusobacterium_percentage),
  Capnocytophaga_percentage=unlist(Capnocytophaga_percentage),
  Streptococcus_percentage =unlist(Streptococcus_percentage)
  
)

rownames(data)<- data$Sample
READ_data<- data[,-1]


####SARC----

otu_table<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Pan-cancer atlas of tumor-resident microbiome, immunity and prognosis/SARC_otutable_bacteria_genus.csv",row.names=1,header=T)

View(otu_table)

otu_table <- otu_table[, !(colnames(otu_table) %in% "cancertype")]

otu_rel_abund <- otu_table / rowSums(otu_table)  

otu_rel_abund[is.na(otu_rel_abund)] <- 0

# 提取 Fusobacterium 的丰度
library(dplyr)
fusobacterium_percentage <- otu_rel_abund %>% select(g__Fusobacterium)
# 提取 Capnocytophaga 的丰度
Capnocytophaga_percentage <- otu_rel_abund %>% select(g__Capnocytophaga)
# 提取 Streptococcus 的丰度
Streptococcus_percentage <- otu_rel_abund %>% select(g__Streptococcus)
data <- data.frame(
  Sample = rownames(otu_rel_abund),
  Fusobacterium_percentage = unlist(fusobacterium_percentage),
  Capnocytophaga_percentage=unlist(Capnocytophaga_percentage),
  Streptococcus_percentage =unlist(Streptococcus_percentage)
  
)

rownames(data)<- data$Sample
SARC_data<- data[,-1]


####SKCM----

otu_table<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Pan-cancer atlas of tumor-resident microbiome, immunity and prognosis/SKCM_otutable_bacteria_genus.csv",row.names=1,header=T)

View(otu_table)

otu_table <- otu_table[, !(colnames(otu_table) %in% "cancertype")]

otu_rel_abund <- otu_table / rowSums(otu_table)  

otu_rel_abund[is.na(otu_rel_abund)] <- 0

# 提取 Fusobacterium 的丰度
library(dplyr)
fusobacterium_percentage <- otu_rel_abund %>% select(g__Fusobacterium)
# 提取 Capnocytophaga 的丰度
Capnocytophaga_percentage <- otu_rel_abund %>% select(g__Capnocytophaga)
# 提取 Streptococcus 的丰度
Streptococcus_percentage <- otu_rel_abund %>% select(g__Streptococcus)
data <- data.frame(
  Sample = rownames(otu_rel_abund),
  Fusobacterium_percentage = unlist(fusobacterium_percentage),
  Capnocytophaga_percentage=unlist(Capnocytophaga_percentage),
  Streptococcus_percentage =unlist(Streptococcus_percentage)
  
)

rownames(data)<- data$Sample
SKCM_data<- data[,-1]



####STAD----

otu_table<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Pan-cancer atlas of tumor-resident microbiome, immunity and prognosis/STAD_otutable_bacteria_genus.csv",row.names=1,header=T)

View(otu_table)

otu_table <- otu_table[, !(colnames(otu_table) %in% "cancertype")]

otu_rel_abund <- otu_table / rowSums(otu_table)  

otu_rel_abund[is.na(otu_rel_abund)] <- 0

# 提取 Fusobacterium 的丰度
library(dplyr)
fusobacterium_percentage <- otu_rel_abund %>% select(g__Fusobacterium)
# 提取 Capnocytophaga 的丰度
Capnocytophaga_percentage <- otu_rel_abund %>% select(g__Capnocytophaga)
# 提取 Streptococcus 的丰度
Streptococcus_percentage <- otu_rel_abund %>% select(g__Streptococcus)
data <- data.frame(
  Sample = rownames(otu_rel_abund),
  Fusobacterium_percentage = unlist(fusobacterium_percentage),
  Capnocytophaga_percentage=unlist(Capnocytophaga_percentage),
  Streptococcus_percentage =unlist(Streptococcus_percentage)
  
)

rownames(data)<- data$Sample
STAD_data<- data[,-1]



####TGCT----

otu_table<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Pan-cancer atlas of tumor-resident microbiome, immunity and prognosis/TGCT_otutable_bacteria_genus.csv",row.names=1,header=T)

View(otu_table)

otu_table <- otu_table[, !(colnames(otu_table) %in% "cancertype")]

otu_rel_abund <- otu_table / rowSums(otu_table)  

otu_rel_abund[is.na(otu_rel_abund)] <- 0

# 提取 Fusobacterium 的丰度
library(dplyr)
fusobacterium_percentage <- otu_rel_abund %>% select(g__Fusobacterium)
# 提取 Capnocytophaga 的丰度
Capnocytophaga_percentage <- otu_rel_abund %>% select(g__Capnocytophaga)
# 提取 Streptococcus 的丰度
Streptococcus_percentage <- otu_rel_abund %>% select(g__Streptococcus)
data <- data.frame(
  Sample = rownames(otu_rel_abund),
  Fusobacterium_percentage = unlist(fusobacterium_percentage),
  Capnocytophaga_percentage=unlist(Capnocytophaga_percentage),
  Streptococcus_percentage =unlist(Streptococcus_percentage)
  
)

rownames(data)<- data$Sample
TGCT_data<- data[,-1]


####THCA----

otu_table<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Pan-cancer atlas of tumor-resident microbiome, immunity and prognosis/THCA_otutable_bacteria_genus.csv",row.names=1,header=T)

View(otu_table)

otu_table <- otu_table[, !(colnames(otu_table) %in% "cancertype")]

otu_rel_abund <- otu_table / rowSums(otu_table)  

otu_rel_abund[is.na(otu_rel_abund)] <- 0

# 提取 Fusobacterium 的丰度
library(dplyr)
fusobacterium_percentage <- otu_rel_abund %>% select(g__Fusobacterium)
# 提取 Capnocytophaga 的丰度
Capnocytophaga_percentage <- otu_rel_abund %>% select(g__Capnocytophaga)
# 提取 Streptococcus 的丰度
Streptococcus_percentage <- otu_rel_abund %>% select(g__Streptococcus)
data <- data.frame(
  Sample = rownames(otu_rel_abund),
  Fusobacterium_percentage = unlist(fusobacterium_percentage),
  Capnocytophaga_percentage=unlist(Capnocytophaga_percentage),
  Streptococcus_percentage =unlist(Streptococcus_percentage)
  
)

rownames(data)<- data$Sample
THCA_data<- data[,-1]



####THYM----

otu_table<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Pan-cancer atlas of tumor-resident microbiome, immunity and prognosis/THYM_otutable_bacteria_genus.csv",row.names=1,header=T)

View(otu_table)

otu_table <- otu_table[, !(colnames(otu_table) %in% "cancertype")]

otu_rel_abund <- otu_table / rowSums(otu_table)  

otu_rel_abund[is.na(otu_rel_abund)] <- 0

# 提取 Fusobacterium 的丰度
library(dplyr)
fusobacterium_percentage <- otu_rel_abund %>% select(g__Fusobacterium)
# 提取 Capnocytophaga 的丰度
Capnocytophaga_percentage <- otu_rel_abund %>% select(g__Capnocytophaga)
# 提取 Streptococcus 的丰度
Streptococcus_percentage <- otu_rel_abund %>% select(g__Streptococcus)
data <- data.frame(
  Sample = rownames(otu_rel_abund),
  Fusobacterium_percentage = unlist(fusobacterium_percentage),
  Capnocytophaga_percentage=unlist(Capnocytophaga_percentage),
  Streptococcus_percentage =unlist(Streptococcus_percentage)
  
)

rownames(data)<- data$Sample
THYM_data<- data[,-1]


####UCEC----

otu_table<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Pan-cancer atlas of tumor-resident microbiome, immunity and prognosis/UCEC_otutable_bacteria_genus.csv",row.names=1,header=T)

View(otu_table)

otu_table <- otu_table[, !(colnames(otu_table) %in% "cancertype")]

otu_rel_abund <- otu_table / rowSums(otu_table)  

otu_rel_abund[is.na(otu_rel_abund)] <- 0

# 提取 Fusobacterium 的丰度
library(dplyr)
fusobacterium_percentage <- otu_rel_abund %>% select(g__Fusobacterium)
# 提取 Capnocytophaga 的丰度
Capnocytophaga_percentage <- otu_rel_abund %>% select(g__Capnocytophaga)
# 提取 Streptococcus 的丰度
Streptococcus_percentage <- otu_rel_abund %>% select(g__Streptococcus)
data <- data.frame(
  Sample = rownames(otu_rel_abund),
  Fusobacterium_percentage = unlist(fusobacterium_percentage),
  Capnocytophaga_percentage=unlist(Capnocytophaga_percentage),
  Streptococcus_percentage =unlist(Streptococcus_percentage)
  
)

rownames(data)<- data$Sample
UCEC_data<- data[,-1]


####UCS----

otu_table<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Pan-cancer atlas of tumor-resident microbiome, immunity and prognosis/UCS_otutable_bacteria_genus.csv",row.names=1,header=T)

View(otu_table)

otu_table <- otu_table[, !(colnames(otu_table) %in% "cancertype")]

otu_rel_abund <- otu_table / rowSums(otu_table)  

otu_rel_abund[is.na(otu_rel_abund)] <- 0

# 提取 Fusobacterium 的丰度
library(dplyr)
fusobacterium_percentage <- otu_rel_abund %>% select(g__Fusobacterium)
# 提取 Capnocytophaga 的丰度
Capnocytophaga_percentage <- otu_rel_abund %>% select(g__Capnocytophaga)
# 提取 Streptococcus 的丰度
Streptococcus_percentage <- otu_rel_abund %>% select(g__Streptococcus)
data <- data.frame(
  Sample = rownames(otu_rel_abund),
  Fusobacterium_percentage = unlist(fusobacterium_percentage),
  Capnocytophaga_percentage=unlist(Capnocytophaga_percentage),
  Streptococcus_percentage =unlist(Streptococcus_percentage)
  
)

rownames(data)<- data$Sample
UCS_data<- data[,-1]
head(UCS_data)


####UVM----

otu_table<- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Pan-cancer atlas of tumor-resident microbiome, immunity and prognosis/UVM_otutable_bacteria_genus.csv",row.names=1,header=T)

View(otu_table)

otu_table <- otu_table[, !(colnames(otu_table) %in% "cancertype")]

otu_rel_abund <- otu_table / rowSums(otu_table)  

otu_rel_abund[is.na(otu_rel_abund)] <- 0

# 提取 Fusobacterium 的丰度
library(dplyr)
fusobacterium_percentage <- otu_rel_abund %>% select(g__Fusobacterium)
# 提取 Capnocytophaga 的丰度
Capnocytophaga_percentage <- otu_rel_abund %>% select(g__Capnocytophaga)
# 提取 Streptococcus 的丰度
Streptococcus_percentage <- otu_rel_abund %>% select(g__Streptococcus)
data <- data.frame(
  Sample = rownames(otu_rel_abund),
  Fusobacterium_percentage = unlist(fusobacterium_percentage),
  Capnocytophaga_percentage=unlist(Capnocytophaga_percentage),
  Streptococcus_percentage =unlist(Streptococcus_percentage)
  
)

rownames(data)<- data$Sample
UVM_data<- data[,-1]

head(UVM_data)


####合并数据-----
library(dplyr)

# 加上来源标签（癌种）
ACC_data$CancerType<-"ACC"
BLCA_data$CancerType<-"BLCA"
BRCA_data$CancerType<-"BRCA"
CESC_data$CancerType<-"CESC"
CHOL_data$CancerType<-"CHOL"
COAD_data$CancerType<-"COAD"
DLBC_data$CancerType<-"DLBC"
ESCA_data$CancerType<-"ESCA"
GBM_data$CancerType<-"GBM"
HNSC_data$CancerType<-"HNSC"
KICH_data$CancerType<-"KICH"
KIRC_data$CancerType<-"KIRC"
KIRP_data$CancerType<-"KIRP"
LAML_data$CancerType<-"LAML"
LGG_data$CancerType<-"LGG"
LIHC_data$CancerType<-"LIHC"
LUAD_data$CancerType<-"LUAD"
LUSC_data$CancerType<-"LUSC"
MESO_data$CancerType<-"MESO"
OV_data$CancerType<-"OV"
PAAD_data$CancerType<-"PAAD"
PCPG_data$CancerType<-"PCPG"
PRAD_data$CancerType<-"PRAD"
READ_data$CancerType<-"READ"
SARC_data$CancerType<-"SARC"
SKCM_data$CancerType<-"SKCM"
STAD_data$CancerType<-"STAD"
TGCT_data$CancerType<-"TGCT"
THCA_data$CancerType<-"THCA"
THYM_data$CancerType<-"THYM"
UCEC_data$CancerType<-"UCEC"
UCS_data$CancerType<-"UCS"
UVM_data$CancerType<-"UVM"

# 然后依次合并
merged_data <- bind_rows(
  ACC_data,
  BLCA_data,
  BRCA_data,
  CESC_data,
  CHOL_data,
  COAD_data,
  DLBC_data,
  ESCA_data,
  GBM_data,
  HNSC_data,
  KICH_data,
  KIRC_data,
  KIRP_data,
  LAML_data,
  LGG_data,
  LIHC_data,
  LUAD_data,
  LUSC_data,
  MESO_data,
  OV_data,
  PAAD_data,
  PCPG_data,
  PRAD_data,
  READ_data,
  SARC_data,
  SKCM_data,
  STAD_data,
  TGCT_data,
  THCA_data,
  THYM_data,
  UCEC_data,
  UCS_data,
  UVM_data)

merged_data

write.csv(merged_data,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Pan-cancer atlas of tumor-resident microbiome, immunity and prognosis/TCGA_microbe_合并.csv")


####TCGA预后分析----

microbe_data<-read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Pan-cancer atlas of tumor-resident microbiome, immunity and prognosis/TCGA_microbe_合并.csv",header = T,row.names = 1)
head(microbe_data)

clin_data<-read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Pan-cancer atlas of tumor-resident microbiome, immunity and prognosis/TCGA_clin_data_合并.csv",header = T,row.names = 1)
head(clin_data)

merged_data <- merge(
  clin_data,
  microbe_data,
  by = "submitter_id",
  all.x = TRUE
)

head(merged_data)


merged_data <- merge(
  clin_data,
  microbe_data,
  by = "submitter_id"
)

head(merged_data)

write.csv(merged_data,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Pan-cancer atlas of tumor-resident microbiome, immunity and prognosis/TCGA_clin_microbe_合并.csv")

merged_data$survival_time <- ifelse(
  !is.na(merged_data$days_to_death),
  merged_data$days_to_death,
  merged_data$days_to_last_follow_up
)

merged_data$event_status <- ifelse(
  merged_data$vital_status == "Dead", 1, 0
)




###########################################

project_ids <- c("TCGA-ACC",
                 "TCGA-BLCA", "TCGA-BRCA", "TCGA-CESC", "TCGA-CHOL", "TCGA-COAD", "TCGA-DLBC", 
                 "TCGA-ESCA", "TCGA-GBM", "TCGA-HNSC", "TCGA-KICH", "TCGA-KIRC", "TCGA-KIRP", 
                 "TCGA-LAML", "TCGA-LGG", "TCGA-LIHC", "TCGA-LUAD", "TCGA-LUSC", "TCGA-MESO", 
                 "TCGA-OV", "TCGA-PAAD", "TCGA-PCPG", "TCGA-PRAD", "TCGA-READ", "TCGA-SARC", 
                 "TCGA-SKCM", "TCGA-STAD", "TCGA-TGCT", "TCGA-THCA", "TCGA-THYM", "TCGA-UCEC", 
                 "TCGA-UCS", "TCGA-UVM"
)

####所有癌症整体分析----

merged_data2<- merged_data
library(survminer)

colnames(merged_data2)

sur.cut <- surv_cutpoint(merged_data2, time = "survival_time",
                         event = "event_status",
                         variables = "Fusobacterium_percentage")
summary(sur.cut)

#使用surv_categorize函数，根据截断值转为二分类变量。

sur.cat <- surv_categorize(sur.cut, labels = c("High", "Low"))


merged_data2$Fusobacterium_cut <-
  ifelse(merged_data2$Fusobacterium_percentage > 0, "High", "Low")

####百分之75的比例分组----

threshold <- quantile(merged_data2$Fusobacterium_percentage, 0.75,na.rm = TRUE)

merged_data2$Fusobacterium_cut <- ifelse(merged_data2$Fusobacterium_percentage <= threshold, "Low", "High")

unique(merged_data2$Fusobacterium_cut)

fit <- survfit(Surv(survival_time, event_status) ~Fusobacterium_cut,  # 创建生存对象 
               data = merged_data2)

pdf("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Pan-cancer atlas of tumor-resident microbiome, immunity and prognosis/预后分析画图/Fusobacterium_COAD_curvep_cut.pdf", width = 6, height =7)
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
           title = "Fusobacterium cancer")# 设置x轴刻度间距
dev.off()



####fusobacterium TCGA-ACC ----

merged_data2<- merged_data %>% filter(project %in% "TCGA-ACC")
library(survminer)

colnames(merged_data2)

sur.cut <- surv_cutpoint(merged_data2, time = "survival_time",
                         event = "event_status",
                         variables = "Fusobacterium_percentage")
summary(sur.cut)

#使用surv_categorize函数，根据截断值转为二分类变量。

sur.cat <- surv_categorize(sur.cut, labels = c("High", "Low"))


merged_data2$Fusobacterium_cut <-
  ifelse(merged_data2$Fusobacterium_percentage > 0.001105118, "High", "Low")

####百分之75的比例分组


threshold <- quantile(merged_data2$Fusobacterium_percentage, 0.75,na.rm = TRUE)

merged_data2$Fusobacterium_cut <- ifelse(merged_data2$Fusobacterium_percentage <= threshold, "Low", "High")

unique(merged_data2$Fusobacterium_cut)

fit <- survfit(Surv(survival_time, event_status) ~Fusobacterium_cut,  # 创建生存对象 
               data = merged_data2)

pdf("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Pan-cancer atlas of tumor-resident microbiome, immunity and prognosis/预后分析画图/Fusobacterium_ACC_curvep_cut.pdf", width = 6, height =7)
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
           title = "Fusobacterium ACC" )  # 设置x轴刻度间距
dev.off()


sur.cut <- surv_cutpoint(merged_data2, time = "survival_time",
                         event = "event_status",
                         variables = "Streptococcus_percentage")
summary(sur.cut)


#使用surv_categorize函数，根据截断值转为二分类变量。
sur.cat <- surv_categorize(sur.cut, labels = c("High", "Low"))


merged_data2$Streptococcus_cut <-
  ifelse(merged_data2$Streptococcus_percentage > 0.002925963, "High", "Low")

threshold <- quantile(merged_data2$Streptococcus_percentage, 0.75,na.rm = TRUE)

merged_data2$Streptococcus_cut <- ifelse(merged_data2$Streptococcus_percentage <= threshold, "Low", "High")

fit <- survfit(Surv(survival_time, event_status) ~Streptococcus_cut,  # 创建生存对象 
               data = merged_data2)

pdf("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Pan-cancer atlas of tumor-resident microbiome, immunity and prognosis/预后分析画图/Streptococcus_ACC_curvep_cut.pdf", width = 6, height =7)
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
           title = "Streptococcus ACC" )  # 设置x轴刻度间距
dev.off()



sur.cut <- surv_cutpoint(merged_data2, time = "survival_time",
                         event = "event_status",
                         variables = "Capnocytophaga_percentage")
summary(sur.cut)


#使用surv_categorize函数，根据截断值转为二分类变量。
sur.cat <- surv_categorize(sur.cut, labels = c("High", "Low"))

cap.cut <- median(merged_data2$Capnocytophaga_percentage, na.rm = TRUE)

merged_data2$Capnocytophaga_cut <-
  ifelse(merged_data2$Capnocytophaga_percentage > 0, "High", "Low")


fit <- survfit(Surv(survival_time, event_status) ~Capnocytophaga_cut,  # 创建生存对象 
               data = merged_data2)

pdf("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Pan-cancer atlas of tumor-resident microbiome, immunity and prognosis/预后分析画图/Capnocytophaga_ACC_curvep_cut.pdf", width = 6, height =7)
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
           title = "Capnocytophaga ACC" )  # 设置x轴刻度间距
dev.off()




merged_data2$Fusobacterium_cut <- factor(
  merged_data2$Fusobacterium_cut,
  levels = c("Low", "High") # 包含所有可能值
)

merged_data2$gender.x <- factor(merged_data2$gender.x)
sum(is.na(merged_data2$gender.x))

names(merged_data2)
unique(merged_data2$years_smoked)
merged_data2$years_smoked_history<-ifelse(merged_data2$years_smoked <= 20, "Low", "High")
merged_data2$portion_weight2 <- log(merged_data2$portion_weight + 1)
merged_data2$years_smoked2 <- log(merged_data2$years_smoked + 1)
merged_data2$years_smoked2[is.na(merged_data2$years_smoked2)] <- mean(merged_data2$years_smoked2, na.rm = TRUE)
merged_data2$Fusobacterium_cut <- relevel(merged_data2$Fusobacterium_cut, ref = "Low")

merged_data2$age_at_index
merged_data2$alcohol_history


cox_model <- coxph(Surv(survival_time, event_status) ~ Fusobacterium_cut+age_at_index+
                     prior_malignancy+
                     portion_weight2+gender.x, data = merged_data2)

summary(cox_model)
f1_Fusobacterium<-ggforest(cox_model,
                           data=merged_data2,
                           main = "Hazard ratio",        # 标题
                           cpositions = c(0.01, 0.1, 0.3), #左边三列的相对距离
                           fontsize = 1, # 字体大小
                           refLabel = "reference", #显示因子的参考水平
                           noDigits = 3)  #小数点位数

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Apan-cancer_analysis_of_the_microbiome/TCGA-BLCA/fusobacterium_BLCA_ggforest_cut_多因素.pdf",f1_Fusobacterium,width =13,height=6)



####fusobacterium TCGA-BLCA ----

merged_data2<- merged_data %>% filter(project %in% "TCGA-BLCA")
library(survminer)

colnames(merged_data2)

sur.cut <- surv_cutpoint(merged_data2, time = "survival_time",
                         event = "event_status",
                         variables = "Fusobacterium_percentage")
summary(sur.cut)

#使用surv_categorize函数，根据截断值转为二分类变量。
sur.cat <- surv_categorize(sur.cut, labels = c("High", "Low"))


merged_data2$Fusobacterium_cut <-
  ifelse(merged_data2$Fusobacterium_percentage > 0, "High", "Low")


threshold <- quantile(merged_data2$Fusobacterium_percentage, 0.75,na.rm = TRUE)

merged_data2$Fusobacterium_cut <- ifelse(merged_data2$Fusobacterium_percentage <=threshold, "Low", "High")


fit <- survfit(Surv(survival_time, event_status) ~Fusobacterium_cut,  # 创建生存对象 
               data = merged_data2)

pdf("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Pan-cancer atlas of tumor-resident microbiome, immunity and prognosis/预后分析画图/Fusobacterium_BLCA_curvep_cut.pdf", width = 6, height =7)
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
           title = "Fusobacterium BLCA" )  # 设置x轴刻度间距
dev.off()

sur.cut <- surv_cutpoint(merged_data2, time = "survival_time",
                         event = "event_status",
                         variables = "Streptococcus_percentage")
summary(sur.cut)


#使用surv_categorize函数，根据截断值转为二分类变量。
sur.cat <- surv_categorize(sur.cut, labels = c("High", "Low"))


merged_data2$Streptococcus_cut <-ifelse(merged_data2$Streptococcus_percentage >0, "High", "Low")

threshold <- quantile(merged_data2$Streptococcus_percentage, 0.75,na.rm = TRUE)

merged_data2$Streptococcus_cut <- ifelse(merged_data2$Streptococcus_percentage <=threshold, "Low", "High")


fit <- survfit(Surv(survival_time, event_status) ~Streptococcus_cut,  # 创建生存对象 
               data = merged_data2)

pdf("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Pan-cancer atlas of tumor-resident microbiome, immunity and prognosis/预后分析画图/Streptococcus_BLCA_curvep_cut.pdf", width = 6, height =7)
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
           title = "Streptococcus BLCA" )  # 设置x轴刻度间距
dev.off()


sur.cut <- surv_cutpoint(merged_data2, time = "survival_time",
                         event = "event_status",
                         variables = "Capnocytophaga_percentage")
summary(sur.cut)


#使用surv_categorize函数，根据截断值转为二分类变量。
sur.cat <- surv_categorize(sur.cut, labels = c("High", "Low"))

cap.cut <- median(merged_data2$Capnocytophaga_percentage, na.rm = TRUE)

merged_data2$Capno_group <- ifelse(
  merged_data2$Capnocytophaga_percentage > cap.cut,
  "High", "Low"
)

table(merged_data2$Capno_group)

merged_data2$Capnocytophaga_cut <-
  ifelse(merged_data2$Capnocytophaga_percentage > 0, "High", "Low")

unique(merged_data2$Capnocytophaga_cut)

fit <- survfit(Surv(survival_time, event_status) ~Capnocytophaga_cut,  # 创建生存对象 
               data = merged_data2)

pdf("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Pan-cancer atlas of tumor-resident microbiome, immunity and prognosis/预后分析画图/Capnocytophaga_BLCA_curvep_cut.pdf", width = 6, height =7)
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
           title = "Capnocytophaga BLCA" )  # 设置x轴刻度间距
dev.off()




####fusobacterium TCGA-BRCA ----

merged_data2<- merged_data %>% filter(project %in% "TCGA-BRCA")
library(survminer)

colnames(merged_data2)

sur.cut <- surv_cutpoint(merged_data2, time = "survival_time",
                         event = "event_status",
                         variables = "Fusobacterium_percentage")
summary(sur.cut)

#使用surv_categorize函数，根据截断值转为二分类变量。
sur.cat <- surv_categorize(sur.cut, labels = c("High", "Low"))


merged_data2$Fusobacterium_cut <-
  ifelse(merged_data2$Fusobacterium_percentage > 0, "High", "Low")


threshold <- quantile(merged_data2$Fusobacterium_percentage, 0.75,na.rm = TRUE)

merged_data2$Fusobacterium_cut <- ifelse(merged_data2$Fusobacterium_percentage <= threshold, "Low", "High")

unique(merged_data2$Fusobacterium_cut)

fit <- survfit(Surv(survival_time, event_status) ~Fusobacterium_cut,  # 创建生存对象 
               data = merged_data2)

pdf("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Pan-cancer atlas of tumor-resident microbiome, immunity and prognosis/预后分析画图/Fusobacterium_BRCA_curvep_cut.pdf", width = 6, height =7)
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
           xlim = c(0, 4000),
           break.x.by = 1000,
           title = "Fusobacterium BRCA" )  # 设置x轴刻度间距
dev.off()

sur.cut <- surv_cutpoint(merged_data2, time = "survival_time",
                         event = "event_status",
                         variables = "Streptococcus_percentage")
summary(sur.cut)


#使用surv_categorize函数，根据截断值转为二分类变量。
sur.cat <- surv_categorize(sur.cut, labels = c("High", "Low"))


merged_data2$Streptococcus_cut <-ifelse(merged_data2$Streptococcus_percentage >0.001069213, "High", "Low")


fit <- survfit(Surv(survival_time, event_status) ~Streptococcus_cut,  # 创建生存对象 
               data = merged_data2)

pdf("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Pan-cancer atlas of tumor-resident microbiome, immunity and prognosis/预后分析画图/Streptococcus_BRCA_curvep_cut.pdf", width = 6, height =7)
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
           title = "Streptococcus BRCA" )  # 设置x轴刻度间距
dev.off()


sur.cut <- surv_cutpoint(merged_data2, time = "survival_time",
                         event = "event_status",
                         variables = "Capnocytophaga_percentage")
summary(sur.cut)


#使用surv_categorize函数，根据截断值转为二分类变量。
sur.cat <- surv_categorize(sur.cut, labels = c("High", "Low"))

cap.cut <- median(merged_data2$Capnocytophaga_percentage, na.rm = TRUE)

merged_data2$Capno_group <- ifelse(
  merged_data2$Capnocytophaga_percentage > cap.cut,
  "High", "Low"
)

table(merged_data2$Capno_group)

merged_data2$Capnocytophaga_cut <-
  ifelse(merged_data2$Capnocytophaga_percentage > 0, "High", "Low")

unique(merged_data2$Capnocytophaga_cut)

fit <- survfit(Surv(survival_time, event_status) ~Capnocytophaga_cut,  # 创建生存对象 
               data = merged_data2)

pdf("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Pan-cancer atlas of tumor-resident microbiome, immunity and prognosis/预后分析画图/Capnocytophaga_BRCA_curvep_cut.pdf", width = 6, height =7)
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
           title = "Capnocytophaga BRCA" )  # 设置x轴刻度间距
dev.off()


####fusobacterium TCGA-CESC ----

merged_data2<- merged_data %>% filter(project %in% "TCGA-CESC")
library(survminer)

colnames(merged_data2)

sur.cut <- surv_cutpoint(merged_data2, time = "survival_time",
                         event = "event_status",
                         variables = "Fusobacterium_percentage")
summary(sur.cut)

#使用surv_categorize函数，根据截断值转为二分类变量。
sur.cat <- surv_categorize(sur.cut, labels = c("High", "Low"))


merged_data2$Fusobacterium_cut <-
  ifelse(merged_data2$Fusobacterium_percentage > 0, "High", "Low")


fit <- survfit(Surv(survival_time, event_status) ~Fusobacterium_cut,  # 创建生存对象 
               data = merged_data2)

pdf("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Pan-cancer atlas of tumor-resident microbiome, immunity and prognosis/预后分析画图/Fusobacterium_CESC_curvep_cut.pdf", width = 6, height =7)
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
           xlim = c(0, 4000),
           break.x.by = 1000,
           title = "Fusobacterium CESC" )  # 设置x轴刻度间距
dev.off()

sur.cut <- surv_cutpoint(merged_data2, time = "survival_time",
                         event = "event_status",
                         variables = "Streptococcus_percentage")
summary(sur.cut)


#使用surv_categorize函数，根据截断值转为二分类变量。
sur.cat <- surv_categorize(sur.cut, labels = c("High", "Low"))


merged_data2$Streptococcus_cut <-ifelse(merged_data2$Streptococcus_percentage >0, "High", "Low")


fit <- survfit(Surv(survival_time, event_status) ~Streptococcus_cut,  # 创建生存对象 
               data = merged_data2)

pdf("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Pan-cancer atlas of tumor-resident microbiome, immunity and prognosis/预后分析画图/Streptococcus_CESC_curvep_cut.pdf", width = 6, height =7)
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
           title = "Streptococcus CESC" )  # 设置x轴刻度间距
dev.off()


sur.cut <- surv_cutpoint(merged_data2, time = "survival_time",
                         event = "event_status",
                         variables = "Capnocytophaga_percentage")
summary(sur.cut)


#使用surv_categorize函数，根据截断值转为二分类变量。
sur.cat <- surv_categorize(sur.cut, labels = c("High", "Low"))

cap.cut <- median(merged_data2$Capnocytophaga_percentage, na.rm = TRUE)

merged_data2$Capno_group <- ifelse(
  merged_data2$Capnocytophaga_percentage > cap.cut,
  "High", "Low"
)

table(merged_data2$Capno_group)

merged_data2$Capnocytophaga_cut <-
  ifelse(merged_data2$Capnocytophaga_percentage >0, "High", "Low")

unique(merged_data2$Capnocytophaga_cut)

fit <- survfit(Surv(survival_time, event_status) ~Capnocytophaga_cut,  # 创建生存对象 
               data = merged_data2)

pdf("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Pan-cancer atlas of tumor-resident microbiome, immunity and prognosis/预后分析画图/Capnocytophaga_CESC_curvep_cut.pdf", width = 6, height =7)
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
           title = "Capnocytophaga CESC" )  # 设置x轴刻度间距
dev.off()





####fusobacterium TCGA-COAD ----

merged_data2<- merged_data %>% filter(project %in% "TCGA-COAD")
library(survminer)

colnames(merged_data2)

sur.cut <- surv_cutpoint(merged_data2, time = "survival_time",
                         event = "event_status",
                         variables = "Fusobacterium_percentage")
summary(sur.cut)

#使用surv_categorize函数，根据截断值转为二分类变量。

sur.cat <- surv_categorize(sur.cut, labels = c("High", "Low"))


merged_data2$Fusobacterium_cut <-
  ifelse(merged_data2$Fusobacterium_percentage > 0.001105118, "High", "Low")

####百分之75的比例分组

threshold <- quantile(merged_data2$Fusobacterium_percentage, 0.75,na.rm = TRUE)

merged_data2$Fusobacterium_cut <- ifelse(merged_data2$Fusobacterium_percentage <= threshold, "Low", "High")

unique(merged_data2$Fusobacterium_cut)

fit <- survfit(Surv(survival_time, event_status) ~Fusobacterium_cut,  # 创建生存对象 
               data = merged_data2)

pdf("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Pan-cancer atlas of tumor-resident microbiome, immunity and prognosis/预后分析画图/Fusobacterium_COAD_curvep_cut.pdf", width = 6, height =7)
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
           title = "Fusobacterium COAD" )  # 设置x轴刻度间距
dev.off()


sur.cut <- surv_cutpoint(merged_data2, time = "survival_time",
                         event = "event_status",
                         variables = "Streptococcus_percentage")
summary(sur.cut)


#使用surv_categorize函数，根据截断值转为二分类变量。
sur.cat <- surv_categorize(sur.cut, labels = c("High", "Low"))


merged_data2$Streptococcus_cut <-
  ifelse(merged_data2$Streptococcus_percentage > 0.002925963, "High", "Low")

threshold <- quantile(merged_data2$Streptococcus_percentage, 0.75,na.rm = TRUE)

merged_data2$Streptococcus_cut <- ifelse(merged_data2$Streptococcus_percentage <= threshold, "Low", "High")

fit <- survfit(Surv(survival_time, event_status) ~Streptococcus_cut,  # 创建生存对象 
               data = merged_data2)

pdf("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Pan-cancer atlas of tumor-resident microbiome, immunity and prognosis/预后分析画图/Streptococcus_COAD_curvep_cut.pdf", width = 6, height =7)
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
           title = "Streptococcus COAD" )  # 设置x轴刻度间距
dev.off()



sur.cut <- surv_cutpoint(merged_data2, time = "survival_time",
                         event = "event_status",
                         variables = "Capnocytophaga_percentage")
summary(sur.cut)


#使用surv_categorize函数，根据截断值转为二分类变量。
sur.cat <- surv_categorize(sur.cut, labels = c("High", "Low"))

cap.cut <- median(merged_data2$Capnocytophaga_percentage, na.rm = TRUE)

merged_data2$Capnocytophaga_cut <-
  ifelse(merged_data2$Capnocytophaga_percentage > 0, "High", "Low")


fit <- survfit(Surv(survival_time, event_status) ~Capnocytophaga_cut,  # 创建生存对象 
               data = merged_data2)

pdf("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Pan-cancer atlas of tumor-resident microbiome, immunity and prognosis/预后分析画图/Capnocytophaga_COAD_curvep_cut.pdf", width = 6, height =7)
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
           title = "Capnocytophaga COAD" )  # 设置x轴刻度间距
dev.off()




merged_data2$Fusobacterium_cut <- factor(
  merged_data2$Fusobacterium_cut,
  levels = c("Low", "High") # 包含所有可能值
)

merged_data2$gender.x <- factor(merged_data2$gender.x)
sum(is.na(merged_data2$gender.x))

names(merged_data2)
unique(merged_data2$years_smoked)
merged_data2$years_smoked_history<-ifelse(merged_data2$years_smoked <= 20, "Low", "High")
merged_data2$portion_weight2 <- log(merged_data2$portion_weight + 1)
merged_data2$years_smoked2 <- log(merged_data2$years_smoked + 1)
merged_data2$years_smoked2[is.na(merged_data2$years_smoked2)] <- mean(merged_data2$years_smoked2, na.rm = TRUE)
merged_data2$Fusobacterium_cut <- relevel(merged_data2$Fusobacterium_cut, ref = "Low")

merged_data2$age_at_index
merged_data2$alcohol_history


cox_model <- coxph(Surv(survival_time, event_status) ~ Fusobacterium_cut+age_at_index+
                     prior_malignancy+
                     portion_weight2+gender.x, data = merged_data2)

summary(cox_model)
f1_Fusobacterium<-ggforest(cox_model,
                           data=merged_data2,
                           main = "Hazard ratio",        # 标题
                           cpositions = c(0.01, 0.1, 0.3), #左边三列的相对距离
                           fontsize = 1, # 字体大小
                           refLabel = "reference", #显示因子的参考水平
                           noDigits = 3)  #小数点位数

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Apan-cancer_analysis_of_the_microbiome/TCGA-BLCA/fusobacterium_BLCA_ggforest_cut_多因素.pdf",f1_Fusobacterium,width =13,height=6)





####fusobacterium TCGA-CHOL ----

merged_data2<- merged_data %>% filter(project %in% "TCGA-CHOL")
library(survminer)

colnames(merged_data2)

sur.cut <- surv_cutpoint(merged_data2, time = "survival_time",
                         event = "event_status",
                         variables = "Fusobacterium_percentage")
summary(sur.cut)

#使用surv_categorize函数，根据截断值转为二分类变量。

sur.cat <- surv_categorize(sur.cut, labels = c("High", "Low"))


merged_data2$Fusobacterium_cut <-
  ifelse(merged_data2$Fusobacterium_percentage > 0.001105118, "High", "Low")

####百分之75的比例分组

threshold <- quantile(merged_data2$Fusobacterium_percentage, 0.75,na.rm = TRUE)

merged_data2$Fusobacterium_cut <- ifelse(merged_data2$Fusobacterium_percentage <= threshold, "Low", "High")

unique(merged_data2$Fusobacterium_cut)

fit <- survfit(Surv(survival_time, event_status) ~Fusobacterium_cut,  # 创建生存对象 
               data = merged_data2)

pdf("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Pan-cancer atlas of tumor-resident microbiome, immunity and prognosis/预后分析画图/Fusobacterium_CHOL_curvep_cut.pdf", width = 6, height =7)
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
           title = "Fusobacterium CHOL" )  # 设置x轴刻度间距
dev.off()


sur.cut <- surv_cutpoint(merged_data2, time = "survival_time",
                         event = "event_status",
                         variables = "Streptococcus_percentage")
summary(sur.cut)


#使用surv_categorize函数，根据截断值转为二分类变量。
sur.cat <- surv_categorize(sur.cut, labels = c("High", "Low"))


merged_data2$Streptococcus_cut <-
  ifelse(merged_data2$Streptococcus_percentage > 0.002925963, "High", "Low")

threshold <- quantile(merged_data2$Streptococcus_percentage, 0.75,na.rm = TRUE)

merged_data2$Streptococcus_cut <- ifelse(merged_data2$Streptococcus_percentage <= threshold, "Low", "High")

fit <- survfit(Surv(survival_time, event_status) ~Streptococcus_cut,  # 创建生存对象 
               data = merged_data2)

pdf("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Pan-cancer atlas of tumor-resident microbiome, immunity and prognosis/预后分析画图/Streptococcus_CHOL_curvep_cut.pdf", width = 6, height =7)
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
           title = "Streptococcus CHOL" )  # 设置x轴刻度间距
dev.off()



sur.cut <- surv_cutpoint(merged_data2, time = "survival_time",
                         event = "event_status",
                         variables = "Capnocytophaga_percentage")
summary(sur.cut)


#使用surv_categorize函数，根据截断值转为二分类变量。
sur.cat <- surv_categorize(sur.cut, labels = c("High", "Low"))

cap.cut <- median(merged_data2$Capnocytophaga_percentage, na.rm = TRUE)

merged_data2$Capnocytophaga_cut <-
  ifelse(merged_data2$Capnocytophaga_percentage > 0, "High", "Low")


fit <- survfit(Surv(survival_time, event_status) ~Capnocytophaga_cut,  # 创建生存对象 
               data = merged_data2)

pdf("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Pan-cancer atlas of tumor-resident microbiome, immunity and prognosis/预后分析画图/Capnocytophaga_CHOL_curvep_cut.pdf", width = 6, height =7)
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
           title = "Capnocytophaga CHOL" )  # 设置x轴刻度间距
dev.off()




merged_data2$Fusobacterium_cut <- factor(
  merged_data2$Fusobacterium_cut,
  levels = c("Low", "High") # 包含所有可能值
)

merged_data2$gender.x <- factor(merged_data2$gender.x)
sum(is.na(merged_data2$gender.x))

names(merged_data2)
unique(merged_data2$years_smoked)
merged_data2$years_smoked_history<-ifelse(merged_data2$years_smoked <= 20, "Low", "High")
merged_data2$portion_weight2 <- log(merged_data2$portion_weight + 1)
merged_data2$years_smoked2 <- log(merged_data2$years_smoked + 1)
merged_data2$years_smoked2[is.na(merged_data2$years_smoked2)] <- mean(merged_data2$years_smoked2, na.rm = TRUE)
merged_data2$Fusobacterium_cut <- relevel(merged_data2$Fusobacterium_cut, ref = "Low")

merged_data2$age_at_index
merged_data2$alcohol_history


cox_model <- coxph(Surv(survival_time, event_status) ~ Fusobacterium_cut+age_at_index+
                     prior_malignancy+
                     portion_weight2+gender.x, data = merged_data2)

summary(cox_model)
f1_Fusobacterium<-ggforest(cox_model,
                           data=merged_data2,
                           main = "Hazard ratio",        # 标题
                           cpositions = c(0.01, 0.1, 0.3), #左边三列的相对距离
                           fontsize = 1, # 字体大小
                           refLabel = "reference", #显示因子的参考水平
                           noDigits = 3)  #小数点位数

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Apan-cancer_analysis_of_the_microbiome/TCGA-BLCA/fusobacterium_BLCA_ggforest_cut_多因素.pdf",f1_Fusobacterium,width =13,height=6)







####fusobacterium TCGA-DLBC ----

merged_data2<- merged_data %>% filter(project %in% "TCGA-DLBC")
library(survminer)

colnames(merged_data2)

sur.cut <- surv_cutpoint(merged_data2, time = "survival_time",
                         event = "event_status",
                         variables = "Fusobacterium_percentage")
summary(sur.cut)

#使用surv_categorize函数，根据截断值转为二分类变量。

sur.cat <- surv_categorize(sur.cut, labels = c("High", "Low"))


merged_data2$Fusobacterium_cut <-
  ifelse(merged_data2$Fusobacterium_percentage > 0.001105118, "High", "Low")

####百分之75的比例分组

threshold <- quantile(merged_data2$Fusobacterium_percentage, 0.75,na.rm = TRUE)

merged_data2$Fusobacterium_cut <- ifelse(merged_data2$Fusobacterium_percentage <= threshold, "Low", "High")

unique(merged_data2$Fusobacterium_cut)

fit <- survfit(Surv(survival_time, event_status) ~Fusobacterium_cut,  # 创建生存对象 
               data = merged_data2)

pdf("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Pan-cancer atlas of tumor-resident microbiome, immunity and prognosis/预后分析画图/Fusobacterium_DLBC_curvep_cut.pdf", width = 6, height =7)
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
           title = "Fusobacterium DLBC" )  # 设置x轴刻度间距
dev.off()


sur.cut <- surv_cutpoint(merged_data2, time = "survival_time",
                         event = "event_status",
                         variables = "Streptococcus_percentage")
summary(sur.cut)


#使用surv_categorize函数，根据截断值转为二分类变量。
sur.cat <- surv_categorize(sur.cut, labels = c("High", "Low"))


merged_data2$Streptococcus_cut <-
  ifelse(merged_data2$Streptococcus_percentage > 0.002925963, "High", "Low")

threshold <- quantile(merged_data2$Streptococcus_percentage, 0.75,na.rm = TRUE)

merged_data2$Streptococcus_cut <- ifelse(merged_data2$Streptococcus_percentage <= threshold, "Low", "High")

fit <- survfit(Surv(survival_time, event_status) ~Streptococcus_cut,  # 创建生存对象 
               data = merged_data2)

pdf("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Pan-cancer atlas of tumor-resident microbiome, immunity and prognosis/预后分析画图/Streptococcus_DLBC_curvep_cut.pdf", width = 6, height =7)
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
           title = "Streptococcus DLBC" )  # 设置x轴刻度间距
dev.off()



sur.cut <- surv_cutpoint(merged_data2, time = "survival_time",
                         event = "event_status",
                         variables = "Capnocytophaga_percentage")
summary(sur.cut)


#使用surv_categorize函数，根据截断值转为二分类变量。
sur.cat <- surv_categorize(sur.cut, labels = c("High", "Low"))

cap.cut <- median(merged_data2$Capnocytophaga_percentage, na.rm = TRUE)

merged_data2$Capnocytophaga_cut <-
  ifelse(merged_data2$Capnocytophaga_percentage > 0, "High", "Low")


fit <- survfit(Surv(survival_time, event_status) ~Capnocytophaga_cut,  # 创建生存对象 
               data = merged_data2)

pdf("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Pan-cancer atlas of tumor-resident microbiome, immunity and prognosis/预后分析画图/Capnocytophaga_DLBC_curvep_cut.pdf", width = 6, height =7)
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
           title = "Capnocytophaga DLBC" )  # 设置x轴刻度间距
dev.off()





####fusobacterium TCGA-ESCA 高的生存率高 ----

merged_data2<- merged_data %>% filter(project %in% "TCGA-ESCA")
library(survminer)

colnames(merged_data2)

sur.cut <- surv_cutpoint(merged_data2, time = "survival_time",
                         event = "event_status",
                         variables = "Fusobacterium_percentage")
summary(sur.cut)

#使用surv_categorize函数，根据截断值转为二分类变量。
sur.cat <- surv_categorize(sur.cut, labels = c("High", "Low"))


merged_data2$Fusobacterium_cut <-
  ifelse(merged_data2$Fusobacterium_percentage > 3.56e-05, "High", "Low")


threshold <- quantile(merged_data2$Fusobacterium_percentage, 0.75,na.rm = TRUE)

merged_data2$Fusobacterium_cut <- ifelse(merged_data2$Fusobacterium_percentage <= threshold, "Low", "High")


fit <- survfit(Surv(survival_time, event_status) ~Fusobacterium_cut,  # 创建生存对象 
               data = merged_data2)

pdf("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Pan-cancer atlas of tumor-resident microbiome, immunity and prognosis/预后分析画图/Fusobacterium_ESCA_curvep_cut.pdf", width = 6, height =7)
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


sur.cut <- surv_cutpoint(merged_data2, time = "survival_time",
                         event = "event_status",
                         variables = "Streptococcus_percentage")
summary(sur.cut)


#使用surv_categorize函数，根据截断值转为二分类变量。
sur.cat <- surv_categorize(sur.cut, labels = c("High", "Low"))


merged_data2$Streptococcus_cut <-
  ifelse(merged_data2$Streptococcus_percentage > 0.000456361, "High", "Low")


fit <- survfit(Surv(survival_time, event_status) ~Streptococcus_cut,  # 创建生存对象 
               data = merged_data2)

pdf("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Pan-cancer atlas of tumor-resident microbiome, immunity and prognosis/预后分析画图/Streptococcus_ESCA_curvep_cut.pdf", width = 6, height =7)
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


sur.cut <- surv_cutpoint(merged_data2, time = "survival_time",
                         event = "event_status",
                         variables = "Capnocytophaga_percentage")
summary(sur.cut)


#使用surv_categorize函数，根据截断值转为二分类变量。
sur.cat <- surv_categorize(sur.cut, labels = c("High", "Low"))


merged_data2$Capnocytophaga_cut <-
  ifelse(merged_data2$Capnocytophaga_percentage > 0, "High", "Low")


fit <- survfit(Surv(survival_time, event_status) ~Capnocytophaga_cut,  # 创建生存对象 
               data = merged_data2)

pdf("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Pan-cancer atlas of tumor-resident microbiome, immunity and prognosis/预后分析画图/Capnocytophaga_ESCA_curvep_cut.pdf", width = 6, height =7)
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
           title = "Capnocytophaga ESCA" )  # 设置x轴刻度间距
dev.off()


####fusobacterium TCGA-GBM----

merged_data2<- merged_data %>% filter(project %in% "TCGA-GBM")
library(survminer)

colnames(merged_data2)

sur.cut <- surv_cutpoint(merged_data2, time = "survival_time",
                         event = "event_status",
                         variables = "Fusobacterium_percentage")
summary(sur.cut)

#使用surv_categorize函数，根据截断值转为二分类变量。
sur.cat <- surv_categorize(sur.cut, labels = c("High", "Low"))


merged_data2$Fusobacterium_cut <-
  ifelse(merged_data2$Fusobacterium_percentage > 3.56e-05, "High", "Low")


threshold <- quantile(merged_data2$Fusobacterium_percentage, 0.75,na.rm = TRUE)

merged_data2$Fusobacterium_cut <- ifelse(merged_data2$Fusobacterium_percentage <= threshold, "Low", "High")


fit <- survfit(Surv(survival_time, event_status) ~Fusobacterium_cut,  # 创建生存对象 
               data = merged_data2)

pdf("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Pan-cancer atlas of tumor-resident microbiome, immunity and prognosis/预后分析画图/Fusobacterium_GBM_curvep_cut.pdf", width = 6, height =7)
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
           title = "Fusobacterium GBM" )  # 设置x轴刻度间距
dev.off()


sur.cut <- surv_cutpoint(merged_data2, time = "survival_time",
                         event = "event_status",
                         variables = "Streptococcus_percentage")
summary(sur.cut)


#使用surv_categorize函数，根据截断值转为二分类变量。
sur.cat <- surv_categorize(sur.cut, labels = c("High", "Low"))


merged_data2$Streptococcus_cut <-
  ifelse(merged_data2$Streptococcus_percentage > 0.000456361, "High", "Low")

threshold <- quantile(merged_data2$Streptococcus_percentage, 0.75,na.rm = TRUE)

merged_data2$Streptococcus_cut <- ifelse(merged_data2$Streptococcus_percentage <= threshold, "Low", "High")


fit <- survfit(Surv(survival_time, event_status) ~Streptococcus_cut,  # 创建生存对象 
               data = merged_data2)

pdf("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Pan-cancer atlas of tumor-resident microbiome, immunity and prognosis/预后分析画图/Streptococcus_GBM_curvep_cut.pdf", width = 6, height =7)
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
           title = "Streptococcus GBM" )  # 设置x轴刻度间距
dev.off()


sur.cut <- surv_cutpoint(merged_data2, time = "survival_time",
                         event = "event_status",
                         variables = "Capnocytophaga_percentage")
summary(sur.cut)


#使用surv_categorize函数，根据截断值转为二分类变量。
sur.cat <- surv_categorize(sur.cut, labels = c("High", "Low"))


merged_data2$Capnocytophaga_cut <-
  ifelse(merged_data2$Capnocytophaga_percentage > 0, "High", "Low")


threshold <- quantile(merged_data2$Capnocytophaga_percentage, 0.75,na.rm = TRUE)

merged_data2$Capnocytophaga_cut <- ifelse(merged_data2$Capnocytophaga_percentage <= threshold, "Low", "High")

fit <- survfit(Surv(survival_time, event_status) ~Capnocytophaga_cut,  # 创建生存对象 
               data = merged_data2)

pdf("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Pan-cancer atlas of tumor-resident microbiome, immunity and prognosis/预后分析画图/Capnocytophaga_GBM_curvep_cut.pdf", width = 6, height =7)
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
           title = "Capnocytophaga GBM" )  # 设置x轴刻度间距
dev.off()


####fusobacterium TCGA-HNSC 没有显著性----

merged_data2<- merged_data %>% filter(project %in% "TCGA-HNSC")
library(survminer)

colnames(merged_data2)

sur.cut <- surv_cutpoint(merged_data2, time = "survival_time",
                         event = "event_status",
                         variables = "Fusobacterium_percentage")
summary(sur.cut)

#使用surv_categorize函数，根据截断值转为二分类变量。
sur.cat <- surv_categorize(sur.cut, labels = c("High", "Low"))


merged_data2$Fusobacterium_cut <-
  ifelse(merged_data2$Fusobacterium_percentage >= 0.006337136, "High", "Low")


threshold <- quantile(merged_data2$Fusobacterium_percentage, 0.75,na.rm = TRUE)

merged_data2$Fusobacterium_cut <- ifelse(merged_data2$Fusobacterium_percentage <= threshold, "Low", "High")


fit <- survfit(Surv(survival_time, event_status) ~Fusobacterium_cut,  # 创建生存对象 
               data = merged_data2)

pdf("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Pan-cancer atlas of tumor-resident microbiome, immunity and prognosis/预后分析画图/Fusobacterium_HNSC_curvep_cut.pdf", width = 6, height =7)
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


merged_data2$Fusobacterium_cut <- factor(
  merged_data2$Fusobacterium_cut,
  levels = c("Low", "High") # 包含所有可能值
)



merged_data2$gender.x <- factor(merged_data2$gender)
sum(is.na(merged_data2$gender))

names(merged_data2)
unique(merged_data2$years_smoked)

merged_data2$age_at_index <- merged_data2$age_at_index 
merged_data2$years_smoked_history<-ifelse(merged_data2$years_smoked <= 20, "Low", "High")
# merged_data2$portion_weight2 <- log(merged_data2$portion_weight + 1)
merged_data2$years_smoked2 <- log(merged_data2$years_smoked + 1)
merged_data2$years_smoked2[is.na(merged_data2$years_smoked2)] <- mean(merged_data2$years_smoked2, na.rm = TRUE)
merged_data2$alcohol_intensity2 <- log(merged_data2$alcohol_intensity + 1)
merged_data2$alcohol_intensity2[is.na(merged_data2$alcohol_intensity2)] <- mean(merged_data2$alcohol_intensity2, na.rm = TRUE)
merged_data2$Fusobacterium_cut <- relevel(merged_data2$Fusobacterium_cut, ref = "Low")

merged_data2$age_at_index
merged_data2$alcohol_history
merged_data2$prior_malignancy

cox_model <- coxph(Surv(survival_time, event_status) ~ Fusobacterium_cut+gender+
                     age_at_index+years_smoked_history+alcohol_intensity2
                     prior_malignancy, data = merged_data2)

colSums(is.na(merged_data2[, c("survival_time", "event_status", 
                               "Fusobacterium_cut", "gender", 
                               "age_at_index", "years_smoked_history", 
                               "alcohol_intensity2", "prior_malignancy")]))


cox_data <- na.omit(merged_data2[, c("survival_time", "event_status", 
                                     "Fusobacterium_cut", "gender", 
                                     "age_at_index", "prior_malignancy")])

merged_data2$Capnocytophaga_group <- relevel(merged_data2$Fusobacterium_cut, ref = "Low")

cox_model <- coxph(Surv(survival_time, event_status) ~ Fusobacterium_cut+gender+
                     age_at_index+
                   prior_malignancy, data = cox_data)


summary(cox_model)
f1_Fusobacterium<-ggforest(cox_model,
                           data=merged_data2,
                           main = "Hazard ratio",        # 标题
                           cpositions = c(0.01, 0.1, 0.3), #左边三列的相对距离
                           fontsize = 1, # 字体大小
                           refLabel = "reference", #显示因子的参考水平
                           noDigits = 3)  #小数点位数

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Apan-cancer_analysis_of_the_microbiome/TCGA-BLCA/fusobacterium_BLCA_ggforest_cut_多因素.pdf",f1_Fusobacterium,width =13,height=6)








sur.cut <- surv_cutpoint(merged_data2, time = "survival_time",
                         event = "event_status",
                         variables = "Streptococcus_percentage")
summary(sur.cut)


#使用surv_categorize函数，根据截断值转为二分类变量。
sur.cat <- surv_categorize(sur.cut, labels = c("High", "Low"))


merged_data2$Streptococcus_cut <-ifelse(merged_data2$Streptococcus_percentage >0, "High", "Low")

threshold <- quantile(merged_data2$Streptococcus_percentage, 0.75,na.rm = TRUE)

merged_data2$Streptococcus_cut <- ifelse(merged_data2$Streptococcus_percentage <=threshold, "Low", "High")


fit <- survfit(Surv(survival_time, event_status) ~Streptococcus_cut,  # 创建生存对象 
               data = merged_data2)

pdf("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Pan-cancer atlas of tumor-resident microbiome, immunity and prognosis/预后分析画图/Streptococcus_HNSC_curvep_cut.pdf", width = 6, height =7)
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


sur.cut <- surv_cutpoint(merged_data2, time = "survival_time",
                         event = "event_status",
                         variables = "Capnocytophaga_percentage")
summary(sur.cut)


#使用surv_categorize函数，根据截断值转为二分类变量。
sur.cat <- surv_categorize(sur.cut, labels = c("High", "Low"))

cap.cut <- median(merged_data2$Capnocytophaga_percentage, na.rm = TRUE)

merged_data2$Capno_group <- ifelse(
  merged_data2$Capnocytophaga_percentage > cap.cut,
  "High", "Low"
)

table(merged_data2$Capno_group)

merged_data2$Capnocytophaga_cut <-
  ifelse(merged_data2$Capnocytophaga_percentage > 0, "High", "Low")

unique(merged_data2$Capnocytophaga_cut)


threshold <- quantile(merged_data2$Capnocytophaga_percentage, 0.75,na.rm = TRUE)

merged_data2$Capnocytophaga_cut <- ifelse(merged_data2$Capnocytophaga_percentage <=threshold, "Low", "High")



fit <- survfit(Surv(survival_time, event_status) ~Capnocytophaga_cut,  # 创建生存对象 
               data = merged_data2)

pdf("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Pan-cancer atlas of tumor-resident microbiome, immunity and prognosis/预后分析画图/Capnocytophaga_HNSC_curvep_cut.pdf", width = 6, height =7)
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
           title = "Capnocytophaga HNSC" )  # 设置x轴刻度间距
dev.off()









####fusobacterium TCGA-KIRC 高丰度高生存率 ----

merged_data2<- merged_data %>% filter(project %in% "TCGA-KIRC")
library(survminer)

colnames(merged_data2)

sur.cut <- surv_cutpoint(merged_data2, time = "survival_time",
                         event = "event_status",
                         variables = "Fusobacterium_percentage")
summary(sur.cut)

#使用surv_categorize函数，根据截断值转为二分类变量。
sur.cat <- surv_categorize(sur.cut, labels = c("High", "Low"))


merged_data2$Fusobacterium_cut <-
  ifelse(merged_data2$Fusobacterium_percentage >= 0.006337136, "High", "Low")


threshold <- quantile(merged_data2$Fusobacterium_percentage, 0.75,na.rm = TRUE)

merged_data2$Fusobacterium_cut <- ifelse(merged_data2$Fusobacterium_percentage <= threshold, "Low", "High")


fit <- survfit(Surv(survival_time, event_status) ~Fusobacterium_cut,  # 创建生存对象 
               data = merged_data2)

pdf("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Pan-cancer atlas of tumor-resident microbiome, immunity and prognosis/预后分析画图/Fusobacterium_KIRC_curvep_cut.pdf", width = 6, height =7)
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
           title = "Fusobacterium KIRC" )  # 设置x轴刻度间距
dev.off()


sur.cut <- surv_cutpoint(merged_data2, time = "survival_time",
                         event = "event_status",
                         variables = "Streptococcus_percentage")
summary(sur.cut)


#使用surv_categorize函数，根据截断值转为二分类变量。
sur.cat <- surv_categorize(sur.cut, labels = c("High", "Low"))


merged_data2$Streptococcus_cut <-ifelse(merged_data2$Streptococcus_percentage >0, "High", "Low")

threshold <- quantile(merged_data2$Streptococcus_percentage, 0.75,na.rm = TRUE)

merged_data2$Streptococcus_cut <- ifelse(merged_data2$Streptococcus_percentage <=threshold, "Low", "High")


fit <- survfit(Surv(survival_time, event_status) ~Streptococcus_cut,  # 创建生存对象 
               data = merged_data2)

pdf("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Pan-cancer atlas of tumor-resident microbiome, immunity and prognosis/预后分析画图/Streptococcus_KIRC_curvep_cut.pdf", width = 6, height =7)
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
           title = "Streptococcus KIRC" )  # 设置x轴刻度间距
dev.off()


sur.cut <- surv_cutpoint(merged_data2, time = "survival_time",
                         event = "event_status",
                         variables = "Capnocytophaga_percentage")
summary(sur.cut)


#使用surv_categorize函数，根据截断值转为二分类变量。
sur.cat <- surv_categorize(sur.cut, labels = c("High", "Low"))

cap.cut <- median(merged_data2$Capnocytophaga_percentage, na.rm = TRUE)

merged_data2$Capno_group <- ifelse(
  merged_data2$Capnocytophaga_percentage > cap.cut,
  "High", "Low"
)

table(merged_data2$Capno_group)

merged_data2$Capnocytophaga_cut <-
  ifelse(merged_data2$Capnocytophaga_percentage > 0, "High", "Low")

unique(merged_data2$Capnocytophaga_cut)


threshold <- quantile(merged_data2$Capnocytophaga_percentage, 0.75,na.rm = TRUE)

merged_data2$Capnocytophaga_cut <- ifelse(merged_data2$Capnocytophaga_percentage <=threshold, "Low", "High")



fit <- survfit(Surv(survival_time, event_status) ~Capnocytophaga_cut,  # 创建生存对象 
               data = merged_data2)

pdf("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Pan-cancer atlas of tumor-resident microbiome, immunity and prognosis/预后分析画图/Capnocytophaga_KIRC_curvep_cut.pdf", width = 6, height =7)
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
           title = "Capnocytophaga KIRC" )  # 设置x轴刻度间距
dev.off()






####fusobacterium TCGA-KIRP 高丰度高生存率 ----

merged_data2<- merged_data %>% filter(project %in% "TCGA-KIRP")
library(survminer)

colnames(merged_data2)

sur.cut <- surv_cutpoint(merged_data2, time = "survival_time",
                         event = "event_status",
                         variables = "Fusobacterium_percentage")
summary(sur.cut)

#使用surv_categorize函数，根据截断值转为二分类变量。
sur.cat <- surv_categorize(sur.cut, labels = c("High", "Low"))


merged_data2$Fusobacterium_cut <-
  ifelse(merged_data2$Fusobacterium_percentage >= 0.006337136, "High", "Low")


threshold <- quantile(merged_data2$Fusobacterium_percentage, 0.75,na.rm = TRUE)

merged_data2$Fusobacterium_cut <- ifelse(merged_data2$Fusobacterium_percentage <= threshold, "Low", "High")


fit <- survfit(Surv(survival_time, event_status) ~Fusobacterium_cut,  # 创建生存对象 
               data = merged_data2)

pdf("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Pan-cancer atlas of tumor-resident microbiome, immunity and prognosis/预后分析画图/Fusobacterium_KIRP_curvep_cut.pdf", width = 6, height =7)
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
           title = "Fusobacterium KIRP" )  # 设置x轴刻度间距
dev.off()


sur.cut <- surv_cutpoint(merged_data2, time = "survival_time",
                         event = "event_status",
                         variables = "Streptococcus_percentage")
summary(sur.cut)


#使用surv_categorize函数，根据截断值转为二分类变量。
sur.cat <- surv_categorize(sur.cut, labels = c("High", "Low"))


merged_data2$Streptococcus_cut <-ifelse(merged_data2$Streptococcus_percentage >0, "High", "Low")

threshold <- quantile(merged_data2$Streptococcus_percentage, 0.75,na.rm = TRUE)

merged_data2$Streptococcus_cut <- ifelse(merged_data2$Streptococcus_percentage <=threshold, "Low", "High")


fit <- survfit(Surv(survival_time, event_status) ~Streptococcus_cut,  # 创建生存对象 
               data = merged_data2)

pdf("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Pan-cancer atlas of tumor-resident microbiome, immunity and prognosis/预后分析画图/Streptococcus_KIRP_curvep_cut.pdf", width = 6, height =7)
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
           title = "Streptococcus KIRP" )  # 设置x轴刻度间距
dev.off()


sur.cut <- surv_cutpoint(merged_data2, time = "survival_time",
                         event = "event_status",
                         variables = "Capnocytophaga_percentage")
summary(sur.cut)


#使用surv_categorize函数，根据截断值转为二分类变量。
sur.cat <- surv_categorize(sur.cut, labels = c("High", "Low"))

cap.cut <- median(merged_data2$Capnocytophaga_percentage, na.rm = TRUE)

merged_data2$Capno_group <- ifelse(
  merged_data2$Capnocytophaga_percentage > cap.cut,
  "High", "Low"
)

table(merged_data2$Capno_group)

merged_data2$Capnocytophaga_cut <-
  ifelse(merged_data2$Capnocytophaga_percentage > 0, "High", "Low")

unique(merged_data2$Capnocytophaga_cut)


threshold <- quantile(merged_data2$Capnocytophaga_percentage, 0.75,na.rm = TRUE)

merged_data2$Capnocytophaga_cut <- ifelse(merged_data2$Capnocytophaga_percentage <=threshold, "Low", "High")



fit <- survfit(Surv(survival_time, event_status) ~Capnocytophaga_cut,  # 创建生存对象 
               data = merged_data2)

pdf("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Pan-cancer atlas of tumor-resident microbiome, immunity and prognosis/预后分析画图/Capnocytophaga_KIRP_curvep_cut.pdf", width = 6, height =7)
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
           title = "Capnocytophaga KIRP" )  # 设置x轴刻度间距
dev.off()









####fusobacterium TCGA-LGG 高丰度高生存率 ----

merged_data2<- merged_data %>% filter(project %in% "TCGA-LGG")
library(survminer)

colnames(merged_data2)

sur.cut <- surv_cutpoint(merged_data2, time = "survival_time",
                         event = "event_status",
                         variables = "Fusobacterium_percentage")
summary(sur.cut)

#使用surv_categorize函数，根据截断值转为二分类变量。
sur.cat <- surv_categorize(sur.cut, labels = c("High", "Low"))


merged_data2$Fusobacterium_cut <-
  ifelse(merged_data2$Fusobacterium_percentage >= 0.006337136, "High", "Low")


threshold <- quantile(merged_data2$Fusobacterium_percentage, 0.75,na.rm = TRUE)

merged_data2$Fusobacterium_cut <- ifelse(merged_data2$Fusobacterium_percentage <= threshold, "Low", "High")


fit <- survfit(Surv(survival_time, event_status) ~Fusobacterium_cut,  # 创建生存对象 
               data = merged_data2)

pdf("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Pan-cancer atlas of tumor-resident microbiome, immunity and prognosis/预后分析画图/Fusobacterium_LGG_curvep_cut.pdf", width = 6, height =7)
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
           break.x.by = 500,
           xlim = c(0, 3000),
           title = "Fusobacterium LGG" )  # 设置x轴刻度间距
dev.off()


sur.cut <- surv_cutpoint(merged_data2, time = "survival_time",
                         event = "event_status",
                         variables = "Streptococcus_percentage")
summary(sur.cut)


#使用surv_categorize函数，根据截断值转为二分类变量。
sur.cat <- surv_categorize(sur.cut, labels = c("High", "Low"))


merged_data2$Streptococcus_cut <-ifelse(merged_data2$Streptococcus_percentage >0, "High", "Low")


threshold <- quantile(merged_data2$Streptococcus_percentage, 0.75,na.rm = TRUE)

merged_data2$Streptococcus_cut <- ifelse(merged_data2$Streptococcus_percentage <=threshold, "Low", "High")


fit <- survfit(Surv(survival_time, event_status) ~Streptococcus_cut,  # 创建生存对象 
               data = merged_data2)

pdf("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Pan-cancer atlas of tumor-resident microbiome, immunity and prognosis/预后分析画图/Streptococcus_LGG_curvep_cut.pdf", width = 6, height =7)
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
           title = "Streptococcus LGG" )  # 设置x轴刻度间距
dev.off()


sur.cut <- surv_cutpoint(merged_data2, time = "survival_time",
                         event = "event_status",
                         variables = "Capnocytophaga_percentage")
summary(sur.cut)


#使用surv_categorize函数，根据截断值转为二分类变量。
sur.cat <- surv_categorize(sur.cut, labels = c("High", "Low"))

cap.cut <- median(merged_data2$Capnocytophaga_percentage, na.rm = TRUE)

merged_data2$Capno_group <- ifelse(
  merged_data2$Capnocytophaga_percentage > cap.cut,
  "High", "Low"
)

table(merged_data2$Capno_group)

merged_data2$Capnocytophaga_cut <-
  ifelse(merged_data2$Capnocytophaga_percentage > 0, "High", "Low")

unique(merged_data2$Capnocytophaga_cut)


threshold <- quantile(merged_data2$Capnocytophaga_percentage, 0.75,na.rm = TRUE)

merged_data2$Capnocytophaga_cut <- ifelse(merged_data2$Capnocytophaga_percentage <=threshold, "Low", "High")



fit <- survfit(Surv(survival_time, event_status) ~Capnocytophaga_cut,  # 创建生存对象 
               data = merged_data2)

pdf("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Pan-cancer atlas of tumor-resident microbiome, immunity and prognosis/预后分析画图/Capnocytophaga_LGG_curvep_cut.pdf", width = 6, height =7)
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
           title = "Capnocytophaga LGG" )  # 设置x轴刻度间距
dev.off()



####fusobacterium TCGA-LIHC ----

merged_data2<- merged_data %>% filter(project %in% "TCGA-LIHC")
library(survminer)

colnames(merged_data2)

sur.cut <- surv_cutpoint(merged_data2, time = "survival_time",
                         event = "event_status",
                         variables = "Fusobacterium_percentage")
summary(sur.cut)

#使用surv_categorize函数，根据截断值转为二分类变量。
sur.cat <- surv_categorize(sur.cut, labels = c("High", "Low"))


merged_data2$Fusobacterium_cut <-
  ifelse(merged_data2$Fusobacterium_percentage > 0, "High", "Low")


threshold <- quantile(merged_data2$Fusobacterium_percentage, 0.75,na.rm = TRUE)

merged_data2$Fusobacterium_cut <- ifelse(merged_data2$Fusobacterium_percentage <=threshold, "Low", "High")


fit <- survfit(Surv(survival_time, event_status) ~Fusobacterium_cut,  # 创建生存对象 
               data = merged_data2)

pdf("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Pan-cancer atlas of tumor-resident microbiome, immunity and prognosis/预后分析画图/Fusobacterium_LIHC_curvep_cut.pdf", width = 6, height =7)
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
           xlim = c(0, 2500),
           break.x.by = 1000,
           title = "Fusobacterium LIHC" )  # 设置x轴刻度间距
dev.off()

sur.cut <- surv_cutpoint(merged_data2, time = "survival_time",
                         event = "event_status",
                         variables = "Streptococcus_percentage")
summary(sur.cut)


#使用surv_categorize函数，根据截断值转为二分类变量。
sur.cat <- surv_categorize(sur.cut, labels = c("High", "Low"))


merged_data2$Streptococcus_cut <-ifelse(merged_data2$Streptococcus_percentage >0.000456361, "High", "Low")


fit <- survfit(Surv(survival_time, event_status) ~Streptococcus_cut,  # 创建生存对象 
               data = merged_data2)

pdf("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Pan-cancer atlas of tumor-resident microbiome, immunity and prognosis/预后分析画图/Streptococcus_LIHC_curvep_cut.pdf", width = 6, height =7)
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
           title = "Streptococcus LIHC" )  # 设置x轴刻度间距
dev.off()


sur.cut <- surv_cutpoint(merged_data2, time = "survival_time",
                         event = "event_status",
                         variables = "Capnocytophaga_percentage")
summary(sur.cut)


#使用surv_categorize函数，根据截断值转为二分类变量。
sur.cat <- surv_categorize(sur.cut, labels = c("High", "Low"))

cap.cut <- median(merged_data2$Capnocytophaga_percentage, na.rm = TRUE)

merged_data2$Capno_group <- ifelse(
  merged_data2$Capnocytophaga_percentage > cap.cut,
  "High", "Low"
)

table(merged_data2$Capno_group)

merged_data2$Capnocytophaga_cut <-
  ifelse(merged_data2$Capnocytophaga_percentage > 0, "High", "Low")

unique(merged_data2$Capnocytophaga_cut)

fit <- survfit(Surv(survival_time, event_status) ~Capnocytophaga_cut,  # 创建生存对象 
               data = merged_data2)

pdf("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Pan-cancer atlas of tumor-resident microbiome, immunity and prognosis/预后分析画图/Capnocytophaga_LIHC_curvep_cut.pdf", width = 6, height =7)
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
           title = "Capnocytophaga LIHC" )  # 设置x轴刻度间距
dev.off()


####fusobacterium TCGA-LUAD ----

merged_data2<- merged_data %>% filter(project %in% "TCGA-LUAD")
library(survminer)

colnames(merged_data2)

sur.cut <- surv_cutpoint(merged_data2, time = "survival_time",
                         event = "event_status",
                         variables = "Fusobacterium_percentage")
summary(sur.cut)

#使用surv_categorize函数，根据截断值转为二分类变量。
sur.cat <- surv_categorize(sur.cut, labels = c("High", "Low"))


merged_data2$Fusobacterium_cut <-
  ifelse(merged_data2$Fusobacterium_percentage > 0, "High", "Low")


threshold <- quantile(merged_data2$Fusobacterium_percentage, 0.75,na.rm = TRUE)

merged_data2$Fusobacterium_cut <- ifelse(merged_data2$Fusobacterium_percentage <= threshold, "Low", "High")


fit <- survfit(Surv(survival_time, event_status) ~Fusobacterium_cut,  # 创建生存对象 
               data = merged_data2)

pdf("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Pan-cancer atlas of tumor-resident microbiome, immunity and prognosis/预后分析画图/Fusobacterium_LUAD_curvep_cut.pdf", width = 6, height =7)
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
           xlim = c(0, 2000),
           break.x.by = 500,
           title = "Fusobacterium LUAD" )  # 设置x轴刻度间距
dev.off()

sur.cut <- surv_cutpoint(merged_data2, time = "survival_time",
                         event = "event_status",
                         variables = "Streptococcus_percentage")
summary(sur.cut)


#使用surv_categorize函数，根据截断值转为二分类变量。
sur.cat <- surv_categorize(sur.cut, labels = c("High", "Low"))
merged_data2$Streptococcus_cut <-ifelse(merged_data2$Streptococcus_percentage >0.000456361, "High", "Low")

fit <- survfit(Surv(survival_time, event_status) ~Streptococcus_cut,  # 创建生存对象 
               data = merged_data2)

pdf("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Pan-cancer atlas of tumor-resident microbiome, immunity and prognosis/预后分析画图/Streptococcus_LUAD_curvep_cut.pdf", width = 6, height =7)
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
           title = "Streptococcus LUAD" )  # 设置x轴刻度间距
dev.off()


sur.cut <- surv_cutpoint(merged_data2, time = "survival_time",
                         event = "event_status",
                         variables = "Capnocytophaga_percentage")
summary(sur.cut)


#使用surv_categorize函数，根据截断值转为二分类变量。
sur.cat <- surv_categorize(sur.cut, labels = c("High", "Low"))

cap.cut <- median(merged_data2$Capnocytophaga_percentage, na.rm = TRUE)

merged_data2$Capno_group <- ifelse(
  merged_data2$Capnocytophaga_percentage > cap.cut,
  "High", "Low"
)

table(merged_data2$Capno_group)

merged_data2$Capnocytophaga_cut <-
  ifelse(merged_data2$Capnocytophaga_percentage > 0.000456361, "High", "Low")

unique(merged_data2$Capnocytophaga_cut)

fit <- survfit(Surv(survival_time, event_status) ~Capnocytophaga_cut,  # 创建生存对象 
               data = merged_data2)

pdf("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Pan-cancer atlas of tumor-resident microbiome, immunity and prognosis/预后分析画图/Capnocytophaga_LUAD_curvep_cut.pdf", width = 6, height =7)
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
           title = "Capnocytophaga LUAD" )  # 设置x轴刻度间距
dev.off()


####fusobacterium TCGA-LUSC ----

merged_data2<- merged_data %>% filter(project %in% "TCGA-LUSC")
library(survminer)

colnames(merged_data2)

sur.cut <- surv_cutpoint(merged_data2, time = "survival_time",
                         event = "event_status",
                         variables = "Fusobacterium_percentage")
summary(sur.cut)

#使用surv_categorize函数，根据截断值转为二分类变量。
sur.cat <- surv_categorize(sur.cut, labels = c("High", "Low"))


merged_data2$Fusobacterium_cut <-
  ifelse(merged_data2$Fusobacterium_percentage > 0, "High", "Low")


fit <- survfit(Surv(survival_time, event_status) ~Fusobacterium_cut,  # 创建生存对象 
               data = merged_data2)

pdf("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Pan-cancer atlas of tumor-resident microbiome, immunity and prognosis/预后分析画图/Fusobacterium_LUSC_curvep_cut.pdf", width = 6, height =7)
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
           xlim = c(0, 4000),
           break.x.by = 1000,
           title = "Fusobacterium LUSC" )  # 设置x轴刻度间距
dev.off()

sur.cut <- surv_cutpoint(merged_data2, time = "survival_time",
                         event = "event_status",
                         variables = "Streptococcus_percentage")
summary(sur.cut)


#使用surv_categorize函数，根据截断值转为二分类变量。
sur.cat <- surv_categorize(sur.cut, labels = c("High", "Low"))


merged_data2$Streptococcus_cut <-ifelse(merged_data2$Streptococcus_percentage >0, "High", "Low")


fit <- survfit(Surv(survival_time, event_status) ~Streptococcus_cut,  # 创建生存对象 
               data = merged_data2)

pdf("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Pan-cancer atlas of tumor-resident microbiome, immunity and prognosis/预后分析画图/Streptococcus_LUSC_curvep_cut.pdf", width = 6, height =7)
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


sur.cut <- surv_cutpoint(merged_data2, time = "survival_time",
                         event = "event_status",
                         variables = "Capnocytophaga_percentage")
summary(sur.cut)


#使用surv_categorize函数，根据截断值转为二分类变量。
sur.cat <- surv_categorize(sur.cut, labels = c("High", "Low"))

cap.cut <- median(merged_data2$Capnocytophaga_percentage, na.rm = TRUE)

merged_data2$Capno_group <- ifelse(
  merged_data2$Capnocytophaga_percentage > cap.cut,
  "High", "Low"
)

table(merged_data2$Capno_group)

merged_data2$Capnocytophaga_cut <-
  ifelse(merged_data2$Capnocytophaga_percentage > 0, "High", "Low")

unique(merged_data2$Capnocytophaga_cut)

fit <- survfit(Surv(survival_time, event_status) ~Capnocytophaga_cut,  # 创建生存对象 
               data = merged_data2)

pdf("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Pan-cancer atlas of tumor-resident microbiome, immunity and prognosis/预后分析画图/Capnocytophaga_LUSC_curvep_cut.pdf", width = 6, height =7)
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
           title = "Capnocytophaga LUSC" )  # 设置x轴刻度间距
dev.off()






####fusobacterium TCGA-OV ----

merged_data2<- merged_data %>% filter(project %in% "TCGA-OV")
library(survminer)

colnames(merged_data2)

sur.cut <- surv_cutpoint(merged_data2, time = "survival_time",
                         event = "event_status",
                         variables = "Fusobacterium_percentage")
summary(sur.cut)

#使用surv_categorize函数，根据截断值转为二分类变量。
sur.cat <- surv_categorize(sur.cut, labels = c("High", "Low"))


merged_data2$Fusobacterium_cut <-
  ifelse(merged_data2$Fusobacterium_percentage > 0, "High", "Low")


threshold <- quantile(merged_data2$Fusobacterium_percentage, 0.75,na.rm = TRUE)

merged_data2$Fusobacterium_cut <- ifelse(merged_data2$Fusobacterium_percentage <= threshold, "Low", "High")


fit <- survfit(Surv(survival_time, event_status) ~Fusobacterium_cut,  # 创建生存对象 
               data = merged_data2)

pdf("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Pan-cancer atlas of tumor-resident microbiome, immunity and prognosis/预后分析画图/Fusobacterium_OV_curvep_cut.pdf", width = 6, height =7)
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
           xlim = c(0, 2000),
           break.x.by = 500,
           title = "Fusobacterium OV" )  # 设置x轴刻度间距
dev.off()

sur.cut <- surv_cutpoint(merged_data2, time = "survival_time",
                         event = "event_status",
                         variables = "Streptococcus_percentage")
summary(sur.cut)


#使用surv_categorize函数，根据截断值转为二分类变量。
sur.cat <- surv_categorize(sur.cut, labels = c("High", "Low"))
merged_data2$Streptococcus_cut <-ifelse(merged_data2$Streptococcus_percentage >0.000456361, "High", "Low")


threshold <- quantile(merged_data2$Streptococcus_percentage, 0.75,na.rm = TRUE)

merged_data2$Streptococcus_cut <- ifelse(merged_data2$Streptococcus_percentage <= threshold, "Low", "High")



fit <- survfit(Surv(survival_time, event_status) ~Streptococcus_cut,  # 创建生存对象 
               data = merged_data2)

pdf("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Pan-cancer atlas of tumor-resident microbiome, immunity and prognosis/预后分析画图/Streptococcus_OV_curvep_cut.pdf", width = 6, height =7)
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
           title = "Streptococcus OV" )  # 设置x轴刻度间距
dev.off()


sur.cut <- surv_cutpoint(merged_data2, time = "survival_time",
                         event = "event_status",
                         variables = "Capnocytophaga_percentage")
summary(sur.cut)


#使用surv_categorize函数，根据截断值转为二分类变量。
sur.cat <- surv_categorize(sur.cut, labels = c("High", "Low"))

cap.cut <- median(merged_data2$Capnocytophaga_percentage, na.rm = TRUE)

merged_data2$Capno_group <- ifelse(
  merged_data2$Capnocytophaga_percentage > cap.cut,
  "High", "Low"
)

table(merged_data2$Capno_group)

merged_data2$Capnocytophaga_cut <-
  ifelse(merged_data2$Capnocytophaga_percentage > 0.000456361, "High", "Low")

unique(merged_data2$Capnocytophaga_cut)


threshold <- quantile(merged_data2$Capnocytophaga_percentage, 0.75,na.rm = TRUE)

merged_data2$Capnocytophaga_cut <- ifelse(merged_data2$Capnocytophaga_percentage <= threshold, "Low", "High")


fit <- survfit(Surv(survival_time, event_status) ~Capnocytophaga_cut,  # 创建生存对象 
               data = merged_data2)

pdf("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Pan-cancer atlas of tumor-resident microbiome, immunity and prognosis/预后分析画图/Capnocytophaga_OV_curvep_cut.pdf", width = 6, height =7)
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
           title = "Capnocytophaga OV" )  # 设置x轴刻度间距
dev.off()

####fusobacterium TCGA-PAAD ----

merged_data2<- merged_data %>% filter(project %in% "TCGA-PAAD")
library(survminer)

colnames(merged_data2)

sur.cut <- surv_cutpoint(merged_data2, time = "survival_time",
                         event = "event_status",
                         variables = "Fusobacterium_percentage")
summary(sur.cut)

#使用surv_categorize函数，根据截断值转为二分类变量。
sur.cat <- surv_categorize(sur.cut, labels = c("High", "Low"))


merged_data2$Fusobacterium_cut <-
  ifelse(merged_data2$Fusobacterium_percentage > 0, "High", "Low")


fit <- survfit(Surv(survival_time, event_status) ~Fusobacterium_cut,  # 创建生存对象 
               data = merged_data2)

pdf("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Pan-cancer atlas of tumor-resident microbiome, immunity and prognosis/预后分析画图/Fusobacterium_PAAD_curvep_cut.pdf", width = 6, height =7)
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
           xlim = c(0, 4000),
           break.x.by = 1000,
           title = "Fusobacterium PAAD" )  # 设置x轴刻度间距
dev.off()

sur.cut <- surv_cutpoint(merged_data2, time = "survival_time",
                         event = "event_status",
                         variables = "Streptococcus_percentage")
summary(sur.cut)


#使用surv_categorize函数，根据截断值转为二分类变量。
sur.cat <- surv_categorize(sur.cut, labels = c("High", "Low"))


merged_data2$Streptococcus_cut <-ifelse(merged_data2$Streptococcus_percentage >0, "High", "Low")


fit <- survfit(Surv(survival_time, event_status) ~Streptococcus_cut,  # 创建生存对象 
               data = merged_data2)

pdf("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Pan-cancer atlas of tumor-resident microbiome, immunity and prognosis/预后分析画图/Streptococcus_PAAD_curvep_cut.pdf", width = 6, height =7)
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
           title = "Streptococcus PAAD" )  # 设置x轴刻度间距
dev.off()


sur.cut <- surv_cutpoint(merged_data2, time = "survival_time",
                         event = "event_status",
                         variables = "Capnocytophaga_percentage")
summary(sur.cut)


#使用surv_categorize函数，根据截断值转为二分类变量。
sur.cat <- surv_categorize(sur.cut, labels = c("High", "Low"))

cap.cut <- median(merged_data2$Capnocytophaga_percentage, na.rm = TRUE)

merged_data2$Capno_group <- ifelse(
  merged_data2$Capnocytophaga_percentage > cap.cut,
  "High", "Low"
)

table(merged_data2$Capno_group)

merged_data2$Capnocytophaga_cut <-
  ifelse(merged_data2$Capnocytophaga_percentage > 0, "High", "Low")

unique(merged_data2$Capnocytophaga_cut)

fit <- survfit(Surv(survival_time, event_status) ~Capnocytophaga_cut,  # 创建生存对象 
               data = merged_data2)

pdf("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Pan-cancer atlas of tumor-resident microbiome, immunity and prognosis/预后分析画图/Capnocytophaga_PAAD_curvep_cut.pdf", width = 6, height =7)
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
           title = "Capnocytophaga PAAD" )  # 设置x轴刻度间距
dev.off()





####fusobacterium TCGA-PCPG ----

merged_data2<- merged_data %>% filter(project %in% "TCGA-PCPG")
library(survminer)

colnames(merged_data2)

sur.cut <- surv_cutpoint(merged_data2, time = "survival_time",
                         event = "event_status",
                         variables = "Fusobacterium_percentage")
summary(sur.cut)

#使用surv_categorize函数，根据截断值转为二分类变量。
sur.cat <- surv_categorize(sur.cut, labels = c("High", "Low"))


merged_data2$Fusobacterium_cut <-
  ifelse(merged_data2$Fusobacterium_percentage > 0, "High", "Low")

threshold <- quantile(merged_data2$Fusobacterium_percentage, 0.75,na.rm = TRUE)

merged_data2$Fusobacterium_cut <- ifelse(merged_data2$Fusobacterium_percentage <= threshold, "Low", "High")


fit <- survfit(Surv(survival_time, event_status) ~Fusobacterium_cut,  # 创建生存对象 
               data = merged_data2)

pdf("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Pan-cancer atlas of tumor-resident microbiome, immunity and prognosis/预后分析画图/Fusobacterium_PCPG_curvep_cut.pdf", width = 6, height =7)
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
           xlim = c(0, 4000),
           break.x.by = 1000,
           title = "Fusobacterium PCPG" )  # 设置x轴刻度间距
dev.off()

sur.cut <- surv_cutpoint(merged_data2, time = "survival_time",
                         event = "event_status",
                         variables = "Streptococcus_percentage")
summary(sur.cut)


#使用surv_categorize函数，根据截断值转为二分类变量。
sur.cat <- surv_categorize(sur.cut, labels = c("High", "Low"))


merged_data2$Streptococcus_cut <-ifelse(merged_data2$Streptococcus_percentage >0, "High", "Low")


threshold <- quantile(merged_data2$Streptococcus_percentage, 0.75,na.rm = TRUE)

merged_data2$Streptococcus_cut <- ifelse(merged_data2$Streptococcus_percentage <= threshold, "Low", "High")

fit <- survfit(Surv(survival_time, event_status) ~Streptococcus_cut,  # 创建生存对象 
               data = merged_data2)

pdf("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Pan-cancer atlas of tumor-resident microbiome, immunity and prognosis/预后分析画图/Streptococcus_PCPG_curvep_cut.pdf", width = 6, height =7)
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
           title = "Streptococcus PCPG" )  # 设置x轴刻度间距
dev.off()


sur.cut <- surv_cutpoint(merged_data2, time = "survival_time",
                         event = "event_status",
                         variables = "Capnocytophaga_percentage")
summary(sur.cut)


#使用surv_categorize函数，根据截断值转为二分类变量。
sur.cat <- surv_categorize(sur.cut, labels = c("High", "Low"))

cap.cut <- median(merged_data2$Capnocytophaga_percentage, na.rm = TRUE)

merged_data2$Capno_group <- ifelse(
  merged_data2$Capnocytophaga_percentage > cap.cut,
  "High", "Low"
)

table(merged_data2$Capno_group)

merged_data2$Capnocytophaga_cut <-
  ifelse(merged_data2$Capnocytophaga_percentage > 0, "High", "Low")

unique(merged_data2$Capnocytophaga_cut)

fit <- survfit(Surv(survival_time, event_status) ~Capnocytophaga_cut,  # 创建生存对象 
               data = merged_data2)

pdf("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Pan-cancer atlas of tumor-resident microbiome, immunity and prognosis/预后分析画图/Capnocytophaga_PCPG_curvep_cut.pdf", width = 6, height =7)
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
           title = "Capnocytophaga PCPG" )  # 设置x轴刻度间距
dev.off()








####fusobacterium TCGA-PRAD ----

merged_data2<- merged_data %>% filter(project %in% "TCGA-PRAD")
library(survminer)

colnames(merged_data2)

sur.cut <- surv_cutpoint(merged_data2, time = "survival_time",
                         event = "event_status",
                         variables = "Fusobacterium_percentage")
summary(sur.cut)

#使用surv_categorize函数，根据截断值转为二分类变量。
sur.cat <- surv_categorize(sur.cut, labels = c("High", "Low"))


merged_data2$Fusobacterium_cut <-
  ifelse(merged_data2$Fusobacterium_percentage > 0, "High", "Low")

threshold <- quantile(merged_data2$Fusobacterium_percentage, 0.75,na.rm = TRUE)

merged_data2$Fusobacterium_cut <- ifelse(merged_data2$Fusobacterium_percentage <= threshold, "Low", "High")


fit <- survfit(Surv(survival_time, event_status) ~Fusobacterium_cut,  # 创建生存对象 
               data = merged_data2)

pdf("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Pan-cancer atlas of tumor-resident microbiome, immunity and prognosis/预后分析画图/Fusobacterium_PRAD_curvep_cut.pdf", width = 6, height =7)
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
           xlim = c(0, 4000),
           break.x.by = 1000,
           title = "Fusobacterium PRAD" )  # 设置x轴刻度间距
dev.off()

sur.cut <- surv_cutpoint(merged_data2, time = "survival_time",
                         event = "event_status",
                         variables = "Streptococcus_percentage")
summary(sur.cut)


#使用surv_categorize函数，根据截断值转为二分类变量。
sur.cat <- surv_categorize(sur.cut, labels = c("High", "Low"))


merged_data2$Streptococcus_cut <-ifelse(merged_data2$Streptococcus_percentage >0, "High", "Low")


threshold <- quantile(merged_data2$Streptococcus_percentage, 0.75,na.rm = TRUE)

merged_data2$Streptococcus_cut <- ifelse(merged_data2$Streptococcus_percentage <= threshold, "Low", "High")

fit <- survfit(Surv(survival_time, event_status) ~Streptococcus_cut,  # 创建生存对象 
               data = merged_data2)

pdf("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Pan-cancer atlas of tumor-resident microbiome, immunity and prognosis/预后分析画图/Streptococcus_PRAD_curvep_cut.pdf", width = 6, height =7)
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
           title = "Streptococcus PRAD" )  # 设置x轴刻度间距
dev.off()


sur.cut <- surv_cutpoint(merged_data2, time = "survival_time",
                         event = "event_status",
                         variables = "Capnocytophaga_percentage")
summary(sur.cut)


#使用surv_categorize函数，根据截断值转为二分类变量。
sur.cat <- surv_categorize(sur.cut, labels = c("High", "Low"))

cap.cut <- median(merged_data2$Capnocytophaga_percentage, na.rm = TRUE)

merged_data2$Capno_group <- ifelse(
  merged_data2$Capnocytophaga_percentage > cap.cut,
  "High", "Low"
)

table(merged_data2$Capno_group)

merged_data2$Capnocytophaga_cut <-
  ifelse(merged_data2$Capnocytophaga_percentage > 0, "High", "Low")

unique(merged_data2$Capnocytophaga_cut)

fit <- survfit(Surv(survival_time, event_status) ~Capnocytophaga_cut,  # 创建生存对象 
               data = merged_data2)

pdf("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Pan-cancer atlas of tumor-resident microbiome, immunity and prognosis/预后分析画图/Capnocytophaga_PRAD_curvep_cut.pdf", width = 6, height =7)
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
           title = "Capnocytophaga PRAD" )  # 设置x轴刻度间距
dev.off()










####fusobacterium TCGA-READ ----

merged_data2<- merged_data %>% filter(project %in% "TCGA-READ")
library(survminer)

colnames(merged_data2)

sur.cut <- surv_cutpoint(merged_data2, time = "survival_time",
                         event = "event_status",
                         variables = "Fusobacterium_percentage")
summary(sur.cut)

#使用surv_categorize函数，根据截断值转为二分类变量。
sur.cat <- surv_categorize(sur.cut, labels = c("High", "Low"))


merged_data2$Fusobacterium_cut <-
  ifelse(merged_data2$Fusobacterium_percentage > 0, "High", "Low")

threshold <- quantile(merged_data2$Fusobacterium_percentage, 0.75,na.rm = TRUE)

merged_data2$Fusobacterium_cut <- ifelse(merged_data2$Fusobacterium_percentage <= threshold, "Low", "High")


fit <- survfit(Surv(survival_time, event_status) ~Fusobacterium_cut,  # 创建生存对象 
               data = merged_data2)

pdf("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Pan-cancer atlas of tumor-resident microbiome, immunity and prognosis/预后分析画图/Fusobacterium_READ_curvep_cut.pdf", width = 6, height =7)
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
           xlim = c(0, 4000),
           break.x.by = 1000,
           title = "Fusobacterium READ" )  # 设置x轴刻度间距
dev.off()

sur.cut <- surv_cutpoint(merged_data2, time = "survival_time",
                         event = "event_status",
                         variables = "Streptococcus_percentage")
summary(sur.cut)


#使用surv_categorize函数，根据截断值转为二分类变量。
sur.cat <- surv_categorize(sur.cut, labels = c("High", "Low"))


merged_data2$Streptococcus_cut <-ifelse(merged_data2$Streptococcus_percentage >0, "High", "Low")


threshold <- quantile(merged_data2$Streptococcus_percentage, 0.75,na.rm = TRUE)

merged_data2$Streptococcus_cut <- ifelse(merged_data2$Streptococcus_percentage <= threshold, "Low", "High")

fit <- survfit(Surv(survival_time, event_status) ~Streptococcus_cut,  # 创建生存对象 
               data = merged_data2)

pdf("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Pan-cancer atlas of tumor-resident microbiome, immunity and prognosis/预后分析画图/Streptococcus_READ_curvep_cut.pdf", width = 6, height =7)
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
           title = "Streptococcus READ" )  # 设置x轴刻度间距
dev.off()


sur.cut <- surv_cutpoint(merged_data2, time = "survival_time",
                         event = "event_status",
                         variables = "Capnocytophaga_percentage")
summary(sur.cut)


#使用surv_categorize函数，根据截断值转为二分类变量。
sur.cat <- surv_categorize(sur.cut, labels = c("High", "Low"))

cap.cut <- median(merged_data2$Capnocytophaga_percentage, na.rm = TRUE)

merged_data2$Capno_group <- ifelse(
  merged_data2$Capnocytophaga_percentage > cap.cut,
  "High", "Low"
)

table(merged_data2$Capno_group)

merged_data2$Capnocytophaga_cut <-
  ifelse(merged_data2$Capnocytophaga_percentage > 0, "High", "Low")

unique(merged_data2$Capnocytophaga_cut)


threshold <- quantile(merged_data2$Capnocytophaga_percentage, 0.75,na.rm = TRUE)

merged_data2$Capnocytophaga_cut <- ifelse(merged_data2$Capnocytophaga_percentage <= threshold, "Low", "High")


fit <- survfit(Surv(survival_time, event_status) ~Capnocytophaga_cut,  # 创建生存对象 
               data = merged_data2)

pdf("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Pan-cancer atlas of tumor-resident microbiome, immunity and prognosis/预后分析画图/Capnocytophaga_READ_curvep_cut.pdf", width = 6, height =7)
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
           title = "Capnocytophaga READ" )  # 设置x轴刻度间距
dev.off()












####fusobacterium TCGA-SARC ----

merged_data2<- merged_data %>% filter(project %in% "TCGA-SARC")
library(survminer)

colnames(merged_data2)

sur.cut <- surv_cutpoint(merged_data2, time = "survival_time",
                         event = "event_status",
                         variables = "Fusobacterium_percentage")
summary(sur.cut)

#使用surv_categorize函数，根据截断值转为二分类变量。
sur.cat <- surv_categorize(sur.cut, labels = c("High", "Low"))


merged_data2$Fusobacterium_cut <-
  ifelse(merged_data2$Fusobacterium_percentage > 0, "High", "Low")

threshold <- quantile(merged_data2$Fusobacterium_percentage, 0.75,na.rm = TRUE)

merged_data2$Fusobacterium_cut <- ifelse(merged_data2$Fusobacterium_percentage <= threshold, "Low", "High")


fit <- survfit(Surv(survival_time, event_status) ~Fusobacterium_cut,  # 创建生存对象 
               data = merged_data2)

pdf("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Pan-cancer atlas of tumor-resident microbiome, immunity and prognosis/预后分析画图/Fusobacterium_SARC_curvep_cut.pdf", width = 6, height =7)
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
           xlim = c(0, 4000),
           break.x.by = 1000,
           title = "Fusobacterium SARC" )  # 设置x轴刻度间距
dev.off()

sur.cut <- surv_cutpoint(merged_data2, time = "survival_time",
                         event = "event_status",
                         variables = "Streptococcus_percentage")
summary(sur.cut)


#使用surv_categorize函数，根据截断值转为二分类变量。
sur.cat <- surv_categorize(sur.cut, labels = c("High", "Low"))


merged_data2$Streptococcus_cut <-ifelse(merged_data2$Streptococcus_percentage >0, "High", "Low")


threshold <- quantile(merged_data2$Streptococcus_percentage, 0.75,na.rm = TRUE)

merged_data2$Streptococcus_cut <- ifelse(merged_data2$Streptococcus_percentage <= threshold, "Low", "High")

fit <- survfit(Surv(survival_time, event_status) ~Streptococcus_cut,  # 创建生存对象 
               data = merged_data2)

pdf("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Pan-cancer atlas of tumor-resident microbiome, immunity and prognosis/预后分析画图/Streptococcus_SARC_curvep_cut.pdf", width = 6, height =7)
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
           title = "Streptococcus SARC" )  # 设置x轴刻度间距
dev.off()


sur.cut <- surv_cutpoint(merged_data2, time = "survival_time",
                         event = "event_status",
                         variables = "Capnocytophaga_percentage")
summary(sur.cut)


#使用surv_categorize函数，根据截断值转为二分类变量。
sur.cat <- surv_categorize(sur.cut, labels = c("High", "Low"))

cap.cut <- median(merged_data2$Capnocytophaga_percentage, na.rm = TRUE)

merged_data2$Capno_group <- ifelse(
  merged_data2$Capnocytophaga_percentage > cap.cut,
  "High", "Low"
)

table(merged_data2$Capno_group)

merged_data2$Capnocytophaga_cut <-
  ifelse(merged_data2$Capnocytophaga_percentage > 0, "High", "Low")

unique(merged_data2$Capnocytophaga_cut)


threshold <- quantile(merged_data2$Capnocytophaga_percentage, 0.75,na.rm = TRUE)

merged_data2$Capnocytophaga_cut <- ifelse(merged_data2$Capnocytophaga_percentage <= threshold, "Low", "High")


fit <- survfit(Surv(survival_time, event_status) ~Capnocytophaga_cut,  # 创建生存对象 
               data = merged_data2)

pdf("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Pan-cancer atlas of tumor-resident microbiome, immunity and prognosis/预后分析画图/Capnocytophaga_SARC_curvep_cut.pdf", width = 6, height =7)
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
           title = "Capnocytophaga SARC" )  # 设置x轴刻度间距
dev.off()











####fusobacterium TCGA-SKCM ----

merged_data2<- merged_data %>% filter(project %in% "TCGA-SKCM")
library(survminer)

colnames(merged_data2)

sur.cut <- surv_cutpoint(merged_data2, time = "survival_time",
                         event = "event_status",
                         variables = "Fusobacterium_percentage")
summary(sur.cut)

#使用surv_categorize函数，根据截断值转为二分类变量。
sur.cat <- surv_categorize(sur.cut, labels = c("High", "Low"))


merged_data2$Fusobacterium_cut <-
  ifelse(merged_data2$Fusobacterium_percentage > 0, "High", "Low")

threshold <- quantile(merged_data2$Fusobacterium_percentage, 0.75,na.rm = TRUE)

merged_data2$Fusobacterium_cut <- ifelse(merged_data2$Fusobacterium_percentage <= threshold, "Low", "High")


fit <- survfit(Surv(survival_time, event_status) ~Fusobacterium_cut,  # 创建生存对象 
               data = merged_data2)

pdf("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Pan-cancer atlas of tumor-resident microbiome, immunity and prognosis/预后分析画图/Fusobacterium_SKCM_curvep_cut.pdf", width = 6, height =7)
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
           xlim = c(0, 4000),
           break.x.by = 1000,
           title = "Fusobacterium SKCM" )  # 设置x轴刻度间距
dev.off()

sur.cut <- surv_cutpoint(merged_data2, time = "survival_time",
                         event = "event_status",
                         variables = "Streptococcus_percentage")
summary(sur.cut)


#使用surv_categorize函数，根据截断值转为二分类变量。
sur.cat <- surv_categorize(sur.cut, labels = c("High", "Low"))


merged_data2$Streptococcus_cut <-ifelse(merged_data2$Streptococcus_percentage >0, "High", "Low")


threshold <- quantile(merged_data2$Streptococcus_percentage, 0.75,na.rm = TRUE)

merged_data2$Streptococcus_cut <- ifelse(merged_data2$Streptococcus_percentage <= threshold, "Low", "High")

fit <- survfit(Surv(survival_time, event_status) ~Streptococcus_cut,  # 创建生存对象 
               data = merged_data2)

pdf("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Pan-cancer atlas of tumor-resident microbiome, immunity and prognosis/预后分析画图/Streptococcus_SKCM_curvep_cut.pdf", width = 6, height =7)
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
           title = "Streptococcus SKCM" )  # 设置x轴刻度间距
dev.off()


sur.cut <- surv_cutpoint(merged_data2, time = "survival_time",
                         event = "event_status",
                         variables = "Capnocytophaga_percentage")
summary(sur.cut)


#使用surv_categorize函数，根据截断值转为二分类变量。
sur.cat <- surv_categorize(sur.cut, labels = c("High", "Low"))

cap.cut <- median(merged_data2$Capnocytophaga_percentage, na.rm = TRUE)

merged_data2$Capno_group <- ifelse(
  merged_data2$Capnocytophaga_percentage > cap.cut,
  "High", "Low"
)

table(merged_data2$Capno_group)

merged_data2$Capnocytophaga_cut <-
  ifelse(merged_data2$Capnocytophaga_percentage > 0, "High", "Low")

unique(merged_data2$Capnocytophaga_cut)


threshold <- quantile(merged_data2$Capnocytophaga_percentage, 0.75,na.rm = TRUE)

merged_data2$Capnocytophaga_cut <- ifelse(merged_data2$Capnocytophaga_percentage <= threshold, "Low", "High")


fit <- survfit(Surv(survival_time, event_status) ~Capnocytophaga_cut,  # 创建生存对象 
               data = merged_data2)

pdf("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Pan-cancer atlas of tumor-resident microbiome, immunity and prognosis/预后分析画图/Capnocytophaga_SKCM_curvep_cut.pdf", width = 6, height =7)
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
           title = "Capnocytophaga SKCM" )  # 设置x轴刻度间距
dev.off()













####fusobacterium TCGA-STAD ----

merged_data2<- merged_data %>% filter(project %in% "TCGA-STAD")
library(survminer)

colnames(merged_data2)

sur.cut <- surv_cutpoint(merged_data2, time = "survival_time",
                         event = "event_status",
                         variables = "Fusobacterium_percentage")
summary(sur.cut)

#使用surv_categorize函数，根据截断值转为二分类变量。
sur.cat <- surv_categorize(sur.cut, labels = c("High", "Low"))


merged_data2$Fusobacterium_cut <-
  ifelse(merged_data2$Fusobacterium_percentage > 0, "High", "Low")


fit <- survfit(Surv(survival_time, event_status) ~Fusobacterium_cut,  # 创建生存对象 
               data = merged_data2)

pdf("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Pan-cancer atlas of tumor-resident microbiome, immunity and prognosis/预后分析画图/Fusobacterium_STAD_curvep_cut.pdf", width = 6, height =7)
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
           xlim = c(0, 4000),
           break.x.by = 1000,
           title = "Fusobacterium STAD" )  # 设置x轴刻度间距
dev.off()

sur.cut <- surv_cutpoint(merged_data2, time = "survival_time",
                         event = "event_status",
                         variables = "Streptococcus_percentage")
summary(sur.cut)


#使用surv_categorize函数，根据截断值转为二分类变量。
sur.cat <- surv_categorize(sur.cut, labels = c("High", "Low"))


merged_data2$Streptococcus_cut <-ifelse(merged_data2$Streptococcus_percentage >0, "High", "Low")


fit <- survfit(Surv(survival_time, event_status) ~Streptococcus_cut,  # 创建生存对象 
               data = merged_data2)

pdf("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Pan-cancer atlas of tumor-resident microbiome, immunity and prognosis/预后分析画图/Streptococcus_STAD_curvep_cut.pdf", width = 6, height =7)
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
           title = "Streptococcus STAD" )  # 设置x轴刻度间距
dev.off()


sur.cut <- surv_cutpoint(merged_data2, time = "survival_time",
                         event = "event_status",
                         variables = "Capnocytophaga_percentage")
summary(sur.cut)


#使用surv_categorize函数，根据截断值转为二分类变量。
sur.cat <- surv_categorize(sur.cut, labels = c("High", "Low"))

cap.cut <- median(merged_data2$Capnocytophaga_percentage, na.rm = TRUE)

merged_data2$Capno_group <- ifelse(
  merged_data2$Capnocytophaga_percentage > cap.cut,
  "High", "Low"
)

table(merged_data2$Capno_group)

merged_data2$Capnocytophaga_cut <-
  ifelse(merged_data2$Capnocytophaga_percentage > 0, "High", "Low")

unique(merged_data2$Capnocytophaga_cut)

fit <- survfit(Surv(survival_time, event_status) ~Capnocytophaga_cut,  # 创建生存对象 
               data = merged_data2)

pdf("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Pan-cancer atlas of tumor-resident microbiome, immunity and prognosis/预后分析画图/Capnocytophaga_STAD_curvep_cut.pdf", width = 6, height =7)
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
           title = "Capnocytophaga STAD" )  # 设置x轴刻度间距
dev.off()






####fusobacterium TCGA-THCA ----

merged_data2<- merged_data %>% filter(project %in% "TCGA-THCA")
library(survminer)

colnames(merged_data2)

sur.cut <- surv_cutpoint(merged_data2, time = "survival_time",
                         event = "event_status",
                         variables = "Fusobacterium_percentage")
summary(sur.cut)

#使用surv_categorize函数，根据截断值转为二分类变量。
sur.cat <- surv_categorize(sur.cut, labels = c("High", "Low"))


merged_data2$Fusobacterium_cut <-
  ifelse(merged_data2$Fusobacterium_percentage > 0, "High", "Low")

threshold <- quantile(merged_data2$Fusobacterium_percentage, 0.75,na.rm = TRUE)

merged_data2$Fusobacterium_cut <- ifelse(merged_data2$Fusobacterium_percentage <= threshold, "Low", "High")


fit <- survfit(Surv(survival_time, event_status) ~Fusobacterium_cut,  # 创建生存对象 
               data = merged_data2)

pdf("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Pan-cancer atlas of tumor-resident microbiome, immunity and prognosis/预后分析画图/Fusobacterium_THCA_curvep_cut.pdf", width = 6, height =7)
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
           xlim = c(0, 4000),
           break.x.by = 1000,
           title = "Fusobacterium THCA" )  # 设置x轴刻度间距
dev.off()

sur.cut <- surv_cutpoint(merged_data2, time = "survival_time",
                         event = "event_status",
                         variables = "Streptococcus_percentage")
summary(sur.cut)


#使用surv_categorize函数，根据截断值转为二分类变量。
sur.cat <- surv_categorize(sur.cut, labels = c("High", "Low"))


merged_data2$Streptococcus_cut <-ifelse(merged_data2$Streptococcus_percentage >0, "High", "Low")


threshold <- quantile(merged_data2$Streptococcus_percentage, 0.75,na.rm = TRUE)

merged_data2$Streptococcus_cut <- ifelse(merged_data2$Streptococcus_percentage <= threshold, "Low", "High")

fit <- survfit(Surv(survival_time, event_status) ~Streptococcus_cut,  # 创建生存对象 
               data = merged_data2)

pdf("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Pan-cancer atlas of tumor-resident microbiome, immunity and prognosis/预后分析画图/Streptococcus_THCA_curvep_cut.pdf", width = 6, height =7)
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
           title = "Streptococcus THCA" )  # 设置x轴刻度间距
dev.off()


sur.cut <- surv_cutpoint(merged_data2, time = "survival_time",
                         event = "event_status",
                         variables = "Capnocytophaga_percentage")
summary(sur.cut)


#使用surv_categorize函数，根据截断值转为二分类变量。
sur.cat <- surv_categorize(sur.cut, labels = c("High", "Low"))

cap.cut <- median(merged_data2$Capnocytophaga_percentage, na.rm = TRUE)

merged_data2$Capno_group <- ifelse(
  merged_data2$Capnocytophaga_percentage > cap.cut,
  "High", "Low"
)

table(merged_data2$Capno_group)

merged_data2$Capnocytophaga_cut <-
  ifelse(merged_data2$Capnocytophaga_percentage > 0, "High", "Low")

unique(merged_data2$Capnocytophaga_cut)


threshold <- quantile(merged_data2$Capnocytophaga_percentage, 0.75,na.rm = TRUE)

merged_data2$Capnocytophaga_cut <- ifelse(merged_data2$Capnocytophaga_percentage <= threshold, "Low", "High")


fit <- survfit(Surv(survival_time, event_status) ~Capnocytophaga_cut,  # 创建生存对象 
               data = merged_data2)

pdf("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Pan-cancer atlas of tumor-resident microbiome, immunity and prognosis/预后分析画图/Capnocytophaga_THCA_curvep_cut.pdf", width = 6, height =7)
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
           title = "Capnocytophaga THCA" )  # 设置x轴刻度间距
dev.off()







####fusobacterium TCGA-THYM ----

merged_data2<- merged_data %>% filter(project %in% "TCGA-THYM")
library(survminer)

colnames(merged_data2)

sur.cut <- surv_cutpoint(merged_data2, time = "survival_time",
                         event = "event_status",
                         variables = "Fusobacterium_percentage")
summary(sur.cut)

#使用surv_categorize函数，根据截断值转为二分类变量。
sur.cat <- surv_categorize(sur.cut, labels = c("High", "Low"))


merged_data2$Fusobacterium_cut <-
  ifelse(merged_data2$Fusobacterium_percentage > 0, "High", "Low")

threshold <- quantile(merged_data2$Fusobacterium_percentage, 0.75,na.rm = TRUE)

merged_data2$Fusobacterium_cut <- ifelse(merged_data2$Fusobacterium_percentage <= threshold, "Low", "High")


fit <- survfit(Surv(survival_time, event_status) ~Fusobacterium_cut,  # 创建生存对象 
               data = merged_data2)

pdf("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Pan-cancer atlas of tumor-resident microbiome, immunity and prognosis/预后分析画图/Fusobacterium_THYM_curvep_cut.pdf", width = 6, height =7)
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
           xlim = c(0, 4000),
           break.x.by = 1000,
           title = "Fusobacterium THYM" )  # 设置x轴刻度间距
dev.off()

sur.cut <- surv_cutpoint(merged_data2, time = "survival_time",
                         event = "event_status",
                         variables = "Streptococcus_percentage")
summary(sur.cut)


#使用surv_categorize函数，根据截断值转为二分类变量。
sur.cat <- surv_categorize(sur.cut, labels = c("High", "Low"))


merged_data2$Streptococcus_cut <-ifelse(merged_data2$Streptococcus_percentage >0, "High", "Low")


threshold <- quantile(merged_data2$Streptococcus_percentage, 0.75,na.rm = TRUE)

merged_data2$Streptococcus_cut <- ifelse(merged_data2$Streptococcus_percentage <= threshold, "Low", "High")

fit <- survfit(Surv(survival_time, event_status) ~Streptococcus_cut,  # 创建生存对象 
               data = merged_data2)

pdf("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Pan-cancer atlas of tumor-resident microbiome, immunity and prognosis/预后分析画图/Streptococcus_THYM_curvep_cut.pdf", width = 6, height =7)
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
           title = "Streptococcus THYM" )  # 设置x轴刻度间距
dev.off()


sur.cut <- surv_cutpoint(merged_data2, time = "survival_time",
                         event = "event_status",
                         variables = "Capnocytophaga_percentage")
summary(sur.cut)


#使用surv_categorize函数，根据截断值转为二分类变量。
sur.cat <- surv_categorize(sur.cut, labels = c("High", "Low"))

cap.cut <- median(merged_data2$Capnocytophaga_percentage, na.rm = TRUE)

merged_data2$Capno_group <- ifelse(
  merged_data2$Capnocytophaga_percentage > cap.cut,
  "High", "Low"
)

table(merged_data2$Capno_group)

merged_data2$Capnocytophaga_cut <-
  ifelse(merged_data2$Capnocytophaga_percentage > 0, "High", "Low")

unique(merged_data2$Capnocytophaga_cut)


threshold <- quantile(merged_data2$Capnocytophaga_percentage, 0.75,na.rm = TRUE)

merged_data2$Capnocytophaga_cut <- ifelse(merged_data2$Capnocytophaga_percentage <= threshold, "Low", "High")


fit <- survfit(Surv(survival_time, event_status) ~Capnocytophaga_cut,  # 创建生存对象 
               data = merged_data2)

pdf("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Pan-cancer atlas of tumor-resident microbiome, immunity and prognosis/预后分析画图/Capnocytophaga_THYM_curvep_cut.pdf", width = 6, height =7)
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
           title = "Capnocytophaga THYM" )  # 设置x轴刻度间距
dev.off()


####fusobacterium TCGA-UCEC ----

merged_data2<- merged_data %>% filter(project %in% "TCGA-UCEC")
library(survminer)

colnames(merged_data2)

sur.cut <- surv_cutpoint(merged_data2, time = "survival_time",
                         event = "event_status",
                         variables = "Fusobacterium_percentage")
summary(sur.cut)

#使用surv_categorize函数，根据截断值转为二分类变量。
sur.cat <- surv_categorize(sur.cut, labels = c("High", "Low"))


merged_data2$Fusobacterium_cut <-
  ifelse(merged_data2$Fusobacterium_percentage > 8.51e-06, "High", "Low")


threshold <- quantile(merged_data2$Fusobacterium_percentage, 0.75,na.rm = TRUE)

merged_data2$Fusobacterium_cut <- ifelse(merged_data2$Fusobacterium_percentage <= threshold, "Low", "High")


fit <- survfit(Surv(survival_time, event_status) ~Fusobacterium_cut,  # 创建生存对象 
               data = merged_data2)

pdf("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Pan-cancer atlas of tumor-resident microbiome, immunity and prognosis/预后分析画图/Fusobacterium_UCEC_curvep_cut.pdf", width = 6, height =7)
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
           xlim = c(0, 4000),
           break.x.by = 1000,
           title = "Fusobacterium UCEC" )  # 设置x轴刻度间距
dev.off()

sur.cut <- surv_cutpoint(merged_data2, time = "survival_time",
                         event = "event_status",
                         variables = "Streptococcus_percentage")
summary(sur.cut)


#使用surv_categorize函数，根据截断值转为二分类变量。
sur.cat <- surv_categorize(sur.cut, labels = c("High", "Low"))


merged_data2$Streptococcus_cut <-ifelse(merged_data2$Streptococcus_percentage >0.000295377, "High", "Low")



threshold <- quantile(merged_data2$Streptococcus_percentage, 0.75,na.rm = TRUE)

merged_data2$Streptococcus_cut <- ifelse(merged_data2$Streptococcus_percentage <= threshold, "Low", "High")



fit <- survfit(Surv(survival_time, event_status) ~Streptococcus_cut,  # 创建生存对象 
               data = merged_data2)

pdf("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Pan-cancer atlas of tumor-resident microbiome, immunity and prognosis/预后分析画图/Streptococcus_UCEC_curvep_cut.pdf", width = 6, height =7)
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
           title = "Streptococcus UCEC" )  # 设置x轴刻度间距
dev.off()


sur.cut <- surv_cutpoint(merged_data2, time = "survival_time",
                         event = "event_status",
                         variables = "Capnocytophaga_percentage")
summary(sur.cut)


#使用surv_categorize函数，根据截断值转为二分类变量。
sur.cat <- surv_categorize(sur.cut, labels = c("High", "Low"))

cap.cut <- median(merged_data2$Capnocytophaga_percentage, na.rm = TRUE)

merged_data2$Capno_group <- ifelse(
  merged_data2$Capnocytophaga_percentage > cap.cut,
  "High", "Low"
)

table(merged_data2$Capno_group)

merged_data2$Capnocytophaga_cut <-
  ifelse(merged_data2$Capnocytophaga_percentage >2.06e-05, "High", "Low")

unique(merged_data2$Capnocytophaga_cut)


threshold <- quantile(merged_data2$Capnocytophaga_percentage, 0.75,na.rm = TRUE)

merged_data2$Capnocytophaga_cut <- ifelse(merged_data2$Capnocytophaga_percentage <= threshold, "Low", "High")


fit <- survfit(Surv(survival_time, event_status) ~Capnocytophaga_cut,  # 创建生存对象 
               data = merged_data2)

pdf("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Pan-cancer atlas of tumor-resident microbiome, immunity and prognosis/预后分析画图/Capnocytophaga_UCEC_curvep_cut.pdf", width = 6, height =7)
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
           title = "Capnocytophaga UCEC" )  # 设置x轴刻度间距
dev.off()




####fusobacterium TCGA-UCS ----

merged_data2<- merged_data %>% filter(project %in% "TCGA-UCS")
library(survminer)

colnames(merged_data2)

sur.cut <- surv_cutpoint(merged_data2, time = "survival_time",
                         event = "event_status",
                         variables = "Fusobacterium_percentage")
summary(sur.cut)

#使用surv_categorize函数，根据截断值转为二分类变量。
sur.cat <- surv_categorize(sur.cut, labels = c("High", "Low"))


merged_data2$Fusobacterium_cut <-
  ifelse(merged_data2$Fusobacterium_percentage > 0.001207527, "High", "Low")


threshold <- quantile(merged_data2$Fusobacterium_percentage, 0.75,na.rm = TRUE)

merged_data2$Fusobacterium_cut <- ifelse(merged_data2$Fusobacterium_percentage <= threshold, "Low", "High")



fit <- survfit(Surv(survival_time, event_status) ~Fusobacterium_cut,  # 创建生存对象 
               data = merged_data2)

pdf("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Pan-cancer atlas of tumor-resident microbiome, immunity and prognosis/预后分析画图/Fusobacterium_UCS_curvep_cut.pdf", width = 6, height =7)
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
           xlim = c(0, 4000),
           break.x.by = 1000,
           title = "Fusobacterium UCS" )  # 设置x轴刻度间距
dev.off()

sur.cut <- surv_cutpoint(merged_data2, time = "survival_time",
                         event = "event_status",
                         variables = "Streptococcus_percentage")
summary(sur.cut)


#使用surv_categorize函数，根据截断值转为二分类变量。
sur.cat <- surv_categorize(sur.cut, labels = c("High", "Low"))


merged_data2$Streptococcus_cut <-ifelse(merged_data2$Streptococcus_percentage >0.000893608, "High", "Low")

threshold <- quantile(merged_data2$Streptococcus_percentage, 0.75,na.rm = TRUE)

merged_data2$Streptococcus_cut <- ifelse(merged_data2$Streptococcus_percentage <= threshold, "Low", "High")


fit <- survfit(Surv(survival_time, event_status) ~Streptococcus_cut,  # 创建生存对象 
               data = merged_data2)

pdf("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Pan-cancer atlas of tumor-resident microbiome, immunity and prognosis/预后分析画图/Streptococcus_UCS_curvep_cut.pdf", width = 6, height =7)
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
           title = "Streptococcus UCS" )  # 设置x轴刻度间距
dev.off()


sur.cut <- surv_cutpoint(merged_data2, time = "survival_time",
                         event = "event_status",
                         variables = "Capnocytophaga_percentage")
summary(sur.cut)


#使用surv_categorize函数，根据截断值转为二分类变量。
sur.cat <- surv_categorize(sur.cut, labels = c("High", "Low"))

cap.cut <- median(merged_data2$Capnocytophaga_percentage, na.rm = TRUE)

merged_data2$Capno_group <- ifelse(
  merged_data2$Capnocytophaga_percentage > cap.cut,
  "High", "Low"
)

table(merged_data2$Capno_group)

merged_data2$Capnocytophaga_cut <-
  ifelse(merged_data2$Capnocytophaga_percentage >0, "High", "Low")

unique(merged_data2$Capnocytophaga_cut)


threshold <- quantile(merged_data2$Capnocytophaga_percentage, 0.75,na.rm = TRUE)

merged_data2$Capnocytophaga_cut <- ifelse(merged_data2$Capnocytophaga_percentage <= threshold, "Low", "High")


fit <- survfit(Surv(survival_time, event_status) ~Capnocytophaga_cut,  # 创建生存对象 
               data = merged_data2)

pdf("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Pan-cancer atlas of tumor-resident microbiome, immunity and prognosis/预后分析画图/Capnocytophaga_UCS_curvep_cut.pdf", width = 6, height =7)
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
           title = "Capnocytophaga UCS" )  # 设置x轴刻度间距
dev.off()



####fusobacterium TCGA-UVM ----

merged_data2<- merged_data %>% filter(project %in% "TCGA-UVM")
library(survminer)

colnames(merged_data2)

sur.cut <- surv_cutpoint(merged_data2, time = "survival_time",
                         event = "event_status",
                         variables = "Fusobacterium_percentage")
summary(sur.cut)

#使用surv_categorize函数，根据截断值转为二分类变量。
sur.cat <- surv_categorize(sur.cut, labels = c("High", "Low"))


merged_data2$Fusobacterium_cut <-
  ifelse(merged_data2$Fusobacterium_percentage > 0, "High", "Low")


threshold <- quantile(merged_data2$Fusobacterium_percentage, 0.75,na.rm = TRUE)

merged_data2$Fusobacterium_cut <- ifelse(merged_data2$Fusobacterium_percentage <= threshold, "Low", "High")


fit <- survfit(Surv(survival_time, event_status) ~Fusobacterium_cut,  # 创建生存对象 
               data = merged_data2)

pdf("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Pan-cancer atlas of tumor-resident microbiome, immunity and prognosis/预后分析画图/Fusobacterium_UVM_curvep_cut.pdf", width = 6, height =7)
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
           xlim = c(0, 2000),
           break.x.by = 500,
           title = "Fusobacterium UVM" )  # 设置x轴刻度间距
dev.off()

sur.cut <- surv_cutpoint(merged_data2, time = "survival_time",
                         event = "event_status",
                         variables = "Streptococcus_percentage")
summary(sur.cut)


#使用surv_categorize函数，根据截断值转为二分类变量。
sur.cat <- surv_categorize(sur.cut, labels = c("High", "Low"))


merged_data2$Streptococcus_cut <-ifelse(merged_data2$Streptococcus_percentage >0, "High", "Low")


threshold <- quantile(merged_data2$Streptococcus_percentage, 0.75,na.rm = TRUE)

merged_data2$Streptococcus_cut <- ifelse(merged_data2$Streptococcus_percentage <= threshold, "Low", "High")



fit <- survfit(Surv(survival_time, event_status) ~Streptococcus_cut,  # 创建生存对象 
               data = merged_data2)

pdf("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Pan-cancer atlas of tumor-resident microbiome, immunity and prognosis/预后分析画图/Streptococcus_UVM_curvep_cut.pdf", width = 6, height =7)
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
           title = "Streptococcus UVM" )  # 设置x轴刻度间距
dev.off()


sur.cut <- surv_cutpoint(merged_data2, time = "survival_time",
                         event = "event_status",
                         variables = "Capnocytophaga_percentage")
summary(sur.cut)


#使用surv_categorize函数，根据截断值转为二分类变量。
sur.cat <- surv_categorize(sur.cut, labels = c("High", "Low"))

cap.cut <- median(merged_data2$Capnocytophaga_percentage, na.rm = TRUE)

merged_data2$Capno_group <- ifelse(
  merged_data2$Capnocytophaga_percentage > cap.cut,
  "High", "Low"
)

table(merged_data2$Capno_group)

merged_data2$Capnocytophaga_cut <-
  ifelse(merged_data2$Capnocytophaga_percentage >0, "High", "Low")

unique(merged_data2$Capnocytophaga_cut)


threshold <- quantile(merged_data2$Capnocytophaga_percentage, 0.75,na.rm = TRUE)

merged_data2$Capnocytophaga_cut <- ifelse(merged_data2$Capnocytophaga_percentage <= threshold, "Low", "High")



fit <- survfit(Surv(survival_time, event_status) ~Capnocytophaga_cut,  # 创建生存对象 
               data = merged_data2)

pdf("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Pan-cancer atlas of tumor-resident microbiome, immunity and prognosis/预后分析画图/Capnocytophaga_UVM_curvep_cut.pdf", width = 6, height =7)
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
           title = "Capnocytophaga UVM" )  # 设置x轴刻度间距
dev.off()




####

data<-read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Pan-cancer atlas of tumor-resident microbiome, immunity and prognosis/TCMbio网站画图/三个物种OS.csv")

library(ggplot2)
library(dplyr)

# 创建数据框
data <- data.frame(
  Bacteria = c("g_Capnocytophaga", "g_Streptococcus", "g_Fusobacterium", 
               "g_Fusobacterium", "g_Capnocytophaga", "g_Streptococcus",
               "g_Fusobacterium", "g_Capnocytophaga", "g_Fusobacterium",
               "g_Capnocytophaga", "g_Capnocytophaga"),
  HR = c(1.457927, 1.056719, 1.122429, 1.178411, 1.146303, 
         1.063195, 1.191763, 1.134414, 1.164162, 0.843395, 1.541184),
  HR.up = c(1.964529, 1.113909, 1.190603, 1.355164, 1.28147, 
            1.101317, 1.362314, 1.278672, 1.348957, 0.991055, 2.229315),
  HR.low = c(1.081965, 1.002466, 1.058158, 1.024711, 1.025393, 
             1.026393, 1.042563, 1.006431, 1.004682, 0.717735, 1.065461),
  P.value = c(0.013224, 0.040214, 0.000124, 0.021318, 0.016355, 
              0.000651, 0.010147, 0.038931, 0.043166, 0.038535, 0.02164),
  CancerType = c("ACC", "BLCA", "CESC", "GBM", "KIRC", 
                 "KIRC", "KIRP", "LUSC", "THCA", "UCEC", "USC")
)

# 创建分组标签
data$Group <- paste(data$Bacteria, data$CancerType, sep = " - ")

# 按HR值排序
data <- data %>% arrange(HR)

# 创建森林图
ggplot(data, aes(x = HR, y = reorder(Group, HR), 
                 xmin = HR.low, xmax = HR.up, 
                 color = Bacteria, shape = HR > 1)) +
  geom_point(size = 3) +
  geom_errorbarh(height = 0.2) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
  scale_x_continuous(trans = "log10", 
                     breaks = c(0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2)) +
  scale_shape_manual(values = c(6, 16), guide = "none") +
  labs(x = "Hazard Ratio (HR) with 95% CI (log scale)", 
       y = "Bacteria - Cancer Type", 
       title = "Association between Bacteria and Cancer Risk",
       subtitle = "HR > 1 indicates increased risk, HR < 1 indicates decreased risk") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        panel.grid.major.y = element_blank())



p<-ggplot(data, aes(x = HR, y = reorder(Group, HR), 
                 xmin = HR.low, xmax = HR.up, 
                 color = Bacteria, shape = HR > 1)) +
  geom_point(size = 3) +  # 增大点的大小
  geom_errorbarh(height = 0.2, size = 1) +  # 加粗误差线
  geom_vline(xintercept = 1, linetype = "dashed", color = "red", size = 1.2) +  # 加粗参考线
  scale_x_continuous(trans = "log10", 
                     breaks = c(0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2)) +
  scale_shape_manual(values = c(6, 16), guide = "none") +
  labs(x = "Hazard Ratio (HR) with 95% CI (log scale)", 
       y = "Bacteria - Cancer Type", 
       title = "Association between Bacteria and Cancer Risk",
       subtitle = "HR > 1 indicates increased risk, HR < 1 indicates decreased risk") +
  theme_minimal(base_size = 12) +  # 增大基础字体大小
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),  # 增大标题
        plot.subtitle = element_text(hjust = 0.5, size = 12),  # 增大副标题
        axis.title = element_text(size = 12, face = "bold"),  # 加粗坐标轴标题
        axis.text = element_text(size = 10),  # 增大坐标轴文字
        legend.title = element_text(face = "bold", size = 10),  # 加粗图例标题
        legend.text = element_text(size = 10),
        panel.grid.major.y = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 4)))  # 增大图例中的点大小

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/Pan-cancer atlas of tumor-resident microbiome, immunity and prognosis/TCMbio网站画图/物种森林图.pdf",p, width=9, height=9)


####微生物丰度和治疗效果----
###

install.packages("readxl")
library(readxl)

Metadata <- read_excel("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/关键微生物与免疫治疗效果关系/s41591-022-01698-2.xlsx", 
                   sheet = "Metadata")

View(Metadata)

# colnames(Metadata)
# rownames(Metadata)<- Metadata$Sample
# 
# Metadata <- as.data.frame(Metadata)  # 转换为 data.frame
# rownames(Metadata) <- Metadata[[1]]  # 设置行名
# Metadata <- Metadata[, -1, drop = FALSE]  # 移除原列

# 检查最终数据结构
str(Metadata)
View(Metadata)


featuretable <- read_excel("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/关键微生物与免疫治疗效果关系/s41591-022-01698-2.xlsx", 
                       sheet = "LKT_PPM")

View(featuretable)
# 将第一列设为行名，并删除原第一列
# featuretable <- as.data.frame(featuretable)  # 转换为 data.frame
# rownames(featuretable) <- featuretable[[1]]  # 设置行名
# featuretable <- featuretable[, -1, drop = FALSE]  # 移除原列

# 检查最终数据结构
str(featuretable)
View(featuretable)

featuretable$row.names


genus<-read_excel("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/关键微生物与免疫治疗效果关系/s41591-022-01698-2.xlsx", 
                   sheet = "LKT_featuretable")

genus$row.names
View(genus)


all(featuretable$row.names %in% genus$row.names)

featuretable_df <- as.data.frame(featuretable)
genus_df <- as.data.frame(genus)

featuretable_df$row.names <- as.character(featuretable_df$row.names)
genus_df$row.names <- as.character(genus_df$row.names)
colnames(featuretable_df)[1] <- "Microbe"
colnames(genus_df)[1] <- "Microbe"

library(dplyr)

combined_data <- left_join(featuretable_df, genus_df, by = "Microbe")

View(combined_data)
combined_data$Genus

colnames(combined_data)


library(dplyr)

# 1. 找到样本列的起始和结束索引（这里假设样本列是第2列到第N列，分类信息从后面开始）
# 你的数据看起来样本列是2到188列，分类信息是189列以后
sample_cols <- 2:187

# 2. 按Genus分组，计算每个样本的属丰度总和
genus_abundance <- combined_data %>%
  group_by(Genus) %>%
  summarise(across(all_of(sample_cols), sum, na.rm = TRUE))

# 3. 计算每个样本的总丰度（所有属加起来）
total_abundance <- colSums(genus_abundance[, -1])

# 4. 计算每个属在每个样本中的相对丰度
genus_rel_abundance <- genus_abundance

# 用 sweep 对每列做除法，计算相对丰度
genus_rel_abundance[, -1] <- sweep(genus_abundance[, -1], 2, total_abundance, FUN = "/")

# 查看结果
head(genus_rel_abundance)

View(genus_rel_abundance)

write.csv(genus_rel_abundance,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/关键微生物与免疫治疗效果关系/s41591-022-01698-2_genus_rel_abundance.csv")

write.csv(combined_data,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/关键微生物与免疫治疗效果关系/s41591-022-01698-2_combined_data.csv")


###筛选Fusobacterium、Capnocytophaga、Streptococcus的丰度

genus_rel_abundance$Genus


target_genera <- c("g__Fusobacterium", "g__Capnocytophaga", "g__Streptococcus")

# 从 genus_rel_abundance 中筛选
selected_genera_abundance <- genus_rel_abundance %>%
  filter(Genus %in% target_genera)

# 查看结果
print(selected_genera_abundance)

colnames(selected_genera_abundance)
Metadata$Sample


library(dplyr)
library(tidyr)

# 1. 把 selected_genera_abundance 转成长格式
long_data <- selected_genera_abundance %>%
  pivot_longer(
    cols = -Genus,
    names_to = "Sample",
    values_to = "RelAbundance"
  )

# 2. 左连接元数据，Metadata 有 Sample 列
long_data_with_meta <- long_data %>%
  left_join(Metadata, by = "Sample")

View(long_data_with_meta)


write.csv(long_data_with_meta,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/关键微生物与免疫治疗效果关系/s41591-022-01698-2_long_data_with_meta.csv")


library(ggplot2)
library(ggpubr)

# 提取要绘图的 Genus 列表
genera_to_plot <- unique(long_data_with_meta$Genus)

# 创建图列表
plot_list <- list()

for (genus in genera_to_plot) {
  
  # 筛选该 Genus 的数据
  df_sub <- long_data_with_meta %>% 
    filter(Genus == genus & Treatment_Response_Current != "N_A")
  
  # 创建柱状图/箱线图 + 显著性标注
  p <- ggplot(df_sub, aes(x = Treatment_Response_Current, y = RelAbundance, fill = Treatment_Response_Current)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.2, size = 1.5, alpha = 0.6) +
    stat_compare_means(method = "kruskal.test", label.y = max(df_sub$RelAbundance, na.rm = TRUE) * 1.05) +  # 非参数检验
    theme_minimal() +
    labs(
      title = paste("Genus:", genus),
      x = "Treatment Response",
      y = "Relative Abundance"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5)
    )
  
  plot_list[[genus]] <- p
}
# 显示一个图，比如 Fusobacterium
print(plot_list[["g__Streptococcus"]])

# 保存所有图
for (genus in names(plot_list)) {
  ggsave(
    filename = paste0("Genus_", genus, "_abundance_by_response.png"),
    plot = plot_list[[genus]],
    width = 7, height = 5, dpi = 300
  )
}
# 分组变量

long_data_with_meta$irAE_status <- ifelse(
  is.na(long_data_with_meta$irAE_first) | long_data_with_meta$irAE_first == "N/A",
  "No_irAE",
  "irAE"
)

library(ggpubr)

library(ggpubr)

p1<- ggboxplot(
  long_data_with_meta %>% filter(Genus == "g__Streptococcus"),
  x = "irAE_status",
  y = "RelAbundance",
  color = "irAE_status",
  add = "jitter",
  size = 1.2      # 加粗箱线图边框
) +
  stat_compare_means(method = "wilcox.test", label = "p.signif") +
  labs(
    title = "Streptococcus Abundance vs irAE Status",
    x = "irAE Status",
    y = "Relative Abundance"
  ) +
  theme_bw(base_size = 14) +  # 整体字体基准变大
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 13)
  )

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/关键微生物与免疫治疗效果关系/s41591-022-01698-2_Streptococcus Abundance vs irAE Status.pdf",p1,width=7,height=5)


colnames(long_data_with_meta)


library(survival)
library(survminer)


strepto_data <- long_data_with_meta %>% 
  filter(Genus == "g__Streptococcus")

# 创建高低组（中位数划分）
strepto_data$Abundance_Group <- ifelse(
  strepto_data$RelAbundance >= median(strepto_data$RelAbundance, na.rm = TRUE),
  "High", "Low"
)


# 替换 "N_A" 为真正的 NA
strepto_data$PFS_days <- as.numeric(ifelse(strepto_data$PFS_days == "N_A", NA, strepto_data$PFS_days))
strepto_data$Progression_Event <- as.numeric(ifelse(strepto_data$Progression_Event == "N_A", NA, strepto_data$Progression_Event))

# 去除含NA的行
strepto_clean <- strepto_data %>% 
  filter(!is.na(PFS_days), !is.na(Progression_Event))

# 建立生存对象
surv_obj <- Surv(time = strepto_clean$PFS_days, event = strepto_clean$Progression_Event)

fit <- survfit(surv_obj ~ Abundance_Group, data = strepto_clean)

ggsurvplot(
  fit,
  data = strepto_clean,
  pval = TRUE,
  conf.int = TRUE,
  risk.table = TRUE,
  surv.median.line = "hv",
  palette = c("#E69F00", "#56B4E9"),
  title = "Streptococcus Abundance vs Progression-Free Survival",
  legend.title = "Abundance Group",
  xlab = "Days",
  ylab = "PFS Probability",
  font.main = c(16, "bold"),
  font.x = c(14), font.y = c(14),
  font.tickslab = 12
)



surv_data <- strepto_data %>%
  mutate(
    RelAbundance = as.numeric(RelAbundance),
    PFS_days = as.numeric(ifelse(PFS_days == "N_A", NA, PFS_days)),
    Progression_Event = as.numeric(ifelse(Progression_Event == "N_A", NA, Progression_Event))
  ) %>%
  filter(!is.na(RelAbundance), !is.na(PFS_days), !is.na(Progression_Event))


library(survminer)

# 自动寻找最佳 cutoff 点
sur.cut <- surv_cutpoint(
  data = surv_data,
  time = "PFS_days",
  event = "Progression_Event",
  variables = "RelAbundance"
)

# 查看建议的 cutoff 值
summary(sur.cut)

# 基于 cutpoint 创建分组
surv_grouped <- surv_categorize(sur.cut)
head(surv_grouped)


# 与原数据合并分组标签
surv_data$Abundance_Group <- surv_grouped$RelAbundance

# 构建生存对象
surv_obj <- Surv(surv_data$PFS_days, surv_data$Progression_Event)

# 拟合生存模型
fit <- survfit(surv_obj ~ Abundance_Group, data = surv_data)

# 画 Kaplan–Meier 曲线
ggsurvplot(
  fit,
  data = surv_data,
  pval = TRUE,
  conf.int = TRUE,
  risk.table = TRUE,
  surv.median.line = "hv",
  title = "RelAbundance vs Progression-Free Survival",
  palette = c("#E69F00", "#56B4E9"),
  legend.title = "Streptococcus Abundance",
  xlab = "Days",
  ylab = "PFS Probability",
  font.main = c(16, "bold"),
  font.x = c(14), font.y = c(14),
  font.tickslab = 12
)


surv_data <- strepto_data %>%
  mutate(
    RelAbundance = as.numeric(RelAbundance),
    Overall_Survival_days = as.numeric(ifelse(Overall_Survival_days == "N_A", NA, Overall_Survival_days)),
    Death_Event = as.numeric(ifelse(Death_Event == "N_A", NA, Death_Event))
  ) %>%
  filter(!is.na(RelAbundance), !is.na(Overall_Survival_days), !is.na(Death_Event))

library(survminer)

sur.cut <- surv_cutpoint(
  data = surv_data,
  time = "Overall_Survival_days",
  event = "Death_Event",
  variables = "RelAbundance"
)

summary(sur.cut)


surv_grouped <- surv_categorize(sur.cut)
surv_data$Abundance_Group <- surv_grouped$RelAbundance

library(survival)
library(survminer)

# 创建 Surv 对象
surv_obj <- Surv(surv_data$Overall_Survival_days, surv_data$Death_Event)

# 拟合生存模型
fit <- survfit(surv_obj ~ Abundance_Group, data = surv_data)

# 绘图
ggsurvplot(
  fit,
  data = surv_data,
  pval = TRUE,
  conf.int = TRUE,
  risk.table = TRUE,
  surv.median.line = "hv",
  title = "Streptococcus RelAbundance vs Overall Survival",
  palette = c("#D55E00", "#009E73"),
  legend.title = "Streptococcus Abundance",
  xlab = "Days",
  ylab = "Overall Survival Probability",
  font.main = c(16, "bold"),
  font.x = c(14), font.y = c(14),
  font.tickslab = 12
)


##### s41591-023-02453-x #####----

library(readxl)

metadata_and_clr_abundances <- read_excel("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/关键微生物与免疫治疗效果关系/s41591-023-02453-x.xlsx", sheet = "metadata_and_clr_abundances")

View(metadata_and_clr_abundances)

colnames(metadata_and_clr_abundances)[1:30]


strain_abundances <-read_excel("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/关键微生物与免疫治疗效果关系/s41591-023-02453-x.xlsx", sheet = "strain_abundances") 


strain_abundances$Fusobacterium_sum <- strain_abundances %>%
  dplyr::select(dplyr::contains("Fusobacterium")) %>%
  rowSums(na.rm = TRUE)

total_abundance <- strain_abundances %>%
  dplyr::select(-lane_id) %>%
  rowSums(na.rm = TRUE)

strain_abundances$Fusobacterium_rel_abundance <- 
  strain_abundances$Fusobacterium_sum / total_abundance



metadata <- read_excel("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/关键微生物与免疫治疗效果关系/s41591-023-02453-x.xlsx", sheet = "patient_metadata")


library(dplyr)

# 合并表格，以 metadata 为主表（保留所有 metadata 样本）
merged_data <- left_join(metadata, strain_abundances, by = "lane_id")

write.csv(merged_data,"D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/关键微生物与免疫治疗效果关系/merged_data.csv")


###########################################################################################################

library(ggplot2)
library(ggpubr)
library(dplyr)
library(readxl)
metadata_and_clr_abundances <- read_excel("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/关键微生物与免疫治疗效果关系/s41591-023-02453-x.xlsx", sheet = "metadata_and_clr_abundances")

# 选取微生物列名
# Fusobacterium
colnames(metadata_and_clr_abundances)[1:30]

library(dplyr)
library(ggplot2)
library(ggpubr)

# 提取含 Fusobacterium 的列
fusobacterium_cols <- grep("Fusobacterium", colnames(metadata_and_clr_abundances), value = TRUE)

# 如果有多个相关列，可先合并求和；这里只用第一个演示
metadata_fuso <- metadata_and_clr_abundances %>%
  select(lane_id, R_vs_PD, PFS,OS,not_progression_free,c__NLR,c__PLAT,c__ALB,c__LDH,all_of(fusobacterium_cols)) %>%
  rename(Fusobacterium_CLR = !!fusobacterium_cols[1])  # 若多个列需先相加再重命名

library(dplyr)
library(tidyr)
table(metadata_fuso$R_vs_PD, useNA = "ifany")

metadata_fuso_clean <- metadata_fuso %>%
  mutate(R_vs_PD = na_if(R_vs_PD, "NA")) %>%
  filter(!is.na(R_vs_PD))

library(ggplot2)
library(ggpubr)


metadata_fuso_clean$R_vs_PD <- factor(metadata_fuso_clean$R_vs_PD, levels = c("TRUE", "FALSE"))

p<- ggplot(metadata_fuso_clean, aes(x = R_vs_PD, y = Fusobacterium_CLR, fill = R_vs_PD)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.6) +
  stat_compare_means(method = "wilcox.test", label = "p.signif") +
  labs(
    title = "Difference of Fusobacterium CLR Abundance between R and PD",
    x = "Group",
    y = "Fusobacterium CLR Abundance"
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),       # 去主网格线
    panel.grid.minor = element_blank(),       # 去次网格线
    axis.title = element_text(face = "bold", size = 14),  # 坐标轴标题加粗字体大一点
    axis.text = element_text(face = "bold", size = 12),   # 坐标轴刻度加粗字体稍大
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),  # 标题加粗居中
    #axis.line = element_line(size = 1, colour = "black"),  # 坐标轴线条加粗
    legend.position = "none"
  )+
  theme(
    panel.border = element_rect(color = "black", size = 1.5, fill = NA)
  )

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/关键微生物与免疫治疗效果关系/s41591-023-02453-x_fusobacterium_cols与治疗后效果的关系.pdf",p, height=5, width=5)


p<- ggplot(metadata_fuso_clean, aes(x = R_vs_PD, y = Fusobacterium_CLR, fill = R_vs_PD)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.6) +
  stat_compare_means(method = "wilcox.test", label = "p.signif") +
  labs(
    title = "Difference of Fusobacterium CLR Abundance between R and PD",
    x = "Group",
    y = "Fusobacterium CLR Abundance"
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),       # 去主网格线
    panel.grid.minor = element_blank(),       # 去次网格线
    axis.title = element_text(face = "bold", size = 14),  # 坐标轴标题加粗字体大一点
    axis.text = element_text(face = "bold", size = 12),   # 坐标轴刻度加粗字体稍大
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),  # 标题加粗居中
    #axis.line = element_line(size = 1, colour = "black"),  # 坐标轴线条加粗
    legend.position = "none"
  )+
  theme(
    panel.border = element_rect(color = "black", size = 1.5, fill = NA)
  )

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/关键微生物与免疫治疗效果关系/s41591-023-02453-x_fusobacterium_cols与治疗后效果的关系.pdf",p, height=5, width=5)



library(survival)
library(survminer)

# 构建生存对象
surv_object <- Surv(time = metadata_fuso$PFS, event = metadata_fuso$not_progression_free)


# 拟合模型或绘制KM曲线
fit <- survfit(surv_object ~ Fusobacterium_CLR, data = metadata_fuso)
ggsurvplot(fit, data = metadata_fuso)

head(metadata_fuso$Fusobacterium_CLR)

metadata_fuso$Fuso_group <- ifelse(
  metadata_fuso$Fusobacterium_CLR >= median(metadata_fuso$Fusobacterium_CLR, na.rm = TRUE),
  "High", "Low"
)

quantiles <- quantile(metadata_fuso$Fusobacterium_CLR, probs = c(0.33, 0.67), na.rm = TRUE)

metadata_fuso$Fuso_group <- with(metadata_fuso, ifelse(
  Fusobacterium_CLR <= quantiles[1], "Low",
  ifelse(Fusobacterium_CLR >= quantiles[2], "High", NA)
))


sur.cut <- surv_cutpoint(metadata_fuso, time = "PFS",
                         event = "not_progression_free",
                         variables = "Fusobacterium_CLR")
summary(sur.cut)


metadata_fuso$Fusobacterium_cut <-
  ifelse(metadata_fuso$Fusobacterium_CLR > -0.281867, "High", "Low")



library(survival)
library(survminer)

# 生存对象
surv_obj <- Surv(time = metadata_fuso$PFS, event = metadata_fuso$not_progression_free)

# 生存模型按组拟合
fit <- survfit(surv_obj ~ Fuso_group, data = metadata_fuso)


fit <- survfit(surv_obj ~ Fusobacterium_cut, data = metadata_fuso)

summary(fit)
# 绘图
fit$surv

pdf("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/关键微生物与免疫治疗效果关系/s41591-023-02453-x_fusobacterium_cols与PFS关系.pdf", width = 6, height =7)
ggsurvplot(
  fit,
  data = metadata_fuso,
  pval = TRUE,
  risk.table = TRUE,
  conf.int = TRUE,
  palette = c("#E69F00", "#56B4E9"),  # 颜色可改
  title = "PFS by Fusobacterium Abundance",
  xlab = "Time (months)",
  ylab = "Progression-Free Survival Probability"
)
dev.off()



library(survival)
library(survminer)

# 生存对象
surv_obj <- Surv(time = metadata_fuso$PFS, event = metadata_fuso$not_progression_free)

# 生存模型按组拟合
fit <- survfit(surv_obj ~ Fuso_group, data = metadata_fuso)


fit <- survfit(surv_obj ~ Fusobacterium_cut, data = metadata_fuso)

# 绘图


pdf("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/关键微生物与免疫治疗效果关系/s41591-023-02453-x_fusobacterium_cols与PFS关系.pdf", width = 6, height =7)
ggsurvplot(
  fit,
  data = metadata_fuso,
  pval = TRUE,
  risk.table = TRUE,
  conf.int = TRUE,
  palette = c("#E69F00", "#56B4E9"),  # 颜色可改
  title = "PFS by Fusobacterium Abundance",
  xlab = "Time (months)",
  ylab = "Progression-Free Survival Probability"
)
dev.off()

####c__NLR----


# 绘图
library(ggplot2)
library(ggpubr)
library(dplyr)

# Step 1: 去除异常值（以 IQR 方法为例）
metadata_filtered <- metadata_fuso_clean %>%
  filter(
    !is.na(c__NLR) & !is.na(Fusobacterium_CLR)
  ) %>%
  filter(
    between(c__NLR, quantile(c__NLR, 0.05), quantile(c__NLR, 0.95)) &  # 去除 NLR 异常值
      between(Fusobacterium_CLR, quantile(Fusobacterium_CLR, 0.05), quantile(Fusobacterium_CLR, 0.95))  # 去除 Fusobacterium CLR 异常值
  )

# Step 2: 绘图
p <- ggplot(metadata_filtered, aes(x = c__NLR, y = Fusobacterium_CLR)) +
  # 散点图
  geom_point(color = "steelblue", size = 2.5, alpha = 0.6) +
  # 加粗的回归线
  geom_smooth(method = "lm", color = "red", se = TRUE, size = 1.2) +
  # 添加相关性系数和显著性
  stat_cor(method = "spearman", 
           label.x.npc = "left", label.y.npc = "top", size = 5,
           aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"))) +
  # 去掉网格背景并加粗字体
  theme_classic(base_size = 14) +  # base_size 控制整体字体大小
  theme(
    axis.title = element_text(face = "bold", size = 16),     # 坐标轴标题加粗
    axis.text = element_text(size = 14),                     # 坐标轴刻度字体
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18)  # 标题居中并加粗
  ) +
  # 标题和标签
  labs(
    title = "Correlation between c__NLR and Fusobacterium_CLR",
    x = "c__NLR", 
    y = "Fusobacterium_CLR"
  )

print(p)

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/关键微生物与免疫治疗效果关系/s41591-023-02453-x_fusobacterium_cols与c__NLR关系.pdf", p, width = 6, height =5)


metadata_filtered <- metadata_fuso_clean %>%
  filter(
    !is.na(c__PLAT) & !is.na(Fusobacterium_CLR)
  ) %>%
  filter(
    between(c__PLAT, quantile(c__PLAT, 0.05), quantile(c__PLAT, 0.95)) &  # 去除 NLR 异常值
      between(Fusobacterium_CLR, quantile(Fusobacterium_CLR, 0.05), quantile(Fusobacterium_CLR, 0.95))  # 去除 Fusobacterium CLR 异常值
  )



p <- ggplot(metadata_filtered, aes(x = c__PLAT, y = Fusobacterium_CLR)) +
  # 散点图
  geom_point(color = "steelblue", size = 2.5, alpha = 0.6) +
  # 加粗的回归线
  geom_smooth(method = "lm", color = "red", se = TRUE, size = 1.2) +
  # 添加相关性系数和显著性
  stat_cor(method = "spearman", 
           label.x.npc = "left", label.y.npc = "top", size = 5,
           aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"))) +
  # 去掉网格背景并加粗字体
  theme_classic(base_size = 14) +  # base_size 控制整体字体大小
  theme(
    axis.title = element_text(face = "bold", size = 16),     # 坐标轴标题加粗
    axis.text = element_text(size = 14),                     # 坐标轴刻度字体
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18)  # 标题居中并加粗
  ) +
  # 标题和标签
  labs(
    title = "Correlation between c__PLAT and Fusobacterium_CLR",
    x = "c__PLAT", 
    y = "Fusobacterium_CLR"
  )

print(p)

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/关键微生物与免疫治疗效果关系/s41591-023-02453-x_fusobacterium_cols与c__PLAT关系.pdf", p, width = 6, height =5)



metadata_filtered <- metadata_fuso_clean %>%
  filter(
    !is.na(c__ALB) & !is.na(Fusobacterium_CLR)
  ) %>%
  filter(
    between(c__ALB, quantile(c__ALB, 0.05), quantile(c__ALB, 0.95)) &  # 去除 NLR 异常值
      between(Fusobacterium_CLR, quantile(Fusobacterium_CLR, 0.05), quantile(Fusobacterium_CLR, 0.95))  # 去除 Fusobacterium CLR 异常值
  )


p <- ggplot(metadata_filtered, aes(x = c__ALB, y = Fusobacterium_CLR)) +
  # 散点图
  geom_point(color = "steelblue", size = 2.5, alpha = 0.6) +
  # 加粗的回归线
  geom_smooth(method = "lm", color = "red", se = TRUE, size = 1.2) +
  # 添加相关性系数和显著性
  stat_cor(method = "spearman", 
           label.x.npc = "left", label.y.npc = "top", size = 5,
           aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"))) +
  # 去掉网格背景并加粗字体
  theme_classic(base_size = 14) +  # base_size 控制整体字体大小
  theme(
    axis.title = element_text(face = "bold", size = 16),     # 坐标轴标题加粗
    axis.text = element_text(size = 14),                     # 坐标轴刻度字体
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18)  # 标题居中并加粗
  ) +
  # 标题和标签
  labs(
    title = "Correlation between c__ALB and Fusobacterium_CLR",
    x = "c__ALB", 
    y = "Fusobacterium_CLR"
  )

print(p)

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/关键微生物与免疫治疗效果关系/s41591-023-02453-x_fusobacterium_cols与c__ALB关系.pdf", p, width = 6, height =5)


metadata_filtered <- metadata_fuso_clean %>%
  filter(
    !is.na(c__LDH) & !is.na(Fusobacterium_CLR)
  ) %>%
  filter(
    between(c__LDH, quantile(c__LDH, 0.05), quantile(c__LDH, 0.95)) &  # 去除 NLR 异常值
      between(Fusobacterium_CLR, quantile(Fusobacterium_CLR, 0.05), quantile(Fusobacterium_CLR, 0.95))  # 去除 Fusobacterium CLR 异常值
  )



p <- ggplot(metadata_filtered, aes(x = c__LDH, y = Fusobacterium_CLR)) +
  # 散点图
  geom_point(color = "steelblue", size = 2.5, alpha = 0.6) +
  # 加粗的回归线
  geom_smooth(method = "lm", color = "red", se = TRUE, size = 1.2) +
  # 添加相关性系数和显著性
  stat_cor(method = "spearman", 
           label.x.npc = "left", label.y.npc = "top", size = 5,
           aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"))) +
  # 去掉网格背景并加粗字体
  theme_classic(base_size = 14) +  # base_size 控制整体字体大小
  theme(
    axis.title = element_text(face = "bold", size = 16),     # 坐标轴标题加粗
    axis.text = element_text(size = 14),                     # 坐标轴刻度字体
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18)  # 标题居中并加粗
  ) +
  # 标题和标签
  labs(
    title = "Correlation between c__LDH and Fusobacterium_CLR",
    x = "c__LDH", 
    y = "Fusobacterium_CLR"
  )

print(p)

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/关键微生物与免疫治疗效果关系/s41591-023-02453-x_fusobacterium_cols与c__LDH关系.pdf", p, width = 6, height =5)



quakeslibrary(tidyverse)
library(ggpubr)

# 示例：选择要分析的物种
data_genus <- read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/关键微生物与免疫治疗效果关系/mrr-4-1-4_KELLY_clinical_metadata_data-genus.csv",header = T ,row.names = 1)

df<- data_genus

saliva_normalised_with_metadata <-read.csv("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/关键微生物与免疫治疗效果关系/mrr-4-1-4_KELLY_saliva_normalised_with_metadata.csv",header = T)

df<-saliva_normalised_with_metadata

colnames(df)

# 转换必要的变量类型
df <- df %>%
  mutate(
    PDL1 = as.factor(PDL1),
    CBR = as.factor(CBR),
    tPFS = as.factor(tPFS)
  )

# 绘图列表
my_colors <- c("1" = "#E64B35", "0" = "#4DBBD5")  # 举例

p1 <- ggboxplot(df, x = "PDL1", y = "Capnocytophaga", fill = "PDL1",
                palette = my_colors, add = "jitter", 
                size = 0.8) +  # 加粗箱线图的边框线条
  stat_compare_means(method = "wilcox.test", size = 5) +  # 设置显著性标注的字体大小
  labs(title = paste("PDL1 vs", "Capnocytophaga")) +
  theme_classic(base_size = 14) +  # 设置整体字体基准大小
  theme(
    axis.title = element_text(size = 16, face = "bold"),     # 坐标轴标题加粗
    axis.text = element_text(size = 14),                     # 坐标轴刻度加大
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),  # 标题加粗并居中
    legend.title = element_text(size = 14, face = "bold"),   # 图例标题加粗
    legend.text = element_text(size = 13)                    # 图例文字
  )

print(p1)

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/关键微生物与免疫治疗效果关系/mrr-4-1-4_Capnocytophaga_PDL1.pdf",p1, width = 6, height =5)

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/关键微生物与免疫治疗效果关系/mrr-4-1-4_saliva_normalised_with_metadata_Capnocytophaga_PDL1.pdf",p1, width = 6, height =5)


p2 <- ggboxplot(df, x = "CBR", y = "Capnocytophaga", fill = "CBR",
                palette = my_colors, add = "jitter", 
                size = 0.8) +  # 加粗箱线图的边框线条
  stat_compare_means(method = "wilcox.test", size = 5) +  # 设置显著性标注的字体大小
  labs(title = paste("CBR vs", "Capnocytophaga")) +
  theme_classic(base_size = 14) +  # 设置整体字体基准大小
  theme(
    axis.title = element_text(size = 16, face = "bold"),     # 坐标轴标题加粗
    axis.text = element_text(size = 14),                     # 坐标轴刻度加大
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),  # 标题加粗并居中
    legend.title = element_text(size = 14, face = "bold"),   # 图例标题加粗
    legend.text = element_text(size = 13)                    # 图例文字
  )

print(p2)

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/关键微生物与免疫治疗效果关系/mrr-4-1-4_Capnocytophaga_CBR.pdf",p2, width = 6, height =5)

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/关键微生物与免疫治疗效果关系/mrr-4-1-4_saliva_normalised_with_metadata_Capnocytophaga_CBR.pdf",p2, width = 6, height =5)



p3 <- ggboxplot(df, x = "tPFS", y = "Capnocytophaga", fill = "tPFS",
                palette = my_colors, add = "jitter", 
                size = 0.8) +  # 加粗箱线图的边框线条
  stat_compare_means(method = "wilcox.test", size = 5) +  # 设置显著性标注的字体大小
  labs(title = paste("tPFS vs", "Capnocytophaga")) +
  theme_classic(base_size = 14) +  # 设置整体字体基准大小
  theme(
    axis.title = element_text(size = 16, face = "bold"),     # 坐标轴标题加粗
    axis.text = element_text(size = 14),                     # 坐标轴刻度加大
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),  # 标题加粗并居中
    legend.title = element_text(size = 14, face = "bold"),   # 图例标题加粗
    legend.text = element_text(size = 13)                    # 图例文字
  )

print(p3)


ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/关键微生物与免疫治疗效果关系/mrr-4-1-4_Capnocytophaga_tPFS.pdf",p3, width = 6, height =5)

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/关键微生物与免疫治疗效果关系/mrr-4-1-4_saliva_normalised_with_metadata_Capnocytophaga_tPFS.pdf",p3, width = 6, height =5)




p4 <- ggboxplot(df, x = "NLR", y = "Capnocytophaga", fill = "NLR",
                palette = my_colors, add = "jitter", 
                size = 0.8) +  # 加粗箱线图的边框线条
  stat_compare_means(method = "wilcox.test", size = 5) +  # 设置显著性标注的字体大小
  labs(title = paste("NLR vs", "Capnocytophaga")) +
  theme_classic(base_size = 14) +  # 设置整体字体基准大小
  theme(
    axis.title = element_text(size = 16, face = "bold"),     # 坐标轴标题加粗
    axis.text = element_text(size = 14),                     # 坐标轴刻度加大
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),  # 标题加粗并居中
    legend.title = element_text(size = 14, face = "bold"),   # 图例标题加粗
    legend.text = element_text(size = 13)                    # 图例文字
  )

print(p4)


ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/关键微生物与免疫治疗效果关系/mrr-4-1-4_Capnocytophaga_NLR.pdf",p4, width = 6, height =5)

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/关键微生物与免疫治疗效果关系/mrr-4-1-4_saliva_normalised_with_metadata_Capnocytophaga_NLR.pdf",p4, width = 6, height =5)



library(ggpubr)

# # 去掉NA
# df_clean <- df[complete.cases(df$RECIST, df$Capnocytophaga), ]
# 
# # 设定明确的分组顺序（因子 levels 必须与你的 comparisons 一致）
# df_clean$RECIST <- factor(df_clean$RECIST, levels = c("CR/PR", "SD", "PD"))

df$RECIST <- factor(df$RECIST, levels = c("SD", "PR", "PD"))

# 设置两两比较组
my_comparisons <- list(
  c("SD", "PR"),
  c("SD", "PD"),
  c("PR", "PD")
)


p <- ggboxplot(df, x = "RECIST", y = "Capnocytophaga", fill = "RECIST", 
               palette = "jco", add = "jitter", size = 1.2) +  # 添加抖动点
  stat_compare_means(comparisons = my_comparisons, 
                     method = "wilcox.test", 
                     label = "p.signif",  # 显著性符号
                     size = 5) +  # 显著性标签字体大小
  labs(
    title = "Capnocytophaga vs RECIST", 
    x = "RECIST Response",
    y = "Capnocytophaga abundance"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    legend.position = "none"
  )

# 显示图形
print(p)

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/关键微生物与免疫治疗效果关系/mrr-4-1-4_saliva_normalised_with_metadata_Capnocytophaga_RECIST.pdf",p1, width = 6, height =5)


# 多图排版
ggarrange(p1, p2, p3, p4,  
          ncol = 2, nrow = 3, 
          labels = c("A", "B", "C", "D", "E"))


####读取tsv---

data <- read.table("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/关键微生物与免疫治疗效果关系/mrr-4-1-4-LefSe analysis KELLY saliva raw reads filtered.tsv", header=T,sep = "\t", stringsAsFactors = FALSE)

df<-data

View(df)

# 加载必要的包
library(ggplot2)
library(ggpubr)

# 假设你的数据框名称为 df，包含微生物丰度（Abundance）和分组信息（CBR）
# 例如：
# df <- data.frame(Abundance = c(...), CBR = c("CB", "no-CB", ...))

View(df)


df_long <- df %>%
  pivot_longer(cols = -CBR,          # 选择除了Species列的所有列
               names_to = "group",         # 将列名转为“CBR”
               values_to = "Abundance") # 将值转为“Abundance”

df_long$group <- gsub("CB\\..*", "CB", df_long$group)


df_capnocytophaga <- df_long %>%
  filter(CBR == "Capnocytophaga")

# 绘制箱形图并加入显著性差异标记
library(ggplot2)
library(ggpubr)

pp<-ggplot(df_capnocytophaga, aes(x = group, y = Abundance, fill = group)) +
  geom_boxplot() +
  stat_compare_means(
    comparisons = list(c("CB", "no.CB")),  # 注意：group 中是 "no.CB" 不是 "no-CB"
    method = "wilcox.test",
    label = "p.signif",
    size = 5
  ) +
  labs(
    title = "Capnocytophaga Abundance Between Groups",
    x = "Group",
    y = "Capnocytophaga Abundance"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5)
  )

ggsave("D:/wangys_uestc_data2/OSCC_data/oscc_process/MPR_NMPR_group/预后分析/关键微生物与免疫治疗效果关系/mrr-4-1-4-Capnocytophaga_Abundance_Between_Groups.pdf", pp,width=7,height = 6 )

