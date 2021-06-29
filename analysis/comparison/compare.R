

# Loading libraries

library(ggplot2)
library(gridExtra)
library(ggpubr)
library(grid)
library(magrittr)
library(dplyr)
library(tidyr)


if (dir.exists("/hpc/cuppen/")){
  
  base_dir <- "/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/"
  devtools::load_all(paste0(base_dir,'/CHORD/processed/scripts_main/mutSigExtractor/'))
  
} else {
  
  library(mutSigExtractor)
}


source("/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/themes.R")


# PCAWG data loading


if (dir.exists("/hpc/cuppen/")){
  pcawg_meta <- readRDS(file = "/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/datasets/pcawg/metadata-per-donor.Rds")
} else {
  pcawg_meta <- readRDS(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/datasets/pcawg/metadata-per-donor.Rds")
}

# head(pcawg_lung_meta)

pcawg_lung_meta <- pcawg_meta[pcawg_meta$organ_system == "LUNG & BRONCHUS",]
<<<<<<< HEAD
pcawg_lung_meta <- pcawg_lung_meta[pcawg_lung_meta$icgc_donor_id != "DO23717",]
=======
>>>>>>> d4b2ba7d869eab81f62853d96bcbc0ac881b0e01



if (dir.exists("/hpc/cuppen/")){
  pcawg_lung_somatic_vcfs <- readRDS(file = "/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/pcawg/list-of-pcawg-lung-somatic-vcfs.rds")
} else {
  pcawg_lung_somatic_vcfs <- readRDS(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/pcawg/toy-list-of-pcawg-lung-somatic-vcfs.rds")
}




if (dir.exists("/hpc/cuppen/")){
  count_df_pcawg <- readRDS(file = "/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/pcawg/count_df_pcawg.rds")
} else {
  count_df_pcawg <- readRDS(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/pcawg/count_df_pcawg.rds")
}


if (dir.exists("/hpc/cuppen/")){
  sv_count_matrix_pcawg <- readRDS(file = "/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/pcawg/simple_sv_count_df_pcawg.rds")
} else {
  sv_count_matrix_pcawg <- readRDS(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/pcawg/simple_sv_count_df_pcawg.rds")
}



if (dir.exists("/hpc/cuppen/")){
  clonality_summary_pcawg <- readRDS(file = "/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/pcawg/clonality_summary_pcawg.rds")
} else {
  clonality_summary_pcawg <- readRDS(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/pcawg/clonality_summary_pcawg.rds")
}

<<<<<<< HEAD
=======

clonality_summary_pcawg$Tot <- rowSums(clonality_summary_pcawg[,3:5], na.rm = T)
clonality_summary_pcawg %<>% gather(key = "Clonality_stat", value = "Counts", 3:5)
clonality_summary_pcawg$Mut_type <- factor(clonality_summary_pcawg$Mut_type, levels = c("All", "SNVs", "MNVs", "Ins", "Del"))
clonality_summary_pcawg$Clonality_stat <- factor(clonality_summary_pcawg$Clonality_stat, levels = c("Clonal", "Non_clonal", "Unknown"))
>>>>>>> d4b2ba7d869eab81f62853d96bcbc0ac881b0e01

clonality_summary_pcawg$Tot <- rowSums(clonality_summary_pcawg[,3:5], na.rm = T)
clonality_summary_pcawg %<>% gather(key = "Clonality_stat", value = "Counts", 3:5)
clonality_summary_pcawg$Mut_type <- factor(clonality_summary_pcawg$Mut_type, levels = c("All", "SNVs", "MNVs", "Ins", "Del"))
clonality_summary_pcawg$Clonality_stat <- factor(clonality_summary_pcawg$Clonality_stat, levels = c("Clonal", "Non_clonal", "Unknown"))


# Ploidy info

# Ploidy info

if (dir.exists("/hpc/cuppen/")){
  pcawg_ploidy_info <- readRDS(file = "/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/pcawg/pcawg-ploidy-info.rds")
} else {
  pcawg_ploidy_info <- readRDS(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/pcawg/pcawg-ploidy-info.rds")
}


pcawg_ploidy_matrix <- matrix(unlist(pcawg_ploidy_info), ncol = ncol(pcawg_ploidy_info[[1]]), byrow = T)
colnames(pcawg_ploidy_matrix) <- colnames(pcawg_ploidy_info[[1]])
pcawg_ploidy_df <- as.data.frame(pcawg_ploidy_matrix)
pcawg_ploidy_df %<>% mutate(sampleId = names(pcawg_ploidy_info))
pcawg_ploidy_df <- pcawg_ploidy_df[,c(25, 1:24)]


for (i in c(2:6, 9:17, 19, 23, 25)) {
  pcawg_ploidy_df[,i] <- as.numeric(pcawg_ploidy_df[,i])
  
}

for (i in c(21, 25)) {
  pcawg_ploidy_df[,i] <- as.integer(pcawg_ploidy_df[,i])
  
}


# HMF data loading



if (dir.exists("/hpc/cuppen/")){
  hmf_meta <- read.csv(file = "/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/datasets/hmf/hmf-lung-cohort.tsv", header = T,
                       sep = "\t", stringsAsFactors = F)
} else {
  hmf_meta <- read.csv(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/datasets/hmf/hmf-lung-cohort.tsv", header = T,
                       sep = "\t", stringsAsFactors = F)
}




if (dir.exists("/hpc/cuppen/")){
  hmf_lung_meta <- read.csv(file = "/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/datasets/hmf/hmf-lung-cohort-final.tsv", header = T,
                            sep = "\t", stringsAsFactors = F)
} else {
  hmf_lung_meta <- read.csv(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/datasets/hmf/hmf-lung-cohort-final.tsv", header = T,
                            sep = "\t", stringsAsFactors = F)
}

hmf_lung_meta_NSC <- hmf_lung_meta[hmf_lung_meta$cancer_type == "Non-Small Cell",]




if (dir.exists("/hpc/cuppen/")){
  hmf_lung_somatic_vcfs <- readRDS(file = "/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/hmf/list-of-hmf-lung-somatic-vcfs.rds")
} else {
  hmf_lung_somatic_vcfs <- readRDS(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/hmf/toy-list-of-hmf-lung-somatic-vcfs.rds")
}



if (dir.exists("/hpc/cuppen/")){
  count_df_hmf <- readRDS(file = "/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/hmf/count_df_hmf.rds")
} else {
  count_df_hmf <- readRDS(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/hmf/count_df_hmf.rds")
}

count_df_hmf_NSC <- count_df_hmf[count_df_hmf$sampleId %in% hmf_lung_meta_NSC$sampleId,]




if (dir.exists("/hpc/cuppen/")){
  sv_count_matrix_hmf <- readRDS(file = "/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/hmf/simple_sv_count_df_hmf.rds")
} else {
  sv_count_matrix_hmf <- readRDS(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/hmf/simple_sv_count_df_hmf.rds")
}


sv_count_matrix_hmf_NSC <- sv_count_matrix_hmf[rownames(sv_count_matrix_hmf) %in% hmf_lung_meta_NSC$sampleId,]



if (dir.exists("/hpc/cuppen/")){
  clonality_summary_hmf <- readRDS(file = "/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/hmf/clonality_summary_hmf.rds")
} else {
  clonality_summary_hmf <- readRDS(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/hmf/clonality_summary_hmf.rds")
}

clonality_summary_hmf_NSC <- clonality_summary_hmf[clonality_summary_hmf$sampleId %in% hmf_lung_meta_NSC$sampleId,]


clonality_summary_hmf_NSC$Tot <- rowSums(clonality_summary_hmf_NSC[,3:5], na.rm = T)
clonality_summary_hmf_NSC %<>% gather(key = "Clonality_stat", value = "Counts", 3:5)
clonality_summary_hmf_NSC$Mut_type <- factor(clonality_summary_hmf_NSC$Mut_type, levels = c("All", "SNVs", "MNVs", "Ins", "Del"))
clonality_summary_hmf_NSC$Clonality_stat <- factor(clonality_summary_hmf_NSC$Clonality_stat, levels = c("Clonal", "Non_clonal", "Unknown"))




# Ploidy info

if (dir.exists("/hpc/cuppen/")){
  hmf_ploidy_info <- readRDS(file = "/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/hmf/hmf-ploidy-info.rds")
} else {
  hmf_ploidy_info <- readRDS(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/hmf/hmf-ploidy-info.rds")
}

hmf_ploidy_matrix <- matrix(unlist(hmf_ploidy_info), ncol = ncol(hmf_ploidy_info[[1]]), byrow = T)
colnames(hmf_ploidy_matrix) <- colnames(hmf_ploidy_info[[1]])
hmf_ploidy_df <- as.data.frame(hmf_ploidy_matrix)
hmf_ploidy_df %<>% mutate(sampleId = names(hmf_ploidy_info))
hmf_ploidy_df <- hmf_ploidy_df[,c(25, 1:24)]


for (i in c(2:6, 9:17, 19, 23, 25)) {
  hmf_ploidy_df[,i] <- as.numeric(hmf_ploidy_df[,i])
  
}

for (i in c(21, 25)) {
  hmf_ploidy_df[,i] <- as.integer(hmf_ploidy_df[,i])
  
}

hmf_ploidy_df_NSC <- hmf_ploidy_df[hmf_ploidy_df$sampleId %in% hmf_lung_meta_NSC$sampleId,]




############################################################################################################################################

# Boxplots for somatic SNVs, MNVs, and small insetions and deletions



options(scipen=999) 


### Using R base boxplot function


boxplot(count_df_pcawg$SNVs, count_df_hmf$SNVs, names = c("Primary Lung Cancer", "Metastatic Lung Cancer"))
hmf_snv_mean <- round(mean(count_df_hmf$SNVs), 2)
pcawg_snv_mean <- round(mean(count_df_pcawg$SNVs), 2)

text(x= 0.7, y= 200000, labels= paste0("Mean value of ", pcawg_snv_mean), font = 2, col = "blue")
text(x= 1.7, y= 200000, labels= paste0("Mean value of ", hmf_snv_mean), font = 2, col = "red")



### Using ggplot2


# prepare the data

SNVs_plot_df1 <- data.frame(label = "Primary Lung Cohort", value = count_df_pcawg$SNVs)
SNVs_plot_df2 <- data.frame(label = "Metastatic Lung Cohort", value = count_df_hmf$SNVs)
SNVs_plot_df <- rbind(SNVs_plot_df1, SNVs_plot_df2)
SNVs_plot_df$label <- factor(SNVs_plot_df$label, levels = c("Primary Lung Cohort", "Metastatic Lung Cohort"))




# plotting

meds_snv <- c(by(SNVs_plot_df$value, SNVs_plot_df$label, median))
meds_snv[[1]] <- -50000
meds_snv[[2]] <- -50000

hmf_snv_mean <- round(mean(count_df_hmf$SNVs), 2)
pcawg_snv_mean <- round(mean(count_df_pcawg$SNVs), 2)


my_comparisons <- list(c("Primary Lung Cohort", "Metastatic Lung Cohort"))


SNVs_plot <- ggplot(data = SNVs_plot_df, aes(x = label, y = value, color = label)) +
  geom_boxplot(outlier.colour="red", outlier.shape=17,
               outlier.size=2, outlier.alpha = 0.5) +
  stat_compare_means(comparisons = my_comparisons, method = "t.test") +
  geom_text(data=data.frame(), aes(x=names(meds_snv), y=meds_snv,
                                   label= c(paste0("Mean of ", pcawg_snv_mean), paste0("Mean of ", hmf_snv_mean))), col = "red", size=4) +
  # scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  scale_color_manual(values = c("blue", "green")) +
  geom_jitter(color="black", size=0.4, alpha=0.5) +
  stat_boxplot(geom ='errorbar') +
  theme_bw() +
  ggtitle("Numbers of somatic SNVs \n") +
  ylab("Number of SNVs \n") +
  xlab("Cohorts \n") +
  theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.5)) +
  theme(axis.title.x = element_text(face = "italic", size = 12)) +
  theme(axis.title.y = element_text(face = "italic", size = 12)) +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))

SNVs_plot


### Using gghalves to create look-alike and consitent figures:

SNVs_plot_df %<>% group_by(label) %>%
  mutate(median = median(value, na.rm = T)) %>%
  mutate(mean = round(mean(value, na.rm = T), digits = 2)) %>%
  ungroup()



SNVs_plot_job <- ggplot(data = SNVs_plot_df, aes(x = label, y = log2(value), color = label, fill = label, label = median)) +
  geomJob_HalfHalfBox +
  ggplot2::scale_color_manual(values = c('#e75480', 'darkblue'), guide = "none") +
  ggplot2::scale_fill_manual(values = c('#e75480', 'darkblue'), guide = "none") +
  # ggplot2::scale_y_continuous(limits = c(20, 90), breaks = seq(0, 90, by = 10)) +
  ggplot2::ggtitle("Numbers of Somatic SNVs \n") +
  ggplot2::labs(x = '\n\n Cohorts', y = 'Number of SNVs (log2)\n\n') +
  theme_Job +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") +
  theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.5))
  # theme(axis.title.x = element_text(face = "italic", size = 12))
  # theme(axis.title.y = element_text(face = "italic", size = 12))
  # theme(plot.margin = unit(c(1,1,1,1), "cm"))
  # geom_text(data=data.frame(), aes(x=names(meds_snv), y=meds_snv,
  #                                  label= c(paste0("Mean = ", pcawg_snv_mean), paste0("Mean = ", hmf_snv_mean))), col = "red", size=4) +





# ================================================================================================

### Using R basic graphics


boxplot(count_df_pcawg$MNVs, count_df_hmf$MNVs, names = c("Primary Lung Cancer", "Metastatic Lung Cancer"))
hmf_mnv_mean <- round(mean(count_df_hmf$MNVs), 2)
pcawg_mnv_mean <- round(mean(count_df_pcawg$MNVs), 2)

text(x= 0.7, y= 4000, labels= paste0("Mean value of ", pcawg_mnv_mean), font = 2, col = "blue")
text(x= 1.7, y= 4000, labels= paste0("Mean value of ", hmf_mnv_mean), font = 2, col = "red")


### Using ggplot2


# prepare the data

MNVs_plot_df1 <- data.frame(label = "Primary Lung Cohort", value = count_df_pcawg$MNVs)
MNVs_plot_df2 <- data.frame(label = "Metastatic Lung Cohort", value = count_df_hmf$MNVs)
MNVs_plot_df <- rbind(MNVs_plot_df1, MNVs_plot_df2)
MNVs_plot_df$label <- factor(MNVs_plot_df$label, levels = c("Primary Lung Cohort", "Metastatic Lung Cohort"))

# plotting

meds_mnv <- c(by(MNVs_plot_df$value, MNVs_plot_df$label, median))
meds_mnv[[1]] <- -2000
meds_mnv[[2]] <- -2000

hmf_mnv_mean <- round(mean(count_df_hmf$MNVs), 2)
pcawg_mnv_mean <- round(mean(count_df_pcawg$MNVs), 2)

MNVs_plot <- ggplot(data = MNVs_plot_df, aes(x = label, y = value, color = label)) +
  geom_boxplot(outlier.colour="red", outlier.shape=17,
               outlier.size=2, outlier.alpha = 0.5) +
  stat_compare_means(comparisons = my_comparisons, method = "t.test") +
  geom_text(data=data.frame(), aes(x=names(meds_mnv), y=meds_mnv, 
                                   label=c(paste0("Mean of ", pcawg_mnv_mean), paste0("Mean of ", hmf_mnv_mean))), col = "red", size=4) +
  scale_color_manual(values = c("blue", "green")) +
  geom_jitter(color="black", size=0.4, alpha=0.5) +
  theme_bw() +
  ggtitle("Numbers of somatic MNVs \n") +
  ylab("Number of MNVs \n") +
  xlab("Cohorts \n") +
  theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.5)) +
  theme(axis.title.x = element_text(face = "italic", size = 12)) +
  theme(axis.title.y = element_text(face = "italic", size = 12)) +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))

MNVs_plot




### Using gghalves to create look-alike and consitent figures:

MNVs_plot_df %<>% group_by(label) %>%
  mutate(median = median(value, na.rm = T)) %>%
  mutate(mean = round(mean(value, na.rm = T), digits = 2)) %>%
  ungroup()



MNVs_plot_job <- ggplot(data = MNVs_plot_df, aes(x = label, y = log2(value), color = label, fill = label, label = median)) +
  geomJob_HalfHalfBox +
  ggplot2::scale_color_manual(values = c('#e75480', 'darkblue'), guide = "none") +
  ggplot2::scale_fill_manual(values = c('#e75480', 'darkblue'), guide = "none") +
  # ggplot2::scale_y_continuous(limits = c(20, 90), breaks = seq(0, 90, by = 10)) +
  ggplot2::ggtitle("Numbers of Somatic MNVs \n") +
  ggplot2::labs(x = '\n\n Cohorts', y = 'Number of MNVs (log2)\n\n') +
  theme_Job +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") +
  theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.5))
# theme(axis.title.x = element_text(face = "italic", size = 12))
# theme(axis.title.y = element_text(face = "italic", size = 12))
# theme(plot.margin = unit(c(1,1,1,1), "cm"))
# geom_text(data=data.frame(), aes(x=names(meds_mnv), y=meds_mnv,
#                                  label= c(paste0("Mean = ", pcawg_mnv_mean), paste0("Mean = ", hmf_mnv_mean))), col = "red", size=4) +






# ====================================================================================

### Using R basic graphics

boxplot(log10(count_df_pcawg$inDEL), log10(count_df_hmf$inDEL), names = c("Primary Lung Cancer", "Metastatic Lung Cancer"),
        main = "Number of small deletions \n \n \n", xlab = "Cohorts", ylab = "Log10 number of deletions", )
hmf_del_mean <- round(mean(count_df_hmf$inDEL), 2)
pcawg_del_mean <- round(mean(count_df_pcawg$inDEL), 2)

text(x= 1, y= 3.65, labels= paste0("Mean value of ", pcawg_del_mean), font = 2, col = "blue")
text(x= 2, y= 4.05, labels= paste0("Mean value of ", hmf_del_mean), font = 2, col = "red")

### Using ggplot2


# prepare the data

inDEL_plot_df1 <- data.frame(label = "Primary Lung Cohort", value = log2(count_df_pcawg$inDEL))
inDEL_plot_df2 <- data.frame(label = "Metastatic Lung Cohort", value = log2(count_df_hmf$inDEL))
inDEL_plot_df <- rbind(inDEL_plot_df1, inDEL_plot_df2)
inDEL_plot_df$label <- factor(inDEL_plot_df$label, levels = c("Primary Lung Cohort", "Metastatic Lung Cohort"))


# plotting

meds_del <- c(by(inDEL_plot_df$value, inDEL_plot_df$label, median))
meds_del[[1]] <- 0
meds_del[[2]] <- 0

hmf_del_mean <- round(mean(count_df_hmf$inDEL), 2)
pcawg_del_mean <- round(mean(count_df_pcawg$inDEL), 2)

inDEL_plot <- ggplot(data = inDEL_plot_df, aes(x = label, y = value, color = label)) +
  geom_boxplot(outlier.colour="red", outlier.shape=17,
               outlier.size=2, outlier.alpha = 0.5) +
  stat_compare_means(method = "t.test") +
  # stat_compare_means(data = inDEL_plot_df, aes(x = label, y = value),comparisons = my_comparisons, method = "t.test", label.y = 2) +
  geom_text(data=data.frame(), aes(x=names(meds_del), y=meds_del, 
                                   label=c(paste0("Mean of ", pcawg_del_mean), paste0("Mean of ", hmf_del_mean))), col = "red", size=4) +
  scale_color_manual(values = c("blue", "green")) +
  geom_jitter(color="black", size=0.4, alpha=0.5) +
  theme_bw() +
  ggtitle("Numbers of somatic small deletions \n") +
  ylab("Number of deletions (log2) \n") +
  xlab("Cohorts \n") +
  theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.5)) +
  theme(axis.title.x = element_text(face = "italic", size = 12)) +
  theme(axis.title.y = element_text(face = "italic", size = 12)) +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))

inDEL_plot




### Using gghalves to create look-alike and consitent figures:

inDEL_plot_df %<>% group_by(label) %>%
  mutate(median = round(median(value, na.rm = T), digits = 1)) %>%
  mutate(mean = round(mean(value, na.rm = T), digits = 2)) %>%
  ungroup()



inDEL_plot_job <- ggplot(data = inDEL_plot_df, aes(x = label, y = log2(value), color = label, fill = label, label = median)) +
  geomJob_HalfHalfBox +
  ggplot2::scale_color_manual(values = c('#e75480', 'darkblue'), guide = "none") +
  ggplot2::scale_fill_manual(values = c('#e75480', 'darkblue'), guide = "none") +
  # ggplot2::scale_y_continuous(limits = c(20, 90), breaks = seq(0, 90, by = 10)) +
  ggplot2::ggtitle("Numbers of Somatic Deletions \n") +
  ggplot2::labs(x = '\n\n Cohorts', y = 'Number of small deletions (log2)\n\n') +
  theme_Job +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") +
  theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.5))
# theme(axis.title.x = element_text(face = "italic", size = 12))
# theme(axis.title.y = element_text(face = "italic", size = 12))
# theme(plot.margin = unit(c(1,1,1,1), "cm"))
# geom_text(data=data.frame(), aes(x=names(meds_mnv), y=meds_mnv,
#                                  label= c(paste0("Mean = ", pcawg_mnv_mean), paste0("Mean = ", hmf_mnv_mean))), col = "red", size=4) +




# =====================================================================================================================================


### Using R basic graphics

boxplot(log10(count_df_pcawg$INdel), log10(count_df_hmf$INdel), names = c("Primary Lung Cancer", "Metastatic Lung Cancer"), main = "Number of small insertions", xlab = "Cohorts", ylab = "Log10 number of insertions")
pcawg_ins_mean <- round(mean(count_df_pcawg$INdel), 2)
hmf_ins_mean <- round(mean(count_df_hmf$INdel), 2)

text(x= 0.7, y= 200000, labels= paste0("Mean value of ", pcawg_ins_mean), font = 2, col = "blue")
text(x= 1.7, y= 200000, labels= paste0("Mean value of ", hmf_ins_mean), font = 2, col = "red")





### Using ggplot2


# prepare the data

INdel_plot_df1 <- data.frame(label = "Primary Lung Cohort", value = log2(count_df_pcawg$INdel))
INdel_plot_df2 <- data.frame(label = "Metastatic Lung Cohort", value = log2(count_df_hmf$INdel))
INdel_plot_df <- rbind(INdel_plot_df1, INdel_plot_df2)
INdel_plot_df$label <- factor(INdel_plot_df$label, levels = c("Primary Lung Cohort", "Metastatic Lung Cohort"))

# plotting

meds_ins <- c(by(INdel_plot_df$value, INdel_plot_df$label, median))
meds_ins[[1]] <- 0
meds_ins[[2]] <- 0

hmf_ins_mean <- round(mean(count_df_hmf$INdel), 2)
pcawg_ins_mean <- round(mean(count_df_pcawg$INdel), 2)

INdel_plot <- ggplot(data = INdel_plot_df, aes(x = label, y = value, color = label)) +
  geom_boxplot(outlier.colour="red", outlier.shape=17,
               outlier.size=2, outlier.alpha = 0.5) +
  stat_compare_means(method = "t.test") +
  geom_text(data=data.frame(), aes(x=names(meds_ins), y=meds_ins, 
                                   label=c(paste0("Mean of ", pcawg_ins_mean), paste0("Mean of ", hmf_ins_mean))), col = "red", size=4) +
  scale_color_manual(values = c("blue", "green")) +
  geom_jitter(color="black", size=0.4, alpha=0.5) +
  theme_bw() +
  ggtitle("Numbers of somatic small insertions \n") +
  ylab("Number of insertions (log2) \n") +
  xlab("Cohorts \n") +
  theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.5)) +
  theme(axis.title.x = element_text(face = "italic", size = 12)) +
  theme(axis.title.y = element_text(face = "italic", size = 12)) +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))

INdel_plot



### Using gghalves to create look-alike and consitent figures:

INdel_plot_df %<>% group_by(label) %>%
  mutate(median = round(median(value, na.rm = T), digits =1)) %>%
  mutate(mean = round(mean(value, na.rm = T), digits = 2)) %>%
  ungroup()



INdel_plot_job <- ggplot(data = INdel_plot_df, aes(x = label, y = log2(value), color = label, fill = label, label = median)) +
  geomJob_HalfHalfBox +
  ggplot2::scale_color_manual(values = c('#e75480', 'darkblue'), guide = "none") +
  ggplot2::scale_fill_manual(values = c('#e75480', 'darkblue'), guide = "none") +
  # ggplot2::scale_y_continuous(limits = c(20, 90), breaks = seq(0, 90, by = 10)) +
  ggplot2::ggtitle("Numbers of Somatic Insertions \n") +
  ggplot2::labs(x = '\n\n Cohorts', y = 'Number of small insertions (log2)\n\n') +
  theme_Job +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") +
  theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.5))
# theme(axis.title.x = element_text(face = "italic", size = 12))
# theme(axis.title.y = element_text(face = "italic", size = 12))
# theme(plot.margin = unit(c(1,1,1,1), "cm"))
# geom_text(data=data.frame(), aes(x=names(meds_mnv), y=meds_mnv,
#                                  label= c(paste0("Mean = ", pcawg_mnv_mean), paste0("Mean = ", hmf_mnv_mean))), col = "red", size=4) +





### combining all boxplots into one figure


tmb_plot_all <- ggarrange(SNVs_plot_job, MNVs_plot_job, inDEL_plot_job, INdel_plot_job, nrow = 2, ncol = 2, common.legend = TRUE, legend = "bottom")

for (i in 1:2){
  if (i == 1) {
    png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/comparison/figs/png/tmb-boxplots-final.png",
        width = 960, heigh = 960)
    print(tmb_plot_all)
    dev.off()
  }
  if (i == 2) {
    pdf(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/comparison/figs/pdf/tmb-boxplots-final.pdf",
        width = 15, heigh = 15, onefile=FALSE)
    print(tmb_plot_all)
    dev.off()
  }
}

############################################################################################################################################


# Comparing the simple SV events between the primary and metastatic datasets in lung cohorts


## Prepare the count dataframes
sv_count_df_pcawg <- as.data.frame(sv_count_matrix_pcawg)
sv_count_df_pcawg %<>% mutate(sampleId = rownames(sv_count_df_pcawg))
rownames(sv_count_df_pcawg) <- 1:nrow(sv_count_df_pcawg)
sv_count_df_pcawg <- sv_count_df_pcawg[,c(6,1:5)]


sv_count_df_hmf_NSC <- as.data.frame(sv_count_matrix_hmf_NSC)
sv_count_df_hmf_NSC %<>% mutate(sampleId = rownames(sv_count_df_hmf_NSC))
rownames(sv_count_df_hmf_NSC) <- 1:nrow(sv_count_df_hmf_NSC)
sv_count_df_hmf_NSC <- sv_count_df_hmf_NSC[,c(6,1:5)]



sv_count_tibb_pcawg <- gather(sv_count_df_pcawg, key = "Simple_SV_Event", value = "Counts", 2:6)
sv_count_tibb_pcawg %<>% mutate(Dataset = "Primary")
sv_count_tibb_hmf_NSC <- gather(sv_count_df_hmf_NSC, key = "Simple_SV_Event", value = "Counts", 2:6)
sv_count_tibb_hmf_NSC %<>% mutate(Dataset = "Metastatic")



sv_count_df_combined <- rbind(sv_count_tibb_pcawg, sv_count_tibb_hmf_NSC)


simple_sv_events <- c("Deletion (> 100kb)", "Deletion (< 100kb)", "Duplication (> 100kb)", "Duplication (< 100kb)", "Inversion")

sv_count_df_combined$Simple_SV_Event <- factor(sv_count_df_combined$Simple_SV_Event, levels = simple_sv_events)
sv_count_df_combined$Dataset <- factor(sv_count_df_combined$Dataset, levels = c("Primary", "Metastatic"))


sv_count_df_combined <- sv_count_df_combined[sv_count_df_combined$Counts != 0,]


<<<<<<< HEAD
stat.test <- sv_count_df_combined %>%
  group_by(Simple_SV_Event) %>%
  wilcox_test(Counts ~ Dataset) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")



stat.test <- stat.test %>%
  add_xy_position(x = "dbs_sig", dodge = 0.8)
stat.test$y.position <- log2(stat.test$y.position)
stat.test$xmin <- c(0.85,1.85,2.85,3.85,4.85)
stat.test$xmax <- c(1.15,2.15,3.15,4.15,5.15)



sv_plot <- sv_count_df_combined %>% ggplot(aes(x = Simple_SV_Event, y = log2(Counts), color = Dataset)) +
=======
sv_plot <- sv_count_df_combined %>% ggplot(aes(x = Simple_SV_Event, y = log2(Counts), fill = Dataset)) +
>>>>>>> d4b2ba7d869eab81f62853d96bcbc0ac881b0e01
  geom_boxplot(outlier.shape = NA) +
  # stat_compare_means(method = "t.test") +
  # geom_text(data=data.frame(), aes(x=names(meds_ins), y=meds_ins, 
  #                                  label=c(paste0("Mean of ", pcawg_ins_mean), paste0("Mean of ", hmf_ins_mean))), col = "red", size=4) +
  # scale_color_manual(values = c("blue", "green")) +
  # geom_jitter(color="black", size=0.4, alpha=0.5) +
  geom_point(position=position_jitterdodge(), size = 0.4, aes(color = Dataset), alpha = 0.5) +
  scale_color_manual(values = c("#ff0101", "#010dff"), guide = "none") +
<<<<<<< HEAD
  stat_pvalue_manual(stat.test,  label = "p.adj.signif", tip.length = 0, hide.ns = F) +
=======
>>>>>>> d4b2ba7d869eab81f62853d96bcbc0ac881b0e01
  theme_bw() +
  ggtitle("SV Event Frequency Comparison \n") +
  ylab("Frequency of Events (log2) \n") +
  ylim(0, 12) +
  xlab("Simple SV Events \n") +
  theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.5)) +
  theme(axis.title.x = element_text(face = "italic", size = 12)) +
  theme(axis.title.y = element_text(face = "italic", size = 12)) +
  labs(fill = "Cancer Type") +
  theme(axis.text.x = element_text(angle = 45, size = 5, vjust = 0.6, face = "bold")) +
<<<<<<< HEAD
  theme(plot.margin = unit(c(1,1,1,1), "cm")) 
=======
  theme(plot.margin = unit(c(1,1,1,1), "cm"))
>>>>>>> d4b2ba7d869eab81f62853d96bcbc0ac881b0e01

sv_plot




for (i in 1:2){
  if (i == 1) {
<<<<<<< HEAD
    png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/comparison/figs/png/sv-boxplots-final.png")
=======
    png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/comparison/figs/png/sv-boxplots.png")
>>>>>>> d4b2ba7d869eab81f62853d96bcbc0ac881b0e01
    print(sv_plot)
    dev.off()
  }
  if (i == 2) {
<<<<<<< HEAD
    pdf(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/comparison/figs/pdf/sv-boxplots-final.pdf")
=======
    pdf(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/comparison/figs/pdf/sv-boxplots.pdf")
>>>>>>> d4b2ba7d869eab81f62853d96bcbc0ac881b0e01
    print(sv_plot)
    dev.off()
  }
}


<<<<<<< HEAD


##### 
# Job theme

sv_count_df_combined %<>% group_by("Simple_SV_Event") %>%
  mutate(median = round(median(Counts, na.rm = T), digits = 1)) %>%
  ungroup


sv_plot_job <- ggplot(data = sv_count_df_combined, aes(x = Dataset, y = log2(Counts), color = Dataset, fill = Dataset)) + facet_wrap(~ Simple_SV_Event) +
  gghalves::geom_half_boxplot(side = 'l', outlier.shape = NA, color = 'black') +
  gghalves::geom_half_point_panel(side = 'r', color = 'black', shape = 21, transformation = ggbeeswarm::position_quasirandom(width = .1, groupOnX = T)) +
  ggplot2::scale_color_manual(values = c('#e75480', 'darkblue'), guide = "none") +
  ggplot2::scale_fill_manual(values = c('#e75480', 'darkblue'), guide = "none") +
  # ggplot2::scale_y_continuous(limits = c(20, 90), breaks = seq(0, 90, by = 10)) +
  ggplot2::ggtitle("Numbers of Somatic Insertions \n") +
  ggplot2::labs(x = '\n\n Cohorts', y = 'Number of small insertions (log2)\n\n') +
  theme_Job +
  # stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") +
  theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.5))
# theme(axis.title.x = element_text(face = "italic", size = 12))
# theme(axis.title.y = element_text(face = "italic", size = 12))
# theme(plot.margin = unit(c(1,1,1,1), "cm"))
# geom_text(data=data.frame(), aes(x=names(meds_mnv), y=meds_mnv,
#                                  label= c(paste0("Mean = ", pcawg_mnv_mean), paste0("Mean = ", hmf_mnv_mean))), col = "red", size=4) +


for (i in 1:2){
  if (i == 1) {
    png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/comparison/figs/png/sv-boxplots-job.png", width = 960, height = 960)
    print(sv_plot_job)
    dev.off()
  }
  if (i == 2) {
    pdf(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/comparison/figs/pdf/sv-boxplots-job.pdf", width = 14, height = 14)
    print(sv_plot_job)
    dev.off()
  }
}


###########################################################################################################################################
# A plot for the correlation between age and TMB


head(pcawg_lung_meta)
pcawg_lung_meta_processed <- pcawg_lung_meta[,c("icgc_donor_id", "donor_sex", "donor_age_at_diagnosis")]
pcawg_lung_meta_processed %<>% mutate(Dataset = "PCAWG") %>% rename(sampleId = icgc_donor_id, donorSex = donor_sex, donorAge = donor_age_at_diagnosis)
head(count_df_pcawg)
all(pcawg_lung_meta_processed$sampleId == count_df_pcawg$sampleId)
pcawg_lung_meta_processed_final <- cbind(pcawg_lung_meta_processed, count_df_pcawg[,2:6])
nrow(pcawg_lung_meta_processed_final)


head(hmf_lung_meta)
# length(as.numeric(sapply(strsplit(hmf_lung_meta$biopsyDate, split = "-"), "[", 1)) - as.numeric(hmf_lung_meta$birthYear))
hmf_lung_meta_NSC %<>% mutate(donor_age_at_biopsy = as.numeric(sapply(strsplit(hmf_lung_meta_NSC$biopsyDate, split = "-"), "[", 1)) - as.numeric(hmf_lung_meta_NSC$birthYear))

hmf_lung_meta_processed <- hmf_lung_meta_NSC[,c("sampleId", "gender", "donor_age_at_biopsy")]
hmf_lung_meta_processed %<>% mutate(Dataset = "HMF") %>% rename(donorSex = gender, donorAge = donor_age_at_biopsy)
all(hmf_lung_meta_processed$sampleId == count_df_hmf_NSC$sampleId)
hmf_lung_meta_processed_final <- cbind(hmf_lung_meta_processed, count_df_hmf_NSC[,2:6])
nrow(hmf_lung_meta_processed_final)


lung_meta_combined <- rbind(pcawg_lung_meta_processed_final, hmf_lung_meta_processed_final)
str(lung_meta_combined)


lung_meta_combined$donorSex[lung_meta_combined$donorSex == "male"] <- "Male"
lung_meta_combined$donorSex[lung_meta_combined$donorSex == "female"] <- "Female"

lung_meta_combined$donorSex <- factor(lung_meta_combined$donorSex, levels = c("Male", "Female"))
lung_meta_combined$Dataset <- factor(lung_meta_combined$Dataset, levels = c("PCAWG", "HMF"))

lung_meta_combined_filtered <- lung_meta_combined[!is.na(lung_meta_combined$donorSex) & !is.na(lung_meta_combined$donorAge), ]

str(lung_meta_combined_filtered)
lung_meta_combined_filtered$all <- log2(lung_meta_combined_filtered$all) 

a <- lm(lung_meta_combined_filtered$all[lung_meta_combined_filtered$donorSex == "Male" & lung_meta_combined_filtered$Dataset == "PCAWG"] ~ lung_meta_combined_filtered$donorAge[lung_meta_combined_filtered$donorSex == "Male" & lung_meta_combined_filtered$Dataset == "PCAWG"])
b <- lm(lung_meta_combined_filtered$all[lung_meta_combined_filtered$donorSex == "Male" & lung_meta_combined_filtered$Dataset == "HMF"] ~ lung_meta_combined_filtered$donorAge[lung_meta_combined_filtered$donorSex == "Male" & lung_meta_combined_filtered$Dataset == "HMF"])
c <- lm(lung_meta_combined_filtered$all[lung_meta_combined_filtered$donorSex == "Female" & lung_meta_combined_filtered$Dataset == "PCAWG"] ~ lung_meta_combined_filtered$donorAge[lung_meta_combined_filtered$donorSex == "Female" & lung_meta_combined_filtered$Dataset == "PCAWG"])
d <- lm(lung_meta_combined_filtered$all[lung_meta_combined_filtered$donorSex == "Female" & lung_meta_combined_filtered$Dataset == "HMF"] ~ lung_meta_combined_filtered$donorAge[lung_meta_combined_filtered$donorSex == "Female" & lung_meta_combined_filtered$Dataset == "HMF"])


eq <- data.frame(donorSex = rep(c("Male", "Female"), each = 2), Dataset = rep(c("PCAWG", "HMF"), times = 2))

eq$Slope <- c(paste0("Slope: ", as.character(round(as.numeric(a$coefficients[2]), 3))),
              paste0("Slope: ", as.character(round(as.numeric(b$coefficients[2]), 3))),
              paste0("Slope: ", as.character(round(as.numeric(c$coefficients[2]), 3))),
              paste0("Slope: ", as.character(round(as.numeric(d$coefficients[2]), 3))))


my_plot <- lung_meta_combined_filtered %>% ggplot(aes(x = donorAge, y = all, color = Dataset)) +
  geom_point(aes(shape = donorSex), size = 2, alpha = 0.75) +
  scale_color_manual(values = c("#ff0101","#0124ff")) +
  gghighlight() +
  geom_text(data=eq,aes(x = 45, y = 9,label=Slope), parse = TRUE, inherit.aes=FALSE, size = 2) +
  facet_wrap(~ donorSex + Dataset) +
  geom_smooth(method = "lm", se=FALSE, aes(color=Dataset), formula = y ~ x) +
  theme_bw() +
  theme(strip.text = element_text(size=10, face="bold", color="black")) +
  theme(strip.background = element_rect(fill="grey", size=1, color="black")) +
  labs(x = "Patient Age", y = "TMB (log2)", shape = "Patient Gender") +
  ggtitle("TMB of PCAWG and HMF Lung Samples by Patient Age \n") +
  theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.35)) +
  theme(plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm"))





for (i in 1:2){
  if (i == 1) {
    png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/comparison/figs/png/age-tmb-final.png")
    print(my_plot)
    dev.off()
  }
  if (i == 2) {
    pdf(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/comparison/figs/pdf/age-tmb-final.pdf")
    print(my_plot)
    dev.off()
  }
}




###########################################################################################################################################
# Clonality

clonality_summary_pcawg$Clona2Nonclonal <- clonality_summary_pcawg$Clonal/clonality_summary_pcawg$Non_clonal
clonality_summary_pcawg_all <- clonality_summary_pcawg[clonality_summary_pcawg$Mut_type == "All",]
clonality_summary_pcawg_all <- clonality_summary_pcawg_all[clonality_summary_pcawg_all$Clonal > 10 & clonality_summary_pcawg_all$Non_clonal > 10,]

clonality_summary_pcawg_all$Dataset <- "PCAWG"




clonality_summary_hmf_NSC$Clona2Nonclonal <- clonality_summary_hmf_NSC$Clonal/clonality_summary_hmf_NSC$Non_clonal
clonality_summary_hmf_all <- clonality_summary_hmf_NSC[clonality_summary_hmf_NSC$Mut_type == "All",]
clonality_summary_hmf_all <- clonality_summary_hmf_all[clonality_summary_hmf_all$Clonal > 10 & clonality_summary_hmf_all$Non_clonal > 10,]

clonality_summary_hmf_all$Dataset <- "HMF"

=======
>>>>>>> d4b2ba7d869eab81f62853d96bcbc0ac881b0e01


##### 
# Job theme

sv_count_df_combined %<>% group_by("Simple_SV_Event") %>%
  mutate(median = round(median(Counts, na.rm = T), digits = 1)) %>%
  ungroup

<<<<<<< HEAD
stat.test <- combined_clonality %>%
  t_test(Clona2Nonclonal ~ Dataset) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")



stat.test <- stat.test %>%
  add_xy_position(x = "dbs_sig", dodge = 0.8)

stat.test$x <- 1.5
stat.test$xmin <- 1
stat.test$xmax <- 2



clonality_plot <-ggplot(combined_clonality, aes(x = Dataset, y = Clona2Nonclonal, color = Dataset)) +
  scale_color_manual(values = c("red", "blue")) +
  geom_boxplot() + 
  
  ylim(c(0,160)) +
  stat_pvalue_manual(stat.test,  label = "p.adj.signif", tip.length = 0, hide.ns = F) +
  theme_classic() +
  ggtitle("Clonality Proportion of Mutations in Lung Cohorts \n \n") +
  theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.5)) +
  labs(x = "\n Datasets", y = "Clonal/Non-clonal Ratio \n", color = "Datasets") +
  geom_hline(yintercept=c(50, 100, 150), linetype="dashed", 
             color = "grey", size=0.5)
=======
>>>>>>> d4b2ba7d869eab81f62853d96bcbc0ac881b0e01

sv_plot_job <- ggplot(data = sv_count_df_combined, aes(x = Dataset, y = log2(Counts), color = Dataset, fill = Dataset)) + facet_wrap(~ Simple_SV_Event) +
  gghalves::geom_half_boxplot(side = 'l', outlier.shape = NA, color = 'black') +
  gghalves::geom_half_point_panel(side = 'r', color = 'black', shape = 21, transformation = ggbeeswarm::position_quasirandom(width = .1, groupOnX = T)) +
  ggplot2::scale_color_manual(values = c('#e75480', 'darkblue'), guide = "none") +
  ggplot2::scale_fill_manual(values = c('#e75480', 'darkblue'), guide = "none") +
  # ggplot2::scale_y_continuous(limits = c(20, 90), breaks = seq(0, 90, by = 10)) +
  ggplot2::ggtitle("Numbers of Somatic Insertions \n") +
  ggplot2::labs(x = '\n\n Cohorts', y = 'Number of small insertions (log2)\n\n') +
  theme_Job +
  # stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") +
  theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.5))
# theme(axis.title.x = element_text(face = "italic", size = 12))
# theme(axis.title.y = element_text(face = "italic", size = 12))
# theme(plot.margin = unit(c(1,1,1,1), "cm"))
# geom_text(data=data.frame(), aes(x=names(meds_mnv), y=meds_mnv,
#                                  label= c(paste0("Mean = ", pcawg_mnv_mean), paste0("Mean = ", hmf_mnv_mean))), col = "red", size=4) +


for (i in 1:2){
  if (i == 1) {
    png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/comparison/figs/png/sv-boxplots-final.png", width = 960, height = 960)
    print(sv_plot_job)
    dev.off()
  }
  if (i == 2) {
    pdf(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/comparison/figs/pdf/sv-boxplots-final.pdf", width = 14, height = 14)
    print(sv_plot_job)
    dev.off()
  }
}


<<<<<<< HEAD
=======
############################################################################################################################################
>>>>>>> d4b2ba7d869eab81f62853d96bcbc0ac881b0e01

# clonality

<<<<<<< HEAD
#### job


combined_clonality %<>% group_by(Dataset) %>%
  mutate(median = round(median(Clona2Nonclonal, na.rm = T), digits =1)) %>%
  # mutate(mean = round(mean(value, na.rm = T), digits = 2)) %>%
  ungroup()

my_comparisons <- list(c("PCAWG", "HMF"))

clonality_plot_job <- ggplot(data = combined_clonality, aes(x = Dataset, y = Clona2Nonclonal, color = Dataset, fill = Dataset, label = median)) +
  geomJob_HalfHalfBox +
  ggplot2::scale_color_manual(values = c('#e75480', 'darkblue'), guide = "none") +
  ggplot2::scale_fill_manual(values = c('#e75480', 'darkblue'), guide = "none") +
  # ggplot2::scale_y_continuous(limits = c(20, 90), breaks = seq(0, 90, by = 10)) +
  ggplot2::ggtitle("Clonality Proportion of Mutations in Lung Cohorts \n \n") +
  ggplot2::labs(x = '\n\n Cohorts', y = 'Clonal/Non-clonal Ratio \n') +
  theme_Job +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") +
  theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.5))
# theme(axis.title.x = element_text(face = "italic", size = 12))
# theme(axis.title.y = element_text(face = "italic", size = 12))
# theme(plot.margin = unit(c(1,1,1,1), "cm"))
# geom_text(data=data.frame(), aes(x=names(meds_mnv), y=meds_mnv,
#                                  label= c(paste0("Mean = ", pcawg_mnv_mean), paste0("Mean = ", hmf_mnv_mean))), col = "red", size=4) +
=======
clonality_summary_pcawg <- clonality_summary_pcawg[clonality_summary_pcawg$Mut_type == "All" & !is.na(clonality_summary_pcawg$Counts),]
nrow(clonality_summary_pcawg)/3
clonality_summary_hmf_NSC

>>>>>>> d4b2ba7d869eab81f62853d96bcbc0ac881b0e01


boxplot(clonality_summary_pcawg$Clonal[clonality_summary_pcawg$Mut_type == "All"], clonality_summary_hmf_NSC$Clonal[clonality_summary_hmf_NSC$Mut_type == "All"])

<<<<<<< HEAD
for (i in 1:2){
  if (i == 1) {
    png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/comparison/figs/png/clonality-plot-final.png")
    print(clonality_plot_job)
    dev.off()
  }
  if (i == 2) {
    pdf(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/comparison/figs/pdf/clonality-plot-final.pdf")
    print(clonality_plot_job)
    dev.off()
  }
}
=======
############################################################################################################################################
>>>>>>> d4b2ba7d869eab81f62853d96bcbc0ac881b0e01

# ploidy

<<<<<<< HEAD
######################################################################################################################################################
# Ploidy


summary(pcawg_ploidy_df$ploidy)
summary(hmf_ploidy_df_NSC$ploidy)

pcawg_ploidy_df$Dataset <- "PCAWG"
hmf_ploidy_df_NSC$Dataset <- "HMF"


com_ploid <- rbind(pcawg_ploidy_df[,c("ploidy", "Dataset")], hmf_ploidy_df_NSC[,c("ploidy", "Dataset")])



com_ploid %<>% group_by(Dataset) %>%
  mutate(median = round(median(ploidy, na.rm = T), digits =1)) %>%
  # mutate(mean = round(mean(value, na.rm = T), digits = 2)) %>%
  ungroup()


com_ploid$Dataset <- factor(com_ploid$Dataset, levels = c("PCAWG", "HMF"))


ploidy_plot_job <- ggplot(data = com_ploid, aes(x = Dataset, y = ploidy, color = Dataset, fill = Dataset, label = median)) +
  geomJob_HalfHalfBox +
  ggplot2::scale_color_manual(values = c('#e75480', 'darkblue'), guide = "none") +
  ggplot2::scale_fill_manual(values = c('#e75480', 'darkblue'), guide = "none") +
  # ggplot2::scale_y_continuous(limits = c(20, 90), breaks = seq(0, 90, by = 10)) +
  ggplot2::ggtitle("Ploidy Level of Mutations in Lung Cohorts \n \n") +
  ggplot2::labs(x = '\n\n Cohorts', y = 'Ploidy \n') +
  theme_Job +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") +
  theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.5))
# theme(axis.title.x = element_text(face = "italic", size = 12))
# theme(axis.title.y = element_text(face = "italic", size = 12))
# theme(plot.margin = unit(c(1,1,1,1), "cm"))
# geom_text(data=data.frame(), aes(x=names(meds_mnv), y=meds_mnv,
#                                  label= c(paste0("Mean = ", pcawg_mnv_mean), paste0("Mean = ", hmf_mnv_mean))), col = "red", size=4) +


for (i in 1:2){
  if (i == 1) {
    png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/comparison/figs/png/ploidy-plot-final.png")
    print(ploidy_plot_job)
    dev.off()
  }
  if (i == 2) {
    pdf(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/comparison/figs/pdf/ploidy-plot-final.pdf")
    print(ploidy_plot_job)
    dev.off()
  }
}
=======




>>>>>>> d4b2ba7d869eab81f62853d96bcbc0ac881b0e01
