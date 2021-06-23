

# libraries

library(ggplot2)
library(gridExtra)
library(ggpubr)
library(grid)
library(magrittr)
library(dplyr)
library(tidyr)
library(gghighlight)
library(rstatix)

if (dir.exists("/hpc/cuppen/")){
  
  base_dir <- "/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/"
  devtools::load_all(paste0(base_dir,'/CHORD/processed/scripts_main/mutSigExtractor/'))
  
} else {
  
  library(mutSigExtractor)
}


# HMF data loading



if (dir.exists("/hpc/cuppen/")){
  hmf_meta <- read.csv(file = "/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/datasets/hmf/metadata_whitelisted.tsv", header = T,
                       sep = "\t", stringsAsFactors = F)
} else {
  hmf_meta <- read.csv(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/datasets/hmf/metadata_whitelisted.tsv", header = T,
                       sep = "\t", stringsAsFactors = F)
}


if (dir.exists("/hpc/cuppen/")){
  hmf_lung_meta <- read.csv(file = "/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/datasets/hmf/hmf-lung-cohort.tsv", header = T,
                       sep = "\t", stringsAsFactors = F)
} else {
  hmf_lung_meta <- read.csv(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/datasets/hmf/hmf-lung-cohort.tsv", header = T,
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

# PCAWG data loading


if (dir.exists("/hpc/cuppen/")){
  pcawg_meta <- readRDS(file = "/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/datasets/pcawg/metadata-per-donor.Rds")
} else {
  pcawg_meta <- readRDS(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/datasets/pcawg/metadata-per-donor.Rds")
}

# head(pcawg_lung_meta)

pcawg_lung_meta <- pcawg_meta[pcawg_meta$organ_system == "LUNG & BRONCHUS",]
pcawg_lung_meta <- pcawg_lung_meta[pcawg_lung_meta$icgc_donor_id != "DO23717",]


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


for (i in c(2:6, 9:17, 19, 23)) {
  pcawg_ploidy_df[,i] <- as.numeric(pcawg_ploidy_df[,i])
  
}


for (i in c(21, 25)) {
  pcawg_ploidy_df[,i] <- as.integer(pcawg_ploidy_df[,i])
  
}

############################################################################################################################################

# Boxplots for somatic SNVs, MNVs, and small insetions and deletions



options(scipen=999) 


### Using R base boxplot function


boxplot(count_df_pcawg$SNVs, count_df_hmf$SNVs, names = c("Primary Lung Cancer", "Metastatic Lung Cancer"))
hmf_snv_mean <- round(mean(count_df_hmf$SNVs), 2)
pcawg_snv_mean <- round(mean(count_df_pcawg$SNVs), 2)

text(x= 0.7, y= 200000, labels= paste0("Mean value of ", pcawg_snv_mean), font = 2, col = "red")
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
  geom_boxplot(outlier.colour="grey", outlier.shape=17,
               outlier.size=2, outlier.alpha = 0.5) +
  stat_compare_means(label.x = 1.3, label.y = 410000, method = "t.test") +
  geom_segment(aes(x = 1, y = 400000, xend = 2, yend = 400000), color = "black", size = 0.05) +
  geom_segment(aes(x = 1, y = 395000, xend = 1, yend = 400000), color = "black", size = 0.05) +
  geom_segment(aes(x = 2, y = 395000, xend = 2, yend = 400000), color = "black", size = 0.05) +
  geom_text(data=data.frame(), aes(x=names(meds_snv), y=meds_snv,
                                   label= c(paste0("Mean of ", pcawg_snv_mean), paste0("Mean of ", hmf_snv_mean))), col = c("red", "blue"), size=4, show.legend = FALSE) +
  # scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  scale_color_manual(name = "Dataset", values = c("red", "blue")) +
  geom_jitter(aes(color = label), size=0.4, alpha=0.5, show.legend = FALSE) +
  stat_boxplot(geom ='errorbar', show.legend = FALSE) +
  theme_bw() +
  ggtitle("Numbers of somatic SNVs \n") +
  ylab("Number of SNVs \n") +
  xlab("\n Cohorts") +
  theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.5)) +
  theme(axis.title.x = element_text(face = "italic", size = 12)) +
  theme(axis.title.y = element_text(face = "italic", size = 12)) +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) 


SNVs_plot


# ================================================================================================

### Using R basic graphics


boxplot(count_df_pcawg$MNVs, count_df_hmf$MNVs, names = c("Primary Lung Cancer", "Metastatic Lung Cancer"))
hmf_mnv_mean <- round(mean(count_df_hmf$MNVs), 2)
pcawg_mnv_mean <- round(mean(count_df_pcawg$MNVs), 2)

text(x= 0.7, y= 4000, labels= paste0("Mean value of ", pcawg_mnv_mean), font = 2, col = "red")
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
  geom_boxplot(outlier.colour="grey", outlier.shape=17,
               outlier.size=2, outlier.alpha = 0.5) +
  stat_compare_means(label.x = 1.3, label.y = 9800, method = "t.test") +
  geom_segment(aes(x = 1, y = 9600, xend = 2, yend = 9600), color = "black", size = 0.05) +
  geom_segment(aes(x = 1, y = 9400, xend = 1, yend = 9600), color = "black", size = 0.05) +
  geom_segment(aes(x = 2, y = 9400, xend = 2, yend = 9600), color = "black", size = 0.05) +
  
  geom_text(data=data.frame(), aes(x=names(meds_mnv), y=meds_mnv, 
                                   label=c(paste0("Mean of ", pcawg_mnv_mean), paste0("Mean of ", hmf_mnv_mean))), col = c("red", "blue"), size=4) +
  scale_color_manual(name = "Dataset", values = c("red", "blue")) +
  geom_jitter(aes(color = label), size=0.4, alpha=0.5) +
  stat_boxplot(geom ='errorbar') +
  theme_bw() +
  ggtitle("Numbers of somatic MNVs \n") +
  ylab("Number of MNVs \n") +
  xlab("\n Cohorts") +
  theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.5)) +
  theme(axis.title.x = element_text(face = "italic", size = 12)) +
  theme(axis.title.y = element_text(face = "italic", size = 12)) +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))

MNVs_plot


# ====================================================================================

### Using R basic graphics

boxplot(log10(count_df_pcawg$inDEL), log10(count_df_hmf$inDEL), names = c("Primary Lung Cancer", "Metastatic Lung Cancer"),
        main = "Number of small deletions \n \n \n", xlab = "Cohorts", ylab = "Log10 number of deletions", )
hmf_del_mean <- round(mean(count_df_hmf$inDEL), 2)
pcawg_del_mean <- round(mean(count_df_pcawg$inDEL), 2)

text(x= 1, y= 3.65, labels= paste0("Mean value of ", pcawg_del_mean), font = 2, col = "red")
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
  geom_boxplot(outlier.colour="grey", outlier.shape=17,
               outlier.size=2, outlier.alpha = 0.5) +
  stat_compare_means(label.x = 1.3, label.y = 18.5, method = "t.test") +
  geom_segment(aes(x = 1, y = 18.2, xend = 2, yend = 18.2), color = "black", size = 0.05) +
  geom_segment(aes(x = 1, y = 17.9, xend = 1, yend = 18.2), color = "black", size = 0.05) +
  geom_segment(aes(x = 2, y = 17.9, xend = 2, yend = 18.2), color = "black", size = 0.05) +
  # stat_compare_means(data = inDEL_plot_df, aes(x = label, y = value),comparisons = my_comparisons, method = "t.test", label.y = 2) +
  geom_text(data=data.frame(), aes(x=names(meds_del), y=meds_del, 
                                   label=c(paste0("Mean of ", pcawg_del_mean), paste0("Mean of ", hmf_del_mean))), col = c("red", "blue"), size=4) +
  scale_color_manual(name = "Dataset", values = c("red", "blue")) +
  geom_jitter(aes(color = label), size=0.4, alpha=0.5) +
  stat_boxplot(geom ='errorbar') +
  theme_bw() +
  ggtitle("Numbers of somatic small deletions \n") +
  ylim(c(0,20)) +
  ylab("Number of deletions (log2)\n") +
  xlab("\n Cohorts") +
  theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.5)) +
  theme(axis.title.x = element_text(face = "italic", size = 12)) +
  theme(axis.title.y = element_text(face = "italic", size = 12)) +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))

inDEL_plot



# =====================================================================================================================================


### Using R basic graphics

boxplot(log10(count_df_pcawg$INdel), log10(count_df_hmf$INdel), names = c("Primary Lung Cancer", "Metastatic Lung Cancer"), main = "Number of small insertions", xlab = "Cohorts", ylab = "Log10 number of insertions")
pcawg_ins_mean <- round(mean(count_df_pcawg$INdel), 2)
hmf_ins_mean <- round(mean(count_df_hmf$INdel), 2)

text(x= 0.7, y= 200000, labels= paste0("Mean value of ", pcawg_ins_mean), font = 2, col = "red")
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
  geom_boxplot(outlier.colour="grey", outlier.shape=17,
               outlier.size=2, outlier.alpha = 0.5) +
  stat_compare_means(label.x = 1.3, label.y = 16, method = "t.test") +
  geom_segment(aes(x = 1, y = 15.7, xend = 2, yend = 15.7), color = "black", size = 0.05) +
  geom_segment(aes(x = 1, y = 15.4, xend = 1, yend = 15.7), color = "black", size = 0.05) +
  geom_segment(aes(x = 2, y = 15.4, xend = 2, yend = 15.7), color = "black", size = 0.05) +
  geom_text(data=data.frame(), aes(x=names(meds_ins), y=meds_ins, 
                                   label=c(paste0("Mean of ", pcawg_ins_mean), paste0("Mean of ", hmf_ins_mean))), col = c("red", "blue"), size=4) +
  scale_color_manual(name = "Dataset", values = c("red", "blue")) +
  geom_jitter(aes(color = label), size=0.4, alpha=0.5) +
  stat_boxplot(geom ='errorbar') +
  theme_bw() +
  ggtitle("Numbers of somatic small insertions \n") +
  ylim(c(0,20)) +
  ylab("Number of insertions (log2) \n") +
  xlab("\n Cohorts") +
  theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.5)) +
  theme(axis.title.x = element_text(face = "italic", size = 12)) +
  theme(axis.title.y = element_text(face = "italic", size = 12)) +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))

INdel_plot



### combining all boxplots into one figure


for (i in 1:2){
  if (i == 1) {
    png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/comparison/figs/png/tmb-boxplots.png", width = 1000, height = 1000)
    print(ggarrange(SNVs_plot, MNVs_plot, inDEL_plot, INdel_plot, nrow = 2, ncol = 2, common.legend = TRUE, legend = "bottom"))
    dev.off()
  }
  if (i == 2) {
    pdf(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/comparison/figs/pdf/tmb-boxplots.pdf", width = 15, heigh = 15, onefile=FALSE)
    print(ggarrange(SNVs_plot, MNVs_plot, inDEL_plot, INdel_plot, nrow = 2, ncol = 2, common.legend = TRUE, legend = "bottom"))
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


sv_count_df_hmf <- as.data.frame(sv_count_matrix_hmf)
sv_count_df_hmf %<>% mutate(sampleId = rownames(sv_count_df_hmf))
rownames(sv_count_df_hmf) <- 1:nrow(sv_count_df_hmf)
sv_count_df_hmf <- sv_count_df_hmf[,c(6,1:5)]



sv_count_tibb_pcawg <- gather(sv_count_df_pcawg, key = "Simple_SV_Event", value = "Counts", 2:6)
sv_count_tibb_pcawg %<>% mutate(Dataset = "Primary")
sv_count_tibb_hmf <- gather(sv_count_df_hmf, key = "Simple_SV_Event", value = "Counts", 2:6)
sv_count_tibb_hmf %<>% mutate(Dataset = "Metastatic")


simple_sv_events <- c("Deletion (> 100kb)", "Deletion (< 100kb)", "Duplication (> 100kb)", "Duplication (< 100kb)", "Inversion")


sv_count_df_combined <- rbind(sv_count_tibb_pcawg, sv_count_tibb_hmf)

sv_count_df_combined$Simple_SV_Event <- factor(sv_count_df_combined$Simple_SV_Event, levels = simple_sv_events)
sv_count_df_combined$Dataset <- factor(sv_count_df_combined$Dataset, levels = c("Primary", "Metastatic"))


sv_count_df_combined <- sv_count_df_combined[sv_count_df_combined$Counts != 0,]


# t.test
head(sv_count_df_combined)

stat.test <- sv_count_df_combined %>%
  group_by(Simple_SV_Event) %>%
  t_test(Counts ~ Dataset) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")



stat.test <- stat.test %>%
  add_xy_position(x = "Simple_SV_Event", dodge = 0.8)

stat.test$y.position <- log2(stat.test$y.position)



sv_plot <- sv_count_df_combined %>% ggplot(aes(x = Simple_SV_Event, y = log2(Counts), color = Dataset)) +
  geom_boxplot(outlier.shape = NA) +
  # stat_compare_means(method = "t.test") +
  # geom_text(data=data.frame(), aes(x=names(meds_ins), y=meds_ins, 
  #                                  label=c(paste0("Mean of ", pcawg_ins_mean), paste0("Mean of ", hmf_ins_mean))), col = "red", size=4) +
  # scale_color_manual(values = c("blue", "green")) +
  # geom_jitter(color="black", size=0.4, alpha=0.5) +
  geom_point(position=position_jitterdodge(), size = 0.4, aes(color = Dataset), alpha = 0.5) +
  scale_color_manual(values = c("#ff0101", "#010dff")) +
  theme_classic() +
  ggtitle("SV Event Frequencies \n") +
  ylab("Frequency of Events (log2) \n") +
  ylim(0, 12) +
  xlab("Simple SV Events \n") +
  theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.5)) +
  theme(axis.title.x = element_text(face = "italic", size = 10)) +
  theme(axis.title.y = element_text(face = "italic", size = 10)) +
  theme(axis.text.x = element_text(angle = 45, size = 5, vjust = 0.6, face = "bold")) +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  stat_pvalue_manual(stat.test,  label = "p.adj.signif", tip.length = 0, hide.ns = F)


sv_plot


# Note that sv-boxplots1.pdf was obtained by looking into ".purple.cnv.somatic.tsv " files while for sv-boxplots2.pdf I used linx.vis_sv_data.tsv files


for (i in 1:2){
  if (i == 1) {
    png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/comparison/figs/png/sv-boxplots2.png")
    print(sv_plot)
    dev.off()
  }
  if (i == 2) {
    pdf(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/comparison/figs/pdf/sv-boxplots2.pdf")
    print(sv_plot)
    dev.off()
  }
}


# Always use Welch's t-test instead of Student's t-test (part. when the sample sizes are drastically unequal)

t.test(sv_count_df_combined$Counts[sv_count_df_combined$Simple_SV_Event == "Deletion (> 100kb)" & sv_count_df_combined$Dataset == "Primary"], sv_count_df_combined$Counts[sv_count_df_combined$Simple_SV_Event == "Deletion (> 100kb)" & sv_count_df_combined$Dataset == "Metastatic"], var.equal = F)
t.test(sv_count_df_combined$Counts[sv_count_df_combined$Simple_SV_Event == "Deletion (< 100kb)" & sv_count_df_combined$Dataset == "Primary"], sv_count_df_combined$Counts[sv_count_df_combined$Simple_SV_Event == "Deletion (< 100kb)" & sv_count_df_combined$Dataset == "Metastatic"], var.equal = F)
t.test(sv_count_df_combined$Counts[sv_count_df_combined$Simple_SV_Event == "Duplication (> 100kb)" & sv_count_df_combined$Dataset == "Primary"], sv_count_df_combined$Counts[sv_count_df_combined$Simple_SV_Event == "Duplication (> 100kb)" & sv_count_df_combined$Dataset == "Metastatic"])
t.test(sv_count_df_combined$Counts[sv_count_df_combined$Simple_SV_Event == "Duplication (< 100kb)" & sv_count_df_combined$Dataset == "Primary"], sv_count_df_combined$Counts[sv_count_df_combined$Simple_SV_Event == "Duplication (< 100kb)" & sv_count_df_combined$Dataset == "Metastatic"])
gg <- t.test(sv_count_df_combined$Counts[sv_count_df_combined$Simple_SV_Event == "Inversion" & sv_count_df_combined$Dataset == "Primary"], sv_count_df_combined$Counts[sv_count_df_combined$Simple_SV_Event == "Inversion" & sv_count_df_combined$Dataset == "Metastatic"])

class(gg)
names(gg)
gg$statistic
gg$p.value








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
hmf_lung_meta %<>% mutate(donor_age_at_biopsy = as.numeric(sapply(strsplit(hmf_lung_meta$biopsyDate, split = "-"), "[", 1)) - as.numeric(hmf_lung_meta$birthYear))

hmf_lung_meta_processed <- hmf_lung_meta[,c("sampleId", "gender", "donor_age_at_biopsy")]
hmf_lung_meta_processed %<>% mutate(Dataset = "HMF") %>% rename(donorSex = gender, donorAge = donor_age_at_biopsy)
all(hmf_lung_meta_processed$sampleId == count_df_hmf$sampleId)
hmf_lung_meta_processed_final <- cbind(hmf_lung_meta_processed, count_df_hmf[,2:6])
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
    png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/comparison/figs/png/age-tmb.png")
    print(my_plot)
    dev.off()
  }
  if (i == 2) {
    pdf(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/comparison/figs/pdf/age-tmb.pdf")
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




clonality_summary_hmf$Clona2Nonclonal <- clonality_summary_hmf$Clonal/clonality_summary_hmf$Non_clonal
clonality_summary_hmf_all <- clonality_summary_hmf[clonality_summary_hmf$Mut_type == "All",]
clonality_summary_hmf_all <- clonality_summary_hmf_all[clonality_summary_hmf_all$Clonal > 10 & clonality_summary_hmf_all$Non_clonal > 10,]

clonality_summary_hmf_all$Dataset <- "HMF"


combined_clonality <- rbind(clonality_summary_pcawg_all, clonality_summary_hmf_all)

combined_clonality$Dataset <- factor(combined_clonality$Dataset, levels = c("PCAWG", "HMF"))


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



for (i in 1:2){
  if (i == 1) {
    png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/comparison/figs/png/clonality_plot.png")
    print(clonality_plot)
    dev.off()
  }
  if (i == 2) {
    pdf(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/comparison/figs/pdf/clonality_plot.pdf")
    print(clonality_plot)
    dev.off()
  }
}


######################################################################################################################################################
# Ploidy


summary(pcawg_ploidy_df$ploidy)
summary(hmf_ploidy_df$ploidy)

boxplot(pcawg_ploidy_df$ploidy, hmf_ploidy_df$ploidy)



for (i in 1:2){
  if (i == 1) {png(filename = "/home/ali313/Desktop/das.png")
    boxplot(pcawg_ploidy_df$ploidy, hmf_ploidy_df$ploidy)
  }
  if (i == 2) {pdf(file = "/home/ali313/Desktop/das.pdf")
    boxplot(pcawg_ploidy_df$ploidy, hmf_ploidy_df$ploidy)
  }
  dev.off()
}


png(filename = "/home/ali313/Desktop/das.png")
pdf(file = "/home/ali313/Desktop/das.pdf")
jpeg()

dev.off()


t.test(pcawg_ploidy_df$ploidy, hmf_ploidy_df$ploidy)
wilcox.test(pcawg_ploidy_df$ploidy, hmf_ploidy_df$ploidy)

table(pcawg_ploidy_df$msStatus)
