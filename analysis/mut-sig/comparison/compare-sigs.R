
# Libraries
library(ggplot2)
library(MutationalPatterns)
library(stringr)
library(RColorBrewer)
library(ggdendro)
library(cowplot)
library(gridExtra)
library(ggpubr)
library(magrittr)
library(dplyr)
library(tidyr)
library(rstatix)


# Directories

wd <- "/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/mut-sig/"
local <- "/home/ali313/Documents/studies/master/umc-project"



# Loading mutSigExtractor package


if (dir.exists("/hpc/cuppen/")){
  
  base_dir <- "/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/"
  devtools::load_all(paste0(base_dir,'/CHORD/processed/scripts_main/mutSigExtractor/'))
  
} else {
  
  library(mutSigExtractor)
}


source("/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/themes.R")


# Reading in pcawg metadata


if (dir.exists("/hpc/cuppen/")){
  pcawg_meta <- readRDS(file = "/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/datasets/pcawg/metadata-per-donor.Rds")
} else {
  pcawg_meta <- readRDS(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/datasets/pcawg/metadata-per-donor.Rds")
}



pcawg_lung_meta <- pcawg_meta[pcawg_meta$organ_system == "LUNG & BRONCHUS",]
pcawg_lung_meta <- pcawg_lung_meta[pcawg_lung_meta$icgc_donor_id != "DO23717",]


# Reading in pcawg mut catolog mat and contrib.




if (dir.exists("/hpc/cuppen/")){
  sig_cont_pcawg <- readRDS(file = "/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/mut-sig/pcawg/list-sig-cont-lung-pcawg.rds")
  sbs_mut_context_mat_pcawg <- readRDS(file = "/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/mut-sig/pcawg/sbs_matrix-sig-context-pcawg.rds")
  dbs_mut_context_mat_pcawg <- readRDS(file = "/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/mut-sig/pcawg/dbs_matrix-sig-context-pcawg.rds")
  id_mut_context_mat_pcawg <- readRDS(file = "/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/mut-sig/pcawg/id_matrix-sig-context-pcawg.rds")
  
} else {
  sig_cont_pcawg <- readRDS(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/mut-sig/pcawg/list-sig-cont-lung-pcawg.rds")
  sbs_mut_context_mat_pcawg <- readRDS(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/mut-sig/pcawg/sbs_matrix-sig-context-pcawg.rds")
  dbs_mut_context_mat_pcawg <- readRDS(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/mut-sig/pcawg/dbs_matrix-sig-context-pcawg.rds")
  id_mut_context_mat_pcawg <- readRDS(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/mut-sig/pcawg/id_matrix-sig-context-pcawg.rds")
  
}


sbs_sig_cont_pcawg_df <- as.data.frame(matrix(unlist(lapply(sig_cont_pcawg, function(x) x[["SBS"]])), ncol = length(lapply(sig_cont_pcawg, function(x) x[["SBS"]])[[1]]), byrow = T))
colnames(sbs_sig_cont_pcawg_df) <- names(lapply(sig_cont_pcawg, function(x) x[["SBS"]])[[1]])
rownames(sbs_sig_cont_pcawg_df) <- names(sig_cont_pcawg)

dbs_sig_cont_pcawg_df <- as.data.frame(matrix(unlist(lapply(sig_cont_pcawg, function(x) x[["DBS"]])), ncol = length(lapply(sig_cont_pcawg, function(x) x[["DBS"]])[[1]]), byrow = T))
colnames(dbs_sig_cont_pcawg_df) <- names(lapply(sig_cont_pcawg, function(x) x[["DBS"]])[[1]])
rownames(dbs_sig_cont_pcawg_df) <- names(sig_cont_pcawg)

id_sig_cont_pcawg_df <- as.data.frame(matrix(unlist(lapply(sig_cont_pcawg, function(x) x[["ID"]])), ncol = length(lapply(sig_cont_pcawg, function(x) x[["ID"]])[[1]]), byrow = T))
colnames(id_sig_cont_pcawg_df) <- names(lapply(sig_cont_pcawg, function(x) x[["ID"]])[[1]])
rownames(id_sig_cont_pcawg_df) <- names(sig_cont_pcawg)




# Reading in hmf metadata


if (dir.exists("/hpc/cuppen/")){
  hmf_meta <- read.csv(file = "/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/datasets/hmf/metadata_whitelisted.tsv", header = T,
                       sep = "\t", stringsAsFactors = F)
} else {
  hmf_meta <- read.csv(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/datasets/hmf/metadata_whitelisted.tsv", header = T,
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


# Reading in hmf mut catolog mat and contrib.


if (dir.exists("/hpc/cuppen/")){
  sig_cont_hmf <- readRDS(file = "/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/mut-sig/hmf/list-sig-cont-lung-hmf.rds")
  sbs_mut_context_mat_hmf <- readRDS(file = paste0(wd, "sbs_matrix-sig-context-hmf.rds"))
  dbs_mut_context_mat_hmf <- readRDS(file = paste0(wd, "dbs_matrix-sig-context-hmf.rds"))
  id_mut_context_mat_hmf <- readRDS(file = paste0(wd, "id_matrix-sig-context-hmf.rds"))
  
} else {
  sig_cont_hmf <- readRDS(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/mut-sig/hmf/list-sig-cont-lung-hmf.rds")
  sbs_mut_context_mat_hmf <- readRDS(file = paste0(local, wd, "hmf/sbs_matrix-sig-context-hmf.rds"))
  dbs_mut_context_mat_hmf <- readRDS(file = paste0(local, wd, "hmf/dbs_matrix-sig-context-hmf.rds"))
  id_mut_context_mat_hmf <- readRDS(file = paste0(local, wd, "hmf/id_matrix-sig-context-hmf.rds"))
  
}

sbs_sig_cont_hmf_df <- as.data.frame(matrix(unlist(lapply(sig_cont_hmf, function(x) x[["SBS"]])), ncol = length(lapply(sig_cont_hmf, function(x) x[["SBS"]])[[1]]), byrow = T))
colnames(sbs_sig_cont_hmf_df) <- names(lapply(sig_cont_hmf, function(x) x[["SBS"]])[[1]])
rownames(sbs_sig_cont_hmf_df) <- names(sig_cont_hmf)

dbs_sig_cont_hmf_df <- as.data.frame(matrix(unlist(lapply(sig_cont_hmf, function(x) x[["DBS"]])), ncol = length(lapply(sig_cont_hmf, function(x) x[["DBS"]])[[1]]), byrow = T))
colnames(dbs_sig_cont_hmf_df) <- names(lapply(sig_cont_hmf, function(x) x[["DBS"]])[[1]])
rownames(dbs_sig_cont_hmf_df) <- names(sig_cont_hmf)

id_sig_cont_hmf_df <- as.data.frame(matrix(unlist(lapply(sig_cont_hmf, function(x) x[["ID"]])), ncol = length(lapply(sig_cont_hmf, function(x) x[["ID"]])[[1]]), byrow = T))
colnames(id_sig_cont_hmf_df) <- names(lapply(sig_cont_hmf, function(x) x[["ID"]])[[1]])
rownames(id_sig_cont_hmf_df) <- names(sig_cont_hmf)


sbs_sig_cont_hmf_df_NSC <- sbs_sig_cont_hmf_df[rownames(sbs_sig_cont_hmf_df) %in% hmf_lung_meta_NSC$sampleId,]
dbs_sig_cont_hmf_df_NSC <- dbs_sig_cont_hmf_df[rownames(dbs_sig_cont_hmf_df) %in% hmf_lung_meta_NSC$sampleId,]
id_sig_cont_hmf_df_NSC <- id_sig_cont_hmf_df[rownames(id_sig_cont_hmf_df) %in% hmf_lung_meta_NSC$sampleId,]



#################################################################################################################################################
# Comparing SBS sigs

## PCAWG

### Absolute
copy_sbs_sig_cont_pcawg_df <- sbs_sig_cont_pcawg_df


copy_sbs_sig_cont_pcawg_df %<>% mutate(sampleId = rownames(copy_sbs_sig_cont_pcawg_df)) %>% relocate(where(is.character), .before = where(is.numeric))
copy_sbs_sig_cont_pcawg_df$sum.of.platinum.sigs.denovo_2_6 <- rowSums(copy_sbs_sig_cont_pcawg_df[,c("SBS35.Platinum.72.denovo_2", "SBS31.Platinum.98.denovo_6")])
copy_sbs_sig_cont_pcawg_df$sum.of.age.sigs.denovo_4_7_10 <- rowSums(copy_sbs_sig_cont_pcawg_df[,c("SBS5.Age.77.denovo_4", "SBS1.Age.97.denovo_7", "SBS40.Age.80.denovo_10")])
copy_sbs_sig_cont_pcawg_df$sum.of.mmrd.sigs.denovo_5_8 <- rowSums(copy_sbs_sig_cont_pcawg_df[,c("SBS44.MMRd.89.denovo_5", "SBS26.MMRd.45.denovo_8")])
copy_sbs_sig_cont_pcawg_df$sum.of.tobacco.sigs.denovo_9_12 <- rowSums(copy_sbs_sig_cont_pcawg_df[,c("SBS4.Tobacco.90.denovo_9", "SBS29.Tobacco_chewing.69.denovo_12")])
copy_sbs_sig_cont_pcawg_df$sum.of.apobec.sigs.denovo_11_13 <- rowSums(copy_sbs_sig_cont_pcawg_df[,c("SBS2.APOBEC.99.denovo_11", "SBS13.APOBEC.83.denovo_13")])
copy_sbs_sig_cont_pcawg_df$sum.of.ros.sigs.denovo_14_15 <- rowSums(copy_sbs_sig_cont_pcawg_df[,c("SBS18.ROS.89.denovo_14", "SBS17b.ROS_5FU.96.denovo_15")])

rownames(copy_sbs_sig_cont_pcawg_df) <- 1:nrow(copy_sbs_sig_cont_pcawg_df)
copy_sbs_sig_cont_pcawg_tibb <- gather(copy_sbs_sig_cont_pcawg_df, key = "SBS_sig", value = "Abs_count", 2:22)
copy_sbs_sig_cont_pcawg_tibb %<>% mutate(Dataset = "Primary")





### Relative

copy2_sbs_sig_cont_pcawg_df <- sbs_sig_cont_pcawg_df


copy2_sbs_sig_cont_pcawg_df$sum.of.platinum.sigs.denovo_2_6 <- rowSums(copy2_sbs_sig_cont_pcawg_df[,c("SBS35.Platinum.72.denovo_2", "SBS31.Platinum.98.denovo_6")])
copy2_sbs_sig_cont_pcawg_df$sum.of.age.sigs.denovo_4_7_10 <- rowSums(copy2_sbs_sig_cont_pcawg_df[,c("SBS5.Age.77.denovo_4", "SBS1.Age.97.denovo_7", "SBS40.Age.80.denovo_10")])
copy2_sbs_sig_cont_pcawg_df$sum.of.mmrd.sigs.denovo_5_8 <- rowSums(copy2_sbs_sig_cont_pcawg_df[,c("SBS44.MMRd.89.denovo_5", "SBS26.MMRd.45.denovo_8")])
copy2_sbs_sig_cont_pcawg_df$sum.of.tobacco.sigs.denovo_9_12 <- rowSums(copy2_sbs_sig_cont_pcawg_df[,c("SBS4.Tobacco.90.denovo_9", "SBS29.Tobacco_chewing.69.denovo_12")])
copy2_sbs_sig_cont_pcawg_df$sum.of.apobec.sigs.denovo_11_13 <- rowSums(copy2_sbs_sig_cont_pcawg_df[,c("SBS2.APOBEC.99.denovo_11", "SBS13.APOBEC.83.denovo_13")])
copy2_sbs_sig_cont_pcawg_df$sum.of.ros.sigs.denovo_14_15 <- rowSums(copy2_sbs_sig_cont_pcawg_df[,c("SBS18.ROS.89.denovo_14", "SBS17b.ROS_5FU.96.denovo_15")])

norm_sbs_sig_cont_pcawg_df <- copy2_sbs_sig_cont_pcawg_df/rowSums(copy2_sbs_sig_cont_pcawg_df)


norm_sbs_sig_cont_pcawg_df %<>% mutate(sampleId = rownames(norm_sbs_sig_cont_pcawg_df)) %>% relocate(where(is.character), .before = where(is.numeric))
rownames(norm_sbs_sig_cont_pcawg_df) <- 1:nrow(norm_sbs_sig_cont_pcawg_df)
norm_sbs_sig_cont_pcawg_tibb <- gather(norm_sbs_sig_cont_pcawg_df, key = "SBS_sig", value = "Rel_count", 2:22)
norm_sbs_sig_cont_pcawg_tibb %<>% mutate(Dataset = "Primary")



## HMF

### Absolute


copy_sbs_sig_cont_hmf_df <- sbs_sig_cont_hmf_df_NSC

copy_sbs_sig_cont_hmf_df %<>% mutate(sampleId = rownames(copy_sbs_sig_cont_hmf_df)) %>% relocate(where(is.character), .before = where(is.numeric))

copy_sbs_sig_cont_hmf_df$sum.of.platinum.sigs.denovo_2_6 <- rowSums(copy_sbs_sig_cont_hmf_df[,c("SBS35.Platinum.72.denovo_2", "SBS31.Platinum.98.denovo_6")])
copy_sbs_sig_cont_hmf_df$sum.of.age.sigs.denovo_4_7_10 <- rowSums(copy_sbs_sig_cont_hmf_df[,c("SBS5.Age.77.denovo_4", "SBS1.Age.97.denovo_7", "SBS40.Age.80.denovo_10")])
copy_sbs_sig_cont_hmf_df$sum.of.mmrd.sigs.denovo_5_8 <- rowSums(copy_sbs_sig_cont_hmf_df[,c("SBS44.MMRd.89.denovo_5", "SBS26.MMRd.45.denovo_8")])
copy_sbs_sig_cont_hmf_df$sum.of.tobacco.sigs.denovo_9_12 <- rowSums(copy_sbs_sig_cont_hmf_df[,c("SBS4.Tobacco.90.denovo_9", "SBS29.Tobacco_chewing.69.denovo_12")])
copy_sbs_sig_cont_hmf_df$sum.of.apobec.sigs.denovo_11_13 <- rowSums(copy_sbs_sig_cont_hmf_df[,c("SBS2.APOBEC.99.denovo_11", "SBS13.APOBEC.83.denovo_13")])
copy_sbs_sig_cont_hmf_df$sum.of.ros.sigs.denovo_14_15 <- rowSums(copy_sbs_sig_cont_hmf_df[,c("SBS18.ROS.89.denovo_14", "SBS17b.ROS_5FU.96.denovo_15")])

rownames(copy_sbs_sig_cont_hmf_df) <- 1:nrow(copy_sbs_sig_cont_hmf_df)
copy_sbs_sig_cont_hmf_tibb <- gather(copy_sbs_sig_cont_hmf_df, key = "SBS_sig", value = "Abs_count", 2:22)
copy_sbs_sig_cont_hmf_tibb %<>% mutate(Dataset = "Metastatic")



### Relative
copy2_sbs_sig_cont_hmf_df <- sbs_sig_cont_hmf_df_NSC



copy2_sbs_sig_cont_hmf_df$sum.of.platinum.sigs.denovo_2_6 <- rowSums(copy2_sbs_sig_cont_hmf_df[,c("SBS35.Platinum.72.denovo_2", "SBS31.Platinum.98.denovo_6")])
copy2_sbs_sig_cont_hmf_df$sum.of.age.sigs.denovo_4_7_10 <- rowSums(copy2_sbs_sig_cont_hmf_df[,c("SBS5.Age.77.denovo_4", "SBS1.Age.97.denovo_7", "SBS40.Age.80.denovo_10")])
copy2_sbs_sig_cont_hmf_df$sum.of.mmrd.sigs.denovo_5_8 <- rowSums(copy2_sbs_sig_cont_hmf_df[,c("SBS44.MMRd.89.denovo_5", "SBS26.MMRd.45.denovo_8")])
copy2_sbs_sig_cont_hmf_df$sum.of.tobacco.sigs.denovo_9_12 <- rowSums(copy2_sbs_sig_cont_hmf_df[,c("SBS4.Tobacco.90.denovo_9", "SBS29.Tobacco_chewing.69.denovo_12")])
copy2_sbs_sig_cont_hmf_df$sum.of.apobec.sigs.denovo_11_13 <- rowSums(copy2_sbs_sig_cont_hmf_df[,c("SBS2.APOBEC.99.denovo_11", "SBS13.APOBEC.83.denovo_13")])
copy2_sbs_sig_cont_hmf_df$sum.of.ros.sigs.denovo_14_15 <- rowSums(copy2_sbs_sig_cont_hmf_df[,c("SBS18.ROS.89.denovo_14", "SBS17b.ROS_5FU.96.denovo_15")])

norm_sbs_sig_cont_hmf_df <- copy2_sbs_sig_cont_hmf_df/rowSums(copy2_sbs_sig_cont_hmf_df)



norm_sbs_sig_cont_hmf_df %<>% mutate(sampleId = rownames(norm_sbs_sig_cont_hmf_df)) %>% relocate(where(is.character), .before = where(is.numeric))
rownames(norm_sbs_sig_cont_hmf_df) <- 1:nrow(norm_sbs_sig_cont_hmf_df)
norm_sbs_sig_cont_hmf_tibb <- gather(norm_sbs_sig_cont_hmf_df, key = "SBS_sig", value = "Rel_count", 2:22)
norm_sbs_sig_cont_hmf_tibb %<>% mutate(Dataset = "Metastatic")





## Combine 


### Absolute



sbs_sig_cont_comb_absolute <- rbind(copy_sbs_sig_cont_pcawg_tibb, copy_sbs_sig_cont_hmf_tibb)
table(sbs_sig_cont_comb_absolute$SBS_sig)
sbs_sig_cont_comb_absolute$SBS_sig <- factor(sbs_sig_cont_comb_absolute$SBS_sig, levels = colnames(copy2_sbs_sig_cont_pcawg_df))
sbs_sig_cont_comb_absolute$Dataset <- factor(sbs_sig_cont_comb_absolute$Dataset, levels = c("Primary", "Metastatic"))




# sum(sbs_sig_cont_comb_absolute$Rel_cont == 0)

sbs_sig_cont_comb_absolute_sep <- sbs_sig_cont_comb_absolute[sbs_sig_cont_comb_absolute$SBS_sig %in% colnames(sbs_sig_cont_pcawg_df),]
sbs_sig_cont_comb_absolute_sep$SBS_sig <- factor(sbs_sig_cont_comb_absolute_sep$SBS_sig, levels = unique(sbs_sig_cont_comb_absolute_sep$SBS_sig))

stat.test <- sbs_sig_cont_comb_absolute_sep %>%
  group_by(SBS_sig) %>%
  wilcox_test(Abs_count ~ Dataset) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")



stat.test <- stat.test %>%
  add_xy_position(x = "SBS_sig", dodge = 0.8)
stat.test$y.position <- stat.test$y.position



sbs_absolute_sep_contrib_plot <-  sbs_sig_cont_comb_absolute_sep %>% ggplot(aes(x = SBS_sig, y = Abs_count, color = Dataset)) +
  geom_boxplot(outlier.shape = NA) +
  # stat_compare_means(method = "t.test") +
  # geom_text(data=data.frame(), aes(x=names(meds_ins), y=meds_ins, 
  #                                  label=c(paste0("Mean of ", pcawg_ins_mean), paste0("Mean of ", hmf_ins_mean))), col = "red", size=4) +
  # scale_color_manual(values = c("blue", "green")) +
  # geom_jitter(color="black", size=0.4, alpha=0.5) +
  geom_point(position=position_jitterdodge(), size = 0.4, aes(color = Dataset), alpha = 0.5) +
  scale_color_manual(values = c("#ff0101", "#010dff")) +
  theme_classic() +
  ggtitle("Single-base Substitution Contribution Comparison \n \n") +
  # ylab("Absolute Contribution (%) \n") +
  # ylim(0, 110) +
  xlab("\n Lung-specific SBS Signatures") +
  theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.4)) +
  theme(axis.title.x = element_text(face = "italic", size = 10)) +
  theme(axis.title.y = element_text(face = "italic", size = 10)) +
  theme(axis.text.x = element_text(angle = 90, size = 5, vjust = 0.6, face = "bold")) +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  scale_y_continuous(name = "Absolute Contribution \n") +
  stat_pvalue_manual(stat.test,  label = "p.adj.signif", tip.length = 0, hide.ns = F)


sbs_absolute_sep_contrib_plot




sbs_sig_cont_comb_absolute_com <- sbs_sig_cont_comb_absolute[sbs_sig_cont_comb_absolute$SBS_sig %notin% colnames(sbs_sig_cont_pcawg_df),]
sbs_sig_cont_comb_absolute_com$SBS_sig <- factor(sbs_sig_cont_comb_absolute_com$SBS_sig, levels = unique(sbs_sig_cont_comb_absolute_com$SBS_sig))


stat.test <- sbs_sig_cont_comb_absolute_com %>%
  group_by(SBS_sig) %>%
  wilcox_test(Abs_count ~ Dataset) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")



stat.test <- stat.test %>%
  add_xy_position(x = "SBS_sig", dodge = 0.8)
stat.test$y.position <- stat.test$y.position



sbs_absolute_com_contrib_plot <-  sbs_sig_cont_comb_absolute_com %>% ggplot(aes(x = SBS_sig, y = Abs_count, color = Dataset)) +
  geom_boxplot(outlier.shape = NA) +
  # stat_compare_means(method = "t.test") +
  # geom_text(data=data.frame(), aes(x=names(meds_ins), y=meds_ins, 
  #                                  label=c(paste0("Mean of ", pcawg_ins_mean), paste0("Mean of ", hmf_ins_mean))), col = "red", size=4) +
  # scale_color_manual(values = c("blue", "green")) +
  # geom_jitter(color="black", size=0.4, alpha=0.5) +
  geom_point(position=position_jitterdodge(), size = 0.4, aes(color = Dataset), alpha = 0.5) +
  scale_color_manual(values = c("#ff0101", "#010dff")) +
  theme_classic() +
  ggtitle("Single-base Substitution Contribution Comparison \n \n") +
  # ylab("Absolute Contribution (%) \n") +
  # ylim(0, 110) +
  xlab("\n Lung-specific SBS Signatures") +
  theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.4)) +
  theme(axis.title.x = element_text(face = "italic", size = 10)) +
  theme(axis.title.y = element_text(face = "italic", size = 10)) +
  theme(axis.text.x = element_text(angle = 90, size = 5, vjust = 0.6, face = "bold")) +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  scale_y_continuous(name = "Absolute Contribution \n") +
  stat_pvalue_manual(stat.test,  label = "p.adj.signif", tip.length = 0, hide.ns = F)


sbs_absolute_com_contrib_plot










for (i in 1:2){
  if (i == 1) {
    png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/mut-sig/comparison/figs/png/sbs-comparison-absolute-contribution-sep-final.png")
    print(sbs_absolute_sep_contrib_plot)
    dev.off()
    
    png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/mut-sig/comparison/figs/png/sbs-comparison-absolute-contribution-comb-final.png")
    print(sbs_absolute_com_contrib_plot)
    dev.off()
    
  }
  if (i == 2) {
    pdf(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/mut-sig/comparison/figs/pdf/sbs-comparison-absolute-contribution-sep-final.pdf")
    print(sbs_absolute_sep_contrib_plot)
    dev.off()
    
    pdf(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/mut-sig/comparison/figs/pdf/sbs-comparison-absolute-contribution-comb-final.pdf")
    print(sbs_absolute_com_contrib_plot)
    dev.off()
    
  }
}






### Relative

sbs_sig_cont_comb_relative <- rbind(norm_sbs_sig_cont_pcawg_tibb, norm_sbs_sig_cont_hmf_tibb)
sbs_sig_cont_comb_relative$SBS_sig <- factor(sbs_sig_cont_comb_relative$SBS_sig, levels = colnames(copy2_sbs_sig_cont_pcawg_df))
sbs_sig_cont_comb_relative$Dataset <- factor(sbs_sig_cont_comb_relative$Dataset, levels = c("Primary", "Metastatic"))


# sum(sbs_sig_cont_comb_relative$Rel_cont == 0)


sbs_sig_cont_comb_relative_sep <- sbs_sig_cont_comb_relative[sbs_sig_cont_comb_relative$SBS_sig %in% colnames(sbs_sig_cont_pcawg_df),]
sbs_sig_cont_comb_relative_sep$SBS_sig <- factor(sbs_sig_cont_comb_relative_sep$SBS_sig, levels = unique(sbs_sig_cont_comb_relative_sep$SBS_sig))




stat.test <- sbs_sig_cont_comb_relative_sep %>%
  group_by(SBS_sig) %>%
  wilcox_test(Rel_count ~ Dataset) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")



stat.test <- stat.test %>%
  add_xy_position(x = "SBS_sig", dodge = 0.8)
stat.test$y.position <- stat.test$y.position*100


sbs_relative_sep_contrib_plot <- sbs_sig_cont_comb_relative_sep %>% ggplot(aes(x = SBS_sig, y = Rel_count*100, color = Dataset)) +
  geom_boxplot(outlier.shape = NA) +
  # stat_compare_means(method = "t.test") +
  # geom_text(data=data.frame(), aes(x=names(meds_ins), y=meds_ins, 
  #                                  label=c(paste0("Mean of ", pcawg_ins_mean), paste0("Mean of ", hmf_ins_mean))), col = "red", size=4) +
  # scale_color_manual(values = c("blue", "green")) +
  # geom_jitter(color="black", size=0.4, alpha=0.5) +
  geom_point(position=position_jitterdodge(), size = 0.4, aes(color = Dataset), alpha = 0.5) +
  scale_color_manual(values = c("#ff0101", "#010dff")) +
  theme_classic() +
  ggtitle("Single-base Substitution Contribution Comparison \n \n") +
  # ylab("Relative Contribution (%) \n") +
  # ylim(0, 110) +
  xlab("\n Lung-specific SBS Signatures") +
  theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.4)) +
  theme(axis.title.x = element_text(face = "italic", size = 10)) +
  theme(axis.title.y = element_text(face = "italic", size = 10)) +
  theme(axis.text.x = element_text(angle = 90, size = 5, vjust = 0.6, face = "bold")) +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  scale_y_continuous(name = "Relative Contribution (%) \n", breaks = c(0,20,40,60,80,100), labels = c(0,20,40,60,80,100)) +
  stat_pvalue_manual(stat.test,  label = "p.adj.signif", tip.length = 0, hide.ns = F)


sbs_relative_sep_contrib_plot





sbs_sig_cont_comb_relative_com <- sbs_sig_cont_comb_relative[sbs_sig_cont_comb_relative$SBS_sig %notin% colnames(sbs_sig_cont_pcawg_df),]
sbs_sig_cont_comb_relative_com$SBS_sig <- factor(sbs_sig_cont_comb_relative_com$SBS_sig, levels = unique(sbs_sig_cont_comb_relative_com$SBS_sig))




stat.test <- sbs_sig_cont_comb_relative_com %>%
  group_by(SBS_sig) %>%
  wilcox_test(Rel_count ~ Dataset) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")



stat.test <- stat.test %>%
  add_xy_position(x = "SBS_sig", dodge = 0.8)
stat.test$y.position <- stat.test$y.position*100


sbs_relative_com_contrib_plot <- sbs_sig_cont_comb_relative_com %>% ggplot(aes(x = SBS_sig, y = Rel_count*100, color = Dataset)) +
  geom_boxplot(outlier.shape = NA) +
  # stat_compare_means(method = "t.test") +
  # geom_text(data=data.frame(), aes(x=names(meds_ins), y=meds_ins, 
  #                                  label=c(paste0("Mean of ", pcawg_ins_mean), paste0("Mean of ", hmf_ins_mean))), col = "red", size=4) +
  # scale_color_manual(values = c("blue", "green")) +
  # geom_jitter(color="black", size=0.4, alpha=0.5) +
  geom_point(position=position_jitterdodge(), size = 0.4, aes(color = Dataset), alpha = 0.5) +
  scale_color_manual(values = c("#ff0101", "#010dff")) +
  theme_classic() +
  ggtitle("Single-base Substitution Contribution Comparison \n \n") +
  # ylab("Relative Contribution (%) \n") +
  # ylim(0, 110) +
  xlab("\n Lung-specific SBS Signatures") +
  theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.4)) +
  theme(axis.title.x = element_text(face = "italic", size = 10)) +
  theme(axis.title.y = element_text(face = "italic", size = 10)) +
  theme(axis.text.x = element_text(angle = 90, size = 5, vjust = 0.6, face = "bold")) +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  scale_y_continuous(name = "Relative Contribution (%) \n", breaks = c(0,20,40,60,80,100), labels = c(0,20,40,60,80,100)) +
  stat_pvalue_manual(stat.test,  label = "p.adj.signif", tip.length = 0, hide.ns = F)


sbs_relative_com_contrib_plot



for (i in 1:2){
  if (i == 1) {
    png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/mut-sig/comparison/figs/png/sbs-comparison-relative-contribution-sep-final.png")
    print(sbs_relative_sep_contrib_plot)
    dev.off()
    
    png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/mut-sig/comparison/figs/png/sbs-comparison-relative-contribution-comb-final.png")
    print(sbs_relative_com_contrib_plot)
    dev.off()
    
  }
  if (i == 2) {
    pdf(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/mut-sig/comparison/figs/pdf/sbs-comparison-relative-contribution-sep-final.pdf")
    print(sbs_relative_sep_contrib_plot)
    dev.off()
    
    pdf(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/mut-sig/comparison/figs/pdf/sbs-comparison-relative-contribution-comb-final.pdf")
    print(sbs_relative_com_contrib_plot)
    dev.off()
    
  }
}






# Making sense of the observation that there are platin-induced mutagenesis but not 5-fu
sum(str_detect(hmf_lung_meta$treatment, pattern = "platin"), na.rm = T)
sum(str_detect(hmf_meta$treatment, pattern = "platin"), na.rm = T)
sum(str_detect(hmf_lung_meta$treatment, pattern = "Fluorouracil"), na.rm = T)
sum(str_detect(hmf_meta$treatment, pattern = "Fluorouracil"), na.rm = T)
hmf_meta$treatment[str_detect(hmf_meta$treatment, pattern = "Fluorouracil")]


#################################################################################################################################################
# Comparing DBS sigs





## PCAWG

### Absolute
copy_dbs_sig_cont_pcawg_df <- dbs_sig_cont_pcawg_df

copy_dbs_sig_cont_pcawg_df %<>% mutate(sampleId = rownames(copy_dbs_sig_cont_pcawg_df)) %>% relocate(where(is.character), .before = where(is.numeric))
rownames(copy_dbs_sig_cont_pcawg_df) <- 1:nrow(copy_dbs_sig_cont_pcawg_df)
copy_dbs_sig_cont_pcawg_tibb <- gather(copy_dbs_sig_cont_pcawg_df, key = "dbs_sig", value = "Abs_count", 2:4)
copy_dbs_sig_cont_pcawg_tibb %<>% mutate(Dataset = "Primary")





### Relative

norm_dbs_sig_cont_pcawg_df <- dbs_sig_cont_pcawg_df/rowSums(dbs_sig_cont_pcawg_df)


norm_dbs_sig_cont_pcawg_df %<>% mutate(sampleId = rownames(norm_dbs_sig_cont_pcawg_df)) %>% relocate(where(is.character), .before = where(is.numeric))
rownames(norm_dbs_sig_cont_pcawg_df) <- 1:nrow(norm_dbs_sig_cont_pcawg_df)
norm_dbs_sig_cont_pcawg_tibb <- gather(norm_dbs_sig_cont_pcawg_df, key = "dbs_sig", value = "Rel_count", 2:4)
norm_dbs_sig_cont_pcawg_tibb %<>% mutate(Dataset = "Primary")



## HMF

### Absolute

copy_dbs_sig_cont_hmf_df <- dbs_sig_cont_hmf_df_NSC

copy_dbs_sig_cont_hmf_df %<>% mutate(sampleId = rownames(copy_dbs_sig_cont_hmf_df)) %>% relocate(where(is.character), .before = where(is.numeric))
rownames(copy_dbs_sig_cont_hmf_df) <- 1:nrow(copy_dbs_sig_cont_hmf_df)
copy_dbs_sig_cont_hmf_tibb <- gather(copy_dbs_sig_cont_hmf_df, key = "dbs_sig", value = "Abs_count", 2:4)
copy_dbs_sig_cont_hmf_tibb %<>% mutate(Dataset = "Metastatic")



### Relative

norm_dbs_sig_cont_hmf_df <- dbs_sig_cont_hmf_df_NSC/rowSums(dbs_sig_cont_hmf_df_NSC)



norm_dbs_sig_cont_hmf_df %<>% mutate(sampleId = rownames(norm_dbs_sig_cont_hmf_df)) %>% relocate(where(is.character), .before = where(is.numeric))
rownames(norm_dbs_sig_cont_hmf_df) <- 1:nrow(norm_dbs_sig_cont_hmf_df)
norm_dbs_sig_cont_hmf_tibb <- gather(norm_dbs_sig_cont_hmf_df, key = "dbs_sig", value = "Rel_count", 2:4)
norm_dbs_sig_cont_hmf_tibb %<>% mutate(Dataset = "Metastatic")





## Combine 


### Absolute





dbs_sig_cont_comb_absolute <- rbind(copy_dbs_sig_cont_pcawg_tibb, copy_dbs_sig_cont_hmf_tibb)

dbs_sig_cont_comb_absolute$dbs_sig <- factor(dbs_sig_cont_comb_absolute$dbs_sig, levels = colnames(dbs_sig_cont_pcawg_df))
dbs_sig_cont_comb_absolute$Dataset <- factor(dbs_sig_cont_comb_absolute$Dataset, levels = c("Primary", "Metastatic"))


# sum(dbs_sig_cont_comb_absolute$Rel_cont == 0)


stat.test <- dbs_sig_cont_comb_absolute %>%
  group_by(dbs_sig) %>%
  wilcox_test(Abs_count ~ Dataset) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")



stat.test <- stat.test %>%
  add_xy_position(x = "dbs_sig", dodge = 0.8)
stat.test$y.position <- stat.test$y.position


dbs_absolute_contrib_plot <- dbs_sig_cont_comb_absolute %>% ggplot(aes(x = dbs_sig, y = Abs_count, color = Dataset)) +
  geom_boxplot(outlier.shape = NA) +
  # stat_compare_means(method = "t.test") +
  # geom_text(data=data.frame(), aes(x=names(meds_ins), y=meds_ins, 
  #                                  label=c(paste0("Mean of ", pcawg_ins_mean), paste0("Mean of ", hmf_ins_mean))), col = "red", size=4) +
  # scale_color_manual(values = c("blue", "green")) +
  # geom_jitter(color="black", size=0.4, alpha=0.5) +
  geom_point(position=position_jitterdodge(), size = 0.4, aes(color = Dataset), alpha = 0.5) +
  scale_color_manual(values = c("#ff0101", "#010dff")) +
  theme_classic() +
  ggtitle("Double-base Substitution Contribution Comparison \n \n") +
  # ylab("Absolute Contribution (%) \n") +
  # ylim(0, 110) +
  xlab("\n Lung-specific DBS Signatures") +
  theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.4)) +
  theme(axis.title.x = element_text(face = "italic", size = 10)) +
  theme(axis.title.y = element_text(face = "italic", size = 10)) +
  theme(axis.text.x = element_text(angle = 90, size = 5, vjust = 0.6, face = "bold")) +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  scale_y_continuous(name = "Absolute Contribution \n") +
  stat_pvalue_manual(stat.test,  label = "p.adj.signif", tip.length = 0, hide.ns = F)


dbs_absolute_contrib_plot


for (i in 1:2){
  if (i == 1) {
    png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/mut-sig/comparison/figs/png/dbs-comparison-absolute-contribution-final.png")
    print(dbs_absolute_contrib_plot)
    dev.off()
  }
  if (i == 2) {
    pdf(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/mut-sig/comparison/figs/pdf/dbs-comparison-absolute-contribution-final.pdf")
    print(dbs_absolute_contrib_plot)
    dev.off()
  }
}








### Relative

dbs_sig_cont_comb_relative <- rbind(norm_dbs_sig_cont_pcawg_tibb, norm_dbs_sig_cont_hmf_tibb)

dbs_sig_cont_comb_relative$dbs_sig <- factor(dbs_sig_cont_comb_relative$dbs_sig, levels = colnames(dbs_sig_cont_pcawg_df))
dbs_sig_cont_comb_relative$Dataset <- factor(dbs_sig_cont_comb_relative$Dataset, levels = c("Primary", "Metastatic"))


# sum(dbs_sig_cont_comb_relative$Rel_cont == 0)


stat.test <- dbs_sig_cont_comb_relative %>%
  group_by(dbs_sig) %>%
  wilcox_test(Rel_count ~ Dataset) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")


  
stat.test <- stat.test %>%
  add_xy_position(x = "dbs_sig", dodge = 0.8)
stat.test$y.position <- stat.test$y.position*100


dbs_relative_contrib_plot <- dbs_sig_cont_comb_relative %>% ggplot(aes(x = dbs_sig, y = Rel_count*100, color = Dataset)) +
  geom_boxplot(outlier.shape = NA) +
  # stat_compare_means(method = "t.test") +
  # geom_text(data=data.frame(), aes(x=names(meds_ins), y=meds_ins, 
  #                                  label=c(paste0("Mean of ", pcawg_ins_mean), paste0("Mean of ", hmf_ins_mean))), col = "red", size=4) +
  # scale_color_manual(values = c("blue", "green")) +
  # geom_jitter(color="black", size=0.4, alpha=0.5) +
  geom_point(position=position_jitterdodge(), size = 0.4, aes(color = Dataset), alpha = 0.5) +
  scale_color_manual(values = c("#ff0101", "#010dff")) +
  theme_classic() +
  ggtitle("Double-base Substitution Contribution Comparison \n \n") +
  # ylab("Relative Contribution (%) \n") +
  # ylim(0, 110) +
  xlab("\n Lung-specific DBS Signatures") +
  theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.4)) +
  theme(axis.title.x = element_text(face = "italic", size = 10)) +
  theme(axis.title.y = element_text(face = "italic", size = 10)) +
  theme(axis.text.x = element_text(angle = 90, size = 5, vjust = 0.6, face = "bold")) +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  scale_y_continuous(name = "Relative Contribution (%) \n", breaks = c(0,20,40,60,80,100), labels = c(0,20,40,60,80,100)) +
  stat_pvalue_manual(stat.test,  label = "p.adj.signif", tip.length = 0, hide.ns = F)


dbs_relative_contrib_plot


for (i in 1:2){
  if (i == 1) {
    png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/mut-sig/comparison/figs/png/dbs-comparison-relative-contribution-final.png")
    print(dbs_relative_contrib_plot)
    dev.off()
  }
  if (i == 2) {
    pdf(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/mut-sig/comparison/figs/pdf/dbs-comparison-relative-contribution-final.pdf")
    print(dbs_relative_contrib_plot)
    dev.off()
  }
}






#################################################################################################################################################
# Comparing ID sigs



## PCAWG

### Absolute
copy_id_sig_cont_pcawg_df <- id_sig_cont_pcawg_df

copy_id_sig_cont_pcawg_df %<>% mutate(sampleId = rownames(copy_id_sig_cont_pcawg_df)) %>% relocate(where(is.character), .before = where(is.numeric))
rownames(copy_id_sig_cont_pcawg_df) <- 1:nrow(copy_id_sig_cont_pcawg_df)
copy_id_sig_cont_pcawg_tibb <- gather(copy_id_sig_cont_pcawg_df, key = "id_sig", value = "Abs_count", 2:7)
copy_id_sig_cont_pcawg_tibb %<>% mutate(Dataset = "Primary")





### Relative

norm_id_sig_cont_pcawg_df <- id_sig_cont_pcawg_df/rowSums(id_sig_cont_pcawg_df)


norm_id_sig_cont_pcawg_df %<>% mutate(sampleId = rownames(norm_id_sig_cont_pcawg_df)) %>% relocate(where(is.character), .before = where(is.numeric))
rownames(norm_id_sig_cont_pcawg_df) <- 1:nrow(norm_id_sig_cont_pcawg_df)
norm_id_sig_cont_pcawg_tibb <- gather(norm_id_sig_cont_pcawg_df, key = "id_sig", value = "Rel_count", 2:7)
norm_id_sig_cont_pcawg_tibb %<>% mutate(Dataset = "Primary")



## HMF

### Absolute


copy_id_sig_cont_hmf_df <- id_sig_cont_hmf_df_NSC

copy_id_sig_cont_hmf_df %<>% mutate(sampleId = rownames(copy_id_sig_cont_hmf_df)) %>% relocate(where(is.character), .before = where(is.numeric))
rownames(copy_id_sig_cont_hmf_df) <- 1:nrow(copy_id_sig_cont_hmf_df)
copy_id_sig_cont_hmf_tibb <- gather(copy_id_sig_cont_hmf_df, key = "id_sig", value = "Abs_count", 2:7)
copy_id_sig_cont_hmf_tibb %<>% mutate(Dataset = "Metastatic")



### Relative

norm_id_sig_cont_hmf_df <- id_sig_cont_hmf_df_NSC/rowSums(id_sig_cont_hmf_df_NSC)



norm_id_sig_cont_hmf_df %<>% mutate(sampleId = rownames(norm_id_sig_cont_hmf_df)) %>% relocate(where(is.character), .before = where(is.numeric))
rownames(norm_id_sig_cont_hmf_df) <- 1:nrow(norm_id_sig_cont_hmf_df)
norm_id_sig_cont_hmf_tibb <- gather(norm_id_sig_cont_hmf_df, key = "id_sig", value = "Rel_count", 2:7)
norm_id_sig_cont_hmf_tibb %<>% mutate(Dataset = "Metastatic")



## Combine 


### Absolute


id_sig_cont_comb_absolute <- rbind(copy_id_sig_cont_pcawg_tibb, copy_id_sig_cont_hmf_tibb)

id_sig_cont_comb_absolute$id_sig <- factor(id_sig_cont_comb_absolute$id_sig, levels = colnames(id_sig_cont_pcawg_df))
id_sig_cont_comb_absolute$Dataset <- factor(id_sig_cont_comb_absolute$Dataset, levels = c("Primary", "Metastatic"))


# sum(id_sig_cont_comb_absolute$Rel_cont == 0)


stat.test <- id_sig_cont_comb_absolute %>%
  group_by(id_sig) %>%
  wilcox_test(Abs_count ~ Dataset) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")



stat.test <- stat.test %>%
  add_xy_position(x = "id_sig", dodge = 0.8)
stat.test$y.position <- stat.test$y.position


id_absolute_contrib_plot <- id_sig_cont_comb_absolute %>% ggplot(aes(x = id_sig, y = Abs_count, color = Dataset)) +
  geom_boxplot(outlier.shape = NA) +
  # stat_compare_means(method = "t.test") +
  # geom_text(data=data.frame(), aes(x=names(meds_ins), y=meds_ins, 
  #                                  label=c(paste0("Mean of ", pcawg_ins_mean), paste0("Mean of ", hmf_ins_mean))), col = "red", size=4) +
  # scale_color_manual(values = c("blue", "green")) +
  # geom_jitter(color="black", size=0.4, alpha=0.5) +
  geom_point(position=position_jitterdodge(), size = 0.4, aes(color = Dataset), alpha = 0.5) +
  scale_color_manual(values = c("#ff0101", "#010dff")) +
  theme_classic() +
  ggtitle("INDEL Contribution Comparison \n \n") +
  # ylab("Absolute Contribution (%) \n") +
  # ylim(0, 110) +
  xlab("\n Lung-specific ID Signatures") +
  theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.4)) +
  theme(axis.title.x = element_text(face = "italic", size = 10)) +
  theme(axis.title.y = element_text(face = "italic", size = 10)) +
  theme(axis.text.x = element_text(angle = 90, size = 5, vjust = 0.6, face = "bold")) +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  scale_y_continuous(name = "Absolute Contribution \n") +
  stat_pvalue_manual(stat.test,  label = "p.adj.signif", tip.length = 0, hide.ns = F)


id_absolute_contrib_plot


for (i in 1:2){
  if (i == 1) {
    png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/mut-sig/comparison/figs/png/id-comparison-absolute-contribution-final.png")
    print(id_absolute_contrib_plot)
    dev.off()
  }
  if (i == 2) {
    pdf(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/mut-sig/comparison/figs/pdf/id-comparison-absolute-contribution-final.pdf")
    print(id_absolute_contrib_plot)
    dev.off()
  }
}








### Relative

id_sig_cont_comb_relative <- rbind(norm_id_sig_cont_pcawg_tibb, norm_id_sig_cont_hmf_tibb)

id_sig_cont_comb_relative$id_sig <- factor(id_sig_cont_comb_relative$id_sig, levels = colnames(id_sig_cont_pcawg_df))
id_sig_cont_comb_relative$Dataset <- factor(id_sig_cont_comb_relative$Dataset, levels = c("Primary", "Metastatic"))


# sum(id_sig_cont_comb_relative$Rel_cont == 0)


stat.test <- id_sig_cont_comb_relative %>%
  group_by(id_sig) %>%
  wilcox_test(Rel_count ~ Dataset) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")



stat.test <- stat.test %>%
  add_xy_position(x = "id_sig", dodge = 0.8)
stat.test$y.position <- stat.test$y.position*100


id_relative_contrib_plot <- id_sig_cont_comb_relative %>% ggplot(aes(x = id_sig, y = Rel_count*100, color = Dataset)) +
  geom_boxplot(outlier.shape = NA) +
  # stat_compare_means(method = "t.test") +
  # geom_text(data=data.frame(), aes(x=names(meds_ins), y=meds_ins, 
  #                                  label=c(paste0("Mean of ", pcawg_ins_mean), paste0("Mean of ", hmf_ins_mean))), col = "red", size=4) +
  # scale_color_manual(values = c("blue", "green")) +
  # geom_jitter(color="black", size=0.4, alpha=0.5) +
  geom_point(position=position_jitterdodge(), size = 0.4, aes(color = Dataset), alpha = 0.5) +
  scale_color_manual(values = c("#ff0101", "#010dff")) +
  theme_classic() +
  ggtitle("INDEL Contribution Comparison \n \n") +
  # ylab("Relative Contribution (%) \n") +
  # ylim(0, 110) +
  xlab("\n Lung-specific ID Signatures") +
  theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.4)) +
  theme(axis.title.x = element_text(face = "italic", size = 10)) +
  theme(axis.title.y = element_text(face = "italic", size = 10)) +
  theme(axis.text.x = element_text(angle = 90, size = 5, vjust = 0.6, face = "bold")) +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  scale_y_continuous(name = "Relative Contribution (%) \n", breaks = c(0,20,40,60,80,100), labels = c(0,20,40,60,80,100)) +
  stat_pvalue_manual(stat.test,  label = "p.adj.signif", tip.length = 0, hide.ns = F)


id_relative_contrib_plot


for (i in 1:2){
  if (i == 1) {
    png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/mut-sig/comparison/figs/png/id-comparison-relative-contribution-final.png")
    print(id_relative_contrib_plot)
    dev.off()
  }
  if (i == 2) {
    pdf(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/mut-sig/comparison/figs/pdf/id-comparison-relative-contribution-final.pdf")
    print(id_relative_contrib_plot)
    dev.off()
  }
}






