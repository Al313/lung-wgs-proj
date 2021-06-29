
# library(vcfR)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(magrittr)
library(ggplot2)


if (dir.exists("/hpc/cuppen/")){
  
  base_dir <- "/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/"
  devtools::load_all(paste0(base_dir,'/CHORD/processed/scripts_main/mutSigExtractor/'))
  
} else {
  
  library(mutSigExtractor)
}



if (dir.exists("/hpc/cuppen/")){
  hmf_lung_meta <- read.csv(file = "/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/datasets/hmf/hmf-lung-cohort-final.tsv", header = T,
                       sep = "\t", stringsAsFactors = F)
} else {
  hmf_lung_meta <- read.csv(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/datasets/hmf/hmf-lung-cohort-final.tsv", header = T,
                      sep = "\t", stringsAsFactors = F)
}

hmf_lung_meta_NSC <- hmf_lung_meta[hmf_lung_meta$cancer_type == "Non-Small Cell",]


# # @0: Getting all hmf lung vcf files in one list
# 
# 
# 
# hmf_lung_somatic_vcfs <- list()
# 
# 
# 
# for (i in 1:nrow(hmf_lung_meta)){
#   print(i)
#   print(hmf_lung_meta$sampleId[i])
#   sample_id <- hmf_lung_meta$sampleId[i]
#   set_name <- hmf_lung_meta$setName[i]
# 
#   vcf_file = T
# 
#   if (dir.exists("/hpc/cuppen/")){
#     if (file.exists(paste0("/hpc/cuppen/shared_resources/HMF_data/DR-104-update3/somatics/",
#                            set_name, "/purple/", sample_id, ".purple.somatic.vcf.gz"))){
# 
#       path <- paste0("/hpc/cuppen/shared_resources/HMF_data/DR-104-update3/somatics/",
#                      set_name, "/purple/", sample_id, ".purple.somatic.vcf.gz")
# 
#       vcf <- variantsFromVcf(vcf.file = path,
#                                               merge.consecutive = T,
#                                               vcf.filter = "PASS",
#                                               vcf.fields = c("CHROM", "POS", "REF", "ALT", "FILTER", "INFO"))
# 
#     } else {
#       vcf_file = F
#     }
#   } else {
#     if (file.exists(paste0("/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/shared_resources/HMF_data/DR-104-update3/somatics/",
#                     set_name, "/purple/", sample_id, ".purple.somatic.vcf.gz"))){
# 
#       path <- paste0("/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/shared_resources/HMF_data/DR-104-update3/somatics/",
#                      set_name, "/purple/", sample_id, ".purple.somatic.vcf.gz")
# 
#       vcf <- variantsFromVcf(vcf.file = path,
#                                               merge.consecutive = T,
#                                               vcf.filter = "PASS",
#                                               vcf.fields = c("CHROM", "POS", "REF", "ALT", "FILTER", "INFO"))
# 
#     } else {
#       vcf_file = F
#     }
# 
#      }
# 
#   if (vcf_file){
# 
#     selelcted_info_fields <- getInfoValues(vcf$info, keys = c("TNC", "SUBCL"))
# 
# 
#     vcf <- cbind(vcf, selelcted_info_fields)
#     vcf$SUBCL <- as.numeric(vcf$SUBCL)
# 
#     # rename that column, then add another column indicating whether this variant is clonal or subclonal or NA
#     vcf <- vcf %>%
#       # rename(SUBCL = "subclonal_likelihood") %>%
#       mutate(clonality = if_else(SUBCL >= 0.8, F, T)) %>%
#       select(-info)
# 
# 
#     colnames(vcf) <- c("CHROM", "POS", "REF", "ALT", "TNC", "SUBCL_SCORE", "CLONALITY")
# 
#     rownames(vcf) <- paste0(vcf$CHROM, ":", vcf$POS, "_", vcf$REF, "/", vcf$ALT)
# 
#     if (any(duplicated(rownames(vcf)))){
#       vcf <- vcf[isUnique(paste0(vcf$CHROM,":", vcf$POS)),]
#         }
# 
#     hmf_lung_somatic_vcfs[[sample_id]] <- vcf
# 
#   } else {
#     hmf_lung_somatic_vcfs[[sample_id]] <- NA
#   }
# }
# 
# 
# 
# if (dir.exists("/hpc/cuppen/")){
#   saveRDS(hmf_lung_somatic_vcfs, file = "/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/hmf/list-of-hmf-lung-somatic-vcfs.rds")
# } else {
#   saveRDS(hmf_lung_somatic_vcfs, file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/hmf/list-of-hmf-lung-somatic-vcfs.rds")
# }
# 
# 
# 
# 
# toy_hmf_lung_somatic_vcfs <- hmf_lung_somatic_vcfs[1:10]
# 
# 
# if (dir.exists("/hpc/cuppen/")){
#   saveRDS(toy_hmf_lung_somatic_vcfs, file = "/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/hmf/toy-list-of-hmf-lung-somatic-vcfs.rds")
# } else {
#   saveRDS(toy_hmf_lung_somatic_vcfs, file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/hmf/toy-list-of-hmf-lung-somatic-vcfs.rds")
# }



if (dir.exists("/hpc/cuppen/")){
  hmf_lung_somatic_vcfs <- readRDS(file = "/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/hmf/list-of-hmf-lung-somatic-vcfs.rds")
} else {
  hmf_lung_somatic_vcfs <- readRDS(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/hmf/toy-list-of-hmf-lung-somatic-vcfs.rds")
}

#========================================================================================================================================================



# # @1:  count matrix for somatic SNVs, MNVs, and indels of hmf lung cohort
# 
# 
# count_matrix_hmf <- matrix(nrow = length(hmf_lung_somatic_vcfs), ncol = 5)
# 
# colnames(count_matrix_hmf) <- c("all", "SNVs", "MNVs", "inDEL", "INdel")
# rownames(count_matrix_hmf) <- names(hmf_lung_somatic_vcfs)
# 
# 
# for (i in 1:length(hmf_lung_somatic_vcfs)){
#   df <- hmf_lung_somatic_vcfs[[i]]
#   count_matrix_hmf[i,1] <- nrow(df)
#   count_matrix_hmf[i,2] <- nrow(df[nchar(df$REF) == 1 & nchar(df$ALT) == 1,])
#   count_matrix_hmf[i,3] <- nrow(df[nchar(df$REF) == nchar(df$ALT) & nchar(df$REF) > 1,])
#   count_matrix_hmf[i,4] <- nrow(df[nchar(df$REF) > nchar(df$ALT) & nchar(df$REF),])
#   count_matrix_hmf[i,5] <- nrow(df[nchar(df$REF) < nchar(df$ALT) & nchar(df$REF),])
# 
# }
# 
# 
# count_df_hmf <- as.data.frame(count_matrix_hmf)
# count_df_hmf["sampleId"] <- rownames(count_df_hmf)
# rownames(count_df_hmf) <- 1:nrow(count_df_hmf)
# count_df_hmf <- count_df_hmf[,c(6,2,3,4,5,1)]
# 
# if (dir.exists("/hpc/cuppen/")){
#   saveRDS(count_df_hmf, file = "/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/hmf/count_df_hmf.rds")
# } else {
#   saveRDS(count_df_hmf, file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/hmf/count_df_hmf.rds")
# }



if (dir.exists("/hpc/cuppen/")){
  count_df_hmf <- readRDS(file = "/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/hmf/count_df_hmf.rds")
} else {
  count_df_hmf <- readRDS(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/hmf/count_df_hmf.rds")
}

count_df_hmf_NSC <- count_df_hmf[count_df_hmf$sampleId %in% hmf_lung_meta_NSC$sampleId,]


#========================================================================================================================================================




# # @2: count matrix for simple structural variants such as deletions, duplications, translocations, and inversions of hmf lung cohort
# 
# simple_sv_events <- c("Deletion (> 100kb)", "Deletion (< 100kb)", "Duplication (> 100kb)", "Duplication (< 100kb)", "Inversion")
# sv_count_matrix_hmf <- matrix(nrow = nrow(hmf_lung_meta), ncol = length(simple_sv_events))
# colnames(sv_count_matrix_hmf) <- simple_sv_events
# rownames(sv_count_matrix_hmf) <- hmf_lung_meta$sampleId
# 
# 
# 
# for (i in 1:nrow(hmf_lung_meta)){
# 
#   print(i)
#   print(hmf_lung_meta$sampleId[i])
#   sample_id <- hmf_lung_meta$sampleId[i]
#   set_name <- hmf_lung_meta$setName[i]
# 
#   if (dir.exists("/hpc/cuppen/")){
#     vcf <- read.csv(file = paste0("/hpc/cuppen/shared_resources/HMF_data/DR-104-update3/somatics/",
#                                   set_name, "/linx14/", sample_id, ".linx.vis_sv_data.tsv"), header = T, sep = "\t", stringsAsFactors = F)
#   } else {
#     vcf <- read.csv(file = paste0("/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/shared_resources/HMF_data/DR-104-update3/somatics/",
#                                   set_name, "/linx14/", sample_id, ".linx.vis_sv_data.tsv"), header = T, sep = "\t", stringsAsFactors = F)
#   }
# 
# 
#   vcf <- mutate(vcf, length = abs(PosEnd - PosStart))
# 
# 
# 
#   vcf1 <- vcf[vcf$ResolvedType == "DEL" & vcf$length > 100000,]
#   sv_count_matrix_hmf[i,1] <- nrow(vcf1)
# 
#   vcf2 <- vcf[vcf$ResolvedType == "DEL" & vcf$length < 100000,]
#   sv_count_matrix_hmf[i,2] <- nrow(vcf2)
# 
#   vcf3 <- vcf[vcf$ResolvedType == "DUP" & vcf$length > 100000,]
#   sv_count_matrix_hmf[i,3] <- nrow(vcf3)
# 
#   vcf4 <- vcf[vcf$ResolvedType == "DUP" & vcf$length < 100000,]
#   sv_count_matrix_hmf[i,4] <- nrow(vcf4)
# 
#   vcf5 <- vcf[vcf$ResolvedType == "INV",]
#   sv_count_matrix_hmf[i,5] <- nrow(vcf5)
# 
# }
# 
# 
# 
# if (dir.exists("/hpc/cuppen/")){
#   saveRDS(sv_count_matrix_hmf, file = "/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/hmf/simple_sv_count_df_hmf.rds")
# } else {
#   saveRDS(sv_count_matrix_hmf, file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/hmf/simple_sv_count_df_hmf.rds")
# }



if (dir.exists("/hpc/cuppen/")){
  sv_count_matrix_hmf <- readRDS(file = "/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/hmf/simple_sv_count_df_hmf.rds")
} else {
  sv_count_matrix_hmf <- readRDS(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/hmf/simple_sv_count_df_hmf.rds")
}

sv_count_matrix_hmf_NSC <- sv_count_matrix_hmf[rownames(sv_count_matrix_hmf) %in% hmf_lung_meta_NSC$sampleId,]


#========================================================================================================================================================




# # @2.1: length matrix for simple structural variants such as deletions, duplications, and inversions of hmf lung cohort
# 
# sv_length_df_hmf <- data.frame()
# sv_length_dfs_hmf <- data.frame()
# 
# 
# 
# for (i in 1:nrow(hmf_lung_meta)){
# 
#   print(i)
#   print(hmf_lung_meta$sampleId[i])
#   sample_id <- hmf_lung_meta$sampleId[i]
#   set_name <- hmf_lung_meta$setName[i]
# 
#   if (dir.exists("/hpc/cuppen/")){
#     vcf <- read.csv(file = paste0("/hpc/cuppen/shared_resources/HMF_data/DR-104-update3/somatics/",
#                                   set_name, "/linx14/", sample_id, ".linx.vis_sv_data.tsv"), header = T, sep = "\t", stringsAsFactors = F)
#   } else {
#     vcf <- read.csv(file = paste0("/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/shared_resources/HMF_data/DR-104-update3/somatics/",
#                                   set_name, "/linx14/", sample_id, ".linx.vis_sv_data.tsv"), header = T, sep = "\t", stringsAsFactors = F)
#   }
# 
#   vcf <- mutate(vcf, length = abs(PosEnd - PosStart))
# 
#   sv_length_df_hmf <- vcf[vcf$ResolvedType %in% c("DEL", "DUP", "INV"), c("SampleId", "ResolvedType", "length")]
#   rownames(sv_length_df_hmf) <- (nrow(sv_length_dfs_hmf) +1 ):(nrow(sv_length_dfs_hmf) + nrow(sv_length_df_hmf))
#   sv_length_dfs_hmf <- rbind(sv_length_dfs_hmf, sv_length_df_hmf)
# 
# }
# 
# 
# 
# if (dir.exists("/hpc/cuppen/")){
#   saveRDS(sv_length_dfs_hmf, file = "/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/hmf/sv_length_dfs_hmf.rds")
# } else {
#   saveRDS(sv_length_dfs_hmf, file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/hmf/sv_length_dfs_hmf.rds")
# }


if (dir.exists("/hpc/cuppen/")){
  sv_length_dfs_hmf <- readRDS(file = "/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/hmf/sv_length_dfs_hmf.rds")
} else {
  sv_length_dfs_hmf <- readRDS(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/hmf/sv_length_dfs_hmf.rds")
}



#========================================================================================================================================================

# # @3: clonality of mutations
# 
# clonality_summary_hmf <- data.frame(sampleId = NA, Mut_type = NA, Clonal = NA, Non_clonal = NA, Unknown = NA)
# 
# 
# for (i in 1:length(hmf_lung_somatic_vcfs)){
# 
#   print(i)
#   j <- (i-1) * 5 + 1
# 
#   clonality_summary_hmf[j:(j+4),1] <- names(hmf_lung_somatic_vcfs)[i]
#   clonality_summary_hmf[j:(j+4),2] <- c("All", "SNVs", "MNVs", "Ins", "Del")
#   for (type in c("SNVs", "MNVs", "Ins", "Del", "All")) {
# 
#     if (type == "All") {
#       tmp_df <- hmf_lung_somatic_vcfs[[i]]
#       clonality_summary_hmf[j,3] <- as.vector(table(tmp_df[,"CLONALITY"],  useNA="always")["TRUE"])
#       clonality_summary_hmf[j,4] <- as.vector(table(tmp_df[,"CLONALITY"],  useNA="always")["FALSE"])
#       clonality_summary_hmf[j,5] <- as.vector(tail(table(tmp_df[,"CLONALITY"],  useNA="always"),1))
#     } else if (type == "SNVs") {
#       tmp_df <- hmf_lung_somatic_vcfs[[i]]
#       clonality_summary_hmf[j+1,3] <- as.vector(table(tmp_df[nchar(tmp_df$REF) == 1 & nchar(tmp_df$ALT) == 1 ,"CLONALITY"],  useNA="always")["TRUE"])
#       clonality_summary_hmf[j+1,4] <- as.vector(table(tmp_df[nchar(tmp_df$REF) == 1 & nchar(tmp_df$ALT) == 1 ,"CLONALITY"],  useNA="always")["FALSE"])
#       clonality_summary_hmf[j+1,5] <- as.vector(tail(table(tmp_df[nchar(tmp_df$REF) == 1 & nchar(tmp_df$ALT) == 1 ,"CLONALITY"],  useNA="always"),1))
#     } else if (type == "MNVs") {
#       tmp_df <- hmf_lung_somatic_vcfs[[i]]
#       clonality_summary_hmf[j+2,3] <- as.vector(table(tmp_df[nchar(tmp_df$REF) == nchar(tmp_df$ALT) & nchar(tmp_df$REF) > 1 ,"CLONALITY"],  useNA="always")["TRUE"])
#       clonality_summary_hmf[j+2,4] <- as.vector(table(tmp_df[nchar(tmp_df$REF) == nchar(tmp_df$ALT) & nchar(tmp_df$REF) > 1 ,"CLONALITY"],  useNA="always")["FALSE"])
#       clonality_summary_hmf[j+2,5] <- as.vector(tail(table(tmp_df[nchar(tmp_df$REF) == nchar(tmp_df$ALT) & nchar(tmp_df$REF) > 1 ,"CLONALITY"],  useNA="always"),1))
#     } else if (type == "Ins") {
#       tmp_df <- hmf_lung_somatic_vcfs[[i]]
#       clonality_summary_hmf[j+3,3] <- as.vector(table(tmp_df[nchar(tmp_df$REF) < nchar(tmp_df$ALT) & nchar(tmp_df$REF) ,"CLONALITY"],  useNA="always")["TRUE"])
#       clonality_summary_hmf[j+3,4] <- as.vector(table(tmp_df[nchar(tmp_df$REF) < nchar(tmp_df$ALT) & nchar(tmp_df$REF) ,"CLONALITY"],  useNA="always")["FALSE"])
#       clonality_summary_hmf[j+3,5] <- as.vector(tail(table(tmp_df[nchar(tmp_df$REF) < nchar(tmp_df$ALT) & nchar(tmp_df$REF) ,"CLONALITY"],  useNA="always"),1))
#     } else if (type == "Del") {
#       tmp_df <- hmf_lung_somatic_vcfs[[i]]
#       clonality_summary_hmf[j+4,3] <- as.vector(table(tmp_df[nchar(tmp_df$REF) > nchar(tmp_df$ALT) & nchar(tmp_df$REF) ,"CLONALITY"],  useNA="always")["TRUE"])
#       clonality_summary_hmf[j+4,4] <- as.vector(table(tmp_df[nchar(tmp_df$REF) > nchar(tmp_df$ALT) & nchar(tmp_df$REF) ,"CLONALITY"],  useNA="always")["FALSE"])
#       clonality_summary_hmf[j+4,5] <- as.vector(tail(table(tmp_df[nchar(tmp_df$REF) > nchar(tmp_df$ALT) & nchar(tmp_df$REF) ,"CLONALITY"],  useNA="always"),1))
#     }
#   }
# }
# 
# 
# if (dir.exists("/hpc/cuppen/")){
#   saveRDS(clonality_summary_hmf, file = "/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/hmf/clonality_summary_hmf.rds")
# } else {
#   saveRDS(clonality_summary_hmf, file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/hmf/clonality_summary_hmf.rds")
# }




if (dir.exists("/hpc/cuppen/")){
  clonality_summary_hmf <- readRDS(file = "/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/hmf/clonality_summary_hmf.rds")
} else {
  clonality_summary_hmf <- readRDS(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/hmf/clonality_summary_hmf.rds")
}

clonality_summary_hmf_NSC <- clonality_summary_hmf[clonality_summary_hmf$sampleId %in% hmf_lung_meta_NSC$sampleId,]


#========================================================================================================================================================

# # @4: ploidy info of samples
# 
# 
# hmf_ploidy_info <- list()
# 
# for (i in 1:nrow(hmf_lung_meta)){
# 
#   print(i)
#   print(hmf_lung_meta$sampleId[i])
#   sample_id <- hmf_lung_meta$sampleId[i]
#   set_name <- hmf_lung_meta$setName[i]
# 
#   if (dir.exists("/hpc/cuppen/")){
#     ploidy_info <- read.csv(file = paste0("/hpc/cuppen/shared_resources/HMF_data/DR-104-update3/somatics/",
#                                   set_name, "/purple/", sample_id, ".purple.purity.tsv"), header = T, sep = "\t", stringsAsFactors = F)
#   } else {
#     ploidy_info <- read.csv(file = paste0("/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/shared_resources/HMF_data/DR-104-update3/somatics/",
#                                   set_name, "/purple/", sample_id, ".purple.purity.tsv"), header = T, sep = "\t", stringsAsFactors = F)
#   }
# 
# 
#   hmf_ploidy_info[[sample_id]] <- ploidy_info
# }
# 
# 
# 
# if (dir.exists("/hpc/cuppen/")){
#   saveRDS(hmf_ploidy_info, file = "/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/hmf/hmf-ploidy-info.rds")
# } else {
#   saveRDS(hmf_ploidy_info, file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/hmf/hmf-ploidy-info.rds")
# }



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


#========================================================================================================================================================


# clustering

head(count_df_hmf_NSC)

count_df_hmf_NSC_norm <- count_df_hmf_NSC[,2:5]#/count_df_hmf_NSC[,6]
rownames(count_df_hmf_NSC_norm) <- count_df_hmf_NSC$sampleId
count_df_hmf_NSC_norm <- t(count_df_hmf_NSC_norm)


sv_count_df_hmf_NSC <- as.data.frame(sv_count_matrix_hmf_NSC)

sv_count_df_hmf_NSC_norm <- sv_count_df_hmf_NSC[,1:5]#/rowSums(sv_count_df_hmf_NSC[,1:5])
rownames(sv_count_df_hmf_NSC_norm) <- rownames(sv_count_df_hmf_NSC)
sv_count_df_hmf_NSC_norm <- t(sv_count_df_hmf_NSC_norm)

ncol(sv_count_df_hmf_NSC_norm)
mut_count_combined <- rbind(count_df_hmf_NSC_norm, sv_count_df_hmf_NSC_norm)
hc.sample.combined <- hclust(dist(t(mut_count_combined)), method = "complete")
sample_order_comb <- colnames(count_df_hmf_NSC_norm)[hc.sample.combined$order]



# Mutation distribution per sample

count_tibb_hmf <- count_df_hmf_NSC %>%
  gather(key = "Mut_type", value = "Count", 2:5)

count_tibb_hmf %<>% select(-all)

count_tibb_hmf$Mut_type <- factor(count_tibb_hmf$Mut_type, levels = rev(c("inDEL", "INdel", "MNVs", "SNVs")))

count_tibb_hmf$sampleId <- factor(count_tibb_hmf$sampleId, levels = sample_order_comb)

cols <- brewer.pal(8, name = "Dark2")
tmb_plot <- count_tibb_hmf %>% ggplot(aes(x = sampleId, y = Count, fill = Mut_type, color = Mut_type)) +
  geom_bar(position="stack", stat="identity", width = 1, size = 0.1) +
  scale_fill_manual(values = cols[c(1,2,7,8)], labels = c("SNVs", "MNVs", "Small Insertions", "Small Deletions")) +
  scale_color_manual(values = rep('black', times =5), guide = "none") +
  theme_bw() +
  theme(
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_line(color = "black", size = 0.5, linetype = "dashed")
  ) +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.text.x = element_blank()) +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  scale_y_continuous(name = "Counts \n",
                     breaks = c(0, 25000, 50000, 75000, 100000, 150000, 200000, 250000, 300000, 350000, 400000),
                     labels = c(0, "25,000", "50,000", "75,000", "100,000", "150,000", "200000", "250000", "300000", "350000", "400000")) +
  labs(x = NULL, y = "Counts", fill = "Mutation Type") +
  ggtitle("Absolute Mutation Counts\n \n") +
  theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.35))
# scale_fill_discrete(labels = c("SNVs", "MNVs", "Small Insertions", "Small Deletions"))





# SV event distributions per sample

sv_count_df_hmf_NSC %<>% mutate(sampleId = rownames(sv_count_df_hmf_NSC))
rownames(sv_count_df_hmf_NSC) <- 1:nrow(sv_count_df_hmf_NSC)

sv_count_tibb_hmf <- sv_count_df_hmf_NSC %>% gather(key = "Mut_type", value = "Count", 1:5)
sv_count_tibb_hmf$Mut_type <- factor(sv_count_tibb_hmf$Mut_type, levels = c("Deletion (< 100kb)", "Deletion (> 100kb)", "Duplication (< 100kb)", "Duplication (> 100kb)", "Inversion"))

sv_count_tibb_hmf$sampleId <- factor(sv_count_tibb_hmf$sampleId, levels = sample_order_comb)


sv_plot <- sv_count_tibb_hmf %>% ggplot(aes(x = sampleId, y = Count, fill = Mut_type, color = Mut_type)) +
  geom_bar(position="stack", stat="identity", width = 1, size = 0.1) +
  scale_fill_manual(values = brewer.pal(5, "Accent")) +
  scale_color_manual(values = rep('black', times =5), guide = "none") +
  theme_bw() +
  theme(
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_line(color = "black", size = 0.5, linetype = "dashed")
  ) +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.text.x = element_blank()) +
  labs(x = "HMF NSCLC Samples (n = 316)", y = "Counts", fill = "SV Type")


# Combine



tmb_plot_grob <- ggplotGrob(tmb_plot)
sv_plot_grob <- ggplotGrob(sv_plot)

hmf_mut_dist <- cowplot::plot_grid(tmb_plot_grob, sv_plot_grob, align = "v", rel_heights = c(1.5,1),nrow = 2, ncol = 1)



# save


for (i in 1:2){
  if (i == 1) {
<<<<<<< HEAD
    png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/hmf/figs/png/mut-count-per-sample-final.png", width = 960)
=======
    png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/hmf/figs/png/mut-count-per-sample-nsc.png", width = 960)
>>>>>>> d4b2ba7d869eab81f62853d96bcbc0ac881b0e01
    print(hmf_mut_dist)
    dev.off()
  }
  if (i == 2) {
<<<<<<< HEAD
    pdf(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/hmf/figs/pdf/mut-count-per-sample-final.pdf", width = 14)
=======
    pdf(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/hmf/figs/pdf/mut-count-per-sample-nsc.pdf", width = 14)
>>>>>>> d4b2ba7d869eab81f62853d96bcbc0ac881b0e01
    print(hmf_mut_dist)
    dev.off()
  }
}




# checking

all(unique(sv_count_tibb_hmf$sampleId) == unique(count_tibb_hmf$sampleId))


count_tibb_hmf[count_tibb_hmf$Count == max(count_tibb_hmf$Count),]
count_tibb_hmf[count_tibb_hmf$sampleId == "DO27747",]

sv_count_tibb_hmf[sv_count_tibb_hmf$Count == max(sv_count_tibb_hmf$Count),]
sv_count_tibb_hmf[sv_count_tibb_hmf$sampleId == "DO25189",]

######################################################################################################################################################
# SV event size distribution

del <- density(sv_length_dfs_hmf[sv_length_dfs_hmf$ResolvedType == "DEL" & sv_length_dfs_hmf$length < 100000,"length"], width = 10000, from = 0)
dup <- density(sv_length_dfs_hmf[sv_length_dfs_hmf$ResolvedType == "DUP" & sv_length_dfs_hmf$length < 100000,"length"], width = 10000, from = 0)
plot(del)
plot(dup)

table(sv_length_dfs_hmf$ResolvedType)


######################################################################################################################################################
# Clonality


# head(clonality_summary_hmf)

clonality_summary_hmf_NSC$Tot <- rowSums(clonality_summary_hmf_NSC[,3:5], na.rm = T)
clonality_summary_hmf_NSC %<>% gather(key = "Clonality_stat", value = "Counts", 3:5)
clonality_summary_hmf_NSC$Mut_type <- factor(clonality_summary_hmf_NSC$Mut_type, levels = c("All", "SNVs", "MNVs", "Ins", "Del"))
clonality_summary_hmf_NSC$Clonality_stat <- factor(clonality_summary_hmf_NSC$Clonality_stat, levels = c("Clonal", "Non_clonal", "Unknown"))


options(scipen=999)



clonality_plot_summary_hmf <- clonality_summary_hmf_NSC %>% ggplot(aes(x = sampleId, y = Counts, fill = Clonality_stat)) + facet_grid(Mut_type ~ ., scales = 'free') +
  geom_bar(position="stack", stat="identity", width = 1) +
  scale_fill_manual(values = c("red", "black", "yellow")) +
  ggtitle("Clonality States of Different Mutation Types \n") +
  labs(x = "HMF NSCLC Samples (n = 316)", y = "Frequency", fill = "Clonality") +
  theme_bw() +
  theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.5)) +
  theme(axis.ticks.x = element_blank()) + theme(axis.text.x = element_blank()) + theme(panel.grid.major.y = element_line(color = "black", size = 0.5, linetype = "dashed"))




# Save

for (i in 1:2){
  if (i == 1) {
<<<<<<< HEAD
    png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/hmf/figs/png/clonality-plot-summary-hmf-final.png", width = 960)
=======
    png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/hmf/figs/png/clonality-plot-summary-hmf-nsc.png", width = 960)
>>>>>>> d4b2ba7d869eab81f62853d96bcbc0ac881b0e01
    print(clonality_plot_summary_hmf)
    dev.off()
  }
  if (i == 2) {
<<<<<<< HEAD
    pdf(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/hmf/figs/pdf/clonality-plot-summary-hmf-final.pdf", width = 14)
=======
    pdf(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/hmf/figs/pdf/clonality-plot-summary-hmf-nsc.pdf", width = 14)
>>>>>>> d4b2ba7d869eab81f62853d96bcbc0ac881b0e01
    print(clonality_plot_summary_hmf)
    dev.off()
  }
}



