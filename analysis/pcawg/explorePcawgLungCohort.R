
# library(vcfR)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(magrittr)


if (dir.exists("/hpc/cuppen/")){
  
  base_dir <- "/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/"
  devtools::load_all(paste0(base_dir,'/CHORD/processed/scripts_main/mutSigExtractor/'))
  
  } else {

  library(mutSigExtractor)
}


if (dir.exists("/hpc/cuppen/")){
  pcawg_meta <- readRDS(file = "/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/datasets/pcawg/metadata-per-donor.Rds")
} else {
  pcawg_meta <- readRDS(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/datasets/pcawg/metadata-per-donor.Rds")
}

pcawg_lung_meta <- pcawg_meta[pcawg_meta$organ_system == "LUNG & BRONCHUS",]
pcawg_lung_meta <- pcawg_lung_meta[pcawg_lung_meta$icgc_donor_id != "DO23717",]


# @0: Getting all pcawg lung vcf files in one list


# pcawg_lung_somatic_vcfs <- list()
# for (i in 1:nrow(pcawg_lung_meta)){
#   print(i)
#   print(pcawg_lung_meta$icgc_donor_id[i])
#   sample_id <- pcawg_lung_meta$icgc_donor_id[i]
# 
#   vcf_file = T
# 
#   if (dir.exists("/hpc/cuppen/")){
#     if (file.exists(paste0("/hpc/cuppen/shared_resources/PCAWG/pipeline5/per-donor/",
#                            sample_id, "-from-jar/purple53/", sample_id, "T.purple.somatic.vcf.gz"))){
# 
#       path <- paste0("/hpc/cuppen/shared_resources/PCAWG/pipeline5/per-donor/",
#                      sample_id, "-from-jar/purple53/", sample_id, "T.purple.somatic.vcf.gz")
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
#     if (file.exists(paste0("/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/shared_resources/PCAWG/pipeline5/per-donor/",
#                            sample_id, "-from-jar/purple53/", sample_id, "T.purple.somatic.vcf.gz"))){
# 
#       path <- paste0("/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/shared_resources/PCAWG/pipeline5/per-donor/",
#                      sample_id, "-from-jar/purple53/", sample_id, "T.purple.somatic.vcf.gz")
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
#     pcawg_lung_somatic_vcfs[[sample_id]] <- vcf
# 
#   } else {
#     pcawg_lung_somatic_vcfs[[sample_id]] <- NA
#   }
# }
# 
# 
# 
# if (dir.exists("/hpc/cuppen/")){
#   saveRDS(pcawg_lung_somatic_vcfs, file = "/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/pcawg/list-of-pcawg-lung-somatic-vcfs.rds")
# } else {
#   saveRDS(pcawg_lung_somatic_vcfs, file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/pcawg/list-of-pcawg-lung-somatic-vcfs.rds")
# }
# 
# 
# 
# 
# toy_pcawg_lung_somatic_vcfs <- pcawg_lung_somatic_vcfs[1:10]
# 
# 
# if (dir.exists("/hpc/cuppen/")){
#   saveRDS(toy_pcawg_lung_somatic_vcfs, file = "/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/pcawg/toy-list-of-pcawg-lung-somatic-vcfs.rds")
# } else {
#   saveRDS(toy_pcawg_lung_somatic_vcfs, file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/pcawg/toy-list-of-pcawg-lung-somatic-vcfs.rds")
# }



if (dir.exists("/hpc/cuppen/")){
  pcawg_lung_somatic_vcfs <- readRDS(file = "/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/pcawg/list-of-pcawg-lung-somatic-vcfs.rds")
} else {
  pcawg_lung_somatic_vcfs <- readRDS(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/pcawg/toy-list-of-pcawg-lung-somatic-vcfs.rds")
}

#========================================================================================================================================================


# @1: Removing PCAWG samples that are not processed by HMF pipeline (but have a directory in per-donor directoy)

# print(names(pcawg_lung_somatic_vcfs[as.vector(is.na(pcawg_lung_somatic_vcfs))]))

# One of the 89 samples did not have purple somatic file (DO23717) this sample will be put aside for the follow-on analysis.

# pcawg_lung_somatic_vcfs <- pcawg_lung_somatic_vcfs[as.vector(!is.na(pcawg_lung_somatic_vcfs))]
# 
# str(pcawg_lung_somatic_vcfs)
# 
# if (dir.exists("/hpc/cuppen/")){
#   saveRDS(pcawg_lung_somatic_vcfs, file = "/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/pcawg/list-of-pcawg-lung-somatic-vcfs.rds")
# } else {
#   saveRDS(pcawg_lung_somatic_vcfs, file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/pcawg/list-of-pcawg-lung-somatic-vcfs.rds")
# }


#========================================================================================================================================================



# @2:  count matrix for somatic SNVs, MNVs, and indels of PCAWG lung cohort


# count_matrix_pcawg <- matrix(nrow = length(pcawg_lung_somatic_vcfs), ncol = 5)
# 
# colnames(count_matrix_pcawg) <- c("all", "SNVs", "MNVs", "inDEL", "INdel")
# rownames(count_matrix_pcawg) <- names(pcawg_lung_somatic_vcfs)
# 
# 
# for (i in 1:length(pcawg_lung_somatic_vcfs)){
#   df <- pcawg_lung_somatic_vcfs[[i]]
#   count_matrix_pcawg[i,1] <- nrow(df)
#   count_matrix_pcawg[i,2] <- nrow(df[nchar(df$REF) == 1 & nchar(df$ALT) == 1,])
#   count_matrix_pcawg[i,3] <- nrow(df[nchar(df$REF) == nchar(df$ALT) & nchar(df$REF) > 1,])
#   count_matrix_pcawg[i,4] <- nrow(df[nchar(df$REF) > nchar(df$ALT) & nchar(df$REF),])
#   count_matrix_pcawg[i,5] <- nrow(df[nchar(df$REF) < nchar(df$ALT) & nchar(df$REF),])
# 
# }
# 
# 
# count_df_pcawg <- as.data.frame(count_matrix_pcawg)
# count_df_pcawg["sampleId"] <- rownames(count_df_pcawg)
# rownames(count_df_pcawg) <- 1:nrow(count_df_pcawg)
# count_df_pcawg <- count_df_pcawg[,c(6,2,3,4,5,1)]
# 
# if (dir.exists("/hpc/cuppen/")){
#   saveRDS(count_df_pcawg, file = "/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/pcawg/count_df_pcawg.rds")
# } else {
#   saveRDS(count_df_pcawg, file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/pcawg/count_df_pcawg.rds")
# }



if (dir.exists("/hpc/cuppen/")){
  count_df_pcawg <- readRDS(file = "/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/pcawg/count_df_pcawg.rds")
} else {
  count_df_pcawg <- readRDS(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/pcawg/count_df_pcawg.rds")
}



#========================================================================================================================================================


# @3: count matrix for simple structural variants such as deletions, duplications, translocations, and inversions of PCAWG lung cohort

# simple_sv_events <- c("Deletion (> 100kb)", "Deletion (< 100kb)", "Duplication (> 100kb)", "Duplication (< 100kb)", "Inversion")
# sv_count_matrix_pcawg <- matrix(nrow = nrow(pcawg_lung_meta), ncol = length(simple_sv_events))
# colnames(sv_count_matrix_pcawg) <- simple_sv_events
# rownames(sv_count_matrix_pcawg) <- pcawg_lung_meta$icgc_donor_id
# 
# 
# 
# for (i in 1:nrow(pcawg_lung_meta)){
# 
#   print(i)
#   print(pcawg_lung_meta$icgc_donor_id[i])
#   if (dir.exists("/hpc/cuppen/")){
#     vcf <- read.csv(file = paste0("/hpc/cuppen/shared_resources/PCAWG/pipeline5/per-donor/",
#                                   pcawg_lung_meta$icgc_donor_id[i], "-from-jar/linx14/", pcawg_lung_meta$icgc_donor_id[i], "T.linx.vis_sv_data.tsv"), header = T, sep = "\t", stringsAsFactors = F)
#   } else {
#     vcf <- read.csv(file = paste0("/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/shared_resources/PCAWG/pipeline5/per-donor/",
#                                   pcawg_lung_meta$icgc_donor_id[i], "-from-jar/linx14/", pcawg_lung_meta$icgc_donor_id[i], "T.linx.vis_sv_data.tsv"), header = T, sep = "\t", stringsAsFactors = F)
#   }
# 
# 
#   vcf <- mutate(vcf, length = abs(PosEnd - PosStart))
# 
# 
# 
#   vcf1 <- vcf[vcf$ResolvedType == "DEL" & vcf$length > 100000,]
#   sv_count_matrix_pcawg[i,1] <- nrow(vcf1)
# 
#   vcf2 <- vcf[vcf$ResolvedType == "DEL" & vcf$length < 100000,]
#   sv_count_matrix_pcawg[i,2] <- nrow(vcf2)
# 
#   vcf3 <- vcf[vcf$ResolvedType == "DUP" & vcf$length > 100000,]
#   sv_count_matrix_pcawg[i,3] <- nrow(vcf3)
# 
#   vcf4 <- vcf[vcf$ResolvedType == "DUP" & vcf$length < 100000,]
#   sv_count_matrix_pcawg[i,4] <- nrow(vcf4)
# 
#   vcf5 <- vcf[vcf$ResolvedType == "INV",]
#   sv_count_matrix_pcawg[i,5] <- nrow(vcf5)
# 
# }
# 
# 
# 
# if (dir.exists("/hpc/cuppen/")){
#   saveRDS(sv_count_matrix_pcawg, file = "/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/pcawg/simple_sv_count_df_pcawg.rds")
# } else {
#   saveRDS(sv_count_matrix_pcawg, file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/pcawg/simple_sv_count_df_pcawg.rds")
# }



if (dir.exists("/hpc/cuppen/")){
  sv_count_matrix_pcawg <- readRDS(file = "/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/pcawg/simple_sv_count_df_pcawg.rds")
} else {
  sv_count_matrix_pcawg <- readRDS(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/pcawg/simple_sv_count_df_pcawg.rds")
}



#========================================================================================================================================================


# @3.1: length matrix for simple structural variants such as deletions, duplications, and inversions of PCAWG lung cohort

# sv_length_df_pcawg <- data.frame()
# sv_length_dfs_pcawg <- data.frame()
# 
# 
# 
# for (i in 1:nrow(pcawg_lung_meta)){
# 
#   print(i)
#   print(pcawg_lung_meta$icgc_donor_id[i])
#   i <- 2
#   
#   if (dir.exists("/hpc/cuppen/")){
#     vcf <- read.csv(file = paste0("/hpc/cuppen/shared_resources/PCAWG/pipeline5/per-donor/",
#                                   pcawg_lung_meta$icgc_donor_id[i], "-from-jar/linx14/", pcawg_lung_meta$icgc_donor_id[i], "T.linx.vis_sv_data.tsv"), header = T, sep = "\t", stringsAsFactors = F)
#   } else {
#     vcf <- read.csv(file = paste0("/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/shared_resources/PCAWG/pipeline5/per-donor/",
#                                   pcawg_lung_meta$icgc_donor_id[i], "-from-jar/linx14/", pcawg_lung_meta$icgc_donor_id[i], "T.linx.vis_sv_data.tsv"), header = T, sep = "\t", stringsAsFactors = F)
#   }
#   
#   vcf <- mutate(vcf, length = abs(PosEnd - PosStart))
#   
#   sv_length_df_pcawg <- vcf[vcf$ResolvedType %in% c("DEL", "DUP", "INV"), c("SampleId", "ResolvedType", "length")]
#   rownames(sv_length_df_pcawg) <- (nrow(sv_length_dfs_pcawg) +1 ):(nrow(sv_length_dfs_pcawg) + nrow(sv_length_df_pcawg))
#   sv_length_dfs_pcawg <- rbind(sv_length_dfs_pcawg, sv_length_df_pcawg)
# 
# }
# 
# 
# 
# if (dir.exists("/hpc/cuppen/")){
#   saveRDS(sv_length_dfs_pcawg, file = "/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/pcawg/sv_length_dfs_pcawg.rds")
# } else {
#   saveRDS(sv_length_dfs_pcawg, file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/pcawg/sv_length_dfs_pcawg.rds")
# }


if (dir.exists("/hpc/cuppen/")){
  sv_length_dfs_pcawg <- readRDS(file = "/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/pcawg/sv_length_dfs_pcawg.rds")
} else {
  sv_length_dfs_pcawg <- readRDS(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/pcawg/sv_length_dfs_pcawg.rds")
}




#========================================================================================================================================================

# @4: clonality of mutations

# clonality_summary_pcawg <- data.frame(sampleId = NA, Mut_type = NA, Clonal = NA, Non_clonal = NA, Unknown = NA)
# 
# 
# for (i in 1:length(pcawg_lung_somatic_vcfs)){
#   
#   print(i)
#   j <- (i-1) * 5 + 1
#   
#   clonality_summary_pcawg[j:(j+4),1] <- names(pcawg_lung_somatic_vcfs)[i]
#   clonality_summary_pcawg[j:(j+4),2] <- c("All", "SNVs", "MNVs", "Ins", "Del")
#   for (type in c("SNVs", "MNVs", "Ins", "Del", "All")) {
#     
#     if (type == "All") {
#       tmp_df <- pcawg_lung_somatic_vcfs[[i]]
#       clonality_summary_pcawg[j,3] <- as.vector(table(tmp_df[,"CLONALITY"],  useNA="always")["TRUE"])
#       clonality_summary_pcawg[j,4] <- as.vector(table(tmp_df[,"CLONALITY"],  useNA="always")["FALSE"])
#       clonality_summary_pcawg[j,5] <- as.vector(tail(table(tmp_df[,"CLONALITY"],  useNA="always"),1))
#     } else if (type == "SNVs") {
#       tmp_df <- pcawg_lung_somatic_vcfs[[i]]
#       clonality_summary_pcawg[j+1,3] <- as.vector(table(tmp_df[nchar(tmp_df$REF) == 1 & nchar(tmp_df$ALT) == 1 ,"CLONALITY"],  useNA="always")["TRUE"])
#       clonality_summary_pcawg[j+1,4] <- as.vector(table(tmp_df[nchar(tmp_df$REF) == 1 & nchar(tmp_df$ALT) == 1 ,"CLONALITY"],  useNA="always")["FALSE"])
#       clonality_summary_pcawg[j+1,5] <- as.vector(tail(table(tmp_df[nchar(tmp_df$REF) == 1 & nchar(tmp_df$ALT) == 1 ,"CLONALITY"],  useNA="always"),1))
#     } else if (type == "MNVs") {
#       tmp_df <- pcawg_lung_somatic_vcfs[[i]]
#       clonality_summary_pcawg[j+2,3] <- as.vector(table(tmp_df[nchar(tmp_df$REF) == nchar(tmp_df$ALT) & nchar(tmp_df$REF) > 1 ,"CLONALITY"],  useNA="always")["TRUE"])
#       clonality_summary_pcawg[j+2,4] <- as.vector(table(tmp_df[nchar(tmp_df$REF) == nchar(tmp_df$ALT) & nchar(tmp_df$REF) > 1 ,"CLONALITY"],  useNA="always")["FALSE"])
#       clonality_summary_pcawg[j+2,5] <- as.vector(tail(table(tmp_df[nchar(tmp_df$REF) == nchar(tmp_df$ALT) & nchar(tmp_df$REF) > 1 ,"CLONALITY"],  useNA="always"),1))
#     } else if (type == "Ins") {
#       tmp_df <- pcawg_lung_somatic_vcfs[[i]]
#       clonality_summary_pcawg[j+3,3] <- as.vector(table(tmp_df[nchar(tmp_df$REF) < nchar(tmp_df$ALT) & nchar(tmp_df$REF) ,"CLONALITY"],  useNA="always")["TRUE"])
#       clonality_summary_pcawg[j+3,4] <- as.vector(table(tmp_df[nchar(tmp_df$REF) < nchar(tmp_df$ALT) & nchar(tmp_df$REF) ,"CLONALITY"],  useNA="always")["FALSE"])
#       clonality_summary_pcawg[j+3,5] <- as.vector(tail(table(tmp_df[nchar(tmp_df$REF) < nchar(tmp_df$ALT) & nchar(tmp_df$REF) ,"CLONALITY"],  useNA="always"),1))
#     } else if (type == "Del") {
#       tmp_df <- pcawg_lung_somatic_vcfs[[i]]
#       clonality_summary_pcawg[j+4,3] <- as.vector(table(tmp_df[nchar(tmp_df$REF) > nchar(tmp_df$ALT) & nchar(tmp_df$REF) ,"CLONALITY"],  useNA="always")["TRUE"])
#       clonality_summary_pcawg[j+4,4] <- as.vector(table(tmp_df[nchar(tmp_df$REF) > nchar(tmp_df$ALT) & nchar(tmp_df$REF) ,"CLONALITY"],  useNA="always")["FALSE"])
#       clonality_summary_pcawg[j+4,5] <- as.vector(tail(table(tmp_df[nchar(tmp_df$REF) > nchar(tmp_df$ALT) & nchar(tmp_df$REF) ,"CLONALITY"],  useNA="always"),1))
#     } 
#   }
# }
# 
# 
# if (dir.exists("/hpc/cuppen/")){
#   saveRDS(clonality_summary_pcawg, file = "/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/pcawg/clonality_summary_pcawg.rds")
# } else {
#   saveRDS(clonality_summary_pcawg, file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/pcawg/clonality_summary_pcawg.rds")
# }




if (dir.exists("/hpc/cuppen/")){
  clonality_summary_pcawg <- readRDS(file = "/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/pcawg/clonality_summary_pcawg.rds")
} else {
  clonality_summary_pcawg <- readRDS(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/pcawg/clonality_summary_pcawg.rds")
}

#========================================================================================================================================================

# @5: ploidy info of samples


# pcawg_ploidy_info <- list()
# 
# for (i in 1:nrow(pcawg_lung_meta)){
# 
#   print(i)
#   print(pcawg_lung_meta$icgc_donor_id[i])
#   
#   if (dir.exists("/hpc/cuppen/")){
#     ploidy_info <- read.csv(file = paste0("/hpc/cuppen/shared_resources/PCAWG/pipeline5/per-donor/",
#                                   pcawg_lung_meta$icgc_donor_id[i], "-from-jar/purple53/", pcawg_lung_meta$icgc_donor_id[i], "T.purple.purity.tsv"), header = T, sep = "\t", stringsAsFactors = F)
#   } else {
#     ploidy_info <- read.csv(file = paste0("/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/shared_resources/PCAWG/pipeline5/per-donor/",
#                                   pcawg_lung_meta$icgc_donor_id[i], "-from-jar/purple53/", pcawg_lung_meta$icgc_donor_id[i], "T.purple.purity.tsv"), header = T, sep = "\t", stringsAsFactors = F)
#   }
# 
#   
#   pcawg_ploidy_info[[pcawg_lung_meta$icgc_donor_id[i]]] <- ploidy_info
# }
# 
# 
# 
# if (dir.exists("/hpc/cuppen/")){
#   saveRDS(pcawg_ploidy_info, file = "/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/pcawg/pcawg-ploidy-info.rds")
# } else {
#   saveRDS(pcawg_ploidy_info, file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/pcawg/pcawg-ploidy-info.rds")
# }




if (dir.exists("/hpc/cuppen/")){
  pcawg_ploidy_info <- readRDS(file = "/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/pcawg/pcawg-ploidy-info.rds")
} else {
  pcawg_ploidy_info <- readRDS(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/pcawg/pcawg-ploidy-info.rds")
}


#========================================================================================================================================================

# clustering

head(count_df_pcawg)

count_df_pcawg_norm <- count_df_pcawg[,2:5]#/count_df_pcawg[,6]
rownames(count_df_pcawg_norm) <- count_df_pcawg$sampleId
count_df_pcawg_norm <- t(count_df_pcawg_norm)


sv_count_df_pcawg <- as.data.frame(sv_count_matrix_pcawg)
head(sv_count_df_pcawg)

sv_count_df_pcawg_norm <- sv_count_df_pcawg[,1:5]#/rowSums(sv_count_df_pcawg[,1:5])

sv_count_df_pcawg_norm <- t(sv_count_df_pcawg_norm)


mut_count_combined <- rbind(count_df_pcawg_norm, sv_count_df_pcawg_norm)
hc.sample.combined <- hclust(dist(t(mut_count_combined)), method = "complete")
sample_order_comb <- colnames(count_df_pcawg_norm)[hc.sample.combined$order]



# Mutation distribution per sample

count_tibb_pcawg <- count_df_pcawg %>%
  gather(key = "Mut_type", value = "Count", 2:5)

count_tibb_pcawg %<>% select(-all)

count_tibb_pcawg$Mut_type <- factor(count_tibb_pcawg$Mut_type, levels = rev(c("inDEL", "INdel", "MNVs", "SNVs")))

count_tibb_pcawg$sampleId <- factor(count_tibb_pcawg$sampleId, levels = sample_order_comb)

cols <- brewer.pal(8, name = "Dark2")
tmb_plot <- count_tibb_pcawg %>% ggplot(aes(x = sampleId, y = Count, fill = Mut_type, color = Mut_type)) +
  geom_bar(position="stack", stat="identity", width = 1, size = 0.1) +
  scale_fill_manual(values = cols[c(1,2,7,8)], labels = c("SNVs", "MNVs", "Small Insertions", "Small Deletions")) +
  scale_color_manual(values = rep('black', times =5), guide=F) +
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
                     breaks = c(0, 25000, 50000, 75000, 100000, 125000, 150000),
                     labels = c(0, "25,000", "50,000", "75,000", "100,000", "125,000", "150,000")) +
  labs(x = NULL, y = "Counts", fill = "Mutation Type") +
  ggtitle("Absolute Mutation Counts in PCAWG Lung Samples \n \n") +
  theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.35))
  # scale_fill_discrete(labels = c("SNVs", "MNVs", "Small Insertions", "Small Deletions"))





# SV event distributions per sample

sv_count_df_pcawg %<>% mutate(sampleId = rownames(sv_count_df_pcawg))
rownames(sv_count_df_pcawg) <- 1:nrow(sv_count_df_pcawg)

sv_count_tibb_pcawg <- sv_count_df_pcawg %>% gather(key = "Mut_type", value = "Count", 1:5)
sv_count_tibb_pcawg$Mut_type <- factor(sv_count_tibb_pcawg$Mut_type, levels = c("Deletion (< 100kb)", "Deletion (> 100kb)", "Duplication (< 100kb)", "Duplication (> 100kb)", "Inversion"))

sv_count_tibb_pcawg$sampleId <- factor(sv_count_tibb_pcawg$sampleId, levels = sample_order_comb)


sv_plot <- sv_count_tibb_pcawg %>% ggplot(aes(x = sampleId, y = Count, fill = Mut_type, color = Mut_type)) +
  geom_bar(position="stack", stat="identity", width = 1, size = 0.1) +
  scale_fill_manual(values = brewer.pal(5, "Accent")) +
  scale_color_manual(values = rep('black', times =5), guide=F) +
  theme_bw() +
  theme(
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_line(color = "black", size = 0.5, linetype = "dashed")
  ) +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.text.x = element_blank()) +
  labs(x = "PCAWG NSCLC Samples (n = 88)", y = "Counts", fill = "SV Type")


# Combine



tmb_plot_grob <- ggplotGrob(tmb_plot)
sv_plot_grob <- ggplotGrob(sv_plot)

pcawg_mut_dist <- cowplot::plot_grid(tmb_plot_grob, sv_plot_grob, align = "v", rel_heights = c(1.5,1),nrow = 2, ncol = 1)



# save

for (i in 1:2){
  if (i == 1) {
    png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/pcawg/figs/png/mut-count-per-sample-final.png")
    print(pcawg_mut_dist)
    dev.off()
  }
  if (i == 2) {
    pdf(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/pcawg/figs/pdf/mut-count-per-sample-final.pdf")
    print(pcawg_mut_dist)
    dev.off()
  }
}




# checking

all(unique(sv_count_tibb_pcawg$sampleId) == unique(count_tibb_pcawg$sampleId))


count_tibb_pcawg[count_tibb_pcawg$Count == max(count_tibb_pcawg$Count),]
count_tibb_pcawg[count_tibb_pcawg$sampleId == "DO27747",]

sv_count_tibb_pcawg[sv_count_tibb_pcawg$Count == max(sv_count_tibb_pcawg$Count),]
sv_count_tibb_pcawg[sv_count_tibb_pcawg$sampleId == "DO25189",]

######################################################################################################################################################
# SV event size distribution

del <- density(sv_length_dfs_pcawg[sv_length_dfs_pcawg$ResolvedType == "DEL" & sv_length_dfs_pcawg$length < 100000,"length"], width = 10000, from = 0)
dup <- density(sv_length_dfs_pcawg[sv_length_dfs_pcawg$ResolvedType == "DUP" & sv_length_dfs_pcawg$length < 100000,"length"], width = 10000, from = 0)
plot(del)
plot(dup)

table(sv_length_dfs_pcawg$ResolvedType)


######################################################################################################################################################
# Clonality


# head(clonality_summary_pcawg)

clonality_summary_pcawg$Tot <- rowSums(clonality_summary_pcawg[,3:5], na.rm = T)
clonality_summary_pcawg %<>% gather(key = "Clonality_stat", value = "Counts", 3:5)
clonality_summary_pcawg$Mut_type <- factor(clonality_summary_pcawg$Mut_type, levels = c("All", "SNVs", "MNVs", "Ins", "Del"))
clonality_summary_pcawg$Clonality_stat <- factor(clonality_summary_pcawg$Clonality_stat, levels = c("Clonal", "Non_clonal", "Unknown"))


options(scipen=999)



clonality_plot_summary_pcawg <- clonality_summary_pcawg %>% ggplot(aes(x = sampleId, y = Counts, fill = Clonality_stat)) + facet_grid(Mut_type ~ ., scales = 'free') +
  geom_bar(position="stack", stat="identity", width = 1) +
  scale_fill_manual(values = c("red", "black", "yellow")) +
  ggtitle("Clonality States of Different Mutation Types \n") +
<<<<<<< HEAD
  labs(x = "PCAWG NSCLC Samples (n = 88)", y = "Frequency", fill = "Clonality") +
=======
  labs(x = "PCAWG Lung Samples", y = "Frequency", fill = "Clonality") +
>>>>>>> d4b2ba7d869eab81f62853d96bcbc0ac881b0e01
  theme_bw() +
  theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.5)) +
  theme(axis.ticks.x = element_blank()) + theme(axis.text.x = element_blank()) + theme(panel.grid.major.y = element_line(color = "black", size = 0.5, linetype = "dashed"))



# Save

for (i in 1:2){
  if (i == 1) {
    png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/pcawg/figs/png/clonality-plot-summary-pcawg-final.png", width = 960)
    print(clonality_plot_summary_pcawg)
    dev.off()
  }
  if (i == 2) {
    pdf(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/analysis/pcawg/figs/pdf/clonality-plot-summary-pcawg-final.pdf", width = 14)
    print(clonality_plot_summary_pcawg)
    dev.off()
  }
}

