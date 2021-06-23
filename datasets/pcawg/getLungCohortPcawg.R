# reading in the required libraries and functions

library("dplyr")
library("readxl")
library("stringr")
`%notin%` <- Negate(`%in%`)


# Let's explore the PCAWG lung cohort and its clinical data a bit

metadata <- read_excel("/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/datasets/pcawg/_pcawg_specimen_histology_August2016_v9.xlsx")
str(metadata)
metadata <- as.data.frame(metadata)

duplicated(metadata$icgc_donor_id)
any(duplicated(metadata$icgc_sample_id))

metadata[duplicated(metadata$icgc_donor_id),"icgc_donor_id"]
metadata[metadata$icgc_donor_id == "DO10172",]




metadata_wgs <- metadata[metadata$specimen_library_strategy %in% c("WGS", "WGS+RNA-Seq"),]


metadata_wgs[duplicated(metadata_wgs$icgc_donor_id),"icgc_donor_id"]
metadata_wgs[metadata_wgs$icgc_donor_id == "DO23543",]

nrow(metadata_wgs[metadata_wgs$organ_system == "LUNG & BRONCHUS",])





clinical_data <- read_excel("/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/datasets/pcawg/_pcawg_donor_clinical_August2016_v9.xlsx")
str(clinical_data)
clinical_data <- as.data.frame(clinical_data)
head(clinical_data)

any(duplicated(clinical_data$icgc_donor_id))

table(clinical_data$tobacco_smoking_history_indicator)




# diplotypes are only available for the ICGC cohort of the PCAWG dataset

# diplotypes <- read.csv(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/annotation-beds/diplotypes_germPonFilt.txt.gz",
#                        sep = "\t", header = T, stringsAsFactors = F)
# 
# 
# pcawg_diplotypes <- diplotypes[str_detect(diplotypes$sample, "^DO", negate = F),]
# length(unique(pcawg_diplotypes$sample))
# 




# This file was on HPC in the Metdata repo

donor <- read.csv(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/datasets/pcawg/donor.tsv",
                  sep = "\t", header = T, stringsAsFactors = F)


# I get the list of HMF-processed PCAWG samples from the "per-donor" repo

per_donor_pcawg <- scan(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/datasets/pcawg/per-donor-pcawg.txt", what = character())

per_donor_pcawg <- sapply(strsplit(per_donor_pcawg, split = "-", ), `[`, 1)
length(per_donor_pcawg)

setdiff(per_donor_pcawg,donor$icgc_donor_id)
# one sample "DO46468" does not have an entry in donor metadata 

setdiff(donor$icgc_donor_id, per_donor_pcawg)
# 142 samples don't have a directory in the per-donor directory


donor_available <- donor[donor$icgc_donor_id %in% per_donor_pcawg,]
nrow(donor_available) # all the samples are present except "DO46468"
head(donor_available)

sum(is.na(donor_available$icgc_donor_id))

for (i in donor_available$icgc_donor_id[1:10]){
  sample_id <- metadata[metadata$icgc_donor_id == i,"icgc_donor_id"]
  if (length(sample_id) > 1){
    print(i)
  }
}

metadata[metadata$icgc_donor_id == "DO23018",]


setdiff(donor_available$icgc_donor_id, metadata_wgs$icgc_donor_id)

cols_from_meta <- data.frame(icgc_donor_id = NA, specimen_donor_treatment_type = NA, histology_tier3 = NA, organ_system = NA, donor_wgs_included_excluded = NA)
colnames(metadata_wgs)
j <- 0
for (i in 1:length(donor_available$icgc_donor_id)){
  # print(donor_available$icgc_donor_id[i])
  sample_id <- metadata_wgs[metadata_wgs$icgc_donor_id == donor_available$icgc_donor_id[i],"icgc_donor_id"]
  
  if (length(sample_id) > 1) {
    j <- j +1
    cols_from_meta[j,1] <- donor_available$icgc_donor_id[i]
    cols_from_meta[j,2:5] <- metadata_wgs[metadata_wgs$icgc_donor_id == donor_available$icgc_donor_id[i],c(26, 16, 12, 27)][1,]
  } else if (length(sample_id) == 1){
    j <- j +1
    cols_from_meta[j,1] <- donor_available$icgc_donor_id[i]
    cols_from_meta[j,2:5] <- metadata_wgs[metadata_wgs$icgc_donor_id == donor_available$icgc_donor_id[i],c(26, 16, 12, 27)]
  }
}

donor_available1 <- left_join(donor_available, cols_from_meta, by = "icgc_donor_id")



setdiff(donor_available$icgc_donor_id, clinical_data$icgc_donor_id)

cols_from_clinical <- data.frame(icgc_donor_id = NA, tobacco_smoking_history_indicator = NA, tobacco_smoking_intensity = NA)

for (i in 1:length(donor_available$icgc_donor_id)){
  # print(donor_available$icgc_donor_id[i])
  sample_id <- clinical_data[clinical_data$icgc_donor_id == donor_available$icgc_donor_id[i],"icgc_donor_id"]
  
  cols_from_clinical[i,1] <- donor_available$icgc_donor_id[i]
  cols_from_clinical[i,2:3] <- clinical_data[clinical_data$icgc_donor_id == donor_available$icgc_donor_id[i],c(14, 15)]

}

donor_available2 <- left_join(donor_available1, cols_from_clinical, by = "icgc_donor_id")
colnames(donor_available2)







saveRDS(donor_available2, file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/datasets/pcawg/metadata-per-donor.Rds")


nrow(metadata_wgs[metadata_wgs$organ_system == "LUNG & BRONCHUS" & metadata_wgs$donor_wgs_included_excluded == "Included",])
nrow(donor_available2[donor_available2$organ_system == "LUNG & BRONCHUS" & donor_available2$donor_wgs_included_excluded == "Included",])


table(metadata_wgs$organ_system)
table(donor_available2$organ_system)

table(metadata_wgs$project_code[metadata_wgs$organ_system == "LUNG & BRONCHUS"])
table(donor_available2$project_code[donor_available2$organ_system == "LUNG & BRONCHUS"])


# Check the integrity of data
included_lung_icgc <- donor_available2$icgc_donor_id[donor_available2$organ_system == "LUNG & BRONCHUS" & donor_available2$donor_wgs_included_excluded == "Included"]
clinical_data[clinical_data$icgc_donor_id %in% included_lung_icgc,]
table(metadata_wgs[metadata_wgs$icgc_donor_id %in% included_lung_icgc,"histology_tier3"])

