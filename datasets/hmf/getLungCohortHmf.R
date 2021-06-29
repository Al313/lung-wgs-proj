# ==========================================================================================================================================
# In this script I want to get the lung samples from HMF data that:
# 1. Have appropriate clinical data available
# 2. Are not secondary or tertiary biopsies of the same patient
# 3. Are not mis-annotated according to Cuplr
# ==========================================================================================================================================

library("dplyr")
library("readxl")
`%notin%` <- Negate(`%in%`)


metadata <- read.csv(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/datasets/hmf/metadata_whitelisted.tsv",
                     header =T, sep = "\t", stringsAsFactors = F)


lung_cohort <- metadata[metadata$primaryTumorLocation == "Lung",]
str(lung_cohort)


clinical_data <- read_excel("/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/datasets/hmf/clinical-data.xlsx")
clinical_data <- as.data.frame(clinical_data)
head(clinical_data)
# clinical_data <- clinical_data[,-c(7,9)]
# colnames(clinical_data)[8] <- "HMF_ID"
str(clinical_data)



# The two following records had a specila character (\n) in their "remarks" field that would cause problems. But I removed those special characters and it works fine now 
# clinical_data[clinical_data$CPCT_sampleID == "DRUP01010116T",]
# clinical_data[clinical_data$CPCT_sampleID == "CPCT02010836T",]

# According to the data collector three samples (DRUP01010010T, CPCT02020294T, CPCT02020639T) did not have optimal clinical records and should be excluded


# clinical_data[clinical_data$CPCT_sampleID == "DRUP01010010T",]
# clinical_data[clinical_data$CPCT_sampleID == "CPCT02020294T",]
# clinical_data[clinical_data$CPCT_sampleID == "CPCT02020639T",]

clinical_data <- clinical_data[clinical_data$CPCT_sampleID %notin% c("DRUP01010010T", "CPCT02020294T", "CPCT02020639T"), ]
nrow(clinical_data)


# According to Cuplr 5 samples (CPCT02011080T, CPCT02040208T, CPCT02040302T, CPCT02190039T, CPCT02040340T) were mis-annotated as lung. We should remove them

# lung_cohort[lung_cohort$sampleId == "CPCT02011080T", ]
# lung_cohort[lung_cohort$sampleId == "CPCT02040208T", ]
# lung_cohort[lung_cohort$sampleId == "CPCT02040302T", ]
# lung_cohort[lung_cohort$sampleId == "CPCT02190039T", ]
# lung_cohort[lung_cohort$sampleId == "CPCT02040340T", ]

lung_cohort <- lung_cohort[lung_cohort$sampleId %notin% c("CPCT02011080T", "CPCT02040208T", "CPCT02040302T", "CPCT02190039T", "CPCT02040340T"),]
nrow(lung_cohort)
head(lung_cohort)


sum(clinical_data$CPCT_sampleID %in% lung_cohort$sampleId)
setdiff(setdiff(clinical_data$CPCT_sampleID, lung_cohort$sampleId), c("CPCT02011080T", "CPCT02040208T", "CPCT02040302T", "CPCT02190039T", "CPCT02040340T"))

clinical_data_available <- clinical_data[clinical_data$CPCT_sampleID %in% lung_cohort$sampleId,]

colnames(clinical_data_available)[6] <- "sampleId"
colnames(clinical_data_available)


lung_metadata <- left_join(lung_cohort, clinical_data_available, by = "sampleId")
str(lung_metadata)
 
# head(lung_metadata)
# nrow(lung_metadata[!is.na(lung_metadata$CPCT_ptID),])
# lung_metadata$X.patientId == lung_metadata$CPCT_ptID

lung_metadata <- lung_metadata[!is.na(lung_metadata$CPCT_ptID),]

# all(lung_metadata$X.patientId == lung_metadata$CPCT_ptID)


table(lung_metadata$primaryTumorSubType)
table(lung_metadata$cancer_type)

# The cancer subtype between the collected clinical data and HMF metadata seemed to be discrepant for some samples. Below I remove those samples to get the final table
nrow(lung_metadata[lung_metadata$primaryTumorSubType == "Non-small cell carcinoma" & lung_metadata$cancer_type == "Non-Small Cell",])
lung_metadata[lung_metadata$primaryTumorSubType == "Non-small cell carcinoma" & lung_metadata$cancer_type == "Non-Small Cell","primaryTumorSubType"] <- "Non-small-cell-carcinoma"
nrow(lung_metadata[lung_metadata$primaryTumorSubType == "Small cell carcinoma" & lung_metadata$cancer_type == "Small Cell",])
lung_metadata[lung_metadata$primaryTumorSubType == "Small cell carcinoma" & lung_metadata$cancer_type == "Small Cell","primaryTumorSubType"] <- "Small-cell-carcinoma"


lung_metadata <- lung_metadata[lung_metadata$primaryTumorSubType == "Non-small-cell-carcinoma" | lung_metadata$primaryTumorSubType == "Small-cell-carcinoma",]
nrow(lung_metadata)

write.table(lung_metadata, file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/datasets/hmf/hmf-lung-cohort.tsv",
          row.names = F, sep = "\t", quote = F)





####### This is the final list of HMF samples that were used for this analysis by Job. It was sent to me in an email on June 25th!

job_data <- readxl::read_xlsx(path = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/datasets/hmf/DR142_Samples.xlsx")
job_sampleIds <- job_data$sampleId

# WIDE01010708T is a special case (it is annotated as a cancer with unknown primary site in hmf metadata but in the clinical data it is lung cancer)
lung_cohort <- metadata[metadata$primaryTumorLocation == "Lung" | metadata$sampleId == 'WIDE01010708T',]

lung_cohort <- lung_cohort[lung_cohort$sampleId %in% job_data$sampleId,]
str(lung_cohort)

# for "CPCT02060086T" "CPCT02060261T" "CPCT02060265T" I don't have updated clinical data so I just try to add them myself
t <- job_sampleIds[job_sampleIds %notin% clinical_data$CPCT_sampleID]

job_data[job_data$sampleId %in% t,]

head(clinical_data)
str(clinical_data)
clinical_data[415:417, "CPCT_sampleID"] <- t
clinical_data[415:417, "pat_sex"] <- as.vector(job_data[job_data$sampleId %in% t,"gender"])
clinical_data[415:417, "cancer_type"] <- as.vector(job_data[job_data$sampleId %in% t,"cancer_type"])

tail(clinical_data)


colnames(clinical_data)[6] <- "sampleId"


clinical_data <- clinical_data[clinical_data$sampleId %in% job_data$sampleId,]

str(clinical_data)

lung_metadata <- left_join(lung_cohort, clinical_data, by = "sampleId")
table(lung_metadata$cancer_type)
table(job_data$cancer_type)

lung_metadata[lung_metadata$cancer_type == "Mix of SCLC & NSCLC","cancer_type"] <- "Small Cell"
lung_metadata[lung_metadata$cancer_type == "Small-Cell","cancer_type"] <- "Small Cell"

write.table(lung_metadata, file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/datasets/hmf/hmf-lung-cohort-final.tsv",
            row.names = F, sep = "\t", quote = F)


######
# Checking with Job's list



lung_metadata <- read.csv(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-lung-proj/datasets/hmf/hmf-lung-cohort.tsv",
                          sep = "\t", header = T, stringsAsFactors = F)
str(lung_metadata)
nrow(lung_metadata[lung_metadata$cancer_type == "Non-Small Cell",])
nrow(lung_metadata[lung_metadata$cancer_type == "Small Cell",])

my_list <- lung_metadata$sampleId[lung_metadata$cancer_type == "Non-Small Cell"]
job_lung_metadata <- read_excel("/home/ali313/Desktop/includedSamples_DR142.xlsx")

str(job_lung_metadata)
nrow(job_lung_metadata[job_lung_metadata$cancer_type == "Non-Small Cell",])
nrow(job_lung_metadata[job_lung_metadata$cancer_type == "Small Cell",])
job_lung_metadata[job_lung_metadata$sampleId == "WIDE01010708T",]
job_list <- job_lung_metadata$sampleId[job_lung_metadata$cancer_type == "Non-Small Cell"]

length(job_list)

setdiff(job_list, my_list)
length(setdiff(my_list, job_list))
setdiff(my_list, job_list)


lung_cohort[lung_cohort$sampleId %in% setdiff(job_list, my_list), ]
metadata[metadata$sampleId == "WIDE01010708T", "primaryTumorLocation"]
lung_cohort[lung_cohort$sampleId == "CPCT02020294T",]
clinical_data[clinical_data$CPCT_sampleID == "WIDE01010708T",]
clinical_data[clinical_data$CPCT_sampleID == "CPCT02020294T",]

job_lung_metadata[job_lung_metadata$sampleId == "WIDE01010708T",]

lung_cohort[lung_cohort$sampleId %in% setdiff(my_list, job_list), ]
metadata[metadata$sampleId %in% setdiff(my_list, job_list),]
metadata[metadata$X.patientId == "CPCT02190020",]

clinical_data[clinical_data$CPCT_sampleID %in% setdiff(my_list, job_list),]
clinical_data[clinical_data$CPCT_sampleID %in% setdiff( job_list, my_list),]
