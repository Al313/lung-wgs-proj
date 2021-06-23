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
