# Lineage tracing of BCRs

# Merge all the BCR files
# Merge the files with metadata from Seurat containing patient vs healthy information
# Split the table into healthy and patient tables
# Reformat to look like filtered_contig_annotations.csv file


# Load libraries
library(Seurat)
library(tidyverse)
library(RCurl)
library(cowplot)
library(future)
library(stringr.tools)
library(Biostrings)


# Parallelize
plan("multicore")

# Change max loading size to 40 GB
options(future.globals.maxSize = 40000 * 1024^2)

# Set working directory
setwd("/home/suk571/Jocelyn_scRNA/KW9079_Farmer/210624_A00873_0376_BHC3KTDSX2_KW9079/cellranger-6.0.1/GRCh38/")

# Load initial seurat file to get metadata
load(file = "results/output_rds/merged_seurat.RData")

# Store metadata in a table
meta <- merged_seurat@meta.data


## Get the B cell barcodes
## Add well number prefixes to match the Seurat object
B_Cell_BC <- list.files("B_VDJ/")
all_clones <- NULL
i <- 1
for (file in B_Cell_BC){
  readFile <- read.csv(file = paste0("B_VDJ/", file), header = TRUE)
  readFile$barcode <- str_prefix(readFile$barcode, paste0("well_", i, "_"))
  all_clones <- rbind(all_clones, readFile)
  i <- i + 1
}
rm(readFile, B_Cell_BC, file, i)

# Use data.table package for efficient calculations. Dplyr takes a long time in my experience.
all_clones <- data.table::as.data.table(all_clones)

# Rename cells column in metadata to match BCR clones table - for merging
meta <- meta %>% dplyr::rename(barcode = cells)

# Merge metadata to BCR clones
merged_bcr <- merge(all_clones, meta)

# Retain only well, sample and patient ID information in BCR clones
merged_bcr <- merged_bcr %>% select(-c(seq_folder, nUMI, nGene, nCount_ADT, nFeature_ADT, log10GenesPerUMI, mitoRatio))

# Split BCR clones table by sample
healthy_BCR <- subset(merged_bcr, sample == "healthy")
patient_BCR <- subset(merged_bcr, sample == "patient")

# Remove sample column
healthy_BCR <- healthy_BCR %>% select(-c(sample))
patient_BCR <- patient_BCR %>% select(-c(sample))

# Remove well information from the barcodes

healthy_BCR$barcode <- str_remove(healthy_BCR$barcode, "well_1_")
healthy_BCR$barcode <- str_remove(healthy_BCR$barcode, "well_2_")
healthy_BCR$barcode <- str_remove(healthy_BCR$barcode, "well_3_")
healthy_BCR$barcode <- str_remove(healthy_BCR$barcode, "well_4_")

patient_BCR$barcode <- str_remove(patient_BCR$barcode, "well_1_")
patient_BCR$barcode <- str_remove(patient_BCR$barcode, "well_2_")
patient_BCR$barcode <- str_remove(patient_BCR$barcode, "well_3_")
patient_BCR$barcode <- str_remove(patient_BCR$barcode, "well_4_")


# Split BCR clones by origin (patient) - 4 technical replicates mixed
# Healthy - hash 1, 3, 5, 6
# Patient - hash 2, 4, 7, 8
healthy_1 <- subset(healthy_BCR, patient_ID == "Hashtag1")
healthy_3 <- subset(healthy_BCR, patient_ID == "Hashtag3")
healthy_5 <- subset(healthy_BCR, patient_ID == "Hashtag5")
healthy_6 <- subset(healthy_BCR, patient_ID == "Hashtag6")

patient_2 <- subset(patient_BCR, patient_ID == "Hashtag2")
patient_4 <- subset(patient_BCR, patient_ID == "Hashtag4")
patient_7 <- subset(patient_BCR, patient_ID == "Hashtag7")
patient_8 <- subset(patient_BCR, patient_ID == "Hashtag8")


# Remove patient ID column

healthy_1 <- healthy_1 %>% select(-c(patient_ID))
healthy_3 <- healthy_3 %>% select(-c(patient_ID))
healthy_5 <- healthy_5 %>% select(-c(patient_ID))
healthy_6 <- healthy_6 %>% select(-c(patient_ID))

patient_2 <- patient_2 %>% select(-c(patient_ID))
patient_4 <- patient_4 %>% select(-c(patient_ID))
patient_7 <- patient_7 %>% select(-c(patient_ID))
patient_8 <- patient_8 %>% select(-c(patient_ID))


# Save the files
write.csv(healthy_1, file = "B_VDJ/healthy/1_filtered_contig_annotations.csv", row.names = FALSE)
write.csv(healthy_3, file = "B_VDJ/healthy/3_filtered_contig_annotations.csv", row.names = FALSE)
write.csv(healthy_5, file = "B_VDJ/healthy/5_filtered_contig_annotations.csv", row.names = FALSE)
write.csv(healthy_6, file = "B_VDJ/healthy/6_filtered_contig_annotations.csv", row.names = FALSE)

write.csv(patient_2, file = "B_VDJ/patient/2_filtered_contig_annotations.csv", row.names = FALSE)
write.csv(patient_4, file = "B_VDJ/patient/4_filtered_contig_annotations.csv", row.names = FALSE)
write.csv(patient_7, file = "B_VDJ/patient/7_filtered_contig_annotations.csv", row.names = FALSE)
write.csv(patient_8, file = "B_VDJ/patient/8_filtered_contig_annotations.csv", row.names = FALSE)



