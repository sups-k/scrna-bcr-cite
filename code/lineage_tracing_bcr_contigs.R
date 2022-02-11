# Data cleanup to perform lineage tracing of BCRs

# Load libraries
library(Seurat)
library(tidyverse)
library(RCurl)
library(cowplot)
library(future)
library(stringr.tools)
library(Biostrings)

# Set working directory
setwd("/home/suk571/Jocelyn_scRNA/KW9079_Farmer/210624_A00873_0376_BHC3KTDSX2_KW9079/cellranger-6.0.1/GRCh38/")



########################################################################
############## READ FASTA FILES AS TABLES FOR EACH WELL ################
########################################################################

# Read FASTA file of well 1 - BRI-1244
fasta_well_1 <- readDNAStringSet("BRI-1244/outs/filtered_contig.fasta")

# Convert to data frame
contig_id <- names(fasta_well_1)
sequence_well <- paste(fasta_well_1)
df_well_1 <- data.frame(contig_id, sequence_well)

# Remove extra variables
rm(fasta_well_1, contig_id, sequence_well)




# Read FASTA file of well 2 - BRI-1247
fasta_well_2 <- readDNAStringSet("BRI-1247/outs/filtered_contig.fasta")
contig_id <- names(fasta_well_2)
sequence_well <- paste(fasta_well_2)
df_well_2 <- data.frame(contig_id, sequence_well)
rm(fasta_well_2, contig_id, sequence_well)

# Read FASTA file of well 3 - BRI-1250
fasta_well_3 <- readDNAStringSet("BRI-1250/outs/filtered_contig.fasta")
contig_id <- names(fasta_well_3)
sequence_well <- paste(fasta_well_3)
df_well_3 <- data.frame(contig_id, sequence_well)
rm(fasta_well_3, contig_id, sequence_well)

# Read FASTA file of well 4 - BRI-1253
fasta_well_4 <- readDNAStringSet("BRI-1253/outs/filtered_contig.fasta")
contig_id <- names(fasta_well_4)
sequence_well <- paste(fasta_well_4)
df_well_4 <- data.frame(contig_id, sequence_well)
rm(fasta_well_4, contig_id, sequence_well)




######################################################################################
############## COMBINE FASTA FILES WITH BCR ANNOTATIONS FOR EACH WELL ################
######################################################################################

## Get the B cell barcodes
bcr_well_1 <- read.csv(file = "B_VDJ/BRI-1242_filtered_contig_annotations.csv", header = TRUE)

# Combine with FASTA sequence
bcr_seq_well_1 <- merge(bcr_well_1, df_well_1)
bcr_seq_well_1$barcode <- str_prefix(bcr_seq_well_1$barcode, "well_1_")
bcr_seq_well_1$contig_id <- str_prefix(bcr_seq_well_1$contig_id, "well_1_") # barcode and contig ID should have the same name

# Remove extra variables
rm(df_well_1, bcr_well_1)



bcr_well_2 <- read.csv(file = "B_VDJ/BRI-1245_filtered_contig_annotations.csv", header = TRUE)
bcr_seq_well_2 <- merge(bcr_well_2, df_well_2)
bcr_seq_well_2$barcode <- str_prefix(bcr_seq_well_2$barcode, "well_2_")
bcr_seq_well_2$contig_id <- str_prefix(bcr_seq_well_2$contig_id, "well_2_")
rm(df_well_2, bcr_well_2)

bcr_well_3 <- read.csv(file = "B_VDJ/BRI-1248_filtered_contig_annotations.csv", header = TRUE)
bcr_seq_well_3 <- merge(bcr_well_3, df_well_3)
bcr_seq_well_3$barcode <- str_prefix(bcr_seq_well_3$barcode, "well_3_")
bcr_seq_well_3$contig_id <- str_prefix(bcr_seq_well_3$contig_id, "well_3_")
rm(df_well_3, bcr_well_3)

bcr_well_4 <- read.csv(file = "B_VDJ/BRI-1251_filtered_contig_annotations.csv", header = TRUE)
bcr_seq_well_4 <- merge(bcr_well_4, df_well_4)
bcr_seq_well_4$barcode <- str_prefix(bcr_seq_well_4$barcode, "well_4_")
bcr_seq_well_4$contig_id <- str_prefix(bcr_seq_well_4$contig_id, "well_4_")
rm(df_well_4, bcr_well_4)






###############################################################
############## READ EACH WELL AS SEURAT OBJECT ################
###############################################################


## Create each individual Seurat object for every sample
for (file in c("BRI-1242/outs/raw_feature_bc_matrix", "BRI-1245/outs/raw_feature_bc_matrix", "BRI-1248/outs/raw_feature_bc_matrix", "BRI-1251/outs/raw_feature_bc_matrix")){
  seurat_data <- Read10X(data.dir = paste0("/home/suk571/Jocelyn_scRNA/KW9079_Farmer/210624_A00873_0376_BHC3KTDSX2_KW9079/cellranger-6.0.1/GRCh38/", file))
  rownames(x = seurat_data[["Antibody Capture"]]) <- gsub(pattern = "_[control_]*TotalSeqB", replacement = "",
                                                          x = rownames(x = seurat_data[["Antibody Capture"]]))
  seurat_obj <- CreateSeuratObject(counts = seurat_data[["Gene Expression"]], 
                                   min.features = 100, 
                                   project = file)
  seurat_obj[["ADT"]] <- CreateAssayObject(seurat_data[["Antibody Capture"]][, colnames(x = seurat_obj)])
  assign(str_split_fixed(file, "/", 2)[1], seurat_obj)
}

# Rename variables
BRI1242 <- `BRI-1242`
BRI1245 <- `BRI-1245`
BRI1248 <- `BRI-1248`
BRI1251 <- `BRI-1251`

## Remove unnecessary variables
rm(`BRI-1242`, `BRI-1245`, `BRI-1248`, `BRI-1251`, seurat_data, seurat_obj, file)

gc()


################################################################################
############## OBTAIN METADATA FOR EACH WELL - ADD PATIENT DATA ################
######################### AND COMBINE WITH BCR TABLES ##########################
################################################################################



##### METADATA OF FIRST WELL #######

## Copy the metadata into a separate data frame to avoid messing up the original
BRI1242_metadata <- BRI1242@meta.data

## Add sample and patient ID columns to the metadata
BRI1242_metadata$sample <- NA
BRI1242_metadata$patient_ID <- NA

## TO FIND OUT PATIENT ID:
# Create a new data frame from the counts matrix of ADT
# Take a column-wise max number - which row did max occur? That row = patient ID
# But exclude 1st 4 rows because they are not patient IDs
BRI1242_subset_counts <- as.data.frame(BRI1242@assays[["ADT"]]@counts)
BRI1242_subset_counts <- BRI1242_subset_counts[-c(1,2,3,4),]

# Subset counts column names == metadata row names. Hence safe to add values sequentially

# Hashtag1 = healthy
# Hashtag3 = healthy
# Hashtag5 = healthy
# Hashtag6 = healthy

# Hashtag2 = patient
# Hashtag4 = patient
# Hashtag7 = patient
# Hashtag8 = patient

## Add patient IDs to the metadata
for (i in 1:ncol(BRI1242_subset_counts)) {
  BRI1242_metadata$patient_ID[i] <- rownames(BRI1242_subset_counts)[which(BRI1242_subset_counts[,i] == max(BRI1242_subset_counts[,i]), arr.ind = TRUE)]
}

## Verify there are no NAs left in the patient IDs
BRI1242_metadata %>% summarise_all(~ sum(is.na(.)))

## Add sample info
BRI1242_metadata$sample[which(str_detect(BRI1242_metadata$patient_ID, "Hashtag2"))] <- "patient"
BRI1242_metadata$sample[which(str_detect(BRI1242_metadata$patient_ID, "Hashtag4"))] <- "patient"
BRI1242_metadata$sample[which(str_detect(BRI1242_metadata$patient_ID, "Hashtag7"))] <- "patient"
BRI1242_metadata$sample[which(str_detect(BRI1242_metadata$patient_ID, "Hashtag8"))] <- "patient"

BRI1242_metadata$sample[which(str_detect(BRI1242_metadata$patient_ID, "Hashtag1"))] <- "healthy"
BRI1242_metadata$sample[which(str_detect(BRI1242_metadata$patient_ID, "Hashtag3"))] <- "healthy"
BRI1242_metadata$sample[which(str_detect(BRI1242_metadata$patient_ID, "Hashtag5"))] <- "healthy"
BRI1242_metadata$sample[which(str_detect(BRI1242_metadata$patient_ID, "Hashtag6"))] <- "healthy"

# Add rownames as barcode column
BRI1242_metadata$barcode <- rownames(BRI1242_metadata)
BRI1242_metadata$barcode <- str_prefix(BRI1242_metadata$barcode, "well_1_")

# Remove extra variables
rm(BRI1242_subset_counts, i, BRI1242)


# Combine with BCR data
well_1 <- merge(BRI1242_metadata, bcr_seq_well_1)
rm(BRI1242_metadata, bcr_seq_well_1)
gc()




##### METADATA OF SECOND WELL #######

## Copy the metadata into a separate data frame to avoid messing up the original
BRI1245_metadata <- BRI1245@meta.data

## Add sample and patient ID columns to the metadata
BRI1245_metadata$sample <- NA
BRI1245_metadata$patient_ID <- NA

## TO FIND OUT PATIENT ID:
BRI1245_subset_counts <- as.data.frame(BRI1245@assays[["ADT"]]@counts)
BRI1245_subset_counts <- BRI1245_subset_counts[-c(1,2,3,4),]

## Add patient IDs to the metadata
for (i in 1:ncol(BRI1245_subset_counts)) {
  BRI1245_metadata$patient_ID[i] <- rownames(BRI1245_subset_counts)[which(BRI1245_subset_counts[,i] == max(BRI1245_subset_counts[,i]), arr.ind = TRUE)]
}

## Verify there are no NAs left in the patient IDs
BRI1245_metadata %>% summarise_all(~ sum(is.na(.)))

## Add sample info
BRI1245_metadata$sample[which(str_detect(BRI1245_metadata$patient_ID, "Hashtag2"))] <- "patient"
BRI1245_metadata$sample[which(str_detect(BRI1245_metadata$patient_ID, "Hashtag4"))] <- "patient"
BRI1245_metadata$sample[which(str_detect(BRI1245_metadata$patient_ID, "Hashtag7"))] <- "patient"
BRI1245_metadata$sample[which(str_detect(BRI1245_metadata$patient_ID, "Hashtag8"))] <- "patient"

BRI1245_metadata$sample[which(str_detect(BRI1245_metadata$patient_ID, "Hashtag1"))] <- "healthy"
BRI1245_metadata$sample[which(str_detect(BRI1245_metadata$patient_ID, "Hashtag3"))] <- "healthy"
BRI1245_metadata$sample[which(str_detect(BRI1245_metadata$patient_ID, "Hashtag5"))] <- "healthy"
BRI1245_metadata$sample[which(str_detect(BRI1245_metadata$patient_ID, "Hashtag6"))] <- "healthy"

# Add rownames as barcode column
BRI1245_metadata$barcode <- rownames(BRI1245_metadata)
BRI1245_metadata$barcode <- str_prefix(BRI1245_metadata$barcode, "well_2_")

# Remove extra variables
rm(BRI1245_subset_counts, i, BRI1245)

# Combine with BCR data
well_2 <- merge(BRI1245_metadata, bcr_seq_well_2)
rm(BRI1245_metadata, bcr_seq_well_2)
gc()





##### METADATA OF THIRD WELL #######

## Copy the metadata into a separate data frame to avoid messing up the original
BRI1248_metadata <- BRI1248@meta.data

## Add sample and patient ID columns to the metadata
BRI1248_metadata$sample <- NA
BRI1248_metadata$patient_ID <- NA

## TO FIND OUT PATIENT ID:
BRI1248_subset_counts <- as.data.frame(BRI1248@assays[["ADT"]]@counts)
BRI1248_subset_counts <- BRI1248_subset_counts[-c(1,2,3,4),]

## Add patient IDs to the metadata
for (i in 1:ncol(BRI1248_subset_counts)) {
  BRI1248_metadata$patient_ID[i] <- rownames(BRI1248_subset_counts)[which(BRI1248_subset_counts[,i] == max(BRI1248_subset_counts[,i]), arr.ind = TRUE)]
}

## Verify there are no NAs left in the patient IDs
BRI1248_metadata %>% summarise_all(~ sum(is.na(.)))

## Add sample info
BRI1248_metadata$sample[which(str_detect(BRI1248_metadata$patient_ID, "Hashtag2"))] <- "patient"
BRI1248_metadata$sample[which(str_detect(BRI1248_metadata$patient_ID, "Hashtag4"))] <- "patient"
BRI1248_metadata$sample[which(str_detect(BRI1248_metadata$patient_ID, "Hashtag7"))] <- "patient"
BRI1248_metadata$sample[which(str_detect(BRI1248_metadata$patient_ID, "Hashtag8"))] <- "patient"

BRI1248_metadata$sample[which(str_detect(BRI1248_metadata$patient_ID, "Hashtag1"))] <- "healthy"
BRI1248_metadata$sample[which(str_detect(BRI1248_metadata$patient_ID, "Hashtag3"))] <- "healthy"
BRI1248_metadata$sample[which(str_detect(BRI1248_metadata$patient_ID, "Hashtag5"))] <- "healthy"
BRI1248_metadata$sample[which(str_detect(BRI1248_metadata$patient_ID, "Hashtag6"))] <- "healthy"

# Add rownames as barcode column
BRI1248_metadata$barcode <- rownames(BRI1248_metadata)
BRI1248_metadata$barcode <- str_prefix(BRI1248_metadata$barcode, "well_3_")

# Remove extra variables
rm(BRI1248_subset_counts, i, BRI1248)

# Combine with BCR data
well_3 <- merge(BRI1248_metadata, bcr_seq_well_3)
rm(BRI1248_metadata, bcr_seq_well_3)
gc()





##### METADATA OF FOURTH (and final) WELL #######

## Copy the metadata into a separate data frame to avoid messing up the original
BRI1251_metadata <- BRI1251@meta.data

## Add sample and patient ID columns to the metadata
BRI1251_metadata$sample <- NA
BRI1251_metadata$patient_ID <- NA

## TO FIND OUT PATIENT ID:
BRI1251_subset_counts <- as.data.frame(BRI1251@assays[["ADT"]]@counts)
BRI1251_subset_counts <- BRI1251_subset_counts[-c(1,2,3,4),]

## Add patient IDs to the metadata
for (i in 1:ncol(BRI1251_subset_counts)) {
  BRI1251_metadata$patient_ID[i] <- rownames(BRI1251_subset_counts)[which(BRI1251_subset_counts[,i] == max(BRI1251_subset_counts[,i]), arr.ind = TRUE)]
}

## Verify there are no NAs left in the patient IDs
BRI1251_metadata %>% summarise_all(~ sum(is.na(.)))

## Add sample info
BRI1251_metadata$sample[which(str_detect(BRI1251_metadata$patient_ID, "Hashtag2"))] <- "patient"
BRI1251_metadata$sample[which(str_detect(BRI1251_metadata$patient_ID, "Hashtag4"))] <- "patient"
BRI1251_metadata$sample[which(str_detect(BRI1251_metadata$patient_ID, "Hashtag7"))] <- "patient"
BRI1251_metadata$sample[which(str_detect(BRI1251_metadata$patient_ID, "Hashtag8"))] <- "patient"

BRI1251_metadata$sample[which(str_detect(BRI1251_metadata$patient_ID, "Hashtag1"))] <- "healthy"
BRI1251_metadata$sample[which(str_detect(BRI1251_metadata$patient_ID, "Hashtag3"))] <- "healthy"
BRI1251_metadata$sample[which(str_detect(BRI1251_metadata$patient_ID, "Hashtag5"))] <- "healthy"
BRI1251_metadata$sample[which(str_detect(BRI1251_metadata$patient_ID, "Hashtag6"))] <- "healthy"

# Add rownames as barcode column
BRI1251_metadata$barcode <- rownames(BRI1251_metadata)
BRI1251_metadata$barcode <- str_prefix(BRI1251_metadata$barcode, "well_4_")

# Remove extra variables
rm(BRI1251_subset_counts, i, BRI1251)

# Combine with BCR data
well_4 <- merge(BRI1251_metadata, bcr_seq_well_4)
rm(BRI1251_metadata, bcr_seq_well_4)
gc()



###################################################################
############## GENERATE FASTA SEQUENCES BY PATIENT ################
###################################################################


###### FOR WELL 1 #########

# Make a copy of the data in well_1
fa_1 <- well_1

# Select only sample, patient ID, contig ID, sequence columns
fa_1 <- fa_1 %>% select(c(sample, patient_ID, contig_id, sequence_well))

# Separate each individual person into a different table
fa_1_healthy_1 <- subset(fa_1, patient_ID == "Hashtag1")
fa_1_healthy_3 <- subset(fa_1, patient_ID == "Hashtag3")
fa_1_healthy_5 <- subset(fa_1, patient_ID == "Hashtag5")
fa_1_healthy_6 <- subset(fa_1, patient_ID == "Hashtag6")

fa_1_patient_2 <- subset(fa_1, patient_ID == "Hashtag2")
fa_1_patient_4 <- subset(fa_1, patient_ID == "Hashtag4")
fa_1_patient_7 <- subset(fa_1, patient_ID == "Hashtag7")
fa_1_patient_8 <- subset(fa_1, patient_ID == "Hashtag8")

# Select only contig ID, sequence columns
fa_1_healthy_1 <- fa_1_healthy_1 %>% select(c(contig_id, sequence_well))
fa_1_healthy_3 <- fa_1_healthy_3 %>% select(c(contig_id, sequence_well))
fa_1_healthy_5 <- fa_1_healthy_5 %>% select(c(contig_id, sequence_well))
fa_1_healthy_6 <- fa_1_healthy_6 %>% select(c(contig_id, sequence_well))

fa_1_patient_2 <- fa_1_patient_2 %>% select(c(contig_id, sequence_well))
fa_1_patient_4 <- fa_1_patient_4 %>% select(c(contig_id, sequence_well))
fa_1_patient_7 <- fa_1_patient_7 %>% select(c(contig_id, sequence_well))
fa_1_patient_8 <- fa_1_patient_8 %>% select(c(contig_id, sequence_well))

# Remove fa_1
rm(fa_1)

###### FOR WELL 2 #########

# Make a copy of the data in well_2
fa_2 <- well_2

# Select only sample, patient ID, contig ID, sequence columns
fa_2 <- fa_2 %>% select(c(sample, patient_ID, contig_id, sequence_well))

# Separate each individual person into a different table
fa_2_healthy_1 <- subset(fa_2, patient_ID == "Hashtag1")
fa_2_healthy_3 <- subset(fa_2, patient_ID == "Hashtag3")
fa_2_healthy_5 <- subset(fa_2, patient_ID == "Hashtag5")
fa_2_healthy_6 <- subset(fa_2, patient_ID == "Hashtag6")

fa_2_patient_2 <- subset(fa_2, patient_ID == "Hashtag2")
fa_2_patient_4 <- subset(fa_2, patient_ID == "Hashtag4")
fa_2_patient_7 <- subset(fa_2, patient_ID == "Hashtag7")
fa_2_patient_8 <- subset(fa_2, patient_ID == "Hashtag8")

# Select only contig ID, sequence columns
fa_2_healthy_1 <- fa_2_healthy_1 %>% select(c(contig_id, sequence_well))
fa_2_healthy_3 <- fa_2_healthy_3 %>% select(c(contig_id, sequence_well))
fa_2_healthy_5 <- fa_2_healthy_5 %>% select(c(contig_id, sequence_well))
fa_2_healthy_6 <- fa_2_healthy_6 %>% select(c(contig_id, sequence_well))

fa_2_patient_2 <- fa_2_patient_2 %>% select(c(contig_id, sequence_well))
fa_2_patient_4 <- fa_2_patient_4 %>% select(c(contig_id, sequence_well))
fa_2_patient_7 <- fa_2_patient_7 %>% select(c(contig_id, sequence_well))
fa_2_patient_8 <- fa_2_patient_8 %>% select(c(contig_id, sequence_well))

# Remove fa_2
rm(fa_2)


###### FOR WELL 3 #########

# Make a copy of the data in well_3
fa_3 <- well_3

# Select only sample, patient ID, contig ID, sequence columns
fa_3 <- fa_3 %>% select(c(sample, patient_ID, contig_id, sequence_well))

# Separate each individual person into a different table
fa_3_healthy_1 <- subset(fa_3, patient_ID == "Hashtag1")
fa_3_healthy_3 <- subset(fa_3, patient_ID == "Hashtag3")
fa_3_healthy_5 <- subset(fa_3, patient_ID == "Hashtag5")
fa_3_healthy_6 <- subset(fa_3, patient_ID == "Hashtag6")

fa_3_patient_2 <- subset(fa_3, patient_ID == "Hashtag2")
fa_3_patient_4 <- subset(fa_3, patient_ID == "Hashtag4")
fa_3_patient_7 <- subset(fa_3, patient_ID == "Hashtag7")
fa_3_patient_8 <- subset(fa_3, patient_ID == "Hashtag8")

# Select only contig ID, sequence columns
fa_3_healthy_1 <- fa_3_healthy_1 %>% select(c(contig_id, sequence_well))
fa_3_healthy_3 <- fa_3_healthy_3 %>% select(c(contig_id, sequence_well))
fa_3_healthy_5 <- fa_3_healthy_5 %>% select(c(contig_id, sequence_well))
fa_3_healthy_6 <- fa_3_healthy_6 %>% select(c(contig_id, sequence_well))

fa_3_patient_2 <- fa_3_patient_2 %>% select(c(contig_id, sequence_well))
fa_3_patient_4 <- fa_3_patient_4 %>% select(c(contig_id, sequence_well))
fa_3_patient_7 <- fa_3_patient_7 %>% select(c(contig_id, sequence_well))
fa_3_patient_8 <- fa_3_patient_8 %>% select(c(contig_id, sequence_well))

# Remove fa_3
rm(fa_3)

###### FOR WELL 4 #########

# Make a copy of the data in well_4
fa_4 <- well_4

# Select only sample, patient ID, contig ID, sequence columns
fa_4 <- fa_4 %>% select(c(sample, patient_ID, contig_id, sequence_well))

# Separate each individual person into a different table
fa_4_healthy_1 <- subset(fa_4, patient_ID == "Hashtag1")
fa_4_healthy_3 <- subset(fa_4, patient_ID == "Hashtag3")
fa_4_healthy_5 <- subset(fa_4, patient_ID == "Hashtag5")
fa_4_healthy_6 <- subset(fa_4, patient_ID == "Hashtag6")

fa_4_patient_2 <- subset(fa_4, patient_ID == "Hashtag2")
fa_4_patient_4 <- subset(fa_4, patient_ID == "Hashtag4")
fa_4_patient_7 <- subset(fa_4, patient_ID == "Hashtag7")
fa_4_patient_8 <- subset(fa_4, patient_ID == "Hashtag8")

# Select only contig ID, sequence columns
fa_4_healthy_1 <- fa_4_healthy_1 %>% select(c(contig_id, sequence_well))
fa_4_healthy_3 <- fa_4_healthy_3 %>% select(c(contig_id, sequence_well))
fa_4_healthy_5 <- fa_4_healthy_5 %>% select(c(contig_id, sequence_well))
fa_4_healthy_6 <- fa_4_healthy_6 %>% select(c(contig_id, sequence_well))

fa_4_patient_2 <- fa_4_patient_2 %>% select(c(contig_id, sequence_well))
fa_4_patient_4 <- fa_4_patient_4 %>% select(c(contig_id, sequence_well))
fa_4_patient_7 <- fa_4_patient_7 %>% select(c(contig_id, sequence_well))
fa_4_patient_8 <- fa_4_patient_8 %>% select(c(contig_id, sequence_well))

# Remove fa_4
rm(fa_4)


###### COMBINING ALL WELLS PER PATIENT ###############

# Healthy 1
healthy_1 <- rbind(fa_1_healthy_1, fa_2_healthy_1)
healthy_1 <- rbind(healthy_1, fa_3_healthy_1)
healthy_1 <- rbind(healthy_1, fa_4_healthy_1)

rm(fa_1_healthy_1, fa_2_healthy_1, fa_3_healthy_1, fa_4_healthy_1)

# Save file
# Convert data frame to FASTA
healthy_1$contig_id <- str_prefix(healthy_1$contig_id, ">")
healthy_1 <- do.call(rbind, lapply(seq(nrow(healthy_1)), function(i) t(healthy_1[i, ])))
write.table(healthy_1, row.names = FALSE, col.names = FALSE, quote = FALSE, file = "B_VDJ/healthy/1_filtered_contig.fasta")

rm(healthy_1)

# Healthy 3
healthy_3 <- rbind(fa_1_healthy_3, fa_2_healthy_3)
healthy_3 <- rbind(healthy_3, fa_3_healthy_3)
healthy_3 <- rbind(healthy_3, fa_4_healthy_3)

rm(fa_1_healthy_3, fa_2_healthy_3, fa_3_healthy_3, fa_4_healthy_3)

# Save file
# Convert data frame to FASTA
healthy_3$contig_id <- str_prefix(healthy_3$contig_id, ">")
healthy_3 <- do.call(rbind, lapply(seq(nrow(healthy_3)), function(i) t(healthy_3[i, ])))
write.table(healthy_3, row.names = FALSE, col.names = FALSE, quote = FALSE, file = "B_VDJ/healthy/3_filtered_contig.fasta")

rm(healthy_3)

# Healthy 5
healthy_5 <- rbind(fa_1_healthy_5, fa_2_healthy_5)
healthy_5 <- rbind(healthy_5, fa_3_healthy_5)
healthy_5 <- rbind(healthy_5, fa_4_healthy_5)

rm(fa_1_healthy_5, fa_2_healthy_5, fa_3_healthy_5, fa_4_healthy_5)

# Save file
# Convert data frame to FASTA
healthy_5$contig_id <- str_prefix(healthy_5$contig_id, ">")
healthy_5 <- do.call(rbind, lapply(seq(nrow(healthy_5)), function(i) t(healthy_5[i, ])))
write.table(healthy_5, row.names = FALSE, col.names = FALSE, quote = FALSE, file = "B_VDJ/healthy/5_filtered_contig.fasta")

rm(healthy_5)


# Healthy 6
healthy_6 <- rbind(fa_1_healthy_6, fa_2_healthy_6)
healthy_6 <- rbind(healthy_6, fa_3_healthy_6)
healthy_6 <- rbind(healthy_6, fa_4_healthy_6)

rm(fa_1_healthy_6, fa_2_healthy_6, fa_3_healthy_6, fa_4_healthy_6)

# Save file
# Convert data frame to FASTA
healthy_6$contig_id <- str_prefix(healthy_6$contig_id, ">")
healthy_6 <- do.call(rbind, lapply(seq(nrow(healthy_6)), function(i) t(healthy_6[i, ])))
write.table(healthy_6, row.names = FALSE, col.names = FALSE, quote = FALSE, file = "B_VDJ/healthy/6_filtered_contig.fasta")

rm(healthy_6)


# Patient 2
patient_2 <- rbind(fa_1_patient_2, fa_2_patient_2)
patient_2 <- rbind(patient_2, fa_3_patient_2)
patient_2 <- rbind(patient_2, fa_4_patient_2)

rm(fa_1_patient_2, fa_2_patient_2, fa_3_patient_2, fa_4_patient_2)

# Save file
# Convert data frame to FASTA
patient_2$contig_id <- str_prefix(patient_2$contig_id, ">")
patient_2 <- do.call(rbind, lapply(seq(nrow(patient_2)), function(i) t(patient_2[i, ])))
write.table(patient_2, row.names = FALSE, col.names = FALSE, quote = FALSE, file = "B_VDJ/patient/2_filtered_contig.fasta")

rm(patient_2)


# Patient 4
patient_4 <- rbind(fa_1_patient_4, fa_2_patient_4)
patient_4 <- rbind(patient_4, fa_3_patient_4)
patient_4 <- rbind(patient_4, fa_4_patient_4)

rm(fa_1_patient_4, fa_2_patient_4, fa_3_patient_4, fa_4_patient_4)

# Save file
# Convert data frame to FASTA
patient_4$contig_id <- str_prefix(patient_4$contig_id, ">")
patient_4 <- do.call(rbind, lapply(seq(nrow(patient_4)), function(i) t(patient_4[i, ])))
write.table(patient_4, row.names = FALSE, col.names = FALSE, quote = FALSE, file = "B_VDJ/patient/4_filtered_contig.fasta")

rm(patient_4)

# Patient 7
patient_7 <- rbind(fa_1_patient_7, fa_2_patient_7)
patient_7 <- rbind(patient_7, fa_3_patient_7)
patient_7 <- rbind(patient_7, fa_4_patient_7)

rm(fa_1_patient_7, fa_2_patient_7, fa_3_patient_7, fa_4_patient_7)

# Save file
# Convert data frame to FASTA
patient_7$contig_id <- str_prefix(patient_7$contig_id, ">")
patient_7 <- do.call(rbind, lapply(seq(nrow(patient_7)), function(i) t(patient_7[i, ])))
write.table(patient_7, row.names = FALSE, col.names = FALSE, quote = FALSE, file = "B_VDJ/patient/7_filtered_contig.fasta")

rm(patient_7)

# Patient 8
patient_8 <- rbind(fa_1_patient_8, fa_2_patient_8)
patient_8 <- rbind(patient_8, fa_3_patient_8)
patient_8 <- rbind(patient_8, fa_4_patient_8)

rm(fa_1_patient_8, fa_2_patient_8, fa_3_patient_8, fa_4_patient_8)

# Save file
# Convert data frame to FASTA
patient_8$contig_id <- str_prefix(patient_8$contig_id, ">")
patient_8 <- do.call(rbind, lapply(seq(nrow(patient_8)), function(i) t(patient_8[i, ])))
write.table(patient_8, row.names = FALSE, col.names = FALSE, quote = FALSE, file = "B_VDJ/patient/8_filtered_contig.fasta")

rm(patient_8)




###################################################################
############## GENERATE BCR ANNOTATIONS BY PATIENT ################
###################################################################


###### FOR WELL 1 #########

# Make a copy of the data in well_1
bcr_1 <- well_1

# Deselect non-BCR columns
bcr_1 <- bcr_1 %>% select(-c(orig.ident, nCount_RNA, nFeature_RNA, nCount_ADT, nFeature_ADT, sequence_well))

# Separate each individual person into a different table
bcr_1_healthy_1 <- subset(bcr_1, patient_ID == "Hashtag1")
bcr_1_healthy_3 <- subset(bcr_1, patient_ID == "Hashtag3")
bcr_1_healthy_5 <- subset(bcr_1, patient_ID == "Hashtag5")
bcr_1_healthy_6 <- subset(bcr_1, patient_ID == "Hashtag6")

bcr_1_patient_2 <- subset(bcr_1, patient_ID == "Hashtag2")
bcr_1_patient_4 <- subset(bcr_1, patient_ID == "Hashtag4")
bcr_1_patient_7 <- subset(bcr_1, patient_ID == "Hashtag7")
bcr_1_patient_8 <- subset(bcr_1, patient_ID == "Hashtag8")

# Drop sample, patient ID columns
bcr_1_healthy_1 <- bcr_1_healthy_1 %>% select(-c(sample, patient_ID))
bcr_1_healthy_3 <- bcr_1_healthy_3 %>% select(-c(sample, patient_ID))
bcr_1_healthy_5 <- bcr_1_healthy_5 %>% select(-c(sample, patient_ID))
bcr_1_healthy_6 <- bcr_1_healthy_6 %>% select(-c(sample, patient_ID))

bcr_1_patient_2 <- bcr_1_patient_2 %>% select(-c(sample, patient_ID))
bcr_1_patient_4 <- bcr_1_patient_4 %>% select(-c(sample, patient_ID))
bcr_1_patient_7 <- bcr_1_patient_7 %>% select(-c(sample, patient_ID))
bcr_1_patient_8 <- bcr_1_patient_8 %>% select(-c(sample, patient_ID))

# Remove bcr_1
rm(bcr_1)

###### FOR WELL 2 #########

# Make a copy of the data in well_2
bcr_2 <- well_2

# Deselect non-BCR columns
bcr_2 <- bcr_2 %>% select(-c(orig.ident, nCount_RNA, nFeature_RNA, nCount_ADT, nFeature_ADT, sequence_well))

# Separate each individual person into a different table
bcr_2_healthy_1 <- subset(bcr_2, patient_ID == "Hashtag1")
bcr_2_healthy_3 <- subset(bcr_2, patient_ID == "Hashtag3")
bcr_2_healthy_5 <- subset(bcr_2, patient_ID == "Hashtag5")
bcr_2_healthy_6 <- subset(bcr_2, patient_ID == "Hashtag6")

bcr_2_patient_2 <- subset(bcr_2, patient_ID == "Hashtag2")
bcr_2_patient_4 <- subset(bcr_2, patient_ID == "Hashtag4")
bcr_2_patient_7 <- subset(bcr_2, patient_ID == "Hashtag7")
bcr_2_patient_8 <- subset(bcr_2, patient_ID == "Hashtag8")

# Drop sample, patient ID columns
bcr_2_healthy_1 <- bcr_2_healthy_1 %>% select(-c(sample, patient_ID))
bcr_2_healthy_3 <- bcr_2_healthy_3 %>% select(-c(sample, patient_ID))
bcr_2_healthy_5 <- bcr_2_healthy_5 %>% select(-c(sample, patient_ID))
bcr_2_healthy_6 <- bcr_2_healthy_6 %>% select(-c(sample, patient_ID))

bcr_2_patient_2 <- bcr_2_patient_2 %>% select(-c(sample, patient_ID))
bcr_2_patient_4 <- bcr_2_patient_4 %>% select(-c(sample, patient_ID))
bcr_2_patient_7 <- bcr_2_patient_7 %>% select(-c(sample, patient_ID))
bcr_2_patient_8 <- bcr_2_patient_8 %>% select(-c(sample, patient_ID))

# Remove bcr_2
rm(bcr_2)

###### FOR WELL 3 #########

# Make a copy of the data in well_3
bcr_3 <- well_3

# Deselect non-BCR columns
bcr_3 <- bcr_3 %>% select(-c(orig.ident, nCount_RNA, nFeature_RNA, nCount_ADT, nFeature_ADT, sequence_well))

# Separate each individual person into a different table
bcr_3_healthy_1 <- subset(bcr_3, patient_ID == "Hashtag1")
bcr_3_healthy_3 <- subset(bcr_3, patient_ID == "Hashtag3")
bcr_3_healthy_5 <- subset(bcr_3, patient_ID == "Hashtag5")
bcr_3_healthy_6 <- subset(bcr_3, patient_ID == "Hashtag6")

bcr_3_patient_2 <- subset(bcr_3, patient_ID == "Hashtag2")
bcr_3_patient_4 <- subset(bcr_3, patient_ID == "Hashtag4")
bcr_3_patient_7 <- subset(bcr_3, patient_ID == "Hashtag7")
bcr_3_patient_8 <- subset(bcr_3, patient_ID == "Hashtag8")

# Drop sample, patient ID columns
bcr_3_healthy_1 <- bcr_3_healthy_1 %>% select(-c(sample, patient_ID))
bcr_3_healthy_3 <- bcr_3_healthy_3 %>% select(-c(sample, patient_ID))
bcr_3_healthy_5 <- bcr_3_healthy_5 %>% select(-c(sample, patient_ID))
bcr_3_healthy_6 <- bcr_3_healthy_6 %>% select(-c(sample, patient_ID))

bcr_3_patient_2 <- bcr_3_patient_2 %>% select(-c(sample, patient_ID))
bcr_3_patient_4 <- bcr_3_patient_4 %>% select(-c(sample, patient_ID))
bcr_3_patient_7 <- bcr_3_patient_7 %>% select(-c(sample, patient_ID))
bcr_3_patient_8 <- bcr_3_patient_8 %>% select(-c(sample, patient_ID))

# Remove bcr_3
rm(bcr_3)


###### FOR WELL 4 #########

# Make a copy of the data in well_4
bcr_4 <- well_4

# Deselect non-BCR columns
bcr_4 <- bcr_4 %>% select(-c(orig.ident, nCount_RNA, nFeature_RNA, nCount_ADT, nFeature_ADT, sequence_well))

# Separate each individual person into a different table
bcr_4_healthy_1 <- subset(bcr_4, patient_ID == "Hashtag1")
bcr_4_healthy_3 <- subset(bcr_4, patient_ID == "Hashtag3")
bcr_4_healthy_5 <- subset(bcr_4, patient_ID == "Hashtag5")
bcr_4_healthy_6 <- subset(bcr_4, patient_ID == "Hashtag6")

bcr_4_patient_2 <- subset(bcr_4, patient_ID == "Hashtag2")
bcr_4_patient_4 <- subset(bcr_4, patient_ID == "Hashtag4")
bcr_4_patient_7 <- subset(bcr_4, patient_ID == "Hashtag7")
bcr_4_patient_8 <- subset(bcr_4, patient_ID == "Hashtag8")

# Drop sample, patient ID columns
bcr_4_healthy_1 <- bcr_4_healthy_1 %>% select(-c(sample, patient_ID))
bcr_4_healthy_3 <- bcr_4_healthy_3 %>% select(-c(sample, patient_ID))
bcr_4_healthy_5 <- bcr_4_healthy_5 %>% select(-c(sample, patient_ID))
bcr_4_healthy_6 <- bcr_4_healthy_6 %>% select(-c(sample, patient_ID))

bcr_4_patient_2 <- bcr_4_patient_2 %>% select(-c(sample, patient_ID))
bcr_4_patient_4 <- bcr_4_patient_4 %>% select(-c(sample, patient_ID))
bcr_4_patient_7 <- bcr_4_patient_7 %>% select(-c(sample, patient_ID))
bcr_4_patient_8 <- bcr_4_patient_8 %>% select(-c(sample, patient_ID))

# Remove bcr_4
rm(bcr_4)



###### COMBINING ALL WELLS PER PATIENT ###############

# Healthy 1
healthy_1 <- rbind(bcr_1_healthy_1, bcr_2_healthy_1)
healthy_1 <- rbind(healthy_1, bcr_3_healthy_1)
healthy_1 <- rbind(healthy_1, bcr_4_healthy_1)
col_order <- c("barcode","is_cell","contig_id","high_confidence","length","chain","v_gene","d_gene","j_gene","c_gene","full_length","productive","fwr1","fwr1_nt","cdr1","cdr1_nt","fwr2","fwr2_nt","cdr2","cdr2_nt","fwr3","fwr3_nt","cdr3","cdr3_nt","fwr4","fwr4_nt","reads","umis","raw_clonotype_id","raw_consensus_id","exact_subclonotype_id")
healthy_1 <- healthy_1[, col_order]

rm(bcr_1_healthy_1, bcr_2_healthy_1, bcr_3_healthy_1, bcr_4_healthy_1)

# Save file
write.csv(healthy_1, file = "B_VDJ/healthy/1_filtered_contig_annotations.csv", row.names = FALSE)

rm(healthy_1)

# Healthy 3
healthy_3 <- rbind(bcr_1_healthy_3, bcr_2_healthy_3)
healthy_3 <- rbind(healthy_3, bcr_3_healthy_3)
healthy_3 <- rbind(healthy_3, bcr_4_healthy_3)
healthy_3 <- healthy_3[, col_order]

rm(bcr_1_healthy_3, bcr_2_healthy_3, bcr_3_healthy_3, bcr_4_healthy_3)

# Save file
write.csv(healthy_3, file = "B_VDJ/healthy/3_filtered_contig_annotations.csv", row.names = FALSE)

rm(healthy_3)

# Healthy 5
healthy_5 <- rbind(bcr_1_healthy_5, bcr_2_healthy_5)
healthy_5 <- rbind(healthy_5, bcr_3_healthy_5)
healthy_5 <- rbind(healthy_5, bcr_4_healthy_5)
healthy_5 <- healthy_5[, col_order]

rm(bcr_1_healthy_5, bcr_2_healthy_5, bcr_3_healthy_5, bcr_4_healthy_5)

# Save file
write.csv(healthy_5, file = "B_VDJ/healthy/5_filtered_contig_annotations.csv", row.names = FALSE)

rm(healthy_5)


# Healthy 6
healthy_6 <- rbind(bcr_1_healthy_6, bcr_2_healthy_6)
healthy_6 <- rbind(healthy_6, bcr_3_healthy_6)
healthy_6 <- rbind(healthy_6, bcr_4_healthy_6)
healthy_6 <- healthy_6[, col_order]

rm(bcr_1_healthy_6, bcr_2_healthy_6, bcr_3_healthy_6, bcr_4_healthy_6)

# Save file
write.csv(healthy_6, file = "B_VDJ/healthy/6_filtered_contig_annotations.csv", row.names = FALSE)

rm(healthy_6)

# Patient 2
patient_2 <- rbind(bcr_1_patient_2, bcr_2_patient_2)
patient_2 <- rbind(patient_2, bcr_3_patient_2)
patient_2 <- rbind(patient_2, bcr_4_patient_2)
patient_2 <- patient_2[, col_order]

rm(bcr_1_patient_2, bcr_2_patient_2, bcr_3_patient_2, bcr_4_patient_2)

# Save file
write.csv(patient_2, file = "B_VDJ/patient/2_filtered_contig_annotations.csv", row.names = FALSE)

rm(patient_2)

# Patient 4
patient_4 <- rbind(bcr_1_patient_4, bcr_2_patient_4)
patient_4 <- rbind(patient_4, bcr_3_patient_4)
patient_4 <- rbind(patient_4, bcr_4_patient_4)
patient_4 <- patient_4[, col_order]

rm(bcr_1_patient_4, bcr_2_patient_4, bcr_3_patient_4, bcr_4_patient_4)

# Save file
write.csv(patient_4, file = "B_VDJ/patient/4_filtered_contig_annotations.csv", row.names = FALSE)

rm(patient_4)

# Patient 7
patient_7 <- rbind(bcr_1_patient_7, bcr_2_patient_7)
patient_7 <- rbind(patient_7, bcr_3_patient_7)
patient_7 <- rbind(patient_7, bcr_4_patient_7)
patient_7 <- patient_7[, col_order]

rm(bcr_1_patient_7, bcr_2_patient_7, bcr_3_patient_7, bcr_4_patient_7)

# Save file
write.csv(patient_7, file = "B_VDJ/patient/7_filtered_contig_annotations.csv", row.names = FALSE)

rm(patient_7)

# Patient 8
patient_8 <- rbind(bcr_1_patient_8, bcr_2_patient_8)
patient_8 <- rbind(patient_8, bcr_3_patient_8)
patient_8 <- rbind(patient_8, bcr_4_patient_8)
patient_8 <- patient_8[, col_order]

rm(bcr_1_patient_8, bcr_2_patient_8, bcr_3_patient_8, bcr_4_patient_8)

# Save file
write.csv(patient_8, file = "B_VDJ/patient/8_filtered_contig_annotations.csv", row.names = FALSE)

rm(patient_8)



########### END OF CODE ##################
rm(well_1, well_2, well_3, well_4, col_order)
gc()
