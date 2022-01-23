## Load libraries
library(Seurat)
library(tidyverse)
library(Matrix)
library(scales)
library(RCurl)
library(stringr.tools)


## Enable parallelisation
future::plan("multisession")
# future::plan("multicore")

## Set working directory
setwd("/home/suk571/Jocelyn_scRNA/KW9079_Farmer/210624_A00873_0376_BHC3KTDSX2_KW9079/cellranger-6.0.1/GRCh38/")

## Create each individual Seurat object for every well
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

################ DE-MULTIPLEXING EACH WELL BY HEALTHY VS PATIENT ########################

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

# Patient ID of the 1st cell
# rownames(BRI1242_subset_counts)[which(BRI1242_subset_counts[,1] == max(BRI1242_subset_counts[,1]), arr.ind = TRUE)]

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

## Add metadata back to Seurat object
BRI1242@meta.data <- BRI1242_metadata

# Remove data frames
rm(BRI1242_metadata, BRI1242_subset_counts, i)



##### METADATA OF SECOND WELL #######

## Copy the metadata into a separate data frame to avoid messing up the original
BRI1245_metadata <- BRI1245@meta.data

## Add sample and patient ID columns to the metadata
BRI1245_metadata$sample <- NA
BRI1245_metadata$patient_ID <- NA

## TO FIND OUT PATIENT ID:
# Create a new data frame from the counts matrix of ADT
# Take a column-wise max number - which row did max occur? That row = patient ID
# But exclude 1st 4 rows because they are not patient IDs
BRI1245_subset_counts <- as.data.frame(BRI1245@assays[["ADT"]]@counts)
BRI1245_subset_counts <- BRI1245_subset_counts[-c(1,2,3,4),]

# Patient ID of the 1st cell
# rownames(BRI1245_subset_counts)[which(BRI1245_subset_counts[,1] == max(BRI1245_subset_counts[,1]), arr.ind = TRUE)]

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

## Add metadata back to Seurat object
BRI1245@meta.data <- BRI1245_metadata

# Remove data frames
rm(BRI1245_metadata, BRI1245_subset_counts, i)


##### METADATA OF THIRD WELL #######

## Copy the metadata into a separate data frame to avoid messing up the original
BRI1248_metadata <- BRI1248@meta.data

## Add sample and patient ID columns to the metadata
BRI1248_metadata$sample <- NA
BRI1248_metadata$patient_ID <- NA

## TO FIND OUT PATIENT ID:
# Create a new data frame from the counts matrix of ADT
# Take a column-wise max number - which row did max occur? That row = patient ID
# But exclude 1st 4 rows because they are not patient IDs
BRI1248_subset_counts <- as.data.frame(BRI1248@assays[["ADT"]]@counts)
BRI1248_subset_counts <- BRI1248_subset_counts[-c(1,2,3,4),]

# Patient ID of the 1st cell
# rownames(BRI1248_subset_counts)[which(BRI1248_subset_counts[,1] == max(BRI1248_subset_counts[,1]), arr.ind = TRUE)]

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

## Add metadata back to Seurat object
BRI1248@meta.data <- BRI1248_metadata

# Remove data frames
rm(BRI1248_metadata, BRI1248_subset_counts, i)

##### METADATA OF FOURTH (and final) WELL #######

## Copy the metadata into a separate data frame to avoid messing up the original
BRI1251_metadata <- BRI1251@meta.data

## Add sample and patient ID columns to the metadata
BRI1251_metadata$sample <- NA
BRI1251_metadata$patient_ID <- NA

## TO FIND OUT PATIENT ID:
# Create a new data frame from the counts matrix of ADT
# Take a column-wise max number - which row did max occur? That row = patient ID
# But exclude 1st 4 rows because they are not patient IDs
BRI1251_subset_counts <- as.data.frame(BRI1251@assays[["ADT"]]@counts)
BRI1251_subset_counts <- BRI1251_subset_counts[-c(1,2,3,4),]

# Patient ID of the 1st cell
# rownames(BRI1251_subset_counts)[which(BRI1251_subset_counts[,1] == max(BRI1251_subset_counts[,1]), arr.ind = TRUE)]

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

## Add metadata back to Seurat object
BRI1251@meta.data <- BRI1251_metadata

# Remove data frames
rm(BRI1251_metadata, BRI1251_subset_counts, i)

############------------------- DE-MUX COMPLETE ----------------################

## Create a merged Seurat object
merged_seurat <- merge(x = BRI1242,
                       y = c(BRI1245, BRI1248, BRI1251),
                       add.cell.id = c("well_1", "well_2", "well_3", "well_4"))


## Calculate number of genes detected per UMI
merged_seurat$log10GenesPerUMI <- log10(merged_seurat$nFeature_RNA) / log10(merged_seurat$nCount_RNA)

## Compute percent mitochondrial ratio
merged_seurat$mitoRatio <- PercentageFeatureSet(object = merged_seurat, pattern = "^MT-")
merged_seurat$mitoRatio <- merged_seurat@meta.data$mitoRatio / 100

## Create metadata dataframe
metadata <- merged_seurat@meta.data

## Add cell IDs to metadata
metadata$cells <- rownames(metadata)

## Create well column
metadata$well <- NA
metadata$well[which(str_detect(metadata$cells, "^well_1_"))] <- "well_1"
metadata$well[which(str_detect(metadata$cells, "^well_2_"))] <- "well_2"
metadata$well[which(str_detect(metadata$cells, "^well_3_"))] <- "well_3"
metadata$well[which(str_detect(metadata$cells, "^well_4_"))] <- "well_4"

## Rename columns
metadata <- metadata %>%
  dplyr::rename(seq_folder = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA)

## Add metadata back to Seurat object
merged_seurat@meta.data <- metadata

## Save as .RData
save(merged_seurat, file="results/output_rds/merged_seurat.RData")

## Remove all other variables - only merged Seurat object is requried now
rm(BRI1242, BRI1245, BRI1248, BRI1251, metadata)
