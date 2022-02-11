## Load libraries
library(Seurat)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)
library(future)
library(stringr.tools)

# For mitochondrial reads
# library(AnnotationHub)
# library(ensembldb)

## Enable parallelisation
plan("multisession")
# plan("multicore")

## Set working directory
setwd("/home/suk571/Jocelyn_scRNA/KW9079_Farmer/210624_A00873_0376_BHC3KTDSX2_KW9079/cellranger-6.0.1/GRCh38/")

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
save(merged_seurat, file="merged_seurat.RData")

## Remove all other variables - only merged Seurat object is requried now
rm(BRI1242, BRI1245, BRI1248, BRI1251, metadata)

## Read the RData file
# load(file = "/home/suk571/Jocelyn_scRNA/KW9079_Farmer/210624_A00873_0376_BHC3KTDSX2_KW9079/cellranger-6.0.1/GRCh38/merged_seurat.RData")


############### ADDING BCR INFORMATION TO METADATA ###################

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

## ******** Below taken from https://www.biostars.org/p/384640/ *******

## In the VDJ data I have, cell barcodes are duplicated in rows because the
## sequencing of heavy and light chain genes creates multiple data points for
## the same cell. To have a more granular detail level in the meta.data slot of
## the Seurat object (such as the actual V/D/J genes sequenced), concatenate all
## the unique data we'd like to keep (which results in the duplication of the cell
## barcodes in the data frame) and collapse the duplicate rows.

# Generate a function that will concatenate unique data entries and collapse duplicate rows
# To do this, I first factorize the data and then get factor levels as unique data points
# Then I paste these data points together separated with "__" to access later on if needed

data_concater <- function(x){
  
  x <- levels(factor(x))
  
  paste(x, collapse = "__")
  
}

# Use data.table package for efficient calculations. Dplyr takes a long time in my experience.

suppressPackageStartupMessages(library(data.table))

all_clones <- as.data.table(all_clones)


# Prepare a progress bar to monitor progress (helpful for large aggregations)

grpn = uniqueN(all_clones$barcode) # uniqueN() is equivalent to
# length(unique(x)) when x is an atomic vector, and nrow(unique(x)) when x is a data.frame or data.table.
pb <- txtProgressBar(min = 0, max = grpn, style = 3)

# This code applies data_concater function per barcode to create a 
# concatenated string with the information we want to keep

all_clones_collapsed <- all_clones[, {setTxtProgressBar(pb,.GRP); lapply(.SD, data_concater)} , by = barcode]

# Assign row names for merging into combined Seurat object

rownames(all_clones_collapsed) <- all_clones_collapsed$barcode

## Add BCR information to metadata of combined Seurat object
merged_seurat <- AddMetaData(merged_seurat, metadata = all_clones_collapsed)

rm(grpn, data_concater, pb, all_clones)

#############----------- BCR INFORMATION ADDED ---------------################

## Save as .RData
save(merged_seurat, file="merged_BCR_seurat.RData")



############### QUALITY CONTROL #################

## We will create our metrics file from the metadata stored in the Seurat object.
metrics <- merged_seurat@meta.data

## We will explore the following metrics through visualizations to decide on which
## cells are low quality and should be removed from the analysis:
#
## -- Cell counts
## -- UMI counts per cell
## -- Genes detected per cell
## -- UMIs vs. genes detected
## -- Mitochondrial counts ratio
## -- Novelty

## ----------------Cell Counts --------------------

## The cell counts are determined by the number of unique cellular barcodes detected.
## Expect the number of unique cellular barcodes to be around 50-60 % of what is loaded.
## After we remove the low quality cells by filtering, we will expect the number
## of cells to be at or a bit below the number of sequenced cells.

## Visualize the number of cell counts per cell
metrics %>% 
  ggplot(aes(x = sample, fill = sample), size = 10) + 
  geom_bar(width = 0.20) + 
  ggtitle("NCells") +
  theme_bw() +
  theme(axis.text.x = element_text(size=rel(3), angle = 45, hjust = 1),
        axis.text.y = element_text(size=rel(4)),
        axis.title = element_text(size=rel(4)),
        plot.title = element_text(size=rel(4), hjust = 0.5, face = "bold"),
        legend.text = element_text(size=20),
        legend.title = element_text(size=20))

## There are 30,000 healthy and 45,000 patient cells.



## --------------- UMI counts (transcripts) per cell --------------------

##The UMI counts per cell should be above 500 as a bare minimum. Although usable,
## it’s still low if between 500-1000 counts. Then the cells probably should have
## been sequenced more deeply.

metrics %>% 
  ggplot(aes(color = sample, x = nUMI, fill = sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  ylab("log10 cell density") +
  geom_vline(xintercept = 500) +
  theme_bw() +
  theme(axis.text.x = element_text(size=30),
        axis.text.y = element_text(size=40),
        axis.title = element_text(size=40),
        plot.title = element_text(size=40, hjust = 0.5, face = "bold"),
        legend.text = element_text(size=20),
        legend.title = element_text(size=20))

## -------------- Genes detected per cell --------------------------------

## We have similar expectations for gene detection as for UMI detection,
## although it may be a bit lower than UMIs. For high quality data, the proportional
## histogram should contain a single large peak that represents cells that were
## encapsulated. If we see a small shoulder to the left of the major peak, or a
## bimodal distribution of the cells, that can indicate a couple of things.
## It might be that there are a set of cells that failed for some reason. It
## could also be that there are biologically different types of cells (i.e.,
## quiescent cell populations, less complex cells of interest), and/or one type
## is much smaller than the other (i.e. cells with high counts may be cells that
## are larger in size). Therefore, this threshold should be assessed with other metrics.


## As a histogram
metrics %>% 
  ggplot(aes(color = sample, x = nGene, fill = sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  geom_vline(xintercept = 300) +
  theme_classic() +
  theme(axis.text.x = element_text(size=30),
        axis.text.y = element_text(size=40),
        axis.title = element_text(size=40),
        plot.title = element_text(size=40, hjust = 0.5, face = "bold"),
        legend.text = element_text(size=20),
        legend.title = element_text(size=20))

## As a bar plot
metrics %>% 
  ggplot(aes(x = sample, y = log10(nGene), fill = sample)) + 
  geom_boxplot() + 
  ggtitle("NCells vs NGenes") +
  theme_bw() +
  theme(axis.text.x = element_text(size=rel(3)),
        axis.text.y = element_text(size=rel(4)),
        axis.title = element_text(size=rel(4)),
        plot.title = element_text(size=rel(4), hjust = 0.5, face = "bold"),
        legend.text = element_text(size=20),
        legend.title = element_text(size=20))

## -------------------- UMIs vs Genes Detected --------------------------------

## Plot the number of genes versus the number of UMIs (unique molecular identifiers) coloured by the fraction of
## mitochondrial reads. Mitochondrial read fractions are only high in particularly
## low count cells with few detected genes (darker colored data points). This could
## be indicative of damaged/dying cells whose cytoplasmic mRNA has leaked out
## through a broken membrane, and thus, only mRNA located in the mitochondria is
## still conserved. These cells are filtered out by our count and gene number thresholds.
## Jointly visualizing the count and gene thresholds shows the joint filtering effect.

## With this plot we also evaluate the slope of the line, and any scatter of data
## points in the bottom right hand quadrant of the plot. These cells have a high
## number of UMIs but only a few number of genes. These could be dying cells, but
## also could represent a population of a low complexity celltype (i.e red blood cells).

## Poor quality cells are likely to have low genes and UMIs per cell. Therefore,
## a poor sample is likely to have cells in the bottom left quadrant of the plot.
## Good cells should exhibit both higher number of genes per cell and higher
## numbers of UMIs. We also expect similar lines with similar slopes for all samples.

metrics %>% 
  ggplot(aes(x = nUMI, y = nGene, color = mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250) +
  facet_wrap(~sample)+
  theme(axis.text.x = element_text(size=rel(3)),
        axis.text.y = element_text(size=rel(4)),
        axis.title = element_text(size=rel(4)),
        plot.title = element_text(size=rel(4), hjust = 0.5, face = "bold"),
        legend.text = element_text(size=20),
        legend.title = element_text(size=20))

## -------------- Mitochondrial counts ratio --------------------------------

## This metric can identify whether there is a large amount of mitochondrial
## contamination from dead or dying cells. Poor quality samples for mitochondrial
## counts would have larger peaks above the 0.2 mitochondrial ratio mark, unless
## it is expected based on sample type.


metrics %>% 
  ggplot(aes(color = sample, x = mitoRatio, fill = sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  geom_vline(xintercept = 0.2) +
  theme_bw() +
  theme(axis.text.x = element_text(size=30),
        axis.text.y = element_text(size=40),
        axis.title = element_text(size=40),
        plot.title = element_text(size=40, hjust = 0.5, face = "bold"),
        legend.title = element_text(size=25),
        legend.text = element_text(size = 25))


## ------------- Novelty/Complexity ---------------------------------

## Outlier cells in samples might be cells that have a less complex RNA species
## than other cells. Sometimes we can detect contamination with low complexity
## cell types like red blood cells via this metric. Generally, we expect the
## novelty score to be above 0.80

metrics %>%
  ggplot(aes(x = log10GenesPerUMI, color = sample, fill = sample)) +
  geom_density(alpha = 0.2) +
  geom_vline(xintercept = 0.8) +
  theme_bw() +
  theme(axis.text.x = element_text(size=30),
        axis.text.y = element_text(size=40),
        axis.title = element_text(size=40),
        plot.title = element_text(size=40, hjust = 0.5, face = "bold"),
        legend.title = element_text(size=25),
        legend.text = element_text(size = 25))










