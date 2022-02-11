#### After running 2_filtering.R #####

# Load libraries
library(Seurat)
library(dplyr)
library(scales)
library(RColorBrewer)
library(AnnotationHub)
library(ensembldb)
library(ineq)
library(vegan)
library(rafalib)
library(parallel)

## Enable parallelisation
future::plan("multisession")

## Set working directory
setwd("/home/suk571/Jocelyn_scRNA/KW9079_Farmer/210624_A00873_0376_BHC3KTDSX2_KW9079/cellranger-6.0.1/GRCh38/")

# Read the RData file
load(file = "results/output_rds/seurat_filtered_after_QC.RData")


########################################
### IDENTIFY REQUESTED GENES TO KEEP ###
########################################
# Ighd, Ighm, Ighg1, Ighg2c, Ighg2b, Ighg3, Igha, Ighe

genes_keep <- c("IGHD", "IGHM", "IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHA1", "IGHE", "IGHA2")
genes_notfound <- setdiff(genes_keep, filtered_seurat_new@assays[["RNA"]]@counts@Dimnames[[1]])
genes_keep <- filtered_seurat_new@assays[["RNA"]]@counts@Dimnames[[1]][filtered_seurat_new@assays[["RNA"]]@counts@Dimnames[[1]] %in% genes_keep]
cat("\nThe following genes will NOT be removed from the data:\n")
cat(genes_keep, "\n")
if(length(genes_notfound) > 0){
  cat("\nWARNING: The following requested genes were NOT FOUND in the data:\n")
  cat(genes_notfound, "\n")
}



###########################################
### SELECTING ONLY PROTEIN-CODING GENES ###
###########################################
# Connect to AnnotationHub
ah <- AnnotationHub()

# Access the Ensembl database for organism
ahDb <- query(ah, pattern = c("Homo sapiens", "EnsDb"), ignore.case = TRUE)
# Acquire the latest annotation files
id <- ahDb %>% mcols() %>% rownames() %>% tail(n = 1)
# Download the appropriate Ensembldb database
edb <- ah[[id]]
annotations <- genes(edb, return.type = "data.frame")
# Select annotations of interest
annotations <- annotations %>%
  dplyr::select(gene_id, gene_name, gene_biotype, seq_name, description, entrezid)


sel <- annotations[match(filtered_seurat_new@assays[["RNA"]]@counts@Dimnames[[1]], annotations$gene_name),]
sel <- sel %>% dplyr::filter(gene_biotype == "protein_coding")
genes_use <- union(as.character(na.omit(sel$gene_name)), genes_keep)
filtered_seurat_new@assays[["RNA"]]@counts <- filtered_seurat_new@assays[["RNA"]]@counts[genes_use,]
filtered_seurat_new@assays[["RNA"]]@data <- filtered_seurat_new@assays[["RNA"]]@data[genes_use,]

#---------

# Create .RData object to load at any time
save(filtered_seurat_new, file="results/output_rds/seurat_filtered_after_QC_only_protein_coding_igh.RData")

