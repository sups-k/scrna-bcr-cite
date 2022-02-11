# Load libraries
library(Seurat)
library(tidyverse)
library(RCurl)
library(cowplot)
library(future)
library(AnnotationHub)
library(ensembldb)
library(multtest) # for clusters
library(metap) # for clusters
library(limma) # for clusters

# Set up multi-core
plan("multicore")

# Set working directory
setwd("/home/suk571/Jocelyn_scRNA/KW9079_Farmer/210624_A00873_0376_BHC3KTDSX2_KW9079/cellranger-6.0.1/GRCh38/")

# Enable storage of large files in R memory
options(future.globals.maxSize = 40000 * 1024^2)

# Load integrated Seurat object
seurat_integrated <- readRDS(file = "results/output_rds/pc10_seurat_integrated.rds")

## Create gene annotations file
# Connect to AnnotationHub
ah <- AnnotationHub()

# Access the Ensembl database for organism
ahDb <- query(ah, 
              pattern = c("Homo sapiens", "EnsDb"), 
              ignore.case = TRUE)

# Acquire the latest annotation files
id <- ahDb %>%
  mcols() %>%
  rownames() %>%
  tail(n = 1)

# Download the appropriate Ensembldb database
edb <- ah[[id]]

# Extract gene-level information from database
annotations <- genes(edb, 
                     return.type = "data.frame")

# Select annotations of interest
annotations <- annotations %>%
  dplyr::select(gene_id, gene_name, seq_name, gene_biotype, description)

# Remove unwanted variables
rm(ah, ahDb, edb, id)

# Select the RNA counts slot to be the default assay for finding markers. Normally it is "integrated"
DefaultAssay(seurat_integrated) <- "RNA"

# Normalize RNA data for visualization purposes, DE analysis, markers
seurat_integrated <- NormalizeData(seurat_integrated, verbose = TRUE)

# Must find conserved markers since we are comparing 2 different conditions
# Create function to get conserved markers for any given cluster
get_conserved <- function(cluster){
  FindConservedMarkers(object = seurat_integrated,
                       ident.1 = cluster,
                       grouping.var = "sample",
                       only.pos = TRUE,
                       assay = "RNA") %>%
    rownames_to_column(var = "gene") %>%
    left_join(y = unique(annotations[, c("gene_name", "description")]),
              by = c("gene" = "gene_name")) %>%
    cbind(cluster_id = cluster, .)
}

# Iterate function across desired clusters - in this case, 28 clusters at 0.8 resolution
Idents(object = seurat_integrated) <- "integrated_snn_res.0.8"
conserved_markers <- map_dfr(0:22, get_conserved)
### In some cases you will have clusters that do not have enough cells for a
### particular group - and your function will fail. For these clusters you will
### need to use FindAllMarkers().

# Extract top 10 markers per cluster
top10 <- conserved_markers %>% 
  mutate(avg_fc = (healthy_avg_log2FC + patient_avg_log2FC) /2) %>% 
  group_by(cluster_id) %>% 
  top_n(n = 10, wt = avg_fc)

# Save the top 10 markers per cluster
write.csv(conserved_markers, file = "results/top10_markers_per_cluster_res0.05.csv", row.names = FALSE)


for (i in seq(from = 0, to = 20, by = 1)) {
  
  # Genes different between conditions
  cells.1 <- colnames(seurat_integrated)[which(seurat_integrated@meta.data$integrated_snn_res.0.4 == i
                                             & seurat_integrated@meta.data$sample == "patient")]
  cells.2 <- colnames(seurat_integrated)[which(seurat_integrated@meta.data$integrated_snn_res.0.4 == i
                                             & seurat_integrated@meta.data$sample == "healthy")]
  cluster <- FindMarkers(seurat_integrated,
                          ident.1 = cells.1,
                          ident.2 = cells.2)
  # Add gene symbols to the DE table
  cluster <- cluster %>%
  rownames_to_column(var = "gene") %>%
  left_join(y = unique(annotations[, c("gene_name", "description")]),
            by = c("gene" = "gene_name"))

  # Reorder columns and sort by padj      
  cluster <- cluster[, c(1, 3:5,2,6:7)]

  cluster <- cluster %>%
  dplyr::arrange(p_val_adj) 

  write.csv(cluster, file = paste0("results/markers/cluster",i,"_P_vs_H_res0.4.csv"), row.names = FALSE)
  rm(cluster, cells.1, cells.2)
  gc()
}
