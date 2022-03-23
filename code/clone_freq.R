# Load libraries
library(Seurat)
library(tidyverse)
library(alakazam)

# Set working directory
setwd("/home/suk571/Jocelyn_scRNA/KW9079_Farmer/210624_A00873_0376_BHC3KTDSX2_KW9079/cellranger-6.0.1/GRCh38/")

# Increase maximum allowed memory to 20GB
options(future.globals.maxSize = 20000 * 1024^2)

# Load integrated Seurat object
seurat_integrated <- readRDS(file = "results/output_rds/sauron_integrated_Jan23_pc43.rds")

# Save metadata in a separate variable
meta <- seurat_integrated@meta.data

# Save cluster information from metadata
meta_clusters <- meta[,c("cells", "integrated_snn_res.0.4", "sample", "patient_ID")]

# Remove large Seurat object
rm(seurat_integrated)
gc()

# Read BCR-seq data
db <- read.delim("~/immcantation-python3.8.12/imm-env/makedb_output/functional/genotype/create_germlines/VDJseq_mutation_quant.tab", stringsAsFactors=F)

# Subset the BCR-seq data to include relevant information only
small_db <- db %>% dplyr::select(SEQUENCE_ID, LOCUS, V_CALL_GENOTYPED, D_CALL, J_CALL, C_CALL, SEQUENCE_IMGT, JUNCTION, JUNCTION_LENGTH, GERMLINE_IMGT, CLONE, MU_FREQ_TOT, CELL)
# Rename 'CELL' column to match 'cells' column in `meta_clusters` for merging
small_db <- small_db %>% dplyr::rename(cells = CELL)

# Merge BCR-seq data with metadata to place each BCR clone in a particular cluster
merged <- merge(meta_clusters, small_db, by = "cells")

# Subset only heavy chains for clones
# heavy <- merged %>% dplyr::filter(LOCUS == "IGH")

# Create a subset of the merged table that only includes information on cell name, cluster, V(D)J, locus and clone
subset <- merged %>% dplyr::select(cells, integrated_snn_res.0.4, sample, patient_ID, CLONE, LOCUS, V_CALL, D_CALL, J_CALL, C_CALL)
# subset <- subset %>% drop_na(CLONE)

# Subset the above table for cluster 0
cluster_0 <- subset %>% dplyr::filter(integrated_snn_res.0.4 == 0)

# List the top 10 most frequently occurring clones in cluster 0
sort(table(cluster_0$CLONE), decreasing = TRUE)[1:10]
# OR use below code
clones <- countClones(cluster_0, clone = "CLONE")

# curve <- estimateAbundance(cluster_0, ci=0.95, nboot=100, clone = "CLONE")
# plot(curve)

# Identify major isotypes (occur in minimum 30 clones) in cluster 0
# Requires that there are no NA values so run `drop_na()`
isotype_test <- alphaDiversity(cluster_0, group="C_CALL", min_q=0, max_q=2, step_q=1, nboot=100, clone="CLONE")
plot(isotype_test)

