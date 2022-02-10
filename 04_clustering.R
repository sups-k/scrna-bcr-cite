# Load libraries
library(Seurat)
library(tidyverse)
library(RCurl)
library(cowplot)
library(future)
library(clustree)

# Set up multi-core in RStudio. For R in shell, use plan("multicore")
plan("multisession")

# Set working directory
setwd("/home/suk571/Jocelyn_scRNA/KW9079_Farmer/210624_A00873_0376_BHC3KTDSX2_KW9079/cellranger-6.0.1/GRCh38/")

# Enable storage of large files in R memory
options(future.globals.maxSize = 40000 * 1024^2)

# Load integrated Seurat object
seurat_integrated <- readRDS(file = "results/output_rds/sauron_integrated_Jan23.rds")

# Identify number of significant harmony embeddings
## Step 1: Determine which embedding exhibits cumulative percent greater than 90% and
## % variation associated with the embedding as less than 5
pct <- seurat_integrated[["pca"]]@stdev / sum(seurat_integrated[["pca"]]@stdev) * 100
cumu <- cumsum(pct)
co1 <- which(cumu > 90 & pct < 5)[1]

## Step 2: Determine the point where the percent change in variation between the
## consecutive PCs is less than 0.1%
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1

## Step 3: Use the maximum of the 2 metrics calculated above
pcs <- max(co1, co2)

## Step 4: Remove unnecessary variables
rm(co1, co2, cumu, pct)

# Perform PCA again with the new number of principal components
seurat_integrated <- RunPCA(object = seurat_integrated, npcs = pcs, verbose = TRUE)

# Determine the K-nearest neighbor graph
seurat_integrated <- FindNeighbors(object = seurat_integrated, reduction = "pca", dims = 1:pcs)

# Determine the clusters for various resolutions (higher the number of cells, more the resolution)                    
seurat_integrated <- FindClusters(object = seurat_integrated,
                                  resolution = seq(0.05, 1, 0.05))

# Calculate UMAP again with new value for number of PCs
seurat_integrated <- RunUMAP(object = seurat_integrated, dims = 1:pcs, reduction = "pca")

# Save the new object
saveRDS(seurat_integrated, "results/output_rds/sauron_integrated_Jan23_pc43.rds")

# Check for over-clustering
clustree(seurat_integrated, prefix = "integrated_snn.res.")
# Based on results, identify which resolution to use for clustering.
