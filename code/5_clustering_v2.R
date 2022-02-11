# Load libraries
library(Seurat)
library(tidyverse)
library(RCurl)
library(cowplot)
library(future)

# Set up multi-core
plan("multisession")

# Set working directory
setwd("/home/suk571/Jocelyn_scRNA/KW9079_Farmer/210624_A00873_0376_BHC3KTDSX2_KW9079/cellranger-6.0.1/GRCh38/")

# Enable storage of large files in R memory
options(future.globals.maxSize = 40000 * 1024^2)

# Load integrated Seurat object
seurat_integrated <- readRDS(file = "results/output_rds/seurat_integrated_v2.rds")

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
                                  resolution = c(0.05, 0.1, 0.4, 0.6, 0.8, 1.0))

# Calculate UMAP again with new value for number of PCs
seurat_integrated <- RunUMAP(object = seurat_integrated, dims = 1:pcs, reduction = "pca")

# Assign resolution of clusters
# Idents(object = seurat_integrated) <- "integrated_snn_res.0.8"

# Plot UMAP
# pdf(file = "/home/suk571/Jocelyn_scRNA/KW9079_Farmer/210624_A00873_0376_BHC3KTDSX2_KW9079/cellranger-6.0.1/GRCh38/results/UMAP_clusters_res_1.6.pdf", width = 10, height = 10)
# DimPlot(seurat_integrated, reduction = "umap", label = TRUE, label.size = 6, split.by = "sample")
# dev.off()

# Save the new object
saveRDS(seurat_integrated, "results/output_rds/pc43_seurat_integrated_v2.rds")

