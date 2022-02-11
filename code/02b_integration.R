library(Seurat)
library(future)
library(harmony)
library(tidyverse)

plan("multicore")
setwd("/home/suk571/Jocelyn_scRNA/KW9079_Farmer/210624_A00873_0376_BHC3KTDSX2_KW9079/cellranger-6.0.1/GRCh38/")
options(future.globals.maxSize = 40000 * 1024^2)

DATA <- readRDS(file = "results/output_rds/sauron_regressed.rds")
# Removing "filtered" column from metadata
DATA@meta.data <- DATA@meta.data[,-50]

# Check whether mitochondiral ratio was regressed
seurat_phase <- RunPCA(object = DATA)
seurat_phase@meta.data$mitoFr <- cut(seurat_phase@meta.data$mitoRatio,
                                     breaks=c(-Inf, 0.029, 0.036, 0.044, Inf),
                                     labels=c("Low","Medium","Medium high", "High"))
DimPlot(seurat_phase, reduction = "pca", group.by= "mitoFr", split.by = "mitoFr")
# Yes it was

# Removing "mitoFr" column from metadata
seurat_phase@meta.data <- seurat_phase@meta.data[,-52]


# Identify number of significant PCs
## Step 1: Determine which PC exhibits cumulative percent greater than 90% and
## % variation associated with the PC as less than 5
pct <- seurat_phase[["pca"]]@stdev / sum(seurat_phase[["pca"]]@stdev) * 100
cumu <- cumsum(pct)
co1 <- which(cumu > 90 & pct < 5)[1]

## Step 2: Determine the point where the percent change in variation between the
## consecutive PCs is less than 0.1%
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1

## Step 3: Use the minimum of the 2 metrics calculated above
pcs <- max(co1, co2)

## Step 4: Remove unnecessary variables
rm(co1, co2, cumu, pct)

# Running PCA again with new number of PCs
seurat_phase <- RunPCA(object = DATA, npcs = pcs)


# Running harmony integration
integ <- RunHarmony(seurat_phase, "patient_ID", plot_convergence = TRUE, assay.use = "SCT")
harmony_embeddings <- Embeddings(integ, "harmony")

# Check integration
DimPlot(object = integ, reduction = "harmony", pt.size = .1, group.by = "patient_ID")

# Integration needed.

integ <- integ %>% 
  RunUMAP(reduction = "harmony", dims = 1:pcs) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:pcs) %>% 
  FindClusters(resolution = 0.5)
integ <- SetIdent(integ, value = integ@meta.data$sample)
DimPlot(integ, reduction = "umap", label = TRUE, label.size = 6, split.by = "sample")
