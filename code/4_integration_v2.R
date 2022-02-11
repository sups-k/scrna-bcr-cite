# Running integration using Harmony because of different patients

library(Seurat)
library(tidyverse)
library(RCurl)
library(cowplot)
library(future)
library(harmony)

plan("multicore")
setwd("/home/suk571/Jocelyn_scRNA/KW9079_Farmer/210624_A00873_0376_BHC3KTDSX2_KW9079/cellranger-6.0.1/GRCh38/")

split_seurat <- readRDS(file = "results/output_rds/split_seurat_reg_mito_patient.rds")

options(future.globals.maxSize = 40000 * 1024^2)
gc()


# Select the most variable features to use for integration
integ_features <- SelectIntegrationFeatures(object.list = split_seurat, 
                                            nfeatures = 3000) 

# Prepare the SCT list object for integration
split_seurat <- PrepSCTIntegration(object.list = split_seurat, 
                                   anchor.features = integ_features)
# Find best buddies - can take a while to run - CCA
integ_anchors <- FindIntegrationAnchors(object.list = split_seurat, 
                                        normalization.method = "SCT", 
                                        anchor.features = integ_features)
# Integrate across conditions
seurat_integrated <- IntegrateData(anchorset = integ_anchors, 
                                   normalization.method = "SCT")

# # Identify number of significant harmony embeddings
# ## Step 1: Determine which embedding exhibits cumulative percent greater than 90% and
# ## % variation associated with the embedding as less than 5
# pct <- seurat_integrated[["pca"]]@stdev / sum(seurat_integrated[["pca"]]@stdev) * 100
# cumu <- cumsum(pct)
# co1 <- which(cumu > 90 & pct < 5)[1]
# 
# ## Step 2: Determine the point where the percent change in variation between the
# ## consecutive PCs is less than 0.1%
# co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
# 
# ## Step 3: Use the minimum of the 2 metrics calculated above
# pcs <- min(co1, co2)
# 
# ## Step 4: Remove unnecessary variables
# rm(co1, co2, cumu, pct)
# 
# # Run PCA
# seurat_integrated <- RunPCA(object = seurat_integrated)
# 
# # Run UMAP
# #seurat_integrated <- RunUMAP(seurat_integrated, dims = 1:40, reduction = "pca")


# Change active.ident from "well" to sample type
seurat_integrated <- SetIdent(seurat_integrated, value = seurat_integrated@meta.data$sample)

# Save integrated Seurat object
saveRDS(seurat_integrated, "results/output_rds/seurat_integrated_v2.rds")
gc()
