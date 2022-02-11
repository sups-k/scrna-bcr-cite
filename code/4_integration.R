library(Seurat)
library(tidyverse)
library(RCurl)
library(cowplot)
library(future)

plan("multicore")
setwd("/home/suk571/Jocelyn_scRNA/KW9079_Farmer/210624_A00873_0376_BHC3KTDSX2_KW9079/cellranger-6.0.1/GRCh38/")

split_seurat <- readRDS(file = "split_seurat.rds")

options(future.globals.maxSize = 40000 * 1024^2)
gc()

# Select the most variable features to use for integration
integ_features <- SelectIntegrationFeatures(object.list = split_seurat, 
                                            nfeatures = 3000) 

# Prepare the SCT list object for integration
split_seurat <- PrepSCTIntegration(object.list = split_seurat, 
                                   anchor.features = integ_features)
# Find best buddies - can take a while to run
integ_anchors <- FindIntegrationAnchors(object.list = split_seurat, 
                                        normalization.method = "SCT", 
                                        anchor.features = integ_features)
# Integrate across conditions
seurat_integrated <- IntegrateData(anchorset = integ_anchors, 
                                   normalization.method = "SCT")

# Run PCA
seurat_integrated <- RunPCA(object = seurat_integrated)

# Run UMAP
#seurat_integrated <- RunUMAP(seurat_integrated, dims = 1:40, reduction = "pca")

# Change active.ident from "well" to sample type
seurat_integrated <- SetIdent(seurat_integrated, value = seurat_integrated@meta.data$sample)

# Save integrated Seurat object
saveRDS(seurat_integrated, "seurat_integrated.rds")
gc()
