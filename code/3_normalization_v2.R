library(Seurat)
library(tidyverse)
library(RCurl)
library(cowplot)
library(future)
plan("multicore")
setwd("/home/suk571/Jocelyn_scRNA/KW9079_Farmer/210624_A00873_0376_BHC3KTDSX2_KW9079/cellranger-6.0.1/GRCh38/")
options(future.globals.maxSize = 20000 * 1024^2)
load(file = "results/output_rds/seurat_filtered_after_QC_only_protein_coding_igh.RData")

# ## Investigating which variables to regress
# seurat_phase <- NormalizeData(filtered_seurat_new)
# seurat_phase <- FindVariableFeatures(seurat_phase, selection.method = "vst", nfeatures = 2000, verbose = TRUE)
# seurat_phase <- ScaleData(seurat_phase)
# seurat_phase <- RunPCA(seurat_phase)

# # Should we regress out by well?
# DimPlot(seurat_phase, reduction = "pca", group.by= "well", split.by = "well")
# # No.
# 
# # Should we regress out by patient ID?
# DimPlot(seurat_phase, reduction = "pca", group.by= "patient_ID", split.by = "patient_ID")
# # YES
# 
# # Adding labels to mitochondrial ratio
# seurat_phase@meta.data$mitoFr <- cut(seurat_phase@meta.data$mitoRatio,
#                                      breaks=c(-Inf, 0.029, 0.036, 0.044, Inf),
#                                      labels=c("Low","Medium","Medium high", "High"))
# # Should we regress out mitochondiral genes?
# DimPlot(seurat_phase, reduction = "pca", group.by= "mitoFr", split.by = "mitoFr")
# # YES

# Split Seurat object into 2 objects - healthy and patient
split_seurat <- SplitObject(filtered_seurat_new, split.by = "sample")
# rm(seurat_phase)
rm(filtered_seurat_new)
gc()

# Perform SCTransform which performs NormalizeData(), FindVariableFeatures(), ScaleData() together - uses MNN
for (i in 1:length(split_seurat)) {
  split_seurat[[i]] <- SCTransform(split_seurat[[i]], vars.to.regress = c("mitoRatio", "patient_ID"))
}

saveRDS(split_seurat, "results/output_rds/split_seurat_reg_mito_patient.rds")
rm(split_seurat)
gc()
