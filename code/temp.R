library(Seurat)
library(tidyverse)
library(RCurl)
library(cowplot)

gc()
future::plan("multicore")
setwd("/home/suk571/Jocelyn_scRNA/KW9079_Farmer/210624_A00873_0376_BHC3KTDSX2_KW9079/cellranger-6.0.1/GRCh38/")

options(future.globals.maxSize = 100000 * 1024^2)

# Loading object
cat("Reading RDS file ...\n")
DATA <- readRDS(file = "results/output_rds/sauron_filtered.rds")

# Split Seurat object into 8 objects - for each sample
cat("Splitting into 8 objects...\n")
split_seurat <- SplitObject(DATA, split.by = "patient_ID")

cat("Removing DATA object...\n")
rm(DATA)
invisible(gc())

# Perform SCTransform which performs NormalizeData(), FindVariableFeatures(), ScaleData() together - uses MNN
for (i in 1:length(split_seurat)) {
  cat(paste0("Performing SCTransform on object ", i, "...\n"))
  split_seurat[[i]] <- SCTransform(split_seurat[[i]], vars.to.regress = c("mitoRatio", "CC.Diff"))
}

cat("Saving regressed object...\n")
saveRDS(split_seurat, "results/output_rds/Jan23_split_seurat_reg_mito_CC.rds")
invisible(gc())


# Select the most variable features to use for integration
cat("Selecting the most variable features to use for integration ...\n")
integ_features <- SelectIntegrationFeatures(object.list = split_seurat, nfeatures = 3000)
# Prepare the SCT list object for integration
cat("Prepare the SCT list object for integration ... \n")
split_seurat <- PrepSCTIntegration(object.list = split_seurat, anchor.features = integ_features)
# Find best buddies - can take a while to run
cat("Finding integration anchors. Takes a long time to run ... \n")
integ_anchors <- FindIntegrationAnchors(object.list = split_seurat, normalization.method = "SCT", anchor.features = integ_features)

# Integrate across conditions
cat("Running integration... \n")
seurat_integrated <- IntegrateData(anchorset = integ_anchors, normalization.method = "SCT")

cat("Removing split_seurat object...\n")
rm(split_seurat)
invisible(gc())

# Save integrated object
cat("Saving integrated object without running PCA... \n")
saveRDS(seurat_integrated, "results/output_rds/Jan23_sauron_integrated_noPCA.rds")

cat("Clearing everything ...\n")
rm(seurat_integrated)
invisible(gc())

cat("END \n")
