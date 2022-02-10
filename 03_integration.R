# Load library
library(Seurat)

gc()
# Use parallel computing on cluster. Do not use in RStudio
future::plan("multicore")

# Set working directory
setwd("/home/suk571/Jocelyn_scRNA/KW9079_Farmer/210624_A00873_0376_BHC3KTDSX2_KW9079/cellranger-6.0.1/GRCh38/")

# Increase maximum allowed object size in R to 100 GB
options(future.globals.maxSize = 100000 * 1024^2)

# Loading object
cat("Reading RDS file ...\n")
split_seurat <- readRDS(file = "results/output_rds/split_sauron_regressed.rds")

# Select the most variable features to use for integration
cat("Selecting the most variable features to use for integration ...\n")
integ_features <- SelectIntegrationFeatures(object.list = split_seurat, nfeatures = 3000)

# Prepare the SCT list object for integration
cat("Prepare the SCT list object for integration ... \n")
split_seurat <- PrepSCTIntegration(object.list = split_seurat, anchor.features = integ_features)

# Find best buddies - can take a while to run
cat("Finding integration anchors. Takes a long time to run ... \n")
integ_anchors <- FindIntegrationAnchors(object.list = split_seurat, normalization.method = "SCT", anchor.features = integ_features)

# Integrate across patients & conditions
cat("Running integration... \n")
seurat_integrated <- IntegrateData(anchorset = integ_anchors, normalization.method = "SCT")

# Run PCA
cat("Running regular PCA... \n")
seurat_integrated <- RunPCA(object = seurat_integrated)

# Save integrated object
cat("Saving integrated object... \n")
saveRDS(seurat_integrated, "results/output_rds/sauron_integrated_Jan23.rds")

cat("Clearing everything ...\n")
rm(split_seurat, seurat_integrated)
gc()

cat("END \n")
