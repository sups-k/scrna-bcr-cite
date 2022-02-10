
# Load library
library(Seurat)

gc()
# Use parallel computing on cluster. Do not use in RStudio
future::plan("multicore")

# Set working directory
setwd("/home/suk571/Jocelyn_scRNA/KW9079_Farmer/210624_A00873_0376_BHC3KTDSX2_KW9079/cellranger-6.0.1/GRCh38/")

# Increase maximum allowed object size in R to 100 GB
options(future.globals.maxSize = 100000 * 1024^2)

# Read Seurat object
cat("Reading RDS file ...\n")
DATA <- readRDS(file = "results/output_rds/split_sauron_regressed.rds")

# Split Seurat object into 8 objects - for each sample
cat("Splitting Seurat object into 8 objects for each sample ...\n")
split_seurat <- SplitObject(filtered_seurat_new, split.by = "patient_ID")
rm(filtered_seurat_new)
gc()

# Perform SCTransform which performs NormalizeData(), FindVariableFeatures(), ScaleData() together - uses MNN
cat("Running SCTransform() to regress mitochondrial ratio. May take a while ...\n")
for (i in 1:length(split_seurat)) {
  split_seurat[[i]] <- SCTransform(split_seurat[[i]], vars.to.regress = "mitoRatio")
}

cat("Saving regressed object ...\n")
saveRDS(split_seurat, "results/output_rds/split_seurat_reg_mito_patient.rds")

cat("Clearing everything ...\n")
rm(split_seurat)
gc()
cat("END \n")
