# Using 2 cores, 80G, it takes 17 min to run SCTransform

library(Seurat)
library(tidyverse)
library(RCurl)
library(cowplot)
library(future)
plan("multicore")
setwd("/home/suk571/Jocelyn_scRNA/KW9079_Farmer/210624_A00873_0376_BHC3KTDSX2_KW9079/cellranger-6.0.1/GRCh38/")
load(file = "seurat_filtered_after_QC.RData")
filtered_seurat_new@meta.data$mitoFr <- cut(filtered_seurat_new@meta.data$mitoRatio, 
                                         breaks=c(-Inf, 0.029, 0.036, 0.044, Inf), 
                                         labels=c("Low","Medium","Medium high", "High"))
split_seurat <- SplitObject(filtered_seurat_new, split.by = "sample")
rm(filtered_seurat_new)
gc()
options(future.globals.maxSize = 20000 * 1024^2)
for (i in 1:length(split_seurat)) {
     split_seurat[[i]] <- SCTransform(split_seurat[[i]], vars.to.regress = c("mitoRatio"))
}

saveRDS(split_seurat, "split_seurat.rds")
