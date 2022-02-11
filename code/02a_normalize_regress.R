library(Seurat)
library(future)

plan("multicore")
setwd("/home/suk571/Jocelyn_scRNA/KW9079_Farmer/210624_A00873_0376_BHC3KTDSX2_KW9079/cellranger-6.0.1/GRCh38/")
options(future.globals.maxSize = 40000 * 1024^2)

DATA <- readRDS(file = "results/output_rds/sauron_filtered.rds")

regressed <- SCTransform(DATA, vars.to.regress = "mitoRatio")

saveRDS(regressed, file = "results/output_rds/sauron_regressed.rds")

rm(DATA, regressed)
invisible(gc())
