library(Seurat)
library(ggplot2)

future::plan("multisession")
setwd("/home/suk571/Jocelyn_scRNA/KW9079_Farmer/210624_A00873_0376_BHC3KTDSX2_KW9079/cellranger-6.0.1/GRCh38/")
options(future.globals.maxSize = 40000 * 1024^2)


DATA <- readRDS(file = "results/output_rds/sauron_integrated_Jan23_pc43.rds")


#######################################
### REMOVE CELLS FROM SEURAT OBJECT ###
#######################################
cat("\n### REMOVING CELLS FROM THE DATA ###\n")

# Determine which cells should be marked as "remove"

meta_cols <- "integrated_snn_res.0.4"
meta_cols <- colnames(DATA@meta.data)[match(meta_cols, casefold(colnames(DATA@meta.data)))]
if (any(is.na(meta_cols))) {
  cat("Could not find the following columns in meta.data:\n")
  cat(unlist(lapply(remove_meta, function(x) x[1]))[is.na(meta_cols)])
  stop('Invalid meta data field.')
}
remove_meta <- c(19, 20)
remove_cells <- matrix(FALSE, nrow=nrow(DATA@meta.data), ncol=length(remove_meta))
for (i in 1:length(remove_meta)) {
  remove_cells[, i] <- DATA@meta.data[[meta_cols]] %in% remove_meta[[i]]
}
remove_cells <- rowSums(remove_cells) > 0

# Plot UMAP showing cells to be removed
DATA@meta.data$remove_cells <- ifelse(remove_cells, "Remove", "Keep")
DimPlot(object=DATA, pt.size=0.1, reduction='umap', group.by='remove_cells', cols=c('grey85','firebrick'))


########################################################
### SAVING LIST OF REMOVED CELLS ###
########################################################
cat("\n### Saving the list of removed cell barcodes ###\n")

write.table(rownames(DATA@meta.data)[remove_cells], file="results/remove_cell_barcodes.txt", row.names=F, col.names=F, quote=F)

rm(DATA, i, meta_cols, remove_cells, remove_meta)
invisible(gc())

########################################################
### LOADING INITIAL SEURAT OBJECT & REMOVING CELLS ###
########################################################

cat("Loading raw Seurat object...")
load(file = "results/output_rds/merged_seurat.RData")
DATA <- merged_seurat
rm(merged_seurat)
invisible(gc())

cat("\nRemoving specified cells from the Seurat object\n")

excl_barcodes <- read.delim("results/remove_cell_barcodes.txt", stringsAsFactors=F, header=F)[[1]]
remove_cells <- colnames(DATA) %in% excl_barcodes
DATA <- subset(DATA, cells=colnames(DATA)[!remove_cells])
cat("Removed", sum(remove_cells), "cells from the Seurat object that were found in the provided list of", length(excl_barcodes), "cell barcodes\n")


cat("\nSaving the RAW Seurat object ...\n")
saveRDS(DATA, file = "results/output_rds/06_removed_clusters19-20_raw.rds" )


#---------
