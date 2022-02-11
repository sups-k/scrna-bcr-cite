library(Seurat)
library(tidyverse)
library(viridisLite)
library(ggplot2)

options(future.globals.maxSize = 100000 * 1024^2)
setwd("~/Jocelyn_scRNA/KW9079_Farmer/210624_A00873_0376_BHC3KTDSX2_KW9079/cellranger-6.0.1/GRCh38/")
DATA <- readRDS("results/output_rds/sauron_integrated_Jan23_pc43.rds")

# Load ChangeO database file
cat('Loading ChangeO database file ...\n')
db <- read.delim("~/immcantation-python3.8.12/imm-env/makedb_output/functional/genotype/create_germlines/VDJseq_mutation_quant.tab", stringsAsFactors=F)

# Get only heavy chains
db <- db %>% dplyr::filter(LOCUS == "IGH")
cell_ids <- db$CELL

# Remove duplicated cell IDs in ChangeO database
dup_ids <- cell_ids[duplicated(cell_ids)]
db <- db[!duplicated(cell_ids), ]
cell_ids <- cell_ids[!duplicated(cell_ids)]

# Add rownames
rownames(db) <- cell_ids

# Input metadata in different table
meta <- DATA@meta.data

# Rename "cells" column
meta <- meta %>% dplyr::rename(CELL = cells)

# Merge metadata and ChangeO database
merged_db <- merge(meta, db, by = "CELL")

# Plot UMAP with mutation frequencies
# define custom plotting function
plotFeat <- function(SeurObj, featName, featMax=Inf, combineMethod='sum', colorPalette=viridis(100)){
  # SeurObj: Seurat Object
  # featName: Column name of metadata to color cells by (NA values will be light gray)
  # featMax: Max value above which feature values will be trimmed
  # combineMethod: If multiple feature names are provided in featName,
  #                'sum': sum the feature values
  #                'sep': plot features in separate plots
  # colorPalette: color palette used to color cells
  
  umap_coords <- SeurObj@reductions$umap@cell.embeddings
  featData <- merged_db[featName]
  nPlots <- 1
  if (length(featName) > 1) {
    if (combineMethod == 'sum') {
      if (length(featName) > 2) {
        featName <- 'Total mutations'
      } else {
        featName <- paste(featName, collapse=' + ')
      }
      featData <- as.data.frame(rowSums(featData)) %>% setNames('Value')
    } else if (combineMethod == 'sep') {
      nPlots <- length(featName)
      if (length(featMax) == 1) {
        featMax <- rep(featMax, nPlots)
      } else if (length(featMax) != nPlots) {
        stop('Number of elements in featMax must be one, or equal to the number of elements in featName.')
      }
    } else {stop('Invalid combineMethod!')}
  }
  
  n_plot_cols <- min(3, nPlots)
  n_plot_rows <- ceiling(nPlots/n_plot_cols)
  par(mfrow=c(n_plot_rows, n_plot_cols), mgp=c(0.5, 0, 0), mar=c(2,2,2,2))
  
  for (i in seq(nPlots)) {
    plot_data <- as.data.frame(cbind(umap_coords, featData[, i]))
    plot_data <- plot_data[order(plot_data[, 3], na.last=NA), ]
    plot_data[plot_data[,3] > featMax[i], 3] <- featMax[i]
    
    colors <- colorPalette[cut(plot_data[,3], breaks=length(colorPalette))]
    plot(umap_coords, xlim=range(umap_coords[,1]), ylim=range(umap_coords[,2]), col='grey85', pch=16, cex=0.3, yaxt='n', xaxt='n')
    par(new=T)
    plot(plot_data[,c(1:2)], xlim=range(umap_coords[,1]), ylim=range(umap_coords[,2]), col=colors, pch=16, cex=0.3, main=featName[i], yaxt='n', xaxt='n')
  }
}

# plot mutation data on UMAP
col_scale <- viridis(100, direction=-1)
png(file.path(dirname(opt$changeo_db_path), 'umap_VDJmutfreq_all.png'), res=300, units='mm', width=120, height=100)
plotFeat(DATA, featName='MU_FREQ_TOT', featMax=0.03, colorPalette=col_scale)
invisible(dev.off())


