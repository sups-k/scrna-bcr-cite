# Load libraries
library(Seurat)
library(tidyverse)
library(RCurl)
library(cowplot)
library(future)
library(AnnotationHub)
library(ensembldb)
library(multtest) # for clusters
library(metap) # for clusters
library(limma) # for clusters

# Set up multi-core
plan("multicore")

# Set working directory
setwd("/home/suk571/Jocelyn_scRNA/KW9079_Farmer/210624_A00873_0376_BHC3KTDSX2_KW9079/cellranger-6.0.1/GRCh38/")

# Enable storage of large files in R memory
options(future.globals.maxSize = 30000 * 1024^2)

# Load integrated Seurat object
seurat_integrated <- readRDS(file = "pc10_integrated_seurat.rds")

# Extract identity and sample information from seurat object to determine the number of cells per cluster per sample
n_cells <- FetchData(seurat_integrated, 
                     vars = c("ident", "sample")) %>%
  dplyr::count(ident, sample) %>%
  tidyr::spread(ident, n)

# View table
View(n_cells)

# Determine metrics to plot present in seurat_integrated@meta.data
metrics <-  c("nUMI", "nGene", "mitoRatio")

FeaturePlot(seurat_integrated, 
            reduction = "umap", 
            features = metrics,
            pt.size = 0.4, 
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE)
## Population 20 has high nUMI and nGene. Opposite for population 1. MitoRatio
## high in population 1.

# Defining the information in the seurat object of interest
columns <- c(paste0("PC_", 1:10),
             "ident",
             "UMAP_1", "UMAP_2")

# Extracting this data from the seurat object
pc_data <- FetchData(seurat_integrated, 
                     vars = columns)

# Adding cluster label to center of cluster on UMAP
umap_label <- FetchData(seurat_integrated, 
                        vars = c("ident", "UMAP_1", "UMAP_2"))  %>%
  group_by(ident) %>%
  summarise(x=mean(UMAP_1), y=mean(UMAP_2))

# Plotting a UMAP plot for each of the PCs
map(paste0("PC_", 1:10), function(pc){
  ggplot(pc_data, 
         aes(UMAP_1, UMAP_2)) +
    geom_point(aes_string(color=pc), 
               alpha = 0.7) +
    scale_color_gradient(low = "grey90", 
                         high = "blue")  +
    geom_text(data=umap_label, 
              aes(label=ident, x, y)) +
    ggtitle(pc)
}) %>% 
  plot_grid(plotlist = .)


## Create gene annotations file
# Connect to AnnotationHub
ah <- AnnotationHub()

# Access the Ensembl database for organism
ahDb <- query(ah, 
              pattern = c("Homo sapiens", "EnsDb"), 
              ignore.case = TRUE)

# Acquire the latest annotation files
id <- ahDb %>%
  mcols() %>%
  rownames() %>%
  tail(n = 1)

# Download the appropriate Ensembldb database
edb <- ah[[id]]

# Extract gene-level information from database
annotations <- genes(edb, 
                     return.type = "data.frame")

# Select annotations of interest
annotations <- annotations %>%
  dplyr::select(gene_id, gene_name, seq_name, gene_biotype, description)

# Remove unwanted variables
rm(ah, ahDb, edb, id)

# Select the RNA counts slot to be the default assay for finding markers. Normally it is "integrated"
DefaultAssay(seurat_integrated) <- "RNA"

# Normalize RNA data for visualization purposes, DE analysis, markers
seurat_integrated <- NormalizeData(seurat_integrated, verbose = TRUE)

# Must find conserved markers since we are comparing 2 different conditions
# Create function to get conserved markers for any given cluster
get_conserved <- function(cluster){
  FindConservedMarkers(object = seurat_integrated,
                       ident.1 = cluster,
                       grouping.var = "sample",
                       only.pos = TRUE,
                       assay = "RNA") %>%
    rownames_to_column(var = "gene") %>%
    left_join(y = unique(annotations[, c("gene_name", "description")]),
              by = c("gene" = "gene_name")) %>%
    cbind(cluster_id = cluster, .)
}

# Iterate function across desired clusters - in this case, 22 clusters at 0.4 resolution
Idents(object = seurat_integrated) <- "integrated_snn_res.0.4"
conserved_markers <- map_dfr(0:22, get_conserved)
### In some cases you will have clusters that do not have enough cells for a
### particular group - and your function will fail. For these clusters you will
### need to use FindAllMarkers().


# Normalize & scale ADT data
seurat_integrated <- NormalizeData(seurat_integrated, assay = "ADT", normalization.method = "CLR")
seurat_integrated <- ScaleData(seurat_integrated, assay = "ADT")

# Now, we will visualize IgD levels for RNA and protein By setting the default assay, we can
# visualize one or the other
DefaultAssay(seurat_integrated) <- "ADT"
p1 <- FeaturePlot(seurat_integrated, "IgD", cols = c("lightgrey", "darkgreen")) + ggtitle("IgD protein")


FeaturePlot(seurat_integrated, 
            reduction = "umap", 
            features = c("IgD", "CD10"), 
            order = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)

DefaultAssay(seurat_integrated) <- "RNA"
p2 <- FeaturePlot(seurat_integrated, "MME") + ggtitle("CD10 RNA")
