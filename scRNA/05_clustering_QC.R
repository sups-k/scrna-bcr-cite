## Run in RStudio


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
options(future.globals.maxSize = 50000 * 1024^2)

# Load integrated Seurat object
seurat_integrated <- readRDS(file = "results/output_rds/sauron_integrated_Jan23_pc43.rds")

# Set resolution based on clustree output from previous code
Idents(object = seurat_integrated) <- "integrated_snn_res.0.4"


########################
## QC FOR CLUSTERING ###
########################

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


# Normalize & scale ADT data
DefaultAssay(seurat_integrated) <- "ADT"
seurat_integrated <- NormalizeData(seurat_integrated, assay = "ADT", normalization.method = "CLR")
seurat_integrated <- ScaleData(seurat_integrated, assay = "ADT")

# Now, we will visualize IgD levels for RNA and protein By setting the default assay, we can
# visualize one or the other

# FeaturePlot(seurat_integrated, 
#            reduction = "umap", 
#            features = c("IgD", "CD10"), 
#            order = TRUE,
#            min.cutoff = 'q10', 
#            label = TRUE)

# p1 <- FeaturePlot(seurat_integrated, "IgD", cols = c("lightgrey", "darkgreen")) + ggtitle("IgD protein")
# DefaultAssay(seurat_integrated) <- "RNA"
# p2 <- FeaturePlot(seurat_integrated, "MME") + ggtitle("CD10 RNA")



# View all CITE-seq markers
FeaturePlot(seurat_integrated,
            reduction = "umap",
            features = c("adt_CD10", "adt_IgD", "adt_CD73", "adt_CD45RB"),
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE)
RidgePlot(seurat_integrated, features = c("adt_CD10", "adt_IgD", "adt_CD73", "adt_CD45RB"), ncol = 2)

# Downsample the clusters to a maximum of 300 cells each (makes the heatmap easier to see for
# small clusters)
cbmc.small <- subset(seurat_integrated, downsample = 300)

# Find protein markers for all clusters, and draw a heatmap
adt.markers <- FindAllMarkers(cbmc.small, assay = "ADT", only.pos = TRUE)
adt.markers <- adt.markers[c(1,2,3,6,10,12,14,16,17,20,21,22),]
DoHeatmap(cbmc.small, features = unique(adt.markers$gene), assay = "ADT", angle = 90) + NoLegend()

# Since we only have 4 markers, instead of doing PCA, we'll just use a standard euclidean
# distance matrix here.
DefaultAssay(seurat_integrated) <- "ADT"
adt.data <- GetAssayData(seurat_integrated, slot = "data")
adt.data <- adt.data[c(1,2,3,4),]
adt.dist <- dist(t(adt.data))
rm(adt.data)
gc()
# Before we recluster the data on ADT levels, we'll stash the RNA cluster IDs for later
seurat_integrated[["rnaClusterID"]] <- Idents(seurat_integrated)

# Now, we rerun UMAP using our distance matrix defined only on ADT (protein) levels.
seurat_integrated[["umap_adt"]] <- RunUMAP(adt.dist, assay = "ADT", reduction.key = "adtUMAP_")
seurat_integrated[["adt_snn"]] <- FindNeighbors(adt.dist)$snn

# saveRDS(seurat_integrated, "results/output_rds/adt_clustering.rds")
seurat_integrated <- FindClusters(seurat_integrated, resolution = 0.4, graph.name = "adt_snn")

# We can compare the RNA and protein clustering, and use this to annotate the protein clustering
# (we could also of course use FindMarkers)
clustering.table <- table(Idents(seurat_integrated), seurat_integrated$rnaClusterID)
clustering.table

# 0 = 0,1,3
# 1 = 0,1
# 2 = 0,1
# 3 = 0,1,3
# 4 = 0,1
# 5 = 4
# 6 = 0,1
# 7 = 2
# 8 = 2
# 9 = 1,2
# 10 = 3,5
# 11 = 3,5
# 12 = 4
# 13 = 4
# 14 = 4

new.cluster.ids <- c("0,1,3", "0,1", "0,1", "0,1,3", "0,1", "4", "0,1", "2", "2", 
                     "1,2", "3,5", "3,5", "4", "4", "4")
names(new.cluster.ids) <- levels(seurat_integrated)
seurat_integrated <- RenameIdents(seurat_integrated, new.cluster.ids)

umap_rnaClusters <- DimPlot(seurat_integrated, reduction = "umap_adt", group.by = "rnaClusterID") + NoLegend()
umap_rnaClusters <- umap_rnaClusters + ggtitle("Clustering based on scRNA-seq") + theme(plot.title = element_text(hjust = 0.5))
umap_rnaClusters <- LabelClusters(plot = umap_rnaClusters, id = "rnaClusterID", size = 4)

umap_adtClusters <- DimPlot(seurat_integrated, reduction = "umap_adt", pt.size = 0.5) + NoLegend()
umap_adtClusters <- umap_adtClusters + ggtitle("Clustering based on ADT signal") + theme(plot.title = element_text(hjust = 0.5))
umap_adtClusters <- LabelClusters(plot = umap_adtClusters, id = "ident", size = 4)

patchwork::wrap_plots(list(umap_rnaClusters, umap_adtClusters), ncol = 2)
