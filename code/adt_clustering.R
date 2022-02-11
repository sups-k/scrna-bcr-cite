# Load libraries
library(Seurat)
library(tidyverse)
library(RCurl)
library(cowplot)

# Set up multi-core
future::plan("multicore")

# Set working directory
setwd("/home/suk571/Jocelyn_scRNA/KW9079_Farmer/210624_A00873_0376_BHC3KTDSX2_KW9079/cellranger-6.0.1/GRCh38/")

# Enable storage of large files in R memory
options(future.globals.maxSize = 400000 * 1024^2)

# Load integrated Seurat object
seurat_integrated <- readRDS(file = "results/output_rds/sauron_integrated_Jan23_pc43.rds")

Idents(object = seurat_integrated) <- "integrated_snn_res.0.4"

# Normalize & scale ADT data
seurat_integrated <- NormalizeData(seurat_integrated, assay = "ADT", normalization.method = "CLR")
seurat_integrated <- ScaleData(seurat_integrated, assay = "ADT")

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
adt.markers <- adt.markers[c(1,2,4,7,9,10,12,15,16,19,20,21),]
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

