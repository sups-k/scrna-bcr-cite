#Load libraries
library(Seurat)
library(tidyverse)
library(RCurl)
library(cowplot)
library(future)

# Set up multi-core
plan("multisession")

# Set working directory
setwd("/home/suk571/Jocelyn_scRNA/KW9079_Farmer/210624_A00873_0376_BHC3KTDSX2_KW9079/cellranger-6.0.1/GRCh38/")

# Enable storage of large files in R memory
options(future.globals.maxSize = 40000 * 1024^2)

# Read Seurat object
seurat_integrated <- readRDS(file = "results/output_rds/pc10_seurat_integrated.rds")

## For differential expression between healthy and patient, refer bottom of code

# Normalize CITE-seq data
seurat_integrated <- NormalizeData(seurat_integrated, assay = "ADT", normalization.method = "CLR")

# Choose resolution 0.8
Idents(object = seurat_integrated) <- "integrated_snn_res.0.8"

# Visualize 2 markers at the same time
FeaturePlot(seurat_integrated, features = c("CR2", "ITGAX"), blend = TRUE, split.by = "sample", min.cutoff = "q10", max.cutoff = "q90")

# Violin plot
VlnPlot(object = seurat_integrated, assay = "ADT", features = c("IgD", "CD73", "CD45RB", "CD10"))

# Dot plot
DotPlot(seurat_integrated, features = c("IgD", "CD73", "CD45RB", "CD10", "CR2", "ITGAX"), split.by = "sample") + RotatedAxis()

# Determine differentiating markers for CD4+ T cell
b_cells <- FindMarkers(seurat_integrated, ident.1 = 5, ident.2 = c(6,17,18,20,26,27))
b_cells <- b_cells %>%
  dplyr::arrange(p_val_adj)

# Rename all clusters
seurat_integrated <- RenameIdents(object = seurat_integrated, 
                                  "0" = "FO",
                                  "1" = "Naive",
                                  "2" = "T3a",
                                  "3" = "FO",
                                  "4" = "IGKV3-20+",
                                  "5" = "T1/T2",
                                  "6" = "IGHV3-23+ TIGIT+ CD27- CCL4+",
                                  "7" = "MZP",
                                  "8" = "FO",
                                  "9" = "FO: IL4R+",
                                  "10" = "T1/T2",
                                  "11" = "Naive: IGHV4-34+ SH2D2A+",
                                  "12" = "Activated Naive",
                                  "13" = "Activated Naive",
                                  "14" = "IgA1+ IgG2+ IGKV3-20+",
                                  "15" = "IgA1+ IgG2+ IGHV3-30+",
                                  "16" = "IgA1+ IgG2+ IGHV3-23+",
                                  "17" = "T1/T2", 
                                  "18" = "Naive: IGHV4-34+ IGKV3-20+ CD27-", 
                                  "19" = "Plasmablasts", 
                                  "20" = "T1/T2",
                                  "21" = "Activated Naive: IGHV4-34+",
                                  "22" = "T3a",
                                  "23" = "T3b",
                                  "24" = "T3a",
                                  "25" = "T3a",
                                  "26" = "T3b",
                                  "27" = "T3b")

# Plot the UMAP of the labelled object
pal <- c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C","#FDBF6F","#FF7F00","#CAB2D6","#6A3D9A","#A6761D","#B15928","black")
levels(seurat_subset_labeled) <- c("T1/T2", "T3a", "T3b", "Naive", "Follicular", "MZP", "Activated Naive", "CD11c+ Double Negative", "CD45RB+ Double Negative", "Class Switched", "Plasma Cells", "CD10+ Switched Memory", "4")
DimPlot(object = seurat_subset_labeled, reduction = "umap", label = FALSE, cols = pal, split.by = "sample")
# DimPlot(seurat_integrated, reduction = "umap", label = TRUE, label.size = 3, repel = TRUE, split.by = "sample")

# Save final R object
write_rds(seurat_integrated, file = "results/output_rds/seurat_labelled.rds")


# Show only some clusters
# seurat_subset_labeled <- subset(seurat_integrated, idents = "Activated Naive")
# Remove some clusters
seurat_subset_labeled <- subset(seurat_integrated, idents = c("Naive",
                                                              "Naive: IGHV4-34+ SH2D2A+",
                                                              "Naive: IGHV4-34+ IGKV3-20+ CD27-",
                                                              "IgA1+ IgG2+ IGKV3-20+",
                                                              "IgA1+ IgG2+ IGHV3-30+",
                                                              "IgA1+ IgG2+ IGHV3-23+",
                                                              "IGKV3-20+",
                                                              "IGHV3-23+ TIGIT+ CD27- CCL4+"),
                                invert = TRUE)

# Re-visualize the clusters
DimPlot(object = seurat_subset_labeled, reduction = "umap", label = TRUE, label.size = 3, repel = TRUE, cols = pal, split.by = "sample")

## For differential expression between healthy and patient
seurat_integrated$celltype.condition <- paste(Idents(seurat_integrated), seurat_integrated$sample, sep = "_")
seurat_integrated$celltype <- Idents(seurat_integrated)
Idents(seurat_integrated) <- "celltype.condition"

# Patient vs healthy - all results are patient w.r.t. healthy
# FindMarkers() tests genes that are detected in a minimum of 10% cells in either of the 2 populations by default
# log2FC = 0.25 by default
# pct.1 = percentage of cells where gene is detected in ident.1
# avg_log2FC: log fold-chage of the average expression between the two groups. Positive values indicate that the gene is more highly expressed in the first group
cluster_13_DEG <- FindMarkers(seurat_integrated, ident.1 = "13_patient", ident.2 = "13_healthy", verbose = TRUE)
write.csv(cluster_13_DEG, file = "results/markers/cluster13_DEG.csv", row.names = FALSE)
# Sort in ascending order of adjusted p-value
cluster_13_DEG <- cluster_13_DEG %>% dplyr::arrange(p_val_adj) %>% 
  dplyr::filter(p_val_adj < 0.1)

# Filter by log2FC >= 1 & adjusted p-value < 10%
cluster_13_DEG <- cluster_13_DEG %>% dplyr::filter(avg_log2FC >= 1 & p_val_adj < 0.1)
cluster_13_DEG <- cluster_13_DEG %>% dplyr::filter(pct.1 > pct.2)


# Number of cells
n_cells <- FetchData(seurat_subset_labeled, 
                     vars = c("ident", "sample")) %>%
  dplyr::count(ident, sample) %>%
  tidyr::spread(ident, n)
rownames(n_cells) <- n_cells$sample
barplot(t(as.matrix(n_cells[,2:14])),
        legend.text = c("T1/T2", "T3a", "T3b", "Naive", "Follicular", "MZP", "Activated Naive", "CD11c+ Double Negative", "CD45RB+ Double Negative", "Class Switched", "Plasma Cells", "CD10+ Switched Memory", "Cluster 4"),
        args.legend = list(x = "topright", inset=c(-0.2,-0.1)),
        main = "Proportion of Cells in Clusters",
        xlab = "Condition",
        ylab = "Number of Cells",
        col = c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C","#FDBF6F","#FF7F00","#CAB2D6","#6A3D9A","#A6761D","#B15928","black")
)


# Add to BCR-seq file
meta <- FetchData(seurat_subset_labeled, vars = c("ident", "sample", "patient_ID"))
meta$cell <- rownames(meta)
db <- read.delim("~/immcantation-python3.8.12/imm-env/makedb_output/functional/genotype/create_germlines/VDJseq_mutation_quant.tab", stringsAsFactors=F)



