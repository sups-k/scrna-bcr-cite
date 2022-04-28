library(Seurat)
library(tidyverse)
library(cowplot)
library(SingleCellExperiment)
library(scater)
library(Matrix.utils)
library(edgeR)
library(dplyr)
library(magrittr)
library(Matrix)
library(purrr)
library(reshape2)
library(S4Vectors)
library(tibble)
library(pheatmap)
library(apeglm)
library(png)
library(DESeq2)
library(RColorBrewer)

# Set working directory
setwd("/home/suk571/Jocelyn_scRNA/KW9079_Farmer/210624_A00873_0376_BHC3KTDSX2_KW9079/cellranger-6.0.1/GRCh38/")

# Enable storage of large files in R memory
options(future.globals.maxSize = 40000 * 1024^2)

# Bring in Seurat object
seurat <- readRDS(file = "results/output_rds/10_pc43_11March.rds")

# Extract raw counts and metadata to create SingleCellExperiment object
counts <- seurat@assays$RNA@counts 

metadata <- seurat@meta.data

Idents(object = seurat) <- "integrated_snn_res.0.4"

# Set up metadata as desired for aggregation and DE analysis
metadata$cluster_id <- factor(seurat@active.ident)
metadata <- metadata[,-c(53,54,55,56,57,58,59,60,61,62,63)] # removing other clustering
metadata$patient_ID <- as.factor(metadata$patient_ID)

# Create single cell experiment object
sce <- SingleCellExperiment(assays = list(counts = counts), 
                            colData = metadata)
rm(seurat)
# Identify groups for aggregation of counts
groups <- colData(sce)[, c("cluster_id", "patient_ID")]

# Named vector of cluster names
kids <- purrr::set_names(levels(sce$cluster_id))
kids

# Total number of clusters
nk <- length(kids)
nk

# Named vector of sample names
sids <- purrr::set_names(levels(sce$patient_ID))

# Total number of samples 
ns <- length(sids)
ns

# Generate sample level metadata

## Determine the number of cells per sample
table(sce$patient_ID)

## Turn named vector into a numeric vector of number of cells per sample
n_cells <- as.numeric(table(sce$patient_ID))

## Determine how to reoder the samples (rows) of the metadata to match the order of sample names in sids vector
m <- match(sids, sce$patient_ID)

## Create the sample level metadata by combining the reordered metadata with the number of cells corresponding to each sample.
ei <- data.frame(colData(sce)[m, ], 
                 n_cells, row.names = NULL) %>% 
  dplyr::select(-"cluster_id")
ei <- ei %>% dplyr::select(sample, patient_ID, n_cells)

# Aggregate the counts per patient_ID and cluster_id

# Subset metadata to only include the cluster and sample IDs to aggregate across
groups <- colData(sce)[, c("cluster_id", "patient_ID")]




# Aggregate across cluster-sample groups
pb <- aggregate.Matrix(t(counts(sce)), 
                       groupings = groups, fun = "sum") 

class(pb)

dim(pb)

pb[1:6, 1:6]

# Not every cluster is present in all samples; create a vector that represents how to split samples
splitf <- sapply(stringr::str_split(rownames(pb), pattern = "_", n = 2), `[`, 1)
# Turn into a list and split the list into components for each cluster and transform, so rows are genes and columns are samples and make rownames as the sample IDs
pb <- split.data.frame(pb, 
                       factor(splitf)) %>%
  lapply(function(u) 
    set_colnames(t(u), 
                 stringr::str_extract(rownames(u), "(?<=_)[:alnum:]+")))

class(pb)

# Explore the different components of list
str(pb)

# Print out the table of cells in each cluster-sample group
options(width = 100)
table(sce$cluster_id, sce$patient_ID)


## STARTING DESEQ2 ANALYSIS ###

# Get sample names for each of the cell type clusters

# prep. data.frame for plotting
get_sample_ids <- function(x){
  pb[[x]] %>%
    colnames()
}

de_samples <- map(1:length(kids), get_sample_ids) %>%
  unlist()

# Get cluster IDs for each of the samples

samples_list <- map(1:length(kids), get_sample_ids)

get_cluster_ids <- function(x){
  rep(names(pb)[x], 
      each = length(samples_list[[x]]))
}

de_cluster_ids <- map(1:length(kids), get_cluster_ids) %>%
  unlist()
# Create a data frame with the sample IDs, cluster IDs and condition

gg_df <- data.frame(cluster_id = de_cluster_ids,
                    patient_ID = de_samples)

gg_df <- dplyr::left_join(gg_df, ei[, c("patient_ID", "sample")]) 


metadata2 <- gg_df %>%
  dplyr::select(cluster_id, patient_ID, sample) 
metadata2$cluster_id <- as.factor(metadata2$cluster_id)
View(metadata2)

# Generate vector of cluster IDs
clusters <- levels(metadata2$cluster_id)
clusters

### RUNNING ANALYSIS ON CLUSTER 0 #####

# Subset the metadata to only cluster 0
cluster_metadata <- metadata2[which(metadata2$cluster_id == clusters[1]), ]
head(cluster_metadata)

# Assign the rownames of the metadata to be the sample IDs
rownames(cluster_metadata) <- cluster_metadata$patient_ID
head(cluster_metadata)

# Subset the counts to only cluster 0
counts <- pb[[clusters[1]]]

cluster_counts <- data.frame(counts[, which(colnames(counts) %in% rownames(cluster_metadata))])

# Check that all of the row names of the metadata are the same and in the same order as the column names of the counts in order to use as input to DESeq2
all(rownames(cluster_metadata) == colnames(cluster_counts))        

# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(cluster_counts, 
                              colData = cluster_metadata, 
                              design = ~ sample)
## Running QC

# Transform counts for data visualization on PCA
rld <- rlog(dds, blind=TRUE)

# QC1: Plot PCA
DESeq2::plotPCA(rld, intgroup = "sample") + geom_text(label = rld@colData@listData[["patient_ID"]])


# Extract the rlog matrix from the object and compute pairwise correlation values for heatmap
rld_mat <- assay(rld)
rld_cor <- cor(rld_mat)

# QC2: Plot heatmap
pheatmap(rld_cor, annotation = cluster_metadata[, c("sample"), drop=F])


## Running DESeq2
# Run DESeq2 differential expression analysis
dds <- DESeq(dds)

# We can check the fit of the model to our data by looking at the plot of dispersion estimates.
# Plot dispersion estimates
plotDispEsts(dds)
# The plot is encouraging, since we expect our dispersions to decrease with increasing mean and follow the line of best fit.

## Comparing patient to healthy
# Output results of Wald test for contrast for patient vs healthy
levels(as.factor(cluster_metadata$sample))[2]
levels(as.factor(cluster_metadata$sample))[1]

contrast <- c("sample", levels(as.factor(cluster_metadata$sample))[2], levels(as.factor(cluster_metadata$sample))[1])

resultsNames(dds)

res <- results(dds, 
               contrast = contrast,
               alpha = 0.05) #FDR cut-off

res <- lfcShrink(dds, 
                 coef =  2,
                 res=res)
# coef = 2 because we want patient vs healthy - index no. of resultsNames(dds)

## Table of results for all genes
# Turn the results object into a tibble for use with tidyverse functions
res_tbl <- res %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>%
  as_tibble()

# Check results output
View(res_tbl)

# Write all results to file
write.csv(res_tbl,
          paste0("results/cluster", clusters[1], "_", levels(as.factor(cluster_metadata$sample))[2], "_vs_", levels(as.factor(cluster_metadata$sample))[1], "_all_genes.csv"),
          quote = FALSE, 
          row.names = FALSE)

# Filter significant genes
# Set thresholds
padj_cutoff <- 0.05

# Subset the significant results
sig_res <- dplyr::filter(res_tbl, padj < padj_cutoff) %>%
  dplyr::arrange(padj)

# Check significant genes output
View(sig_res)

# Write significant results to file
write.csv(sig_res,
          paste0("results/cluster", clusters[1], "_", levels(as.factor(cluster_metadata$sample))[2], "_vs_", levels(as.factor(cluster_metadata$sample))[1], "_sig_genes.csv"),
          quote = FALSE, 
          row.names = FALSE)




## Scatterplot of normalized expression of top 20 most significant genes

## ggplot of top genes
normalized_counts <- counts(dds, 
                            normalized = TRUE)

## Order results by padj values
top20_sig_genes <- sig_res %>%
  dplyr::arrange(padj) %>%
  dplyr::pull(gene) %>%
  head(n=20)


top20_sig_norm <- data.frame(normalized_counts) %>%
  rownames_to_column(var = "gene") %>%
  dplyr::filter(gene %in% top20_sig_genes)

gathered_top20_sig <- top20_sig_norm %>%
  gather(colnames(top20_sig_norm)[2:length(colnames(top20_sig_norm))], key = "samplename", value = "normalized_counts")

gathered_top20_sig <- dplyr::inner_join(ei[, c("patient_ID", "sample" )], gathered_top20_sig, by = c("patient_ID" = "samplename"))

## plot using ggplot2
ggplot(gathered_top20_sig) +
  geom_point(aes(x = gene, 
                 y = normalized_counts, 
                 color = sample), 
             position=position_jitter(w=0.1,h=0)) +
  scale_y_log10() +
  xlab("Genes") +
  ylab("log10 Normalized Counts") +
  ggtitle("Top 20 Significant DE Genes") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  theme(plot.title = element_text(hjust = 0.5))


## Heatmap of all significant genes

# Extract normalized counts for only the significant genes
sig_norm <- data.frame(normalized_counts) %>%
  rownames_to_column(var = "gene") %>%
  dplyr::filter(gene %in% sig_res$gene)

# Set a color palette
heat_colors <- rev(brewer.pal(6, "RdYlBu"))

# Run pheatmap using the metadata data frame for the annotation
pheatmap(sig_norm[ , 2:length(colnames(sig_norm))], 
         color = heat_colors, 
         cluster_rows = T, 
         show_rownames = F,
         annotation = cluster_metadata[, c("sample", "cluster_id")], 
         border_color = NA, 
         fontsize = 10, 
         scale = "row", 
         fontsize_row = 10, 
         height = 20)


## Volcano plot of results

## Obtain logical vector where TRUE values denote padj values < 0.05 and fold change > 1.5 in either direction
res_table_thres <- res_tbl %>% 
  dplyr::mutate(threshold = padj < 0.05 & abs(log2FoldChange) >= 1.5)

## Volcano plot
ggplot(res_table_thres) +
  geom_point(aes(x = log2FoldChange, y = -log10(padj), colour = threshold)) +
  ggtitle("Volcano plot of patient cluster 0 cells relative to healthy") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  scale_y_continuous(limits = c(0,50)) +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))           


