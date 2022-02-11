## Load libraries
library(Seurat)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)
library(future)
library(stringr.tools)

## Enable parallelisation
plan("multisession")
#plan("multicore")

## Set working directory
setwd("/home/suk571/Jocelyn_scRNA/KW9079_Farmer/210624_A00873_0376_BHC3KTDSX2_KW9079/cellranger-6.0.1/GRCh38/")

# Read the RData file
load(file = "merged_BCR_seurat.RData")

# Cell-level filtering
## Now that we have visualized the various metrics, we can decide on the
## thresholds to use to remove the low quality. Often the recommendations
## mentioned earlier are a rough guideline, but the specific experiment needs
## to inform the exact thresholds chosen. We will use the following thresholds:
  
## nUMI >= 500
## nGene >= 250
## log10GenesPerUMI > 0.8
## mitoRatio < 0.2

## According to the QC plots of the raw data, Sam Kazer believes that the
## cutoffs are "reasonable and liberal enough to allow B cells that may have
## lower RNA content  or higher mitochondria read alignments."
## "I would start to consider relaxing cutoffs when you see a large proportion of
## cells near your cutoff parameters."

# Filter out low quality cells using selected thresholds - these will change with experiment
filtered_seurat <- subset(x = merged_seurat, 
                          subset= (nUMI >= 500) & 
                            (nGene >= 250) & 
                            (log10GenesPerUMI > 0.80) & 
                            (mitoRatio < 0.20))


# Gene-level filtering
## Within our data we will have many genes with zero counts. These genes can
## dramatically reduce the average expression for a cell and so we will remove
## them from our data. We will start by identifying which genes have a zero
## count in each cell:

# Extract counts
counts <- GetAssayData(object = filtered_seurat, slot = "counts")

# Output a logical matrix specifying for each gene on whether or not there are more than zero counts per cell
nonzero <- counts > 0


## Now, we will perform some filtering by prevalence. If a gene is only expressed
## in a handful of cells, it is not particularly meaningful as it still brings
## down the averages for all other cells it is not expressed in. For our data we
## choose to **keep only genes which are expressed in 10 or more cells**.
## By using this filter, genes which have zero counts in all cells will
## effectively be removed.


# Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
keep_genes <- Matrix::rowSums(nonzero) >= 10

# Only keeping those genes expressed in more than 10 cells
filtered_counts <- counts[keep_genes, ]


## Finally, take those filtered counts and create a new Seurat object for downstream analysis.

rm(merged_seurat)
gc()

# Reassign to filtered Seurat object
filtered_seurat_new <- filtered_seurat
filtered_seurat_new <- CreateSeuratObject(filtered_counts, meta.data = filtered_seurat@meta.data)
filtered_seurat_new[["ADT"]] <- filtered_seurat@assays$ADT

# Save filtered subset to new metadata
metadata_clean <- filtered_seurat_new@meta.data

# Remove old data
rm(counts, filtered_counts, nonzero, keep_genes, filtered_seurat)
gc()



## Re-assess QC Metrics

## ----------------Cell Counts --------------------

metadata_clean %>% 
  ggplot(aes(x=sample, fill=sample), size = 10) + 
  geom_bar(width = 0.20) + 
  ggtitle("NCells") +
  theme_bw() +
  theme(axis.text.x = element_text(size=rel(3), angle = 45, hjust = 1),
        axis.text.y = element_text(size=rel(4)),
        axis.title = element_text(size=rel(4)),
        plot.title = element_text(size=rel(4), hjust = 0.5, face = "bold"),
        legend.text = element_text(size=20),
        legend.title = element_text(size=20))


## --------------- UMI counts (transcripts) per cell --------------------

metadata_clean %>% 
  ggplot(aes(color=sample, x=nUMI, fill=sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  ylab("log10 cell density") +
  geom_vline(xintercept = 500) +
  theme_bw() +
  theme(axis.text.x = element_text(size=30),
        axis.text.y = element_text(size=40),
        axis.title = element_text(size=40),
        plot.title = element_text(size=40, hjust = 0.5, face = "bold"),
        legend.text = element_text(size=20),
        legend.title = element_text(size=20))

## -------------- Genes detected per cell --------------------------------


## As a histogram
metadata_clean %>% 
  ggplot(aes(color=sample, x=nGene, fill=sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  geom_vline(xintercept = 300) +
  theme_classic() +
  theme(axis.text.x = element_text(size=30),
        axis.text.y = element_text(size=40),
        axis.title = element_text(size=40),
        plot.title = element_text(size=40, hjust = 0.5, face = "bold"),
        legend.text = element_text(size=20),
        legend.title = element_text(size=20))

## As a bar plot
metadata_clean %>% 
  ggplot(aes(x=sample, y=log10(nGene), fill=sample)) + 
  geom_boxplot() + 
  ggtitle("NCells vs NGenes") +
  theme_bw() +
  theme(axis.text.x = element_text(size=rel(3)),
        axis.text.y = element_text(size=rel(4)),
        axis.title = element_text(size=rel(4)),
        plot.title = element_text(size=rel(4), hjust = 0.5, face = "bold"),
        legend.text = element_text(size=20),
        legend.title = element_text(size=20))

## -------------------- UMIs vs Genes Detected --------------------------------

metadata_clean %>% 
  ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250) +
  facet_wrap(~sample)+
  theme(axis.text.x = element_text(size=rel(3)),
        axis.text.y = element_text(size=rel(4)),
        axis.title = element_text(size=rel(4)),
        plot.title = element_text(size=rel(4), hjust = 0.5, face = "bold"),
        legend.text = element_text(size=20),
        legend.title = element_text(size=20))

## -------------- Mitochondrial counts ratio --------------------------------

metadata_clean %>% 
  ggplot(aes(color=sample, x=mitoRatio, fill=sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  geom_vline(xintercept = 0.2) +
  theme_bw() +
  theme(axis.text.x = element_text(size=30),
        axis.text.y = element_text(size=40),
        axis.title = element_text(size=40),
        plot.title = element_text(size=40, hjust = 0.5, face = "bold"),
        legend.title = element_text(size=25),
        legend.text = element_text(size = 25))



## ------------- Novelty/Complexity ---------------------------------

metadata_clean %>%
  ggplot(aes(x=log10GenesPerUMI, color = sample, fill=sample)) +
  geom_density(alpha = 0.2) +
  geom_vline(xintercept = 0.8) +
  theme_bw() +
  theme(axis.text.x = element_text(size=30),
        axis.text.y = element_text(size=40),
        axis.title = element_text(size=40),
        plot.title = element_text(size=40, hjust = 0.5, face = "bold"),
        legend.title = element_text(size=25),
        legend.text = element_text(size = 25))


# Create .RData object to load at any time
save(filtered_seurat_new, file="seurat_filtered_after_QC.RData")


### Run 2_filtering_v2.R before 3_normalization.R #####