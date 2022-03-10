# Load libraries
library(Seurat)
library(tidyverse)
library(AnnotationHub)
library(ensembldb)

future::plan("multisession")
# future::plan("multicore")

# Load the Seurat object
load(file = "/home/suk571/Jocelyn_scRNA/KW9079_Farmer/210624_A00873_0376_BHC3KTDSX2_KW9079/cellranger-6.0.1/GRCh38/results/output_rds/merged_seurat.RData")
merged_seurat <- subset(merged_seurat, cells=colnames(merged_seurat) [ Matrix::colSums(merged_seurat) > 0 ],
                        features=rownames(merged_seurat) [ Matrix::rowSums(merged_seurat) > 0 ] )

cat("The total dimensions of your dataset is: ",dim(merged_seurat),"\n")
# The total dimensions of the dataset is: 26122 x 75789
# Now 26010 x 75683 after removing T/NK & macrophages

# Calculate data diversity indices of gene expression
indexes <- t(apply(merged_seurat@assays[["RNA"]]@counts,2,function(x) {
  c(vegan::diversity(x,index = "simpson"),
    vegan::diversity(x,index = "invsimpson"),
    vegan::diversity(x,index = "shannon"),
    ineq::Gini(x)) }))
merged_seurat$simp_index <- indexes[,1]
merged_seurat$invsimp_index <- indexes[,2]
merged_seurat$shan_index <- indexes[,3]
merged_seurat$gini_index <- indexes[,4]
invisible(gc())

#############################################
### CALCULATE PERCENTAGE OF GENE FAMILIES ###
#############################################
cat("\nCalculating percentage of mitocondrial/ribosomal genes ...\n")
Gene.groups <- substring(rownames(x = merged_seurat@assays[["RNA"]]@counts),1,3)
seq_depth <- Matrix::colSums(merged_seurat@assays[["RNA"]]@counts)
temp <- rowsum(as.matrix(merged_seurat@assays[["RNA"]]@counts),Gene.groups)
perc <- sort(apply( t(temp) / seq_depth,2,median) ,decreasing = T)*100
tot <- sort(rowSums(temp)/sum(temp),decreasing = T)*100

#Compute the relative expression of each gene per cell
rel_expression <- Matrix::t( Matrix::t(merged_seurat@assays[["RNA"]]@counts) / (Matrix::colSums(merged_seurat@assays[["RNA"]]@counts)) ) * 100
most_expressed <- sort(apply(rel_expression,1,mean),T)[1:100] / ncol(merged_seurat)

png(filename = "~/Jocelyn_scRNA/KW9079_Farmer/210624_A00873_0376_BHC3KTDSX2_KW9079/cellranger-6.0.1/GRCh38/QC_raw_data/gene_family proportions.png",width = 600*3,height = 4*600,res = 150)
rafalib::mypar(4,1,mar=c(5,5,2,1))
boxplot( as.matrix(Matrix::t(rel_expression[names(most_expressed),])),cex=.1,outline=T,las=2,main="% total count per cell",col=scales::hue_pal()(100))
boxplot( (t(temp)/seq_depth) [,names(perc)[1:100]]*100,outline=T,las=2,main="% reads per cell",col=scales::hue_pal()(100))
boxplot(t(temp)[,names(perc)[1:100]], outline=T,las=2,main="reads per cell",col=scales::hue_pal()(100) )
barplot(tot[names(tot)[1:100]],las=2,xaxs="i",main="Total % reads (all cells)",col=scales::hue_pal()(100))
invisible(dev.off())

plot_gene_family = c("RPS", "RPL", "mito", "HB")
for(i in unique( c("rpl","rps","hb","mito",unlist(strsplit(casefold(plot_gene_family),","))))){
  cat(i,"\t")
  family.genes <- rownames(merged_seurat@assays[["RNA"]]@counts)[grep(pattern = paste0("^",ifelse(i=="mito","mt-",i)), x = casefold(rownames(merged_seurat@assays[["RNA"]]@counts)), value = F)]
  if(length(family.genes)>1){merged_seurat <- PercentageFeatureSet(merged_seurat,features = family.genes,assay = "RNA",col.name = paste0("perc_",ifelse(i=="mt-","mito",i)) )}
}

rm("temp","perc","tot","Gene.groups","i","indexes")
invisible(gc())
#---------



############################################
### CALCULATING GENE BIOTYPE PERCENTAGES ###
############################################


cat("\nCalculating gene biotype percentages ...\n")
# Connect to AnnotationHub
ah <- AnnotationHub()

# Access the Ensembl database for organism
ahDb <- query(ah, pattern = c("Homo sapiens", "EnsDb"), ignore.case = TRUE)
# Acquire the latest annotation files
id <- ahDb %>% mcols() %>% rownames() %>% tail(n = 1)
# Download the appropriate Ensembldb database
edb <- ah[[id]]
annotations <- genes(edb, return.type = "data.frame")
# Select annotations of interest
annotations <- annotations %>%
  dplyr::select(gene_id, gene_name, gene_biotype, seq_name, description, entrezid)

output_path <- "~/Jocelyn_scRNA/KW9079_Farmer/210624_A00873_0376_BHC3KTDSX2_KW9079/cellranger-6.0.1/GRCh38/QC_raw_data"

item <- annotations[match(rownames(merged_seurat@assays[["RNA"]]@counts) , annotations[,2]), "gene_biotype"]
item[is.na(item)] <- "unknown"

png(filename = paste0(output_path,"/gene_biotype_proportions.png"),width = 600*3,height = 600,res = 150)
rafalib::mypar(1,3,mar=c(4,2,2,1))
pie(sort(table(item),decreasing = T), clockwise = T,col = scales::hue_pal()(length(unique(item))))
title("before filtering")
par(mar=c(10,2,2,1))
barplot(sort(table(item),decreasing = T),las=2,xaxs="i",main="Total reads (all cells)",col=scales::hue_pal()(100))

temp <- rowsum(as.matrix(merged_seurat@assays[["RNA"]]@counts),group=item)
o <- order(apply(temp,1,median),decreasing = T)
boxplot( (t(temp)/Matrix::colSums(merged_seurat@assays[["RNA"]]@counts))[,o]*100,outline=F,las=2,main="% reads per cell",col=scales::hue_pal()(100))
invisible(dev.off())

aaa <- setNames(as.data.frame(((t(temp)/Matrix::colSums(merged_seurat@assays[["RNA"]]@counts))[,o]*100)[,names(sort(table(item),decreasing = T))]),paste0("perc_",names(sort(table(item),decreasing = T))))
merged_seurat@meta.data <- merged_seurat@meta.data[,!(colnames(merged_seurat@meta.data) %in% colnames(aaa))]
merged_seurat@meta.data <- cbind(merged_seurat@meta.data,aaa)

#---------

##########################
### CELL CYCLE SCORING ###
##########################
cat("\nLog Normalizing counts for cell cycle scoring...\n")
merged_seurat <- NormalizeData(object = merged_seurat, scale.factor = 1000)

cat("\nPredicting cell cycle scores with Seurat ...\n")
s.genes <- Seurat::cc.genes$s.genes
g2m.genes <- Seurat::cc.genes$g2m.genes

merged_seurat <- CellCycleScoring(object = merged_seurat, s.features = s.genes, g2m.features = g2m.genes)
merged_seurat$G1.Score <- 1 - ( merged_seurat$S.Score + merged_seurat$G2M.Score )
merged_seurat$CC.Diff <- merged_seurat$S.Score - merged_seurat$G2M.Score
#---------

###############
### PLOT QC ###
###############
cat("\nPlotting QC metrics ...\n")

feats <- colnames(merged_seurat@meta.data) [ grepl("nFeature|nCount|_index|[.]Score",colnames(merged_seurat@meta.data) ) ]
feats <- c(feats,"perc_mito" ,"perc_rps","perc_rpl","perc_hb", "perc_protein_coding" ,"perc_lincRNA","perc_snRNA","perc_miRNA","perc_processed_pseudogene",
             "perc_unknown","perc_Chr_1","perc_Chr_X","perc_Chr_Y","perc_Chr_MT")
feats <- feats[feats %in% colnames(merged_seurat@meta.data)]

png(filename = paste0(output_path,"/QC_patient_ALL.png"),width = 1200*(length(unique(merged_seurat@meta.data[,"patient_ID"]))/2+1),height = 700*ceiling(length(feats)/5),res = 200)
print(VlnPlot(object = merged_seurat, features  = feats, ncol = 5, group.by = "patient_ID", pt.size = .1, assay = "RNA"))
invisible(dev.off())
#---------



######################
### CELL FILTERING ###
######################
cat("\nFiltering low quality cells ...\n")
NF <-  merged_seurat@meta.data [ grepl("nFeature",colnames(merged_seurat@meta.data)) ][,1]
NC <-  merged_seurat@meta.data [ grepl("nCount",colnames(merged_seurat@meta.data)) ][,1]

mito_range <- c(0,25)
ribo_range <- c(0,25)


Ts <- data.frame(
  MitoT = between(merged_seurat$perc_mito, mito_range[1], mito_range[2]),
  RpsT = between(merged_seurat$perc_rps, ribo_range[1], ribo_range[2]),
  RplT = between(merged_seurat$perc_rpl, ribo_range[1], ribo_range[2]),
  nUMIT = between(NF,quantile(NF,probs = c(0.005)),quantile(NF,probs = c(0.995))),
  nCountT = between(NC,quantile(NC,probs = c(0.005)),quantile(NC,probs = c(0.995))),
  GiniT = between(merged_seurat$gini_index,0.8,1),
  SimpT = between(merged_seurat$simp_index,0.8,1),
  protein_codingT = between(merged_seurat$perc_protein_coding,50,100),
  row.names = rownames(merged_seurat@meta.data) )
print(head(Ts,90))

dim(merged_seurat)
cell_use <- rownames(Ts)[ rowSums(!Ts) == 0 ]
merged_seurat$filtered <- (rowSums(Ts) == 0)
#---------

rm(aaa, ah, ahDb, edb, rel_expression, temp, Ts, s.genes, ribo_range, seq_depth,
   o, NC, NF, item, g2m.genes, id, family.genes, mito_range, most_expressed)
invisible(gc())

########################################
### IDENTIFY REQUESTED GENES TO KEEP ###
########################################
# Ighd, Ighm, Ighg1, Ighg2c, Ighg2b, Ighg3, Igha, Ighe

genes_keep <- c("IGHD", "IGHM", "IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHA1", "IGHE", "IGHA2")
genes_notfound <- setdiff(genes_keep, merged_seurat@assays[["RNA"]]@counts@Dimnames[[1]])
genes_keep <- merged_seurat@assays[["RNA"]]@counts@Dimnames[[1]][merged_seurat@assays[["RNA"]]@counts@Dimnames[[1]] %in% genes_keep]
cat("\nThe following genes will NOT be removed from the data:\n")
cat(genes_keep, "\n")
if(length(genes_notfound) > 0){
  cat("\nWARNING: The following requested genes were NOT FOUND in the data:\n")
  cat(genes_notfound, "\n")
}


###########################################
### SELECTING ONLY PROTEIN-CODING GENES ###
###########################################

cat("\nSelect only the protein-coding genes ...\n")
sel <- annotations[match(merged_seurat@assays[["RNA"]]@counts@Dimnames[[1]], annotations$gene_name),]
sel <- sel %>% dplyr::filter(gene_biotype == "protein_coding")
genes_use <- union(as.character(na.omit(sel$gene_name)), genes_keep)
merged_seurat@assays[["RNA"]]@counts <- merged_seurat@assays[["RNA"]]@counts[genes_use,]
merged_seurat@assays[["RNA"]]@data <- merged_seurat@assays[["RNA"]]@data[genes_use,]

#---------


#############################################
### REMOVING SELECTED GENES FROM THE DATA ###
#############################################
cat("\nRemoving mitochondrial genes from the data ...\n")
genes_use <- rownames(merged_seurat@assays[["RNA"]]@counts)[!grepl(gsub(",","|", sub('mito','mt-',casefold("mito")) ), casefold(rownames(merged_seurat@assays[["RNA"]]@counts)))]
genes_use <- union(genes_use, genes_keep)
merged_seurat@assays[["RNA"]]@counts <- merged_seurat@assays[["RNA"]]@counts[genes_use,]
merged_seurat@assays[["RNA"]]@data <- merged_seurat@assays[["RNA"]]@data[genes_use,]
a <- merged_seurat
load(file = "/home/suk571/Jocelyn_scRNA/KW9079_Farmer/210624_A00873_0376_BHC3KTDSX2_KW9079/cellranger-6.0.1/GRCh38/results/output_rds/merged_seurat.RData")
a[["ADT"]] <- merged_seurat@assays$ADT
#---------


##########################
### FILTERING & SAVING COUNTS ###
##########################
DATA <- CreateSeuratObject(counts = a@assays[["RNA"]]@counts[, cell_use], assay = "RNA", meta.data = a@meta.data[cell_use,], min.cells = 5, min.features = 200)
b <- CreateSeuratObject(counts = a@assays[["ADT"]]@counts[, colnames(DATA)], assay = "ADT", meta.data = a@meta.data[colnames(DATA),])

DATA[["ADT"]] <- b@assays$ADT

rm(a, b, merged_seurat, annotations)
rm(sel, cell_use, genes_keep, genes_notfound, genes_use)
invisible(gc())

# Filtering low-quality cells
filtered_seurat <- subset(x = DATA, subset= log10GenesPerUMI > 0.80)
filtered_seurat <- subset(x = filtered_seurat, subset= nCount_RNA >= 500)

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
# Reassign to filtered Seurat object
filtered_seurat_new <- filtered_seurat
filtered_seurat_new <- CreateSeuratObject(filtered_counts, meta.data = filtered_seurat@meta.data)
filtered_seurat_new[["ADT"]] <- filtered_seurat@assays$ADT

saveRDS(filtered_seurat_new, file = "~/Jocelyn_scRNA/KW9079_Farmer/210624_A00873_0376_BHC3KTDSX2_KW9079/cellranger-6.0.1/GRCh38/results/output_rds/sauron_filtered.rds")

#---------

rm(counts, DATA, filtered_counts, filtered_seurat, nonzero, keep_genes)
invisible(gc())

###############
### PLOT QC ###
###############
cat("\nPlotting QC metrics ...\n")
## Normalization for visualizing
# Change scale factor from default to 1000 (total counts per cell) for cell-level normalization
filtered_seurat_new <- NormalizeData(object = filtered_seurat_new, scale.factor = 1000, normalization.method = "LogNormalize", assay = "RNA")
# Normalizing across cells (not features)
# filtered_seurat_new <- NormalizeData(object = filtered_seurat_new, normalization.method = "CLR", margin = 2, assay = "ADT")

metadata_clean <- filtered_seurat_new@meta.data


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

## As a box plot
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
  geom_vline(xintercept = 0.25) +
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



png(filename = paste0(output_path,"/QC_patient_FILTERED.png"),width = 1200*(length(unique(filtered_seurat_new@meta.data[,"patient_ID"]))/2+1),height = 700*ceiling(length(feats)/5),res = 200)
print(VlnPlot(object = filtered_seurat_new, features  = feats, ncol = 5, group.by = "patient_ID", pt.size = .1, assay = "RNA"))
invisible(dev.off())

rm(filtered_seurat_new, metadata_clean, feats, output_path)
invisible(gc())
