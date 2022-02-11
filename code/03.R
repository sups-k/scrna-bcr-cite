library(Seurat) 
library(dplyr)
library(scales)
library(RColorBrewer)
library(biomaRt)
library(igraph)
#library(sva)
library(rafalib)
library(parallel)
library(scran)
#library(scater)
library(dbscan)

future::plan("multisession")
options(future.globals.maxSize = 100000 * 1024^2)
out_path <- "~/Jocelyn_scRNA/KW9079_Farmer/210624_A00873_0376_BHC3KTDSX2_KW9079/cellranger-6.0.1/GRCh38"

#############################
### LOAD Seurat.v3 OBJECT ###
#############################
cat("\n### LOADING Seurat.v3 OBJECT ###\n")
DATA <- readRDS(file = paste0(out_path, "/results/output_rds/12Jan_seurat_object.rds"))
DATA@active.assay <- "RNA"
#---------


###################################
### SELECT CELLS FROM A CLUSTER ###
###################################
# No cluster selected. Use all cells
cells_use <- rownames(DATA@meta.data)
#---------


###################
### Running PCA ###
###################
cat("\n### Running PCA ###\n")
if(!dir.exists(paste0(out_path,"/results/pca_plots"))){dir.create(paste0(out_path,"/results/pca_plots"),recursive = T)}
DATA <- RunPCA(DATA, do.print = F, assay = "RNA",npcs = 100)
write.csv(DATA@reductions$pca@cell.embeddings, paste0(out_path,"/results/pca_plots/PCA_coordinates.csv"))
write.csv(DATA@reductions$pca@feature.loadings, paste0(out_path,"/results/pca_plots/PCA_feature_loadings.csv"))

ggplot2::ggsave(PCHeatmap(DATA,ncol=5,dims=1:10),filename = paste0("PCA_heatmap.png"), path = paste0(out_path,"/results/pca_plots"), dpi = 300,units = "mm",width = 150*5,height = 150*4.5,limitsize = FALSE)

var_expl <- (DATA@reductions$pca@stdev^2)/sum(DATA@reductions$pca@stdev^2)
top_PCs <- 50

png(filename = paste0(out_path,"/results/pca_plots/Variance_explained_PC.png"),width = 1500,height =1200,res = 200)
plot( var_expl*100,yaxs="i",bg="grey",pch=21,type="l",ylab="% Variance",xlab="PCs",main="all cells",las=1,ylim=c(0,1.2*max(var_expl*100)))
points( var_expl*100,bg=c(rep("orange",top_PCs),rep("grey",100 - top_PCs)),pch=21)
invisible(dev.off())

# Identify number of significant PCs
## Step 1: Determine which PC exhibits cumulative percent greater than 90% and
## % variation associated with the PC as less than 5
pct <- DATA[["pca"]]@stdev / sum(DATA[["pca"]]@stdev) * 100
cumu <- cumsum(pct)
co1 <- which(cumu > 90 & pct < 5)[1]

## Step 2: Determine the point where the percent change in variation between the
## consecutive PCs is less than 0.1%
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1

## Step 3: Use the minimum of the 2 metrics calculated above
top_PCs <- min(co1, co2)

## Step 4: Remove unnecessary variables
rm(co1, co2, cumu, pct)

# Perform PCA again with the new number of principal components
DATA <- RunPCA(object = DATA, npcs = top_PCs, verbose = TRUE, assay = "RNA")

#---------


#################################################
### Parse dimensionality reduction parameters ###
#################################################
dim_reduct_params <- c('umap, n.neighbors=50, min.dist=0.1, spread=5, repulsion.strength=0.5, n.epochs=500, learning.rate=0.5, negative.sample.rate=7, metric="euclidean", seed.use=42;',
                       'umap10, n.neighbors=50, min.dist=0.1, spread=5, repulsion.strength=0.5, n.epochs=500, learning.rate=0.5, negative.sample.rate=7, metric="euclidean", seed.use=42')
param_set <- paste0(sub(",", "=list(", casefold(unlist(strsplit(dim_reduct_params, ";")))), ")", collapse=",")
dim_reduct_params <- eval(parse(text=paste0("list(", param_set, ")")))
#---------


####################
### Running UMAP ###
####################

if(!dir.exists(paste0(out_path,"/results/umap_plots"))){dir.create(paste0(out_path,"/results/umap_plots"),recursive = T)}

umap_params <- list(dims = 1:top_PCs,
                        n.neighbors = 10,
                        spread = 0.3,
                        repulsion.strength = 1,
                        learning.rate = 1,
                        min.dist= 0.001,
                        verbose = T,
                        num_threads=0,
                        n.epochs = 200,
                        metric = "euclidean",
                        negative.sample.rate = 5L,
                        seed.use = 42)
# Overwrite defaults with specified parameters
umap_params <- modifyList(umap_params, dim_reduct_params$umap)
    
ttt <- Sys.time()
DATA <- RunUMAP(object=DATA, reduction=casefold("mnn"), n.components=2,
                    dims=umap_params$dims, n.neighbors=umap_params$n.neighbors, learning.rate=umap_params$learning.rate, spread=umap_params$spread,
                    repulsion.strength=umap_params$repulsion.strength, min.dist=umap_params$min.dist, verbose=umap_params$verbose,
                    num_threads=umap_params$num_threads, n.epochs=umap_params$n.epochs, metric=umap_params$metric, seed.use=umap_params$seed.use,
                    negative.sample.rate=umap_params$negative.sample.rate)
cat("UMAP_2dimensions ran in ",difftime(Sys.time(), ttt, units='mins'),"\n")
invisible(gc())
    
umap10_params <- list(dims = 1:top_PCs,
                          n.neighbors = 15,
                          spread = 1,
                          repulsion.strength = 1,
                          learning.rate = 0.2,
                          min.dist= 0.1,
                          verbose = T,
                          num_threads=0,
                          n.epochs = 200,
                          metric = "euclidean",
                          negative.sample.rate = 5L,
                          seed.use = 42)
# Overwrite defaults with specified parameters
umap10_params <- modifyList(umap10_params, dim_reduct_params$umap10)
    
ttt <- Sys.time()
DATA <- RunUMAP(object=DATA, reduction=casefold("mnn"), n.components=10, reduction.name="umap10", reduction.key="umap10_",
                    dims=umap10_params$dims, n.neighbors=umap10_params$n.neighbors, spread=umap10_params$spread, repulsion.strength=umap10_params$repulsion.strength,
                    learning.rate=umap10_params$learning.rate, min.dist= umap10_params$min.dist, verbose=umap10_params$verbose, num_threads=umap10_params$num_threads,
                    n.epochs=umap10_params$n.epochs, metric=umap10_params$metric, seed.use=umap10_params$seed.use, negative.sample.rate=umap10_params$negative.sample.rate)
cat("UMAP_10dimensions ran in ",difftime(Sys.time(), ttt, units='mins'))
invisible(gc())

write.csv(DATA@reductions$umap@cell.embeddings, paste0(out_path,"/results/umap_plots/UMAP_coordinates.csv"))
write.csv(DATA@reductions$umap10@cell.embeddings, paste0(out_path,"/results/umap_plots/UMAP10_coordinates.csv"))
invisible(gc())

#---------

###################
### Running SNN ###
###################
cat("\n### Running SNN on ",casefold("mnn")," ###\n")
DATA <- FindNeighbors(DATA, assay = "RNA", graph.name="SNN", prune.SNN = .2,k.param = 30,force.recalc = T,reduction = casefold("mnn"), dims = 1:top_PCs )
g <- graph_from_adjacency_matrix(DATA@graphs$SNN,weighted = T,diag=F)
g <- simplify(g)
saveRDS(DATA@graphs$SNN, file = paste0(out_path,"/results/SNN_Graph.rds") )

for(i in c("pca", "umap")){
  png(filename = paste0(out_path,"/results/",i,"_plots","/",i,"_plot_with_SNN_overlay.png"),width = 200,height =205,res = 600,units = "mm")
  plot.igraph(x = g, layout = DATA@reductions[[i]]@cell.embeddings[,1:2], edge.width = E(graph = g)$weight/4, vertex.label = NA,
              edge.color = colorRampPalette(c("grey90","black"))(50)[round(E(graph = g)$weight/2*49+1)],edge.arrow.size=E(graph = g)$weight*0,
              vertex.size = 1,vertex.frame.color=hue_pal(l=50, c=80)(length(levels(factor(DATA$orig.ident))))[factor(DATA$orig.ident)],
              vertex.color = hue_pal()(length(unique(DATA$orig.ident)))[factor(DATA$orig.ident)])
  invisible(dev.off())
}
rm(g); invisible(gc())
#---------


#########################################
### Plotting Dimensionality Reduction ###
#########################################
mtdt <- colnames(DATA@meta.data) [ grepl("nFeature|nCount|_index|[.]Score",colnames(DATA@meta.data) ) ]
mtdt <- c(mtdt, "perc_mito" ,"perc_rps","perc_rpl","perc_hb", "perc_protein_coding" ,"perc_lincRNA","perc_snRNA","perc_miRNA","perc_processed_pseudogene",
          "perc_unknown","perc_Chr_1","perc_Chr_X","perc_Chr_Y","perc_Chr_MT")
mtdt <- mtdt[mtdt %in% colnames(DATA@meta.data)]
j <- c("sample", "patient_ID")

col_scale <- c("grey85","navy")
for(i in c("pca", "umap")){
  temp <- FeaturePlot(object = DATA, features = mtdt, cols = col_scale,pt.size = .5,reduction = i,ncol = 5,dims = 1:2)
  ggplot2::ggsave(temp,filename = paste0(i,"_metadata_dim1_dim2.png"), path = paste0(out_path,"/results/",i,"_plots"), dpi = 300,units = "mm",width = 170*5,height = 150*ceiling(length(mtdt)/5),limitsize = FALSE )
  
  temp2 <- DimPlot(DATA,dims = 1:2,reduction = i,group.by = j,pt.size = .3,ncol = 5)
  ggplot2::ggsave(temp2,filename = paste0(i,"_metadata_factors_dim1_dim2.png"), path = paste0(out_path,"/results/",i,"_plots"), dpi = 300,units = "mm",width = 170*5,height = 150*ceiling(length(j)/5),limitsize = FALSE )
  
  if(i == "pca"){
    temp <- FeaturePlot(object = DATA, features = mtdt, cols = col_scale,pt.size = .5,reduction = i,ncol = 5,dims = 3:4)
    ggplot2::ggsave(temp,filename = paste0(i,"_metadata_dim3_dim4.png"), path = paste0(out_path,"/results/",i,"_plots"), dpi = 300,units = "mm",width = 170*5,height = 150*ceiling(length(mtdt)/5),limitsize = FALSE )
    
    temp2 <- DimPlot(DATA,dims = 3:4,reduction = i,group.by = j,pt.size = .3,ncol = 5)
    ggplot2::ggsave(temp2,filename = paste0(i,"_metadata_factors_dim3_dim4.png"), path = paste0(out_path,"/results/",i,"_plots"), dpi = 300,units = "mm",width = 170*5,height = 150*ceiling(length(j)/5),limitsize = FALSE )
  } }
rm(temp,temp2); invisible(gc())
#---------


################################
### Clustering using Louvain ###
################################

if(!dir.exists(paste0(out_path,"/results/clustering"))){dir.create(paste0(out_path,"/results/clustering"))}
#remove old clustering with this name
DATA@meta.data <- DATA@meta.data[ , ! grepl( "louvain_", colnames(DATA@meta.data) ) ]
  
DATA <- FindClusters(object = DATA, reduction.type = "pca", dims.use = 1:top_PCs, resolution = seq(.05,2,by=.05), verbose = T,graph.name = "SNN",algorithm = 1)
colnames(DATA@meta.data) <- sub("SNN_res.","louvain_",colnames(DATA@meta.data))
 
for(i in c("pca","umap")){
  s <- colnames(DATA@meta.data)[grep("louvain_",colnames(DATA@meta.data))]
  plot_list <- lapply(s,function(j){ DimPlot(DATA,dims = 1:2,reduction = i,group.by = j, pt.size = .3,ncol = 5 ,label = T) +
        ggplot2::theme(legend.position = "none") + ggplot2::ggtitle(label =j) })
  p <- cowplot::plot_grid(plotlist = plot_list, ncol=5)
  ggplot2::ggsave(p,filename = paste0("clustering_louvain_",i,".png"), path = paste0(out_path,"/results/clustering"), dpi = 300,
                    units = "mm",width = 150*5,height = 140*ceiling(length(plot_list)/5),limitsize = FALSE )
}
  

rm(temp2); invisible(gc())
#---------

#######################
### K-MEANS PARTITIONING on UMAP ###
#######################

# Not running hclust because size of distance matrix exceeds maximum (65536)
if(!dir.exists(paste0(out_path,"/results/clustering"))){dir.create(paste0(out_path,"/results/clustering"))}
#remove old clustering with this name
DATA@meta.data <- DATA@meta.data[ , ! grepl( "kmeans_", colnames(DATA@meta.data) ) ]

cat("\n### Clustering with K-means on UMAP-10dims ###\n")

mincells <- 15
ideal <- round(ncol(DATA) / mincells)
clcl<- kmeans(DATA@reductions$umap10@cell.embeddings,centers = ideal,iter.max = 50)
DATA <- AddMetaData(DATA,metadata = setNames(clcl$cluster,colnames(DATA)), paste0("kmeans_",ideal))

temp <- rowsum(DATA@reductions$umap10@cell.embeddings,clcl$cluster) / as.vector(table(clcl$cluster))
cors <- cor(t(temp))

cat("\nMerging clusters ...\n")
merge_par <- seq(.8,.99,.02)

for(j in merge_par){
  if( j > min(cors) ){
    tcors <- (cors > j)*1
    cell_clust <- clcl$cluster
    clust <- rownames(tcors)
    
    for( i in clust){
      sel <- rownames(tcors)[ tcors[i,] > 0 ]
      cell_clust[cell_clust %in% sel] <- sel[1]
    }
    DATA <- AddMetaData(object = DATA, metadata = cell_clust, col.name = paste0("kmeans_merged_",j))
    
    temp <- UMAPPlot(object = DATA, group.by=paste0("kmeans_merged_",j), pt.size = .5) +
      ggplot2::theme(legend.position = "none")
    ggplot2::ggsave(temp,filename = paste0("results/clustering/UMAP_kmeans_merged_",j,".png"), path = out_path, dpi = 300,units = "mm",width = 170,height = 150 )
  }
}


for(i in c("pca", "umap")){
  s <- colnames(DATA@meta.data)[grep("kmeans_",colnames(DATA@meta.data))]
  plot_list <- lapply(s,function(j){ DimPlot(DATA,dims = 1:2,reduction = i,group.by = j, pt.size = .3,ncol = 5 ,label = T) +
      ggplot2::theme(legend.position = "none") + ggplot2::ggtitle(label = j) })
  p <- cowplot::plot_grid(plotlist = plot_list, ncol=5)
  ggplot2::ggsave(p,filename = paste0("clustering_kmeans_",i,".png"), path = paste0(out_path,"/results/clustering"), dpi = 300,
                  units = "mm",width = 150*5,height = 140*ceiling(length(plot_list)/5),limitsize = FALSE )
}

  

rm(temp); invisible(gc())
#---------


###################################
### SAVING RAW Seurat.v3 OBJECT ###
###################################
cat("\n### Saving the Seurat object ###\n")
saveRDS(DATA, file = paste0(out_path,"/results/output_rds/18Jan_sauron_seurat_object.rds") )
write.csv2(DATA@meta.data,paste0(out_path,"/results/output_rds/Metadata_with_clustering.csv"),row.names = T)
#---------
