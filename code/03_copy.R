library(Seurat) 
library(dplyr)
library(scales)
library(RColorBrewer)
library(rafalib)
library(parallel)
library(scran)
library(dbscan)

# future::plan("multicore")
options(future.globals.maxSize = 100000 * 1024^2)
out_path <- "~/Jocelyn_scRNA/KW9079_Farmer/210624_A00873_0376_BHC3KTDSX2_KW9079/cellranger-6.0.1/GRCh38"
setwd("~/Jocelyn_scRNA/KW9079_Farmer/210624_A00873_0376_BHC3KTDSX2_KW9079/cellranger-6.0.1/GRCh38")

#############################
### LOAD Seurat.v3 OBJECT ###
#############################
cat("\n### LOADING Seurat.v3 OBJECT ###\n")
DATA <- readRDS(file = paste0(out_path, "/results/output_rds/17Janclusters.rds"))

#######################
### HBDSCAN on UMAP ###
#######################

# Not running hclust because size of distance matrix exceeds maximum (65536)

cat("\n### Clustering with HDBSCAN on UMAP-10dims ###\n")

  
for(k in seq(5,100,by=2)){
  cat(k,"\t")
  clusters <- hdbscan(DATA@reductions$umap10@cell.embeddings, minPts = k)
  names(clusters$cluster) <- rownames(DATA@meta.data)
  DATA <- AddMetaData(object = DATA, metadata = clusters$cluster, col.name = paste0("hdbscan_",k))
}
  
for(i in c("pca", "umap")){
s <- colnames(DATA@meta.data)[grep("hdbscan_",colnames(DATA@meta.data))]
plot_list <- lapply(s,function(j){ DimPlot(DATA,dims = 1:2,reduction = i,group.by = j, pt.size = .3,ncol = 5 ,label = T) +
      ggplot2::theme(legend.position = "none") + ggplot2::ggtitle(label = j) })
p <- cowplot::plot_grid(plotlist = plot_list, ncol=5)
ggplot2::ggsave(p,filename = paste0("clustering_hdbscan_",i,".png"), path = paste0(out_path,"/results/clustering"), dpi = 300,
                  units = "mm",width = 150*5,height = 140*ceiling(length(plot_list)/5),limitsize = FALSE )
}


rm(temp2); invisible(gc())
#---------

###################################
### SAVING Seurat.v3 OBJECT ###
###################################
cat("\n### Saving the Seurat object ###\n")
saveRDS(DATA, file = paste0(out_path,"/results/output_rds/Jan18seurat_object.rds") )
write.csv2(DATA@meta.data,paste0(out_path,"/results/output_rds/Metadata_with_clustering.csv"),row.names = T)
