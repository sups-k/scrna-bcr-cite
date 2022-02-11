library(Seurat)
library(dplyr)
library(scales)
library(RColorBrewer)
library(igraph)
#library(sva)
library(rafalib)
library(parallel)
library(batchelor)
#library(scran)
#library(scater)

# var_to_regress='nCount_RNA,nFeature_RNA,nCount_ADT,perc_mito'
future::plan("multisession")
options(future.globals.maxSize = 100000 * 1024^2)



# Loading object
DATA <- readRDS(file = "results/output_rds/sauron_filtered.rds")

########################################
### SCALE DATA AND REGRESS VARIABLES ###
########################################
vars_regress <- c("nCount_RNA", "nFeature_RNA", "mitoRatio")

fast_ScaleData <- function(DATA, assay="RNA",vars.to.regress=NULL,do.scale=T,do.center=T,scale.max=NULL){
  require(parallel)
  #Check spelling to see if the variables to regress are in the dataset
  vars.to.regress <- vars.to.regress[vars.to.regress %in% colnames(DATA@meta.data)]
  cat("The following variables were found in your data and will be used to regress out the counts\n")
  print(vars.to.regress)
  
  #future::plan("multisession")
  
  #Regress factors, if any
  if(!is.null(vars.to.regress)){
    cat("Your data contains ",ncol(eval(parse(text=paste0("DATA@assays$",assay,"@data")))),"samples and",nrow(eval(parse(text=paste0("DATA@assays$",assay,"@data"))))," features\n")
    cat("fast_ScaleData regression started running, please wait ...\n")
    
    l1 <- apply(eval(parse(text=paste0("DATA@assays$",assay,"@data"))),1,function (x){
      model <- paste0("x ~ ", paste(paste0("DATA$",vars.to.regress),collapse = "+"))
      m <- glm(eval(parse(text=model)))
      return(scale(m$residuals,T,T))
    })
    
    l2 <- apply(l1,2,function (x){  return(scale(x,do.scale,do.center)) })
    
  }
  
  invisible(gc())
  rownames(l2) <- colnames(DATA@assays[[assay]]@data)
  rm(l1)
  invisible(gc())
  
  #Assign it back to the DATA object
  DATA@assays[[assay]]@scale.data <-  as.matrix(t(l2))
  return(DATA)
}


DATA <- fast_ScaleData(DATA, vars.to.regress=vars_regress, assay="RNA")



####################################
### INTEGRATE DATASETS USING MNN ###
####################################

cat("\n### INTEGRATING DATASETS USING MNN ###\n")

  
# Defining batch variables
cat("\nCreating dataset list\n")
batch <- as.character(factor(DATA@meta.data[,"patient_ID"]))
DATA.list <- SplitObject(DATA, split.by = "patient_ID")


compute_hvgs <- function(object,VAR_choice,output_path,assay="RNA", nSeurat=3000){

  if(!dir.exists(output_path)){dir.create(output_path,recursive = T)}
  perc1 <- rowSums(as.matrix(object@assays[[assay]]@data) > 0) / ncol(object@assays[[assay]]@data)
  perc <- names(perc1)[ (perc1 < 0.95) & (perc1 > 1/ncol(object@assays[[assay]]@data) ) ]
    
  #########################################################
  ### Running SEURAT method for variable gene selection ###
  #########################################################
  cat("\nCalculating highly variable genes with Seurat ...\n")
  y_cut <- 1
    
    
  #Defining the variable genes based on the mean gene expression above the 5% quantile and the dispersion above 2.
  object <- FindVariableFeatures(object = object,nfeatures = nSeurat)
  m <- max(quantile(object@assays[[assay]]@meta.features$vst.mean,probs = c(.025)) , 0.01)
    
  object@assays[[assay]]@var.features <- object@assays[[assay]]@var.features[object@assays[[assay]]@var.features %in% perc]
  object@assays[[assay]]@meta.features$use <- rownames(object@assays[[assay]]@meta.features) %in% object@assays[[assay]]@var.features
  write.csv2(object@assays[[assay]]@meta.features, paste0(output_path,"/HVG_info_seurat.csv"))
    
  png(filename = paste0(output_path,"/Var_vst_exp_disp_gene_selection_seurat.png"),width = 1000,height = 1050,res = 200)
  plot(log2(object@assays[[assay]]@meta.features$vst.variance.expected),object@assays[[assay]]@meta.features$vst.variance.standardized,cex=.1,main="HVG selection",
         col=ifelse(rownames(object@assays[[assay]]@meta.features)%in% object@assays[[assay]]@var.features,"red","black" ),ylab="vst.variable",xlab="log2(var.expected)")
  abline(v=log2(m),h=y_cut,lty=2,col="grey20",lwd=1)
  invisible(dev.off())
    
  png(filename = paste0(output_path,"/Var_vst_mean_disp_gene_selection_seurat.png"),width = 1000,height = 1050,res = 200)
  plot(log2(object@assays[[assay]]@meta.features$vst.mean),object@assays[[assay]]@meta.features$vst.variance.standardized,cex=.1,main="HVG selection",
         col=ifelse(rownames(object@assays[[assay]]@meta.features)%in% object@assays[[assay]]@var.features,"red","black" ),ylab="vst.variable",xlab="log2(mean expression)")
  abline(v=log2(m),h=y_cut,lty=2,col="grey20",lwd=1)
  invisible(dev.off())
    #---------
  return(object)
}



# Define HVGs per dataset
cat("\nComputing HVGs\n")
VAR_choice <- "seurat"
for (i in 1:length(DATA.list)) {
  cat("\nProcessing dataset: ",names(DATA.list)[i]," \n")
  DATA.list[[i]] <- compute_hvgs(DATA.list[[i]],VAR_choice,paste0(out_path,"/results/var_genes_",names(DATA.list)[i]),assay = "RNA", nSeurat=3000)
}
    
# Select the most informative genes that are shared across all datasets:
cat("\nComputing HVGs\n")
universe <- unique(unlist(lapply(DATA.list,function(x){x@assays[["RNA"]]@var.features})))
head(universe, 50)
cat("\n",length(universe)," genes found as variable within datasets\n")
    
#Separating batch matricies
DATA.list <- lapply(DATA.list, function(x){ x@assays[["RNA"]]@scale.data[universe,]} )
myinput <- list()
    
myinput[["k"]] <- 20

myinput[["approximate"]] <-  TRUE
myinput[["d"]] <-  51

    
#Applying MNN correction on scaled counts
cat("\nApplying MNN correction on scaled counts\n")
out <- do.call(fastMNN, args=c(DATA.list, myinput))
cat("\nMNN computation done\n")

# Corrected values for use in clustering, etc.
out2 <- t(reducedDim(out))
colnames(out2) <- unlist(lapply(DATA.list,function(x){colnames(x)}))
out2 <- out2[,colnames(DATA)]
rownames(out2) <- paste0("dim",1:myinput$d)
DATA@reductions[["mnn"]] <- CreateDimReducObject(embeddings = t(out2),key = "MNN_",assay = "RNA")
rm(out, out2, myinput)
invisible(gc())


#---------

###########################
### FIND VARIABLE GENES ###
###########################
if(DefaultAssay(DATA) == "RNA"){
  output_path <- paste0(out_path,"results/All_datasets_together")
  DATA <- compute_hvgs(DATA, VAR_choice, output_path, assay="RNA", nSeurat=3000)
}
#---------

###################################
### SAVING Seurat.v3 OBJECT ###
###################################
cat("\n### Saving Seurat object ###\n")
saveRDS(DATA, file = paste0(out_path,"/12Jan_seurat_object.rds") )
#---------

