library(dowser)
library(ggtree)
library(tidyverse)
library(alakazam)
library(shazam)
future::plan("multisession")
results_db <- readChangeoDb("~/immcantation-python3.8.12/imm-env/makedb_output/AIRR/functional/filtered_from_Seurat/seq_not_combined/genotype/light_cluster_clones/1_germline_germ-pass.tsv")
clones <- formatClones(results_db,
                       v_call = "v_call_genotyped",
                       nproc = 4,
                       traits = c("c_call", "integrated_snn_res.0.4"))

trees <-  getTrees(clones, build = "dnapars",
                   exec = "~/immcantation-python3.8.12/imm-env/bin/bin/phylip-3.697/exe/dnapars",
                   nproc = 4)
pal <- c("#8E0152","#C51B7D","#DE77AE","#F1B6DA","#FDE0EF",
         "#E6F5D0","#B8E186","#7FBC41","#4D9221","#276419",
         "#67001F","#B2182B","#D6604D","#F4A582", "#D1E5F0","#92C5DE",
         "#4393C3","#2166AC","#053061")
plotTrees(trees)[[19]] +
  geom_tippoint(aes(colour = factor(integrated_snn_res.0.4), size=60)) +
  geom_tiplab(aes(label = c_call), offset = 0.002) +
  scale_colour_manual(values = pal)

