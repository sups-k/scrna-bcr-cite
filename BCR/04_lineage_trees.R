library(dowser)
library(ggtree)
library(tidyverse)
library(alakazam)
library(shazam)

future::plan("multicore")
# Read reconstructed germlines
results_db <- readChangeoDb("~/immcantation-python3.8.12/imm-env/makedb_output/functional/filtered_from_Seurat/genotype/create_germlines/all_ph_genotyped_germ-pass.tab")

# Check output column
#results_db$GERMLINE_IMGT_D_MASK[1]
results_db$germline_alignment_d_mask[1]

# Read Seurat metadata
meta <- read.csv(file = "~/Jocelyn_scRNA/KW9079_Farmer/210624_A00873_0376_BHC3KTDSX2_KW9079/cellranger-6.0.1/GRCh38/results/10_pc43_11March_METADATA.csv") %>%
  as_tibble()
# Rename cells column
meta <- meta %>% rename(CELL = cells)

# Merge metadata to BCR information
merged_db <- merge(meta, results_db, by = "CELL")
# Select only cluster ID
merged_db <- merged_db %>% select(colnames(results_db), integrated_snn_res.0.4)

# Make clone objects with aligned, processed sequences
# Collapse identical sequences unless differ by trait
# Add up duplicate_count column for collapsed sequences
# Store isotype, cluster ID
# Discard clones with < 5 distinct sequences

clones <- formatClones(merged_db,
                       id = "SEQUENCE_ID",
                       seq = "SEQUENCE_IMGT",
                       clone = "CLONE",
                       germ = "GERMLINE_IMGT_D_MASK",
                       v_call = "V_CALL",
                       j_call = "J_CALL",
                       junc_len = "JUNCTION_LENGTH",
                       nproc = 6,
                       cell = "CELL",
                       locus = "LOCUS",
                       traits = c("C_CALL", "integrated_snn_res.0.4"),
                       num_fields = "UMICOUNT",
                       minseq = 5)

# B cell specific maximum likelihood with IgPhyML
trees <-  getTrees(clones, build = "igphyml",
                 exec = "~/immcantation-python3.8.12/imm-env/bin/igphyml/src/igphyml",
                 nproc = 6)
save.image(file = "~/immcantation-python3.8.12/imm-env/lineage_trees.RData")
saveRDS(trees, file = "~/immcantation-python3.8.12/imm-env/lineage_trees.rds")
# # Plot all trees
# plots <- plotTrees(trees, tips = "isotype", tipsize = 2)
# 
# # Plot the largest tree
# plots[[1]]
# 
# # Save PDF of all trees
# treesToPDF(plots, file = "final_data_trees.pdf", nrow = 2, ncol = 2)
# 
# 
# # Scale branches to mutations rather than mutations/site
# trees <-  scaleBranches(trees)
# 
# # Make fancy plot of second largest tree
# plotTrees(trees, scale = 5)[[2]] +
#   geom_tippoint(aes(colour = cluster_ID)) +
#   geom_tiplab(aes(label = isotype), offset = 0.002) +
#   scale_colour_manual(values = col_anno)
# 
# 
# 
# 
# # Reconstruct sequences at internal nodes of the tree
# 
# # Collapse nodes with identical sequences
# trees <- collapseNodes(trees)
# 
# # node_nums=TRUE labels each internal node
# p <- plotTrees(trees, node_nums = TRUE, labelsize = 6, scale = 5)[[2]] +
#   geom_tippoint(aes(colour = cluster_ID)) +
#   geom_tiplab(aes(label = isotype), offset = 0.002) +
#   scale_colour_manual(values = col_anno)
# 
# # Get sequence at node 26
# getSeq(trees, clone = trees$clone_id[2], node = 26)
