library(dowser)
library(ggtree)

# Read in IMGT data
references <- readIMGT(dir = "~/immcantation-python3.8.12/imm-env/imgt/human/vdj")

# Reconstruct germlines
results_db_2 <- createGermlines(results_db, references)

# Check output column
results_db_2$germline_alignment_d_mask[1]

# Make cline objects with aligned, processed sequences
# Collapse identical sequences unless differ by trait
# Add up duplicate_count column for collapsed sequences
# Store isotype, cluster ID
# Discard clones with < 5 distinct sequences

clones <- formatClones(results_db_2,
                       traits = c("isotype", "cluster_ID"),
                       num_fields = c("duplicate_count"),
                       minseq = 5)

# B cell specific maximum likelihood with IgPhyML
trees = getTrees(clones, build = "igphyml",
                 exec = "src/igphyml",
                 nproc = 2)

# Plot all trees
plots <- plotTrees(trees, tips = "isotype", tipsize = 2)

# Plot the largest tree
plots[[1]]

# Save PDF of all trees
treesToPDF(plots, file = "final_data_trees.pdf", nrow = 2, ncol = 2)


# Scale branches to mutations rather than mutations/site
trees <-  scaleBranches(trees)

# Make fancy plot of second largest tree
plotTrees(trees, scale = 5)[[2]] +
  geom_tippoint(aes(colour = cluster_ID)) +
  geom_tiplab(aes(label = isotype), offset = 0.002) +
  scale_colour_manual(values = col_anno)




# Reconstruct sequences at internal nodes of the tree

# Collapse nodes with identical sequences
trees <- collapseNodes(trees)

# node_nums=TRUE labels each internal node
p <- plotTrees(trees, node_nums = TRUE, labelsize = 6, scale = 5)[[2]] +
  geom_tippoint(aes(colour = cluster_ID)) +
  geom_tiplab(aes(label = isotype), offset = 0.002) +
  scale_colour_manual(values = col_anno)

# Get sequence at node 26
getSeq(trees, clone = trees$clone_id[2], node = 26)