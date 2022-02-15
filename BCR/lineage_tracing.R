library(alakazam)
library(igraph)
library(dplyr)


# specify directory containing genotyped VDJ files
geno_dir <- "~/immcantation-python3.8.12/imm-env/makedb_output/functional/genotype/create_germlines/"

# get list of available files
geno_files <- dir(geno_dir, 'germ-pass[.]tab', full.names=T)

# specify which file to load
g_file <- geno_files[2]

# load the Change-O database file with germline sequence information (*_germ-pass.tab file)
db <- readChangeoDb(g_file)

# select a desired clone
head(sort(table(db$CLONE), decreasing=T))  # view most abundant clone IDs
sub_db <- subset(db, CLONE == 3202)

# create ChangeOclone object for clone
clone <- makeChangeoClone(sub_db,
                          id = "SEQUENCE_ID",
                          seq = "SEQUENCE_IMGT",
                          germ = "GERMLINE_IMGT",
                          v_call = "V_CALL",
                          j_call = "J_CALL",
                          junc_len = "JUNCTION_LENGTH",
                          clone = "CLONE",
                          text_fields=c("CELL","C_CALL"),
                          num_field="UMICOUNT")

# Run PHYLIP and parse output
# Need to download and install PHYLIP from http://evolution.genetics.washington.edu/phylip/getme-new1.html
phylip_exec <- "~/immcantation-python3.8.12/imm-env/bin/phylip-3.697/exe/dnapars"  # path to dnapars executable within PHYLIP package
g <- buildPhylipLineage(clone,
                        phylip_exec,
                        rm_temp=TRUE,
                        verbose=TRUE)

# Modify graph and plot attributes
V(g)$color <- "steelblue"
V(g)$color[V(g)$name == "Germline"] <- "black"
V(g)$color[grepl("Inferred", V(g)$name)] <- "white"
V(g)$label <- V(g)$C_CALL
E(g)$label <- ""

# Remove large default margins
# par(mar=c(0, 0, 0, 0) + 0.1)

# Plot graph
plot(g, layout=layout_as_tree, edge.arrow.mode=0, vertex.frame.color="black",
     vertex.label.color="black", vertex.size=5, vertex.label.family='sans', vertex.label.cex=0.7)

# Add legend
legend("topleft", c("Germline", "Inferred", "Sample"), 
       fill=c("black", "white", "steelblue"), cex=0.75)
