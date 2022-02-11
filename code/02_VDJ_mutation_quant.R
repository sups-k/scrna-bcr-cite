# Using SHazaM to find mutational load

# Loading libraries
library(alakazam)
library(shazam)
library(tigger)
library(dplyr)
library(ggplot2)

# Get list of available files and iterate through each file
geno_files <- dir("/home/suk571/immcantation-python3.8.12/imm-env/makedb_output/functional/genotype/create_germlines", "germ-pass[.]tab", full.names=T)

db_full <- NULL

for (g_file in geno_files) {
  
  cat('Processing file:', basename(g_file), '...\n')
  
  # Load the Change-O database file with germline sequence information (*_germ-pass.tab file)
  db <- readChangeoDb(g_file)
  
  # Calculate mutation counts and frequencies
  # Mutaion counts per individual regions
  db <- observedMutations(db,
                          sequenceColumn="SEQUENCE_IMGT",
                          germlineColumn="GERMLINE_IMGT_D_MASK",
                          frequency=F,
                          regionDefinition=IMGT_V_BY_REGIONS,
                          combine=F)
  # Total mutation counts
  db <- observedMutations(db,
                          sequenceColumn="SEQUENCE_IMGT",
                          germlineColumn="GERMLINE_IMGT_D_MASK",
                          frequency=F,
                          combine=T)
  colnames(db)[colnames(db) == 'MU_COUNT'] <- 'MU_COUNT_TOT'
  
  # Mutation frequency per individual regions
  db <- observedMutations(db,
                          sequenceColumn="SEQUENCE_IMGT",
                          germlineColumn="GERMLINE_IMGT_D_MASK",
                          frequency=T,
                          regionDefinition=IMGT_V_BY_REGIONS,
                          combine=F)
  # Total mutation frequency
  db <- observedMutations(db,
                          sequenceColumn="SEQUENCE_IMGT",
                          germlineColumn="GERMLINE_IMGT_D_MASK",
                          frequency=T,
                          combine=T)
  colnames(db)[colnames(db) == 'MU_FREQ'] <- 'MU_FREQ_TOT'
  
  # Combine data with full database
  if (is.null(db_full)) {
    db_full <- db
  } else {
    db_full <- rbind(db_full, db)
  }
  
}


# Write the merged ChangeO database containing mutation data to a file
out_file <- "VDJseq_mutation_quant.tab"

cat('\nWriting merged ChangeO database file:', out_file, '...\n')

writeChangeoDb(db_full, file.path("/home/suk571/immcantation-python3.8.12/imm-env/makedb_output/functional/genotype/create_germlines", out_file))

cat('Done!\n\n')
gc()

