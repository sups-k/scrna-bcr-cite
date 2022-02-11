# Adapted from https://kleinstein.bitbucket.io/tutorials/intro-lab/index.html
# and https://github.com/angelettilab/scMouseBcellFlu/blob/master/scripts/VDJ_analysis/01_VDJ_genotype_and_threshold.R

# Load libraries
library(alakazam)
library(shazam)
library(tigger)

patient_ID <- c("Hashtag1", "Hashtag2", "Hashtag3", "Hashtag4", "Hashtag5", "Hashtag6", "Hashtag7", "Hashtag8")

#####################################################
### FIND NOVEL V SEQ ALLELES AND INFER GENOTYPE ###
#####################################################

# V(D)J assignment is a key step in analyzing repertoires, and it is done by
# matching the sequences against a database of known V(D)J alleles. However,
# current databases are incomplete, and this process will fail for sequences
# that utilize previously undetected alleles. Some assignments may also occur to
# genes which are not carried by the individual. The TIgGER R package infers
# subject-specific V genotypes, including novel alleles, and then uses these
# results to improve the V gene annotations.

# Load V-segment germline sequences
ighv <- readIgFasta("~/immcantation-python3.8.12/imm-env/imgt/human/vdj/imgt_human_IGHV.fasta")

# Import ChangeO-formatted sequence database files (heavy and light chains)
seqdb <- NULL

# For healthy 1
tmp <- readChangeoDb("~/immcantation-python3.8.12/imm-env/makedb_output/functional/heavy/1_heavy_parse-select.tab")
tmp$PATIENT_ID <- patient_ID[1]
seqdb <- rbind(seqdb, tmp)
# For patient 2
tmp <- readChangeoDb("~/immcantation-python3.8.12/imm-env/makedb_output/functional/heavy/2_heavy_parse-select.tab")
tmp$PATIENT_ID <- patient_ID[2]
seqdb <- rbind(seqdb, tmp)
# For healthy 3
tmp <- readChangeoDb("~/immcantation-python3.8.12/imm-env/makedb_output/functional/heavy/3_heavy_parse-select.tab")
tmp$PATIENT_ID <- patient_ID[3]
seqdb <- rbind(seqdb, tmp)
# For patient 4
tmp <- readChangeoDb("~/immcantation-python3.8.12/imm-env/makedb_output/functional/heavy/4_heavy_parse-select.tab")
tmp$PATIENT_ID <- patient_ID[4]
seqdb <- rbind(seqdb, tmp)
# For healthy 5
tmp <- readChangeoDb("~/immcantation-python3.8.12/imm-env/makedb_output/functional/heavy/5_heavy_parse-select.tab")
tmp$PATIENT_ID <- patient_ID[5]
seqdb <- rbind(seqdb, tmp)
# For healthy 6
tmp <- readChangeoDb("~/immcantation-python3.8.12/imm-env/makedb_output/functional/heavy/6_heavy_parse-select.tab")
tmp$PATIENT_ID <- patient_ID[6]
seqdb <- rbind(seqdb, tmp)
# For patient 7
tmp <- readChangeoDb("~/immcantation-python3.8.12/imm-env/makedb_output/functional/heavy/7_heavy_parse-select.tab")
tmp$PATIENT_ID <- patient_ID[7]
seqdb <- rbind(seqdb, tmp)
# For patient 8
tmp <- readChangeoDb("~/immcantation-python3.8.12/imm-env/makedb_output/functional/heavy/8_heavy_parse-select.tab")
tmp$PATIENT_ID <- patient_ID[8]
seqdb <- rbind(seqdb, tmp)

# Infer genotype (performed separately for each patient)
rm(tmp)
gc()

# For each V gene allele, findNovelAlleles analyzes the sequences that were
# assigned the allele, and evaluates the apparent mutation frequency at each
# position as a function of sequence-wide mutation counts. Positions that
# contain polymorphisms (rather than somatic hypermutations) will exhibit a high
# apparent mutation frequency even when the sequence-wide mutation count is low.

genotype_dir <- "~/immcantation-python3.8.12/imm-env/makedb_output/functional/"


for(i in patient_ID){
  seqdb_patient <- seqdb[seqdb$PATIENT_ID %in% i, ]
  novel_rows <- NULL
  nv <- NULL
  nv <- findNovelAlleles(data = seqdb_patient,
                         germline_db = ighv,
                         seq = "SEQUENCE_IMGT",
                         v_call = "V_CALL",
                         j_call = "J_CALL",
                         junction = "JUNCTION",
                         junction_length = "JUNCTION_LENGTH")
  
  # Show novel alleles
  novel_rows <- selectNovel(nv)
  # Visualize first novel allele if it exists
  if (!is.null(novel_rows) && (nrow(novel_rows) > 0)){
    pdf(file = paste0(genotype_dir, "/novel_alleles/", i, "_novel_alleles.pdf"), height=10, width=10)
    plotNovel(data = seqdb_patient,
              novel_row = novel_rows[1, ],
              v_call = "V_CALL",
              j_call = "J_CALL",
              seq = "SEQUENCE_IMGT",
              junction = "JUNCTION",
              junction_length = "JUNCTION_LENGTH")
    invisible(dev.off())
  }
  
  
  # Infer patient genotype
  gt_pat <- inferGenotype(data = seqdb_patient,
                          germline_db = ighv,
                          novel = nv,
                          v_call = "V_CALL",
                          seq = "SEQUENCE_IMGT")
  
  pdf(file = paste0(genotype_dir, "genotype/IGHV_genotype_plot_", i, ".pdf"), height=10, width=10)
  plotGenotype(gt_pat, gene_sort = "position", text_size = 8)
  invisible(dev.off())
  
  # Convert genotype table to vector of nucleotide sequences
  gtseq_pat <- genotypeFasta(genotype = gt_pat, germline_db = ighv, novel=nv)
  writeFasta(gtseq_pat, paste0(genotype_dir, "genotype/IGHV_genotype_", i, ".fasta"))
  
  # Correct allele calls based on the personalized genotype
  seqdb_patient <- reassignAlleles(data = seqdb_patient,
                                   genotype_db = gtseq_pat,
                                   v_call = "V_CALL",
                                   seq = "SEQUENCE_IMGT")
  
  # Save data
  writeChangeoDb(seqdb_patient, paste0(genotype_dir, "genotype/IGHV-genotyped_", i, ".tab"))
  
}


### 1 novel allele in healthy 1
# Warning messages:
# In max(FRACTION) : no non-missing arguments to max; returning -Inf
### No novel alleles in patient 2
### No novel alleles in healthy 3
### No novel alleles in patient 4
### No novel alleles in healthy 5
### No novel alleles in healthy 6
# Warning messages:
# 1: In max(FRACTION) : no non-missing arguments to max; returning -Inf
# 2: In max(FRACTION) : no non-missing arguments to max; returning -Inf
# 3: In max(FRACTION) : no non-missing arguments to max; returning -Inf
# 4: In max(FRACTION) : no non-missing arguments to max; returning -Inf
### No novel alleles in patient 7
### No novel alleles found for patient 8







####################################################
####### CLONAL DIVERSITY ANALYSIS USING CDR3 #######
####################################################

# It practice, we first split sequences into groups that share the same V and J
# gene assignments, and that have the same junction (or, equivalently CDR3)
# length. This is based on the assumption that members of a clone will
# necessarily share all of these properties. distToNearest performs this
# grouping step and then, for each group, it counts the number of mismatches in
# the junction region between all pairs of sequences in the group, and returns
# the smallest non-zero value for each sequence. At the end of this step, a new
# column, dist_nearest, which contains the distances to the closest
# (non-identical) sequence in the group, will be added to seqdb_patient.

# findThreshold uses the distribution of distances calculated in the previous
# step to determine an appropriate threshold for the dataset. This can be done
# using either a density or mixture based method. The function plot can be used
# to visualize the distance-to-nearest distribution and the threshold.

predicted_thresholds <- data.frame(patient_ID = patient_ID, threshold=rep(0.1, length(patient_ID)))

for(i in patient_ID){
  seqdb_patient <- seqdb[seqdb$PATIENT_ID %in% i, ]
  # Calculate distances to nearest neighbors
  dist_ham <- distToNearest(db = seqdb_patient,
                           sequenceColumn = "JUNCTION",
                           vCallColumn="V_CALL",
                           jCallColumn = "J_CALL",
                           model="ham",
                           normalize="len",
                           nproc=1)

  # Perform automatic threshold estimation
  output <- findThreshold(distances = dist_ham$dist_nearest, method = "density")

  # Save threshold values
  if (is.null(output) || is.na(output@threshold)) {  
    cat('\n\tThreshold estimation failed. Reverting to default threshold: 0.1', '\n')
  } else {
    predicted_thresholds$threshold[predicted_thresholds$patient_ID == i] <- output@threshold
  }
}

# Save predicted thresholds
write.csv(predicted_thresholds, file = "~/immcantation-python3.8.12/imm-env/makedb_output/functional/predicted_thresholds.csv", quote = FALSE, row.names = FALSE)


### USING SCOPER ####
# In some datasets, the distance-to-nearest distribution may not be bimodal and
# findThreshold may fail to determine the distance threshold. In these cases,
# spectral clustering with an adaptive threshold to determine the local sequence
# neighborhood may be used. This is be done with functions from the SCOPer R package.


library(scoper)

# Initialize predicted thresholds table
predicted_thresholds <- data.frame(patient_ID = patient_ID, threshold=rep(0.1, length(patient_ID)))

for(i in patient_ID){
  seqdb_patient <- seqdb[seqdb$PATIENT_ID %in% i, ]

  # Running spectralClones since findThreshold gave unimodal distribution
  results <- spectralClones(db = seqdb_patient,
                 method = "novj",
                 germline = "GERMLINE_IMGT",
                 sequence = "SEQUENCE_IMGT",
                 junction = "JUNCTION",
                 v_call = "V_CALL",
                 j_call = "J_CALL",
                 locus = "LOCUS",
                 cell_id = "SEQUENCE_ID",
                 split_light = FALSE,
                 only_heavy = TRUE)
  
  # To get threshold, plot the results. Use the effective threshold value.
  # plot(results)
  
  # Add threshold to table - change for each patient
  predicted_thresholds$threshold[predicted_thresholds$patient_ID == i] <- results@eff_threshold
}

# Save predicted thresholds
write.csv(predicted_thresholds, file = "~/immcantation-python3.8.12/imm-env/makedb_output/functional/predicted_thresholds.csv", quote = FALSE, row.names = FALSE)
