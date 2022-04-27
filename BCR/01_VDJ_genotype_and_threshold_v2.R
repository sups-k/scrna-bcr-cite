library(alakazam)
library(shazam)
library(tigger)
library(scoper)
library(tidyverse)

future::plan("multisession")

# Load V-segment germline sequences
ighv <- readIgFasta("~/immcantation-python3.8.12/imm-env/imgt/human/vdj/heavy/imgt_human_IGHV.fasta")
# Import ChangeO-formatted sequence database files (heavy and light chains)
seqdb <- NULL
patient_ID <- c("Hashtag1", "Hashtag2", "Hashtag3", "Hashtag4", "Hashtag5", "Hashtag6", "Hashtag7", "Hashtag8")

# For healthy 1
tmp <- readChangeoDb("~/immcantation-python3.8.12/imm-env/makedb_output/AIRR/functional/filtered_from_Seurat/heavy/1_heavy_parse-select.tab")
novel_rows <- NULL
nv <- NULL
genotype_dir <- "~/immcantation-python3.8.12/imm-env/makedb_output/AIRR/functional/filtered_from_Seurat/"

nv <- findNovelAlleles(data = tmp,
                       germline_db = ighv,
                       seq = "SEQUENCE_IMGT",
                       v_call = "V_CALL",
                       j_call = "J_CALL",
                       junction = "JUNCTION",
                       junction_length = "JUNCTION_LENGTH")

novel_rows <- selectNovel(nv)

gt_pat <- inferGenotype(data = tmp,
                        germline_db = ighv,
                        novel = nv,
                        v_call = "V_CALL",
                        seq = "SEQUENCE_IMGT")

# Convert genotype table to vector of nucleotide sequences
gtseq_pat <- genotypeFasta(genotype = gt_pat, germline_db = ighv, novel=nv)
writeFasta(gtseq_pat, paste0(genotype_dir, "genotype/IGHV_genotype_Hashtag1.fasta"))

seqdb_patient <- reassignAlleles(data = tmp,
                                 genotype_db = gtseq_pat,
                                 v_call = "V_CALL",
                                 seq = "SEQUENCE_IMGT")

writeChangeoDb(seqdb_patient, paste0(genotype_dir, "genotype/IGHV-genotyped_Hashtag1.tab"))

rm(tmp, novel_rows, nv, gtseq_pat, gt_pat)

### I know that no other sample has novel alleles so I skipped the above step for the others ###

# Removing V_CALL and designating it to v_call_genotyped
seqdb_patient$patient_id <- patient_ID[1]
seqdb <- rbind(seqdb, seqdb_patient)
seqdb <- seqdb %>% select(-V_CALL)
seqdb <- seqdb %>% rename(V_CALL = v_call_genotyped)


# For patient 2
tmp <- readChangeoDb("~/immcantation-python3.8.12/imm-env/makedb_output/AIRR/functional/filtered_from_Seurat/heavy/2_heavy_parse-select.tab")
tmp$patient_id <- patient_ID[2]
seqdb <- rbind(seqdb, tmp)
# For healthy 3
tmp <- readChangeoDb("~/immcantation-python3.8.12/imm-env/makedb_output/AIRR/functional/filtered_from_Seurat/heavy/3_heavy_parse-select.tab")
tmp$patient_id <- patient_ID[3]
seqdb <- rbind(seqdb, tmp)
# For patient 4
tmp <- readChangeoDb("~/immcantation-python3.8.12/imm-env/makedb_output/AIRR/functional/filtered_from_Seurat/heavy/4_heavy_parse-select.tab")
tmp$patient_id <- patient_ID[4]
seqdb <- rbind(seqdb, tmp)
# For healthy 5
tmp <- readChangeoDb("~/immcantation-python3.8.12/imm-env/makedb_output/AIRR/functional/filtered_from_Seurat/heavy/5_heavy_parse-select.tab")
tmp$patient_id <- patient_ID[5]
seqdb <- rbind(seqdb, tmp)
# For healthy 6
tmp <- readChangeoDb("~/immcantation-python3.8.12/imm-env/makedb_output/AIRR/functional/filtered_from_Seurat/heavy/6_heavy_parse-select.tab")
tmp$patient_id <- patient_ID[6]
seqdb <- rbind(seqdb, tmp)
# For patient 7
tmp <- readChangeoDb("~/immcantation-python3.8.12/imm-env/makedb_output/AIRR/functional/filtered_from_Seurat/heavy/7_heavy_parse-select.tab")
tmp$patient_id <- patient_ID[7]
seqdb <- rbind(seqdb, tmp)
# For patient 8
tmp <- readChangeoDb("~/immcantation-python3.8.12/imm-env/makedb_output/AIRR/functional/filtered_from_Seurat/heavy/8_heavy_parse-select.tab")
tmp$patient_id <- patient_ID[8]
seqdb <- rbind(seqdb, tmp)


rm(tmp, seqdb_patient)
gc()

# Calculate threshold for CDR3 distance
dist_ham <- distToNearest(db = seqdb,
sequenceColumn = "JUNCTION",
vCallColumn = "V_CALL",
jCallColumn = "J_CALL",
model="ham",
normalize="len",
nproc=2)
output <- findThreshold(distances = dist_ham$dist_nearest, method = "density")

# Threshold was NA because distribution is unimodal.

# Define clones using SCOPer instead
results <- spectralClones(db = seqdb,
                          method = "vj", # Group cells by common VH, JH, junction length
                          germline = "GERMLINE_IMGT",
                          sequence = "SEQUENCE_IMGT",
                          junction = "JUNCTION",
                          v_call = "V_CALL",
                          j_call = "J_CALL",
                          clone = "CLONE",
                          locus = "LOCUS",
                          cell_id = "SEQUENCE_ID",
#split_light = FALSE,
#only_heavy = TRUE,
                          nproc = 2)

# View and note down threshold - 0.27
View(results@eff_threshold)

# Get results data frame
results_db <- as.data.frame(results)
writeChangeoDb(results_db, paste0(genotype_dir, "genotype/all-samples_IGHV-genotyped.tab"))

db <- seqdb %>% rename(sequence_id = SEQUENCE_ID, sequence = SEQUENCE_INPUT, productive = FUNCTIONAL,
                       vj_in_frame = IN_FRAME, stop_codon = STOP, locus = LOCUS,
                       v_call = V_CALL, d_call = D_CALL, j_call = J_CALL,
                       sequence_alignment = SEQUENCE_IMGT, v_sequence_start = V_SEQ_START,
                       v_germline_start = V_GERM_START_IMGT, np1_length = NP1_LENGTH,
                       d_sequence_start = D_SEQ_START, d_germline_start = D_GERM_START,
                       np2_length = NP2_LENGTH, j_sequence_start = J_SEQ_START,
                       j_sequence_end = J_SEQ_LENGTH, j_germline_start = J_GERM_START,
                       junction = JUNCTION, junction_length = JUNCTION_LENGTH,
                       germline_alignment = GERMLINE_IMGT, v_score = V_SCORE,
                       v_identity = V_IDENTITY, v_support = V_EVALUE, v_cigar = V_CIGAR,
                       d_score = D_SCORE, d_identity = D_IDENTITY, d_support = D_EVALUE,
                       d_cigar = D_CIGAR, j_score = J_SCORE, j_identity = J_IDENTITY,
                       j_support = J_EVALUE, j_cigar = J_CIGAR, fwr1 = FWR1_IMGT,
                       fwr2 = FWR2_IMGT, fwr3 = FWR3_IMGT, fwr4 = FWR4_IMGT,
                       cdr1 = CDR1_IMGT, cdr2 = CDR2_IMGT, cdr3 = CDR3_IMGT,
                       cell_id = CELL, c_call = C_CALL, consensus_count = CONSCOUNT,
                       patient_id = PATIENT_ID, v_sequence_end = V_SEQ_LENGTH,
                       ) %>% 
  select(-c(MUTATED_INVARIANT, INDELS, SEQUENCE_VDJ, V_GERM_START_VDJ, V_GERM_LENGTH_VDJ,
            V_GERM_LENGTH_IMGT, D_SEQ_LENGTH, D_GERM_LENGTH, J_GERM_LENGTH, UMICOUNT,
            V_CALL_10X, D_CALL_10X, J_CALL_10X, JUNCTION_10X, JUNCTION_10X_AA))
