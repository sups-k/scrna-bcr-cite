library(alakazam)
library(shazam)
library(tigger)
library(scoper)
library(future)

plan("multicore")

genotype_dir <- "~/immcantation-python3.8.12/imm-env/makedb_output/functional/"
patient_ID <- c("Hashtag1", "Hashtag2", "Hashtag3", "Hashtag4", "Hashtag5", "Hashtag6", "Hashtag7", "Hashtag8")
predicted_thresholds <- data.frame(patient_ID = patient_ID, threshold=rep(0.1, length(patient_ID)))
seqdb <- NULL
igv <- readIgFasta("~/immcantation-python3.8.12/imm-env/imgt/human/vdj/imgt_human_IGV.fasta")

# For healthy 1
tmp <- readChangeoDb("~/immcantation-python3.8.12/imm-env/makedb_output/functional/1_functional_parse-select.tab")
tmp$PATIENT_ID <- patient_ID[1]
seqdb <- rbind(seqdb, tmp)
cat('\n\tSuccessfully read file 1', '\n')
# For patient 2
tmp <- readChangeoDb("~/immcantation-python3.8.12/imm-env/makedb_output/functional/2_functional_parse-select.tab")
tmp$PATIENT_ID <- patient_ID[2]
seqdb <- rbind(seqdb, tmp)
cat('\n\tSuccessfully read file 2', '\n')
# For healthy 3
tmp <- readChangeoDb("~/immcantation-python3.8.12/imm-env/makedb_output/functional/3_functional_parse-select.tab")
tmp$PATIENT_ID <- patient_ID[3]
seqdb <- rbind(seqdb, tmp)
cat('\n\tSuccessfully read file 3', '\n')
# For patient 4
tmp <- readChangeoDb("~/immcantation-python3.8.12/imm-env/makedb_output/functional/4_functional_parse-select.tab")
tmp$PATIENT_ID <- patient_ID[4]
seqdb <- rbind(seqdb, tmp)
cat('\n\tSuccessfully read file 4', '\n')
# For healthy 5
tmp <- readChangeoDb("~/immcantation-python3.8.12/imm-env/makedb_output/functional/5_functional_parse-select.tab")
tmp$PATIENT_ID <- patient_ID[5]
seqdb <- rbind(seqdb, tmp)
cat('\n\tSuccessfully read file 5', '\n')
# For healthy 6
tmp <- readChangeoDb("~/immcantation-python3.8.12/imm-env/makedb_output/functional/6_functional_parse-select.tab")
tmp$PATIENT_ID <- patient_ID[6]
seqdb <- rbind(seqdb, tmp)
cat('\n\tSuccessfully read file 6', '\n')
# For patient 7
tmp <- readChangeoDb("~/immcantation-python3.8.12/imm-env/makedb_output/functional/7_functional_parse-select.tab")
tmp$PATIENT_ID <- patient_ID[7]
seqdb <- rbind(seqdb, tmp)
cat('\n\tSuccessfully read file 7', '\n')
# For patient 8
tmp <- readChangeoDb("~/immcantation-python3.8.12/imm-env/makedb_output/functional/8_functional_parse-select.tab")
tmp$PATIENT_ID <- patient_ID[8]
seqdb <- rbind(seqdb, tmp)
cat('\n\tSuccessfully read file 8', '\n')

rm(tmp)
gc()

for(i in patient_ID){
  seqdb_patient <- seqdb[seqdb$PATIENT_ID %in% i, ]
  
  novel_rows <- NULL
  nv <- NULL
  nv <- findNovelAlleles(data = seqdb_patient,
                         germline_db = igv,
                         seq = "SEQUENCE_IMGT",
                         v_call = "V_CALL",
                         j_call = "J_CALL",
                         junction = "JUNCTION",
                         junction_length = "JUNCTION_LENGTH",
                         nproc = 4)
  
  
  # Infer patient genotype
  gt_pat <- inferGenotype(data = seqdb_patient,
                          germline_db = igv,
                          novel = nv,
                          v_call = "V_CALL",
                          seq = "SEQUENCE_IMGT")
  
  # Convert genotype table to vector of nucleotide sequences
  gtseq_pat <- genotypeFasta(genotype = gt_pat, germline_db = igv, novel=nv)

  # Correct allele calls based on the personalized genotype
  seqdb_patient <- reassignAlleles(data = seqdb_patient,
                                   genotype_db = gtseq_pat,
                                   v_call = "V_CALL",
                                   seq = "SEQUENCE_IMGT")
  
  # Running spectralClones since findThreshold gave unimodal distribution
  results <- spectralClones(db = seqdb_patient,
                            method = "novj",
                            germline = "GERMLINE_IMGT",
                            sequence = "SEQUENCE_IMGT",
                            junction = "JUNCTION",
                            v_call = "V_CALL",
                            j_call = "J_CALL",
                            clone = "CLONE",
                            locus = "LOCUS",
                            cell_id = "SEQUENCE_ID",
                            split_light = TRUE,
                            only_heavy = FALSE,
                            nproc = 4)

  # Get results data frame
  results_db <- as.data.frame(results)

  # Add threshold to table - change for each patient
  predicted_thresholds$threshold[predicted_thresholds$patient_ID == i] <- results@eff_threshold
  
  # Save data
  writeChangeoDb(results_db, paste0(genotype_dir, "genotype/", i, "_IGV-genotyped.tab"))
}

write.csv(predicted_thresholds, file = "~/immcantation-python3.8.12/imm-env/makedb_output/functional/predicted_thresholds.csv", quote = FALSE, row.names = FALSE)

