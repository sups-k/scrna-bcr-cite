library(tidyverse)

# Read metadata from Seurat object as tibble
meta <- read.csv(file = "../../KW9079_Farmer/210624_A00873_0376_BHC3KTDSX2_KW9079/cellranger-6.0.1/GRCh38/results/10_pc43_11March_METADATA.csv") %>% 
  as_tibble()

# Read functional makedb output
h1 <- read_delim(file = "../../../../../../../../../../home/suk571/immcantation-python3.8.12/imm-env/makedb_output/functional/1_functional_parse-select.tab")

# Subset metadata for sample read above & match column name of cell to makedb output
h1_m <- meta %>%
  filter(patient_ID == "Hashtag1") %>%
  rename(CELL = cells)

# Merge both files
merged_h1 <- merge(h1_m, h1, by = "CELL")

# Keep original column names
merged_h1 <- merged_h1 %>% select(colnames(h1))

# Save to functional
write_tsv(merged_h1, file = "../../../../../../../../../../home/suk571/immcantation-python3.8.12/imm-env/makedb_output/functional/filtered_from_Seurat/1_functional_parse-select.tab")
rm(h1, h1_m, merged_h1)