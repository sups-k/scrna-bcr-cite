# Based on the function `occupiedscRepertoire()` from the `scRepertoire` package

# Load library
library(tidyverse)

# Read BCR-seq data
bcr <- read_csv(file = "~/Downloads/Pillai_B1/scRNA-PI3Kdelta/BCR-seq/bcr_annot.csv") %>%
  as_tibble() %>%
  dplyr::filter(integrated_snn_res.0.4 <= 18)
# Removing clusters 19 and 20 because they are contamination


# Function to create a table for the stacked bar chart. It will indicate
# how many clones fall into these categories for each cluster: hyper-expanded,
# large, medium, small, and single (not expanded).
# In the function, x refers to the patient ID.

occupied_clonal_repertoire <- function(x, ID){
  # Empty table to be filled with clone type information for each cluster
  db_full <- NULL
  
  # Clone Types:
  # Hyper-expanded is 100 < freq <= 500
  # Large is 20 < freq <= 100
  # Medium is 5 < freq <= 20
  # Small is 1 < freq <= 5
  # Single is freq = 1
  cols <- c(Single = NA_integer_, Small = NA_integer_, Medium = NA_integer_, Large = NA_integer_, Hyperexpanded = NA_integer_)
  
  
  # Create table for bar chart
  for (i in c(0:18)) {
    # Split by cluster
    clu <- x %>% filter(integrated_snn_res.0.4 == i)
    
    # Frequency of occurrence of clones
    freq <- stack(sort(table(drop_na(clu)$CLONE), decreasing = TRUE)) %>%
      dplyr::rename(frequency = values, clone = ind)
    
    # Add column indicating the clone type
    freq <- freq %>% dplyr::mutate(type = case_when(frequency == 1 ~ "Single",
                                                    frequency <= 5 ~ "Small",
                                                    frequency <= 20 ~ "Medium",
                                                    frequency <= 100 ~ "Large",
                                                    frequency <= 500 ~ "Hyperexpanded"))
    
    # Create new table containing the frequency of each clone type
    temp <- stack(sort(table(freq$type), decreasing = TRUE))
    temp <- t(temp) # Transpose the table
    colnames(temp) <- temp[2,] # Make the last row the column name for the table
    rownames(temp) <- NULL # Remove row names
    temp <- as_tibble(temp) # Turn into tibble for removing the last row
    temp <- temp[1,] # Remove last row since it contains column names
    
    # To ensure number of columns stays the same throughout, add a missing column
    temp <- temp %>% add_column(!!!cols[!names(cols) %in% names(.)])
    # Replace NA with 0
    temp <- temp %>% 
      dplyr::mutate(across(where(anyNA), ~ replace_na(., 0)))
    
    # Make all values in the column numeric. They are in character form.
    temp$Single <- as.numeric(temp$Single)
    temp$Small <- as.numeric(temp$Small)
    temp$Medium <- as.numeric(temp$Medium)
    temp$Large <- as.numeric(temp$Large)
    temp$Hyperexpanded <- as.numeric(temp$Hyperexpanded)
    
    # Add cluster identity
    temp <- temp %>% dplyr::mutate(cluster = i)
    
    # Combine data with full database
    if (is.null(db_full)) {
      db_full <- temp
    } else {
      db_full <- rbind(db_full, temp)
    } 
  }
  
  # Remove unnecessary variables
  rm(x, clu, freq, temp, i)
  
  # Create stacked bar plot
  db_full <- as.data.frame(db_full)
  rownames(db_full) <- db_full$cluster
  barplot(t(as.matrix(db_full[,1:5])),
          legend.text = c("Single", "Small", "Medium", "Large", "Hyperexpanded"),
          args.legend = list(x = "topright", inset=c(-0.2,-0.1)),
          main = ID,
          xlab = "Cluster",
          ylab = "Frequency",
          col = c("#253494", "#2C7FB8", "#41B6C4", "#A1DAB4", "#FFFFCC"))
}


# Split by patient - can plot this only on a patient-by-patient basis
a <- bcr %>% filter(patient_ID == "Hashtag8")
occupied_clonal_repertoire(a, "Hashtag 8: Patient")

