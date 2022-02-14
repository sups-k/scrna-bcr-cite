# load library
library(tidyverse)
library(ggplot2)

# Read BCR-seq data
bcr <- read_csv(file = "~/Downloads/Pillai_B1/scRNA-PI3Kdelta/BCR-seq/bcr_annot.csv")

# Get healthy control data
healthy <- bcr %>% filter(sample == "healthy")
# Get patient data
patient <- bcr %>% filter(sample == "patient")

# Function: Make donut chart of isotypes in a cluster
donut_isotype <- function(input, title){
  # Get isotypes in descending order of frequency
  data <- stack(sort(table(drop_na(input)$C_CALL), decreasing = TRUE))
  
  # Reorder the columns
  data <- data %>% select(ind, values)
  
  # Rename the columns
  data <- data %>% rename(Isotype = ind, Frequency = values)
  
  # Change "Isotype" to character from factor
  data$Isotype <- as.character(data$Isotype)
  
  # Compute percentages
  data$fraction = data$Frequency / sum(data$Frequency)
  
  # Compute the cumulative percentages (top of each rectangle)
  data$ymax = cumsum(data$fraction)
  
  # Compute the bottom of each rectangle
  data$ymin = c(0, head(data$ymax, n=-1))
  
  # Make the donut chart
  plot <- ggplot(data, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Isotype)) +
    geom_rect() +
    scale_fill_brewer(palette = "RdYlBu") + # must be a color blind friendly palette from RColorBrewer
    scale_color_brewer(palette = "RdYlBu") +
    coord_polar(theta="y") +
    xlim(c(-1, 4)) +
    theme_void() +
    theme(legend.position = "right") +
    ggtitle(title)
  
  # Save the plot
  ggsave(filename = paste0("~/Downloads/Pillai_B1/scRNA-PI3Kdelta/", str_replace_all(string = title, pattern = " ", replacement = "_"), ".pdf"), plot = plot, device = "pdf")
}

# Function: Stacked bar plot of isotypes
barplot_isotype <- function(healthy, patient, title){
  # Get isotypes in descending order of frequency
  data1 <- stack(sort(table(drop_na(healthy)$C_CALL), decreasing = TRUE))
  data2 <- stack(sort(table(drop_na(patient)$C_CALL), decreasing = TRUE))
  
  # Reorder the columns
  data1 <- data1 %>% select(ind, values)
  data2 <- data2 %>% select(ind, values)
  
  # Rename the columns
  data1 <- data1 %>% rename(Isotype = ind, Frequency = values)
  data2 <- data2 %>% rename(Isotype = ind, Frequency = values)
  
  # Change "Isotype" to character from factor
  data1$Isotype <- as.character(data1$Isotype)
  data2$Isotype <- as.character(data2$Isotype)
  
  # Add "condition" column
  data1$Condition <- rep("Healthy", nrow(data1))
  data2$Condition <- rep("Patient", nrow(data2))
  
  # Merge both data
  data <- merge(data1, data2, all = TRUE)
  
  # Make the stacked bar chart
  plot <- ggplot(data, aes(fill=Isotype, y=Frequency, x=Condition)) + 
    geom_bar(position="stack", stat="identity") +
    scale_fill_brewer(palette = "RdYlBu") +
    ggtitle(title) +
    theme_bw() +
    xlab("")
  
  # Save it
  ggsave(filename = paste0("~/Downloads/Pillai_B1/scRNA-PI3Kdelta/barplot_isotypes_", str_replace_all(string = title, pattern = " ", replacement = "_"), ".pdf"), plot = plot, device = "pdf")
}


# Plot stacked bar plots for all cluster isotypes
for(i in 0:18){
  h <- healthy %>% dplyr::filter(integrated_snn_res.0.4 == i)
  p <- patient %>% dplyr::filter(integrated_snn_res.0.4 == i)
  title <- paste0("Cluster ", i)
  barplot_isotype(h, p, title)
}

# Plot and save donut charts for all cluster isotypes
for (i in 0:18) {
  h <- healthy %>% dplyr::filter(integrated_snn_res.0.4 == i)
  p <- patient %>% dplyr::filter(integrated_snn_res.0.4 == i)
  donut_isotype(h, paste0("Isotypes of Cluster ", i, ": Healthy"))
  donut_isotype(p, paste0("Isotypes of Cluster ", i, ": Patient"))
}


# Function: Stacked bar plot of top 10 heavy-chain clones
barplot_clones <- function(healthy, patient, title){
  # Get isotypes in descending order of frequency
  data1 <- stack(sort(table(drop_na(healthy)$CLONE), decreasing = TRUE)[1:10])
  data2 <- stack(sort(table(drop_na(patient)$CLONE), decreasing = TRUE)[1:10])
  
  # Reorder the columns
  data1 <- data1 %>% select(ind, values)
  data2 <- data2 %>% select(ind, values)
  
  # Rename the columns
  data1 <- data1 %>% rename(Clone = ind, Frequency = values)
  data2 <- data2 %>% rename(Clone = ind, Frequency = values)
  
  # Change "Clone" to character from factor
  data1$Clone <- as.character(data1$Clone)
  data2$Clone <- as.character(data2$Clone)
  
  # Add "condition" column
  data1$Condition <- rep("Healthy", nrow(data1))
  data2$Condition <- rep("Patient", nrow(data2))
  
  # Merge both data
  data <- merge(data1, data2, all = TRUE)
  
  # Make the stacked bar chart
  plot <- ggplot(data, aes(fill=Clone, y=Frequency, x=Condition)) + 
    geom_bar(position="stack", stat="identity") +
    scale_fill_brewer(palette = "RdYlBu") +
    ggtitle(title) +
    theme_bw() +
    xlab("")
  
  # Save it
  ggsave(filename = paste0("~/Downloads/Pillai_B1/scRNA-PI3Kdelta/barplot_clones_", str_replace_all(string = title, pattern = " ", replacement = "_"), ".pdf"), plot = plot, device = "pdf")
}





