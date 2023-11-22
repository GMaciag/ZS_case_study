
# LOAD LIBRARIES ----------------------------------------------------------

# Load tidyverse for handling files
library(tidyverse)
# Load ggplot2 for plotting
library(ggplot2)
# Load pheatmap to plot heatmaps
library(pheatmap)
# The following setting is important, do not omit
options(stringsAsFactors = FALSE)

# LOAD DATA ---------------------------------------------------------------

# Load the data.
HpaConsensus <- read_tsv("../resources/rna_tissue_consensus.tsv")

# PREPARE THE DATA --------------------------------------------------------

# Print the data table to the console
HpaConsensus
# Remove the unnecessary data field
HpaConsensus <- HpaConsensus[, -1]
# Rename a column to remove whitespace
HpaConsensus <- HpaConsensus %>% rename(Gene_name = `Gene name`)
# Reshape from long to wide
HpaConsensus_wide <- HpaConsensus %>% 
  pivot_wider(names_from = Tissue, values_from = nTPM, values_fn = list)
# Remove rows containing NA values
HpaConsensus_wide_noNA <- HpaConsensus_wide %>% 
  drop_na()
# Sanity check
'IGF2R' %in% HpaConsensus_wide_noNA$Gene_name





