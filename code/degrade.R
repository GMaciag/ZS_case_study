
# LOAD LIBRARIES ----------------------------------------------------------

# Load tidyverse for handling files
library(tidyverse)
# Load pheatmap to plot heatmaps
library(pheatmap)
# The following setting is important, do not omit
options(stringsAsFactors = FALSE)

# LOAD DATA ---------------------------------------------------------------

# Load the data.
HpaConsensus <- read_tsv("../resources/rna_tissue_consensus.tsv")

# EXPLORE THE DATA --------------------------------------------------------

# Print the data table to the console
HpaConsensus
# Remove the unnecessary data field
HpaConsensus <- HpaConsensus[, -1]
