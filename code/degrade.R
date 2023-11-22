
# LOAD LIBRARIES ----------------------------------------------------------

# Load libraries for handling files
library(tidyverse)
library(reshape2)
# Load libraries for plotting
library(ggplot2)
library(ggpubr)
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
# Transform the data from long to wide format
HpaConsensus_wide <- dcast(HpaConsensus, Gene_name ~ Tissue, value.var = "nTPM", fun.aggregate = mean, na.rm = TRUE)
# Remove rows containing NA values
HpaConsensus_wide_noNA <- HpaConsensus_wide %>% 
  drop_na()
# Sanity check
'IGF2R' %in% HpaConsensus_wide_noNA$Gene_name

# EXPLORE IGF2R -----------------------------------------------------------

# Plot barplot of IGF2R expression across tissues
ggplot(data=HpaConsensus[HpaConsensus$Gene_name=="IGF2R",], aes(x=reorder(Tissue, -nTPM), y=nTPM)) +
  geom_bar(stat="identity", fill="steelblue") +
  ylab("Normalized transcript expression") +
  xlab("Tissue") +
  theme_pubr() +
  guides(x =  guide_axis(angle = 90))

ggsave("plots/IGF2R_expression.pdf", device = "pdf")


