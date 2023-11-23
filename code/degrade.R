
# LOAD LIBRARIES ----------------------------------------------------------

# Load libraries for handling files
library(tidyverse)
library(reshape2)
# Load libraries for plotting
library(ggplot2)
library(ggpubr)
library(ComplexHeatmap)
# Load libraries for clustering and visualisation
library(cluster)    
library(factoextra) 
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
HpaConsensus_df <- dcast(HpaConsensus, Gene_name ~ Tissue, value.var = "nTPM", fun.aggregate = mean, na.rm = TRUE)
# Remove rows containing NA values
HpaConsensus_df <- HpaConsensus_df %>% 
  drop_na()
# Remove lowly expressed genes
HpaConsensus_df <- HpaConsensus_df %>% 
  .[-which(rowSums(.[sapply(., is.numeric)]) < 100),]
# Sanity check
'IGF2R' %in% HpaConsensus_df$Gene_name
# Set gene names to rownames, as a lot of functions expects a matrix
rownames(HpaConsensus_df) <- HpaConsensus_df$Gene_name
HpaConsensus_df <- subset(HpaConsensus_df, select=-Gene_name)

# EXPLORE IGF2R -----------------------------------------------------------

# Plot barplot of IGF2R expression across tissues
ggplot(data=HpaConsensus[HpaConsensus$Gene_name=="IGF2R",], aes(x=reorder(Tissue, -nTPM), y=nTPM)) +
  geom_bar(stat="identity", fill="steelblue") +
  ylab("Normalized transcript expression") +
  xlab("Tissue") +
  theme_pubr() +
  guides(x =  guide_axis(angle = 90))
ggsave("plots/IGF2R_expression.pdf", device = "pdf")

# CALCULATE CORRELATION -----------------------------------------------------------------

# Save expression of the gene of interest to a separate variable
Igf2r_expression <- HpaConsensus_df["IGF2R",]
# Calculate the correlation vector between IGF2R expression across tissues and all other genes in the dataset
corr_vector <- Igf2r_expression %>% 
  t() %>% 
  cor(., t(HpaConsensus_df[rownames(HpaConsensus_df)!="IGF2R",]), method = "spearman") %>% # Don't correlate IGF2R with itself
  t() %>% as.data.frame()
# Order the correlation vector by value
corr_vector <- corr_vector[order(-corr_vector$IGF2R, decreasing = FALSE), , drop = FALSE]
# Print top 10 correlated genes
head(corr_vector)

# SUBSET AND SCALE -----------------------------------------------------------------

# Keep only genes with high (>0.7) correlation
corr_vector_high <- corr_vector[corr_vector>0.7, ,drop=FALSE]
# Save the highly correlated genes genes to a csv file
write.csv(corr_vector_high, "output/corr_spear_07.csv", row.names=TRUE)
# Select the most highly correlated genes (including IGF2R)
selected_genes <- c("IGF2R", rownames(corr_vector_high))
# Subset the main dataset based on high correlation with IGF2R
HpaConsensus_df_subset <- as.matrix(HpaConsensus_df[rownames(HpaConsensus_df) %in% selected_genes,])
scaled_mat = t(scale(t(HpaConsensus_df_subset), center=TRUE))

# K-MEANS CLUSTERING ------------------------------------------------------

# Compute gap statistic to check for optimal number of clusters
set.seed(2023) # Set the seed before each run to create reproducible results
gap_stat <- clusGap(scaled_mat, FUN = kmeans, nstart = 25,
                    K.max = 10, B = 50)
# Visualise the gap statstic results
set.seed(2023) 
fviz_gap_stat(gap_stat) # k = 4
ggsave("plots/gap_stat.pdf", device = "pdf")
# Compute k-means clustering with the optimal k = 4
set.seed(2023)
k4 <- kmeans(scaled_mat, centers = 4, nstart = 25)
# Visualise the clustering results on a PCA plot
set.seed(2023)
fviz_cluster(k4, data = scaled_mat)
ggsave("plots/clustering_k4.pdf", device = "pdf")
# Extract the number of the cluster containg IGF2R
igf2r_clust <- k4$cluster[[which(names(k4$cluster)=='IGF2R')]]
# Get the list of genes in the same cluster as IGF2R
target_genes <- names(k4$cluster[k4$cluster==1])
# Save the list of target genes
write.csv(target_genes, "output/target_genes.csv", row.names=TRUE)




