
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
# Load package for GO and KEGG analysis
library(limma)
library(org.Hs.eg.db)
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
HpaConsensus <- HpaConsensus %>% dplyr::rename(Gene_name = `Gene name`)
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

# Scale and center
HpaConsensus_scaled_mat = HpaConsensus_df %>% 
  t() %>% 
  scale(center = TRUE) %>% 
  t()

# EXPLORE IGF2R -----------------------------------------------------------

# Plot barplot of IGF2R expression across tissues
ggplot(data=HpaConsensus[HpaConsensus$Gene_name=="IGF2R",], aes(x=reorder(Tissue, -nTPM), y=nTPM)) +
  geom_bar(stat="identity", fill="steelblue") +
  ylab("Normalized transcript expression") +
  xlab("Tissue") +
  theme_pubr() +
  guides(x =  guide_axis(angle = 90))
ggsave("output_plots/IGF2R_expression.pdf", device = "pdf")

# CALCULATE CORRELATION AND SUBSET DATASET -----------------------------------------------------------------

# Save expression of the gene of interest to a separate variable
Igf2r_expression <- HpaConsensus_scaled_mat["IGF2R",]
# Calculate the correlation vector between IGF2R expression across tissues and all other genes in the dataset
corr_vector <- Igf2r_expression %>% 
  cor(., t(HpaConsensus_scaled_mat[rownames(HpaConsensus_scaled_mat)!="IGF2R",]), method = "spearman") %>% # Don't correlate IGF2R with itself
  t() %>% as.data.frame() %>% setNames("IGF2R correlation")
# Order the correlation vector by value
corr_vector <- corr_vector[order(-corr_vector$IGF2R, decreasing = FALSE), , drop = FALSE]
# Print top 10 correlated genes
head(corr_vector)

# Keep only genes with high (>0.7) correlation
corr_vector_high <- corr_vector[corr_vector>0.7, ,drop=FALSE]
# Save the highly correlated genes genes to a csv file
write.csv(corr_vector_high, "output_files/corr_spear_07.csv", row.names=TRUE)
# Select the most highly correlated genes (including IGF2R)
select_genes_corr_high <- c("IGF2R", rownames(corr_vector_high))
# Subset the main dataset based on high correlation with IGF2R
HpaConsensus_scaled_mat_corr_high <- HpaConsensus_scaled_mat[rownames(HpaConsensus_scaled_mat) %in% select_genes_corr_high,]

# K-MEANS CLUSTERING ------------------------------------------------------

# Compute gap statistic to check for optimal number of clusters
set.seed(2023) # Set the seed before each run to create reproducible results
gap_stat <- clusGap(HpaConsensus_scaled_mat_corr_high, FUN = kmeans, nstart = 25,
                    K.max = 10, B = 50)
# Visualise the gap statstic results
set.seed(2023) 
fviz_gap_stat(gap_stat) # k = 4
ggsave("output_plots/gap_stat.pdf", device = "pdf")

# Compute k-means clustering with the optimal k = 4
set.seed(2023)
k4 <- kmeans(HpaConsensus_scaled_mat_corr_high, centers = 4, nstart = 25)
# Visualise the clustering results on a PCA plot
set.seed(2023)
fviz_cluster(k4, data = HpaConsensus_scaled_mat_corr_high)
ggsave("output_plots/clustering_k4.pdf", device = "pdf")

# Extract the number of the cluster containing IGF2R
igf2r_clust <- k4$cluster[[which(names(k4$cluster)=='IGF2R')]]
# Get the list of genes in the same cluster as IGF2R
target_genes <- names(k4$cluster[k4$cluster==1])
# Save the list of target genes
target_genes %>% 
  as.data.frame() %>% 
  setNames("IGF2R target genes") %>% 
  write.csv("output_files/target_genes.csv", row.names=TRUE)

# HEATMAP -----------------------------------------------------------------

# Subset the scaled matrix to include target genes and IGF2R
target_genes_scaled_mat <- HpaConsensus_scaled_mat_corr_high[rownames(HpaConsensus_scaled_mat_corr_high) %in% c("IGF2R", target_genes),]

# Create grouping of rows to highlight IGF2R on the heatmap
which_row_igf2r = which(grepl("IGF2R", rownames(target_genes_scaled_mat))) # Which row is IGF2R in
# Create a split of scaled_mat_target_genes rows that singles out IGF2R
split = data.frame(x = c(rep("A", which_row_igf2r - 1), "B",
                         rep("C", nrow(target_genes_scaled_mat) - which_row_igf2r)))

# Plot a heatmap of the target genes and cluster tissues 
set.seed(2023)
ht <- Heatmap(target_genes_scaled_mat, 
              show_row_names = TRUE,
              cluster_rows = FALSE,
              row_split = split,
              row_title = NULL,
              clustering_distance_columns='spearman',
              name = "Z-score", 
              column_km = 4, # k-means clustering 
              column_km_repeats=50)
draw(ht)
# Save the heatmap to file
pdf("output_plots/target_genes_heatmap.pdf", width=9, height=7)
draw(ht)
dev.off()

# GO AND KEGG ANALYSIS ----------------------------------------------------

# Take all of the genes in the original dataset as the universe for enrichment testing
universe <- rownames(HpaConsensus_df)
# Convert gene symbols to Entrez Gene IDs
target_genes_entrez <- mapIds(org.Hs.eg.db, keys = target_genes, keytype = "SYMBOL", column = "ENTREZID")
universe <- mapIds(org.Hs.eg.db, keys = universe, keytype = "SYMBOL", column = "ENTREZID")
# Remove any NA values (genes not found in the database)
target_genes_entrez <- target_genes_entrez[!is.na(target_genes_entrez)]
universe <- universe[!is.na(universe)]

# Run GO enrichment analysis
GOtest <- goana(de=target_genes_entrez, universe=universe, species="Hs")
# Summarise GO analysis results
GOtopten <- GOtest %>% 
  rownames_to_column('ID') %>% 
  filter(N<2000) %>% 
  column_to_rownames('ID') %>% 
  topGO(ontology="BP", number=10)
# Save top 10 GO results to file
write.csv(file = "output_files/top10_GO.csv", x = GOtopten, quote = FALSE)

# Run KEGG pathway analysis
KEGGtest <- kegga(de = target_genes_entrez, universe = universe, species = "Hs")
# Summarise GO analysis results
KEGGtopten <- KEGGtest %>% 
  rownames_to_column('ID') %>% 
  filter(N<2000) %>% 
  column_to_rownames('ID') %>% 
  topKEGG(number=10)
# Save top 10 GO results to file
write.csv(file = "output_files/top10_KEGG.csv", x = KEGGtopten, quote = FALSE)

# BONUS QUESTION 1 --------------------------------------------------------

# Does your analysis give any evidence to whether or not we can expect to degrade any of these proteins? 
relevant_genes <- c(
  'PCSK9',
  'TARDBP',
  'UCP2',
  'DCN',
  'APOD'
)
# Are those genes among the 75 highly correlated genes?
relevant_genes %in% rownames(corr_vector_high) # No

# Subset the scaled matrix to include the relevant genes and IGF2R
relevant_genes_scaled_mat <- HpaConsensus_scaled_mat[rownames(HpaConsensus_scaled_mat) %in% c("IGF2R", relevant_genes),]
# Plot a heatmap of the target genes and cluster tissues 
set.seed(2023)
ht <- Heatmap(relevant_genes_scaled_mat, 
              show_row_names = TRUE,
              cluster_rows = FALSE,
              clustering_distance_columns='spearman',
              name = "Z-score")
draw(ht)
# Save the heatmap to file
pdf("output_plots/relevant_genes_heatmap_bonus1.pdf", width=9, height=5)
draw(ht)
dev.off()
