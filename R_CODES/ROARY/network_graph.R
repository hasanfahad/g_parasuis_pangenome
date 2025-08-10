# Libaries
library(ggplot2)
library(tidyverse)
library(reshape2)
library(igraph)
library(ggraph)
library(pheatmap)

# Load data
gpa <- read.csv("./data/gene_presence_absence.csv", stringsAsFactors = FALSE)

# Extract strain column
strain_cols <- 15:ncol(gpa)
strain_names <- colnames(gpa)[strain_cols]
gene_names <- gpa$Gene

# Creating a binary matrix
binary_matrix <- matrix(0, nrow = nrow(gpa), ncol = length(strain_names))

for (i in 1:length(strain_names)) {
  binary_matrix[, i] <- ifelse(gpa[, strain_cols[i]] != "" & !is.na(gpa[, strain_cols[i]]), 1, 0)
}
colnames(binary_matrix) <- strain_names
rownames(binary_matrix) <- gene_names

# Jaccard Index
calculate_jaccard <- function(x, y) {
  intersection <- sum(x & y)
  union <- sum(x | y)
  return(intersection / union)
}

# Create similarity index between genomes
n_strains <- length(strain_names)
similarity_matrix <- matrix(0, nrow = n_strains, ncol = n_strains)
rownames(similarity_matrix) <- strain_names
colnames(similarity_matrix) <- strain_names

for (i in 1:n_strains) {
  for (j in i:n_strains) {
    sim <- calculate_jaccard(binary_matrix[, i], binary_matrix[, j])
    similarity_matrix[i, j] <- sim
    similarity_matrix[j, i] <- sim
  }
}

#-------------------Network Graph-----------------------------------------------

# Filter edges by similarity threshold
threshold <- 0.7

edge_list <- data.frame()

for (i in 1:(n_strains-1)) {
  for (j in (i+1):n_strains) {
    if (similarity_matrix[i, j] >= threshold) {
      edge_list <- rbind(edge_list, data.frame(
        from = strain_names[i],
        to = strain_names[j],
        weight = similarity_matrix[i, j]
      ))
    }
  }
}

# Create igraph object
if (nrow(edge_list) > 0) {
  g <- graph_from_data_frame(edge_list, directed = FALSE, vertices = strain_names)
} else {
  g <- make_empty_graph(n = n_strains)
  V(g)$name <- strain_names
}

# Calculate genome specific statistics
genome_stats <- data.frame(
  strain = strain_names,
  total_genes = colSums(binary_matrix),
  unique_genes = apply(binary_matrix, 2, function(x) sum(x == 1 & rowSums(binary_matrix) == 1)),
  core_genes = apply(binary_matrix, 2, function(x) sum(x == 1 & rowSums(binary_matrix) == n_strains))
)

# Add attribute to graph
V(g)$total_genes <- genome_stats$total_genes[match(V(g)$name, genome_stats$strain)]
V(g)$unique_genes <- genome_stats$unique_genes[match(V(g)$name, genome_stats$strain)]

# Network plot
p1 <- ggraph(g, layout = "fr") +  # Fruchterman-Reingold layout
  geom_edge_link(aes(width = weight), alpha = 0.6, color = "gray70") +
  geom_node_point(aes(size = total_genes, color = unique_genes), alpha = 0.8) +
  #geom_node_text(aes(label = name), size = 2, repel = TRUE) +
  scale_size_continuous(name = "Total Genes", range = c(2, 8)) +
  scale_color_gradient(name = "Unique Genes", low = "blue", high = "red") +
  scale_edge_width(name = "Similarity", range = c(0.5, 3)) +
  theme_graph() +
  labs(title = "Pangenome Network Graph",
       subtitle = paste("Similarity threshold:", threshold)) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))
p1

ggsave(file ='./result/network_plot.png', p1, height = 6, width = 8, dpi = 600)

# Circular plot
p2 <- ggraph(g, layout = "circle") +
  geom_edge_arc(aes(width = weight), alpha = 0.3, color = "gray70") +
  geom_node_point(aes(size = total_genes, color = unique_genes), alpha = 0.8) +
  #geom_node_text(aes(label = name), size = 2, repel = TRUE) +
  scale_size_continuous(name = "Total Genes", range = c(2, 8)) +
  scale_color_gradient(name = "Unique Genes", low = "blue", high = "red") +
  scale_edge_width(name = "Similarity", range = c(0.5, 3)) +
  theme_graph() +
  labs(title = "Pangenome Network Graph (Circular Layout)",
       subtitle = paste("Similarity threshold:", threshold)) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

p2

ggsave(file ='./result/circular_plot.png', p2, height = 8, width = 8, dpi = 600)

# Heatmap

heatmap_data <- melt(similarity_matrix)
colnames(heatmap_data) <- c("Genome1", "Genome2", "Similarity")

p3 <- ggplot(heatmap_data, aes(x = Genome1, y = Genome2, fill = Similarity)) +
  geom_tile(color = "white", size = 0.1) +
  scale_fill_viridis_c(name = "Jaccard\nSimilarity", option = "plasma") +
  labs(
    title = "Genome Similarity Heatmap",
    subtitle = paste("Based on", nrow(gpa), "gene families"),
    x = "Genomes",
    y = "Genomes"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
    axis.text.y = element_text(size = 6),
    plot.title = element_text(hjust = 0.5, size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 10),
    panel.grid = element_blank()
  ) +
  coord_fixed()

p3
