# Load library
library(pheatmap)
library(ComplexHeatmap)
library(RColorBrewer)
library(tidyverse)
library(reshape2)
library(viridis)

# Load data
ani_table <- read_tsv('./data/all_vs_all_fastani_output.txt', col_names = FALSE)
colnames(ani_table) <- c("query", "reference", "ani", "frags_mapped", "frags_total")
ani_table %>% View()

# QC
ani_table$ani <- as.numeric(trimws(ani_table$ani))
ani_table$frags_mapped <- as.numeric(trimws(ani_table$frags_mapped))
ani_table$frags_total <- as.numeric(trimws(ani_table$frags_total))

ani_table$query <- gsub('assembly_list/', '', ani_table$query)
ani_table$query <- gsub('_genomic.fna', '', ani_table$query)
ani_table$query <- gsub('.fna', '', ani_table$query)
ani_table$reference <- gsub('assembly_list/', '', ani_table$reference)
ani_table$reference <- gsub('_genomic.fna', '', ani_table$reference)
ani_table$reference <- gsub('.fna', '', ani_table$reference)

# Get all unique strain values
all_strain <- sort(unique(c(ani_table$query, ani_table$reference)))

# Get a square matrix
ani_matrix <- matrix(NA, nrow = length(all_strain), ncol = length(all_strain), dimnames = list(all_strain, all_strain))

# Fill the matrix
for (i in 1:nrow(ani_table)) {
  row <- ani_table[i, ]
  ani_matrix[row$query, row$reference] <- row$ani
  ani_matrix[row$reference, row$query] <- row$ani  # symmetric fill
}

ani_matrix %>% View()

# Hierarchical clustering on both axes
row_dend <- hclust(dist(ani_matrix))
col_dend <- hclust(dist(t(ani_matrix)))

# Order indices
row_order <- row_dend$order
col_order <- col_dend$order

ani_matrix <- ani_matrix[row_order, col_order]

# Create heatmap dataset
ani_long <- melt(ani_matrix, varnames = c("query", "reference"), value.name = "ANI")

# Plot heatmap
ani_plot <- ggplot(ani_long, aes(x = query, y = reference, fill = ANI)) +
  geom_tile(color = "white") +  # only colored tiles, no text
  scale_fill_viridis(
    option = "plasma",
    name = "ANI (%)",
    limits = c(96, 100)
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "right"
  )

ani_plot

ggsave(file = './result/ani_plot.png', ani_plot, height = 6, width = 8, dpi = 300)






