# Libaries
library(ggplot2)
library(tidyverse)
library(reshape2)

pan <- read.csv("./data/gene_presence_absence.csv", header = TRUE)
pan

# Roary Metadata
strain_cols <- 15:ncol(pan)
strain_names <- colnames(pan)[strain_cols]
n_strains <- length(strain_names)
cat(paste("Found", n_strains, "strains\n"))

# Binary Matrix
binary_matrix <- matrix(0, nrow = nrow(pan), ncol = n_strains)

for (i in 1:n_strains) {
  binary_matrix[, i] <- ifelse(pan[, strain_cols[i]] != "" & !is.na(pan[, strain_cols[i]]), 1, 0)
}
colnames(binary_matrix) <- strain_names
binary_matrix %>% View()

# Fucntion to calculate Pangenome
calc_pancore_growth <- function(strain_order) {
  pan_sizes <- numeric(n_strains)
  core_sizes <- numeric(n_strains)
  
  for (i in 1:n_strains) {
    # Get subset of strains up to position i
    current_strains <- strain_order[1:i]
    current_matrix <- binary_matrix[, current_strains, drop = FALSE]
    
    # Pangenome: genes present in at least one strain
    pan_sizes[i] <- sum(rowSums(current_matrix) > 0)
    
    # Core genome: genes present in ALL strains
    core_sizes[i] <- sum(rowSums(current_matrix) == i)
    
  }
  
  return(list(pan = pan_sizes, core = core_sizes))
}

# Calculate growth curve
n_shuffles <- 100
all_pan <- matrix(0, nrow = n_strains, ncol = n_shuffles)
all_core <- matrix(0, nrow = n_strains, ncol = n_shuffles)

for (shuffle in 1:n_shuffles) {
  strain_order <- sample(strain_names)
  growth <- calc_pancore_growth(strain_order)
  
  all_pan[, shuffle] <- growth$pan
  all_core[, shuffle] <- growth$core
}

# Take means across shuffles
growth <- list(
  pan = rowMeans(all_pan),
  core = rowMeans(all_core)
)

# Prepare data for plotting
plot_data <- data.frame(
  Genomes = 1:n_strains,
  Pangenome = growth$pan,
  Core = growth$core
)

# Pangenome Growth Curve
p1 <- ggplot(plot_data, aes(x = Genomes, y = Pangenome)) +
  geom_line(color = "steelblue", size = 2) +
  geom_point(color = "steelblue", size = 2.5, alpha = 0.7) +
  labs(
    x = "Number of Genomes",
    y = "Pangenome Size (gene families)",
    title = "Pangenome Growth Curve"
  ) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
p1

# Core Genome Concentration Curve
p2 <- ggplot(plot_data, aes(x = Genomes, y = Core)) +
  geom_line(color = "firebrick", size = 1) +
  geom_point(color = "firebrick", size = 1.5, alpha = 0.7) +
  labs(
    x = "Number of Genomes",
    y = "Core Genome Size (gene families)",
    title = "Core Genome Contraction Curve"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

p2

# Combine plot
plot_data_long <- melt(plot_data, id.vars = "Genomes",
                       variable.name = "Type", value.name = "Size")

p3 <- ggplot(plot_data_long, aes(x = Genomes, y = Size, color = Type)) +
  geom_line(size = 1) +
  geom_point(size = 1.5, alpha = 0.7) +
  scale_color_manual(values = c("Pangenome" = "steelblue", "Core" = "firebrick"),
                     labels = c("Core Genome", "Pangenome")) +
  labs(
    x = "Number of Genomes",
    y = "Gene Families",
    title = "Pangenome vs Core Genome Dynamics",
    color = ""
  ) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),       # Centers and enlarges the title
        axis.title.x = element_text(size = 14, face = "bold"),                  # Enlarges x-axis label
        axis.title.y = element_text(size = 14, face = "bold"),                  # Enlarges y-axis label
        axis.text.x = element_text(size = 15),                                  # Enlarges x-axis tick labels
        axis.text.y = element_text(size = 15),
        legend.position = "top")
p3

ggsave(file ='./result/pan_v_core.png', p3, height = 6, width = 8, dpi = 600)
