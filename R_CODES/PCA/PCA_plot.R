# Load library
library(ggplot2)
library(colorspace)
library(tidyverse)
library(factoextra)

# Load data
data <- read.csv('./data/gene_presence_absence.csv', header =TRUE, na.strings = c('', 'NA'))
data %>% View()

# Remove metadata columns to get the binary matrix
columns_to_remove <- c("Gene","Non.unique.Gene.name","No..isolates","No..sequences","Avg.sequences.per.isolate",
                       "Genome.Fragment","Order.within.Fragment", "Accessory.Fragment","Accessory.Order.with.Fragment",
                       "QC","Min.group.size.nuc","Max.group.size.nuc","Avg.group.size.nuc")
pa_data <- data[ , !(names(data) %in% columns_to_remove) ]
pa_data %>% View()

# Remove gene annotation columns if any are still present
rownames(pa_data) <- data$Gene  # Use gene names as row names if desired
# Change NAs to 0 (absent), other values to 1 (present)
pa_bin <- as.matrix(pa_data)
pa_bin[is.na(pa_bin)] <- 0
pa_bin[pa_bin != 0] <- 1
# Some values may be factors/strings—convert to numeric
pa_bin <- apply(pa_bin, 2, as.numeric)
rownames(pa_bin) <- data$Gene

# Transpose if your matrix is genes x samples, as PCA needs samples x features
pca_input <- t(pa_bin)
pca_res <- prcomp(pca_input, center = TRUE, scale. = FALSE)

plot(pca_res$x[,1], pca_res$x[,2], xlab = "PC1", ylab = "PC2",
     main = "PCA of Roary Gene Presence/Absence Matrix")
text(pca_res$x[,1], pca_res$x[,2], labels = rownames(pca_res$x), cex=0.7, pos=3)

# Prepare data for ggplot
df <- as.data.frame(pca_res$x)
df$Strain <- rownames(df)

pca_ggplot <- ggplot(df, aes(x = PC1, y = PC2, label = Strain)) +
  geom_point(size = 2.5, alpha = 0.3, color = 'purple') +
  labs(
    title = "PCA of Gene Presence/Absence",
    x = paste0("PC1 (", round(summary(pca_res)$importance[2, 1] * 100, 1), "% variance)"),
    y = paste0("PC2 (", round(summary(pca_res)$importance[2, 2] * 100, 1), "% variance)"),
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),       # Centers and enlarges the title
    axis.title.x = element_text(size = 14, face = "bold"),                  # Enlarges x-axis label
    axis.title.y = element_text(size = 14, face = "bold"),                  # Enlarges y-axis label
    axis.text.x = element_text(size = 15),                                  # Enlarges x-axis tick labels
    axis.text.y = element_text(size = 15)                                   # Enlarges y-axis tick labels
  ) +
  scale_color_brewer(palette = "Set1")

ggsave(file = './result/pca_plot.png', pca_ggplot, height = 6, width = 6, dpi = 600)
