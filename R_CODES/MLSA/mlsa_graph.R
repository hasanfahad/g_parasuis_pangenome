# Load library
library(ape)
library(ggtree)
library(tidyverse)
library(ggplot2)
library(rentrez)

tree <- read.tree('./data/mlsa_supermatrix.fas.treefile')
tree$tip.label <- gsub("\\.\\d+$", "", tree$tip.label)
tree$tip.label %>% as_data_frame() %>% View()

p <- ggtree(tree, layout = "circular", branch.length = "none", size = 0.6, color = "#2c3e50") +
  theme_tree2() +
  ggtitle("Circular MLSA Phylogenetic Tree")

p

# Phyla
accessions <- tree$tip.label
accessions

lapply(accessions, function(acc) {
  summary <- entrez_summary(db="assembly", id=acc)
  taxid <- summary$taxid
  taxo_info <- entrez_summary(db="taxonomy", id=taxid)
  taxo_info$scientificname
  # For extended taxonomy, use entrez_fetch for full lineage
})
