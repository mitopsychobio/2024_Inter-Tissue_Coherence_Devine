library(igraph)
library(corrplot)
library(grDevices)
library(tidyverse)
library(naniar)
library(circlize)
library(ggfortify)
library(psych)


# run correlation matrix
cor_matrix <- cor(mito_nDNA_percent, use = "pairwise.complete.obs", method = "spearman")

# count sample size of each tissue pair
n_tissue_cor <- pairwiseCount(mito_nDNA_percent) %>%
  as.data.frame()

# apply sample size cutoff
n_tissue_cor <- n_tissue_cor %>%
  replace_with_na_if(.predicate = is.double, condition = ~.x <10) %>%
  `rownames<-`(colnames(n_tissue_cor))


tissue_cor_cutoff <- cor_matrix
tissue_cor_cutoff[is.na(n_tissue_cor)] <- NA


# Set a threshold for correlations to consider an edge
threshold <- 0.2

# create adjacency matrix 
adj_matrix <- abs(tissue_cor_cutoff) > threshold

# set NAs to False
adj_matrix[is.na(adj_matrix)] <- FALSE

g <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected", diag = FALSE)


edge_list <- get.edgelist(g)
weights <- tissue_cor_cutoff[as.matrix(edge_list)]


E(g)$weight <- abs(weights)


group_assignments <- list(
  'CNS' = c("Brain - Caudate (basal ganglia)", "Brain - Anterior cingulate cortex (BA24)", "Brain - Putamen (basal ganglia)", "Brain - Nucleus accumbens (basal ganglia)", 
            "Brain - Amygdala", "Brain - Substantia nigra", "Brain - Hippocampus", "Brain - Hypothalamus","Brain - Cortex", "Brain - Frontal Cortex (BA9)",
            "Brain - Cerebellum", "Brain - Cerebellar Hemisphere", "Brain - Spinal cord (cervical c-1)"),
  'Anabolic' = c("Kidney - Cortex","Liver"),
  'Contractile' = c("Heart - Atrial Appendage", "Muscle - Skeletal", "Heart - Left Ventricle"),
  'Reproductive' = c("Testis", "Vagina", "Uterus", "Ovary", "Prostate"),
  'Digestive' = c("Esophagus - Muscularis","Esophagus - Gastroesophageal Junction", "Colon - Transverse", "Stomach", "Colon - Sigmoid", "Esophagus - Mucosa"),
  'Secretory' = c("Pancreas","Adrenal Gland", "Pituitary", "Thyroid", "Adipose - Subcutaneous","Adipose - Visceral (Omentum)", "Minor Salivary Gland"),
  'Other' = c("Skin - Not Sun Exposed (Suprapubic)", "Skin - Sun Exposed (Lower leg)", "Whole Blood",
              "Lung", "Artery - Tibial", "Artery - Coronary", "Artery - Aorta", "Breast - Mammary Tissue", "Nerve - Tibial")
)

# Create a mapping of column names to group names
column_names <- colnames(mito_nDNA_percent)
groups <- rep(NA, length(column_names))  

for (group in names(group_assignments)) {
  columns_in_group <- group_assignments[[group]]
  groups[column_names %in% columns_in_group] <- group
}


default_group <- "Ungrouped"
groups[is.na(groups)] <- default_group


vertex_names <- V(g)$name
groups <- groups[match(vertex_names, column_names)]


group_colors <- c(
  "CNS" = "#20A600", 
  "Anabolic" = "#EA4497",
  "Contractile" = "#F3AC03", 
  "Reproductive" = "#E6020E",
  "Digestive" = "#204679",
  "Other" = "#EDE279", 
  "Secretory" = "#969696",
  "Ungrouped" = "#BBBBBB"
)

group_colors_transparent <- sapply(group_colors, function(color) adjustcolor(color, alpha.f = 0.8))

# Assign colors to nodes
V(g)$color <- group_colors_transparent[groups]

# Label nodes
V(g)$label <- vertex_names
V(g)$label.cex <- 0.3  
V(g)$label.font <- 2   

# edge width proportional to the strength of correlation
E(g)$width <- E(g)$weight * 3

# node sizes proportional to degree
V(g)$size <- igraph::degree(g) + 4

set.seed(123)
layout <- layout_with_fr(g, weights = E(g)$weight, niter = 500)

plot(g, layout = layout, vertex.frame.color = "white", edge.color = "darkgrey", 
     vertex.label.color = "black", edge.width = E(g)$width, 
     #vertex.label = NA
)