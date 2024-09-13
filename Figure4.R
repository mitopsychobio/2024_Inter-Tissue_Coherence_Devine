library(tidyverse)
library(ComplexHeatmap)
library(naniar)
library(circlize)
library(ggfortify)
library(dendextend)
library(psych)



# run correlation matrix on mito_nDNA% (from Figure 3 R file)
cor_matrix <- cor(mito_nDNA_percent, use = "pairwise.complete.obs", method = "spearman")

# count sample size of each tissue pair
n_tissue_cor <- pairwiseCount(mito_nDNA_percent) %>%
  as.data.frame()

# apply sample size cutoff
n_tissue_cor <- n_tissue_cor %>%
  replace_with_na_if(.predicate = is.double, condition = ~.x <10) %>%
  `rownames<-`(colnames(n_tissue_cor))


tissue_cor_cutoff <- cor_matrix

names_order <- c("Brain - Caudate (basal ganglia)", "Brain - Anterior cingulate cortex (BA24)", "Brain - Putamen (basal ganglia)", "Brain - Nucleus accumbens (basal ganglia)", 
                 "Brain - Amygdala", "Brain - Substantia nigra", "Brain - Hippocampus", "Brain - Hypothalamus","Brain - Cortex", "Brain - Frontal Cortex (BA9)",
                 "Brain - Cerebellum", "Brain - Cerebellar Hemisphere", "Brain - Spinal cord (cervical c-1)", "Pituitary", "Kidney - Cortex","Liver", "Esophagus - Muscularis","Esophagus - Gastroesophageal Junction", 
                 "Thyroid", "Breast - Mammary Tissue", "Colon - Sigmoid",  "Adipose - Subcutaneous", "Nerve - Tibial",
                 "Lung", "Artery - Tibial", "Artery - Coronary", "Artery - Aorta", "Ovary", "Prostate", "Adipose - Visceral (Omentum)", 
                 "Heart - Atrial Appendage", "Muscle - Skeletal", "Colon - Transverse", "Heart - Left Ventricle", "Stomach",
                 "Minor Salivary Gland", "Skin - Not Sun Exposed (Suprapubic)", "Skin - Sun Exposed (Lower leg)", "Whole Blood", "Uterus",
                 "Testis", "Vagina", "Esophagus - Mucosa", "Pancreas","Adrenal Gland")

tissue_cor_cutoff[is.na(n_tissue_cor)] <- NA
tissue_cor_cutoff <- tissue_cor_cutoff[names_order, names_order]


mask_value <- min(cor_matrix, na.rm = T) - 1
tissue_cor_cutoff[upper.tri(tissue_cor_cutoff, diag = F)] <- mask_value

color_fun <- colorRamp2(c(mask_value, -1, 
                          0, 1), c("white", "blue", "white", "red")
)

# heatmap of inter-tissue correlation matrix
Heatmap(tissue_cor_cutoff, name = "Spearman r", row_names_gp = gpar(fontsize = 9), column_names_gp = gpar(fontsize = 9), cluster_columns = F,
        cluster_rows = F, show_column_names = T, na_col = "grey", row_names_side = "left", column_names_rot = 45, col = color_fun
)



# brain-brain distribution
brain <- mito_nDNA_percent %>%
  select(`Brain - Cortex`,
         `Brain - Frontal Cortex (BA9)`, 
         `Brain - Hypothalamus`, 
         `Brain - Hippocampus`, 
         `Brain - Anterior cingulate cortex (BA24)`, 
         `Brain - Amygdala`,
         `Brain - Cerebellum`, `Brain - Cerebellar Hemisphere`, 
         `Brain - Substantia nigra`, 
         `Brain - Caudate (basal ganglia)`, 
         `Brain - Putamen (basal ganglia)`, 
         `Brain - Nucleus accumbens (basal ganglia)`
  ) 

cor_brain <- cor(brain, use = "pairwise.complete.obs", method = "spearman")

cor_brain <- cor_brain[lower.tri(cor_brain)]
cor_brain_vector <- c(cor_brain)

median(cor_brain_vector, na.rm = TRUE)



# body-body distribution
body <- mito_nDNA_percent %>%
  select(-`Brain - Amygdala`, -`Brain - Anterior cingulate cortex (BA24)`, -`Brain - Caudate (basal ganglia)`, -`Brain - Cerebellar Hemisphere`,
         -`Brain - Cerebellum`, -`Brain - Cortex`, -`Brain - Frontal Cortex (BA9)`, -`Brain - Hippocampus`, -`Brain - Hypothalamus`,
         -`Brain - Nucleus accumbens (basal ganglia)`, -`Brain - Putamen (basal ganglia)`, -`Brain - Substantia nigra`)

cor_body <- cor(body, use = "pairwise.complete.obs", method = "spearman")

n_cor_body <- pairwiseCount(body)
n_cor_body <- n_cor_body %>%
  as.data.frame()

n_cor_body <- n_cor_body %>%
  replace_with_na_if(.predicate = is.double, condition = ~.x <10) %>%
  `rownames<-`(colnames(n_cor_body))


cor_body_cutoff <- cor_body

cor_body_cutoff[is.na(n_cor_body)] <- NA
cor_body_cutoff <- cor_body_cutoff[lower.tri(cor_body_cutoff)]
cor_body_cutoff_vector <- c(cor_body_cutoff)
cor_body_cutoff_vector <- cor_body_cutoff_vector[!is.na(cor_body_cutoff_vector)]
median(cor_body_cutoff_vector, na.rm = TRUE)


# brain-body distribtion
tissue_cor_cutoff <- cor_matrix
tissue_cor_cutoff[is.na(n_tissue_cor)] <- NA
tissue_cor_cutoff <- tissue_cor_cutoff[names_order, names_order]

cor_brain_body <- tissue_cor_cutoff[13:nrow(tissue_cor_cutoff),1:12]
cor_brain_body_vector <- c(cor_brain_body)
median(cor_brain_body_vector, na.rm = TRUE)




# frequency distribution of brain-brain, body-body, brain-body
p <- ggplot() + geom_histogram(aes(cor_body_cutoff_vector), binwidth = 0.05, alpha = 0.4, fill="#0D6153", color = "black", size = 0.2) +
  geom_histogram(aes(cor_brain_body_vector), binwidth = 0.05, fill = "#e0b61d", alpha = 0.4, color = "black", size = 0.2) + 
  geom_histogram(aes(cor_brain_vector), fill = "#9f3864", binwidth = 0.05, alpha = 0.4, color = "black", size = 0.2) + 
  labs(x = "Spearman r", y = "Count") + 
  xlim(-0.8, 0.8) + theme_bw() +
  theme(
    axis.title = element_text(size = 12), 
    axis.text = element_text(size =10),
    plot.margin = margin(t=10, r= 10, b=20, l=10)
  ) 
plotly::ggplotly(p)



wilcox.test(cor_body_cutoff_vector, mu = 0)
wilcox.test(cor_brain_body_vector, mu = 0)
wilcox.test(cor_brain_vector, mu = 0)


### biplots

p <- ggplot() + aes(x = mito_nDNA_percent$`Adrenal Gland`, y = mito_nDNA_percent$Liver) + geom_point(alpha = 0.3, size = 3, color = "#0D6153") + 
  labs(x = "Adrenal Gland", y = "Liver") + theme_bw() + geom_smooth(method = "lm", se = T, color = "#0D6153", fill = "#0D6153", alpha = 0.15) +theme(
    axis.title = element_text(size = 12), 
    axis.text = element_text(size =12),
    plot.margin = margin(t=10, r= 10, b=20, l=10),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()
  ) 
plotly::ggplotly(p)
cor(mito_nDNA_percent$`Adrenal Gland`, mito_nDNA_percent$Liver, use = "pairwise.complete.obs", method = "spearman")
cor.test(mito_nDNA_percent$`Adrenal Gland`, mito_nDNA_percent$Liver, method = "spearman")



p <- ggplot() + aes(x = mito_nDNA_percent$`Brain - Nucleus accumbens (basal ganglia)`, y = mito_nDNA_percent$`Brain - Substantia nigra`) + geom_point(alpha = 0.3, size = 3, color = "#9f3864") +
  labs(x = "NAc", y = "SN") + theme_bw() + geom_smooth(method = "lm", se = T, color = "#9f3864", fill = "#9f3864", alpha = 0.15) +theme(
    axis.title = element_text(size = 12), 
    axis.text = element_text(size =12),
    plot.margin = margin(t=10, r= 10, b=20, l=10),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()
  ) 
plotly::ggplotly(p)
cor(mito_nDNA_percent$`Brain - Substantia nigra`, mito_nDNA_percent$`Brain - Nucleus accumbens (basal ganglia)`, use = "pairwise.complete.obs", method = "spearman")
cor.test(mito_nDNA_percent$`Brain - Substantia nigra`, mito_nDNA_percent$`Brain - Nucleus accumbens (basal ganglia)`, method = "spearman")


p <- ggplot() + aes(x = mito_nDNA_percent$`Adrenal Gland`, y = mito_nDNA_percent$`Brain - Amygdala`) + geom_point(alpha = 0.3, size = 3, color = "#e0b61d") + 
  labs(x = "Adrenal Gland", y = "Amygdala") + geom_smooth(method = "lm", se = T, color = "#e0b61d", fill = "#e0b61d", alpha = 0.15) + theme_bw() +theme(
    axis.title = element_text(size = 12), 
    axis.text = element_text(size =12),
    plot.margin = margin(t=10, r= 10, b=20, l=10),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()
  ) 
plotly::ggplotly(p)
cor(mito_nDNA_percent$`Adrenal Gland`, mito_nDNA_percent$`Brain - Amygdala`, use = "pairwise.complete.obs", method = "spearman")
cor.test(mito_nDNA_percent$`Brain - Amygdala`, mito_nDNA_percent$`Adrenal Gland`, method = "spearman")
