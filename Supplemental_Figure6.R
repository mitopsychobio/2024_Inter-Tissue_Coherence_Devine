library(tidyverse)
library(ComplexHeatmap)
library(naniar)
library(circlize)
library(ggfortify)
library(psych)


mtDNAcn <- readxl::read_xlsx("pnas.2402291121.sd01.xlsx", sheet = 1)

mtDNAcn <- mtDNAcn[,11:ncol(mtDNAcn)] %>%
  mutate(SUBJID = substring(mtDNAcn$`Sample ID for data sharing and public release`,1,10)) %>%
  as.data.frame()

mtDNAcn$`mtCN Per Cell` <- as.numeric(mtDNAcn$`mtCN Per Cell`)

mtDNAcn <- mtDNAcn %>%
  mutate(SUBJID_tmp = case_when(
    substr(SUBJID, nchar(SUBJID), nchar(SUBJID)) == "-" ~ substr(SUBJID, 1, nchar(SUBJID)-1),
    .default = SUBJID
  )) %>%
  dplyr::select(-SUBJID) %>%
  rename(SUBJID = SUBJID_tmp) %>%
  as.data.frame() %>%
  select(SUBJID, `mtCN Per Cell`, `Tissue Site Detail`) %>%
  group_by(SUBJID,`Tissue Site Detail`) %>%
  unique() %>%
  mutate(avg_exprs = mean(`mtCN Per Cell`, na.rm = T)) %>%
  select(SUBJID, avg_exprs, `Tissue Site Detail`) %>%
  unique() %>%
  pivot_wider(names_from = `Tissue Site Detail`, values_from = avg_exprs) %>%
  column_to_rownames("SUBJID")

mtDNAcn <- mtDNAcn %>%
  select(-`Cervix - Ectocervix`, -`Cervix - Endocervix`, -`Cells - EBV-transformed lymphocytes`, - `Fallopian Tube`, -Bladder,
         -Spleen, -`Small Intestine - Terminal Ileum`)

cor_mtDNAcn <- cor(mtDNAcn, use = "pairwise.complete.obs", method = "spearman")


n_tissue_cor <- pairwiseCount(mtDNAcn) %>%
  as.data.frame()

n_tissue_cor <- n_tissue_cor %>%
  replace_with_na_if(.predicate = is.double, condition = ~.x <10) %>%
  `rownames<-`(colnames(n_tissue_cor))

tissue_cor_cutoff <- cor_mtDNAcn

# same order as transcript-based analysis
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

mask_value <- min(cor_mtDNAcn, na.rm = T) - 1
tissue_cor_cutoff[upper.tri(tissue_cor_cutoff, diag = F)] <- mask_value

color_fun <- colorRamp2(c(mask_value, -1, 
                          0, 1), c("white", "blue", "white", "red")
)


Heatmap(tissue_cor_cutoff, name = "Spearman r", row_names_gp = gpar(fontsize = 9), column_names_gp = gpar(fontsize = 9), cluster_columns = F,
        cluster_rows = F, show_column_names = T, na_col = "grey", row_names_side = "left", column_names_rot = 45, col = color_fun
)
tissue_cor_cutoff <- tissue_cor_cutoff[lower.tri(tissue_cor_cutoff)]
cor_vector_cutoff <- c(tissue_cor_cutoff)
cor_vector_cutoff <- cor_vector_cutoff[!is.na(cor_vector_cutoff)]

p <- ggplot() + geom_histogram(aes(cor_vector_cutoff), binwidth = 0.05) + labs(x = "spearman r", y = "number of pairs of tissues") + xlim(-1, 1) + theme_bw() +
  theme(
    axis.title = element_text(size = 18), 
    axis.text = element_text(size =14),
    plot.margin = margin(t=10, r= 10, b=20, l=10)
  ) 
plotly::ggplotly(p)

median(cor_vector_cutoff, na.rm = TRUE)



# brain-brain distribution
brain <- mtDNAcn %>%
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
body <- mtDNAcn %>%
  select(-`Brain - Amygdala`, -`Brain - Anterior cingulate cortex (BA24)`, -`Brain - Caudate (basal ganglia)`, -`Brain - Cerebellar Hemisphere`,
         -`Brain - Cerebellum`, -`Brain - Cortex`, -`Brain - Frontal Cortex (BA9)`, -`Brain - Hippocampus`, -`Brain - Hypothalamus`,
         -`Brain - Nucleus accumbens (basal ganglia)`, -`Brain - Putamen (basal ganglia)`, -`Brain - Substantia nigra`)

cor_body <- cor(body, use = "pairwise.complete.obs", method = "spearman")

n_cor_body <- pairwiseCount(body)

n_cor_body <- n_cor_body %>%
  as.data.frame %>%
  replace_with_na_if(.predicate = is.double, condition = ~.x <10) %>%
  `rownames<-`(colnames(n_cor_body))

cor_body_cutoff <- cor_body
cor_body_cutoff[is.na(n_cor_body)] <- NA
cor_body_cutoff <- cor_body_cutoff[lower.tri(cor_body_cutoff)]
cor_body_cutoff_vector <- c(cor_body_cutoff)
cor_body_cutoff_vector <- cor_body_cutoff_vector[!is.na(cor_body_cutoff_vector)]
median(cor_body_cutoff_vector, na.rm = TRUE)


# brain-body distribution
tissue_cor_cutoff <- cor_mtDNAcn
tissue_cor_cutoff[is.na(n_tissue_cor)] <- NA
tissue_cor_cutoff <- tissue_cor_cutoff[names_order, names_order]

cor_brain_body <- tissue_cor_cutoff[13:nrow(tissue_cor_cutoff),1:12]
cor_brain_body_vector <- c(cor_brain_body)
median(cor_brain_body, na.rm = TRUE)




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


