library(tidyverse)
library(ComplexHeatmap)
library(naniar)
library(circlize)
library(ggfortify)


# mtDNA% vs mito-nDNA% correlation within each tissue

correlation_with_pvalues <- sapply(names(mtDNA_percent), function(col_name) {
  cor_result <- cor.test(mtDNA_percent[[col_name]], mito_nDNA_percent[[col_name]])
  list(column = col_name, correlation = cor(mtDNA_percent[[col_name]], mito_nDNA_percent[[col_name]], method = "spearman", use = "pairwise.complete.obs"), p_value = cor_result$p.value)
}, simplify = FALSE)

cor_mtDNA_nDNA <- as.tibble(do.call(rbind, correlation_with_pvalues)) %>%
  column_to_rownames("column") 

cor_mtDNA_nDNA <- cor_mtDNA_nDNA %>%
  mutate(across(where(is.list), ~ map_dbl(., as.numeric))) %>%
  arrange(desc(correlation))

cor_mtDNA_nDNA <- cor_mtDNA_nDNA %>%
  rownames_to_column("Tissue")



#proliferation index
exprs_proliferative <- exprs %>%
  as.data.frame() %>%
  filter(Gene %in% c("MKI67", "RRM2","TOP2A"))


exprs_proliferative <- exprs_proliferative %>%
  as.data.frame() %>%
  column_to_rownames("Gene") %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("X")


exprs_proliferative <- exprs_proliferative %>%
  pivot_longer(colnames(exprs_proliferative[2:ncol(exprs_proliferative)])) %>%
  group_by(X) %>%
  mutate(summed_exprs = sum(value)) %>%
  select(X, summed_exprs) %>%
  ungroup() %>%
  unique()


exprs_proliferative <- exprs_proliferative %>%
  full_join(annotations_sample_attributes, by = "X") %>%
  select(summed_exprs, SMTSD, SUBJID, SMRIN) %>%
  filter(SMRIN>=5.5) %>%
  dplyr::select(-SMRIN) %>%
  unique() %>%
  pivot_wider(names_from = SMTSD, values_from = summed_exprs) %>%
  column_to_rownames("SUBJID")

#ensure dataframes in same order
exprs_proliferative <- exprs_proliferative[rownames(exprs_total_nDNA), colnames(exprs_total_nDNA)]

normalised_proliferative <- (exprs_proliferative/exprs_total_nDNA)*100



proliferative_ranked <- normalised_proliferative %>%
  colMeans(na.rm = T) %>%
  as.matrix() %>%
  as.data.frame() %>%
  arrange(desc(V1)) %>%
  rename("proliferation_score" = V1) %>%
  rownames_to_column("Tissue")

proliferative_ranked <- proliferative_ranked %>%
  filter(! proliferative_ranked$Tissue %in% c("Cervix - Ectocervix", "Cervix - Endocervix", "Fallopian Tube", "Cells - Cultured fibroblasts",
                                              "Bladder", "Cells - EBV-transformed lymphocytes",
                                              "Cells - EBV-transformed lymphocytes", "Kidney - Medulla", "Small Intestine - Terminal Ileum", "Spleen"
  ))


proliferative_ranked <- proliferative_ranked%>%
  mutate(category = case_when(
    grepl("Brain", proliferative_ranked$Tissue) ~ "CNS",
    grepl("Muscle", Tissue) ~ "Contractile",
    grepl("Heart", Tissue) ~ "Contractile",
    grepl("Liver", Tissue) ~ "Anabolic",
    grepl("Kidney", Tissue) ~ "Anabolic",
    grepl("Testis", Tissue) ~ "Reproductive",
    grepl("Vagina", Tissue) ~ "Reproductive",
    grepl("Ovary", Tissue) ~ "Reproductive",
    grepl("Uterus", Tissue) ~ "Reproductive",
    grepl("Prostate", Tissue) ~ "Reproductive",
    grepl("Colon", Tissue) ~ "Digestive",
    grepl("Esophagus", Tissue) ~ "Digestive",
    grepl("Stomach", Tissue) ~ "Digestive",
    grepl("Adipose", Tissue) ~ "Secretory",
    grepl("Salivary", Tissue) ~ "Secretory",
    grepl("Pituitary", Tissue) ~ "Secretory",
    grepl("Thyroid", Tissue) ~ "Secretory",
    grepl("Adrenal", Tissue) ~ "Secretory",
    TRUE ~ "Other"
  ))

category_colors <- c("CNS" = "#20A600", "Anabolic" = "#EA4497", "Contractile" = "#F3AC03", "Reproductive" = "#E6020E","Digestive" = "#204679",
                     "Other" = "#EDE279", "Secretory" = "#969696")





joined_df <- cor_mtDNA_nDNA %>%
  full_join(proliferative_ranked, by = "Tissue")



p <- ggplot(joined_df) + aes(x = log10(proliferation_score), y = correlation, color = category, label = Tissue) + geom_point(size = 3, alpha = 0.7) + 
  labs(x = "log10 Proliferation score", y = "mtDNA% and mito nDNA% (spearman r)") + theme_bw()+ geom_smooth(method = "lm", se = T, color = "#999999", fill = "#999999", alpha = 0.15) + 
  scale_color_manual(values = category_colors) + labs(color = "Tissues") +theme(
    axis.title = element_text(size = 12), 
    axis.text = element_text(size =12),
    plot.margin = margin(t=10, r= 10, b=20, l=10),
    legend.position = "none",
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()
  )
plotly::ggplotly(p)
cor(joined_df$correlation, joined_df$proliferation_score, use = "pairwise.complete.obs", method = "spearman")
cor.test(joined_df$correlation, joined_df$proliferation_score, use = "pairwise.complete.obs", method = "spearman")
