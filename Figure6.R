library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(ggfortify)


# PGC1-alpha

# filter exprs dataframe (from Figure 3 R script) for PGC1alpha
exprs_pgc1a <- exprs %>%
  as.data.frame() %>%
  filter(Gene %in% c("PPARGC1A"))


exprs_pgc1a <- exprs_pgc1a %>%
  column_to_rownames("Gene") %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("X")


exprs_pgc1a <- exprs_pgc1a %>%
  pivot_longer(colnames(exprs_pgc1a[2:ncol(exprs_pgc1a)])) %>%
  group_by(X) %>%
  mutate(summed_exprs = sum(value)) %>%
  select(X, summed_exprs) %>%
  ungroup() %>%
  unique()

# join with sample annotations (from Figure 3 R script)
exprs_pgc1a <- exprs_pgc1a %>%
  full_join(annotations_sample_attributes, by = "X") %>%
  select(summed_exprs, SMTSD, SUBJID, SMRIN) %>%
  filter(SMRIN>=5.5) %>%
  select(-SMRIN) %>%
  unique() %>%
  pivot_wider(names_from = SMTSD, values_from = summed_exprs) %>%
  column_to_rownames("SUBJID")

# express pgc1alpha transcript abundance as percentage of all nuclear transcripts
pgc1a_percentage <- (exprs_pgc1a/exprs_total_nDNA)*100
rm(exprs_pgc1a)

# select 45 tissues of interest
tissue_cor_pgc1a <- pgc1a_percentage %>%
  select( - `Cervix - Ectocervix`, - `Cervix - Endocervix`, - `Fallopian Tube`, -`Cells - Cultured fibroblasts`,
          - Bladder, -`Cells - EBV-transformed lymphocytes`, -`Small Intestine - Terminal Ileum`, - Spleen,
          -`Cells - EBV-transformed lymphocytes`, - `Kidney - Medulla`)

rm(pgc1a_percentage)



# correlation of pgc1a with mito nDNA%
correlation_with_pvalues <- sapply(names(tissue_cor_pgc1a), function(col_name) {
  cor_result <- cor.test(tissue_cor_pgc1a[[col_name]], mito_nDNA_percent[[col_name]])
  list(column = col_name, correlation = cor(tissue_cor_pgc1a[[col_name]], mito_nDNA_percent[[col_name]], method = "spearman", use = "pairwise.complete.obs"), p_value = cor_result$p.value)
}, simplify = FALSE)


cor_df <- as.tibble(do.call(rbind, correlation_with_pvalues)) %>%
  column_to_rownames("column") 

cor_df_pgc1a <- cor_df %>%
  mutate(across(where(is.list), ~ map_dbl(., as.numeric))) %>%
  arrange(desc(correlation)) 

Heatmap(cor_df_pgc1a[1], name =  "Spearman r", cluster_rows = F, cluster_columns = F, cell_fun = function(j, i, x, y, width, height, fill) {grid.text(sprintf("%.2f", cor_df_pgc1a[i, j]), x, y, gp = gpar(fontsize = 7))},
        width = ncol(cor_df)*unit(8, "mm"), row_names_gp = gpar(fontsize = 6), col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")), show_column_names = F, row_names_side = "left")


# ISR

# filter exprs dataframe (from Figure 3 R script) for ISR genes
exprs_isr <- exprs %>%
  as.data.frame() %>%
  filter(Gene %in% c("ATF4", "ATF5", "DDIT3", "GDF15"))


exprs_isr <- exprs_isr %>%
  column_to_rownames("Gene") %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("X")

# sum ISR transcripts in each sample
exprs_isr <- exprs_isr %>%
  pivot_longer(colnames(exprs_isr[2:ncol(exprs_isr)])) %>%
  group_by(X) %>%
  mutate(avg_exprs = sum(value)) %>%
  select(X, avg_exprs) %>%
  ungroup() %>%
  unique()

# join with sample annotations
exprs_isr <- exprs_isr %>%
  full_join(annotations_sample_attributes, by = "X") %>%
  select(avg_exprs, SMTSD, SUBJID, SMRIN) %>%
  filter(SMRIN>=5.5) %>%
  select(-SMRIN) %>%
  unique() %>%
  pivot_wider(names_from = SMTSD, values_from = avg_exprs) %>%
  column_to_rownames("SUBJID")


# express ISR transcript abundance as a percentage of all nuclear transcripts
isr_percentage <- exprs_isr/exprs_total_nDNA


# filter for tissues of interest
tissue_cor_isr <- isr_percentage %>%
  select( - `Cervix - Ectocervix`, - `Cervix - Endocervix`, - `Fallopian Tube`, -`Cells - Cultured fibroblasts`,
          - Bladder, -`Cells - EBV-transformed lymphocytes`, -`Small Intestine - Terminal Ileum`, - Spleen,
          -`Cells - EBV-transformed lymphocytes`, - `Kidney - Medulla`)


# correlation of ISR with mito nDNA%

correlation_with_pvalues_isr <- sapply(names(tissue_cor_isr), function(col_name) {
  cor_result <- cor.test(tissue_cor_isr[[col_name]], mito_nDNA_percent[[col_name]])
  list(column = col_name, correlation = cor(tissue_cor_isr[[col_name]], mito_nDNA_percent[[col_name]], method = "spearman", use = "pairwise.complete.obs"), p_value = cor_result$p.value)
}, simplify = FALSE)

cor_df <- as.tibble(do.call(rbind, correlation_with_pvalues_isr)) %>%
  column_to_rownames("column") 

cor_df <- cor_df %>%
  mutate(across(where(is.list), ~ map_dbl(., as.numeric))) %>%
  arrange(desc(correlation))

Heatmap(cor_df[1], name =  "Spearman r", cluster_rows = F, cluster_columns = F, cell_fun = function(j, i, x, y, width, height, fill) {grid.text(sprintf("%.2f", cor_df[i, j]), x, y, gp = gpar(fontsize = 7))},
        width = ncol(cor_df)*unit(9, "mm"), row_names_gp = gpar(fontsize = 6), col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")), show_column_names = F, row_names_side = "left")
