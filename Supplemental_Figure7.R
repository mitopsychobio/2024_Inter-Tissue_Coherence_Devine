library(tidyverse)
library(ComplexHeatmap)
library(naniar)
library(psych)

proteomics_data <- readxl::read_xlsx("NIHMS1624446-supplement-2.xlsx", sheet = 5)

# sum of all mitochondrial proteins in each sample
mito_proteins <- proteomics_data[1:nrow(proteomics_data), 2:ncol(proteomics_data)] %>%
  t() %>%
  as.data.frame() %>%
  filter(!`V3` == "reference")

colnames(mito_proteins) <- mito_proteins[1,]

colnames(mito_proteins)[1:3] <- c("Run", "Tissue", "Sample_id")

mito_proteins <- mito_proteins[2:nrow(mito_proteins),] %>%
  pivot_longer(cols = colnames(mito_proteins[4:ncol(mito_proteins)]))

# filter for mito genes using mito_ensembls dataframe (from Figure 3 R file) 
mito_ensembls <- unique(match_by_ensembl$Ensembl_gtex_iso)
cleaned_ensembls <- sub("\\..*", "", mito_ensembls)

mito_proteins <- mito_proteins %>%
  filter(name %in% cleaned_ensembls)

mito_proteins$value <- as.numeric(mito_proteins$value)

mito_proteins <- mito_proteins %>%
  group_by(Run, Tissue, Sample_id) %>%
  mutate(summed_exprs = sum(value, na.rm = T)) %>%
  ungroup() %>%
  dplyr::select("Run", "Tissue", "Sample_id", "summed_exprs") %>%
  unique() %>%
  dplyr::select(-Run) %>%
  group_by(Tissue, Sample_id) %>%
  mutate(summed_exprs_updated = mean(summed_exprs, na.rm = T)) %>%
  ungroup() %>%
  dplyr::select(-summed_exprs) %>%
  unique()

mito_proteins$Sample_id <- substring(mito_proteins$Sample_id,1,10)

mito_proteins <- mito_proteins %>%
  pivot_wider(names_from = Tissue, values_from = summed_exprs_updated) %>%
  column_to_rownames("Sample_id") 

# filter for tissues of interest
mito_proteins <- mito_proteins %>%
  select(`Muscle - Skeletal`, `Esophagus - Gastroesophageal Junction`, `Artery - Aorta`,
         Thyroid, `Heart - Atrial Appendage`, Stomach, `Colon - Sigmoid`
  )


# Sample ids for joining
sample_ids_df <- proteomics_data[1:nrow(proteomics_data), 2:ncol(proteomics_data)] %>%
  t() %>%
  as.data.frame() %>%
  filter(!`V3` == "reference")

colnames(sample_ids_df) <- sample_ids_df[1,]

colnames(sample_ids_df)[1:3] <- c("Run", "Tissue", "Sample_id")



# sum of all proteins in each sample
temp <- proteomics_data[1:nrow(proteomics_data), 2:ncol(proteomics_data)] %>%
  as.data.frame()

colnames(temp) <- temp[3,]

temp <- temp[4:nrow(temp), 2:ncol(temp)]

temp <- as.data.frame(lapply(temp, as.numeric))
sums_df <- temp %>%
  map_dbl(~sum(.,na.rm=T)) %>%
  as.data.frame()

sums_df <- sums_df %>%
  rownames_to_column("Tissue") %>%
  rename("summed_exprs" = ".")

all_proteins <- sums_df %>%
  filter(!case_when(
    str_detect(Tissue, "reference") ~ TRUE,
    TRUE ~ FALSE
  ))

rm(temp, sums_df)


all_proteins <- all_proteins %>%
  mutate(Tissue = case_when(
    str_ends(Tissue, ".1") ~ str_sub(Tissue, 1, -3),
    str_ends(Tissue, ".2") ~ str_sub(Tissue, 1, -3),
    TRUE ~ Tissue)) %>%
  group_by(Tissue) %>%
  mutate(exprs = mean(summed_exprs)) %>%
  select(-summed_exprs) %>%
  unique() %>%
  rename("Sample_id" = Tissue) %>%
  mutate(Sample_id = gsub("\\.", "-", Sample_id)) %>%
  full_join(sample_ids_df, by = "Sample_id") %>%
  select(Sample_id, Tissue, exprs) %>%
  unique() %>%
  na.omit()

all_proteins$Sample_id <- substring(all_proteins$Sample_id,1,10)

all_proteins <- all_proteins %>%
  pivot_wider(names_from = Tissue, values_from = exprs) %>%
  column_to_rownames("Sample_id") %>%
  select(`Muscle - Skeletal`, `Esophagus - Gastroesophageal Junction`, `Artery - Aorta`,
         Thyroid, `Heart - Atrial Appendage`, Stomach, `Colon - Sigmoid`)




### percentage of each sample composed of mitochondrial proteins
mito_percent <- (mito_proteins/all_proteins)*100

h <- cor(mito_percent, method = "spearman", use = "pairwise.complete.obs")

n_proteomics_cor <- pairwiseCount(mito_proteins)

n_proteomics_cor <- n_proteomics_cor %>%
  as.data.frame() %>%
  replace_with_na_if(.predicate = is.double, condition = ~.x <8) %>%
  `rownames<-`(colnames(n_proteomics_cor))


proteomics_cor_cutoff <- h

proteomics_cor_cutoff[is.na(n_proteomics_cor)] <- NA

mask_value <- min(proteomics_cor_cutoff, na.rm = T) - 1
proteomics_cor_cutoff[upper.tri(proteomics_cor_cutoff, diag = F)] <- mask_value

color_fun <- colorRamp2(c(mask_value, -1, 
                          0, 1), c("white", "blue", "white", "red")
)



Heatmap(proteomics_cor_cutoff, name = "Spearman r", show_column_dend = F, cluster_rows = F, cluster_columns = F,
        col = color_fun, row_names_side = "left", row_names_gp = gpar(fontsize = 14), column_names_gp = gpar(fontsize = 14),
        column_names_rot = 45, show_heatmap_legend = F, row_names_max_width = unit(10, "cm"), 
        column_names_max_height = unit(10, "cm")
)


proteomics_cor_cutoff_vector <- proteomics_cor_cutoff[lower.tri(proteomics_cor_cutoff)]
proteomics_cor_cutoff_vector <- proteomics_cor_cutoff_vector[!is.na(proteomics_cor_cutoff_vector)]

p <- ggplot() + geom_histogram(aes(proteomics_cor_cutoff_vector), binwidth = 0.05) + labs(x = "Spearman r", y = "Count") + xlim(-1, 1) + theme_bw() +theme(
  axis.title = element_text(size = 16), 
  axis.text = element_text(size =16),
  plot.margin = margin(t=10, r= 10, b=20, l=10),
  panel.grid.minor = element_blank(),
  panel.grid.major = element_blank()
) 
plotly::ggplotly(p)

median(proteomics_cor_cutoff_vector, na.rm = T)

