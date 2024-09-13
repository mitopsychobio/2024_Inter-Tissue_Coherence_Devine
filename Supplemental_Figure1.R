library(tidyverse)
library(ComplexHeatmap)
library(circlize)


cohort1_data <- readxl::read_xlsx("Supplemental_File_1.xlsx", sheet = 1) %>%
  column_to_rownames("Animal ID")
cor_matrix <- cor(cohort1_data, use = "pairwise.complete.obs", method = "spearman")


mask_value <- min(cor_matrix, na.rm = T) - 1
cor_matrix[upper.tri(cor_matrix, diag = F)] <- mask_value

color_fun <- colorRamp2(c(mask_value, -1, 
                          0, 1), c("white", "blue", "white", "red")
)


Heatmap(cor_matrix, name = "Spearman r", row_names_gp = gpar(fontsize = 9), column_names_gp = gpar(fontsize = 9), cluster_columns = F,
        cluster_rows = F, show_column_names = F, na_col = "grey", row_names_side = "left", column_names_rot = 45, col = color_fun, show_row_names = F
)



#CI
CI_cor <- cohort1_data %>%
  select(matches("CI_"))

CI_cor_matrix <- cor(CI_cor, use = "pairwise.complete.obs", method = "spearman")
CI_cor_vector <- CI_cor_matrix[lower.tri(CI_cor_matrix)]
CI_cor_vector <- c(CI_cor_vector)
CI_cor_vector <- CI_cor_vector[!is.na(CI_cor_vector)]


#CII
CII_cor <- cohort1_data %>%
  select(matches("_CII"))

CII_cor_matrix <- cor(CII_cor, use = "pairwise.complete.obs", method = "spearman")
CII_cor_vector <- CII_cor_matrix[lower.tri(CII_cor_matrix)]
CII_cor_vector <- c(CII_cor_vector)
CII_cor_vector <- CII_cor_vector[!is.na(CII_cor_vector)]


#CIV
CIV_cor <- cohort1_data %>%
  select(matches("_CIV"))

CIV_cor_matrix <- cor(CIV_cor, use = "pairwise.complete.obs", method = "spearman")
CIV_cor_vector <- CIV_cor_matrix[lower.tri(CIV_cor_matrix)]
CIV_cor_vector <- c(CIV_cor_vector)
CIV_cor_vector <- CIV_cor_vector[!is.na(CIV_cor_vector)]



#CS
CS_cor <- cohort1_data %>%
  select(matches("_CS"))

CS_cor_matrix <- cor(CS_cor, use = "pairwise.complete.obs", method = "spearman")
CS_cor_vector <- CS_cor_matrix[lower.tri(CS_cor_matrix)]
CS_cor_vector <- c(CS_cor_vector)
CS_cor_vector <- CS_cor_vector[!is.na(CS_cor_vector)]



#mtDNA
mtDNA_cor <- cohort1_data %>%
  select(matches("_mtDNA"))

mtDNA_cor_matrix <- cor(mtDNA_cor, use = "pairwise.complete.obs", method = "spearman")
mtDNA_cor_vector <- mtDNA_cor_matrix[lower.tri(mtDNA_cor_matrix)]
mtDNA_cor_vector <- c(mtDNA_cor_vector)
mtDNA_cor_vector <- mtDNA_cor_vector[!is.na(mtDNA_cor_vector)]



# frequency distribution of correlation coefficients
freq_dist_vector <- append(mtDNA_cor_vector, c(CS_cor_vector,CI_cor_vector, CII_cor_vector, CIV_cor_vector))


p <- ggplot() + geom_histogram(aes(freq_dist_vector), binwidth = 0.05) + labs(x = "Spearman r", y = "Count") + xlim(-1, 1) + theme_bw() +theme(
  axis.title = element_text(size = 12), 
  axis.text = element_text(size =12),
  plot.margin = margin(t=10, r= 10, b=20, l=10),
  panel.grid.minor = element_blank(),
  panel.grid.major = element_blank()
) 
plotly::ggplotly(p)


median(freq_dist_vector, na.rm = TRUE)
wilcox.test(freq_dist_vector, mu = 0)





p <- ggplot() + aes(x = cohort1_data$Hippocampus_CS, y = cohort1_data$Liver_CS) + geom_point(alpha = 0.15, size = 3) + 
  labs(x = "Hippocampus CS Activity", y = "Liver CS Activity") + theme_bw() + geom_smooth(method = "lm", se = T, color = "#999999", fill = "#999999", alpha = 0.15) +
  theme(
    axis.title = element_text(size = 12), 
    axis.text = element_text(size =12),
    plot.margin = margin(t=10, r= 10, b=20, l=10),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()
  ) 
plotly::ggplotly(p)
cor(cohort1_data$Liver_CS, cohort1_data$Hippocampus_CS, use = "pairwise.complete.obs", method = "spearman")
cor.test(cohort1_data$Hippocampus_CS, cohort1_data$Liver_CS, method = "spearman")


p <- ggplot() + aes(x = cohort1_data$BrownFat_CIV, y = cohort1_data$Liver_CIV) + geom_point(alpha = 0.15, size = 3) + 
  labs(x = "Brown Fat CIV Activity", y = "Liver CIV Activity") + theme_bw() + geom_smooth(method = "lm", se = T, color = "#999999", fill = "#999999", alpha = 0.15) +
  theme(
    axis.title = element_text(size = 12), 
    axis.text = element_text(size =12),
    plot.margin = margin(t=10, r= 10, b=20, l=10),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()
  ) 
plotly::ggplotly(p)
cor(cohort1_data$Liver_CIV, cohort1_data$BrownFat_CIV, use = "pairwise.complete.obs", method = "spearman")
cor.test(cohort1_data$BrownFat_CIV, cohort1_data$Liver_CIV, method = "spearman")



p <- ggplot() + aes(x = cohort1_data$Hippocampus_CS, y = cohort1_data$Muscle_CS) + geom_point(alpha = 0.15, size = 3) + 
  labs(x = "Hippocampus CS Activity", y = "Muscle CS Activity") + theme_bw() + geom_smooth(method = "lm", se = T, color = "#999999", fill = "#999999", alpha = 0.15) +
  theme(
    axis.title = element_text(size = 12), 
    axis.text = element_text(size =12),
    plot.margin = margin(t=10, r= 10, b=20, l=10),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()
  ) 
plotly::ggplotly(p)
cor(cohort1_data$Hippocampus_CS, cohort1_data$Muscle_CS, use = "pairwise.complete.obs", method = "spearman")
cor.test(cohort1_data$Hippocampus_CS, cohort1_data$Muscle_CS, method = "spearman")

