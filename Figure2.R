library(tidyverse)
library(ComplexHeatmap)
library(naniar)


cohort2_data <- readxl::read_xlsx("Supplemental_File_1.xlsx", sheet = 2) %>%
  column_to_rownames("Animal ID:")

cor_matrix <- cor(cohort2_data, use = "pairwise.complete.obs", method = "spearman")


Heatmap(cor_matrix, name = "Spearman r", row_names_gp = gpar(fontsize = 6), column_names_gp = gpar(fontsize = 6), cluster_columns = F, cluster_rows = F,
        show_row_names = F, show_column_names = F)




#CI
CI_cor <- cohort2_data %>%
  select(matches("CI_"))

CI_cor_matrix <- cor(CI_cor, use = "pairwise.complete.obs", method = "spearman")

CI_brain <- CI_cor_matrix[1:17, 1:17]
CI_brain_vector <- CI_brain[lower.tri(CI_brain)]
CI_brain_vector <- c(CI_brain_vector)
CI_brain_vector <- CI_brain_vector[!is.na(CI_brain_vector)]

CI_body <- CI_cor_matrix[18:22,18:22]
CI_body_vector <- CI_body[lower.tri(CI_body)]
CI_body_vector <- c(CI_body_vector)
CI_body_vector <- CI_body_vector[!is.na(CI_body_vector)]

CI_brain_body <- CI_cor_matrix[18:22,1:17]
CI_brain_body_vector <- c(CI_brain_body)
CI_brain_body_vector <- CI_brain_body_vector[!is.na(CI_brain_body_vector)]

#CII
CII_cor <- cohort2_data %>%
  select(matches("CII_"))

CII_cor_matrix <- cor(CII_cor, use = "pairwise.complete.obs", method = "spearman")

CII_brain <- CII_cor_matrix[1:17, 1:17]
CII_brain_vector <- CII_brain[lower.tri(CII_brain)]
CII_brain_vector <- c(CII_brain_vector)
CII_brain_vector <- CII_brain_vector[!is.na(CII_brain_vector)]

CII_body <- CII_cor_matrix[18:22,18:22]
CII_body_vector <- CII_body[lower.tri(CII_body)]
CII_body_vector <- c(CII_body_vector)
CII_body_vector <- CII_body_vector[!is.na(CII_body_vector)]

CII_brain_body <- CII_cor_matrix[18:22,1:17]
CII_brain_body_vector <- c(CII_brain_body)
CII_brain_body_vector <- CII_brain_body_vector[!is.na(CII_brain_body_vector)]


#CIV
CIV_cor <- cohort2_data %>%
  select(matches("CIV_"))

CIV_cor_matrix <- cor(CIV_cor, use = "pairwise.complete.obs", method = "spearman")

CIV_brain <- CIV_cor_matrix[1:17, 1:17]
CIV_brain_vector <- CIV_brain[lower.tri(CIV_brain)]
CIV_brain_vector <- c(CIV_brain_vector)
CIV_brain_vector <- CIV_brain_vector[!is.na(CIV_brain_vector)]

CIV_body <- CIV_cor_matrix[18:22,18:22]
CIV_body_vector <- CIV_body[lower.tri(CIV_body)]
CIV_body_vector <- c(CIV_body_vector)
CIV_body_vector <- CIV_body_vector[!is.na(CIV_body_vector)]

CIV_brain_body <- CIV_cor_matrix[18:22,1:17]
CIV_brain_body_vector <- c(CIV_brain_body)
CIV_brain_body_vector <- CIV_brain_body_vector[!is.na(CIV_brain_body_vector)]

#CS
CS_cor <- cohort2_data %>%
  select(matches("CS_"))

CS_cor_matrix <- cor(CS_cor, use = "pairwise.complete.obs", method = "spearman")

CS_brain <- CS_cor_matrix[1:17, 1:17]
CS_brain_vector <- CS_brain[lower.tri(CS_brain)]
CS_brain_vector <- c(CS_brain_vector)
CS_brain_vector <- CS_brain_vector[!is.na(CS_brain_vector)]

CS_body <- CS_cor_matrix[18:22,18:22]
CS_body_vector <- CS_body[lower.tri(CS_body)]
CS_body_vector <- c(CS_body_vector)
CS_body_vector <- CS_body_vector[!is.na(CS_body_vector)]

CS_brain_body <- CS_cor_matrix[18:22,1:17]
CS_brain_body_vector <- c(CS_brain_body)
CS_brain_body_vector <- CS_brain_body_vector[!is.na(CS_brain_body_vector)]

#mtDNA
mtDNA_cor <- cohort2_data %>%
  select(matches("mtDNA_"))

mtDNA_cor_matrix <- cor(mtDNA_cor, use = "pairwise.complete.obs", method = "spearman")

mtDNA_brain <- mtDNA_cor_matrix[1:17, 1:17]
mtDNA_brain_vector <- mtDNA_brain[lower.tri(mtDNA_brain)]
mtDNA_brain_vector <- c(mtDNA_brain_vector)
mtDNA_brain_vector <- mtDNA_brain_vector[!is.na(mtDNA_brain_vector)]

mtDNA_body <- mtDNA_cor_matrix[18:22,18:22]
mtDNA_body_vector <- mtDNA_body[lower.tri(mtDNA_body)]
mtDNA_body_vector <- c(mtDNA_body_vector)
mtDNA_body_vector <- mtDNA_body_vector[!is.na(mtDNA_body_vector)]

mtDNA_brain_body <- mtDNA_cor_matrix[18:22,1:17]
mtDNA_brain_body_vector <- c(mtDNA_brain_body)
mtDNA_brain_body_vector <- mtDNA_brain_body_vector[!is.na(mtDNA_brain_body_vector)]


# frequency distributions
freq_dist_brain <- append(CI_brain_vector, c(CII_brain_vector, CIV_brain_vector, CS_brain_vector, mtDNA_brain_vector))
freq_dist_body <- append(CI_body_vector, c(CII_body_vector, CIV_body_vector, CS_body_vector, mtDNA_body_vector))
freq_dist_brain_body <- append(CI_brain_body_vector, c(CII_brain_body_vector, CIV_brain_body_vector, CS_brain_body_vector, mtDNA_brain_body_vector))


p <- ggplot()+ geom_histogram(aes(freq_dist_brain), fill = "#9f3864", binwidth = 0.05, alpha = 0.3, color = "black", size = 0.2)  +
  geom_histogram(aes(freq_dist_brain_body), binwidth = 0.05, fill = "#e0b61d", alpha = 0.45, color = "black", size = 0.2)+ 
  geom_histogram(aes(freq_dist_body), binwidth = 0.05, alpha = 0.4, fill="#0D6153", color = "black", size = 0.2)  + 
  labs(x = "Spearman r", y = "Count") + 
  xlim(-0.8, 0.8) + theme_bw() +
  theme(
    axis.title = element_text(size = 12), 
    axis.text = element_text(size =10),
    plot.margin = margin(t=10, r= 10, b=20, l=10),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) 
plotly::ggplotly(p)


median(freq_dist_brain, na.rm = T)
median(freq_dist_body, na.rm = T)
median(freq_dist_brain_body, na.rm = T)

wilcox.test(freq_dist_brain, mu = 0)
wilcox.test(freq_dist_body, mu = 0)
wilcox.test(freq_dist_brain_body, mu = 0)
