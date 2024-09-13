library(tidyverse)
library(factoextra)
library(ComplexHeatmap)


exprs_sub <- mito_nDNA_percent %>%
  dplyr::select(`Heart - Atrial Appendage`,`Muscle - Skeletal`, `Brain - Cortex`, `Adipose - Subcutaneous`) %>%
  na.omit() 


ratio_df <- exprs_sub %>%
  mutate(
    "Adipose - Subcutaneous/Muscle - Skeletal" = exprs_sub$`Adipose - Subcutaneous`/exprs_sub$`Muscle - Skeletal`,
    "Adipose - Subcutaneous/Heart - Atrial Appendage" = exprs_sub$`Adipose - Subcutaneous`/exprs_sub$`Heart - Atrial Appendage`,
    "Adipose - Subcutaneous/Brain" = exprs_sub$`Adipose - Subcutaneous`/exprs_sub$`Brain - Cortex`,
    "Muscle - Skeletal/Adipose - Subcutaneous" = exprs_sub$`Muscle - Skeletal`/exprs_sub$`Adipose - Subcutaneous`,
    "Muscle - Skeletal/Heart - Atrial Appendage" = exprs_sub$`Muscle - Skeletal`/exprs_sub$`Heart - Atrial Appendage`,
    "Muscle - Skeletal/Brian - Cortex" = exprs_sub$`Muscle - Skeletal`/exprs_sub$`Brain - Cortex`,
    "Heart - Atrial Appendage/Adipose - Subcutaneous" = exprs_sub$`Heart - Atrial Appendage`/exprs_sub$`Adipose - Subcutaneous`,
    "Heart - Atrial Appendage/Muscle - Skeletal" = exprs_sub$`Heart - Atrial Appendage`/exprs_sub$`Muscle - Skeletal`,
    "Heart - Atrial Appendage/Brain - Cortex" = exprs_sub$`Heart - Atrial Appendage`/exprs_sub$`Brain - Cortex`,
    "Brain - Cortex/Adipose - Subcutaneous" = exprs_sub$`Brain - Cortex`/exprs_sub$`Adipose - Subcutaneous`,
    "Brain - Cortex/Muscle - Skeletal" = exprs_sub$`Brain - Cortex`/exprs_sub$`Muscle - Skeletal`,
    "Brain - Cortex/Heart - Atrial Appendage" = exprs_sub$`Brain - Cortex`/exprs_sub$`Heart - Atrial Appendage`
  ) %>%
  dplyr::select(-`Adipose - Subcutaneous`, -`Muscle - Skeletal`, -`Brain - Cortex`, - `Heart - Atrial Appendage`)




set.seed(123)
kmeans_res <- kmeans(scale(ratio_df), centers = 3, nstart = 50)

clusters <- kmeans_res$cluster

p <- fviz_cluster(kmeans_res, data = ratio_df, geom = "point",
                  ellipse.type = "convex", ggtheme = theme_bw(),
                  palette = c("1" = "#A48AD3",
                              "2" = "#1CC5FE",
                              "3" = "#6FC7CF"), 
                  pointsize = 3,
                  show.clust.cent = FALSE,
                  main = NULL) +theme(
                    plot.margin = margin(t=10, r= 10, b=20, l=10),
                    panel.grid.minor = element_blank(),
                    panel.grid.major = element_blank(),
                    legend.position = "none", 
                    axis.title = element_text(size = 14),      
                    axis.text = element_text(size = 12)) 


plotly::ggplotly(p)





exprs_sub <- mito_nDNA_percent %>%
  dplyr::select(`Heart - Atrial Appendage`,`Muscle - Skeletal`, `Brain - Cortex`, `Adipose - Subcutaneous`) %>%
  na.omit() %>%
  mutate("cluster" = clusters) %>%
  arrange(cluster)

cluster1_nDNA <- exprs_sub %>%
  filter(cluster == "1")

cluster2_nDNA <- exprs_sub %>%
  filter(cluster == "2") 

cluster3_nDNA <- exprs_sub %>%
  filter(cluster == "3")




# Clinical Phenotypic patterns

cluster1_subjects <- annotations_subject_phenotypes %>%
  filter(SUBJID %in% rownames(cluster1_nDNA)) %>%
  select(MHCOPD, MHCVD, MHHRTATT, MHHRTDISB, MHHTN, MHLVRDIS, MHPNMIAB, MHPNMIAB,
         MHSTRDLT, MHT1D, MHT2D, DTHVNT, MHSMKSTS, MHDRNKSTS, MHALS, MHALZDMT, MHALZHMR, MHARTHTS,
         MHASCITES, MHASTHMA, MHBCTINF, MHBLDDND, MHBLDOCNT, MHCANCER5, MHCANCERC, MHCANCERNM,
         MHCLLULTS, MHCLRD, MHCOCAINE5, MHCOUGHU, MHDLYSIS, MHDMNTIA, MHDPRSSN,MHENCEPHA, MHFLU,
         MHFNGINF, MHFVRU, MHGNRR12M, MHHEPBCT, MHHEPCCT, MHHEROIN, MHHGH, MHHIVCT, MHHMPHLIA, MHMENINA,
         MHMS, MHNEPH, MHNPHYS4W, MHOPPINF, MHORGNTP, MHOSTMYLTS, MHPLLABS, MHPRKNSN,MHRA,
         MHRNLFLR, MHTBHX, MHUREMIA, MHWKNSSU, TRORGNS) %>%
  mutate(MHSMKSTS = ifelse(is.na(MHSMKSTS), NA, ifelse(MHSMKSTS == "Yes", 1, 0))) %>%
  mutate(MHDRNKSTS = ifelse(is.na(MHDRNKSTS), NA, ifelse(MHDRNKSTS == "Yes", 1, 0)))


cluster2_subjects <- annotations_subject_phenotypes %>%
  filter(SUBJID %in% rownames(cluster2_nDNA)) %>%
  select(MHCOPD, MHCVD, MHHRTATT, MHHRTDISB, MHHTN, MHLVRDIS, MHPNMIAB, MHPNMIAB,
         MHSTRDLT, MHT1D, MHT2D, DTHVNT, MHSMKSTS, MHDRNKSTS, MHALS, MHALZDMT, MHALZHMR, MHARTHTS,
         MHASCITES, MHASTHMA, MHBCTINF, MHBLDDND, MHBLDOCNT, MHCANCER5, MHCANCERC, MHCANCERNM,
         MHCLLULTS, MHCLRD, MHCOCAINE5, MHCOUGHU, MHDLYSIS, MHDMNTIA, MHDPRSSN,MHENCEPHA, MHFLU,
         MHFNGINF, MHFVRU, MHGNRR12M, MHHEPBCT, MHHEPCCT, MHHEROIN, MHHGH, MHHIVCT, MHHMPHLIA, MHMENINA,
         MHMS, MHNEPH, MHNPHYS4W, MHOPPINF, MHORGNTP, MHOSTMYLTS, MHPLLABS, MHPRKNSN,MHRA,
         MHRNLFLR, MHTBHX, MHUREMIA, MHWKNSSU, TRORGNS) %>%
  mutate(MHSMKSTS = ifelse(is.na(MHSMKSTS), NA, ifelse(MHSMKSTS == "Yes", 1, 0))) %>%
  mutate(MHDRNKSTS = ifelse(is.na(MHDRNKSTS), NA, ifelse(MHDRNKSTS == "Yes", 1, 0)))

cluster3_subjects <- annotations_subject_phenotypes %>%
  filter(SUBJID %in% rownames(cluster3_nDNA)) %>%
  select(MHCOPD, MHCVD, MHHRTATT, MHHRTDISB, MHHTN, MHLVRDIS, MHPNMIAB, MHPNMIAB,
         MHSTRDLT, MHT1D, MHT2D, DTHVNT, MHSMKSTS, MHDRNKSTS, MHALS, MHALZDMT, MHALZHMR, MHARTHTS,
         MHASCITES, MHASTHMA, MHBCTINF, MHBLDDND, MHBLDOCNT, MHCANCER5, MHCANCERC, MHCANCERNM,
         MHCLLULTS, MHCLRD, MHCOCAINE5, MHCOUGHU, MHDLYSIS, MHDMNTIA, MHDPRSSN,MHENCEPHA, MHFLU,
         MHFNGINF, MHFVRU, MHGNRR12M, MHHEPBCT, MHHEPCCT, MHHEROIN, MHHGH, MHHIVCT, MHHMPHLIA, MHMENINA,
         MHMS, MHNEPH, MHNPHYS4W, MHOPPINF, MHORGNTP, MHOSTMYLTS, MHPLLABS, MHPRKNSN,MHRA,
         MHRNLFLR, MHTBHX, MHUREMIA, MHWKNSSU, TRORGNS) %>%
  mutate(MHSMKSTS = ifelse(is.na(MHSMKSTS), NA, ifelse(MHSMKSTS == "Yes", 1, 0))) %>%
  mutate(MHDRNKSTS = ifelse(is.na(MHDRNKSTS), NA, ifelse(MHDRNKSTS == "Yes", 1, 0)))


count_1 <- function(df) {
  result <- sapply(df, function(x) { sum(x == 1, na.rm = T)})
  
  return(result)
}

counted_cluster1 <- as.data.frame(count_1(cluster1_subjects))
counted_cluster2 <- as.data.frame(count_1(cluster2_subjects))
counted_cluster3<- as.data.frame(count_1(cluster3_subjects))

combined_counts <- cbind(counted_cluster1, counted_cluster2, counted_cluster3) %>%
  rename("Cluster 1" = `count_1(cluster1_subjects)`,
         "Cluster 2" = `count_1(cluster2_subjects)`,
         "Cluster 3" = `count_1(cluster3_subjects)`)

remove_rows_with_low_sum <- function(combined_counts, threshold = 10) {
  row_sums <- rowSums(combined_counts)
  df_filtered <- combined_counts[row_sums >= threshold, ]
  return(df_filtered)
}

# Apply the function to df_counts
df_counts_filtered <- remove_rows_with_low_sum(combined_counts)


sample_size <- c("Cluster 1" = length(rownames(cluster1_nDNA)), "Cluster 2" = length(rownames(cluster2_nDNA)), 
                 "Cluster 3" = length(rownames(cluster3_nDNA)))

# Function to divide each column by its corresponding divisor
divide_by <- function(df, sample_size) {
  df %>%
    mutate(across(all_of(names(sample_size)), ~ . / sample_size[cur_column()]))
}


combined_counted_percentage <- divide_by(df_counts_filtered, sample_size) %>%
  t() %>%
  as.data.frame() 


col_fun <- colorRamp2(c(-2, 0, 2), c("#311432", "#FFFFFF", "#008000"))

h <- Heatmap(scale(combined_counted_percentage), cluster_rows = F, col = col_fun,
             row_names_side = "left",       
             column_names_rot = 45,         
             show_column_dend = FALSE, 
             show_row_names = TRUE,
             show_column_names = T,
             heatmap_legend_param = list(title = "z-score Frequency %", title_gp = gpar(fontsize = 12, col = "black"), 
                                         title_position = "leftcenter-rot"))


draw(h)
