library(tidyverse)



gtex <- read.delim(gzfile("bulk-gex_v8_rna-seq_GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz"), skip = 2)

annotations_sample_attributes <- read_tsv("GTEx_Analysis_2017-06-05_v8_Annotations_GTEx_Analysis_2017-06-05_v8_Annotations_SampleAttributesDS.tsv")
annotations_subject_phenotypes <- read_tsv("GTEx_Analysis_2017-06-05_v8_Annotations_GTEx_Analysis_2017-06-05_v8_Annotations_SubjectPhenotypesDS.tsv")


annotations_sample_attributes <- annotations_sample_attributes %>%
  filter(SMAFRZE == "RNASEQ") %>%
  mutate(SUBJID = substring(SAMPID,1,10))

annotations_sample_attributes <- annotations_sample_attributes %>%
  mutate(SUBJID_tmp = case_when(
    substr(SUBJID, nchar(SUBJID), nchar(SUBJID)) == "-" ~ substr(SUBJID, 1, nchar(SUBJID)-1),
    .default = SUBJID
  )) %>%
  dplyr::select(-SUBJID) %>%
  rename(SUBJID = SUBJID_tmp)


annotations_sample_attributes <- annotations_sample_attributes %>%
  rename("X" = SAMPID)
annotations_sample_attributes$X <- gsub("\\-", ".", annotations_sample_attributes$X)

exprs <- gtex %>%
  dplyr::select(-Name) %>%
  rename("Gene" = Description) 

#filter out genes with 0 for every value
exprs <- exprs %>%
  filter(rowSums(across(everything(), ~ . == 0)) != (ncol(exprs)-1))


#read in Mitocarta genes
mitocarta <- readxl::read_xls(here::here("HumanMitoCarta3_0.xls"), sheet = 2) 
mitogenes <- unique(mitocarta$Symbol)
# Mitocarta symbol, synonyms and ensembl
mitogenes_df <- mitocarta %>%
  dplyr::select(Symbol, Synonyms, EnsemblGeneID_mapping_version_20200130) %>%
  dplyr::rename(Symbol_MC = Symbol, Synonyms_MC = Synonyms, Ensembl_MC = EnsemblGeneID_mapping_version_20200130) %>%
  mutate(Ensembl_MC = gsub("|", " ", Ensembl_MC, fixed=TRUE)) %>%
  separate(Ensembl_MC, into = c("a", "b", "c", "d", "e", "f", "g", "h", "i", "j"), sep = " ") %>%
  pivot_longer(cols = c("a", "b", "c", "d", "e", "f", "g", "h", "i", "j"), values_to = "Ensembl_MC") %>%
  dplyr::select(-name) %>%
  filter(!is.na(Ensembl_MC)) %>%
  mutate(Synonyms_MC = gsub("|", " ", Synonyms_MC, fixed=TRUE)) %>%
  separate(Synonyms_MC, into = as.character(seq(1:15)), sep = " ") %>%
  pivot_longer(cols = as.character(seq(1:15)), values_to = "Synonyms_MC") %>%
  dplyr::select(-name) %>%
  filter(!is.na(Synonyms_MC)) %>%
  mutate(merge = Ensembl_MC) 

geneset <- gtex %>%
  select(Name, Description) 

gtex_genes <- geneset %>%
  mutate(name2 = sub("\\..*", "",Name)) %>%
  dplyr::rename(Ensembl_gtex_iso = Name, 
                Ensembl_gtex = name2,
                Symbol_gtex = Description) %>%
  dplyr::mutate(merge = Ensembl_gtex) 

match_by_ensembl <- full_join(gtex_genes, mitogenes_df, by = "merge") %>%
  filter(!is.na(Ensembl_MC)) %>%
  filter(!is.na(Ensembl_gtex)) %>%
  unique() 


genes_not_found <- mitogenes[which(!mitogenes %in% match_by_ensembl$Symbol_MC)]


mitocarta_genes <- unique(match_by_ensembl$Symbol_gtex)

mtDNA_genes <- c("MT-ATP6", "MT-ATP8", "MT-CO1", "MT-CO2", "MT-CO3", "MT-CYB",
                 "MT-ND1", "MT-ND2", "MT-ND3", "MT-ND4", "MT-ND4L", "MT-ND5",
                 "MT-ND6", "MT-RNR1", "MT-RNR2", "MT-TA", "MT-TC", "MT-TD",
                 "MT-TE", "MT-TF", "MT-TG", "MT-TH", "MT-TI", "MT-TK", "MT-TL1", "MT-TL2",
                 "MT-TM", "MT-TN", "MT-TP", "MT-TQ", "MT-TR", "MT-TS1", "MT-TS2",
                 "MT-TT", "MT-TV", "MT-TW", "MT-TY")


# filter for nuclear mitocarta genes
exprs_nDNA <- exprs %>%
  as.data.frame() %>%
  filter(Gene %in% mitocarta_genes) %>%
  filter(!Gene %in% mtDNA_genes)


exprs_nDNA <- exprs_nDNA %>%
  column_to_rownames("Gene") %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("X")



# sum all nuclear mitochondrial genes in each sample
exprs_nDNA <- exprs_nDNA %>%
  pivot_longer(colnames(exprs_nDNA[2:ncol(exprs_nDNA)])) %>%
  group_by(X) %>%
  mutate(summed_exprs = sum(value)) %>%
  dplyr::select(X, summed_exprs) %>%
  ungroup() %>%
  unique()

# join with sample attributes data and pivot wider
exprs_nDNA <- exprs_nDNA %>%
  full_join(annotations_sample_attributes, by = "X") %>%
  dplyr::select(summed_exprs, SMTSD, SUBJID, SMRIN) %>%
  filter(SMRIN >= 5.5) %>%
  dplyr::select(-SMRIN) %>%
  pivot_wider(names_from = SMTSD, values_from = summed_exprs) %>%
  column_to_rownames("SUBJID")



# all total nDNA
exprs_total_nDNA <- exprs %>%
  filter(!Gene %in% mtDNA_genes) #remove mtDNA transcripts from each sample

# sum all remaining nuclear transcipts
exprs_total_nDNA <- exprs_total_nDNA %>%
  dplyr::select(-Gene)%>%
  map_dbl(~sum(.,na.rm=T)) %>%
  as.data.frame()


exprs_total_nDNA <- exprs_total_nDNA %>%
  rownames_to_column("X") %>%
  rename("summed_exprs" = ".")

exprs_total_nDNA <- exprs_total_nDNA %>%
  full_join(annotations_sample_attributes, by = "X") %>%
  dplyr::select(summed_exprs, SMTSD, SUBJID, SMRIN) %>%
  filter(SMRIN>=5.5) %>%
  dplyr::select(-SMRIN) %>%
  unique() %>%
  pivot_wider(names_from = SMTSD, values_from = summed_exprs) %>%
  column_to_rownames("SUBJID")

#ensure dataframes in same order
exprs_nDNA <- exprs_nDNA[rownames(exprs_total_nDNA), colnames(exprs_total_nDNA)]

mito_nDNA_percent <- (exprs_nDNA/exprs_total_nDNA)*100

# filter for tissues of interest
mito_nDNA_percent <- mito_nDNA_percent %>%
  dplyr::select( - `Cervix - Ectocervix`, - `Cervix - Endocervix`, - `Fallopian Tube`, -`Cells - Cultured fibroblasts`,
                 - Bladder, -`Cells - EBV-transformed lymphocytes`, - `Small Intestine - Terminal Ileum`, -Spleen,
                 -`Cells - EBV-transformed lymphocytes`, - `Kidney - Medulla`
  ) 





# all transcripts

exprs_all_transcripts <- exprs

exprs_all_transcripts <- exprs_all_transcripts %>%
  as.data.frame() %>%
  dplyr::select(-Gene) %>%
  map_dbl(~sum(.,na.rm=T)) %>%
  as.data.frame()

exprs_all_transcripts <- exprs_all_transcripts %>%
  rownames_to_column("X") %>%
  rename("summed_exprs" = ".")

exprs_all_transcripts <- exprs_all_transcripts %>%
  full_join(annotations_sample_attributes, by = "X") %>%
  dplyr::select(summed_exprs, SMTSD, SUBJID, SMRIN) %>%
  filter(SMRIN>=5.5) %>%
  dplyr::select(-SMRIN) %>%
  unique() %>%
  pivot_wider(names_from = SMTSD, values_from = summed_exprs) %>%
  column_to_rownames("SUBJID")


# mtDNA transcripts

exprs_mtDNA <- exprs %>%
  filter(Gene %in% mtDNA_genes) %>%
  column_to_rownames("Gene") %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("X")


exprs_mtDNA <- exprs_mtDNA %>%
  pivot_longer(colnames(exprs_mtDNA[2:ncol(exprs_mtDNA)])) %>%
  group_by(X) %>%
  mutate(summed_exprs = sum(value)) %>%
  dplyr::select(X, summed_exprs) %>%
  ungroup() %>%
  unique()


exprs_mtDNA <- exprs_mtDNA %>%
  full_join(annotations_sample_attributes, by = "X") %>%
  dplyr::select(summed_exprs, SMTSD, SUBJID, SMRIN) %>%
  filter(SMRIN>=5.5) %>%
  dplyr::select(-SMRIN) %>%
  unique() %>%
  pivot_wider(names_from = SMTSD, values_from = summed_exprs) %>%
  column_to_rownames("SUBJID")

#ensure dataframes in same order
exprs_mtDNA <- exprs_mtDNA[rownames(exprs_all_transcripts), colnames(exprs_all_transcripts)]

mtDNA_percent <- (exprs_mtDNA/exprs_all_transcripts)*100

#filter for tissues of interest
mtDNA_percent <- mtDNA_percent %>%
  dplyr::select( - `Cervix - Ectocervix`, - `Cervix - Endocervix`, - `Fallopian Tube`, -`Cells - Cultured fibroblasts`,
                 - Bladder, -`Cells - EBV-transformed lymphocytes`, - `Small Intestine - Terminal Ileum`, -Spleen,
                 -`Cells - EBV-transformed lymphocytes`, - `Kidney - Medulla`
  ) 
