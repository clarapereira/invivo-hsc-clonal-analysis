library(tidyverse)

source("/Users/clarapereira/Dropbox/Boston_partners/HSC_clones_shared_folder/repository/random_hsc_functions.R")

# ==============================================================================================================================
# user-defined variables: 
# ==============================================================================================================================
dataset <- "all_replace_with_downsampled"
min_coverage <- 10
coverage_threshold <- min_coverage


# ==============================================================================================================================
# import data: 
# ==============================================================================================================================
aici_table_long.path <- paste0("../tables/",dataset,"/AICI_table_thrcov",min_coverage,dataset,"_longer.tsv")
aici_table_long <- data.table::fread(aici_table_long.path)
head(aici_table_long)

# ==============================================================================================================================
# add all metadata, including: 
# - study metadata 
# - imprinted genes db
# - biomart db
# - expression data 
# - filter by genes by expression (choose threshold)
# ==============================================================================================================================

aici_table_long.expressed <- aici_table_long %>% 
  addAllMetaWrapper(
    biomart_path = "../tables/genes_biomaRt.tsv", 
    norm_path = "../abundance_edgeR/B_vs_T_v2/mean_abundance.tsv", 
    cpm_threshold = 10
  ) %>% 
  # keep only chr of interest
  dplyr::filter(!grepl("^chrJ|^chrG|^chrMT|^chrY", chr)) ; head(aici_table_long.expressed)

# ==============================================================================================================================
# now, add also LOH: 
# ==============================================================================================================================
aici_table_long.expressed.lohwes <- aici_table_long.expressed %>% 
  addLohFromWes(
    loh_path = "../tables/wes/LOH_wes_light_thrcov_90.tsv"
  )
aici_table_long.expressed.lohwes.lohrna <- aici_table_long.expressed.lohwes %>% 
  addLohFromE6Rnaseq(
    min  = 0.25,  
    max = 0.75
  )
