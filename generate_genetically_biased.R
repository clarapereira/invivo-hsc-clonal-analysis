library(tidyverse)

source("/Users/clarapereira/Dropbox/Boston_partners/HSC_clones_shared_folder/repository/random_hsc_functions.R")


# ==============================================================================================================================
# user-defined variables: 
# ==============================================================================================================================
dataset <- "all_replace_with_downsampled"
min_coverage <- 10
coverage_threshold <- min_coverage
#
biascutof <- 10 
min_impr  <- 0.1 
max_impr <- 0.9
output_folder.lists <- "../tables/imprinted/gene_lists/"
output_folder.dfs <- "../imprinted/tables/"
# ==============================================================================================================================
# import data: 
# ==============================================================================================================================
aici_table_long.path <- paste0("../tables/",dataset,"/AICI_table_thrcov",min_coverage,dataset,"_longer.tsv")
aici_table_long <- data.table::fread(aici_table_long.path)
head(aici_table_long)

# ==============================================================================================================================
# add all metadata: 
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
# classify genetically biased, biascutof: (output is wider)
# ==============================================================================================================================
aici_table_long.expressed.biased.biascutof <- aici_table_long.expressed %>% 
  getGeneticallyBiased(min_impr  = min_impr, max_impr = max_impr); head(aici_table_long.expressed.biased)
aici_table_long.expressed.biased.biascutof %>% 
  write_delim(paste0(output_folder.dfs, "genetically_biased_biascutof",biascutof,"_dfwider.tsv"), delim = "\t")


potentially_imprinted.biascutof <- aici_table_long.expressed.biased %>%
  dplyr::select(
    ID, 
    imprinted_Bcells, 
    imprinted_Tcells
    ) %>% 
  filter(
    imprinted_Bcells == "genetically_biased" | imprinted_Tcells == "genetically_biased"
    ); head(potentially_imprinted.biascutof10)
# how many?
potentially_imprinted.biascutof %>% nrow()
potentially_imprinted.biascutof %>% 
  write_delim(paste0(output_folder.lists, "genetically_biased_biascutof",biascutof,"_genelist.tsv"), delim = "\t")




