#! #!/usr/bin/Rscript
#

library(tidyverse)
source("R/random_hsc_functions.R")
source("src/rme_abelson_aistdev.R")



# ==============================================================================================================================
# user-defined variables: 
# ==============================================================================================================================
dataset <- "all_replace_with_downsampled"
min_coverage <- 10
coverage_threshold <- min_coverage

out_path <- "../coordination/abelson/"


# ==============================================================================================================================
# create dir
# ==============================================================================================================================
dir.create(file.path(out_path), recursive = TRUE, showWarnings = FALSE)


# ==============================================================================================================================
# import abl data: 
# ==============================================================================================================================
abl_aici_table_long.path <- paste0("../tables/abelson/AICI_table_abl_fitCovThr30_thrcov_",min_coverage,"_longer.tsv")
abl_aici_table_long <- data.table::fread(abl_aici_table_long.path)
head(abl_aici_table_long)

# ==============================================================================================================================
# add all metadata: 
# ==============================================================================================================================

aici_table_long.expressed <- abl_aici_table_long %>% 
  addAllMetaWrapper(
    biomart_path = "../tables/genes_biomaRt.tsv", 
    norm_path = "../abundance_edgeR/abelson2/mean_abundance.tsv", 
    imprinted_path = "../tables/gtf.mm10.v68.geneimprint.tucci2919.annotated.tsv",
    cpm_threshold = 10
  ) %>% 
  # keep only chr of interest
  dplyr::filter(!grepl("^chrJ|^chrG|chrMT|chrY", chr)) ; head(aici_table_long.expressed)


# add LOH info: 
ableson_loh <- data.table::fread("../tables/abelson/LOH_abl_thrcov_10_longer.tsv")
# remove LOH: 
aici_table_long.expressed.noloh <- aici_table_long.expressed %>% 
  left_join(
    ableson_loh,
    by = c("ID", "sample")
  ) %>% 
  filter(
    loh != "TRUE"
  )


# ==============================================================================================================================
# filter our rme_abelson_aistdev genes, 
# classify direction,
# export coordination table
# ==============================================================================================================================
aici_table_long.expressed.noloh.coord <- aici_table_long.expressed.noloh %>% #head()
  filter( gene_name %in% rme_abelson_aistdev) %>% 
  #filter(clonality == "monoclonal") %>% 
  #filter(cell_type == "Bcells") %>% 
  distinct() %>% 
  mutate(direction = case_when(
    AI > 0.5 ~ "M",
    AI < 0.5 ~ "P",
    T ~ "no_bias"
  )) %>% 
  select(gene_name, sample, direction) %>% 
  pivot_wider(
    id_cols = c(sample),
    names_from = gene_name,
    names_sep = "_",
    names_repair = "check_unique",
    values_from = direction) 


aici_table_long.expressed.noloh.coord %>% 
  write_delim(
    paste0(out_path, "coordination_abelson.tsv"), 
    delim = "\t"
  )


# by chromosome: 

my_chr <- paste0("chr", c(1:19, "X"))

for (i in 1:length(my_chr)){
  
  bychr <- aici_table_long.expressed.noloh %>% #head()
    dplyr::filter(seqnames == my_chr[i]) %>%
    filter( gene_name %in% rme_abelson_aistdev) %>% 
    #filter(clonality == "monoclonal") %>% 
    #filter(cell_type == "Bcells") %>% 
    distinct() %>% 
    mutate(direction = case_when(
      AI > 0.5 ~ "M",
      AI < 0.5 ~ "P",
      T ~ "no_bias"
    )) %>% 
    select(gene_name, sample, direction) %>% 
    pivot_wider(
      id_cols = c(sample),
      names_from = gene_name,
      names_sep = "_",
      names_repair = "check_unique",
      values_from = direction) 
  
  bychr %>% 
    write_delim(
      paste0(out_path, "coordination_abelson",my_chr[i],".tsv"), 
      delim = "\t"
    )
  
  message(my_chr[i])
  print(bychr)
}







