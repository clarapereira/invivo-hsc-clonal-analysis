#! #!/usr/bin/Rscript
#

library(tidyverse)
source("R/random_hsc_functions.R")
source("src/rme_hsc_bcell.R")



# ==============================================================================================================================
# user-defined variables: 
# ==============================================================================================================================
dataset <- "all_replace_with_downsampled"
min_coverage <- 10
coverage_threshold <- min_coverage

out_path <- "../coordination/hsc_expanded/"


# ==============================================================================================================================
# create dir
# ==============================================================================================================================
dir.create(file.path(out_path), recursive = TRUE, showWarnings = FALSE)


# ==============================================================================================================================
# import data: 
# ==============================================================================================================================
aici_table_long.path <- paste0("../tables/",dataset,"/AICI_table_thrcov",min_coverage,dataset,"_longer.tsv")
aici_table_long <- data.table::fread(aici_table_long.path)
head(aici_table_long)

# ==============================================================================================================================
# add all metadata, remove all potential LOH: 
# ==============================================================================================================================
aici_table_long.expressed <- aici_table_long %>% 
  addAllMetaWrapper(
    biomart_path = "../tables/genes_biomaRt.tsv", 
    bedfile = "../asynchronous_replication_overlap/coordinates_mm10_exon.bed",
    norm_path = "../abundance_edgeR/B_vs_T_v2/mean_abundance.tsv", 
    imprinted_path = "../tables/gtf.mm10.v68.geneimprint.tucci2919.annotated.tsv",
    cpm_threshold = 10
  ) %>% 
  # keep only chr of interest
  dplyr::filter(!grepl("^chrJ|^chrG|^chrMT|^chrY", chr)) ; head(aici_table_long.expressed)

aici_table_long.expressed.lohwes <- aici_table_long.expressed %>% 
  addLohFromWes(
    loh_path = "../tables/wes/LOH_wes_light_thrcov_90.tsv"
  )
aici_table_long.expressed.lohwes.lohrna <- aici_table_long.expressed.lohwes %>% 
  addLohFromE6Rnaseq(
    min  = 0.25,  
    max = 0.75
  )
# remove LOH: 
aici_table_long.expressed.noloh <- aici_table_long.expressed.lohwes.lohrna %>% #nrow()
  #dplyr::filter(LOH_ANYsamples == "not_loh" & LOH_E6_RNA == "not_loh")
  dplyr::filter(LOH_ANYsamples == "not_loh" & LOH_E6_RNA == "not_loh")
aici_table_long.expressed.noloh %>% nrow()


# ==============================================================================================================================
# filter our rme_hsc_bcell genes, 
# classify direction,
# export coordination table
# ==============================================================================================================================

aici_table_long.expressed.noloh.coord <- aici_table_long.expressed.noloh %>% #head()
  filter( gene_name %in% rme_hsc_bcell) %>% 
  filter(clonality == "monoclonal") %>% 
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
    paste0(out_path, "coordination_monoclonal_BT.tsv"), 
    delim = "\t"
  )


# by chromosome: 

my_chr <- paste0("chr", c(1:19, "X"))

for (i in 1:length(my_chr)){
  
  bychr <- aici_table_long.expressed.noloh %>% #head()
    dplyr::filter(seqnames == my_chr[i]) %>%
    filter( gene_name %in% rme_hsc_bcell) %>% 
    filter(clonality == "monoclonal") %>% 
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
  
  message(my_chr[i])
  print(bychr)
}






aici_table_wider.stats_poly %>% 
  filter(chr!="chrX") %>% 
  filter(Std_AI_abl_clone > 0.15) %>% 
  select(ID, gene) %>% 
  distinct() %>% 
  #nrow()
  #View() 
  select(gene) %>% 
  write_delim()





