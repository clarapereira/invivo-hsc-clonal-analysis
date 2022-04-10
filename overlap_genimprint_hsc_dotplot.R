# this script overlaps our data with geneimprint dataset (imprinted list of genes)
# it filters our genes for the expressed ones (but we can change that by varying the var cpm_threshold - default is cpm_threshold = 10)
# outputs a dot plot with the overlaps and AI values per gene in our data


library(tidyverse)

source("R/random_hsc_functions.R")

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
# 1. add all metadata, including: 
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
# 2. add also LOH: 
# ==============================================================================================================================
# from WES: 
aici_table_long.expressed.lohwes <- aici_table_long.expressed %>% 
  addLohFromWes(
    loh_path = "../tables/wes/LOH_wes_light_thrcov_90.tsv"
  )
# # from E6: 
# aici_table_long.expressed.lohwes.lohrna <- aici_table_long.expressed.lohwes %>% 
#   addLohFromE6Rnaseq(
#     min  = 0.25,  
#     max = 0.75
#  )

# remove LOH: 
aici_table_long.expressed.noloh <- aici_table_long.expressed.lohwes %>% #nrow()
  dplyr::filter(LOH_ANYsamples == "not_loh")#; aici_table_long.expressed.noloh %>% nrow()

# how many genes are expressed in our dataset?
aici_table_long.expressed %>% 
  select(ID) %>% 
  distinct() %>% 
  nrow()

# and after removing LOH?
aici_table_long.expressed.noloh %>% 
  select(ID) %>% 
  distinct() %>% 
  nrow()

# how many of these are imprinted?
# all database: 
aici_table_long.expressed.noloh %>% 
  select(ID, gene, imprinted_status) %>% 
  distinct() %>% 
  group_by(imprinted_status) %>% 
  summarise(count = n())
# in B cells: 
aici_table_long.expressed.noloh %>% 
  filter(cell_type == "Bcells") %>% 
  select(ID, gene, imprinted_status) %>% 
  distinct() %>% 
  group_by(imprinted_status) %>% 
  summarise(count = n())
# in T cells: 
aici_table_long.expressed.noloh %>% 
  filter(cell_type == "Tcells") %>% 
  select(ID, gene, imprinted_status) %>% 
  distinct() %>% 
  group_by(imprinted_status) %>% 
  summarise(count = n())

# how do these imprinted look like in our data?
# with dots
aici_table_long.expressed.noloh %>% 
  filter(imprinted_status == "Imprinted") %>% 
  ggplot(aes(x=gene, y = AI, alpha = abundance, color = cell_type)) + 
  theme_light() +
  #geom_violin() +
  geom_point() +
  geom_jitter() +
  coord_flip() +
  facet_wrap(vars(clonality)) +
  #facet_grid(vars(cell_type), vars(clonality)) + 
  ggtitle(
    "Geneimprint 'Imprinted' genes expressed in our samples",
    subtitle = "(Lymplocyte populations from single-cell expanded HSCs)"
  )

ggsave(
  "../imprinted/plots/imprinted_genes_expressed_in_our_samples.pdf",
  width = 6,
  height = 5
)