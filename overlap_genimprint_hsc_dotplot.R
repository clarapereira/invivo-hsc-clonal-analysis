#!/usr/bin/Rscript
# this script overlaps our data with geneimprint dataset (imprinted list of genes)
# it filters our genes for the expressed ones (but we can change that by varying the var cpm_threshold - default is cpm_threshold = 10)
# outputs a dot plot with the overlaps and AI values per gene in our data


library(tidyverse)

source("R/random_hsc_functions.R")

# ==============================================================================================================================
# user-defined variables: 
# ==============================================================================================================================
dataset <- "all_replace_with_downsampled"
min_coverage <- 10    # coverage threshold
cpm_threshold <- 10   # abundance threshold (is different from coverage)
biomart_path <- "../tables/genes_biomaRt.tsv"   # already contains geneimprint data
bedfile <- "../asynchronous_replication_overlap/coordinates_mm10_exon.bed"
abundance_path <- "../abundance_edgeR/B_vs_T_v2/mean_abundance.tsv"
loh_wes_path <- "../tables/wes/LOH_wes_light_thrcov_90.tsv"
output_plot_path <- "../imprinted/plots"


# ==============================================================================================================================
# import data: 
# ==============================================================================================================================
aici_table_long.path <- paste0("../tables/",dataset,"/AICI_table_thrcov",min_coverage,dataset,"_longer.tsv")
aici_table_long <- data.table::fread(aici_table_long.path)
head(aici_table_long)

# ==============================================================================================================================
# 1. add all metadata, including: 
# ==============================================================================================================================
# - study metadata 
# - imprinted genes db
# - biomart db
# - expression data 
# - filter by genes by expression (choose threshold)
aici_table_long.expressed <- aici_table_long %>% 
  addAllMetaWrapper(
    biomart_path = biomart_path, 
    bedfile = bedfile,
    norm_path = abundance_path, 
    cpm_threshold = cpm_threshold
  ) %>% 
  # keep only chr of interest
  dplyr::filter(!grepl("^chrJ|^chrG|^chrMT|^chrY", chr)) ; head(aici_table_long.expressed)


# # retrieve as many imprinted genes as possible: 
# geneimprint.annotated <- data.table::fread("../tables/geneimprint.annotated.tsv") 
# 
# 
# aici_table_long.expressed <- aici_table_long.expressed.badImprint %>% 
#   dplyr::rename(imprinted_status.old = imprinted_status) %>% 
#   left_join(
#     geneimprint.annotated %>% dplyr::select(ensembl_id, imprinted_status),
#     by = c("ID" = "ensembl_id")
#   ) 
# aici_table_long.expressed %>% 
#   filter(imprinted_status.old == "Imprinted" | imprinted_status == "Imprinted" ) %>% 
#   dplyr::select(ID, gene, gene_name, imprinted_status.old, imprinted_status) %>% 
#   dplyr::distinct() #%>% 
#   #dplyr::group_by(imprinted_status.old) %>% 
#   #dplyr::group_by(imprinted_status) %>%
#   #dplyr::summarise(n())


aici_table_long.expressed %>% 
  dplyr::select(ID, imprinted_status, gene_name, gene) %>% 
  dplyr::distinct() %>% 
  dplyr::filter(is.na(gene_name))

# from_gtf %>% filter(gene_id == "ENSMUSG00000023795")
# from_gtf %>% filter(gene_name == "Pisd-ps2")
# aici_table_long.expressed %>% filter(gene == "Pisd-ps2") %>% distinct()
# ==============================================================================================================================
# 2. add also LOH: 
# ==============================================================================================================================
# 2.1 add from WES: 
aici_table_long.expressed.lohwes <- aici_table_long.expressed %>% 
  addLohFromWes(
    loh_path = loh_wes_path
  )
# # from E6:
# aici_table_long.expressed.lohwes.lohrna <- aici_table_long.expressed.lohwes %>%
#   addLohFromE6Rnaseq(
#     min  = 0.25,
#     max = 0.75
#  )

# 2.2 remove LOH: 
aici_table_long.expressed.noloh <- aici_table_long.expressed.lohwes %>% #nrow()
  dplyr::filter(LOH_ANYsamples == "not_loh")#; aici_table_long.expressed.noloh %>% nrow()

# =======================================
# Output some reports
# =======================================
message(
  "Total genes expressed in the dataset: ",
  aici_table_long.expressed %>% 
    dplyr::select(ID) %>% 
    dplyr::distinct() %>% 
    nrow(),
  "\nAfter LOH removal, remained: ",
  aici_table_long.expressed.noloh %>% 
    dplyr::select(ID) %>% 
    dplyr::distinct() %>% 
    nrow()
)

message(
  "How many of these are imprinted?"
)
# all database: 
aici_table_long.expressed.noloh %>% 
  dplyr::select(ID, gene, imprinted_status) %>% 
  dplyr::distinct() %>% 
  dplyr::group_by(imprinted_status) %>% 
  dplyr::summarise(count = n())
# in B cells: 
message(
  "How many of these are imprinted in B cells?"
)
aici_table_long.expressed.noloh %>% 
  dplyr::filter(cell_type == "Bcells") %>% 
  dplyr::select(ID, gene, imprinted_status) %>% 
  dplyr::distinct() %>% 
  dplyr::group_by(imprinted_status) %>% 
  dplyr::summarise(count = n())
# in T cells: 
message(
  "How many of these are imprinted in T cells?"
)
aici_table_long.expressed.noloh %>% 
  dplyr::filter(cell_type == "Tcells") %>% 
  dplyr::select(ID, gene, imprinted_status) %>% 
  dplyr::distinct() %>% 
  dplyr::group_by(imprinted_status) %>% 
  dplyr::summarise(count = n())

# ==============================================================================================================================
# 3. Plot
# ==============================================================================================================================
# how do these "imprinted" genes look like in our data?
# with dots

#tissueColors <- c("#FD6467", "#5B1A18")
#tissueColors <- c("pink", "brown")

# Box plot


aici_table_long.expressed.noloh %>% 
  dplyr::filter(imprinted_status == "Imprinted") %>% 
  ggplot(aes(x=gene_name, y = AI, alpha = abundance, color = cell_type)) + 
  theme_light() +
  #geom_violin() +
  geom_point() +
  geom_jitter() +
  coord_flip() +
  #scale_color_manual(values=tissueColors) +
  facet_wrap(vars(clonality)) +
  #facet_grid(vars(cell_type), vars(clonality)) + 
  ggtitle(
    "'Imprinted' genes expressed in B cells and T cells",
    subtitle = "(lymphocyte populations expanded from 1 HSCs in vivo)"
  )

ggsave(
  paste0(output_plot_path, "/imprinted_genes_expressed_cpm",cpm_threshold,"_in_hsc_samples.pdf"),
  width = 7,
  height = 9
)

message(
  "Plot saved as: ", paste0(output_plot_path, "/imprinted_genes_expressed_cpm",cpm_threshold,"_in_hsc_samples.pdf") 
)
