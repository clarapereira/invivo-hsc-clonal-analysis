library(tidyverse)
library(GenomicRanges)

source("R/random_hsc_functions.R")

# ==============================================================================================================================
# user-defined variables: 
# ==============================================================================================================================
dataset <- "all_replace_with_downsampled"
min_coverage <- 10
coverage_threshold <- min_coverage

plot_out_path <- "../asynchronous_replication_overlap/plots/abelson/"

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
aici_table_long.expressed.noloh <- aici_table_long.expressed %>% 
  left_join(
    ableson_loh,
    by = c("ID", "sample")
  ) %>% 
  filter(
    loh != "TRUE"
  )
head(aici_table_long.expressed.noloh)
nrow(aici_table_long.expressed.noloh)
nrow(aici_table_long.expressed)

# how many genes are expressed in abls dataset?
aici_table_long.expressed %>% 
 select(ID) %>% 
 distinct() %>% 
 nrow()

# and after loh removal?
aici_table_long.expressed.noloh %>% 
  select(ID) %>% 
  distinct() %>% 
  nrow()

# how many of these are imprinted?
aici_table_long.expressed.noloh %>% 
 select(ID, gene, imprinted_status) %>% 
 distinct() %>% 
 group_by(imprinted_status) %>% 
 summarise(count = n())
# how do these imprinted look like in abelson data?
# with dots
aici_table_long.expressed.noloh %>% 
 filter(imprinted_status == "Imprinted") %>% 
 ggplot(aes(x=gene, y = AI, alpha = abundance, color = sample)) +
 theme_light() +
 #geom_violin() +
 geom_point() +
 geom_jitter() +
 coord_flip() +
  ggtitle(
    "Imprinted genes expressed in abelson clones"
  )

ggsave(
  "../imprinted/plots/imprinted_genes_expressed_in_abelson.pdf",
  width = 6,
  height = 10
)

# ==============================================================================================================================
# Import Blumenfeld regions, after liftover mm9 --> mm10 (UCSC table browser, 5 April 2022)
# ==============================================================================================================================
blumenfeld_bed <- "/Users/clarapereira/Downloads/13.Literature_Mar2022/Chromosomal_coordination/liftover/hglft_genome_51b0_c6c2f0.bed"
get_lifted <- data.table::fread(blumenfeld_bed) %>% 
 dplyr::rename(
  "chromosome" = "V1",
  "start" = "V2",
  "end" = "V3",
  "assynch" = "V4",
  "regionID" = "V5"
 )

# ==============================================================================================================================
# get coordinates from our annotation 
# ==============================================================================================================================

bed = "/Users/clarapereira/Dropbox/Boston_partners/HSC_clones_shared_folder/asynchronous_replication_overlap/coordinates_mm10_exon.bed"

coord_from_gtf <- data.table::fread(bed) %>% 
 # keep only autosomes and X chromosome: 
 dplyr::filter(!grepl("^J|^G|MT|Y", V1)) %>% 
 dplyr::mutate(V1b = paste0("chr", V1)) %>% 
 dplyr::select(V1b, V2,V3,V4,V5, V7) %>% 
 dplyr::rename(
  "seqnames" = "V1b",
  "start" = "V2",
  "end" = "V3",
  "strand" = "V7",
  "gene_id" = "V4",
  "gene_name" = "V5"
 ); head(coord_from_gtf)

# ==============================================================================================================================
# Make the overlaps: 
# ==============================================================================================================================
# 1. transform our dataframes into granges objects: 
coord_from_gtf.granges <- makeGRangesFromDataFrame(
 coord_from_gtf, 
 keep.extra.columns=T
 )
get_lifted.granges <- makeGRangesFromDataFrame(
 get_lifted, 
 keep.extra.columns=T
 )
# find overlaps function for Dave Tang
intersect_bed <- function(a, b){
 library(GenomicRanges)
 my_hit <- findOverlaps(a, b)
 my_df <- cbind(as.data.frame(a[queryHits(my_hit)]),
   as.data.frame(b[subjectHits(my_hit)]))
}

# ==============================================================================================================================
# How many of our annotation version genes are in the Blumenfeld regions?
# ==============================================================================================================================
all_genes_in_blumenfeld <- intersect_bed(
 coord_from_gtf.granges, 
 get_lifted.granges
 ) ; all_genes_in_blumenfeld %>% head()

colnames(all_genes_in_blumenfeld) <- c("seqnames", "start", "end" , "width", "strand",  "gene_id",
   "gene_name", "Chromosome", "start", "end",  "width", "strand",  
   "assynch",  "regionID" )
assy_genes <- all_genes_in_blumenfeld %>% 
 dplyr::select(seqnames, gene_id, gene_name, assynch, regionID) %>% 
 dplyr::distinct(); head(assy_genes)

assy_genes %>% 
 select(gene_name) %>% 
 distinct() %>% 
 nrow()

# ==============================================================================================================================
# How many of these genes are expressed in our abelson dataset?
# ==============================================================================================================================
assy_genes.expressed.data <- assy_genes %>% 
 left_join(
  aici_table_long.expressed.noloh %>% select(ID) %>% distinct() %>% mutate(expressed = "yes"),
  by = c("gene_id" = "ID")
 ); head(assy_genes.expressed.data)
#
assy_genes.expressed.data %>% 
 na.omit() %>% 
 nrow()

# ==============================================================================================================================
# How do they look like?
# ==============================================================================================================================
assy_genes.expressed.data.abelson <- assy_genes.expressed.data %>% select(-gene_name, -seqnames) %>% 
 na.omit() %>% 
 left_join(
  aici_table_long.expressed.noloh, 
  by = c("gene_id" = "ID" )
 ); head(assy_genes.expressed.data.abelson)

assy_genes.expressed.data.abelson %>% 
 dplyr::filter(grepl("^chr11|^chr18", 
                     seqnames)) %>% 
 #dplyr::filter(seqnames == "chr11") %>% 
 ggplot(aes(x=gene_name, y = AI, alpha = abundance, color = sample)) + 
  scale_y_continuous(limits = c(0, 1))  +
  theme_light() +
  geom_point() +
  geom_jitter() +
  coord_flip() +
  ggtitle("chr11 & chr18 abelson expressed genes",
          subtitle = "within Blumenfeld et. al regions")

ggsave(
  paste0(plot_out_path, "/blumenfeld_overlap_chr11chr18genes_abelson.pdf"),
  width = 6,
  height = 5
)

my_chr <- paste0("chr", c(1:19, "X"))

for (i in 1:length(my_chr)){
  
  assy_genes.expressed.data.abelson %>% #head()
    dplyr::filter(seqnames == my_chr[i]) %>% 
    ggplot(aes(x=gene_name, y = AI, alpha = abundance, color = sample)) +
    scale_y_continuous(limits = c(0, 1))  +
    theme_light() +
    geom_point() +
    geom_jitter() +
    coord_flip() +
    #facet_wrap(vars(clonality)) +
    ggtitle(paste0(my_chr[i], " abelson expressed genes"),
            subtitle = "within Blumenfeld et. al regions"
            )
  ggsave(
    paste0(paste0(plot_out_path, "/blumenfeld_overlap_",my_chr[i],"genes.pdf")),
    width = 6,
    height = 5
  )
}

 
 
 
 
 
 
 