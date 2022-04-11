library(tidyverse)
library(GenomicRanges)

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
  dplyr::filter(LOH_ANYsamples == "not_loh" & LOH_E6_RNA == "not_loh")
aici_table_long.expressed.noloh %>% nrow()

# how many genes are expressed in our dataset?
aici_table_long.expressed %>% 
  dplyr::select(ID) %>% 
  dplyr::distinct() %>% 
  nrow()

# and after removing LOH?
aici_table_long.expressed.noloh %>% 
  dplyr::select(ID) %>% 
  dplyr::distinct() %>% 
  nrow()

# how many of these are imprinted?
aici_table_long.expressed.noloh %>% 
  dplyr::select(ID, gene, imprinted_status) %>% 
  dplyr::distinct() %>% 
  dplyr::group_by(imprinted_status) %>% 
  dplyr::summarise(count = n())

# how do these imprinted look like in our data?
aici_table_long.expressed.noloh %>% #head()
  #filter(clonality == "non_monoclonal") %>% 
  dplyr::filter(imprinted_status == "Imprinted") %>% 
  #na.omit() %>% 
  ggplot(aes(x=gene, y = AI)) +
  theme_light() +
  geom_violin() +
  coord_flip() +
  facet_wrap(vars(clonality))
# with dots
aici_table_long.expressed.noloh %>% 
  dplyr::filter(imprinted_status == "Imprinted") %>% 
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
# get cordinates from our annotation (via Invoke a System Command)
# ==============================================================================================================================
# 1. the following command will filter only CDS, and extract a bed file from our GTF file; it will save the bed file, and we will open it afterwards
#    the command is : awk 'BEGIN{OFS="\t";} $3=="CDS" {print $1,$4-1,$5,$10,$16,$6,$7}' "/Users/clarapereira/Dropbox/Barreto_lab/reference/release-68/Mus_musculus.GRCm38.68.gtf" | sed 's/"//g' | sed 's/;//g' | sortBed > "/Users/clarapereira/Dropbox/Boston_partners/HSC_clones_shared_folder/asynchronous_replication_overlap/coordinates_mm10.bed"
input_gtf <- "/Users/clarapereira/Dropbox/Barreto_lab/reference/release-68/Mus_musculus.GRCm38.68.gtf"
output_bed = "../asynchronous_replication_overlap/coordinates_mm10_exon.bed"
system(
  #"awk 'BEGIN{OFS=\"\t\";} $3==\"CDS\" {print $1,$4-1,$5,$10,$16,$6,$7}' \"/Users/clarapereira/Dropbox/Barreto_lab/reference/release-68/Mus_musculus.GRCm38.68.gtf\" | sed 's/\"//g' | sed 's/;//g' | sortBed > \"/Users/clarapereira/Dropbox/Boston_partners/HSC_clones_shared_folder/asynchronous_replication_overlap/coordinates_mm10.bed\"" # this works
  paste0("awk 'BEGIN{OFS=\"\t\";} $3==\"exon\" {print $1,$4-1,$5,$10,$16,$6,$7}' ",input_gtf," | sed 's/\"//g' | sed 's/;//g' | sortBed > ",output_bed )
  )
# 2. now we can get the data we' are interested in're looking for: 
coord_from_gtf <- data.table::fread(output_bed) %>% 
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
  )



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
  my_df  <- cbind(as.data.frame(a[queryHits(my_hit)]),
                  as.data.frame(b[subjectHits(my_hit)]))
}

# ==============================================================================================================================
# How many of our annotation version genes are in the Blumenfeld regions?
# ==============================================================================================================================
all_genes_in_blumenfeld <- intersect_bed(
  coord_from_gtf.granges, 
  get_lifted.granges
  ) ; all_genes_in_blumenfeld %>% head()

colnames(all_genes_in_blumenfeld) <- c("seqnames",  "start",     "end" ,      "width",     "strand",    "gene_id",
                         "gene_name", "Chromosome",  "start",     "end",       "width",     "strand",   
                         "assynch",   "regionID" )
assy_genes <- all_genes_in_blumenfeld %>% 
  dplyr::select(seqnames, gene_id, gene_name, assynch, regionID) %>% 
  dplyr::distinct(); head(assy_genes)

assy_genes %>% 
  dplyr::select(gene_name) %>% 
  dplyr::distinct() %>% 
  nrow()

# ==============================================================================================================================
# How many of these genes are expressed in our dataset?
# ==============================================================================================================================
assy_genes.expressed.data <- assy_genes %>% 
  dplyr::left_join(
    aici_table_long.expressed.noloh %>% dplyr::select(ID) %>% dplyr::distinct() %>% dplyr::mutate(expressed = "yes"),
    by = c("gene_id" = "ID")
  ); head(assy_genes.expressed.data)
#
assy_genes.expressed.data %>% 
  na.omit() %>% 
  nrow()

# ==============================================================================================================================
# Is any of these genes one of our RME from the B-cell clonal HSC analysis?
# ==============================================================================================================================
rme_hsc <- c(
  "Aldh4a1",
  "Apba3",
  "Dpp7",
  "Fam32a",
  "Flnb",
  "Gnpda2",
  "Igkv6-25",
  "Kdm8",
  "Pacsin1",
  "Pkp3",
  "Plekha8",
  "Slc38a7",
  "Trim5",
  "Zfp873"
)

assy_genes.expressed.data %>% 
  dplyr::filter(gene_name %in% rme_hsc)

assy_genes %>% 
  dplyr::filter(gene_name == "Vmn1r211")

# ==============================================================================================================================
# How many of the genes in the Blumenfeld regions are imprinted according with geneimprint?
# ==============================================================================================================================

assy_genes.imprint.data <- assy_genes.expressed.data %>% #nrow()
  dplyr::rename(ID = gene_id) %>% 
  addBiomartMeta(biomart_path = "../tables/genes_biomaRt.tsv"); head(assy_genes.imprint.data)
#assy_genes.imprint.data %>% filter(is.na(gene))
assy_genes.imprint.data %>% 
  group_by(imprinted_status) %>% 
  summarise(count = n())


assy_genes.imprint.data %>% filter(imprinted_status == "Imprinted")


# ==============================================================================================================================
# How do all genes look like, chromosome by chromosome
# ==============================================================================================================================
assy_genes.expressed.data.hsc <- assy_genes.expressed.data %>% 
  na.omit() %>% 
  left_join(
    aici_table_long.expressed.noloh, 
    by = c("gene_id" = "ID" )
  ); head(assy_genes.expressed.data.hsc)

my_chr <- paste0("chr", c(1:19, "X"))

for (i in 1:length(my_chr)){
  
  assy_genes.expressed.data.hsc %>% #filter(is.na(gene))
    #dplyr::filter(grepl("^chr11|^chr18|^chr19", seqnames)) %>% 
    dplyr::filter(seqnames == my_chr[i]) %>% 
    ggplot(aes(x=gene_name, y = AI, alpha = abundance, color = cell_type)) +
    theme_light() +
    #geom_violin() +
    geom_point() +
    geom_jitter() +
    coord_flip() +
    facet_wrap(vars(clonality)) +
    #facet_grid(vars(cell_type), vars(clonality))
    #facet_grid( vars(seqnames),vars(clonality))
    ggtitle(
      paste0(my_chr[i], " expressed genes within Blumenfeld et. al regions",
      subtitle = "(Lymplocyte populations from single-cell expanded HSCs)")
    )
  ggsave(
    paste0("../asynchronous_replication_overlap/plots/blumenfeld_overlap_",my_chr[i],"genes.pdf"),
    width = 6,
    height = 5
  )
}


# How do "our" "biased" genes overlap with the blumenfeld regions?
genetically_biased <- data.table::fread("/Users/clarapereira/Dropbox/Boston_partners/HSC_clones_shared_folder/tables/imprinted/gene_lists/genetically_biased_cutof10.tsv")
aici_table_long.expressed.genetically_biased <- aici_table_long.expressed.noloh %>% 
  left_join(
    genetically_biased,
    by = "ID"
  ); head(aici_table_long.expressed.genetically_biased )


assy_genes.expressed.genetically_biased <- assy_genes.expressed.data %>% 
  na.omit() %>% 
  left_join(
    aici_table_long.expressed.genetically_biased, 
    by = c("gene_id" = "ID" )
  ) %>% 
  filter(imprinted_Bcells == "genetically_biased" | imprinted_Tcells == "genetically_biased"); head(assy_genes.expressed.genetically_biased)

assy_genes.expressed.genetically_biased %>% 
  #filter(imprinted_status == "Imprinted") %>% 
  ggplot(aes(x=gene_name, y = AI, alpha = abundance, color = cell_type)) +
  theme_light() +
  #geom_violin() +
  geom_point() +
  geom_jitter() +
  coord_flip() +
  facet_wrap(vars(clonality)) +
  #facet_grid(vars(cell_type), vars(clonality))
  ggtitle(
    "Genetically biased genes (.1) within Blumenfeld et. al regions",
    subtitle = "(Lymplocyte populations from single-cell expanded HSCs)"
    )

ggsave(
  "../asynchronous_replication_overlap/plots/blumenfeld_overlap_genetically_biased_cutof10.pdf",
  width = 6,
  height = 3
    )
