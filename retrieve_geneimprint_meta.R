#!/usr/bin/Rscript
# retrieve as many imprinted genes as possible

library(tidyverse)
source("R/random_hsc_functions.R")

genes_biomaRt <- data.table::fread("../tables/genes_biomaRt.tsv") %>% dplyr::select(ensembl_gene_id, external_gene_name, chr, imprinted.status)
genes_biomaRt <- dplyr::rename(genes_biomaRt, ID=ensembl_gene_id, gene=external_gene_name, imprinted_status=imprinted.status)
genes_biomaRt$ID <- as.character(genes_biomaRt$ID)
genes_biomaRt %>% 
  filter(!str_detect(chr, "PATCH")) %>% 
  dplyr::select(ID, imprinted_status) %>% 
  dplyr::distinct() %>% 
  dplyr::group_by(imprinted_status) %>% 
  dplyr::summarise(n())

imprinted_biomart <- genes_biomaRt %>% 
  filter(!str_detect(chr, "PATCH")) %>% 
  filter(imprinted_status == "Imprinted") #%>% View()

# get all gene names from GTF: 
from_gtf <- getGtfCoordinates(bedfile = "../asynchronous_replication_overlap/coordinates_mm10_exon.bed"); head(from_gtf)

# get all gene names from GTF: 
geneimprint <- readxl::read_excel("/Users/clarapereira/Dropbox/Boston_partners/BarretoLab/imprinted_loci/genimprint_semicurated.xlsx"); head(geneimprint )

#merge everything
geneimprint.annotated <- geneimprint %>% #colnames()
  left_join(
    from_gtf %>% filter(!str_detect(seqnames, "PATCH")) %>% dplyr::select(gene_id,  gene_name, strand) %>% dplyr::distinct(),
    by = c("Gene" = "gene_name")
  ) %>% #filter(!is.na(gene_id)) %>% dplyr::select(gene_id) %>% distinct() %>% nrow()
  left_join(
    genes_biomaRt %>% filter(!str_detect(chr, "PATCH")) ,
    by = c("Gene" = "gene")
  ) %>% #filter(!is.na(ID)) %>% dplyr::select(ID) %>% distinct() %>% nrow()
  #filter(is.na(gene_id) | is.na(ID)) %>% #View()
  mutate(
    ensembl_id = case_when(
      !is.na(ID) ~ ID,
      !is.na(gene_id) ~ gene_id
    )
  ) %>% 
  dplyr::select(-imprinted_status, -gene_id, -ID) %>% 
  dplyr::rename("imprinted_status" = "genimprint.com_status")
  
geneimprint.annotated %>% 
  filter(is.na(ensembl_id))
geneimprint.annotated %>%   #View()
  dplyr::select(imprinted_status, ensembl_id) %>% 
  dplyr::group_by(imprinted_status) %>% 
  dplyr::summarise(n())


write_delim(geneimprint.annotated, "../tables/geneimprint.annotated.tsv", delim = "\t")

from_gtf %>% filter(gene_name =="Peg3as")

