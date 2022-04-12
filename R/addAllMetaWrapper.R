
addAllMetaWrapper <- function(df_long = aici_table_long,  biomart_path = "../tables/genes_biomaRt.tsv", bedfile = "../asynchronous_replication_overlap/coordinates_mm10_exon.bed", norm_path = "../abundance_edgeR/B_vs_T_v2/mean_abundance.tsv", cpm_threshold = 10){
  #' @family this is a project-specific function (hsc rme project);
  #' 
  # ==============================================================================================================================
  # add metadata: 
  # ==============================================================================================================================
  df_long_meta <- df_long %>% addHscMetadata()
  # ==============================================================================================================================
  # add biomart and imprinting data: 
  # ==============================================================================================================================
  df_long_meta_biomart <- df_long_meta %>% 
    addBiomartMeta(
      biomart_path = biomart_path
    )
  
  # ==============================================================================================================================
  # add gene_names from GTF file
  # ==============================================================================================================================
  from_gtf <- getGtfCoordinates(bedfile = bedfile)
  df_long_meta_biomart_gtf <- df_long_meta_biomart %>% 
    dplyr::left_join(
      from_gtf %>% dplyr::select(gene_id,  gene_name, strand, seqnames) %>% dplyr::distinct(),
      by = c("ID" = "gene_id")
    ) %>% 
    mutate(
      chr = case_when(
        !is.na(seqnames) ~ seqnames,
        !is.na(chr) ~ chr,
        T ~ seqnames
      )
    )
  
  # ==============================================================================================================================
  # correct genemprint annotation
  # ==============================================================================================================================
  #geneimprint.annotated <- data.table::fread("../tables/geneimprint.annotated.tsv") 
  imprinted.annotated <- data.table::fread("../tables/gtf.mm10.v68.geneimprint.tucci2919.annotated.tsv") 
  df_long_meta_biomart_gtf_imprintscorr <-  df_long_meta_biomart_gtf %>% 
    dplyr::rename(imprinted_status.old = imprinted_status) %>% 
    left_join(
      imprinted.annotated %>% dplyr::select(gene_id, imprinted_status, expressed_allele) %>% rename("ensembl_id" = "gene_id"),
      by = c("ID" = "ensembl_id")
    ) 
  df_long_meta_biomart_gtf_imprintscorr$imprinted_status  <- df_long_meta_biomart_gtf_imprintscorr$imprinted_status  %>% replace_na("ND")
  
  # ==============================================================================================================================
  # add normalized gene expression values and filter for CPM>10
  # ==============================================================================================================================
  df_long_meta_biomart_gtf.expressed <-  df_long_meta_biomart_gtf_imprintscorr %>% 
    addNormalizedAndFilter(
      norm_path = norm_path, 
      cpm_threshold = cpm_threshold
    )
}
