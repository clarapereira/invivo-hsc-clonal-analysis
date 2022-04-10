
addAllMetaWrapper <- function(df_long = aici_table_long,  biomart_path = "../tables/genes_biomaRt.tsv", norm_path = "../abundance_edgeR/B_vs_T_v2/mean_abundance.tsv", cpm_threshold = 10){
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
  # add normalized gene expression values and filter for CPM>10
  # ==============================================================================================================================
  df_long_meta_biomart.expressed <-  df_long_meta_biomart %>% 
    addNormalizedAndFilter(
      norm_path = norm_path, 
      cpm_threshold = cpm_threshold
    )
}
