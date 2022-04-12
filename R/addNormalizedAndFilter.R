


# ==============================================================================================================================
# add normalized gene expression values and filter for CPM>10
# ==============================================================================================================================

addNormalizedAndFilter <- function(df_long = aici_table_long_meta_biomart, norm_path = "../abundance_edgeR/B_vs_T_v2/mean_abundance.tsv", cpm_threshold = 10){
  #' @family this is a project-specific function (hsc rme project); 
  #' @description adds normalized gene expression values and filter for CPM > threshold 
  #' @df_long, a dataframe in long format containing a "ID","sample" columns
  #' @norm_path path to file with EDGR normalised data, contains "ID","sample", "mean_abundance", "low_expr_genes", "std_abundance" columns
  #' @example aici_table_long.expressed <- aici_table_long_meta_biomart %>% addNormalizedAndFilter(norm_path = "../abundance_edgeR/B_vs_T_v2/mean_abundance.tsv", cpm_threshold = 10)
  #' 
  abundance <- data.table::fread(norm_path)
  abundance$mean_abundance <- as.numeric(round(abundance$mean_abundance, 2))
  abundance$std_abundance <- as.numeric(round(abundance$std_abundance, 2))
  #
  df_long.expressed <- df_long %>% 
    left_join(
      abundance, 
      by = c("ID","sample") 
      ) %>%
    filter(mean_abundance > cpm_threshold) %>%
    dplyr::rename(abundance = "mean_abundance") %>%
    dplyr::select(-low_expr_genes, -std_abundance)
  return(df_long.expressed)
}

