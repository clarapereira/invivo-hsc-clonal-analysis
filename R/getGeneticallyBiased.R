

getGeneticallyBiased <- function(df_long = aici_table_long, min_impr  = 0.1, max_impr = 0.9){
  #' @family this is a project-specific function (hsc rme project);
  #' Identify genetically biased genesin our samples
  #' "para cada gene vi se em todas as amostras poly e mono com desvio, se o desvio AI sempre AI > 0.9 ou sempre AI < 0.1"
  #' here, I'm using all data (not only the significant)
  #' 
  
  # 1. make wider df
  df_wider <- df_long %>% 
    #filter(chr != "chrX")%>%
    #filter(abundance > 10)%>%
    #filter(sample %in% c("control_T", "E13.1_T", "E13.2_T", "E13.24_T", "E13.29_T", "control_B", "E13.1_B", "E13.2_B", "E15.2_B", "E13.24_B", "E13.29_B", "E15.10_B") ) %>% 
    dplyr::select(ID, gene, chr, sample, AI)%>%
    pivot_wider(
      id_cols = c(ID, gene, chr),
      names_from = sample,
      names_sep = "_",
      names_repair = "check_unique",
      values_from = `AI`)
  
  genetically_biased <- df_wider %>%
    dplyr::mutate(
      imprinted_Bcells = if_else(
        (control_B < min_impr & E13.1_B < min_impr & E13.2_B < min_impr & E15.2_B < min_impr & E13.24_B < min_impr & E13.29_B < min_impr & E15.10_B < min_impr & E6.1_B < min_impr & E6.2_B < min_impr & E6.42_B < min_impr & E6.43_B < min_impr) |
          (control_B > max_impr & E13.1_B > max_impr & E13.2_B > max_impr & E15.2_B > max_impr & E13.24_B > max_impr & E13.29_B > max_impr & E15.10_B > max_impr & E6.1_B > max_impr & E6.2_B > max_impr & E6.42_B > max_impr & E6.43_B > max_impr), 
        "genetically_biased", 
        "not_genetically_biased")
      ) %>%
    dplyr::mutate(
      imprinted_Tcells = if_else(
        (control_T < min_impr & E13.1_T < min_impr & E13.2_T < min_impr & E13.24_T < min_impr & E13.29_T< min_impr) |
          (control_T > max_impr & E13.1_T > max_impr & E13.2_T > max_impr & E13.24_T > max_impr & E13.29_T > max_impr), 
        "genetically_biased", 
        "not_genetically_biased")
      )
  return(genetically_biased)
}
