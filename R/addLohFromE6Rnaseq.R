

addLohFromE6Rnaseq <- function(df_long = aici_table_long, min  = 0.25,  max = 0.75){
  #' Flag possible LOH genes within experiment 6 (E6) using RNA-seq data using threshold given by user
  #'
  #'
  df.E6_wider <- df_long %>% 
    filter(chr != "chrX")%>%
    #filter(abundance > 10)%>%
    filter(sample %in% c("E6.1_B", "E6.2_B", "E6.42_B", "E6.43_B") ) %>% 
    dplyr::select(ID, gene, chr, sample, sumCOV, matCOV, AI)%>%
    pivot_wider(
      id_cols = c(ID, gene, chr),
      names_from = sample,
      names_sep = "_",
      names_repair = "check_unique",
      values_from =`sumCOV`:`AI`)
  
  df.E6_wider.loh_rna <-  df.E6_wider %>%
    dplyr::mutate(
      LOH_E6_RNA = if_else(
        (AI_E6.1_B < min & AI_E6.2_B < min & AI_E6.42_B < min & AI_E6.43_B < min) |
          (AI_E6.1_B > max & AI_E6.2_B > max & AI_E6.42_B > max & AI_E6.43_B > max), 
        "TRUE", 
        "not_loh")
      ) %>%
    dplyr::select(ID, LOH_E6_RNA)
  
  df_long.loh_rna <- df_long %>% 
    left_join(
      df.E6_wider.loh_rna, 
      by = "ID"
      )
  df_long.loh_rna$LOH_E6_RNA  <- df_long.loh_rna$LOH_E6_RNA %>% replace_na("not_loh")
  
  return(df_long.loh_rna)
}
