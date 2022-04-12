

getGeneticallyBiased <- function(df_long = aici_table_long.expressed, min_impr  = 0.1, max_impr = 0.9){
  #' @family this is a project-specific function (hsc rme project);
  #' Identify genetically biased genesin our samples
  #' From OLD: "para cada gene vi se em todas as amostras poly e mono com desvio, se o desvio AI sempre AI > 0.9 ou sempre AI < 0.1"
  #' Note: I'm including the controls here!!. To avoid that, add "starts_with("E") to the selecting function across()
  #' here, I'm using all data (not only the significant)
  #' Criteria for inclusion:
  #' - at least 4 out of the 5 Tcell samples must have a parentally-specific biased AI value --> "genetically_biased"
  #' - at least 6 out of the 11 Bcell samples must have a parentally-specific biased AI value --> "genetically_biased"
  #' - if less than 4 Tcell or less than 6 Bcells have a parentally-specific biased AI value --> "possibly_biased"
  #' - none have a value --> "no_data"
  #' 
  df_long$AI <- as.numeric(round(df_long$AI, 2))
  # 1. recode data
  df_recode.wider <- df_long %>% 
    mutate(
      bias = case_when(
        is.na(AI) ~ "missing",
        AI > max_impr ~ "maternal",
        AI < min_impr ~ "paternal",
        T ~ "no_imprint"
      )
    ) %>% 
  # 2. make wider df
    dplyr::select( ID, gene, chr, sample, AI, bias) %>%
    pivot_wider(
      id_cols = c(ID, gene, chr),
      names_from = sample,
      names_sep = "_",
      names_repair = "check_unique",
      values_from = `AI`) %>% 
   # mutate(
   #   across(everything(), ~replace_na(.x, "missing"))
   # ) %>% 
    mutate(
      
      T_total_maternal = rowSums(across(ends_with("_T")) > max_impr, na.rm = TRUE),
      T_total_paternal = rowSums(across(ends_with("_T")) < min_impr, na.rm = TRUE),
      T_total_notbiased = rowSums(across(ends_with("_T")) > min_impr & across(ends_with("_T")) < max_impr, na.rm = TRUE),
      T_total_missing = rowSums(is.na(across(ends_with("_T"))) ),
      T_total_samples = T_total_paternal + T_total_maternal + T_total_notbiased + T_total_missing ,
      
      B_total_paternal = rowSums(across(ends_with("_B")) < min_impr, na.rm = TRUE),
      B_total_maternal = rowSums(across(ends_with("_B")) > max_impr, na.rm = TRUE),
      B_total_notbiased = rowSums(across(ends_with("_B")) > min_impr & across(ends_with("_B")) < max_impr, na.rm = TRUE),
      B_total_missing = rowSums(is.na(across(ends_with("_B"))) ),
      B_total_samples = B_total_paternal + B_total_maternal + B_total_notbiased + B_total_missing 

      ) %>% 
    mutate(
      imprinted_Tcells = case_when(
        T_total_samples == T_total_paternal | T_total_samples == T_total_maternal ~ "genetically_biased",
        ( (T_total_samples == T_total_paternal + T_total_missing) & T_total_missing == 1 ) | (( T_total_samples == T_total_maternal + T_total_missing ) & T_total_missing == 1 ) ~ "genetically_biased",
        T_total_samples == T_total_missing ~ "no_data",
        T_total_notbiased > 0 ~ "not_genetically_biased",
        (T_total_samples == T_total_paternal + T_total_missing) & T_total_missing > 1 | ( T_total_samples == T_total_maternal + T_total_missing ) & T_total_missing > 1 ~ "possibly_biased",
        T_total_paternal > 0 & T_total_maternal > 0 ~ "possibly_random",
        T ~ "possibly_random"
      ),
      imprinted_Bcells = case_when(
        B_total_samples == B_total_paternal | B_total_samples == B_total_maternal ~ "genetically_biased",
        ( (B_total_samples == B_total_paternal + B_total_missing) & B_total_missing == 1 ) | (( B_total_samples == B_total_maternal + B_total_missing ) & B_total_missing == 1) ~ "genetically_biased",
        ( (B_total_samples == B_total_paternal + B_total_missing) & B_total_missing == 2 ) | (( B_total_samples == B_total_maternal + B_total_missing ) & B_total_missing == 2) ~ "genetically_biased",
        ( (B_total_samples == B_total_paternal + B_total_missing) & B_total_missing == 3 ) | (( B_total_samples == B_total_maternal + B_total_missing ) & B_total_missing == 3) ~ "genetically_biased",
        ( (B_total_samples == B_total_paternal + B_total_missing) & B_total_missing == 4 ) | (( B_total_samples == B_total_maternal + B_total_missing ) & B_total_missing == 4) ~ "genetically_biased",
        ( (B_total_samples == B_total_paternal + B_total_missing) & B_total_missing == 5 ) | (( B_total_samples == B_total_maternal + B_total_missing ) & B_total_missing == 5) ~ "genetically_biased",
        B_total_samples == B_total_missing ~ "no_data",
        B_total_notbiased > 0 ~ "not_genetically_biased",
        ( (B_total_samples == B_total_paternal + B_total_missing) & B_total_missing > 5 ) | (( B_total_samples == B_total_maternal + B_total_missing ) & B_total_missing > 5 ) ~ "possibly_biased",
        B_total_paternal > 0 & B_total_maternal > 0 ~ "possibly_random",
        T ~ "possibly_random"
    )) %>% 
    select(-starts_with("T_"), -starts_with("B_"))
  
  
  return(df_recode.wider)
}



