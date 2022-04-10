

addLohFromWes <- function(df_long = aici_table_long, loh_path = "../tables/wes/LOH_wes_light_thrcov_90.tsv"){
  #'
  #'
  #'
  #'
  loh_wes <- data.table::fread(loh_path)%>%
    dplyr::select(ID, LOH_ANYsamples)
  
  df_long.loh_wes <- df_long %>% 
    left_join(
      loh_wes, 
      by = "ID"
      )
  df_long.loh_wes$LOH_ANYsamples  <- df_long.loh_wes$LOH_ANYsamples %>% replace_na("not_loh")
  
  return(df_long.loh_wes)
}