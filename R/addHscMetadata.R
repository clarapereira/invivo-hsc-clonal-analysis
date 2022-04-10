
addHscMetadata <- function(df_long){
  #' @family this is a project-specific function
  #' @description adds dataset metadata to data table in long format, left_joins by "sample"
  #' @df_long, a dataframe in long format containing a "sample" column
  #' @example aici_table_long_meta <- aici_table_long %>% addHscMetadata()

  sample <- c("control_T", "E13.1_T", "E13.2_T", "E13.24_T", "E13.29_T", "control_B", "E13.1_B", "E13.2_B", "E15.2_B", "E13.24_B", "E13.29_B", "E15.10_B", "E6.1_B", "E6.2_B", "E6.42_B", "E6.43_B")	
  cell_type <- c("Tcells", "Tcells", "Tcells", "Tcells", "Tcells", "Bcells", "Bcells","Bcells","Bcells","Bcells","Bcells","Bcells","Bcells","Bcells","Bcells","Bcells")
  clonality <- c ("non_monoclonal", "non_monoclonal" , "non_monoclonal", "monoclonal", "monoclonal", "non_monoclonal", "non_monoclonal", "non_monoclonal", "non_monoclonal", "monoclonal", "monoclonal", "monoclonal", "non_monoclonal", "non_monoclonal", "monoclonal", "monoclonal")
  
  metadata <- data.frame(sample, clonality, cell_type)
  
  df_long_meta <- left_join(df_long, metadata, by = "sample")
  head(df_long_meta)
  return(df_long_meta)
}