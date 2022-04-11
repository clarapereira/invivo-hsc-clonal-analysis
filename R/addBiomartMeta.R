
addBiomartMeta <- function(df_long = aici_table_long_meta, biomart_path = "../tables/genes_biomaRt.tsv"){
  #' @family this is a project-specific function (hsc rme project); this function is useful because its used through the project many times
  #' @description adds data previously extracted from biomart and other databases, left_joins by "ID"
  #' @biomart_path, path to table previously generated
  #' @df_long, a dataframe in long format containing a "ID" column
  #' @example aici_table_long_meta_biomart <- aici_table_long_meta %>% addBiomartMeta(biomart_path = "../tables/genes_biomaRt.tsv")
  #' 
  genes_biomaRt <- data.table::fread(biomart_path) %>% 
    dplyr::select(
      ensembl_gene_id, 
      external_gene_name, 
      chr, 
      imprinted.status
      )
  #
  genes_biomaRt <- dplyr::rename(
    genes_biomaRt, 
    ID = ensembl_gene_id, 
    gene = external_gene_name, 
    imprinted_status = imprinted.status
    )
  #
  genes_biomaRt$ID <- as.character(genes_biomaRt$ID)
  #
  df_long_biomart <- df_long %>% 
    left_join(genes_biomaRt, by = "ID") 
  #df_long_biomart %>% filter(is.na(gene)) %>% select(ID, gene) %>% distinct() # about 936 gene names missing???
  #
  head(df_long_biomart)
  #
  return(df_long_biomart)
}
