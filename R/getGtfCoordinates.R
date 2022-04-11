#'
#' bedfile, a file with the columns: V1-V5 in the following order: chr, start, end, strand, gene_id, gene_name
#' opens the file and extracts all chr except ones starting with letters: "^J|^G|MT|Y"
#' adds "chr" string to chr name/number
#'
#'

getGtfCoordinates <- function(bedfile = "../asynchronous_replication_overlap/coordinates_mm10_exon.bed"){
  require(tidyverse)
  coord_from_gtf <- data.table::fread(bedfile) %>% 
    # keep only autosomes and X chromosome: 
    dplyr::filter(!grepl("^J|^G|MT|Y", V1)) %>% 
    dplyr::mutate(V1b = paste0("chr", V1)) %>% 
    dplyr::select(V1b, V2,V3,V4,V5, V7) %>% 
    dplyr::rename(
      "seqnames" = "V1b",
      "start" = "V2",
      "end" = "V3",
      "strand" = "V7",
      "gene_id" = "V4",
      "gene_name" = "V5"
    )
  return(coord_from_gtf)
}

