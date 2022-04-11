

# this acript exposrts a bed file using bash / system commands
# also exposrts a 


# ==============================================================================================================================
# get cordinates from our annotation (via Invoke a System Command)
# ==============================================================================================================================
# 1. the following command will filter only CDS, and extract a bed file from our GTF file; it will save the bed file, and we will open it afterwards
#    the command is : awk 'BEGIN{OFS="\t";} $3=="CDS" {print $1,$4-1,$5,$10,$16,$6,$7}' "/Users/clarapereira/Dropbox/Barreto_lab/reference/release-68/Mus_musculus.GRCm38.68.gtf" | sed 's/"//g' | sed 's/;//g' | sortBed > "/Users/clarapereira/Dropbox/Boston_partners/HSC_clones_shared_folder/asynchronous_replication_overlap/coordinates_mm10.bed"
input_gtf <- "/Users/clarapereira/Dropbox/Barreto_lab/reference/release-68/Mus_musculus.GRCm38.68.gtf"
output_bed.cds = "../asynchronous_replication_overlap/coordinates_mm10_cds.bed"
output_bed.exon = "../asynchronous_replication_overlap/coordinates_mm10_exon.bed"
system(
  #"awk 'BEGIN{OFS=\"\t\";} $3==\"CDS\" {print $1,$4-1,$5,$10,$16,$6,$7}' \"/Users/clarapereira/Dropbox/Barreto_lab/reference/release-68/Mus_musculus.GRCm38.68.gtf\" | sed 's/\"//g' | sed 's/;//g' | sortBed > \"/Users/clarapereira/Dropbox/Boston_partners/HSC_clones_shared_folder/asynchronous_replication_overlap/coordinates_mm10.bed\"" # this works
  paste0("awk 'BEGIN{OFS=\"\t\";} $3==\"CDS\" {print $1,$4-1,$5,$10,$16,$6,$7}' ",input_gtf," | sed 's/\"//g' | sed 's/;//g' | sortBed > ",output_bed.cds )
)
system(
  paste0("awk 'BEGIN{OFS=\"\t\";} $3==\"exon\" {print $1,$4-1,$5,$10,$16,$6,$7}' ",input_gtf," | sed 's/\"//g' | sed 's/;//g' | sortBed > ",output_bed.exon )
)

# 2. now we can get the data we' are interested in're looking for: 
coord_from_gtf <- data.table::fread(output_bed.exon) %>% 
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

# aici_table_long.expressed.noloh %>% 
#   left_join(
#     coord_from_gtf, 
#     by = c("ID" = "gene_id")
#   ) %>% 
#   filter(is.na(gene))
