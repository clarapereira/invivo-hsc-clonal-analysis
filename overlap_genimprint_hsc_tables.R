#!/usr/bin/Rscript
# this script overlaps our data with geneimprint dataset (imprinted list of genes)
# it filters our genes for the expressed ones (but we can change that by varying the var cpm_threshold - default is cpm_threshold = 10)


library(tidyverse)
source("R/random_hsc_functions.R")

# some stats on our data: 
# Number imprinted genes recoverd from Tucci (2019) --> 167 
# Total number imprinted genes expressed (TPM >10) from Tucci (2019) --> 66
# Number imprinted genes recoverd from geneimprint --> 83 
# Total number imprinted genes expressed from geneimprint --> 22


# ==============================================================================================================================
# user-defined variables: 
# ==============================================================================================================================
dataset <- "all_replace_with_downsampled"
min_coverage <- 10    # coverage threshold
cpm_threshold <- 10   # abundance threshold (is different from coverage) 
biomart_path <- "../tables/genes_biomaRt.tsv"   # already contains geneimprint data
bedfile <- "../asynchronous_replication_overlap/coordinates_mm10_exon.bed"
abundance_path <- "../abundance_edgeR/B_vs_T_v2/mean_abundance.tsv"
loh_wes_path <- "../tables/wes/LOH_wes_light_thrcov_90.tsv"
#
biascutof <- 15
min_impr  <- 0.15
max_impr <- 0.85
output_plot_path <- "../imprinted/plots"
output_folder.lists <- "../tables/imprinted/gene_lists/"
output_folder.dfs <- "../imprinted/tables/"

# ==============================================================================================================================
# import data, only B cells: 
# ==============================================================================================================================
aici_table_long.path <- paste0("../tables/",dataset,"/AICI_table_thrcov",min_coverage,dataset,"_longer.tsv")
aici_table_long <- data.table::fread(aici_table_long.path)
head(aici_table_long)

# ==============================================================================================================================
# 1. add all metadata, including: 
# ==============================================================================================================================
# - study metadata 
# - imprinted genes db
# - biomart db
# - expression data 
# - filter by genes by expression (choose threshold)
aici_table_long.expressed <- aici_table_long %>% 
  addAllMetaWrapper(
    biomart_path = biomart_path, 
    bedfile = bedfile,
    norm_path = abundance_path, 
    cpm_threshold = cpm_threshold
  ) %>% 
  # keep only chr of interest
  dplyr::filter(!grepl("^chrJ|^chrG|^chrMT|^chrY", chr)) ; head(aici_table_long.expressed)

# export stats on imprinted genes: 
aici_table_long.expressed %>% 
  select(ID, imprinted_status) %>% 
  distinct() %>% 
  group_by(imprinted_status) %>% 
  summarise(count = n())

# ==============================================================================================================================
# 2. add also LOH data: 
# ==============================================================================================================================
# 2.1 add from WES: 
aici_table_long.expressed.lohwes <- aici_table_long.expressed %>% 
  addLohFromWes(
    loh_path = loh_wes_path
  )
# 2.2 from E6:
aici_table_long.expressed.lohwes.lohrna <- aici_table_long.expressed.lohwes %>%
  addLohFromE6Rnaseq(
    min  = 0.25,
    max = 0.75
 )

# 2.3 Remove LOH:
aici_table_long.expressed.lohwes.lohrna.nolohwes <- aici_table_long.expressed.lohwes.lohrna %>% 
  filter(LOH_ANYsamples != TRUE )# %>% 
  #filter(LOH_E6_RNA != TRUE & cell_type == "Bcells")

# ==============================================================================================================================
# 3. classify genetically biased with biascutof, and add to the previous:
# ==============================================================================================================================
# 3.1 classify:
aici_table_long.expressed.biased.biascutof <- aici_table_long.expressed.lohwes.lohrna.nolohwes %>% 
  getGeneticallyBiased(min_impr  = min_impr, max_impr = max_impr); head(aici_table_long.expressed.biased.biascutof)
#aici_table_long.expressed.biased.biascutof %>% View()
 # filter(gene == "Igf2r") %>% View()

  
# 3.2 add classification:
aici_table_long.expressed.lohwes.lohrna.biased.biascutof <- aici_table_long.expressed.lohwes.lohrna %>% 
  left_join(
    aici_table_long.expressed.biased.biascutof %>% dplyr::select(ID, imprinted_Bcells, imprinted_Tcells),
    by = "ID"
    ); head(aici_table_long.expressed.lohwes.lohrna.biased.biascutof)

aici_table_long.expressed.lohwes.lohrna.biased.biascutof %>% 
  filter(gene == "Igf2r")

# ==============================================================================================================================
# 4. Filter B cells table and make it wider: 
# ==============================================================================================================================

makeTableWStats <- function(tissue = "Bcells", imprint_column = "imprinted_Bcells"){
  
  aici_table.wider <- aici_table_long.expressed.lohwes.lohrna.biased.biascutof  %>% 
    filter(cell_type == tissue ) %>%
    filter(abundance > 10) %>%             # redundant
    dplyr::filter(.data[[imprint_column]] == "genetically_biased" | imprinted_status == "Imprinted") %>% 
    #filter(str_detect(imprinted_Bcells, "genetically_biased|possibly_biased") | imprinted_status == "Imprinted") %>% 
    filter(LOH_ANYsamples != "TRUE")%>%
    dplyr::select(ID, gene, gene_name, chr, sample, AI, .data[[imprint_column]], imprinted_status) %>%
    distinct() %>%
    pivot_wider(
      id_cols = c(ID, gene, gene_name, chr, .data[[imprint_column]], imprinted_status),
      names_from = sample,
      names_sep = "_",
      names_repair = "check_unique",
      values_from = `AI`)
  
  stat <- aici_table_long.expressed.lohwes.lohrna.biased.biascutof %>% 
    filter(abundance > 10 ) %>%
    filter(cell_type == tissue) %>%
    dplyr::filter(.data[[imprint_column]] == "genetically_biased" | imprinted_status == "Imprinted") %>% 
    #dplyr::filter( (str_detect(.data[[imprint_column]], "genetically_biased") | str_detect(imprinted_status,"Imprinted") )) %>% 
    #filter(str_detect(imprinted_Bcells, "genetically_biased|possibly_biased") | imprinted_status == "Imprinted") %>% 
    filter(LOH_ANYsamples != "TRUE" ) %>%
    group_by(ID)%>%
    summarise(Mean_AI = round(mean(AI, na.rm = TRUE), 2), Std_AI = round(sd(AI, na.rm = TRUE), 2)) %>%
    ungroup()
  
  stat_abundance <- aici_table_long.expressed.lohwes.lohrna.biased.biascutof %>% 
    filter(abundance > 10 ) %>%
    filter(cell_type == tissue ) %>%
    dplyr::filter(.data[[imprint_column]] == "genetically_biased" | imprinted_status == "Imprinted") %>% 
    #filter(str_detect(imprinted_Bcells, "genetically_biased|possibly_biased") | imprinted_status == "Imprinted") %>% 
    filter(LOH_ANYsamples != "TRUE") %>%
    group_by(ID) %>%
    summarise(Mean_abundance = round(mean(abundance, na.rm = TRUE), 2), Std_abundance = round(sd(abundance, na.rm = TRUE), 2) ) %>%
    ungroup()
  
  table <- aici_table.wider %>% 
    left_join(stat, 
              by = "ID") %>% 
    left_join(stat_abundance, 
              by = "ID") 
  
  return(table)
}

B_table <- makeTableWStats(tissue = "Bcells", imprint_column = "imprinted_Bcells")
T_table <- makeTableWStats(tissue = "Tcells", imprint_column = "imprinted_Tcells")

B_table.plus <- B_table %>%  
  left_join(
    T_table %>% filter(imprinted_status == "Imprinted") %>% mutate(imprinted_B_T = "yes") %>% select(ID, imprinted_B_T),
    by = "ID"
    ) %>% 
  left_join(
    T_table %>% filter(imprinted_Tcells == "genetically_biased") %>% mutate(biased_B_T = "yes") %>% select(ID, biased_B_T),
    by = "ID"
  ) 

B_table %>% select(ID, imprinted_status) %>% distinct() %>% group_by(imprinted_status) %>% summarise(count = n()) %>% 
  write_delim(paste0("../tables/imprinted/stats/stats_imprinted_B_biascutoff",biascutof,"_TMM",cpm_threshold,".tsv"), delim = "\t")
B_table %>% select(ID, imprinted_Bcells) %>% distinct() %>% group_by(imprinted_Bcells) %>% summarise(count = n()) %>% 
  write_delim(paste0("../tables/imprinted/stats/stats_genetically_biased_B_biascutoff",biascutof,"_TMM",cpm_threshold,".tsv"), delim = "\t")


# Order by chromosome: 
chrOrder <-c(paste0("chr", (1:22)),"chrX")
biasedOrder1 <- c("genetically_biased", "possibly_biased", "no_data", "not_genetically_biased", "possibly_random")
imprintedOrder1 <- c("ND", "Not Imprinted", "Provisional Data", "Imprinted")

B_table.plus$chr <- factor(B_table.plus$chr, chrOrder, ordered=TRUE)
B_table.plus$imprinted_Bcells <- factor(B_table.plus$imprinted_Bcells, biasedOrder1, ordered=TRUE)
B_table.plus$imprinted_status <- factor(B_table.plus$imprinted_status, imprintedOrder1, ordered=TRUE)

B_table.plus %>% arrange(imprinted_Bcells, imprinted_status, chr) %>% head()
B_table.plus %>% nrow()
#B_table %>% arrange(chr) %>% View()

write_delim(B_table.plus %>% 
              dplyr::arrange(imprinted_Bcells, imprinted_status, chr) , paste0("../tables/imprinted/imprinted_B_biascutoff",biascutof,"_TMM",cpm_threshold,".tsv"), delim = "\t")

###################################################################################################################################
#### NOW EXPORT TABLE IN PDF: 
###################################################################################################################################
library(gt)
library(glue)


colnames(B_table) <- c(
  "ID", "Gene.vOld", "Gene", "Chr", "imprinted_Bcells", "imprinted_status",
  "Control", "E13.1", "E13.2", "E15.2", "E13.24",
  "E13.29", "E15.10", "E6.1", "E6.2", "E6.42", "E6.43",     
  "Mean(AI)", "Std(AI)", "Mean", "Std")

B_table.tidy <- B_table %>% 
  dplyr::mutate_if(is.numeric, ~round(., 2)) %>% 
  dplyr::mutate(
    Biased = case_when(
      imprinted_Bcells == "genetically_biased" ~ "yes",
      imprinted_Bcells == "not_genetically_biased" ~ "no",
      imprinted_Bcells == "possibly_biased" ~ "maybe",
      imprinted_Bcells == "possibly_random" ~ "no",
      imprinted_Bcells == "no_data" ~ "-"
    ),
    Imprinted = case_when(
      imprinted_status == "ND" ~ "-",
      imprinted_status == "Imprinted" ~ "yes",
      imprinted_status == "Not Imprinted" ~ "no", 
      imprinted_status == "Provisional Data" ~ "no", 
      T ~ imprinted_status
    )
  ) 

biasedOrder <- c("yes", "maybe", "-", "no")
#biasedOrder <- c("Biased", "Not biased", "No criteria for classification")
B_table.tidy$Biased <- factor(B_table.tidy$Biased, biasedOrder, ordered=TRUE)

B_table.image <- B_table.tidy %>% 
  dplyr::arrange(Biased, Imprinted, Chr) %>% 
  dplyr::select(Gene, Chr, Biased, Imprinted, -imprinted_status, -ID,-imprinted_Bcells, Control:Std ) %>% 
  #gt(rowname_col = "Imprinted", groupname_col = "Biased") %>% 
  gt() %>% 
  cols_align(
    align = "center",
    columns = c(Chr:Std)
  ) %>% 
  tab_spanner(label = "AI", columns = matches("Control|^E")) %>%
  tab_spanner(label = "Abundance", columns = matches("Mean$|Std$")) %>%
  tab_header(
    title = md("Supplementary Table X. List of genetically biased and imprinted genes in B samples."),
    subtitle = md("")
  ) %>%
  tab_source_note(
    md("Genes expressed in our samples (>10 TMM-normalized counts) were classified as “genetically biased” if AI values were less than 0.1 or higher than 0.9 in at least 11 samples in B cells (in blue). \n
    A list of imprinted genes was retrieved from https://www.geneimprint.com/site/genes-by-species (in orange). Out of the 126 imprinted genes listed, we could identify 82 in our annotation file \n
    (ftp://ftp.ensembl.org/pub/release-68/gtf/Mus_musculus.GRCm38.68.gtf), but only 19 were “expressed” in our samples (i.e., TMM-normalized counts >10). Of the 19 imprinted genes detected \n
    in our samples, only one (Zrsr1, in bold) had been detected as “genetically biased” in our samples. The abundance mean (mean of TMM-normalized counts) of all analyzed genes is 115.03 ± 477.90 \n
    in B cells. Low expressed genes (<10 TMM-normalized counts) and genes suspect of LOH were removed from study."))
B_table.image
#library(webshot)
gtsave(B_table.image, paste0("../tables/imprinted/imprinted_B_biascutoff",biascutof,"_TMM",cpm_threshold,".pdf"))
########################################################################################################################


##################################
#### T cells: 
##################################

T_table %>% select(ID, imprinted_status) %>% distinct() %>% group_by(imprinted_status) %>% summarise(count = n()) %>% 
  write_delim(paste0("../tables/imprinted/stats/stats_imprinted_T_biascutoff",biascutof,"_TMM",cpm_threshold,".tsv"), delim = "\t")
T_table %>% select(ID, imprinted_Tcells) %>% distinct() %>% group_by(imprinted_Tcells) %>% summarise(count = n()) %>% 
  write_delim(paste0("../tables/imprinted/stats/stats_genetically_biased_T_biascutoff",biascutof,"_TMM",cpm_threshold,".tsv"), delim = "\t")

T_table.plus <- T_table %>%  
  left_join(
    B_table %>% filter(imprinted_status == "Imprinted") %>% mutate(imprinted_B_T = "yes") %>% select(ID, imprinted_B_T),
    by = "ID"
  ) %>% 
  left_join(
    B_table %>% filter(imprinted_Bcells == "genetically_biased") %>% mutate(biased_B_T = "yes") %>% select(ID, biased_B_T),
    by = "ID"
  ) 

# Order

T_table.plus$chr <- factor(T_table$chr, chrOrder, ordered=TRUE)
T_table.plus$imprinted_Tcells <- factor(T_table.plus$imprinted_Tcells, biasedOrder1, ordered=TRUE)
T_table.plus$imprinted_status <- factor(T_table.plus$imprinted_status, imprintedOrder1, ordered=TRUE)


T_table.plus %>% 
  dplyr::arrange(imprinted_Tcells, imprinted_status, chr) %>% head()
T_table %>% nrow()
#B_table %>% arrange(chr) %>% View()

write_delim(T_table.plus %>% 
              dplyr::arrange(imprinted_Tcells, imprinted_status, chr), paste0("../tables/imprinted/imprinted_T_biascutoff",biascutof,"_TMM",cpm_threshold,".tsv"), delim = "\t")


###################################################################################################################################
#### NOW EXPORT TABLE IN PDF: 
###################################################################################################################################

colnames(T_table) <- c(
  "ID", "Gene.vOld", "Gene", "Chr", "imprinted_Tcells", "imprinted_status",
  "Control", "E13.1", "E13.2", "E13.24", "E13.29",
  "Mean(AI)", "Std(AI)", "Mean", "Std")

T_table.tidy <- T_table %>% 
  dplyr::mutate_if(is.numeric, ~round(., 2)) %>% 
  dplyr::mutate(
    Biased = case_when(
      imprinted_Tcells == "genetically_biased" ~ "yes",
      imprinted_Tcells == "not_genetically_biased" ~ "no",
      imprinted_Tcells == "possibly_biased" ~ "maybe",
      imprinted_Tcells == "possibly_random" ~ "no",
      imprinted_Tcells == "no_data" ~ "-"
    ),
    Imprinted = case_when(
      imprinted_status == "ND" ~ "-",
      imprinted_status == "Imprinted" ~ "yes",
      imprinted_status == "Not Imprinted" ~ "no", 
      imprinted_status == "Provisional Data" ~ "no", 
      T ~ imprinted_status
    )
  ) 

biasedOrder <- c("yes", "maybe", "-", "no")
#biasedOrder <- c("Biased", "Not biased", "No criteria for classification")
T_table.tidy$Biased <- factor(T_table.tidy$Biased, biasedOrder, ordered=TRUE)

T_table.image <- T_table.tidy %>% 
  dplyr::arrange(Biased, Imprinted, Chr) %>% 
  dplyr::select(Gene, Chr, Biased, Imprinted, -imprinted_status, -ID,-imprinted_Tcells, Control:Std ) %>% 
  #gt(rowname_col = "Imprinted", groupname_col = "Biased") %>% 
  gt() %>% 
  cols_align(
    align = "center",
    columns = c(Chr:Std)
  ) %>% 
  tab_spanner(label = "AI", columns = matches("Control|E13.1|E13.2|E13.24|E13.29")) %>%
  tab_spanner(label = "Abundance", columns = matches("Mean$|Std$")) %>%
  tab_header(
    title = md("Supplementary Table X. List of genetically biased and imprinted genes in T samples."),
    subtitle = md("")
  ) %>%
  tab_source_note(
    md("Genes expressed in our samples (>10 TMM-normalized counts) were classified as “genetically biased” if AI values were less than 0.1 or higher than 0.9 in at least 5 samples for T cells (in blue). \n
    A list of imprinted genes was retrieved from https://www.geneimprint.com/site/genes-by-species (in orange). Out of the 126 imprinted genes listed, we could identify 82 in our annotation file \n
    (ftp://ftp.ensembl.org/pub/release-68/gtf/Mus_musculus.GRCm38.68.gtf), but only 21 were “expressed” in our samples (i.e., TMM-normalized counts >10). Of the 21 imprinted genes detected in \n
    our samples, only two (Igf2r and Zrsr1, in bold) had been detected as “genetically biased” in our samples. The abundance mean (mean of TMM-normalized counts) of all analyzed genes is 101.65 ± 218.74 \n
    in T cells. Low expressed genes (<10 TMM-normalized counts) and genes suspect of LOH were removed from study."))

T_table.image

gtsave(T_table.image, paste0("../tables/imprinted/imprinted_T_biascutoff",biascutof,"_TMM",cpm_threshold,".pdf"))
###################################################################################################################################












