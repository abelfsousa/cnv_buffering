# Abel Sousa 2018

# PhD project
# Study of Control of Protein Abundance using Multi-omics Data from Cancer Samples



# --



suppressMessages(library(tidyverse))

options(bitmapType = "cairo")





# import protein attenuation data
protein_attenuation <- read_tsv("./files/protein_attenuation.txt")


# import genome-wide protein associations with CNV and RNA
# filtered by genomic co-localization
protein_associations_cnv_rna <- read_tsv("./files/signf_protein_associations_cnv_mrna_filtered_close_controlling_genomeWide.txt") %>%
  mutate(pair = paste(controlling, controlled, sep="_")) %>%
  dplyr::select(pair, everything())


# import protein list
protein_list_genome_wide <- read_tsv("./files/protein_list_control_status_genome_wide.txt")



# import E1 E2 E3 ubiquitination lists
e1_ubi <- read_tsv("./data/go/E1.ubi.GO.0004839.amigo2.txt", col_names = F)
e2_ubi <- read_tsv("./data/go/E2.ubi.GO.0061631.amigo2.txt", col_names = F)
e3_ubi <- read_tsv("./data/go/E3.ubi.GO.0061630.amigo2.txt", col_names = F)


length( intersect(e1_ubi$X2, protein_list_genome_wide %>% filter(control_status == "controlling") %>% pull(gene) ) )

length( intersect(e2_ubi$X2, protein_list_genome_wide %>% filter(control_status == "controlling") %>% pull(gene) ) )

length( intersect(e3_ubi$X2, protein_list_genome_wide %>% filter(control_status == "controlling") %>% pull(gene) ) )
