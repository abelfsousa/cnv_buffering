suppressMessages(library(tidyverse))


# import protein attenuation data
protein_attenuation <- read_tsv("./files/protein_attenuation.txt")


#new
protein_associations_cnv_rna_gw1 <- read_tsv("./files/signf_protein_associations_cnv_mrna_filtered_close_controlling_genomeWide.txt") %>%
  mutate(pair = paste(controlling, controlled, sep="_")) %>%
  dplyr::select(pair, everything())


protein_associations_cnv_rna_gw1_all <- read_tsv("./files/signf_protein_associations_cnv_mrna_genomeWide.txt") %>%
  mutate(pair = paste(controlling, controlled, sep="_")) %>%
  dplyr::select(pair, everything())


protein_associations_cnv <- readRDS("./files/protein_associations_tcga_ccle_cnv_reg_genome_wide.rds")
protein_associations_rna <- readRDS("./files/protein_associations_tcga_ccle_rna_reg_genome_wide.rds")



#old
protein_associations_cnv_rna_gw2 <- read_tsv("./files/older_genome_wide/signf_protein_associations_cnv_mrna_filtered_close_controlling_genomeWide.txt") %>%
  mutate(pair = paste(controlling, controlled, sep="_")) %>%
  dplyr::select(pair, everything())


protein_associations_cnv2 <- readRDS("./files/older_genome_wide/protein_associations_tcga_ccle_cnv_reg_genome_wide.rds")
protein_associations_rna2 <- readRDS("./files/older_genome_wide/protein_associations_tcga_ccle_rna_reg_genome_wide.rds")




specific_older <- protein_associations_cnv_rna_gw2 %>%
  filter(!pair %in% protein_associations_cnv_rna_gw1$pair)

protein_associations_cnv_rna_gw1_all %>% filter(pair %in% specific_older$pair)

protein_associations_cnv %>% inner_join(specific_older[, c("controlling", "controlled")], by=c("controlling", "controlled"))

specific_older2 <- specific_older[, c("controlling", "controlled")] %>% anti_join(protein_associations_cnv, by=c("controlling", "controlled")) %>%
  mutate(pair = paste(controlling, controlled, sep="_"))
