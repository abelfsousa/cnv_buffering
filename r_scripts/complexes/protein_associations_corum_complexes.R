# Abel Sousa 2018

# PhD project
# Study of Control of Protein Abundance using Multi-omics Data from Cancer Samples


# -- Number of protein associations found in corum complexes



suppressMessages(library(tidyverse))
options(bitmapType = "cairo")




# load significative protein associations with CNV and RNA
protein_associations_rna_cnv <- read_tsv("./files/signf_protein_associations_cnv_mrna_filtered_close_controlling.txt") %>%
  mutate(pair = paste(controlling, controlled, sep="_")) %>%
  dplyr::select(pair, everything())



# load corum data
corum <- read_tsv("./data/corum/coreComplexes_29_05_2018.txt")



# define all protein pairs by complex
corum_complexes <- corum %>%
	filter(Organism == "Human") %>%
	dplyr::select(ComplexName, `subunits(Gene name)`) %>%
	rename(GeneName = `subunits(Gene name)`) %>%
	group_by(ComplexName) %>%
	summarise(GeneName = paste(GeneName, collapse=";")) %>%
	ungroup() %>%
	mutate(GeneName = str_split(GeneName, ";\\s?")) %>%
	unnest() %>%
	filter(GeneName != "") %>%
	distinct() %>%
  group_by(ComplexName) %>%
	mutate(pair = list(tidyr::crossing(GeneName, GeneName))) %>%
	ungroup() %>%
	unnest() %>%
  dplyr::select(-GeneName) %>%
  distinct() %>%
	filter(GeneName1 != GeneName2) %>%
  mutate(pair = paste(GeneName1, GeneName2, sep="_")) %>%
  dplyr::select(-c(GeneName1, GeneName2))




#
complexes_associations <- corum_complexes %>%
  mutate(associations = pair %in% protein_associations_rna_cnv$pair) %>%
  group_by(ComplexName) %>%
  summarise(associations = sum(associations)) %>%
  arrange(desc(associations))
write.table(complexes_associations %>% filter(associations >= 4), "./files/corum_complexes_associations_4.txt", sep = "\t", quote = F, row.names=F)



complexes_associations_list <- corum_complexes %>%
  mutate(associations = pair %in% protein_associations_rna_cnv$pair) %>%
  filter(associations) %>%
  dplyr::select(-associations) %>%
  group_by(ComplexName) %>%
  mutate(associations = n()) %>%
  ungroup() %>%
  arrange(desc(associations), ComplexName, pair)
write.table(complexes_associations_list %>% filter(associations >= 4), "./files/corum_complexes_associations_4_list.txt", sep = "\t", quote = F, row.names=F)
