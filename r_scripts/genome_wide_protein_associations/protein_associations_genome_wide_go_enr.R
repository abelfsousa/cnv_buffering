# Abel Sousa 2018

# PhD project
# Study of Control of Protein Abundance using Multi-omics Data from Cancer Samples



# -- Enrichment on genome-wide protein associations



suppressMessages(library(tidyverse))
suppressMessages(library(clusterProfiler))
suppressMessages(library(org.Hs.eg.db))
suppressMessages(library(igraph))
suppressMessages(library(gridExtra))

options(bitmapType = "cairo")
set.seed(123)



# import protein attenuation data
protein_attenuation <- read_tsv("./files/protein_attenuation.txt")




# import physical protein associations with CNV and RNA
# filtered by genomic co-localization
protein_associations_cnv_rna_physical <- read_tsv("./files/signf_protein_associations_cnv_mrna_filtered_close_controlling.txt") %>%
  mutate(pair = paste(controlling, controlled, sep="_")) %>%
  dplyr::select(pair, everything())




# import genome-wide protein associations with CNV and RNA
# filtered by genomic co-localization
protein_associations_cnv_rna_gw <- read_tsv("./files/signf_protein_associations_cnv_mrna_filtered_close_controlling_genomeWide_fdr10.txt") %>%
  mutate(pair = paste(controlling, controlled, sep="_")) %>%
  dplyr::select(pair, everything()) %>%
  filter(!pair %in% protein_associations_cnv_rna_physical$pair) %>%
  #filter(controlled %in% (protein_attenuation %>% filter(class == "high-attenuated-protein") %>% pull(gene))) %>%
  filter(beta_cnv < 0, beta_rna < 0)




# number of times each protein is found and control status
protein_list <- protein_associations_cnv_rna_gw %>%
  dplyr::select(controlling, controlled) %>%
  dplyr::combine() %>%
  as_tibble() %>%
  dplyr::rename(gene = value) %>%
  group_by(gene) %>%
  dplyr::summarise(interactions = n()) %>%
  arrange(desc(interactions)) %>%
  dplyr::mutate(
    control_status = if_else(gene %in% setdiff(protein_associations_cnv_rna_gw$controlling, protein_associations_cnv_rna_gw$controlled), "controlling",
    if_else(gene %in% setdiff(protein_associations_cnv_rna_gw$controlled, protein_associations_cnv_rna_gw$controlling), "controlled", "both"))) %>%
  arrange(control_status)
write.table(protein_list, "./files/protein_list_control_status_genome_wide.txt", sep="\t", quote=F, row.names=F)



# import E1 E2 E3 ligases lists
e1_ubi <- read_tsv("./data/go/E1.ubi.GO.0004839.amigo2.txt", col_names = F) %>%
  dplyr::select(gene=X2) %>%
  mutate(ligase = "E1", term="ligase")

e2_ubi <- read_tsv("./data/go/E2.ubi.GO.0061631.amigo2.txt", col_names = F) %>%
  dplyr::select(gene=X2) %>%
  mutate(ligase = "E2", term="ligase") %>%
  filter(!gene %in% c("ube2n-ube2v2_human", "ube2n-ube2v1_human", "TMEM189-UBE2V1"))

e3_ubi <- read_tsv("./data/go/E3.ubi.GO.0061630.amigo2.txt", col_names = F) %>%
  dplyr::select(gene=X2) %>%
  mutate(ligase = "E3", term="ligase") %>%
  filter(gene != "brcc_human")

ligases <- bind_rows(e1_ubi, e2_ubi, e3_ubi)



# import all go terms
all_go_terms <- readRDS("./data/go/all_go_terms_symbols.rds")
all_go_terms <- tibble(term = names(all_go_terms), gene = all_go_terms) %>%
  unnest(gene)

# import KEGG pathways
kegg_pathways <- readRDS("./data/kegg/kegg_pathways_symbols.rds")
kegg_pathways <- tibble(term = names(kegg_pathways), gene = kegg_pathways) %>%
  unnest(gene)

degradation_terms <- "LIGASE|UBIQUITIN|AUTOPHAGY|ERAD|PROTEASOME|PROTEIN_TRANSPORT"

degradation_go_terms <- all_go_terms  %>%
  filter(str_detect(term, degradation_terms))


degradation_kegg_pathways <- kegg_pathways  %>%
  filter(str_detect(term, degradation_terms))


terms <- bind_rows(ligases[, c("term", "gene")], degradation_go_terms, degradation_kegg_pathways, protein_attenuation %>% dplyr::select(gene) %>% mutate(term = "ALL"))



# ENRICHMENT ANALYSIS


#controlling proteins
protein_set_controlling <- protein_list %>% filter(control_status == "controlling") %>% pull(gene)
universe_set <- protein_attenuation %>% filter(class == "non-attenuated") %>% pull(gene)
#universe_set <- protein_attenuation %>% pull(gene)


controlling_enr_degradation <- enricher(
  gene=protein_set_controlling,
  universe = universe_set,
  TERM2GENE = terms,
  pvalueCutoff = 1,
  qvalueCutoff = 1,
  pAdjustMethod = "BH",
  minGSSize = 10,
  maxGSSize = 10000)
write.table(
  as.data.frame(controlling_enr_degradation@result %>% dplyr::select(ID, pvalue, Count) %>% arrange(pvalue) %>% head(10)),
  "./files/degradation_terms_all_controlled.txt", sep="\t", quote=F, row.names=F)
controlling_enr_degradation@result <- controlling_enr_degradation@result %>% arrange(pvalue) %>% head(10)
controlling_enr_degradation_dot <- dotplot(controlling_enr_degradation, x = "geneRatio", colorBy = "pvalue", font.size = 10)
ggsave(filename="controlling_enrichment_degradation.png", plot=controlling_enr_degradation_dot, path="./plots/protein_associations_genome_wide/", width = 10, height = 4)
unlink("controlling_enrichment_degradation.png")



# GO over-representation test
protein_set_controlling_overRepGO <- enrichGO(
  protein_set_controlling,
  OrgDb = org.Hs.eg.db,
  universe = universe_set,
  keyType = "SYMBOL",
  ont = "ALL",
  pvalueCutoff = 1,
  qvalueCutoff = 1,
  pAdjustMethod = "BH",
  minGSSize = 10,
  maxGSSize = 500)
protein_set_controlling_overRepGO@result <- protein_set_controlling_overRepGO@result %>% arrange(p.adjust) %>% head
write.table(protein_set_controlling_overRepGO@result %>% dplyr::select(Description, ONTOLOGY, pvalue, p.adjust, Count, geneID),
  "./files/protein_set_controlling_over_representation_GO.txt", sep="\t", quote=F, row.names=F)
protein_set_controlling_overRepGO_dot <- dotplot(protein_set_controlling_overRepGO, x = "geneRatio", colorBy = "p.adjust", font.size = 15)
ggsave(filename="protein_set_controlling_over_representation_GO.png", plot=protein_set_controlling_overRepGO_dot, path="./plots/protein_associations_genome_wide/", width = 10, height = 4)
unlink("protein_set_controlling_over_representation_GO.png")


protein_set_controlling_overRepKEGG <- enrichKEGG(
  bitr(protein_set_controlling, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db", drop=TRUE)[,c(2)],
  organism = "hsa",
  keyType = "ncbi-geneid",
  pvalueCutoff = 0.2,
  qvalueCutoff = 0.2,
  pAdjustMethod = "BH",
  minGSSize = 10,
  maxGSSize = 500,
  use_internal_data = FALSE)



#controlled proteins
protein_set_controlled <- protein_list %>% filter(control_status == "controlled") %>% pull(gene)

protein_set_controlled_overRepGO <- enrichGO(
  protein_set_controlled,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "ALL",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  minGSSize = 10,
  maxGSSize = 500,
  universe = universe_set)
#write.table(as.data.frame(protein_set_controlled_overRepGO@result), "./files/protein_set_controlled_over_representation_GO.txt", sep="\t", quote=F, row.names=F)

protein_set_controlled_overRepGO_dot <- dotplot(protein_set_controlled_overRepGO, split="ONTOLOGY", x = "geneRatio", colorBy = "p.adjust", font.size = 10)
ggsave(filename="protein_set_controlled_over_representation_GO.png", plot=protein_set_controlled_overRepGO_dot, path="./plots/protein_associations_genome_wide/", width = 10, height = 15)
unlink("protein_set_controlled_over_representation_GO.png")

protein_set_controlled_overRepKEGG <- enrichKEGG(
  bitr(protein_set_controlled, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db", drop=TRUE)[,c(2)],
  organism = "hsa",
  keyType = "ncbi-geneid",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  minGSSize = 10,
  maxGSSize = 500,
  use_internal_data = FALSE)
