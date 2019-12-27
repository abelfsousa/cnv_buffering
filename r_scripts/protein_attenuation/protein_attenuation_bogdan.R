suppressMessages(library(tidyverse))
suppressMessages(library(clusterProfiler))
suppressMessages(library(org.Hs.eg.db))
suppressMessages(library(data.table))

options(bitmapType = "cairo")


# map of uniprot/gene symbol/protein length
uniprot2gene <- read_tsv("./data/uniprot/uniprot_gene_plength_19_06_2018.tab") %>%
    dplyr::rename(UNIPROT=Entry, SYMBOL=`Gene names  (primary )`, LENGTH=Length) %>%
    dplyr::select(-Status, -Organism)


# corr(CNV, RNA) and corr(CNV, Protein)
protein_attenuation <- read_tsv("./files/protein_attenuation.txt") %>%
	inner_join(uniprot2gene[, c("UNIPROT", "SYMBOL")], by=c("gene" = "SYMBOL")) %>%
	dplyr::select(gene_id = gene, uniprot_id = UNIPROT, everything()) %>%
	group_by(gene_id) %>%
	filter(n() == 1) %>%
	ungroup()
write.table(protein_attenuation, "./files/protein_attenuation_bogdan.txt", sep="\t", quote=F, row.names=F)







# import not-assessed genes from Bogdan
not_assessed <- fread("./data/bogdan/gtexhigh_protein_attenuation_all.txt")
not_assessed <- not_assessed[, -c(1,7)] %>%
  as.tibble() %>%
  filter(attenuation == "not-assesed") %>%
  pull(gene_name) %>%
  unique()




# GO over-representation test
not_assessed_overRepGO <- enrichGO(
  not_assessed,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "ALL",
  pvalueCutoff = 0.1,
  qvalueCutoff = 0.1,
  pAdjustMethod = "BH",
  minGSSize = 10,
  maxGSSize = 500)
write.table(as.data.frame(not_assessed_overRepGO@result), "./files/not_assessed_over_representation_GO.txt", sep="\t", quote=F, row.names=F)

not_assessed_overRepGO_dot <- dotplot(not_assessed_overRepGO, x = "geneRatio", colorBy = "p.adjust", split = "ONTOLOGY", font.size = 15)
ggsave(filename="not_assessed_over_representation_GO.png", plot=not_assessed_overRepGO_dot, path="./plots/bogdan/", width = 10, height = 5)
unlink("not_assessed_over_representation_GO.png")
