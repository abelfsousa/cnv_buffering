# Abel Sousa 2018

# PhD project
# Study of Control of Protein Abundance using Multi-omics Data from Cancer Samples



# -- Test essentiality of protein association pairs found with CNV across cell lines




suppressMessages(library(tidyverse))
#suppressMessages(library(CePa))
suppressMessages(library(ggpubr))
suppressMessages(library(gplots))

options(bitmapType = "cairo")




# significative protein associations with CNV and RNA
signf_protein_associations_cnv_rna <- read_tsv("./files/signf_protein_associations_cnv_mrna_filtered_close_controlling.txt")





# load Achilles CRISPR Ceres data
achilles_CRISPR <- read_csv("/nfs/research1/beltrao/memon/databases/cell_lines/crispr/ceres-gene-effects.csv") %>%
	dplyr::rename(gene=X1) %>%
	gather(key="cell_line", value="essentiality_score", -gene)
#17670 genes

achilles_CRISPR_sd <- achilles_CRISPR %>% pull(essentiality_score) %>% sd(na.rm=T)

achilles_CRISPR <- achilles_CRISPR %>%
	group_by(gene) %>%
	filter(sum(essentiality_score < -achilles_CRISPR_sd, na.rm = TRUE) > n()*0.05) %>%
	ungroup()
achilles_CRISPR$gene %>% unique %>% length
#5532 genes




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




# corum protein complexes from protein pairs with gene essentiality scores
ppairs_crispr_complexes <- signf_protein_associations_cnv_rna %>%
  mutate(pair = paste(controlling, controlled, sep="_")) %>%
  dplyr::select(pair, controlling, controlled) %>%
	filter(controlling %in% achilles_CRISPR$gene & controlled %in% achilles_CRISPR$gene) %>%
  inner_join(corum_complexes, by = "pair") %>%
  dplyr::select(ComplexName, everything()) %>%
  group_by(ComplexName) %>%
  summarise(counts = n()) %>%
  arrange(desc(counts))





# essentiality scores for genes from protein pairs
ppairs_crispr <- signf_protein_associations_cnv_rna %>%
  mutate(pair = paste(controlling, controlled, sep="_")) %>%
  dplyr::select(pair, controlling, controlled) %>%
	filter(controlling %in% achilles_CRISPR$gene & controlled %in% achilles_CRISPR$gene) %>%
  inner_join(corum_complexes, by = "pair") %>%
  dplyr::select(ComplexName, everything()) %>%
  group_by(ComplexName) %>%
  mutate(counts = n()) %>%
  arrange(desc(counts), ComplexName) %>%
  mutate(gene = list(unique(c(controlled, controlling)))) %>%
  ungroup() %>%
  dplyr::select(-c(pair, controlling, controlled, counts)) %>%
  unnest() %>%
  distinct() %>%
  inner_join(achilles_CRISPR, by="gene") %>%
  spread(key="cell_line", value = "essentiality_score")




pdf(file="./plots/complexes/complexes_essentiality_crispr_row.pdf", w=10, h=4)
#par(oma=c(2,0,0,3))


for(complex in ppairs_crispr_complexes %>% filter(counts >= 4) %>% pull(ComplexName)){

example <- ppairs_crispr %>%
  filter(ComplexName == complex) %>%
  dplyr::select(-c(ComplexName)) %>%
  as.data.frame()

rownames(example) <- example$gene
example <- example[, -c(1)]


heatmap.2( x=as.matrix(example),
  Rowv=T,
	Colv=T,
	dendrogram="both",
	scale = "row",
	trace = "none",
	#col = greenred,
	col = colorpanel(50,"green","black","red"),
	main = complex,
	xlab = "",
	ylab = "",
  labCol = NA,
	key.title=NA,
	cexCol = 1,
	cexRow = 0.6)

}

dev.off()






save(list=ls(), file="./r_workspaces/complexes_essentiality.RData")
