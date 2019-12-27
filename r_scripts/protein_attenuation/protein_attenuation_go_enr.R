# Abel Sousa 2018

# PhD project
# Study of Control of Protein Abundance using Multi-omics Data from Cancer Samples



# -- GO enrichment of attenuated proteins


suppressMessages(library(tidyverse))
suppressMessages(library(clusterProfiler))
suppressMessages(library(org.Hs.eg.db))
suppressMessages(library(igraph))
suppressMessages(library(viridis))

options(bitmapType = "cairo")



# import protein attenuation data
protein_attenuation <- read_tsv("./files/protein_attenuation.txt")


non_attenuated <- protein_attenuation %>% filter(class == "non-attenuated") %>% pull(gene)
low_attenuated <- protein_attenuation %>% filter(class == "low-attenuated-protein") %>% pull(gene)
high_attenuated <- protein_attenuation %>% filter(class == "high-attenuated-protein") %>% pull(gene)
attenuated <- protein_attenuation %>% filter(class == "low-attenuated-protein" | class == "high-attenuated-protein") %>% pull(gene)



# GO over-representation test
go_test <- function(gene_set, universe, ontology, p.adj){
  enr <- enrichGO(
    gene = gene_set,
    universe = universe,
    ont = ontology,
    OrgDb = org.Hs.eg.db,
    keyType = "SYMBOL",
    pvalueCutoff = p.adj,
    qvalueCutoff = p.adj,
    pAdjustMethod = "BH",
    minGSSize = 5,
    maxGSSize = 500)
}



attenuated_GOenr <- go_test(attenuated, protein_attenuation$gene, "CC", 0.05)

attenuated_GOenr_dotplot <- dotplot(attenuated_GOenr, x = "geneRatio", colorBy = "p.adjust", showCategory = 20, font.size = 25, title = "")
ggsave(filename="attenuated_GOenr_dotplot.png", plot=attenuated_GOenr_dotplot, path="./plots/protein_attenuation_tcga_ccle/", width = 20, height = 15)
unlink("attenuated_GOenr_dotplot.png")


non_attenuated_GOenr <- go_test(non_attenuated, protein_attenuation$gene, "CC", 0.05)

non_attenuated_GOenr_dotplot <- dotplot(non_attenuated_GOenr, x = "geneRatio", colorBy = "p.adjust", showCategory = 20, font.size = 25, title = "")
ggsave(filename="non_attenuated_GOenr_dotplot.png", plot=non_attenuated_GOenr_dotplot, path="./plots/protein_attenuation_tcga_ccle/", width = 20, height = 15)
unlink("non_attenuated_GOenr_dotplot.png")





# load corum data
corum <- read_tsv("./data/corum/coreComplexes_29_05_2018.txt")


overlap_complexes <- function(complexes){
  x <- rep(1, length(complexes))

  for(i in 1:(length(complexes)-1)){
    for(j in (i+1):length(complexes)){
      if(x[i] == 1 & x[j] == 1){
        a <- complexes[[i]]
        b <- complexes[[j]]
        jac <- (length(intersect(a, b)))/(length(union(a, b)))
        if(jac >= 0.9){
          x[j] = 0
        }
      }
    }
  }
  return(x)
}


# protein subunits by corum complex
corum_complexes <- corum %>%
	filter(Organism == "Human") %>%
	dplyr::select(ComplexName, `subunits(Gene name)`) %>%
	dplyr::rename(GeneName = `subunits(Gene name)`) %>%
	group_by(ComplexName) %>%
	summarise(GeneName = paste(GeneName, collapse=";")) %>%
	ungroup() %>%
	mutate(GeneName = str_split(GeneName, ";\\s?")) %>%
	unnest() %>%
	filter(GeneName != "") %>%
	distinct() %>%
  group_by(ComplexName) %>%
  summarise(GeneName = list(GeneName)) %>%
  ungroup() %>%
  mutate(keep = overlap_complexes(GeneName)) %>%
  filter(keep == 1) %>%
  dplyr::select(-keep) %>%
  unnest(GeneName)


corum_test <- function(gene_set, universe, terms, p_adj, q_value){
  enr <- enricher(
    gene = gene_set,
    universe = universe,
    TERM2GENE = terms,
    pvalueCutoff = p_adj,
    qvalueCutoff = q_value,
    pAdjustMethod = "BH",
    minGSSize = 5,
    maxGSSize = 500)
}

attenuated_CORUMenr <- corum_test(attenuated, protein_attenuation$gene, corum_complexes, 1, 1)
attenuated_CORUMenr_tab <- attenuated_CORUMenr@result

non_attenuated_CORUMenr <- corum_test(non_attenuated, protein_attenuation$gene, corum_complexes, 1, 1)
non_attenuated_CORUMenr_tab <- non_attenuated_CORUMenr@result


corum_enr_tab <- attenuated_CORUMenr_tab %>%
  dplyr::select(ID, p.adjust_att=p.adjust) %>%
  inner_join(non_attenuated_CORUMenr_tab[, c("ID", "p.adjust")], by = "ID") %>%
  dplyr::rename(p.adjust_non=p.adjust) %>%
  gather(key = "p.adjust", value="value", -ID) %>%
  ggplot(mapping = aes(x = ID, y = value, color = p.adjust)) +
  geom_point() +
  geom_hline(yintercept=0.05, linetype="dashed", color = "black") +
  theme_classic() +
  theme(axis.title=element_text(colour="black", size=15),
        axis.text.y=element_text(colour="black", size=12),
        legend.title=element_text(colour="black", size=15),
        legend.text=element_text(colour="black", size=12),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  scale_x_discrete(name = "CORUM complexes") +
  scale_y_continuous(name = "Adjusted P-value") +
  scale_color_discrete(name = "Attenuation", labels = c("attenuated", "non-attenuated"))
ggsave(filename="corum_enrichment_attenuation.png", plot=corum_enr_tab, width=7, height=5, path = "./plots/protein_attenuation_tcga_ccle/")
unlink("corum_enrichment_attenuation.png")


attenuated_CORUMenr2 <- attenuated_CORUMenr
attenuated_CORUMenr2@result <- attenuated_CORUMenr2@result[attenuated_CORUMenr2@result$ID %in% (corum_enr_tab$data %>% filter(p.adjust == "p.adjust_att") %>% filter(value < 0.05) %>% pull(ID)), ]

attenuated_CORUMenr2_dotplot <- dotplot(attenuated_CORUMenr2, x = "geneRatio", colorBy = "p.adjust", showCategory = 23, font.size = 25, title = "")
ggsave(filename="corum_enrichment_attenuated_dotplot.png", plot=attenuated_CORUMenr2_dotplot, path="./plots/protein_attenuation_tcga_ccle/", width = 16, height = 12)
unlink("corum_enrichment_attenuated_dotplot.png")


attenuated_CORUMenr_barplot <- attenuated_CORUMenr_tab %>%
  as.tibble() %>%
  dplyr::select(ID, p_adj=p.adjust) %>%
  filter(p_adj < 0.05) %>%
  mutate(p_adj = -log10(p_adj)) %>%
  arrange(p_adj) %>%
  filter(!str_detect(ID, "ribosomal subunit, cytoplasmic")) %>%
  filter(!str_detect(ID, "ribosomal subunit, mitochondrial")) %>%
  filter(!str_detect(ID, "C complex spliceosome")) %>%
  filter(!str_detect(ID, "20S proteasome")) %>%
  filter(!str_detect(ID, "NuA4/Tip60 HAT complex")) %>%
  filter(!str_detect(ID, "ARC-L complex")) %>%
  filter(!ID == "STAGA complex") %>%
  mutate_if(is.character, as.factor) %>%
  mutate(ID = fct_inorder(ID)) %>%
  ggplot(mapping = aes(y = p_adj, x=ID, fill = p_adj)) +
    geom_bar(stat="identity") +
    theme_classic() +
    theme(
      axis.title.x=element_text(colour="black", size=22),
      axis.title.y=element_blank(),
      axis.text.y=element_text(colour="black", size=22),
      axis.text.x=element_text(colour="black", size=20),
      plot.title = element_blank(),
      legend.text=element_text(size=16),
      legend.title=element_text(size=18)) +
    coord_flip() +
    scale_fill_viridis(option="D", name="", guide=F) +
    scale_y_continuous(name = "Adjusted P-value (-log10)")
ggsave(filename="corum_enrichment_attenuated_barplot.png", plot=attenuated_CORUMenr_barplot, path="./plots/protein_attenuation_tcga_ccle/", width = 14, height = 9)
ggsave(filename="corum_enrichment_attenuated_barplot.pdf", plot=attenuated_CORUMenr_barplot, path="./plots/protein_attenuation_tcga_ccle/", width = 14, height = 9)
unlink("corum_enrichment_attenuated_barplot.png")
unlink("corum_enrichment_attenuated_barplot.pdf")


save(list=ls(), file="./r_workspaces/protein_attenuation_go_enr.RData")
