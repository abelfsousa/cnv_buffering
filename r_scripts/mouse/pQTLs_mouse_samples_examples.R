# Abel Sousa 2018

# PhD project
# Study of Control of Protein Abundance using Multi-omics Data from Cancer Samples

# -- Analysis of protein associations using normal samples from mouse

# -- Examples of pQTLs and eQTLs affecting protein and RNA expression in cis/trans





suppressMessages(library(tidyverse))
suppressMessages(library(ggpubr))
suppressMessages(library(gginnards))

options(bitmapType = "cairo")
set.seed(123)


load("./data/ChickMungeretal2016/ChickMungeretal2016_DiversityOutbred.Rdata")
load("./r_workspaces/protein_associations_QTLs_mouse_samples.RData")

chr15_genotypes <- read_tsv("./data/ChickMungeretal2016/genotypes_chr_15/Svenson-183_Svenson_DO-MegaMUGA-calls.txt") %>%
  gather(key="sample", value="genotype", -c(marker, chr, pos))


mice_pQTLs <- mice_pQTLs %>%
  mutate(pQTL_LocusID = as.numeric(str_split(pQTL_LocusID, "_", simplify=T)[,2])) %>%
  dplyr::select(-Gene_mouse)

#mice_pQTLs %>% filter(pQTL_Chromosome == 15) %>% inner_join(chr15_genotypes[, c("pos", "sample", "genotype")], by=c("pQTL_LocusID" = "pos")
#none pQTL position is mapping to the genotype positions



# load gene annotation
rna_IDs <- read_tsv("./data/ChickMungeretal2016/DO_mice_rna_IDs.txt")
protein_IDs <- read_tsv("./data/ChickMungeretal2016/DO_mice_protein_IDs.txt")



expr.protein.192 <- t(expr.protein.192)
expr.protein.192 <- cbind.data.frame(ID = rownames(expr.protein.192), expr.protein.192)

expr.protein.192 <- expr.protein.192 %>%
  as.tibble() %>%
  gather(key="sample", value="protein", -c(ID)) %>%
  inner_join(protein_IDs, by=c("ID" = "Protein_ID")) %>%
  dplyr::select(Gene_Symbol, sample, protein)


expr.rna.192 <- t(expr.rna.192)
expr.rna.192 <- cbind.data.frame(ID = rownames(expr.rna.192), expr.rna.192)

expr.rna.192 <- expr.rna.192 %>%
  as.tibble() %>%
  gather(key="sample", value="rna", -c(ID)) %>%
  inner_join(rna_IDs, by=c("ID" = "Gene_ID")) %>%
  dplyr::select(Gene_Symbol, sample, rna)



expr.rna.protein <- inner_join(expr.rna.192, expr.protein.192, by=c("Gene_Symbol", "sample"))


positions <- chr15_genotypes %>%
  filter(pos >= 82021830 & pos <= 82994447) %>%
  pull(pos) %>%
  unique()

for (i in positions){
print(i)
xrcc6_xrcc5 <- expr.rna.protein %>%
  filter(Gene_Symbol == "Xrcc6" | Gene_Symbol == "Xrcc5") %>%
  #inner_join(chr15_genotypes %>% filter(pos == "82375773") %>% dplyr::select(sample, genotype), by="sample") %>%
  inner_join(chr15_genotypes %>% filter(pos == i) %>% dplyr::select(sample, genotype), by="sample") %>%
  #filter(genotype != "--") %>%
  mutate(control_status = if_else(Gene_Symbol == "Xrcc5", "controlled", "controlling")) %>%
  gather(key="data", value="expression", -c(Gene_Symbol, sample, genotype, control_status)) %>%
  ggplot(mapping = aes(x = genotype, y = expression)) +
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(width = 0.1, alpha = 0.1) +
  #stat_compare_means(comparisons = list( c(1, 2), c(1, 3), c(2, 3) )) +
  facet_grid(control_status ~ data) +
  theme(axis.title.y=element_text(colour="black", size=15),
    axis.title.x=element_text(colour="black", size=15),
    axis.text.y=element_text(colour="black", size=12),
    axis.text.x=element_text(colour="black", size=12),
    strip.text = element_text(size=15)) +
  labs(x = "", y = "")
#xrcc6_xrcc5$layers[[which_layers(xrcc6_xrcc5, "GeomSignif")]]$aes_params$textsize <- 2
ggsave(filename=paste0("xrcc6_pQTL_xrcc5_expression_", i, ".png"), plot=xrcc6_xrcc5, width=5, height=5, path = "./plots/mouse_QTLs/")
unlink(paste0("xrcc6_pQTL_xrcc5_expression_", i, ".png"))
}
