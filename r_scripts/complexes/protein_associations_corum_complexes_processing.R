# Abel Sousa 2018

# PhD project
# Study of Control of Protein Abundance using Multi-omics Data from Cancer Samples


# -- Description of protein subunits inside of CORUM complexes


suppressMessages(library(tidyverse))
suppressMessages(library(ggpubr))
suppressMessages(library(rjson))
options(bitmapType = "cairo")


# export FREESASA=/homes/abelsousa/freesasa-master/src/freesasa


# number of protein interactions and control status
protein_list <- read_tsv("./files/protein_list_control_status.txt")


# import protein attenuation state
protein_attenuation <- read_tsv("./files/protein_attenuation.txt") %>%
  mutate_at(vars(contains("class")), as.factor) %>%
  mutate_if(is.factor, fct_infreq)


# map of uniprot/gene symbol/protein length
uniprot2gene <- read_tsv("./data/uniprot/uniprot_gene_plength_19_06_2018.tab") %>%
    dplyr::rename(UNIPROT=Entry, SYMBOL=`Gene names  (primary )`, LENGTH=Length) %>%
    dplyr::select(-Status, -Organism)



get_areas <- function(complex, json_file){
  sasa <- fromJSON(file = json_file)

  chains <- sapply(sasa$results[[1]]$structure[[1]]$chains, function(x) x$label)
  area_incomplex <- sapply(sasa$results[[1]]$structure[[1]]$chains, function(x) x$area$total)
  area_outcomplex <- sapply(sasa$results, function(x) x$structure[[1]]$chains[[1]]$area$total)
  area_outcomplex <- area_outcomplex[-c(1)]

  res <- tibble(complex, chains, area_incomplex, area_outcomplex, perc_area_incomplex = ((area_outcomplex-area_incomplex)/area_outcomplex)*100)

  return(res)
}



# COP9_signalosome
# $FREESASA --format=json --chain-groups=A+B+C+D+E+F+G+H 4d10_half_dimer.pdb > 4d10_half_dimer_acc_areas.json

# load chain gene map
cop9_chain_gene <- read_tsv("./data/pdbs/COP9_signalosome/freesasa/4d10_chain_gene.txt")

cop9 <- get_areas("COP9_signalosome", "./data/pdbs/COP9_signalosome/freesasa/4d10_half_dimer_acc_areas.json") %>%
  inner_join(cop9_chain_gene, by=c("chains" = "chain")) %>%
  inner_join(uniprot2gene[, c("SYMBOL", "LENGTH")], by=c("gene" = "SYMBOL")) %>%
  inner_join(protein_attenuation[, c("gene", "attenuation", "class")], by="gene") %>%
  left_join(protein_list[, c("gene", "control_status")], by="gene") %>%
  dplyr::select(complex, chains, gene, length=LENGTH, attenuation, class, control_status, area_incomplex, area_outcomplex, perc_area_incomplex)



# Spliceosome (total) with RNA molecules
# $FREESASA --format=json --chain-groups=A+B+C+D+E+F+G+H+I+J+K+L+M+N+O+P+Q+R+S+T+U+V+W+X+Y+Z+a+b+c+d+f+e+g+h+i+j+k+l+m+n+o+p+q+r+s+t+u+v+w+x 5xjc.pdb > 5xjc_acc_area.json

# load chain gene map
spliceosome_chain_gene <- read_tsv("./data/pdbs/spliceosome/freesasa/5xjc_chain_gene.txt")

spliceosome1 <- get_areas("Spliceosome", "./data/pdbs/spliceosome/freesasa/5xjc_acc_area.json") %>%
  inner_join(spliceosome_chain_gene, by=c("chains" = "chain")) %>%
  filter(!is.na(gene)) %>%
  inner_join(uniprot2gene[, c("SYMBOL", "LENGTH")], by=c("gene" = "SYMBOL")) %>%
  inner_join(protein_attenuation[, c("gene", "attenuation", "class")], by="gene") %>%
  left_join(protein_list[, c("gene", "control_status")], by="gene") %>%
  dplyr::select(complex, chains, gene, length=LENGTH, attenuation, class, control_status, area_incomplex, area_outcomplex, perc_area_incomplex) %>%
  group_by(gene) %>%
  summarise(complex = unique(complex), chains = paste(chains, collapse=";"), length = unique(length), attenuation = unique(attenuation), class = unique(class), control_status = unique(control_status), area_incomplex = mean(area_incomplex), area_outcomplex = mean(area_outcomplex), perc_area_incomplex = mean(perc_area_incomplex)) %>%
  dplyr::select(complex, chains, gene, everything())


# Spliceosome (total) without RNA molecules
# $FREESASA --format=json --chain-groups=A+C+D+E+I+J+K+L+M+N+O+P+Q+R+S+T+U+V+W+X+Y+Z+a+b+c+d+f+e+g+h+i+j+k+l+m+n+o+p+q+r+s+t+u+v+w+x 5xjc_just_proteins.pdb > 5xjc_just_proteins_acc_area.json
spliceosome2 <- get_areas("Spliceosome", "./data/pdbs/spliceosome/freesasa/5xjc_just_proteins_acc_area.json") %>%
  inner_join(spliceosome_chain_gene, by=c("chains" = "chain")) %>%
  inner_join(uniprot2gene[, c("SYMBOL", "LENGTH")], by=c("gene" = "SYMBOL")) %>%
  inner_join(protein_attenuation[, c("gene", "attenuation", "class")], by="gene") %>%
  left_join(protein_list[, c("gene", "control_status")], by="gene") %>%
  dplyr::select(complex, chains, gene, length=LENGTH, attenuation, class, control_status, area_incomplex, area_outcomplex, perc_area_incomplex) %>%
  group_by(gene) %>%
  summarise(complex = unique(complex), chains = paste(chains, collapse=";"), length = unique(length), attenuation = unique(attenuation), class = unique(class), control_status = unique(control_status), area_incomplex = mean(area_incomplex), area_outcomplex = mean(area_outcomplex), perc_area_incomplex = mean(perc_area_incomplex)) %>%
  dplyr::select(complex, chains, gene, everything())



# RNA polymerase II
# map orthologous genes between bos taurus (cow) and homo sapiens

# load human cow orthologous
ortho_human_cow <-  read_tsv("./data/ensembl/mart_export-human-cow-orthologous.txt") %>%
  filter(Cow_homology_type == "ortholog_one2one")

# load chain gene map
RNApolII_chain_gene <- read_tsv("./data/pdbs/RNA_polymerase_II/freesasa/5flm_chain_gene.txt")

# without DNA/RNA molecules
# $FREESASA --format=json --chain-groups=A+B+C+D+E+F+G+H+I+J+K+L 5flm_just_proteins.pdb > 5flm_just_proteins_acc_area.json
RNApolII1 <- get_areas("RNA_polymerase_II", "./data/pdbs/RNA_polymerase_II/freesasa/5flm_just_proteins_acc_area.json") %>%
  inner_join(RNApolII_chain_gene, by=c("chains" = "chain")) %>%
  inner_join(ortho_human_cow[, c("Gene_name", "Cow_gene_name")], by=c("gene" = "Cow_gene_name")) %>%
  dplyr::select(-gene) %>%
  dplyr::rename(gene=Gene_name) %>%
  inner_join(uniprot2gene[, c("SYMBOL", "LENGTH")], by=c("gene" = "SYMBOL")) %>%
  inner_join(protein_attenuation[, c("gene", "attenuation", "class")], by="gene") %>%
  left_join(protein_list[, c("gene", "control_status")], by="gene") %>%
  dplyr::select(complex, chains, gene, length=LENGTH, attenuation, class, control_status, area_incomplex, area_outcomplex, perc_area_incomplex)


# with DNA/RNA molecules
# $FREESASA --format=json --chain-groups=A+B+C+D+E+F+G+H+I+J+K+L+N+P+T 5flm.pdb > 5flm_acc_area.json
RNApolII2 <- get_areas("RNA_polymerase_II", "./data/pdbs/RNA_polymerase_II/freesasa/5flm_acc_area.json") %>%
  inner_join(RNApolII_chain_gene, by=c("chains" = "chain")) %>%
  filter(!is.na(gene)) %>%
  inner_join(ortho_human_cow[, c("Gene_name", "Cow_gene_name")], by=c("gene" = "Cow_gene_name")) %>%
  dplyr::select(-gene) %>%
  dplyr::rename(gene=Gene_name) %>%
  inner_join(uniprot2gene[, c("SYMBOL", "LENGTH")], by=c("gene" = "SYMBOL")) %>%
  inner_join(protein_attenuation[, c("gene", "attenuation", "class")], by="gene") %>%
  left_join(protein_list[, c("gene", "control_status")], by="gene") %>%
  dplyr::select(complex, chains, gene, length=LENGTH, attenuation, class, control_status, area_incomplex, area_outcomplex, perc_area_incomplex)



# eIF3
# map orthologous genes between oryctolagus cuniculus (rabbit) and homo sapiens

# load human rabbit orthologous
ortho_human_rabbit <-  read_tsv("./data/ensembl/mart_export-human-rabbit-orthologous.txt") %>%
  filter(Rabbit_homology_type == "ortholog_one2one")

# load chain gene map
eif3_chain_gene <- read_tsv("./data/pdbs/eIF3/freesasa/5a5t_chain_gene.txt")

# $FREESASA --format=json --chain-groups=A+C+E+F+H+K+L+M 5a5t.pdb > 5a5t_acc_areas.json
eif3 <- get_areas("eIF3", "./data/pdbs/eIF3/freesasa/5a5t_acc_areas.json") %>%
  inner_join(eif3_chain_gene, by=c("chains" = "chain")) %>%
  inner_join(ortho_human_rabbit[, c("Gene_name", "Rabbit_gene_name")], by=c("gene" = "Rabbit_gene_name")) %>%
  dplyr::select(-gene) %>%
  dplyr::rename(gene=Gene_name) %>%
  inner_join(uniprot2gene[, c("SYMBOL", "LENGTH")], by=c("gene" = "SYMBOL")) %>%
  inner_join(protein_attenuation[, c("gene", "attenuation", "class")], by="gene") %>%
  left_join(protein_list[, c("gene", "control_status")], by="gene") %>%
  dplyr::select(complex, chains, gene, length=LENGTH, attenuation, class, control_status, area_incomplex, area_outcomplex, perc_area_incomplex)




# V-ATPase
# map orthologous genes between s cerevisiae and homo sapiens

# load human s cerevisiae orthologous
ortho_human_scerevisiae <-  read_tsv("./data/ensembl/mart_export-human-s_cerevisiae-orthologous.txt") %>%
  filter(Saccharomyces_cerevisiae_homology_type == "ortholog_one2one")

# load chain gene map
VATPase_chain_gene <- read_tsv("./data/pdbs/V_ATPase/freesasa/3j9t_chain_gene.txt")


# $FREESASA --format=json --chain-groups=M+N+A+B+C+D+E+F+Q+L+K+P+b+O+H+G+J+I+Y+R+U+V+T+W+S+X+Z+a 3j9t.pdb > 3j9t_acc_areas.json
VATPase <- get_areas("V-ATPase", "./data/pdbs/V_ATPase/freesasa/3j9t_acc_areas.json") %>%
  inner_join(VATPase_chain_gene, by=c("chains" = "chain")) %>%
  inner_join(ortho_human_scerevisiae[, c("Gene_name", "Saccharomyces_cerevisiae_gene_name")], by=c("gene" = "Saccharomyces_cerevisiae_gene_name")) %>%
  dplyr::select(-gene) %>%
  dplyr::rename(gene=Gene_name) %>%
  inner_join(uniprot2gene[, c("SYMBOL", "LENGTH")], by=c("gene" = "SYMBOL")) %>%
  inner_join(protein_attenuation[, c("gene", "attenuation", "class")], by="gene") %>%
  left_join(protein_list[, c("gene", "control_status")], by="gene") %>%
  dplyr::select(complex, chains, gene, length=LENGTH, attenuation, class, control_status, area_incomplex, area_outcomplex, perc_area_incomplex) %>%
  group_by(gene) %>%
  summarise(complex = unique(complex), chains = paste(chains, collapse=";"), length = unique(length), attenuation = unique(attenuation), class = unique(class), control_status = unique(control_status), area_incomplex = mean(area_incomplex), area_outcomplex = mean(area_outcomplex), perc_area_incomplex = mean(perc_area_incomplex)) %>%
  dplyr::select(complex, chains, gene, everything())




# CCT complex
# map orthologous genes between s cerevisiae and homo sapiens

# load chain gene map
CCT_chain_gene <- read_tsv("./data/pdbs/CCT/freesasa/4v94_half_dimer_chain_gene.txt")


# $FREESASA --format=json --chain-groups=A+B+C+D+E+F+G+H+I+J+K+L+M+N+O+P 4v94_half_dimer.pdb > 4v94_half_dimer_acc_area.json
CCT <- get_areas("CCT", "./data/pdbs/CCT/freesasa/4v94_half_dimer_acc_area.json") %>%
  inner_join(CCT_chain_gene, by=c("chains" = "chain")) %>%
  inner_join(ortho_human_scerevisiae[, c("Gene_name", "Saccharomyces_cerevisiae_gene_name")], by=c("gene" = "Saccharomyces_cerevisiae_gene_name")) %>%
  dplyr::select(-gene) %>%
  dplyr::rename(gene=Gene_name) %>%
  inner_join(uniprot2gene[, c("SYMBOL", "LENGTH")], by=c("gene" = "SYMBOL")) %>%
  inner_join(protein_attenuation[, c("gene", "attenuation", "class")], by="gene") %>%
  left_join(protein_list[, c("gene", "control_status")], by="gene") %>%
  dplyr::select(complex, chains, gene, length=LENGTH, attenuation, class, control_status, area_incomplex, area_outcomplex, perc_area_incomplex) %>%
  group_by(gene) %>%
  summarise(complex = unique(complex), chains = paste(chains, collapse=";"), length = unique(length), attenuation = unique(attenuation), class = unique(class), control_status = unique(control_status), area_incomplex = mean(area_incomplex), area_outcomplex = mean(area_outcomplex), perc_area_incomplex = mean(perc_area_incomplex)) %>%
  dplyr::select(complex, chains, gene, everything())



# merge complexes
complexes <- bind_rows(cop9, spliceosome2, RNApolII1, eif3, VATPase, CCT)
write.table(complexes, "./files/complexes_sasa.txt", sep="\t", row.names = F, quote=F)




# import protein associations with CNV and RNA
# filtered by genomic co-localization
protein_associations_cnv_mrna <- read_tsv("./files/signf_protein_associations_cnv_mrna_filtered_close_controlling.txt") %>%
  mutate(pair = paste(controlling, controlled, sep="_")) %>%
  dplyr::select(pair, everything())



#controlled proteins
controlled_sasa <- complexes %>%
  filter(control_status == "controlled") %>%
  inner_join(protein_associations_cnv_mrna %>% dplyr::select(controlled, beta_cnv, fdr_cnv), by=c("gene"="controlled")) %>%
  mutate(beta_cnv = abs(beta_cnv), fdr_cnv = -log10(fdr_cnv)) %>%
  group_by(complex, chains, gene, area_incomplex, area_outcomplex, perc_area_incomplex) %>%
  summarise(beta_cnv = median(beta_cnv), fdr_cnv = median(fdr_cnv)) %>%
  gather(key="stat", value="value", -c("gene", "complex", "chains", "area_incomplex", "area_outcomplex", "perc_area_incomplex")) %>%
  gather(key="area_type", value="area_value", -c("gene", "complex", "chains", "stat", "value")) %>%
  ggplot(mapping = aes(x = value, y = area_value)) +
    geom_point() +
    geom_smooth(method=lm, color="blue", se=T) +
    facet_grid(area_type ~ stat, scales = "free") +
    stat_cor() +
    theme_classic() +
    theme(axis.title.y=element_text(colour="black", size=15),
          axis.title.x=element_text(colour="black", size=15),
          axis.text=element_text(colour="black", size=12),
          strip.text = element_text(size=12),
          strip.background = element_blank(),
          plot.title=element_text(size=12)) +
    labs(x = "Model statistics", y = "SASA", title="")
ggsave(filename="controlled_cor_sasa_stats.png", plot=controlled_sasa, path = "./plots/complexes/")
unlink("controlled_cor_sasa_stats.png")


#controlling proteins
controlling_sasa <- complexes %>%
  filter(control_status == "controlling") %>%
  inner_join(protein_associations_cnv_mrna %>% dplyr::select(controlling, beta_cnv, fdr_cnv), by=c("gene"="controlling")) %>%
  mutate(beta_cnv = abs(beta_cnv), fdr_cnv = -log10(fdr_cnv)) %>%
  group_by(complex, chains, gene, area_incomplex, area_outcomplex, perc_area_incomplex) %>%
  summarise(beta_cnv = median(beta_cnv), fdr_cnv = median(fdr_cnv)) %>%
  gather(key="stat", value="value", -c("gene", "complex", "chains", "area_incomplex", "area_outcomplex", "perc_area_incomplex")) %>%
  gather(key="area_type", value="area_value", -c("gene", "complex", "chains", "stat", "value")) %>%
  ggplot(mapping = aes(x = value, y = area_value)) +
    geom_point() +
    geom_smooth(method=lm, color="blue", se=T) +
    facet_grid(area_type ~ stat, scales = "free") +
    stat_cor() +
    theme_classic() +
    theme(axis.title.y=element_text(colour="black", size=15),
          axis.title.x=element_text(colour="black", size=15),
          axis.text=element_text(colour="black", size=12),
          strip.text = element_text(size=12),
          strip.background = element_blank(),
          plot.title=element_text(size=12)) +
    labs(x = "Model statistics", y = "SASA", title="")
ggsave(filename="controlling_cor_sasa_stats.png", plot=controlling_sasa, path = "./plots/complexes/")
unlink("controlling_cor_sasa_stats.png")


# control status by percentage of area inside complex
control_status_sasa <- ggplot(data = complexes %>% filter(!is.na(control_status)), mapping = aes(x = control_status, y = perc_area_incomplex, color=control_status)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.4, colour="grey") +
  theme_classic() +
  stat_compare_means(comparisons = list( c("both", "controlled"), c("controlled", "controlling"), c("both", "controlling") )) +
  theme(axis.title=element_text(colour="black", size=15),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    axis.text.y=element_text(colour="black", size=12),
    plot.title=element_text(size=15)) +
  scale_color_discrete(name = "") +
  labs(x = "Control status", y = "Percentage of are inside complex", title="")
ggsave(filename="complexes_control_status_sasa.png", plot = control_status_sasa, path = "./plots/complexes/")
unlink("complexes_control_status_sasa.png")



#attenuation
attenuation_sasa1 <- ggplot(data = complexes, mapping = aes(x = attenuation, y = perc_area_incomplex)) +
    geom_point() +
    geom_smooth(method=lm , color="blue", se=T) +
    stat_cor() +
    theme_classic() +
    facet_wrap( ~ complex, scales = "free") +
    theme(axis.title=element_text(colour="black", size=15),
        axis.text=element_text(colour="black", size=12),
        strip.text.x = element_text(size=15),
        strip.background = element_blank(),
        plot.title=element_text(size=15)) +
    labs(x = "Protein attenuation", y = "Percentage of are inside complex", title="")
ggsave(filename="attenuation_cor_sasa1.png", plot = attenuation_sasa1, path = "./plots/complexes/")
unlink("attenuation_cor_sasa1.png")


attenuation_sasa2 <- complexes %>%
  dplyr::select(-c(length, class, control_status)) %>%
  gather(key="area_type", value="area_value", -c("gene", "complex", "chains", "attenuation")) %>%
  ggplot(mapping = aes(x = attenuation, y = area_value)) +
    geom_point() +
    geom_smooth(method=lm, color="blue", se=T) +
    facet_grid(area_type ~ complex, scales = "free") +
    stat_cor() +
    theme_classic() +
    theme(axis.title.y=element_text(colour="black", size=15),
          axis.title.x=element_text(colour="black", size=15),
          axis.text=element_text(colour="black", size=12),
          strip.text = element_text(size=12),
          strip.background = element_blank(),
          plot.title=element_text(size=12)) +
    labs(x = "Protein attenuation", y = "SASA", title="")
ggsave(filename="attenuation_cor_sasa2.png", plot=attenuation_sasa2, height=10, width=12, path = "./plots/complexes/")
unlink("attenuation_cor_sasa2.png")



attenuation_sasa3 <- ggplot(data = complexes, mapping = aes(x = attenuation, y = perc_area_incomplex)) +
    geom_point(aes(color = complex)) +
    geom_smooth(method=lm , color="blue", se=T) +
    stat_cor() +
    theme_classic() +
    theme(axis.title=element_text(colour="black", size=15),
        axis.text=element_text(colour="black", size=12),
        plot.title=element_text(size=15)) +
    labs(x = "Protein attenuation", y = "Percentage of are inside complex", title="")
ggsave(filename="attenuation_cor_sasa3.png", plot = attenuation_sasa3, path = "./plots/complexes/")
unlink("attenuation_cor_sasa3.png")


attenuation_sasa4 <- ggplot(data = complexes, mapping = aes(x = attenuation, y = area_incomplex)) +
    geom_point(aes(color = complex)) +
    geom_smooth(method=lm , color="blue", se=T) +
    stat_cor() +
    theme_classic() +
    theme(axis.title=element_text(colour="black", size=15),
        axis.text=element_text(colour="black", size=12),
        plot.title=element_text(size=15)) +
    labs(x = "Protein attenuation", y = "Solvent accessible area inside complex", title="")
ggsave(filename="attenuation_cor_sasa4.png", plot = attenuation_sasa4, path = "./plots/complexes/")
unlink("attenuation_cor_sasa4.png")


attenuation_sasa5 <- ggplot(data = complexes, mapping = aes(x = attenuation, y = area_outcomplex)) +
    geom_point(aes(color = complex)) +
    geom_smooth(method=lm , color="blue", se=T) +
    stat_cor() +
    theme_classic() +
    theme(axis.title=element_text(colour="black", size=15),
        axis.text=element_text(colour="black", size=12),
        plot.title=element_text(size=15)) +
    labs(x = "Protein attenuation", y = "Solvent accessible area outside complex", title="")
ggsave(filename="attenuation_cor_sasa5.png", plot = attenuation_sasa5, path = "./plots/complexes/")
unlink("attenuation_cor_sasa5.png")



attenuation_sasa6 <- ggplot(data = complexes, mapping = aes(x = class, y = perc_area_incomplex, color=class)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.4, colour="grey") +
  theme_classic() +
  theme(axis.title=element_text(colour="black", size=15),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    axis.text.y=element_text(colour="black", size=12),
    plot.title=element_text(size=15)) +
  scale_color_discrete(name = "") +
  labs(x = "Attenuation state", y = "Percentage of are inside complex", title="")
ggsave(filename="complexes_attenuation_sasa.png", plot = attenuation_sasa6, path = "./plots/complexes/")
unlink("complexes_attenuation_sasa.png")



cop9_sasa <- ggplot(data = complexes %>% filter(complex == "COP9_signalosome"), mapping = aes(x = attenuation, y = perc_area_incomplex)) +
    geom_point(aes(color=class)) +
    geom_smooth(method=lm , color="black", fill = "grey", se=T, alpha=0.2) +
    stat_cor(size=3, label.y=50) +
    #annotate("text", label="r = 0.47", x = 0.3, y = 45, size=3) +
    theme_classic() +
    theme(axis.title=element_text(colour="black", size=10),
        axis.text=element_text(colour="black", size=8),
        legend.title=element_text(colour="black", size=10),
        legend.text=element_text(colour="black", size=8),
        plot.title=element_blank(),
        aspect.ratio = 1,
        legend.position="bottom") +
    scale_color_manual(values=c("#9ecae1", "#3182bd"), labels = c("lowly-attenuated", "highly-attenuated"), name = "Attenuation") +
    guides(color = guide_legend(title.position = "top")) +
    labs(x = "Protein attenuation", y = "Area inside complex (%)", title="")
ggsave(filename="cor_attenuation_sasa_cop9.png", plot = cop9_sasa, path = "./plots/complexes/", width=3, height=3)
ggsave(filename="cor_attenuation_sasa_cop9.pdf", plot = cop9_sasa, path = "./plots/complexes/", width=3, height=3)
unlink("cor_attenuation_sasa_cop9.png")
unlink("cor_attenuation_sasa_cop9.pdf")


save(list=ls(), file="./r_workspaces/protein_associations_corum_complexes_processing.RData")
