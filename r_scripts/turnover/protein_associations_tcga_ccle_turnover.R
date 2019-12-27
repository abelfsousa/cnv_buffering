# Abel Sousa 2018

# PhD project
# Study of Control of Protein Abundance using Multi-omics Data from Cancer Samples


# -- Analysis of protein turnover


suppressMessages(library(tidyverse))
suppressMessages(library(ggpubr))
options(bitmapType = "cairo")



# load datasets


# number of times each protein is found and control status
protein_list <- read_tsv("./files/protein_list_control_status.txt")



# import protein attenuation state
protein_attenuation <- read_tsv("./files/protein_attenuation.txt")



# load protein half-lifes (log10)
half_lifes <- read_tsv("./files/protein_half_lifes.txt") %>%
	gather(key="sample", value="log10_half_life", -gene)





# half-life for controlling/controlled/both protein
protein_list_half_life <- protein_list[, -c(2)] %>%
	inner_join(half_lifes, by="gene")



protein_list_half_life_boxpl <- ggplot(data=protein_list_half_life, mapping=aes(x=control_status, y=log10_half_life, fill=control_status)) +
    geom_boxplot(outlier.size=0.5) +
    #geom_jitter(mapping=aes(x=control_status, y=log10_half_life), fill="black", width=0.1, alpha=0.2) +
    stat_compare_means(comparisons = list( c(1, 2), c(2, 3), c(1, 3) )) +
    facet_grid(. ~ sample) +
    theme(axis.title.y=element_text(colour="black", size=12),
          axis.title.x=element_text(colour="black", size=12),
          axis.text.y=element_text(colour="black", size=10),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) +
    labs(x = "Cell type", y = "Protein Half-life (log10)", fill="", title = "")
ggsave("protein_control_status_half_life.png", plot=protein_list_half_life_boxpl, path="./plots/protein_associations_turnover/")
unlink("protein_control_status_half_life.png")



protein_list_half_life_boxpl2 <- ggplot(data=protein_list_half_life %>% group_by(gene) %>% summarise(control_status = unique(control_status), log10_half_life = median(log10_half_life,  na.rm = T)), mapping=aes(x=control_status, y=log10_half_life, fill=control_status)) +
    geom_boxplot(outlier.size=0.5, outlier.shape=NA) +
    geom_jitter(mapping=aes(x=control_status, y=log10_half_life), fill="black", width=0.1, alpha=0.2) +
    stat_compare_means(comparisons = list( c(1, 2), c(2, 3), c(1, 3) )) +
    theme(axis.title.y=element_text(colour="black", size=12),
          axis.title.x=element_text(colour="black", size=12),
          axis.text.y=element_text(colour="black", size=10),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) +
    labs(x = "Control status", y = "Protein Half-life (log10)", fill="", title = "")
ggsave("protein_control_status_half_life2.png", plot=protein_list_half_life_boxpl2, path="./plots/protein_associations_turnover/")
unlink("protein_control_status_half_life2.png")



# half-life by attenuation state
protein_attenuation_half_life <- protein_attenuation[, c("gene", "attenuation", "class")] %>%
  inner_join(half_lifes, by="gene") %>%
  mutate_at(vars(contains("class")), as.factor) %>%
  mutate_if(is.factor, fct_infreq)



protein_attenuation_half_life_boxpl <- ggplot(data=protein_attenuation_half_life, mapping=aes(x=class, y=log10_half_life, fill=class)) +
    geom_boxplot(outlier.size=0.5) +
    #geom_jitter(mapping=aes(x=class, y=log10_half_life), fill="black", width=0.1, alpha=0.1) +
    stat_compare_means(comparisons = list( c(1, 2), c(2, 3), c(1, 3) )) +
    facet_grid(. ~ sample) +
    theme(axis.title.y=element_text(colour="black", size=12),
          axis.title.x=element_text(colour="black", size=12),
          axis.text.y=element_text(colour="black", size=10),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) +
    labs(x = "Attenuation state", y = "Protein Half-life (log10)", fill="", title = "")
ggsave("protein_attenuation_state_half_life.png", plot=protein_attenuation_half_life_boxpl, width=8, path="./plots/protein_associations_turnover/")
unlink("protein_attenuation_state_half_life.png")



protein_attenuation_half_life_boxpl2 <- ggplot(data=protein_attenuation_half_life %>% group_by(gene) %>% summarise(class = unique(class), log10_half_life = median(log10_half_life, na.rm = T)), mapping=aes(x=class, y=log10_half_life, fill=class)) +
    geom_boxplot(outlier.size=0.5, outlier.shape=NA) +
    geom_jitter(mapping=aes(x=class, y=log10_half_life), fill="black", width=0.1, alpha=0.1) +
    stat_compare_means(comparisons = list( c(1, 2), c(2, 3), c(1, 3) )) +
    theme(axis.title.y=element_text(colour="black", size=12),
          axis.title.x=element_text(colour="black", size=12),
          axis.text.y=element_text(colour="black", size=10),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) +
    labs(x = "Attenuation state", y = "Protein Half-life (log10)", fill="", title = "")
ggsave("protein_attenuation_state_half_life2.png", plot=protein_attenuation_half_life_boxpl2, path="./plots/protein_associations_turnover/")
unlink("protein_attenuation_state_half_life2.png")






# significative protein associations with CNV and RNA
signf_protein_associations_cnv_rna <- read_tsv("./files/signf_protein_associations_cnv_mrna_filtered_close_controlling.txt")


# load corum protein pairs
corum_pairs <- read_tsv("./files/corum_protein_pairs.txt")



# select 500 random corum pairs
corum_pairs_random <- corum_pairs[sample(nrow(corum_pairs), 500), ]



# select genes involved in protein interactions
interactive_genes <- read_tsv("./files/biogrid_corum_ER_pairs.txt") %>%
  dplyr::select(1,2) %>%
  unite(interactor.A, interactor.B) %>%
  distinct() %>%
  pull(interactor.A)
length(interactive_genes)
#16583



# get 500 random gene pairs
gene_pairs_random <- data.frame(
  interactor.A = sample(interactive_genes, 500),
  interactor.B = sample(interactive_genes, 500),
  stringsAsFactors = F) %>%
  as.tibble() %>%
  distinct() %>%
  filter(interactor.A != interactor.B) %>%
  mutate(experiment.type = "random")



Reduce(intersect, list(
  paste(signf_protein_associations_cnv_rna$controlling, signf_protein_associations_cnv_rna$controlled, sep="_"),
  paste(corum_pairs_random$interactor.A, corum_pairs_random$interactor.B, sep="_"),
  paste(gene_pairs_random$interactor.A, gene_pairs_random$interactor.B, sep="_"))) %>%
  length()





# correlations of protein half-lifes between protein pairs


# association pairs with CNV
ppairs_half_life_cor <- signf_protein_associations_cnv_rna %>%
  filter(controlling %in% half_lifes$gene & controlled %in% half_lifes$gene) %>%
  dplyr::select(controlling, controlled) %>%
  left_join(half_lifes, by=c("controlling" = "gene")) %>%
  dplyr::rename(controlling_halflife = log10_half_life) %>%
  left_join(half_lifes, by=c("controlled" = "gene", "sample")) %>%
  dplyr::rename(controlled_halflife = log10_half_life) %>%
  group_by(controlling, controlled) %>%
  filter(sum(is.na(controlling_halflife)) < 1 & sum(is.na(controlled_halflife)) < 1) %>%
  do(broom::tidy(cor.test(.$controlling_halflife, .$controlled_halflife, method = "pearson"))) %>%
  dplyr::select(controlling, controlled, estimate, p.value) %>%
  dplyr::rename(r_half_life = estimate) %>%
  ungroup() %>%
  gather(key="statistic", value="value", -c(controlling, controlled)) %>%
  arrange(controlling, controlled) %>%
  mutate(protein_pair = "association pairs with CNV") %>%
  dplyr::select(protein_pair, everything())



# random protein pairs
random_half_life_cor <- gene_pairs_random %>%
  filter(interactor.A %in% half_lifes$gene & interactor.B %in% half_lifes$gene) %>%
  dplyr::select(interactor.A, interactor.B) %>%
  left_join(half_lifes, by=c("interactor.A" = "gene")) %>%
  dplyr::rename(interactor.A_halflife = log10_half_life) %>%
  left_join(half_lifes, by=c("interactor.B" = "gene", "sample")) %>%
  dplyr::rename(interactor.B_halflife = log10_half_life) %>%
  group_by(interactor.A, interactor.B) %>%
  filter(sum(is.na(interactor.A_halflife)) < 1 & sum(is.na(interactor.B_halflife)) < 1) %>%
  do(broom::tidy(cor.test(.$interactor.A_halflife, .$interactor.B_halflife, method = "pearson"))) %>%
  dplyr::select(interactor.A, interactor.B, estimate, p.value) %>%
  dplyr::rename(r_half_life = estimate) %>%
  ungroup() %>%
  gather(key="statistic", value="value", -c(interactor.A, interactor.B)) %>%
  arrange(interactor.A, interactor.B) %>%
  mutate(protein_pair = "random protein pairs") %>%
  dplyr::select(protein_pair, everything())



# random corum pairs
corum_half_life_cor <- corum_pairs_random %>%
  filter(interactor.A %in% half_lifes$gene & interactor.B %in% half_lifes$gene) %>%
  dplyr::select(interactor.A, interactor.B) %>%
  left_join(half_lifes, by=c("interactor.A" = "gene")) %>%
  dplyr::rename(interactor.A_halflife = log10_half_life) %>%
  left_join(half_lifes, by=c("interactor.B" = "gene", "sample")) %>%
  dplyr::rename(interactor.B_halflife = log10_half_life) %>%
  group_by(interactor.A, interactor.B) %>%
  filter(sum(is.na(interactor.A_halflife)) < 1 & sum(is.na(interactor.B_halflife)) < 1) %>%
  do(broom::tidy(cor.test(.$interactor.A_halflife, .$interactor.B_halflife, method = "pearson"))) %>%
  dplyr::select(interactor.A, interactor.B, estimate, p.value) %>%
  dplyr::rename(r_half_life = estimate) %>%
  ungroup() %>%
  gather(key="statistic", value="value", -c(interactor.A, interactor.B)) %>%
  arrange(interactor.A, interactor.B) %>%
  mutate(protein_pair = "corum random pairs") %>%
  dplyr::select(protein_pair, everything())





pairs_half_lifes <- bind_rows(ppairs_half_life_cor[, c(1,4,5)], random_half_life_cor[, c(1,4,5)], corum_half_life_cor[, c(1,4,5)])



pairs_half_lifes_boxpl <- ggplot(data=pairs_half_lifes[pairs_half_lifes$statistic == "r_half_life", ], mapping=aes(x=protein_pair, y=value, fill=protein_pair)) +
    geom_boxplot(outlier.shape=NA) +
    geom_jitter(fill="black", width=0.1, alpha=0.1) +
    stat_compare_means(comparisons = list( c(1, 2), c(2, 3), c(1, 3)) ) +
    theme(axis.title.y=element_text(colour="black", size=12),
          axis.title.x=element_text(colour="black", size=12),
          axis.text.y=element_text(colour="black", size=10),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          plot.title=element_text(colour="black", size=15)) +
    labs(x = "Type of protein pair", y = "Correlation of protein-pairs half-lifes (Pearson r)", fill="", title = "Half-life correlation for protein pairs")
ggsave("pairs_half_lifes.png", plot=pairs_half_lifes_boxpl, path="./plots/protein_associations_turnover/")
unlink("pairs_half_lifes.png")












save(list=ls(), file="./r_workspaces/protein_associations_tcga_ccle_turnover.RData")
