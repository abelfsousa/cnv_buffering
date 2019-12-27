# Abel Sousa 2018

# PhD project
# Study of Control of Protein Abundance using Multi-omics Data from Cancer Samples

# -- Number of GTEx eQTLs by control status and tissue

# -- Gene expression variance by control status and tissue




suppressMessages(library(tidyverse))
suppressMessages(library(ggpubr))
suppressMessages(library(data.table))

options(bitmapType = "cairo")



# import number of protein interactions and control status
protein_list <- read_tsv("./files/protein_list_control_status.txt")



# import protein attenuation state
protein_attenuation <- read_tsv("./files/protein_attenuation.txt") %>%
  mutate_at(vars(contains("class")), as.factor) %>%
  mutate_if(is.factor, fct_infreq)



# import GTEx eQTLs for all tissue types
gtex_eqtls <- fread("./data/gtex/all_tissues.v6p.signif_snpgene_pairs.txt")
gtex_eqtls <- gtex_eqtls %>%
  as.tibble() %>%
  dplyr::rename(tissue = Adipose_Subcutaneous)



# load GENCODE v19 annotation
gene_annot_v19 <- read_tsv("./data/gtex/geneAnnot.gencode.v19.txt")

gtex_eqtls <- gtex_eqtls %>%
  inner_join(gene_annot_v19, by=c("gene_id" = "geneID")) %>%
  dplyr::select(-c(geneType, chrom)) %>%
  dplyr::rename(gene_name = geneName) %>%
  dplyr::select(variant_id, gene_id, gene_name, everything())



gene_list_control1 <- gtex_eqtls %>%
  dplyr::select(variant_id, gene_name, tissue) %>%
  inner_join(protein_list %>% dplyr::select(gene, control_status), by=c("gene_name" = "gene")) %>%
  group_by(tissue) %>%
  mutate(eQTLs_tissue = n()) %>%
  mutate(genes_tissue = length(unique(gene_name))) %>%
  group_by(tissue, control_status) %>%
  summarise(eQTLs_perc = n()/unique(eQTLs_tissue)) %>%
  ungroup() %>%
  ggplot(mapping = aes(x=tissue, y = eQTLs_perc, group=control_status, fill=control_status)) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_classic() +
		theme(axis.title = element_text(colour="black", size=15),
    		axis.text.y = element_text(colour="black", size=13),
    		axis.text.x = element_text(colour="black", size=7, angle=45, hjust=1)) +
  labs(x = "Tissue", y = "Number of eQTLs", title = "")
ggsave(filename="gtex_tissues_eQTLs_control_status.png", plot=gene_list_control1, path = "./plots/protein_associations_normal")
unlink("gtex_tissues_eQTLs_control_status.png")


gene_list_attenuation1 <- gtex_eqtls %>%
  dplyr::select(variant_id, gene_name, tissue) %>%
  inner_join(protein_attenuation %>% dplyr::select(gene, class), by=c("gene_name" = "gene")) %>%
  group_by(tissue) %>%
  mutate(eQTLs_tissue = n()) %>%
  mutate(genes_tissue = length(unique(gene_name))) %>%
  group_by(tissue, class) %>%
  summarise(eQTLs_perc = n()/unique(eQTLs_tissue)) %>%
  ungroup() %>%
  ggplot(mapping = aes(x=tissue, y = eQTLs_perc, group=class, fill=class)) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_classic() +
		theme(axis.title = element_text(colour="black", size=15),
    		axis.text.y = element_text(colour="black", size=13),
    		axis.text.x = element_text(colour="black", size=7, angle=45, hjust=1)) +
  labs(x = "Tissue", y = "Number of eQTLs", title = "")
ggsave(filename="gtex_tissues_eQTLs_attenuation_state.png", plot=gene_list_attenuation1, path = "./plots/protein_associations_normal")
unlink("gtex_tissues_eQTLs_attenuation_state.png")


gene_list_control2 <- gtex_eqtls %>%
  dplyr::select(variant_id, gene_name, tissue) %>%
  inner_join(protein_list %>% dplyr::select(gene, control_status), by=c("gene_name" = "gene")) %>%
  group_by(gene_name) %>%
  summarise(control_status = unique(control_status), n_eQTLs = n()) %>%
  inner_join(protein_list %>% group_by(control_status) %>% summarise(n_genes = n()), by=c("control_status")) %>%
  group_by(control_status) %>%
  summarise(eQTLs_perc = n()/unique(n_genes)) %>%
  ggplot(mapping = aes(x=control_status, y=eQTLs_perc, fill=control_status)) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_classic() +
		theme(axis.title = element_text(colour="black", size=15),
    		axis.text.y = element_text(colour="black", size=13),
    		axis.text.x = element_text(colour="black", size=13)) +
    scale_y_continuous(limits = c(0,1)) +
  labs(x = "Control status", y = "Percentage of eQTLs", title = "")
ggsave(filename="gtex_tissues_eQTLs_control_status2.png", plot=gene_list_control2, path = "./plots/protein_associations_normal")
unlink("gtex_tissues_eQTLs_control_status2.png")


gene_list_attenuation2 <- gtex_eqtls %>%
  dplyr::select(variant_id, gene_name, tissue) %>%
  inner_join(protein_attenuation %>% dplyr::select(gene, class), by=c("gene_name" = "gene")) %>%
  group_by(gene_name) %>%
  summarise(class = unique(class), n_eQTLs = n()) %>%
  inner_join(protein_attenuation %>% group_by(class) %>% summarise(n_genes = n()), by=c("class")) %>%
  group_by(class) %>%
  summarise(eQTLs_perc = n()/unique(n_genes)) %>%
  ggplot(mapping = aes(x=class, y=eQTLs_perc, fill=class)) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_classic() +
		theme(axis.title = element_text(colour="black", size=15),
    		axis.text.y = element_text(colour="black", size=13),
    		axis.text.x = element_blank()) +
    scale_y_continuous(limits = c(0,1)) +
  labs(x = "Attenuation state", y = "Percentage of eQTLs", title = "")
ggsave(filename="gtex_tissues_eQTLs_attenuation_state2.png", plot=gene_list_attenuation2, path = "./plots/protein_associations_normal")
unlink("gtex_tissues_eQTLs_attenuation_state2.png")


gene_list_control3 <- gtex_eqtls %>%
  dplyr::select(variant_id, gene_name, tissue) %>%
  inner_join(protein_list %>% dplyr::select(gene, control_status), by=c("gene_name" = "gene")) %>%
  group_by(control_status) %>%
  summarise(n_eQTLs = n()) %>%
  inner_join(protein_list %>% group_by(control_status) %>% summarise(n_genes = n()), by=c("control_status")) %>%
  mutate(eQTLs_perc = n_eQTLs/n_genes) %>%
  ggplot(mapping = aes(x=control_status, y=eQTLs_perc, fill=control_status)) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_classic() +
		theme(axis.title = element_text(colour="black", size=15),
    		axis.text.y = element_text(colour="black", size=13),
    		axis.text.x = element_text(colour="black", size=13)) +
  labs(x = "Control status", y = "Nº eQTLs per gene", title = "")
ggsave(filename="gtex_tissues_eQTLs_control_status3.png", plot=gene_list_control3, path = "./plots/protein_associations_normal")
unlink("gtex_tissues_eQTLs_control_status3.png")


gene_list_attenuation3 <- gtex_eqtls %>%
  dplyr::select(variant_id, gene_name, tissue) %>%
  inner_join(protein_attenuation %>% dplyr::select(gene, class), by=c("gene_name" = "gene")) %>%
  group_by(class) %>%
  summarise(n_eQTLs = n()) %>%
  inner_join(protein_attenuation %>% group_by(class) %>% summarise(n_genes = n()), by=c("class")) %>%
  mutate(eQTLs_perc = n_eQTLs/n_genes) %>%
  ggplot(mapping = aes(x=class, y=eQTLs_perc, fill=class)) +
    theme_classic() +
    geom_bar(stat = "identity", position = "dodge") +
		theme(axis.title = element_text(colour="black", size=15),
    		axis.text.y = element_text(colour="black", size=13),
    		axis.text.x = element_blank()) +
  labs(x = "Attenuation state", y = "Nº eQTLs per gene", title = "")
ggsave(filename="gtex_tissues_eQTLs_attenuation_state3.png", plot=gene_list_attenuation3, path = "./plots/protein_associations_normal")
unlink("gtex_tissues_eQTLs_attenuation_state3.png")


gene_list_control4 <- gtex_eqtls %>%
  dplyr::select(variant_id, gene_name, tissue) %>%
  inner_join(protein_list %>% dplyr::select(gene, control_status), by=c("gene_name" = "gene")) %>%
  group_by(gene_name) %>%
  summarise(control_status = unique(control_status), n_eQTLs = n()) %>%
  ggplot(mapping = aes(x=control_status, y=n_eQTLs, fill=control_status)) +
    geom_boxplot() +
    theme_classic() +
    stat_compare_means(comparisons = list( c(1, 2), c(2, 3), c(1, 3) )) +
		theme(axis.title = element_text(colour="black", size=15),
    		axis.text.y = element_text(colour="black", size=13),
    		axis.text.x = element_text(colour="black", size=13)) +
    scale_y_continuous(trans = "log2") +
  labs(x = "Control status", y = "Number of eQTLs (log2)", title = "")
ggsave(filename="gtex_tissues_eQTLs_control_status4.png", plot=gene_list_control4, path = "./plots/protein_associations_normal")
unlink("gtex_tissues_eQTLs_control_status4.png")


gene_list_attenuation4 <- gtex_eqtls %>%
  dplyr::select(variant_id, gene_name, tissue) %>%
  inner_join(protein_attenuation %>% dplyr::select(gene, class), by=c("gene_name" = "gene")) %>%
  group_by(gene_name) %>%
  summarise(class = unique(class), n_eQTLs = n()) %>%
  ggplot(mapping = aes(x=class, y=n_eQTLs, fill=class)) +
    geom_boxplot() +
    theme_classic() +
    stat_compare_means(comparisons = list( c(1, 2), c(2, 3), c(1, 3), c(1,4), c(2,4), c(3,4) )) +
		theme(axis.title = element_text(colour="black", size=15),
    		axis.text.y = element_text(colour="black", size=13),
    		axis.text.x = element_blank()) +
    scale_y_continuous(trans = "log2") +
    labs(x = "Attenuation state", y = "Number of eQTLs (log2)", title = "")
ggsave(filename="gtex_tissues_eQTLs_attenuation_state4.png", plot=gene_list_attenuation4, path = "./plots/protein_associations_normal")
unlink("gtex_tissues_eQTLs_attenuation_state4.png")
