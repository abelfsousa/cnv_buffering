# Abel Sousa 2018

# PhD project
# Study of Control of Protein Abundance using Multi-omics Data from Cancer Samples



# -- Test essentiality of protein association pairs found with CNV across cell lines




suppressMessages(library(tidyverse))
suppressMessages(library(CePa))
suppressMessages(library(ggpubr))
options(bitmapType = "cairo")
set.seed(123)




# significative protein associations with CNV and RNA
signf_protein_associations_cnv_rna <- read_tsv("./files/signf_protein_associations_cnv_mrna_filtered_close_controlling.txt")



# number of times each protein is found and control status
protein_list <- read_tsv("./files/protein_list_control_status.txt")



# protein attenuation
protein_attenuation <- read_tsv("./files/protein_attenuation.txt")



# load DRIVE RNAi data
drive_RNAi_ataris <- readRDS("/nfs/research1/beltrao/memon/databases/cell_lines/rnai_datasets/DRIVE/DRIVE_ATARiS_data.rds")
drive_RNAi_ataris <- drive_RNAi_ataris %>%
	mutate(gene=rownames(drive_RNAi_ataris)) %>%
	as.tibble() %>%
	gather(key="cell_line", value="essentiality_score", -gene)
drive_RNAi_ataris$gene %>% unique %>% length
#6557 genes


drive_RNAi_rsa <- readRDS("/nfs/research1/beltrao/memon/databases/cell_lines/rnai_datasets/DRIVE/DRIVE_RSA_data.rds")
drive_RNAi_rsa <- drive_RNAi_rsa %>%
	mutate(gene=rownames(drive_RNAi_rsa)) %>%
	as.tibble() %>%
	gather(key="cell_line", value="essentiality_score", -gene) %>%
	group_by(gene) %>%
	#filter by essentiality score (lower than -2 in at least 5% of the cell lines)
	filter(sum(essentiality_score <= -2, na.rm = TRUE) > n()*0.05) %>%
	ungroup()
drive_RNAi_rsa$gene %>% unique %>% length
#3554 genes




#achilles_RNAi <- read.gct("/nfs/research1/beltrao/memon/databases/cell_lines/rnai_datasets/Achilles_v2.20.2_GeneSolutions.gct")
#achilles_RNAi <- achilles_RNAi %>%
#	as.data.frame() %>%
#	mutate(gene=rownames(achilles_RNAi)) %>%
#	as.tibble() %>%
#	gather(key="cell_line", value="essentiality_score", -gene)


# load Achilles RNAi DEMETER data
achilles_RNAi <- read_csv("./data/Achilles/ExpandedGeneZSolsCleaned.csv") %>%
	dplyr::rename(gene=X1) %>%
	gather(key="cell_line", value="essentiality_score", -gene)
#17098 genes

achilles_RNAi_sd <- achilles_RNAi %>% pull(essentiality_score) %>% sd(na.rm=T)

#filter by essentiality score (lower than the sd in at least 5% of the cell lines)
achilles_RNAi <- achilles_RNAi %>%
	group_by(gene) %>%
	filter(sum(essentiality_score < -achilles_RNAi_sd, na.rm = TRUE) > n()*0.05) %>%
	ungroup()
achilles_RNAi$gene %>% unique %>% length
#8989 genes





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





# gene essentiality distribution

# rnai_drive_ataris_boxpl <- ggplot(data=drive_RNAi_ataris, mapping=aes(x=cell_line, y=essentiality_score)) +
# 	geom_boxplot(outlier.size=0.1) +
# 	theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
# ggsave("rnai_drive_ataris.png", plot=rnai_drive_ataris_boxpl, path="./plots/protein_associations_essentiality/", width=20)
# unlink("rnai_drive_ataris.png")


# rnai_drive_rsa_boxpl <- ggplot(data=drive_RNAi_rsa, mapping=aes(x=cell_line, y=essentiality_score)) +
# 	geom_boxplot(outlier.size=0.1) +
# 	theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
# ggsave("rnai_drive_rsa.png", plot=rnai_drive_rsa_boxpl, path="./plots/protein_associations_essentiality/", width=20)
# unlink("rnai_drive_rsa.png")


# rnai_achilles_boxpl <- ggplot(data=achilles_RNAi, mapping=aes(x=cell_line, y=essentiality_score)) +
# 	geom_boxplot(outlier.size=0.1) +
# 	theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
# ggsave("rnai_achilles.png", plot=rnai_achilles_boxpl, path="./plots/protein_associations_essentiality/", width=20)
# unlink("rnai_achilles.png")


# crispr_achilles_boxpl <- ggplot(data=achilles_CRISPR, mapping=aes(x=cell_line, y=essentiality_score)) +
# 	geom_boxplot(outlier.size=0.1) +
# 	theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
# ggsave("crispr_achilles.png", plot=crispr_achilles_boxpl, path="./plots/protein_associations_essentiality/", width=20)
# unlink("crispr_achilles.png")




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
#16583 genes




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



# correlations between protein association pairs

ppairs_essentiality_drive_rsa <- signf_protein_associations_cnv_rna %>%
	#filter(abs(beta_cnv) > 0.4) %>%
	#filter(grepl("corum", experiment.type, fixed=T)) %>%
	filter(controlling %in% drive_RNAi_rsa$gene & controlled %in% drive_RNAi_rsa$gene) %>%
	dplyr::select(controlling, controlled) %>%
	left_join(drive_RNAi_rsa, by=c("controlling" = "gene")) %>%
	dplyr::rename(controlling_essent = essentiality_score) %>%
	left_join(drive_RNAi_rsa, by=c("controlled" = "gene", "cell_line")) %>%
	dplyr::rename(controlled_essent = essentiality_score) %>%
	group_by(controlling, controlled) %>%
	do(broom::tidy(cor.test(.$controlling_essent, .$controlled_essent, method = "pearson"))) %>%
	dplyr::select(controlling, controlled, estimate) %>%
	dplyr::rename(r_drive_rnai_rsa = estimate) %>%
	ungroup()


ppairs_essentiality_achilles_rnai <- signf_protein_associations_cnv_rna %>%
	#filter(abs(beta_cnv) > 0.4) %>%
	#filter(grepl("corum", experiment.type, fixed=T)) %>%
	filter(controlling %in% achilles_RNAi$gene & controlled %in% achilles_RNAi$gene) %>%
	dplyr::select(controlling, controlled) %>%
	left_join(achilles_RNAi, by=c("controlling" = "gene")) %>%
	dplyr::rename(controlling_essent = essentiality_score) %>%
	left_join(achilles_RNAi, by=c("controlled" = "gene", "cell_line")) %>%
	dplyr::rename(controlled_essent = essentiality_score) %>%
	na.exclude() %>%
	group_by(controlling, controlled) %>%
	do(broom::tidy(cor.test(.$controlling_essent, .$controlled_essent, method = "pearson"))) %>%
	dplyr::select(controlling, controlled, estimate) %>%
	dplyr::rename(r_achilles_rnai = estimate) %>%
	ungroup()


ppairs_essentiality_achilles_crispr <- signf_protein_associations_cnv_rna %>%
	#filter(abs(beta_cnv) > 0.4) %>%
	#filter(grepl("corum", experiment.type, fixed=T)) %>%
	filter(controlling %in% achilles_CRISPR$gene & controlled %in% achilles_CRISPR$gene) %>%
	dplyr::select(controlling, controlled) %>%
	left_join(achilles_CRISPR, by=c("controlling" = "gene")) %>%
	dplyr::rename(controlling_essent = essentiality_score) %>%
	left_join(achilles_CRISPR, by=c("controlled" = "gene", "cell_line")) %>%
	dplyr::rename(controlled_essent = essentiality_score) %>%
	group_by(controlling, controlled) %>%
	do(broom::tidy(cor.test(.$controlling_essent, .$controlled_essent, method = "pearson"))) %>%
	dplyr::select(controlling, controlled, estimate) %>%
	dplyr::rename(r_achilles_crispr = estimate) %>%
	ungroup()


ppairs_essentiality <- full_join(ppairs_essentiality_drive_rsa, ppairs_essentiality_achilles_rnai, by=c("controlling", "controlled")) %>%
	full_join(ppairs_essentiality_achilles_crispr, by=c("controlling", "controlled")) %>%
	gather(key="data", value="pearson_r", -c(controlling, controlled)) %>%
	arrange(controlling, controlled) %>%
	mutate(protein_pair = "association pairs with CNV") %>%
	dplyr::select(protein_pair, everything(), -controlling, -controlled)




# correlations between random pairs

random_essentiality_drive_rsa <- gene_pairs_random %>%
	filter(interactor.A %in% drive_RNAi_rsa$gene & interactor.B %in% drive_RNAi_rsa$gene) %>%
	left_join(drive_RNAi_rsa, by=c("interactor.A" = "gene")) %>%
	dplyr::rename(interactor.A_essent = essentiality_score) %>%
	left_join(drive_RNAi_rsa, by=c("interactor.B" = "gene", "cell_line")) %>%
	dplyr::rename(interactor.B_essent = essentiality_score) %>%
	group_by(interactor.A, interactor.B) %>%
	do(broom::tidy(cor.test(.$interactor.A_essent, .$interactor.B_essent, method = "pearson"))) %>%
	dplyr::select(interactor.A, interactor.B, estimate) %>%
	dplyr::rename(r_drive_rnai_rsa = estimate) %>%
	ungroup()


random_essentiality_achilles_rnai <- gene_pairs_random %>%
	filter(interactor.A %in% achilles_RNAi$gene & interactor.B %in% achilles_RNAi$gene) %>%
	dplyr::select(interactor.A, interactor.B) %>%
	left_join(achilles_RNAi, by=c("interactor.A" = "gene")) %>%
	dplyr::rename(interactor.A_essent = essentiality_score) %>%
	left_join(achilles_RNAi, by=c("interactor.B" = "gene", "cell_line")) %>%
	dplyr::rename(interactor.B_essent = essentiality_score) %>%
	na.exclude() %>%
	group_by(interactor.A, interactor.B) %>%
	do(broom::tidy(cor.test(.$interactor.A_essent, .$interactor.B_essent, method = "pearson"))) %>%
	dplyr::select(interactor.A, interactor.B, estimate) %>%
	dplyr::rename(r_achilles_rnai = estimate) %>%
	ungroup()


random_essentiality_achilles_crispr <- gene_pairs_random %>%
	filter(interactor.A %in% achilles_CRISPR$gene & interactor.B %in% achilles_CRISPR$gene) %>%
	dplyr::select(interactor.A, interactor.B) %>%
	left_join(achilles_CRISPR, by=c("interactor.A" = "gene")) %>%
	dplyr::rename(interactor.A_essent = essentiality_score) %>%
	left_join(achilles_CRISPR, by=c("interactor.B" = "gene", "cell_line")) %>%
	dplyr::rename(interactor.B_essent = essentiality_score) %>%
	group_by(interactor.A, interactor.B) %>%
	do(broom::tidy(cor.test(.$interactor.A_essent, .$interactor.B_essent, method = "pearson"))) %>%
	dplyr::select(interactor.A, interactor.B, estimate) %>%
	dplyr::rename(r_achilles_crispr = estimate) %>%
	ungroup()


random_essentiality <- full_join(random_essentiality_drive_rsa, random_essentiality_achilles_rnai, by=c("interactor.A", "interactor.B")) %>%
	full_join(random_essentiality_achilles_crispr, by=c("interactor.A", "interactor.B")) %>%
	gather(key="data", value="pearson_r", -c(interactor.A, interactor.B)) %>%
	arrange(interactor.A, interactor.B) %>%
	mutate(protein_pair = "random protein pairs") %>%
	dplyr::select(protein_pair, everything(), -interactor.A, -interactor.B)




# correlations between random corum pairs

corum_essentiality_drive_rsa <- corum_pairs_random %>%
	filter(interactor.A %in% drive_RNAi_rsa$gene & interactor.B %in% drive_RNAi_rsa$gene) %>%
	left_join(drive_RNAi_rsa, by=c("interactor.A" = "gene")) %>%
	dplyr::rename(interactor.A_essent = essentiality_score) %>%
	left_join(drive_RNAi_rsa, by=c("interactor.B" = "gene", "cell_line")) %>%
	dplyr::rename(interactor.B_essent = essentiality_score) %>%
	group_by(interactor.A, interactor.B) %>%
	do(broom::tidy(cor.test(.$interactor.A_essent, .$interactor.B_essent, method = "pearson"))) %>%
	dplyr::select(interactor.A, interactor.B, estimate) %>%
	dplyr::rename(r_drive_rnai_rsa = estimate) %>%
	ungroup()


corum_essentiality_achilles_rnai <- corum_pairs_random %>%
	filter(interactor.A %in% achilles_RNAi$gene & interactor.B %in% achilles_RNAi$gene) %>%
	dplyr::select(interactor.A, interactor.B) %>%
	left_join(achilles_RNAi, by=c("interactor.A" = "gene")) %>%
	dplyr::rename(interactor.A_essent = essentiality_score) %>%
	left_join(achilles_RNAi, by=c("interactor.B" = "gene", "cell_line")) %>%
	dplyr::rename(interactor.B_essent = essentiality_score) %>%
	na.exclude() %>%
	group_by(interactor.A, interactor.B) %>%
	do(broom::tidy(cor.test(.$interactor.A_essent, .$interactor.B_essent, method = "pearson"))) %>%
	dplyr::select(interactor.A, interactor.B, estimate) %>%
	dplyr::rename(r_achilles_rnai = estimate) %>%
	ungroup()


corum_essentiality_achilles_crispr <- corum_pairs_random %>%
	filter(interactor.A %in% achilles_CRISPR$gene & interactor.B %in% achilles_CRISPR$gene) %>%
	dplyr::select(interactor.A, interactor.B) %>%
	left_join(achilles_CRISPR, by=c("interactor.A" = "gene")) %>%
	dplyr::rename(interactor.A_essent = essentiality_score) %>%
	left_join(achilles_CRISPR, by=c("interactor.B" = "gene", "cell_line")) %>%
	dplyr::rename(interactor.B_essent = essentiality_score) %>%
	group_by(interactor.A, interactor.B) %>%
	do(broom::tidy(cor.test(.$interactor.A_essent, .$interactor.B_essent, method = "pearson"))) %>%
	dplyr::select(interactor.A, interactor.B, estimate) %>%
	dplyr::rename(r_achilles_crispr = estimate) %>%
	ungroup()


corum_essentiality <- full_join(corum_essentiality_drive_rsa, corum_essentiality_achilles_rnai, by=c("interactor.A", "interactor.B")) %>%
	full_join(corum_essentiality_achilles_crispr, by=c("interactor.A", "interactor.B")) %>%
	gather(key="data", value="pearson_r", -c(interactor.A, interactor.B)) %>%
	arrange(interactor.A, interactor.B) %>%
	mutate(protein_pair = "corum random pairs") %>%
	dplyr::select(protein_pair, everything(), -interactor.A, -interactor.B)



essentiality <- bind_rows(ppairs_essentiality, random_essentiality, corum_essentiality)



essentiality_boxpl <- ggplot(data=essentiality, mapping=aes(x=protein_pair, y=pearson_r, fill=protein_pair)) +
    geom_boxplot(outlier.shape=NA) +
    geom_jitter(fill="black", width=0.1, alpha=0.1) +
    stat_compare_means(comparisons = list( c(1, 2), c(2, 3), c(1, 3)) ) +
    facet_grid(. ~ data) +
    theme(axis.title.y=element_text(colour="black", size=15),
          axis.title.x=element_text(colour="black", size=15),
          axis.text.y=element_text(colour="black", size=12),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          plot.title=element_text(colour="black", size=15),
          strip.text.x=element_text(size=15),
          legend.text=element_text(colour="black", size=12)) +
    labs(x = "Type of protein pair", y = "Correlation of protein-pairs essentiality (Pearson r)", fill="", title = "Correlation of essentiality profiles for gene pairs")
ggsave("protein_pairs_cnv_essentiality.png", plot=essentiality_boxpl, width=9, path="./plots/protein_associations_essentiality/")
unlink("protein_pairs_cnv_essentiality.png")





# gene essentiality by controlling status for RNAi and Crispr studies

protein_essentiality_drive_rnai_rsa <- protein_list %>%
	dplyr::select(-interactions) %>%
	inner_join(drive_RNAi_rsa, by="gene") %>%
	group_by(gene) %>%
	summarise(control_status = unique(control_status), rnai_drive_RSA = median(essentiality_score, na.rm = TRUE))

protein_essentiality_achilles_rnai <- protein_list %>%
	dplyr::select(-interactions) %>%
	inner_join(achilles_RNAi, by="gene") %>%
	group_by(gene) %>%
	summarise(control_status = unique(control_status), rnai_achilles = median(essentiality_score, na.rm = TRUE))

protein_essentiality_achilles_crispr <- protein_list %>%
	dplyr::select(-interactions) %>%
	inner_join(achilles_CRISPR, by="gene") %>%
	group_by(gene) %>%
	summarise(control_status = unique(control_status), crispr_ceres = median(essentiality_score, na.rm = TRUE))



protein_list_essentiality <- full_join(protein_essentiality_drive_rnai_rsa, protein_essentiality_achilles_rnai, by=c("gene", "control_status")) %>%
	full_join(protein_essentiality_achilles_crispr, by=c("gene", "control_status")) %>%
	gather(key="study", value="median_essentiality_score", -c(gene, control_status))


protein_list_essentiality_drive_rnai_rsa_boxpl <- ggplot(data=protein_list_essentiality[protein_list_essentiality$study == "rnai_drive_RSA", ], mapping=aes(x=control_status, y=median_essentiality_score, fill=control_status)) +
    geom_boxplot(outlier.shape=NA) +
    geom_jitter(fill="black", width=0.1, alpha=0.2) +
    stat_compare_means(comparisons = list( c(1, 2), c(2, 3), c(1, 3)) ) +
    theme(axis.title.y=element_text(colour="black", size=12),
          axis.title.x=element_text(colour="black", size=12),
          axis.text.y=element_text(colour="black", size=10),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) +
    labs(x = "Control status", y = "Median gene essentiality across cell lines", fill="", title = "DRIVE RNAi")
ggsave("protein_list_essentiality_drive_rnai_rsa.png", plot=protein_list_essentiality_drive_rnai_rsa_boxpl, path="./plots/protein_associations_essentiality/", width=4)
unlink("protein_list_essentiality_drive_rnai_rsa.png")

protein_list_essentiality_achilles_rnai_boxpl <- ggplot(data=protein_list_essentiality[protein_list_essentiality$study == "rnai_achilles", ], mapping=aes(x=control_status, y=median_essentiality_score, fill=control_status)) +
    geom_boxplot(outlier.shape=NA) +
    geom_jitter(fill="black", width=0.1, alpha=0.2) +
    stat_compare_means(comparisons = list( c(1, 2), c(2, 3), c(1, 3)) ) +
    theme(axis.title.y=element_text(colour="black", size=12),
          axis.title.x=element_text(colour="black", size=12),
          axis.text.y=element_text(colour="black", size=10),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) +
    labs(x = "Control status", y = "Median gene essentiality across cell lines", fill="", title = "Achilles RNAi")
ggsave("protein_list_essentiality_achilles_rnai.png", plot=protein_list_essentiality_achilles_rnai_boxpl, path="./plots/protein_associations_essentiality/", width=4)
unlink("protein_list_essentiality_achilles_rnai.png")

protein_list_essentiality_achilles_crispr_boxpl <- ggplot(data=protein_list_essentiality[protein_list_essentiality$study == "crispr_ceres", ], mapping=aes(x=control_status, y=median_essentiality_score, fill=control_status)) +
    geom_boxplot(outlier.shape=NA) +
    geom_jitter(fill="black", width=0.1, alpha=0.2) +
    stat_compare_means(comparisons = list( c(1, 2), c(2, 3), c(1, 3)) ) +
    theme(axis.title.y=element_text(colour="black", size=12),
          axis.title.x=element_text(colour="black", size=12),
          axis.text.y=element_text(colour="black", size=10),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) +
    labs(x = "Control status", y = "Median gene essentiality across cell lines", fill="", title = "Achilles CRISPR")
ggsave("protein_list_essentiality_achilles_crispr.png", plot=protein_list_essentiality_achilles_crispr_boxpl, path="./plots/protein_associations_essentiality/", width=4)
unlink("protein_list_essentiality_achilles_crispr.png")






# gene essentiality by attenuation state for RNAi and Crispr studies

protein_attenuation_drive_rnai_rsa <- protein_attenuation %>%
	dplyr::select(gene, class) %>%
	inner_join(drive_RNAi_rsa, by="gene") %>%
	group_by(gene) %>%
	summarise(class = unique(class), rnai_drive_RSA = median(essentiality_score, na.rm = TRUE))

protein_attenuation_achilles_rnai <- protein_attenuation %>%
	dplyr::select(gene, class) %>%
	inner_join(achilles_RNAi, by="gene") %>%
	group_by(gene) %>%
	summarise(class = unique(class), rnai_achilles = median(essentiality_score, na.rm = TRUE))

protein_attenuation_achilles_crispr <- protein_attenuation %>%
	dplyr::select(gene, class) %>%
	inner_join(achilles_CRISPR, by="gene") %>%
	group_by(gene) %>%
	summarise(class = unique(class), crispr_ceres = median(essentiality_score, na.rm = TRUE))



protein_attenuation_median_essentiality <- full_join(protein_attenuation_drive_rnai_rsa, protein_attenuation_achilles_rnai, by=c("gene", "class")) %>%
	full_join(protein_attenuation_achilles_crispr, by=c("gene", "class")) %>%
	gather(key="study", value="median_essentiality_score", -c(gene, class)) %>%
	mutate_at(vars(contains("class")), as.factor) %>%
	mutate_if(is.factor, fct_infreq)


protein_attenuation_median_essentiality_drive_rnai_rsa_boxpl <- ggplot(data=protein_attenuation_median_essentiality[protein_attenuation_median_essentiality$study == "rnai_drive_RSA", ], mapping=aes(x=class, y=median_essentiality_score, fill=class)) +
    geom_boxplot(outlier.shape=NA) +
    geom_jitter(fill="black", width=0.1, alpha=0.2) +
    stat_compare_means(comparisons = list( c(1, 2), c(2, 3), c(1, 3)) ) +
    theme(axis.title.y=element_text(colour="black", size=12),
          axis.title.x=element_text(colour="black", size=12),
          axis.text.y=element_text(colour="black", size=10),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) +
    labs(x = "Attenuation state", y = "Median gene essentiality across cell lines", fill="", title = "DRIVE RNAi")
ggsave("protein_attenuation_median_essentiality_drive_rnai_rsa.png", plot=protein_attenuation_median_essentiality_drive_rnai_rsa_boxpl, path="./plots/protein_associations_essentiality/", width=4)
unlink("protein_attenuation_median_essentiality_drive_rnai_rsa.png")

protein_attenuation_median_essentiality_achilles_rnai_boxpl <- ggplot(data=protein_attenuation_median_essentiality[protein_attenuation_median_essentiality$study == "rnai_achilles", ], mapping=aes(x=class, y=median_essentiality_score, fill=class)) +
    geom_boxplot(outlier.shape=NA) +
    geom_jitter(fill="black", width=0.1, alpha=0.2) +
    stat_compare_means(comparisons = list( c(1, 2), c(2, 3), c(1, 3)) ) +
    theme(axis.title.y=element_text(colour="black", size=12),
          axis.title.x=element_text(colour="black", size=12),
          axis.text.y=element_text(colour="black", size=10),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) +
    labs(x = "Attenuation state", y = "Median gene essentiality across cell lines", fill="", title = "Achilles RNAi")
ggsave("protein_attenuation_median_essentiality_achilles_rnai.png", plot=protein_attenuation_median_essentiality_achilles_rnai_boxpl, path="./plots/protein_associations_essentiality/", width=4)
unlink("protein_attenuation_median_essentiality_achilles_rnai.png")

protein_attenuation_median_essentiality_achilles_crispr_boxpl <- ggplot(data=na.omit(protein_attenuation_median_essentiality[protein_attenuation_median_essentiality$study == "crispr_ceres", ]), mapping=aes(x=class, y=median_essentiality_score, fill=class)) +
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(color="black", width=0.1, alpha=0.05) +
  stat_compare_means(comparisons = list( c(1, 3)) ) +
	theme_classic() +
  theme(axis.title.y=element_text(colour="black", size=15),
		axis.title.x=element_text(colour="black", size=15),
    axis.text.y=element_text(colour="black", size=12),
		axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
		plot.title=element_blank(),
		legend.text=element_text(colour="black", size=12),
		legend.title=element_text(colour="black", size=15)) +
	scale_fill_brewer(type = "seq", palette = "Blues", labels = c("non-attenuated", "lowly-attenuated", "highly-attenuated"), name = "Attenuation") +
  labs(x = "Attenuation", y = "Median CRISPR-Cas9 essentiality score")
ggsave("protein_attenuation_median_essentiality_achilles_crispr.png", plot=protein_attenuation_median_essentiality_achilles_crispr_boxpl, path="./plots/protein_associations_essentiality/", width=5, height=6)
ggsave("protein_attenuation_median_essentiality_achilles_crispr.pdf", plot=protein_attenuation_median_essentiality_achilles_crispr_boxpl, path="./plots/protein_associations_essentiality/", width=5, height=6)
unlink("protein_attenuation_median_essentiality_achilles_crispr.png")
unlink("protein_attenuation_median_essentiality_achilles_crispr.pdf")












# number of cell lines in which a gene is essential by controlling status for RNAi and Crispr studies

number_cells_essential_achilles_rnai <- protein_list %>%
	dplyr::select(-interactions) %>%
	inner_join(achilles_RNAi %>% group_by(gene) %>% summarise(n_cell_lines = sum(essentiality_score <= -3, na.rm=T)) %>% filter(n_cell_lines!=0) %>% ungroup(), by="gene") %>%
	mutate(project = "achilles_rnai")


number_cells_essential_achilles_crispr <- protein_list %>%
	dplyr::select(-interactions) %>%
	inner_join(achilles_CRISPR %>% group_by(gene) %>% summarise(n_cell_lines = sum(essentiality_score <= -1, na.rm=T)) %>% filter(n_cell_lines!=0) %>% ungroup(), by="gene") %>%
	mutate(project = "achilles_crispr")


number_cells_essential_drive_rnai <- protein_list %>%
	dplyr::select(-interactions) %>%
	inner_join(drive_RNAi_rsa %>% group_by(gene) %>% summarise(n_cell_lines = sum(essentiality_score <= -3, na.rm=T)) %>% filter(n_cell_lines!=0) %>% ungroup(), by="gene") %>%
	mutate(project = "drive_rnai")



number_cells_essential <- bind_rows(number_cells_essential_achilles_rnai, number_cells_essential_achilles_crispr, number_cells_essential_drive_rnai)




number_cells_essential_achilles_rnai_boxpl <- ggplot(data=number_cells_essential_achilles_rnai, mapping=aes(x=control_status, y=n_cell_lines, fill=control_status)) +
    geom_boxplot(outlier.shape=NA) +
    geom_jitter(fill="black", width=0.1, alpha=0.2) +
    stat_compare_means(comparisons = list( c(1, 2), c(2, 3), c(1, 3)) ) +
    theme(axis.title.y=element_text(colour="black", size=15),
          axis.title.x=element_text(colour="black", size=15),
          axis.text.y=element_text(colour="black", size=13),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) +
    labs(x = "Control status", y = "Number of cell lines in which the gene is essential", fill="", title = "Achilles RNAi")
ggsave("number_cells_essential_achilles_rnai.png", plot=number_cells_essential_achilles_rnai_boxpl, path="./plots/protein_associations_essentiality/")
unlink("number_cells_essential_achilles_rnai.png")




number_cells_essential_drive_rnai_boxpl <- ggplot(data=number_cells_essential_drive_rnai, mapping=aes(x=control_status, y=n_cell_lines, fill=control_status)) +
    geom_boxplot(outlier.shape=NA) +
    geom_jitter(fill="black", width=0.1, alpha=0.2) +
    stat_compare_means(comparisons = list( c(1, 2), c(2, 3), c(1, 3)) ) +
    theme(axis.title.y=element_text(colour="black", size=15),
          axis.title.x=element_text(colour="black", size=15),
          axis.text.y=element_text(colour="black", size=13),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) +
    labs(x = "Control status", y = "Number of cell lines in which the gene is essential", fill="", title = "DRIVE RNAi")
ggsave("number_cells_essential_drive_rnai.png", plot=number_cells_essential_drive_rnai_boxpl, path="./plots/protein_associations_essentiality/")
unlink("number_cells_essential_drive_rnai.png")




number_cells_essential_achilles_crispr_boxpl <- ggplot(data=number_cells_essential_achilles_crispr, mapping=aes(x=control_status, y=n_cell_lines, fill=control_status)) +
    geom_boxplot(outlier.shape=NA) +
    geom_jitter(fill="black", width=0.1, alpha=0.2) +
    stat_compare_means(comparisons = list( c(1, 2), c(2, 3), c(1, 3)) ) +
    theme(axis.title.y=element_text(colour="black", size=15),
          axis.title.x=element_text(colour="black", size=15),
          axis.text.y=element_text(colour="black", size=13),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) +
    labs(x = "Control status", y = "Number of cell lines in which the gene is essential", fill="", title = "Achilles CRISPR")
ggsave("number_cells_essential_achilles_crispr.png", plot=number_cells_essential_achilles_crispr_boxpl, path="./plots/protein_associations_essentiality/")
unlink("number_cells_essential_achilles_crispr.png")





number_cells_essential_boxpl <- ggplot(data=number_cells_essential, mapping=aes(x=control_status, y=n_cell_lines, fill=control_status)) +
    geom_boxplot(outlier.shape=NA) +
    geom_jitter(fill="black", width=0.1, alpha=0.2) +
    stat_compare_means(comparisons = list( c(1, 2), c(2, 3), c(1, 3)) ) +
    facet_grid(. ~ project) +
    theme(axis.title.y=element_text(colour="black", size=15),
          axis.title.x=element_text(colour="black", size=15),
          axis.text.y=element_text(colour="black", size=13),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          strip.text.x=element_text(size=15)) +
    labs(x = "Control status", y = "Number of cell lines in which the gene is essential", fill="", title = "")
ggsave("number_cells_essential_achilles_by_project.png", plot=number_cells_essential_boxpl, path="./plots/protein_associations_essentiality/")
unlink("number_cells_essential_achilles_by_project.png")






# number of cell lines in which a gene is essential by attenuation state for RNAi and Crispr studies

attenuation_number_cells_essential_achilles_rnai <- protein_attenuation %>%
	dplyr::select(gene, class) %>%
	inner_join(achilles_RNAi %>% group_by(gene) %>% summarise(n_cell_lines = sum(essentiality_score <= -3, na.rm=T)) %>% filter(n_cell_lines!=0) %>% ungroup(), by="gene") %>%
	mutate(project = "achilles_rnai")


attenuation_number_cells_essential_achilles_crispr <- protein_attenuation %>%
	dplyr::select(gene, class) %>%
	inner_join(achilles_CRISPR %>% group_by(gene) %>% summarise(n_cell_lines = sum(essentiality_score <= -1, na.rm=T)) %>% filter(n_cell_lines!=0) %>% ungroup(), by="gene") %>%
	mutate(project = "achilles_crispr")


attenuation_number_cells_essential_drive_rnai <- protein_attenuation %>%
	dplyr::select(gene, class) %>%
	inner_join(drive_RNAi_rsa %>% group_by(gene) %>% summarise(n_cell_lines = sum(essentiality_score <= -3, na.rm=T)) %>% filter(n_cell_lines!=0) %>% ungroup(), by="gene") %>%
	mutate(project = "drive_rnai")



attenuation_number_cells_essential <- bind_rows(attenuation_number_cells_essential_achilles_rnai, attenuation_number_cells_essential_achilles_crispr, attenuation_number_cells_essential_drive_rnai) %>%
	mutate_at(vars(contains("class")), as.factor) %>%
	mutate_if(is.factor, fct_infreq)



attenuation_number_cells_essential_boxpl <- ggplot(data=attenuation_number_cells_essential %>% filter(project == "achilles_crispr"), mapping=aes(x=class, y=n_cell_lines, fill=class)) +
  geom_boxplot() +
  #geom_jitter(fill="black", width=0.1, alpha=0.2) +
  #stat_compare_means(comparisons = list( c(1, 2), c(2, 3), c(1, 3)) ) +
  #facet_grid(. ~ project) +
	theme_classic() +
  theme(axis.title=element_text(colour="black", size=12),
    axis.text.y=element_text(colour="black", size=10),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    plot.title = element_blank()) +
	scale_fill_brewer(type = "seq", palette = "Oranges", labels = c("non-attenuated", "lowly-attenuated", "highly-attenuated"), name = "") +
  labs(x = "Attenuation", y = "Number of cell lines in which the gene is essential")
ggsave("attenuation_number_cells_essential_achilles_by_project.png", plot=attenuation_number_cells_essential_boxpl, path="./plots/protein_associations_essentiality/", width=5, height=5)
unlink("attenuation_number_cells_essential_achilles_by_project.png")




# cell lines CNV data
ccle_cnv <- read_tsv("./data/ccle/data_CNA.txt") %>%
	dplyr::rename(gene = Hugo_Symbol) %>%
	dplyr::select(-c(Entrez_Gene_Id, Cytoband)) %>%
	gather(key = "cell_line", value = "cnv_gistic2", -gene) %>%
	filter(cell_line %in% achilles_CRISPR$cell_line) %>%
	group_by(gene) %>%
	summarise(amplified = sum(cnv_gistic2 > 0)) %>%
	ungroup()



attenuation_number_cells_amplified <- protein_attenuation %>%
	dplyr::select(gene, class) %>%
	inner_join(ccle_cnv, by="gene")


attenuation_number_cells_amplified_boxpl <- ggplot(data=attenuation_number_cells_amplified, mapping=aes(x=class, y=amplified, fill=class)) +
	geom_boxplot() +
	#geom_jitter(fill="black", width=0.1, alpha=0.2) +
	#stat_compare_means(comparisons = list( c(1, 2), c(2, 3), c(1, 3)) ) +
	theme_classic() +
	theme(axis.title=element_text(colour="black", size=12),
		axis.text.y=element_text(colour="black", size=10),
	  axis.text.x=element_blank(),
	  axis.ticks.x=element_blank(),
	  plot.title = element_blank()) +
	scale_fill_brewer(type = "seq", palette = "Oranges", labels = c("non-attenuated", "lowly-attenuated", "highly-attenuated"), name = "") +
	labs(x = "Attenuation", y = "Number of cell lines in which the gene is amplified")
ggsave("attenuation_number_cells_amplified.png", plot=attenuation_number_cells_amplified_boxpl, path="./plots/protein_associations_essentiality/", width=5, height=5)
unlink("attenuation_number_cells_amplified.png")




save(list=ls(), file="./r_workspaces/protein_associations_tcga_ccle_essentiality.RData")
