# Abel Sousa 2018

# PhD project
# Study of Control of Protein Abundance using Multi-omics Data from Cancer Samples

# -- Analysis of protein associations using normal samples of GTEx and Human Proteome Map

# -- Correlations of protein/RNA expression across normal tissues




suppressMessages(library(tidyverse))
suppressMessages(library(ggpubr))
suppressMessages(library(biomaRt))
suppressMessages(library(limma))

options(bitmapType = "cairo")
set.seed(123)



# non-significative protein associations with CNV
non_signf_protein_associations_cnv <- read_tsv("./files/protein_associations_tcga_ccle_cnv_reg.txt") %>%
    filter(fdr > 0.05)


# non-significative protein associations with RNA
non_signf_protein_associations_rna <- read_tsv("./files/protein_associations_tcga_ccle_mRNA_reg.txt") %>%
    filter(fdr > 0.05)


# non-significative protein associations with CNV and RNA
non_signf_protein_associations_cnv_rna <- inner_join(
	non_signf_protein_associations_cnv %>% dplyr::select(controlling, controlled, beta_cnv, fdr_cnv = fdr),
	non_signf_protein_associations_rna %>% dplyr::select(controlling, controlled, beta_rna, fdr_rna = fdr),
	by=c("controlling", "controlled"))


# load significative protein associations with CNV and RNA
protein_associations_rna_cnv <- read_tsv("./files/signf_protein_associations_cnv_mrna_filtered_close_controlling.txt")


# load protein attenuation
protein_attenuation <- read_tsv("./files/protein_attenuation.txt")



# load GTEx median expression per tissue
gtex_tissue <- read_tsv("./data/gtex/GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_median_rpkm.gct", skip=2) %>%
	dplyr::select(
		Gene_Name = `Description`,
		frontal_cortex = `Brain - Frontal Cortex (BA9)`,
		spinal_cord = `Brain - Spinal cord (cervical c-1)`,
		`Heart - Atrial Appendage`, `Heart - Left Ventricle`,
		liver=`Liver`, ovary=`Ovary`, testis=`Testis`, lung=`Lung`,
		adrenal_gland = `Adrenal Gland`, pancreas=`Pancreas`,
		kidney = `Kidney - Cortex`,
		`Esophagus - Gastroesophageal Junction`, `Esophagus - Mucosa`, `Esophagus - Muscularis`,
		`Colon - Sigmoid`, `Colon - Transverse`,
		urinary_bladder = `Bladder`, prostate_gland = `Prostate`) %>%
	mutate(
		heart = rowMeans(dplyr::select(., `Heart - Atrial Appendage`, `Heart - Left Ventricle`)),
		esophagus = rowMeans(dplyr::select(., `Esophagus - Gastroesophageal Junction`, `Esophagus - Mucosa`, `Esophagus - Muscularis`)),
		colon = rowMeans(dplyr::select(., `Colon - Sigmoid`, `Colon - Transverse`))) %>%
	dplyr::select(-c(`Heart - Atrial Appendage`, `Heart - Left Ventricle`, `Esophagus - Gastroesophageal Junction`, `Esophagus - Mucosa`, `Esophagus - Muscularis`,`Colon - Sigmoid`, `Colon - Transverse`))



# load gene annotation
ensembl <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", version=76, GRCh=37)
ensembl_v76 <- getBM(attributes=c('ensembl_gene_id', 'hgnc_symbol'), filters = "", values = "", mart = ensembl) %>%
  as.tibble() %>%
  filter(hgnc_symbol != "")



# load Human Proteome Map expression per tissue
# re-processed and with more proteins than expression atlas
hpm_tissue <- read_tsv("./data/human_proteome_map/EnsembleGenePPKM_V3_Pandey.txt", na = c("x")) %>%
  inner_join(ensembl_v76, by=c("GENE" = "ensembl_gene_id")) %>%
  dplyr::select(-c(GENE)) %>%
  dplyr::select(Gene_Name = hgnc_symbol, frontal_cortex = frontalcortex, spinal_cord = spinalcord,
  liver, ovary, testis, lung, adrenal_gland, pancreas, kidney, urinary_bladder = urinarybladder,
  prostate_gland = prostate, heart, esophagus, colon)




# melt datasets
gtex_tissue <- gtex_tissue %>%
  gather(key="Tissue", value="RNA", -Gene_Name) %>%
  #mutate(RNA = log2(RNA + 1)) %>%
  group_by(Gene_Name, Tissue) %>%
  summarise(RNA = mean(RNA)) %>%
  group_by(Gene_Name) %>%
  filter( sum(RNA != 0) >= 10 ) %>%
  mutate(RNA = scale(RNA)[,1]) %>%
  ungroup()



hpm_tissue <- hpm_tissue %>%
  gather(key="Tissue", value="Protein", -Gene_Name) %>%
  #mutate(Protein = log2(Protein)) %>%
  group_by(Gene_Name, Tissue) %>%
  summarise(Protein = mean(Protein, na.rm=T)) %>%
  group_by(Gene_Name) %>%
  filter( sum(!is.na(Protein)) >= 10 ) %>%
  mutate(Protein = scale(Protein)[,1]) %>%
  ungroup()




# quantile-quantile normalization
gtex_tissue <- as.data.frame(spread(gtex_tissue, key = "Tissue", value = "RNA"))
rownames(gtex_tissue) <- gtex_tissue$Gene_Name
gtex_tissue <- gtex_tissue[, -c(1)]

gtex_tissue <- normalizeQuantiles( as.matrix(gtex_tissue) )
gtex_tissue <- cbind(data.frame(Gene_Name=rownames(gtex_tissue), stringsAsFactors=F), as.data.frame(gtex_tissue))


gtex_tissue <- gtex_tissue %>%
  as.tibble() %>%
  gather(key = "Tissue", value = "RNA", -Gene_Name)




hpm_tissue <- as.data.frame(spread(hpm_tissue, key = "Tissue", value = "Protein"))
rownames(hpm_tissue) <- hpm_tissue$Gene_Name
hpm_tissue <- hpm_tissue[, -c(1)]

hpm_tissue <- normalizeQuantiles( as.matrix(hpm_tissue) )
hpm_tissue <- cbind(data.frame(Gene_Name=rownames(hpm_tissue), stringsAsFactors=F), as.data.frame(hpm_tissue))


hpm_tissue <- hpm_tissue %>%
  as.tibble() %>%
  gather(key = "Tissue", value = "Protein", -Gene_Name)



# boxplot of RNA distribution
gtex_tissue_boxpl <- gtex_tissue %>%
#mutate(RNA = log2(RNA + 1)) %>%
ggplot(mapping = aes(x = Tissue, y = RNA, fill = Tissue)) +
  geom_boxplot(outlier.size = 0.5) +
  scale_fill_discrete(guide=F) +
  theme(axis.title = element_text(colour="black", size=15),
    axis.text.y = element_text(colour="black", size=13),
    axis.text.x = element_text(colour="black", size=13, angle=45, hjust=1)) +
  labs(x = "Tissue", y = "GTEx median RPKM per gene (zscore)")
ggsave(filename="gtex_rna_expression2.png", plot=gtex_tissue_boxpl, path = "./plots/pre_processing")
unlink("gtex_rna_expression2.png")



# boxplot of Protein distribution
# HPM protein atlas (reprocessed)
hpm_tissue_boxpl <- hpm_tissue %>%
  ggplot(mapping = aes(x = Tissue, y = Protein, fill = Tissue)) +
  geom_boxplot(outlier.size = 0.5) +
  scale_fill_discrete(guide=F) +
  theme(axis.title = element_text(colour="black", size=15),
    axis.text.y = element_text(colour="black", size=13),
    axis.text.x = element_text(colour="black", size=13, angle=45, hjust=1)) +
  labs(x = "Tissue", y = "HPM expression atlas reprocessed (zscore)")
ggsave(filename="hpm_atl_reprocessed_protein_expression.png", plot=hpm_tissue_boxpl, path = "./plots/pre_processing")
unlink("hpm_atl_reprocessed_protein_expression.png")




tissue_expression_gtex_hpm <- inner_join(gtex_tissue, hpm_tissue, by=c("Gene_Name", "Tissue"))




# boxplot of RNA and Protein distribution
tissue_protein_rna_boxpl <- tissue_expression_gtex_hpm %>%
  #mutate(Protein = log2(Protein), RNA = log2(RNA)) %>%
	gather(key = "var", value="value", -c(Gene_Name, Tissue)) %>%
	ggplot(mapping = aes(x = Tissue, y = value, fill = var)) +
	geom_boxplot(outlier.size = 0.5) +
	scale_fill_discrete(name="") +
    theme(axis.title = element_text(colour="black", size=15),
    	axis.text.y = element_text(colour="black", size=13),
    	axis.text.x = element_text(colour="black", size=13, angle=45, hjust=1)) +
    labs(x = "Tissue", y = "Protein and RNA expression (zscore)")
ggsave(filename="gtex_hpm_reprocessed_expression_distribution.png", plot=tissue_protein_rna_boxpl, path = "./plots/pre_processing")
unlink("gtex_hpm_reprocessed_expression_distribution.png")



# hpm correlations with GTEx
# across genes
hpm_gtex_cor <- tissue_expression_gtex_hpm %>%
	group_by(Tissue) %>%
	#summarise(counts = n())
	do(broom::tidy(cor.test(.$Protein, .$RNA, method = "pearson"))) %>%
	ungroup()


hpm_gtex_cor_bar <- hpm_gtex_cor %>%
	ggplot(mapping = aes(x = Tissue, y = estimate)) +
		geom_bar(stat = "identity", fill="white", color="black") +
		geom_hline(yintercept = 0.5, colour = "black", linetype="dashed") +
    geom_hline(yintercept = 0.4, colour = "black", linetype="dashed") +
		theme(axis.title = element_text(colour="black", size=15),
    		axis.text.y = element_text(colour="black", size=13),
    		axis.text.x = element_text(colour="black", size=13, angle=45, hjust=1)) +
    	labs(x = "Tissue", y = "Pearson correlation", title = paste("Correlation of Expression Atlas Human Proteome re-processed with GTEx", paste(length(unique(tissue_expression_gtex_hpm$Gene_Name)), "genes", sep=" "), sep="\n"))
ggsave(filename="hpm_reprocessed_gtex_correlation.png", plot=hpm_gtex_cor_bar, path = "./plots/pre_processing")
unlink("hpm_reprocessed_gtex_correlation.png")



# hpm correlations with GTEx
# across tissues
hpm_gtex_cor2 <- tissue_expression_gtex_hpm %>%
	group_by(Gene_Name) %>%
	#summarise(counts = n())
	do(broom::tidy(cor.test(.$Protein, .$RNA, method = "pearson"))) %>%
	ungroup()


hpm_gtex_cor_boxplot <- hpm_gtex_cor2 %>%
	ggplot(mapping = aes(x = factor(0), y = estimate, color = Correlation)) +
		geom_boxplot(fill="white", color="black") +
		theme(axis.title = element_text(colour="black", size=15),
    		axis.text.y = element_text(colour="black", size=13),
    		axis.text.x = element_text(colour="black", size=13)) +
    	labs(x = "", y = "Pearson correlation distribution", title = paste("Correlation of Expression Atlas Human Proteome re-processed with GTEx", paste(length(unique(tissue_expression_gtex_hpm$Gene_Name)), "genes", sep=" "), sep="\n"))
ggsave(filename="hpm_reprocessed_gtex_correlation2.png", plot=hpm_gtex_cor_boxplot, path = "./plots/pre_processing")
unlink("hpm_reprocessed_gtex_correlation2.png")



# hpm correlations with GTEx
# across tissues
# by attenuation state
hpm_gtex_cor_att <- tissue_expression_gtex_hpm %>%
	inner_join(protein_attenuation[, c("gene", "class")], by=c("Gene_Name" = "gene")) %>%
	group_by(class, Gene_Name) %>%
	#summarise(counts = n())
	do(broom::tidy(cor.test(.$Protein, .$RNA, method = "pearson"))) %>%
	ungroup() %>%
	mutate_at(vars(contains("class")), as.factor) %>%
  mutate_if(is.factor, fct_infreq)



hpm_gtex_cor_att_boxplot <- ggplot(data = hpm_gtex_cor_att, mapping = aes(x = class, y = estimate, color=class)) +
		geom_boxplot(fill="white") +
		scale_color_discrete(name="") +
		theme(axis.title = element_text(colour="black", size=15),
    		axis.text.y = element_text(colour="black", size=13),
    		axis.text.x = element_blank(),
    		axis.ticks.x = element_blank()) +
    	labs(x = "Attenuation state", y = "Pearson correlation distribution", title = paste("Human Proteome Map and GTEx correlation by attenuation state", paste(length(unique(tissue_expression_gtex_hpm$Gene_Name)), "genes", sep=" "), sep="\n"))
ggsave(filename="hpm_reprocessed_gtex_correlation_attenuation.png", plot=hpm_gtex_cor_att_boxplot, path = "./plots/pre_processing")
unlink("hpm_reprocessed_gtex_correlation_attenuation.png")









# merge with highly-significative protein associations
# FDR CNV < 1e-3
p_associations_hsignf <- protein_associations_rna_cnv %>%
	filter(fdr_cnv < 1e-3) %>%
	dplyr::select(controlling, controlled) %>%
	inner_join(tissue_expression_gtex_hpm, by=c("controlling" = "Gene_Name")) %>%
	dplyr::rename(Protein_controlling = Protein, RNA_controlling = RNA) %>%
	inner_join(tissue_expression_gtex_hpm, by=c("controlled" = "Gene_Name", "Tissue")) %>%
	dplyr::rename(Protein_controlled = Protein, RNA_controlled = RNA) %>%
	mutate(signf = "highly-significant")



# merge with significative protein associations
# 5e-2 >= FDR CNV >= 1e-3
p_associations_signf <- protein_associations_rna_cnv %>%
	filter(fdr_cnv >= 1e-3) %>%
	dplyr::select(controlling, controlled) %>%
  inner_join(tissue_expression_gtex_hpm, by=c("controlling" = "Gene_Name")) %>%
	dplyr::rename(Protein_controlling = Protein, RNA_controlling = RNA) %>%
	inner_join(tissue_expression_gtex_hpm, by=c("controlled" = "Gene_Name", "Tissue")) %>%
	dplyr::rename(Protein_controlled = Protein, RNA_controlled = RNA) %>%
	mutate(signf = "significant")



# merge with non-significative protein associations
# FDR CNV > 5e-2
# sample 516 associations
p_associations_nonsignf <- non_signf_protein_associations_cnv_rna %>%
	#dplyr::sample_n(size = 516) %>%
	dplyr::select(controlling, controlled) %>%
  inner_join(tissue_expression_gtex_hpm, by=c("controlling" = "Gene_Name")) %>%
	dplyr::rename(Protein_controlling = Protein, RNA_controlling = RNA) %>%
	inner_join(tissue_expression_gtex_hpm, by=c("controlled" = "Gene_Name", "Tissue")) %>%
	dplyr::rename(Protein_controlled = Protein, RNA_controlled = RNA) %>%
	mutate(signf = "non-significant")


p_associations <- bind_rows(p_associations_hsignf, p_associations_signf, p_associations_nonsignf)



# compute correlations

#corr(protein controlling, protein controlled)
p_associations_1 <- p_associations %>%
    group_by(signf, controlling, controlled) %>%
    do(broom::tidy(cor.test(.$Protein_controlling, .$Protein_controlled, method = "pearson"))) %>%
    ungroup() %>%
    mutate(correlation = "cor_Pcontrolling_Pcontrolled")

#corr(rna controlled, protein controlled)
p_associations_2 <- p_associations %>%
    group_by(signf, controlling, controlled) %>%
    do(broom::tidy(cor.test(.$RNA_controlled, .$Protein_controlled, method = "pearson"))) %>%
    ungroup() %>%
    mutate(correlation = "cor_RNAcontrolled_Pcontrolled")

#corr(rna controlling, rna controlled)
p_associations_3 <- p_associations %>%
    group_by(signf, controlling, controlled) %>%
    do(broom::tidy(cor.test(.$RNA_controlling, .$RNA_controlled, method = "pearson"))) %>%
    ungroup() %>%
    mutate(correlation = "cor_RNAcontrolling_RNAcontrolled")

#corr(rna controlling, protein controlled)
p_associations_4 <- p_associations %>%
  group_by(signf, controlling, controlled) %>%
  do(broom::tidy(cor.test(.$RNA_controlling, .$Protein_controlled, method = "pearson"))) %>%
  ungroup() %>%
  mutate(correlation = "cor_RNAcontrolling_Proteincontrolled")


# correlation for the association pairs
p_associations_5 <- bind_rows(p_associations_1, p_associations_2, p_associations_3) %>%
  dplyr::select(-c(method, statistic, p.value, parameter, conf.low, conf.high, alternative)) %>%
  spread(key=correlation, value=estimate) %>%
  filter(cor_RNAcontrolling_RNAcontrolled > 0 & cor_RNAcontrolling_RNAcontrolled < 0.4) %>%
  gather(key="correlation", value="estimate", -c(signf, controlling, controlled)) %>%
	mutate_at(vars(matches("signf|correlation")), as.factor) %>%
  mutate(signf = fct_relevel(signf, c("non-significant", "significant", "highly-significant"))) %>%
  mutate(correlation = fct_relevel(correlation, c("cor_Pcontrolling_Pcontrolled", "cor_RNAcontrolling_RNAcontrolled", "cor_RNAcontrolled_Pcontrolled"))) %>%
	ggplot(mapping = aes(x = correlation, y = estimate, colour = correlation)) +
	geom_boxplot(outlier.shape = NA) +
  #geom_jitter(width = 0.2, alpha = 0.4, colour="grey") +
  #stat_compare_means(comparisons = list( c(1, 2), c(1, 3), c(2, 3))) +
  facet_grid(. ~ signf, scales = "free", labeller=labeller(signf = c("highly-significant" = "Highly significant\nFDR CNV < 1e-3", "significant" = "Significant\n5e-2 >= FDR CNV >= 1e-3", "non-significant" = "Non-significant\nFDR CNV > 5e-2"))) +
  theme_classic() +
  scale_color_discrete(name="", label=c("Protein controlling vs Protein controlled", "RNA controlling vs RNA controlled", "RNA controlled vs Protein controlled")) +
  theme(axis.title.y=element_text(colour="black", size=15),
    axis.title.x=element_text(colour="black", size=15),
    axis.text.y=element_text(colour="black", size=12),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    strip.text = element_text(size=11),
    strip.background = element_blank(),
    plot.title=element_text(size=15),
    legend.text=element_text(size=10),
    panel.spacing = unit(2, "lines")) +
  labs(x = "Comparison", y = "Pearson correlation distribution", title="Correlation for the association pairs")
ggsave(filename="gtex_hpm_reprocessed_control_status_correlation.png", plot=p_associations_5, width = 10, path = "./plots/protein_associations_normal/")
unlink("gtex_hpm_reprocessed_control_status_correlation.png")


# correlation difference for the association pairs
p_associations_6 <- bind_rows(p_associations_1, p_associations_2) %>%
	dplyr::select(signf, correlation, controlling, controlled, estimate) %>%
	spread(key="correlation", value="estimate") %>%
	mutate(diff_cor = cor_Pcontrolling_Pcontrolled - cor_RNAcontrolled_Pcontrolled) %>%
	mutate_at(vars(contains("signf")), as.factor) %>%
    mutate(signf = fct_relevel(signf, c("highly-significant", "significant", "non-significant"))) %>%
	ggplot(mapping = aes(x = signf, y = diff_cor, colour = signf)) +
	geom_boxplot(outlier.shape = NA) +
    #geom_jitter(width = 0.2, alpha = 0.4, colour="grey") +
    stat_compare_means(comparisons = list( c(1, 2), c(1, 3), c(2, 3))) +
    theme_classic() +
    scale_color_discrete(name="", label=c("Highly significant\nFDR CNV < 1e-3", "Significant\n5e-2 >= FDR CNV >= 1e-3", "Non-significant\nFDR CNV > 5e-2")) +
    theme(axis.title.y=element_text(colour="black", size=15),
          axis.title.x=element_text(colour="black", size=15),
          axis.text.y=element_text(colour="black", size=12),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          strip.text = element_text(size=12),
          strip.background = element_blank(),
          plot.title=element_text(size=12),
          legend.text=element_text(size=10),
          legend.key.size = unit(2, "lines")) +
    labs(x = "Significance", y = "Difference of pearson correlation", title="Correlation difference for the association pairs\ncorr(controlling protein, controlled protein) - corr(controlled RNA, controlled protein)")
ggsave(filename="gtex_hpm_reprocessed_control_status_correlation_diff.png", plot=p_associations_6, width = 8, path = "./plots/protein_associations_normal/")
unlink("gtex_hpm_reprocessed_control_status_correlation_diff.png")


p_associations_7 <- bind_rows(p_associations_1, p_associations_2) %>%
	dplyr::select(signf, correlation, controlling, controlled, estimate) %>%
	spread(key="correlation", value="estimate") %>%
	mutate(diff_cor = cor_Pcontrolling_Pcontrolled - cor_RNAcontrolled_Pcontrolled) %>%
	group_by(controlling, controlled) %>%
	inner_join(bind_rows(protein_associations_rna_cnv, non_signf_protein_associations_cnv_rna) %>% dplyr::select(controlling, controlled, fdr_cnv) %>% mutate(fdr_cnv = -log10(fdr_cnv)), by=c("controlling", "controlled")) %>%
	ungroup() %>%
	mutate_at(vars(contains("signf")), as.factor) %>%
    mutate(signf = fct_relevel(signf, c("highly-significant", "significant", "non-significant"))) %>%
	ggplot(mapping = aes(x = fdr_cnv, y = diff_cor)) +
	geom_point() +
    geom_smooth(method=lm, color="blue", se=T) +
    #facet_grid(. ~ signf, scales = "free", labeller=labeller(signf = c("highly-significant" = "Highly significant\nFDR CNV < 1e-3", "significant" = "Significant\n5e-2 >= FDR CNV >= 1e-3", "non-significant" = "Non-significant\nFDR CNV > 5e-2"))) +
    stat_cor(label.y = 1.5) +
    theme_classic() +
    theme(axis.title.y=element_text(colour="black", size=15),
          axis.title.x=element_text(colour="black", size=15),
          axis.text=element_text(colour="black", size=12),
          strip.text = element_text(size=12),
          strip.background = element_blank(),
          plot.title=element_text(size=12)) +
    labs(x = "FDR CNV model (-log10)", y = "Difference of pearson correlation", title="Correlation between the correlation difference for the association pairs and the -log10 FDR")
ggsave(filename="gtex_hpm_reprocessed_control_status_cor_correlation_diff_fdr.png", plot=p_associations_7, width = 10, height=5, path = "./plots/protein_associations_normal/")
unlink("gtex_hpm_reprocessed_control_status_cor_correlation_diff_fdr.png")


p_associations_8 <- bind_rows(p_associations_1, p_associations_2) %>%
	dplyr::select(signf, correlation, controlling, controlled, estimate) %>%
	spread(key="correlation", value="estimate") %>%
	mutate(diff_cor = cor_Pcontrolling_Pcontrolled - cor_RNAcontrolled_Pcontrolled) %>%
	group_by(controlling, controlled) %>%
	inner_join(bind_rows(protein_associations_rna_cnv, non_signf_protein_associations_cnv_rna) %>% dplyr::select(controlling, controlled, fdr_cnv) %>% mutate(fdr_cnv = -log10(fdr_cnv)), by=c("controlling", "controlled")) %>%
	ungroup() %>%
	mutate_at(vars(contains("signf")), as.factor) %>%
    mutate(signf = fct_relevel(signf, c("highly-significant", "significant", "non-significant"))) %>%
	ggplot(mapping = aes(x = fdr_cnv, y = cor_Pcontrolling_Pcontrolled)) +
	geom_point() +
    geom_smooth(method=lm, color="blue", se=T) +
    #facet_grid(. ~ signf, scales = "free", labeller=labeller(signf = c("highly-significant" = "Highly significant\nFDR CNV < 1e-3", "significant" = "Significant\n5e-2 >= FDR CNV >= 1e-3", "non-significant" = "Non-significant\nFDR CNV > 5e-2"))) +
    stat_cor(label.y = 1.5) +
    theme_classic() +
    theme(axis.title.y=element_text(colour="black", size=13),
          axis.title.x=element_text(colour="black", size=15),
          axis.text=element_text(colour="black", size=12),
          strip.text = element_text(size=12),
          strip.background = element_blank(),
          plot.title=element_text(size=12)) +
    labs(x = "FDR CNV model (-log10)", y = "Correlation between controlling and controlled protein", title="Correlation between the correlation of the controlling protein and controlled protein and the -log10 FDR")
ggsave(filename="gtex_hpm_reprocessed_control_status_cor_correlation_controlling_controlled_fdr.png", plot=p_associations_8, width = 10, height=5, path = "./plots/protein_associations_normal/")
unlink("gtex_hpm_reprocessed_control_status_cor_correlation_controlling_controlled_fdr.png")





save(list=ls(), file="./r_workspaces/normal_samples_gtex_hpm_atl_reprocessed.RData")
