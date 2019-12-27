# Abel Sousa 2018

# PhD project
# Study of Control of Protein Abundance using Multi-omics Data from Cancer Samples

# -- Analysis of protein associations using normal samples of GTEx and Human Proteome Map

# -- Correlations of protein/RNA expression across normal tissues




suppressMessages(library(tidyverse))
suppressMessages(library(ggpubr))
suppressMessages(library(limma))
suppressMessages(library(viridis))

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



# load Human Proteome Map expression per tissue
# downloaded from HPM website
hpm_tissue1 <- read_csv("./data/human_proteome_map/HPM_gene_level_epxression_matrix_Kim_et_al_052914.csv") %>%
	dplyr::select(matches("Gene|Adult")) %>%
	dplyr::rename_at(vars(contains("Adult")), funs( stringr::str_replace_all(., " ", "_") ) ) %>%
	dplyr::rename_at(vars(contains("Adult")), funs( stringr::str_replace_all(., "Adult_", "") ) ) %>%
  dplyr::rename(
  	Gene_Name = Gene, frontal_cortex = Frontal_Cortex, spinal_cord = Spinal_Cord, retina = Retina,
  	heart = Heart, liver = Liver, ovary = Ovary, testis = Testis, lung = Lung, adrenal_gland = Adrenal,
  	gall_bladder = Gallbladder, pancreas = Pancreas, kidney = Kidney, esophagus = Esophagus, colon = Colon,
  	rectum = Rectum, urinary_bladder = Urinary_Bladder, prostate_gland = Prostate)



# load Human Proteome Map expression per tissue
# re-processed by expression atlas
hpm_tissue2 <- read_tsv("./data/human_proteome_map/E-PROT-1-query-results.tsv", skip=4) %>%
	dplyr::select(-`Gene ID`) %>%
	dplyr::rename_all(funs ( stringr::str_replace_all(., " ", "_") ) )




# melt datasets
gtex_tissue <- gtex_tissue %>%
  gather(key="Tissue", value="RNA", -Gene_Name)

hpm_tissue1 <- hpm_tissue1 %>%
  gather(key="Tissue", value="Protein_hpm", -Gene_Name)

hpm_tissue2 <- hpm_tissue2 %>%
  gather(key="Tissue", value="Protein_atl", -Gene_Name)

tissues <- gtex_tissue$Tissue %>% intersect(., hpm_tissue1$Tissue) %>% intersect(., hpm_tissue2$Tissue)
#tissues <- tissues[ tissues != "testis"]
print(tissues)



# filter datasets
gtex_tissue <- gtex_tissue %>%
  group_by(Tissue) %>%
  filter( Tissue %in% tissues ) %>%
  ungroup() %>%
  #mutate(RNA = log2(RNA + 1)) %>%
  group_by(Gene_Name, Tissue) %>%
  summarise(RNA = mean(RNA)) %>%
  group_by(Gene_Name) %>%
  filter( sum(RNA != 0) >= 10 ) %>%
  mutate(RNA = scale(RNA)[,1]) %>%
  ungroup()


hpm_tissue1 <- hpm_tissue1 %>%
  filter( Gene_Name %in% unique(hpm_tissue2$Gene_Name) ) %>%
  group_by(Tissue) %>%
  filter( Tissue %in% tissues ) %>%
  ungroup() %>%
  #mutate(Protein_hpm = log2(Protein_hpm + 1)) %>%
  group_by(Gene_Name) %>%
  filter( sum(Protein_hpm != 0) >= 10 ) %>%
  mutate(Protein_hpm = scale(Protein_hpm)[,1]) %>%
  ungroup()


hpm_tissue2 <- hpm_tissue2 %>%
  group_by(Tissue) %>%
  filter( Tissue %in% tissues ) %>%
  ungroup() %>%
  #mutate(Protein_atl = log2(Protein_atl)) %>%
  group_by(Gene_Name, Tissue) %>%
  summarise(Protein_atl = mean(Protein_atl, na.rm = T)) %>%
  group_by(Gene_Name) %>%
  filter( sum(!is.na(Protein_atl)) >= 5 ) %>%
  mutate(Protein_atl = scale(Protein_atl)[,1]) %>%
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




hpm_tissue1 <- as.data.frame(spread(hpm_tissue1, key = "Tissue", value = "Protein_hpm"))
rownames(hpm_tissue1) <- hpm_tissue1$Gene_Name
hpm_tissue1 <- hpm_tissue1[, -c(1)]

hpm_tissue1 <- normalizeQuantiles( as.matrix(hpm_tissue1) )
hpm_tissue1 <- cbind(data.frame(Gene_Name=rownames(hpm_tissue1), stringsAsFactors=F), as.data.frame(hpm_tissue1))


hpm_tissue1 <- hpm_tissue1 %>%
  as.tibble() %>%
  gather(key = "Tissue", value = "Protein_hpm", -Gene_Name)



hpm_tissue2 <- as.data.frame(spread(hpm_tissue2, key = "Tissue", value = "Protein_atl"))
rownames(hpm_tissue2) <- hpm_tissue2$Gene_Name
hpm_tissue2 <- hpm_tissue2[, -c(1)]

hpm_tissue2 <- normalizeQuantiles( as.matrix(hpm_tissue2) )
hpm_tissue2 <- cbind(data.frame(Gene_Name=rownames(hpm_tissue2), stringsAsFactors=F), as.data.frame(hpm_tissue2))


hpm_tissue2 <- hpm_tissue2 %>%
  as.tibble() %>%
  gather(key = "Tissue", value = "Protein_atl", -Gene_Name)




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
ggsave(filename="gtex_rna_expression.png", plot=gtex_tissue_boxpl, path = "./plots/pre_processing")
unlink("gtex_rna_expression.png")



# boxplot of Protein distribution
# HPM
hpm_tissue1_boxpl <- hpm_tissue1 %>%
  #mutate(Protein_hpm = log2(Protein_hpm + 1)) %>%
  ggplot(mapping = aes(x = Tissue, y = Protein_hpm, fill = Tissue)) +
  geom_boxplot(outlier.size = 0.5) +
  scale_fill_discrete(guide=F) +
    theme(axis.title = element_text(colour="black", size=15),
      axis.text.y = element_text(colour="black", size=13),
      axis.text.x = element_text(colour="black", size=13, angle=45, hjust=1)) +
    labs(x = "Tissue", y = "HPM median counts per gene (zscore)")
ggsave(filename="hpm_protein_expression.png", plot=hpm_tissue1_boxpl, path = "./plots/pre_processing")
unlink("hpm_protein_expression.png")



# boxplot of Protein distribution
# HPM protein atlas
hpm_tissue2_boxpl <- hpm_tissue2 %>%
  #mutate(Protein_atl = log2(Protein_atl)) %>%
  ggplot(mapping = aes(x = Tissue, y = Protein_atl, fill = Tissue)) +
  geom_boxplot(outlier.size = 0.5) +
  scale_fill_discrete(guide=F) +
    theme(axis.title = element_text(colour="black", size=15),
      axis.text.y = element_text(colour="black", size=13),
      axis.text.x = element_text(colour="black", size=13, angle=45, hjust=1)) +
    labs(x = "Tissue", y = "HPM expression atlas (zscore)")
ggsave(filename="hpm_atl_protein_expression.png", plot=hpm_tissue2_boxpl, path = "./plots/pre_processing")
unlink("hpm_atl_protein_expression.png")



# hpm correlations
# across genes
hpm_cor <- inner_join(hpm_tissue1, hpm_tissue2, by=c("Gene_Name", "Tissue")) %>%
	group_by(Tissue) %>%
	#summarise(counts = n())
	do(broom::tidy(cor.test(.$Protein_hpm, .$Protein_atl, method = "pearson"))) %>%
	ggplot(mapping = aes(x = Tissue, y = estimate, label = p.value)) +
		geom_bar(stat = "identity", fill="white", colour="grey") +
		geom_text() +
		geom_hline(yintercept = 0.5, colour = "black", linetype="dashed") +
		theme(axis.title = element_text(colour="black", size=15),
	    	axis.text.y = element_text(colour="black", size=13),
	    	axis.text.x = element_text(colour="black", size=13, angle=45, hjust=1)) +
    	labs(x = "Tissue", y = "Pearson correlation", title = "Correlation between Expression Atlas and Human Proteome Map\n 3858 genes")
ggsave(filename="hpm_correlation1.png", plot=hpm_cor, path = "./plots/pre_processing")
unlink("hpm_correlation1.png")


# hpm correlations
# across tissues
hpm_cor2 <- inner_join(hpm_tissue1, hpm_tissue2, by=c("Gene_Name", "Tissue")) %>%
	group_by(Gene_Name) %>%
	#summarise(counts = n())
	do(broom::tidy(cor.test(.$Protein_hpm, .$Protein_atl, method = "pearson"))) %>%
	ggplot(mapping = aes(x = factor(0), y = estimate)) +
		geom_boxplot(fill="white") +
		scale_x_discrete(label = "Correlation HPM vs Expression Atlas") +
		theme(axis.title = element_text(colour="black", size=15),
	    	axis.text.y = element_text(colour="black", size=13),
	    	axis.text.x = element_text(colour="black", size=13)) +
    	labs(x = "", y = "Pearson correlation", title = "Correlation between Expression Atlas and Human Proteome Map\n(3856 genes across 14 tissues)")
ggsave(filename="hpm_correlation2.png", plot=hpm_cor2, path = "./plots/pre_processing")
unlink("hpm_correlation2.png")


# hpm correlations with GTEx
# across genes
hpm_gtex_cor1 <- inner_join(hpm_tissue1, hpm_tissue2, by=c("Gene_Name", "Tissue")) %>%
	inner_join(gtex_tissue, by=c("Gene_Name", "Tissue")) %>%
	group_by(Tissue) %>%
	#summarise(counts = n())
	do(broom::tidy(cor.test(.$Protein_hpm, .$RNA, method = "pearson"))) %>%
	ungroup()

hpm_gtex_cor2 <- inner_join(hpm_tissue1, hpm_tissue2, by=c("Gene_Name", "Tissue")) %>%
	inner_join(gtex_tissue, by=c("Gene_Name", "Tissue")) %>%
	group_by(Tissue) %>%
	#summarise(counts = n())
	do(broom::tidy(cor.test(.$Protein_atl, .$RNA, method = "pearson"))) %>%
	ungroup()


hpm_gtex_cor_bar <- inner_join(hpm_gtex_cor1[, c("Tissue", "estimate", "p.value")], hpm_gtex_cor2[, c("Tissue", "estimate", "p.value")], by = c("Tissue")) %>%
	dplyr::rename(cor_hpm_gtex = estimate.x, pvalue_hpm_gtex = p.value.x, cor_atl_gtex = estimate.y, pvalue_atl_gtex = p.value.y) %>%
	gather(key="Correlation", value="Value", -c(Tissue, pvalue_hpm_gtex, pvalue_atl_gtex)) %>%
	ggplot(mapping = aes(x = Tissue, y = Value, color = Correlation)) +
		geom_bar(stat = "identity", position = "dodge", fill="white") +
		geom_hline(yintercept = 0.5, colour = "black", linetype="dashed") +
		theme(axis.title = element_text(colour="black", size=15),
    		axis.text.y = element_text(colour="black", size=13),
    		axis.text.x = element_text(colour="black", size=13, angle=45, hjust=1)) +
    	labs(x = "Tissue", y = "Pearson correlation", title = "Correlation of Expression Atlas and Human Proteome Map with GTEx")
ggsave(filename="hpm_gtex_correlation.png", plot=hpm_gtex_cor_bar, path = "./plots/pre_processing")
unlink("hpm_gtex_correlation.png")



# hpm correlations with GTEx
# across tissues
hpm_gtex_cor3 <- inner_join(hpm_tissue1, hpm_tissue2, by=c("Gene_Name", "Tissue")) %>%
	inner_join(gtex_tissue, by=c("Gene_Name", "Tissue")) %>%
	group_by(Gene_Name) %>%
	#summarise(counts = n())
	do(broom::tidy(cor.test(.$Protein_hpm, .$RNA, method = "pearson"))) %>%
	ungroup()

hpm_gtex_cor4 <- inner_join(hpm_tissue1, hpm_tissue2, by=c("Gene_Name", "Tissue")) %>%
	inner_join(gtex_tissue, by=c("Gene_Name", "Tissue")) %>%
	group_by(Gene_Name) %>%
	#summarise(counts = n())
	do(broom::tidy(cor.test(.$Protein_atl, .$RNA, method = "pearson"))) %>%
	ungroup()


hpm_gtex_cor_boxplot <- inner_join(hpm_gtex_cor3[, c("Gene_Name", "estimate")], hpm_gtex_cor4[, c("Gene_Name", "estimate")], by = c("Gene_Name")) %>%
	dplyr::rename(cor_hpm_gtex = estimate.x, cor_atl_gtex = estimate.y) %>%
	gather(key="Correlation", value="Value", -Gene_Name) %>%
	ggplot(mapping = aes(x = Correlation, y = Value, color = Correlation)) +
		geom_boxplot(fill="white") +
		theme(axis.title = element_text(colour="black", size=15),
    		axis.text.y = element_text(colour="black", size=13),
    		axis.text.x = element_text(colour="black", size=13)) +
    	labs(x = "Comparison", y = "Pearson correlation distribution", title = "Correlation of Expression Atlas and Human Proteome Map with GTEx\n(3826 genes across 14 tissues)")
ggsave(filename="hpm_gtex_correlation2.png", plot=hpm_gtex_cor_boxplot, path = "./plots/pre_processing")
unlink("hpm_gtex_correlation2.png")



# hpm correlations with GTEx
# across tissues
# by attenuation state
hpm_gtex_cor_att <- inner_join(hpm_tissue1, gtex_tissue, by=c("Gene_Name", "Tissue")) %>%
	inner_join(protein_attenuation[, c("gene", "class")], by=c("Gene_Name" = "gene")) %>%
	group_by(class, Gene_Name) %>%
	#summarise(counts = n())
	do(broom::tidy(cor.test(.$Protein_hpm, .$RNA, method = "pearson"))) %>%
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
    	labs(x = "Attenuation state", y = "Pearson correlation distribution", title = "Human Proteome Map and GTEx correlation by attenuation state\n(4462 genes across 14 tissues)")
ggsave(filename="hpm_gtex_correlation_attenuation.png", plot=hpm_gtex_cor_att_boxplot, path = "./plots/pre_processing")
unlink("hpm_gtex_correlation_attenuation.png")








tissue_expression_gtex_hpm <- inner_join(hpm_tissue1, gtex_tissue, by=c("Gene_Name", "Tissue"))





# merge with highly-significative protein associations
# FDR CNV < 1e-3
p_associations_hsignf <- protein_associations_rna_cnv %>%
	filter(fdr_cnv < 1e-3) %>%
	dplyr::select(controlling, controlled) %>%
	inner_join(tissue_expression_gtex_hpm, by=c("controlling" = "Gene_Name")) %>%
	dplyr::rename(Protein_controlling = Protein_hpm, RNA_controlling = RNA) %>%
	inner_join(tissue_expression_gtex_hpm, by=c("controlled" = "Gene_Name", "Tissue")) %>%
	dplyr::rename(Protein_controlled = Protein_hpm, RNA_controlled = RNA) %>%
	mutate(signf = "highly-significant")



# merge with significative protein associations
# 5e-2 >= FDR CNV >= 1e-3
p_associations_signf <- protein_associations_rna_cnv %>%
	filter(fdr_cnv >= 1e-3) %>%
	dplyr::select(controlling, controlled) %>%
	inner_join(tissue_expression_gtex_hpm, by=c("controlling" = "Gene_Name")) %>%
	dplyr::rename(Protein_controlling = Protein_hpm, RNA_controlling = RNA) %>%
	inner_join(tissue_expression_gtex_hpm, by=c("controlled" = "Gene_Name", "Tissue")) %>%
	dplyr::rename(Protein_controlled = Protein_hpm, RNA_controlled = RNA) %>%
	mutate(signf = "significant")



# merge with non-significative protein associations
# FDR CNV > 5e-2
# sample 516 associations
p_associations_nonsignf <- non_signf_protein_associations_cnv_rna %>%
	#dplyr::sample_n(size = 516) %>%
	dplyr::select(controlling, controlled) %>%
	inner_join(tissue_expression_gtex_hpm, by=c("controlling" = "Gene_Name")) %>%
	dplyr::rename(Protein_controlling = Protein_hpm, RNA_controlling = RNA) %>%
	inner_join(tissue_expression_gtex_hpm, by=c("controlled" = "Gene_Name", "Tissue")) %>%
	dplyr::rename(Protein_controlled = Protein_hpm, RNA_controlled = RNA) %>%
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

#corr(rna controlling, protein controlling)
p_associations_4 <- p_associations %>%
  group_by(signf, controlling, controlled) %>%
  do(broom::tidy(cor.test(.$RNA_controlling, .$Protein_controlling, method = "pearson"))) %>%
  ungroup() %>%
  mutate(correlation = "cor_RNAcontrolling_Proteincontrolling")

#corr(rna controlled, protein controlling)
p_associations_5 <- p_associations %>%
  group_by(signf, controlling, controlled) %>%
  do(broom::tidy(cor.test(.$RNA_controlled, .$Protein_controlling, method = "pearson"))) %>%
  ungroup() %>%
  mutate(correlation = "cor_RNAcontrolled_Proteincontrolling")


# correlation for the association pairs
gtex_hpm_control_status_cor <- bind_rows(p_associations_1, p_associations_2, p_associations_3) %>%
  dplyr::select(signf, correlation, controlling, controlled, estimate) %>%
  spread(key=correlation, value=estimate) %>%
  filter(cor_RNAcontrolling_RNAcontrolled > 0 & cor_RNAcontrolling_RNAcontrolled < 0.4) %>%
  gather(key="correlation", value="estimate", -c(signf, controlling, controlled)) %>%
	mutate_at(vars(matches("signf|correlation")), as.factor) %>%
  mutate(signf = fct_relevel(signf, c("non-significant", "significant", "highly-significant"))) %>%
  mutate(correlation = fct_relevel(correlation, c("cor_Pcontrolling_Pcontrolled", "cor_RNAcontrolling_RNAcontrolled", "cor_RNAcontrolled_Pcontrolled"))) %>%
	ggplot(mapping = aes(x = correlation, y = estimate, fill = correlation)) +
	geom_boxplot() +
  facet_wrap( ~ signf, strip.position = "bottom", labeller=labeller(signf = c("highly-significant" = "Highly-significant", "significant" = "Significant", "non-significant" = "Non-significant"))) +
  theme_classic() +
  theme(axis.title.y=element_text(colour="black", size=15),
    axis.title.x=element_blank(),
    axis.text.y=element_text(colour="black", size=12),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    strip.text = element_text(size=14),
    strip.background = element_blank(),
    plot.title=element_blank(),
    legend.text=element_text(size=12),
    legend.title = element_blank(),
    #legend.box.background = element_rect(colour = "black"),
    legend.position=c(0.65, 0.12)) +
  #scale_fill_discrete(name="Correlation", label=c("corr(controlling, controlled)\nProtein vs Protein", "corr(controlling, controlled)\nRNA vs RNA", "corr(controlled, controlled)\nProtein vs RNA")) +
  scale_fill_viridis(discrete=TRUE, option = "E", alpha=0.8, label=c("corr(controlling protein, controlled protein)", "corr(controlling RNA, controlled RNA)", "corr(controlled protein, controlled RNA)")) +
  scale_y_continuous(name = "Pearson r")
  #guides(fill = guide_legend(title.position = "top"))
ggsave(filename="gtex_hpm_control_status_correlation.png", plot=gtex_hpm_control_status_cor, height=5, width=6, path = "./plots/protein_associations_normal/")
ggsave(filename="gtex_hpm_control_status_correlation.pdf", plot=gtex_hpm_control_status_cor, height=5, width=6, path = "./plots/protein_associations_normal/")
unlink("gtex_hpm_control_status_correlation.png")
unlink("gtex_hpm_control_status_correlation.pdf")


gtex_hpm_control_status_cor2 <- bind_rows(p_associations_1, p_associations_2, p_associations_3, p_associations_4, p_associations_5) %>%
  dplyr::select(-c(method, statistic, p.value, parameter, conf.low, conf.high, alternative)) %>%
  spread(key=correlation, value=estimate) %>%
  gather(key="correlation", value="estimate", -c(signf, controlling, controlled)) %>%
	mutate_at(vars(matches("signf|correlation")), as.factor) %>%
  mutate(signf = fct_relevel(signf, c("non-significant", "significant", "highly-significant"))) %>%
  #mutate(correlation = fct_relevel(correlation, c("cor_Pcontrolling_Pcontrolled", "cor_RNAcontrolling_RNAcontrolled", "cor_RNAcontrolled_Pcontrolled"))) %>%
	ggplot(mapping = aes(x = correlation, y = estimate, fill = correlation)) +
	geom_boxplot() +
  #geom_jitter(width = 0.1, alpha = 0.1, color="black") +
  #stat_compare_means(comparisons = list( c(1, 2), c(1, 3), c(2, 3))) +
  facet_wrap( ~ signf, strip.position = "bottom", labeller=labeller(signf = c("highly-significant" = "Highly-significant", "significant" = "Significant", "non-significant" = "Non-significant"))) +
  theme_classic() +
  theme(axis.title.y=element_text(colour="black", size=15),
    axis.title.x=element_blank(),
    axis.text.y=element_text(colour="black", size=12),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    strip.text = element_text(size=14),
    strip.background = element_blank(),
    plot.title=element_blank(),
    legend.text=element_text(size=10),
    legend.title=element_text(size=12),
    legend.position="right") +
  #scale_fill_discrete(name="Correlation", label=c("corr(controlling, controlled)\nProtein vs Protein", "corr(controlling, controlled)\nRNA vs RNA", "corr(controlled, controlled)\nProtein vs RNA")) +
  #scale_fill_brewer(name="Correlation", label=c("controlling protein vs.\ncontrolled protein", "controlling RNA vs.\ncontrolled RNA", "controlled protein vs.\ncontrolled RNA"), type="qual", palette="Set2") +
  scale_y_continuous(name = "Pearson r")
ggsave(filename="gtex_hpm_control_status_correlation2.png", plot=gtex_hpm_control_status_cor2, height=5, width=10, path = "./plots/protein_associations_normal/")
unlink("gtex_hpm_control_status_correlation2.png")



# correlation difference for the association pairs
gtex_hpm_control_status_cor_diff <- bind_rows(p_associations_1, p_associations_2) %>%
	dplyr::select(signf, correlation, controlling, controlled, estimate) %>%
	spread(key="correlation", value="estimate") %>%
	mutate(diff_cor = cor_Pcontrolling_Pcontrolled - cor_RNAcontrolled_Pcontrolled) %>%
	mutate_at(vars(contains("signf")), as.factor) %>%
  mutate(signf = fct_relevel(signf, c("non-significant", "significant", "highly-significant"))) %>%
	ggplot(mapping = aes(x = signf, y = diff_cor, fill = signf)) +
	geom_boxplot() +
    #geom_jitter(width = 0.1, alpha = 0.1, colour="black") +
    #stat_compare_means(comparisons = list( c(1, 2), c(1, 3), c(2, 3))) +
    theme_classic() +
    scale_fill_brewer(name="", label=c("Non-significant", "Significant", "Highly-significant"), type="seq", palette="Greys") +
    theme(axis.title.y=element_text(colour="black", size=15),
          axis.title.x=element_text(colour="black", size=15),
          axis.text.y=element_text(colour="black", size=12),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          plot.title=element_blank(),
          legend.text=element_text(size=10),
          legend.title=element_text(size=12),
          legend.key.size = unit(2, "lines")) +
    labs(x = "Significance", y = "Difference of Pearson r")
ggsave(filename="gtex_hpm_control_status_correlation_diff.png", plot=gtex_hpm_control_status_cor_diff, height=5, width = 4, path = "./plots/protein_associations_normal/")
unlink("gtex_hpm_control_status_correlation_diff.png")


cor_fdr_gtex_hpm_control_status_cor_diff <- bind_rows(p_associations_1, p_associations_2) %>%
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
ggsave(filename="gtex_hpm_control_status_cor_correlation_diff_fdr.png", plot=cor_fdr_gtex_hpm_control_status_cor_diff, width = 10, height=5, path = "./plots/protein_associations_normal/")
unlink("gtex_hpm_control_status_cor_correlation_diff_fdr.png")


cor_fdr_gtex_hpm_control_status_cor <- bind_rows(p_associations_1, p_associations_2) %>%
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
ggsave(filename="gtex_hpm_control_status_cor_correlation_controlling_controlled_fdr.png", plot=cor_fdr_gtex_hpm_control_status_cor, width = 10, height=5, path = "./plots/protein_associations_normal/")
unlink("gtex_hpm_control_status_cor_correlation_controlling_controlled_fdr.png")





save(list=ls(), file="./r_workspaces/normal_samples_gtex_hpm.RData")
