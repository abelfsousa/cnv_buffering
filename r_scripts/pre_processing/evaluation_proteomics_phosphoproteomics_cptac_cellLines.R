# Abel Sousa 2018

# PhD project
# Study of Control of Protein Abundance using Multi-omics Data from Cancer Samples


# -- Evaluation of proteomics and phosphoproteomics datasets




suppressMessages(library(tidyverse))
options(bitmapType = "cairo")
set.seed(123)



# load proteomics datasets
proteomics <- read_tsv("./files/proteomics_cptac_cellLines.txt") %>%
	gather(key="sample", value="log2FC", -gene) %>%
	mutate(data = if_else(grepl("TCGA", sample, fixed=T), "tcga", "cell_lines")) %>%
	group_by(gene) %>%
	filter( sum(!is.na(log2FC)) >= sum(data=="tcga")*0.5 & sum(!is.na(log2FC)) >= sum(data=="cell_lines")*0.5 ) %>%
	ungroup()


proteomicsQ <- read_tsv("./files/proteomicsQ_cptac_cellLines.txt") %>%
	gather(key="sample", value="log2FC", -gene) %>%
	mutate(data = if_else(grepl("TCGA", sample, fixed=T), "tcga", "cell_lines")) %>%
	group_by(gene) %>%
	filter( sum(!is.na(log2FC)) >= sum(data=="tcga")*0.5 & sum(!is.na(log2FC)) >= sum(data=="cell_lines")*0.5 ) %>%
	ungroup()



# load phosphoproteomics datasets
phosphoproteomics <- read_tsv("./files/phosphoproteomics_cptac_cellLines.txt") %>%
	gather(key="sample", value="log2FC", -phospho_site) %>%
	mutate(data = if_else(grepl("TCGA", sample, fixed=T), "tcga", "cell_lines")) %>%
	group_by(phospho_site) %>%
	filter( sum(!is.na(log2FC)) >= sum(data=="tcga")*0.5 & sum(!is.na(log2FC)) >= sum(data=="cell_lines")*0.5 ) %>%
	ungroup() %>%
	separate(phospho_site, c("gene", "phospho_site"), sep="_")


phosphoproteomicsQ <- read_tsv("./files/phosphoproteomicsQ_cptac_cellLines.txt") %>%
	gather(key="sample", value="log2FC", -phospho_site) %>%
	mutate(data = if_else(grepl("TCGA", sample, fixed=T), "tcga", "cell_lines")) %>%
	group_by(phospho_site) %>%
	filter( sum(!is.na(log2FC)) >= sum(data=="tcga")*0.5 & sum(!is.na(log2FC)) >= sum(data=="cell_lines")*0.5 ) %>%
	ungroup() %>%
	separate(phospho_site, c("gene", "phospho_site"), sep="_")







proteomics_phospho <- inner_join(proteomics, phosphoproteomics, by=c("gene", "sample")) %>%
	rename(log2FC.prot = log2FC.x, log2FC.phospho = log2FC.y) %>%
	dplyr::select(gene, phospho_site, sample, log2FC.prot, log2FC.phospho) %>%
	group_by(gene, phospho_site) %>%
	do(broom::tidy(cor.test(.$log2FC.prot, .$log2FC.phospho))) %>%
	ungroup()

proteomicsQ_phosphoQ <- inner_join(proteomicsQ, phosphoproteomicsQ, by=c("gene", "sample")) %>%
	rename(log2FC.prot = log2FC.x, log2FC.phospho = log2FC.y) %>%
	dplyr::select(gene, phospho_site, sample, log2FC.prot, log2FC.phospho) %>%
	group_by(gene, phospho_site) %>%
	do(broom::tidy(cor.test(.$log2FC.prot, .$log2FC.phospho))) %>%
	ungroup()


proteomics_phospho_comp <- inner_join(proteomics_phospho, proteomicsQ_phosphoQ, by=c("gene", "phospho_site")) %>%
	dplyr::select(gene, phospho_site, estimate.x, estimate.y, p.value.x, p.value.y) %>%
	rename(estimate.proteomics_phospho = estimate.x, estimate.proteomicsQ_phosphoQ = estimate.y, p.value.proteomics_phospho = p.value.x, p.value.proteomicsQ_phosphoQ = p.value.y) %>%
	mutate(p.value.proteomics_phospho = -log10(p.value.proteomics_phospho), p.value.proteomicsQ_phosphoQ = -log10(p.value.proteomicsQ_phosphoQ))


proteomics_phospho_comp2 <- proteomics_phospho_comp[, c("gene", "phospho_site", "estimate.proteomics_phospho", "p.value.proteomics_phospho")] %>%
	gather(key="stat", value="value", -c(gene, phospho_site)) %>%
	mutate(var = if_else(stat == "estimate.proteomics_phospho", "pearson_r", "p_value"))

proteomics_phospho_comp3 <- proteomics_phospho_comp[, c("gene", "phospho_site", "estimate.proteomicsQ_phosphoQ", "p.value.proteomicsQ_phosphoQ")] %>%
	gather(key="stat", value="value", -c(gene, phospho_site)) %>%
	mutate(var = if_else(stat == "estimate.proteomicsQ_phosphoQ", "pearson_r", "p_value"))


proteomics_phospho_comp4 <- rbind(proteomics_phospho_comp2, proteomics_phospho_comp3)


proteomics_phospho_comp_r_boxpl <- ggplot(data = proteomics_phospho_comp4[proteomics_phospho_comp4$var == "pearson_r", ], mapping = aes(x = stat, y = value, fill=stat)) +
    geom_boxplot() +
    theme(axis.title.y=element_text(colour="black", size=12),
          axis.title.x=element_text(colour="black", size=12),
          axis.text.y=element_text(colour="black", size=10),
          axis.text.x=element_blank()) +
    labs(x = "", y = "pearson r")
ggsave(filename="proteomics_phospho_comp_r_boxpl.png", path = "./plots/", plot=proteomics_phospho_comp_r_boxpl)
unlink("proteomics_phospho_comp_r_boxpl.png")



proteomics_phospho_comp_p_boxpl <- ggplot(data = proteomics_phospho_comp4[proteomics_phospho_comp4$var == "p_value", ], mapping = aes(x = stat, y = value, fill=stat)) +
    geom_boxplot() +
    theme(axis.title.y=element_text(colour="black", size=12),
          axis.title.x=element_text(colour="black", size=12),
          axis.text.y=element_text(colour="black", size=10),
          axis.text.x=element_blank()) +
    labs(x = "", y = "-log10 p_value")
ggsave(filename="proteomics_phospho_comp_p_boxpl.png", path = "./plots/", plot=proteomics_phospho_comp_p_boxpl)
unlink("proteomics_phospho_comp_p_boxpl.png")









save(list=ls(), file="./r_workspaces/evaluation_proteomics_phosphoproteomics_cptac_cellLines.RData")





