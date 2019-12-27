# Abel Sousa 2018

# PhD project
# Study of Control of Protein Abundance using Multi-omics Data from Cancer Samples




suppressMessages(library(tidyverse))
suppressMessages(library(RColorBrewer))
suppressMessages(library(ggpubr))
suppressMessages(library(grid))
suppressMessages(library(gridExtra))
suppressMessages(library(cowplot))
options(bitmapType = "cairo")



# import protein associations with CNV and RNA
# filtered by genomic co-localization
protein_associations_cnv_mrna <- read_tsv("./files/signf_protein_associations_cnv_mrna_filtered_close_controlling.txt") %>%
  mutate(pair = paste(controlling, controlled, sep="_")) %>%
  dplyr::select(pair, everything())



# import genome-wide protein associations with CNV and RNA
# filtered by genomic co-localization
protein_associations_cnv_mrna_gw <- read_tsv("./files/signf_protein_associations_cnv_mrna_filtered_close_controlling_genomeWide_fdr10.txt") %>%
  mutate(pair = paste(controlling, controlled, sep="_")) %>%
  dplyr::select(pair, everything())



# load common proteomics rna cnv metadata samples
common_samples <- read_tsv("./files/cptac_cellLines_samples_prot_rna_cnv.txt", col_names = FALSE)


# load metadata
metadata <- read_tsv("./files/metadata_cptac_cellLines.txt") %>%
  arrange(sample)


cnv <- read_tsv("./files/cnv_tcga_cellLines.txt") %>%
    gather(key="sample", value="cnv_gistic2", -gene) %>%
    arrange(sample) %>%
    filter(sample %in% common_samples$X1)


proteomics <- read_tsv("./files/proteomicsQ_cptac_cellLines.txt") %>%
    gather(key="sample", value="prot_log2FC", -gene) %>%
    arrange(sample) %>%
    filter(sample %in% common_samples$X1) %>%
    group_by(gene) %>%
    filter( sum(!is.na(prot_log2FC)) > n()*0.25 ) %>%
    mutate(prot_log2FC = scale(prot_log2FC)[,1]) %>%
    ungroup()


rna <- read_tsv("./files/rna_tcga_cellLines.txt") %>%
    gather(key="sample", value="rna_log2CPM", -gene) %>%
    arrange(sample) %>%
    filter(sample %in% common_samples$X1) %>%
    group_by(gene) %>%
    mutate(rna_log2CPM = scale(rna_log2CPM)[,1]) %>%
    ungroup()


linear_model <- function(data){
  data <- as.data.frame(data)
  rownames(data) <- data$sample

  covars <- data[, -c(1)]

  # remove nas
  covars <- na.exclude(covars)


  # remove covariates with less than two levels
  covars <- covars[, apply(covars, 2, function(x) length(levels(as.factor(x))) >= 2), drop=F]


  # fit the model
  lreg_f <- as.formula( paste( "proteomics_controlled", "~", paste(colnames(covars[, -match("proteomics_controlled", colnames(covars)), drop=F]), collapse="+") ) )
  lreg <- lm(lreg_f, data = covars)


  # model residuals
  resd <- as.tibble(data.frame(sample = names(residuals(lreg)), proteomics_controlled_ = unname(residuals(lreg)), stringsAsFactors = F))

  return(resd)

}


data <- protein_associations_cnv_mrna %>%
	dplyr::select(controlling, controlled) %>%
  inner_join(proteomics, by=c("controlled" = "gene")) %>%
  dplyr::rename(proteomics_controlled = prot_log2FC) %>%
  inner_join(rna, by=c("controlled" = "gene", "sample")) %>%
  dplyr::rename(rna_controlled = rna_log2CPM) %>%
  inner_join(metadata, by = "sample") %>%
  #regress-out rna from protein abundance
  nest(-c(controlling, controlled)) %>%
  mutate(proteomics_controlled_ = purrr::map(data, ~ linear_model(.))) %>%
  #rowwise() %>%
  #mutate(proteomics_controlled_ = list(linear_model(data))) %>%
  dplyr::select(-data) %>%
  unnest(proteomics_controlled_) %>%
  #inner_join(proteomics, by=c("controlled" = "gene", "sample")) %>%
  #dplyr::rename(proteomics_controlled = prot_log2FC) %>%
  #inner_join(rna, by=c("controlled" = "gene", "sample")) %>%
  #dplyr::rename(rna_controlled = rna_log2CPM) %>%
	inner_join(cnv, by=c("controlling" = "gene", "sample")) %>%
	dplyr::rename(cnv_controlling = cnv_gistic2) %>%
  inner_join(proteomics, by=c("controlling" = "gene", "sample")) %>%
  dplyr::rename(proteomics_controlling = prot_log2FC) %>%
	mutate(pair = paste(controlling, controlled, sep="_")) %>%
	dplyr::select(pair, controlling, controlled, sample, everything()) %>%
	mutate(cnv_controlling = as.factor(cnv_controlling)) %>%
  inner_join(metadata[, c("sample", "batch")], by = "sample") %>%
  mutate_at(vars(contains("batch")), as.factor) %>%
  mutate(batch = fct_relevel(batch, "LAW", "LPK", "TCGA_BRCA", "RMLT", "TCGA_COREAD", "TCGA_OV"))


data2 <- protein_associations_cnv_mrna_gw %>%
  dplyr::select(controlling, controlled) %>%
  inner_join(proteomics, by=c("controlled" = "gene")) %>%
  dplyr::rename(proteomics_controlled = prot_log2FC) %>%
  inner_join(rna, by=c("controlled" = "gene", "sample")) %>%
  dplyr::rename(rna_controlled = rna_log2CPM) %>%
  inner_join(metadata, by = "sample") %>%
  #regress-out rna from protein abundance
  nest(-c(controlling, controlled)) %>%
  mutate(proteomics_controlled_ = purrr::map(data, ~ linear_model(.))) %>%
  #rowwise() %>%
  #mutate(proteomics_controlled_ = list(linear_model(data))) %>%
  dplyr::select(-data) %>%
  unnest(proteomics_controlled_) %>%
  #inner_join(proteomics, by=c("controlled" = "gene", "sample")) %>%
  #dplyr::rename(proteomics_controlled = prot_log2FC) %>%
  #inner_join(rna, by=c("controlled" = "gene", "sample")) %>%
  #dplyr::rename(rna_controlled = rna_log2CPM) %>%
  inner_join(cnv, by=c("controlling" = "gene", "sample")) %>%
  dplyr::rename(cnv_controlling = cnv_gistic2) %>%
  inner_join(proteomics, by=c("controlling" = "gene", "sample")) %>%
  dplyr::rename(proteomics_controlling = prot_log2FC) %>%
  mutate(pair = paste(controlling, controlled, sep="_")) %>%
	dplyr::select(pair, controlling, controlled, sample, everything()) %>%
  mutate(cnv_controlling = as.factor(cnv_controlling)) %>%
  inner_join(metadata[, c("sample", "batch")], by = "sample") %>%
  mutate_at(vars(contains("batch")), as.factor) %>%
  mutate(batch = fct_relevel(batch, "LAW", "LPK", "TCGA_BRCA", "RMLT", "TCGA_COREAD", "TCGA_OV"))



#COG3 CNV -> COG2 Prot'
#data %>% filter(controlling == "COG3", controlled == "COG2") %>% group_by(controlling, controlled) %>% summarise(cor = cor.test(as.numeric(as.character(cnv_controlling)), proteomics_controlled_)$estimate)
plot1 <- ggplot(data = data %>% filter(controlling == "COG3", controlled == "COG2"), mapping = aes(x = cnv_controlling, y = proteomics_controlled_)) +
    geom_boxplot(outlier.shape=NA, aes(fill=cnv_controlling)) +
    geom_jitter(aes(fill=cnv_controlling), pch=21, width=0.1, alpha=0.5, color = "black") +
    annotate("text", label=paste("r = ", "0.37, ", "p < ", "3.79e-13"), x = 1.5, y = 2) +
    theme_classic() +
    theme(axis.title.y=element_text(colour="black", size=14),
          axis.title.x=element_text(colour="black", size=14),
          axis.text.y=element_text(colour="black", size=12),
          axis.text.x=element_text(colour="black", size=12),
          aspect.ratio=1) +
    scale_x_discrete(name = "COG3 CNV (gistic2)") +
    scale_y_continuous(name = "COG2 proteomics residuals (log2FC)") +
    scale_fill_brewer(type = "qual", palette = "GnBu", guide=F)
ggsave(filename="cog3_cnv_cog2_prot.png", plot=plot1, width=5, height=5, path = "./plots/protein_pairs_examples/")
unlink("cog3_cnv_cog2_prot.png")


#COG3 Prot -> COG2 Prot'
plot2 <- ggplot(data = data %>% filter(controlling == "COG3", controlled == "COG2"), mapping = aes(x = proteomics_controlling, y = proteomics_controlled_)) +
    geom_point(aes(color=batch)) +
    geom_smooth(method=lm , color="black", fill = "grey", se=T, size=0.5) +
    stat_cor() +
    theme_classic() +
    theme(axis.title.y=element_text(colour="black", size=14),
          axis.title.x=element_text(colour="black", size=14),
          axis.text.y=element_text(colour="black", size=12),
          axis.text.x=element_text(colour="black", size=12),
          aspect.ratio=1,
          legend.position="bottom") +
    scale_x_continuous(name = "COG3 proteomics (log2FC)") +
    scale_y_continuous(name = "COG2 proteomics residuals (log2FC)")
ggsave(filename="cog3_prot_cog2_prot.png", plot=plot2, width=5, height=5, path = "./plots/protein_pairs_examples/")
unlink("cog3_prot_cog2_prot.png")



#TRMT6 CNV -> TRMT61A Prot'
#data %>% filter(controlling == "TRMT6", controlled == "TRMT61A") %>% group_by(controlling, controlled) %>% summarise(cor = cor.test(as.numeric(as.character(cnv_controlling)), proteomics_controlled_)$estimate)
plot3 <- ggplot(data = data %>% filter(controlling == "TRMT6", controlled == "TRMT61A"), mapping = aes(x = cnv_controlling, y = proteomics_controlled_)) +
    geom_boxplot(outlier.shape=NA, aes(fill=cnv_controlling)) +
    geom_jitter(width=0.1, alpha=0.5, pch=21, aes(fill=cnv_controlling), color = "black") +
    annotate("text", label=paste("r = ", "0.45, ", "p < ", "1.09e-15"), x = 1.5, y = 2) +
    theme_classic() +
    theme(axis.title.y=element_text(colour="black", size=14),
      axis.title.x=element_text(colour="black", size=14),
      axis.text.y=element_text(colour="black", size=12),
      axis.text.x=element_text(colour="black", size=12),
      aspect.ratio = 1) +
    scale_x_discrete(name = "TRMT6 CNV (gistic2)") +
    scale_y_continuous(name = "TRMT61A proteomics residuals (log2FC)", limits=c(-2.5, 2.5)) +
    scale_fill_brewer(type = "qual", palette = "GnBu", guide=F)
ggsave(filename="trmt6_cnv_trmt61A_prot.png", plot=plot3, width=5, height=5, path = "./plots/protein_pairs_examples/")
unlink("trmt6_cnv_trmt61A_prot.png")



#TRMT6 Prot -> TRMT61A Prot'
plot4 <- ggplot(data = data %>% filter(controlling == "TRMT6", controlled == "TRMT61A"), mapping = aes(x = proteomics_controlling, y = proteomics_controlled_)) +
    geom_point(aes(color=batch)) +
    geom_smooth(method=lm , color="black", fill = "grey", se=T, size=0.5) +
    stat_cor() +
    theme_classic() +
    theme(axis.title.y=element_text(colour="black", size=14),
          axis.title.x=element_text(colour="black", size=14),
          axis.text.y=element_text(colour="black", size=12),
          axis.text.x=element_text(colour="black", size=12),
          aspect.ratio=1,
          legend.position="bottom") +
    scale_x_continuous(name = "TRMT6 proteomics (log2FC)") +
    scale_y_continuous(name = "TRMT61A proteomics residuals (log2FC)")
ggsave(filename="trmt6_prot_trmt61A_prot.png", plot=plot4, width=5, height=5, path = "./plots/protein_pairs_examples/")
unlink("trmt6_prot_trmt61A_prot.png")



#IDH3A CNV -> IDH3B Prot'
#data %>% filter(controlling == "IDH3A", controlled == "IDH3B") %>% group_by(controlling, controlled) %>% summarise(cor = cor.test(as.numeric(as.character(cnv_controlling)), proteomics_controlled_)$estimate)
plot5 <- ggplot(data = data %>% filter(controlling == "IDH3A", controlled == "IDH3B"), mapping = aes(x = cnv_controlling, y = proteomics_controlled_)) +
    geom_boxplot(outlier.shape=NA, aes(fill=cnv_controlling)) +
    geom_jitter(width=0.1, alpha=0.5, pch=21, aes(fill=cnv_controlling), color = "black") +
    annotate("text", label=paste("r = ", "0.39, ", "p < ", "1.57e-14"), x = 1.5, y = 2) +
    theme_classic() +
    theme(axis.title.y=element_text(colour="black", size=14),
          axis.title.x=element_text(colour="black", size=14),
          axis.text.y=element_text(colour="black", size=12),
          axis.text.x=element_text(colour="black", size=12),
          aspect.ratio = 1) +
    scale_x_discrete(name = "IDH3A CNV (gistic2)") +
    scale_y_continuous(name = "IDH3B proteomics residuals (log2FC)") +
    scale_fill_brewer(type = "qual", palette = "GnBu", guide=F)
ggsave(filename="idh3a_cnv_idh3b_prot.png", plot=plot5, width=5, height=5, path = "./plots/protein_pairs_examples/")
unlink("idh3a_cnv_idh3b_prot.png")



#IDH3A Prot -> IDH3B Prot'
plot6 <- ggplot(data = data %>% filter(controlling == "IDH3A", controlled == "IDH3B"), mapping = aes(x = proteomics_controlling, y = proteomics_controlled_)) +
    geom_point(aes(color=batch)) +
    geom_smooth(method=lm , color="black", fill = "grey", se=T, size=0.5) +
    stat_cor() +
    theme_classic() +
    theme(axis.title.y=element_text(colour="black", size=14),
          axis.title.x=element_text(colour="black", size=14),
          axis.text.y=element_text(colour="black", size=12),
          axis.text.x=element_text(colour="black", size=12),
          aspect.ratio=1,
          legend.position="bottom") +
    scale_x_continuous(name = "IDH3A proteomics (log2FC)") +
    scale_y_continuous(name = "IDH3B proteomics residuals (log2FC)")
ggsave(filename="idh3a_prot_idh3b_prot.png", plot=plot6, width=5, height=5, path = "./plots/protein_pairs_examples/")
unlink("idh3a_prot_idh3b_prot.png")


#TCP1 CNV -> CCT3/5/7/8 Prot'
plot7_stats <- data %>%
  filter(controlling == "TCP1" & controlled == "CCT3" | controlled == "CCT7" | controlled == "CCT5" | controlled == "CCT8") %>%
  group_by(controlling, controlled) %>%
  do(broom::tidy(cor.test(as.numeric(as.character(.$cnv_controlling)), .$proteomics_controlled_, method = "pearson"))) %>%
  ungroup() %>%
  dplyr::select(controlling, controlled, estimate, p.value) %>%
  mutate(text = paste("R = ", round(estimate, 2), ", p = ", signif(p.value, 3))) %>%
  mutate(x = rep(2.5, 4), y = rep(4, 4))

plot7 <- ggplot(data = data %>% filter(controlling == "TCP1" & controlled == "CCT3" | controlled == "CCT7" | controlled == "CCT5" | controlled == "CCT8"), mapping = aes(x = cnv_controlling, y = proteomics_controlled_)) +
    geom_boxplot(outlier.shape=NA, aes(fill=cnv_controlling)) +
    geom_jitter(width=0.1, alpha=0.5, pch=21, aes(fill=cnv_controlling), color = "black") +
    geom_text(data = plot7_stats, mapping = aes(x = x, y = y, label = text), size = 2.5) +
    facet_wrap( ~ controlled, scales = "free") +
    theme_classic() +
    theme(axis.title.y=element_text(colour="black", size=12),
          axis.title.x=element_text(colour="black", size=12),
          axis.text.y=element_text(colour="black", size=10),
          axis.text.x=element_text(colour="black", size=10),
          strip.text.x = element_text(size=12),
          strip.background = element_blank(),
          aspect.ratio=1) +
    scale_x_discrete(name = "TCP1 CNV") +
    scale_y_continuous(name = "Protein residuals (log2FC)") +
    scale_fill_brewer(type = "seq", palette = "YlOrBr", guide = F)
ggsave(filename="tcp1_cnv_cct3_5_7_8_prot.png", plot=plot7, width=4, height=4, path = "./plots/protein_pairs_examples/")
ggsave(filename="tcp1_cnv_cct3_5_7_8_prot.pdf", plot=plot7, width=4, height=4, path = "./plots/protein_pairs_examples/")
unlink("tcp1_cnv_cct3_5_7_8_prot.png")
unlink("tcp1_cnv_cct3_5_7_8_prot.pdf")



#TCP1 Prot -> CCT3/5/7/8 Prot'
plot8 <- ggplot(data = data %>% filter(controlling == "TCP1" & controlled == "CCT3" | controlled == "CCT7" | controlled == "CCT5" | controlled == "CCT8"), mapping = aes(x = proteomics_controlling, y = proteomics_controlled_)) +
    geom_point(aes(color=batch)) +
    geom_smooth(method=lm , color="black", fill = "grey", se=T, size=0.5) +
    stat_cor(size=2.5, label.x=-2, label.y=4) +
    facet_wrap( ~ controlled, scales = "free") +
    theme_classic() +
    theme(axis.title.y=element_text(colour="black", size=12),
          axis.title.x=element_text(colour="black", size=12),
          axis.text.y=element_text(colour="black", size=10),
          axis.text.x=element_text(colour="black", size=10),
          strip.text.x = element_text(size=12),
          strip.background = element_blank(),
          aspect.ratio=1,
          legend.position="bottom") +
    scale_color_brewer(type = "div", palette = "RdYlGn", labels=c("Lawrence et al", "Lapek et al", "TCGA BRCA", "Roumeliotis et al", "TCGA COREAD", "TCGA HGSC"), name="Samples", guide=F, direction=-1) +
    scale_x_continuous(name = "TCP1 protein (log2FC)") +
    scale_y_continuous(name = "Protein residuals (log2FC)")
ggsave(filename="tcp1_prot_cct3_5_7_8_prot.png", plot=plot8, width=4, height=4, path = "./plots/protein_pairs_examples/")
ggsave(filename="tcp1_prot_cct3_5_7_8_prot.pdf", plot=plot8, width=4, height=4, path = "./plots/protein_pairs_examples/")
unlink("tcp1_prot_cct3_5_7_8_prot.png")
unlink("tcp1_prot_cct3_5_7_8_prot.pdf")



#SMARCA2 CNV -> SMARCA4 Prot'
#data %>% filter(controlling == "SMARCA2", controlled == "SMARCA4") %>% group_by(controlling, controlled) %>% summarise(cor = cor.test(as.numeric(as.character(cnv_controlling)), proteomics_controlled_)$estimate)
plot9 <- ggplot(data = data %>% filter(controlling == "SMARCA2", controlled == "SMARCA4"), mapping = aes(x = cnv_controlling, y = proteomics_controlled_)) +
    geom_boxplot(outlier.shape=NA, aes(fill=cnv_controlling)) +
    geom_jitter(width=0.1, alpha=0.5, pch=21, aes(fill=cnv_controlling), color = "black") +
    annotate("text", label=paste("r = ", "-0.29, ", "p < ", "5e-08"), x = 1.5, y = 3) +
    theme_classic() +
    theme(axis.title.y=element_text(colour="black", size=14),
          axis.title.x=element_text(colour="black", size=14),
          axis.text.y=element_text(colour="black", size=12),
          axis.text.x=element_text(colour="black", size=12),
          aspect.ratio=1) +
    scale_x_discrete(name = "SMARCA2 CNV (gistic2)") +
    scale_y_continuous(name = "SMARCA4 proteomics residuals (log2FC)") +
    scale_fill_brewer(type = "qual", palette = "GnBu", guide=F)
ggsave(filename="smarca2_cnv_smarca4_prot.png", plot=plot9, width=5, height=5, path = "./plots/protein_pairs_examples/")
unlink("smarca2_cnv_smarca4_prot.png")



#SMARCA2 Prot -> SMARCA4 Prot'
plot10 <- ggplot(data = data %>% filter(controlling == "SMARCA2", controlled == "SMARCA4"), mapping = aes(x = proteomics_controlling, y = proteomics_controlled_)) +
    geom_point(aes(color=batch)) +
    geom_smooth(method=lm , color="black", fill = "grey", se=T, size=0.5) +
    stat_cor() +
    theme_classic() +
    theme(axis.title.y=element_text(colour="black", size=14),
          axis.title.x=element_text(colour="black", size=14),
          axis.text.y=element_text(colour="black", size=12),
          axis.text.x=element_text(colour="black", size=12),
          aspect.ratio=1,
          legend.position="bottom") +
    scale_x_continuous(name = "SMARCA2 proteomics (log2FC)") +
    scale_y_continuous(name = "SMARCA4 proteomics residuals (log2FC)")
ggsave(filename="smarca2_prot_smarca4_prot.png", plot=plot10, width=5, height=5, path = "./plots/protein_pairs_examples/")
unlink("smarca2_prot_smarca4_prot.png")




#EIF4A3 CNV -> FYTTD1 Prot'
#data %>% filter(controlling == "EIF4A3", controlled == "FYTTD1") %>% group_by(controlling, controlled) %>% summarise(cor = cor.test(as.numeric(as.character(cnv_controlling)), proteomics_controlled_)$estimate)
plot11 <- ggplot(data = data %>% filter(controlling == "EIF4A3", controlled == "FYTTD1"), mapping = aes(x = cnv_controlling, y = proteomics_controlled_)) +
    geom_boxplot(outlier.shape=NA, aes(fill=cnv_controlling)) +
    geom_jitter(width=0.1, alpha=0.5, pch=21, aes(fill=cnv_controlling), color = "black") +
    annotate("text", label=paste("r = ", "-0.35, ", "p < ", "1.42e-11"), x = 1.5, y = 3) +
    theme_classic() +
    theme(axis.title.y=element_text(colour="black", size=14),
          axis.title.x=element_text(colour="black", size=14),
          axis.text.y=element_text(colour="black", size=12),
          axis.text.x=element_text(colour="black", size=12),
          aspect.ratio=1) +
    scale_x_discrete(name = "EIF4A3 CNV (gistic2)") +
    scale_y_continuous(name = "FYTTD1 proteomics residuals (log2FC)") +
    scale_fill_brewer(type = "qual", palette = "GnBu", guide=F)
ggsave(filename="eif4a3_cnv_fyttd1_prot.png", plot=plot11, width=5, height=5, path = "./plots/protein_pairs_examples/")
unlink("eif4a3_cnv_fyttd1_prot.png")



#EIF4A3 Prot -> FYTTD1 Prot'
plot12 <- ggplot(data = data %>% filter(controlling == "EIF4A3", controlled == "FYTTD1"), mapping = aes(x = proteomics_controlling, y = proteomics_controlled_)) +
    geom_point(aes(color=batch)) +
    geom_smooth(method=lm , color="black", fill = "grey", se=T, size=0.5) +
    stat_cor() +
    theme_classic() +
    theme(axis.title.y=element_text(colour="black", size=14),
          axis.title.x=element_text(colour="black", size=14),
          axis.text.y=element_text(colour="black", size=12),
          axis.text.x=element_text(colour="black", size=12),
          aspect.ratio=1,
          legend.position="bottom") +
    scale_x_continuous(name = "EIF4A3 proteomics (log2FC)") +
    scale_y_continuous(name = "FYTTD1 proteomics residuals (log2FC)")
ggsave(filename="eif4a3_prot_fyttd1_prot.png", plot=plot12, width=5, height=5, path = "./plots/protein_pairs_examples/")
unlink("eif4a3_prot_fyttd1_prot.png")



#DDX5 CNV -> DDX17 Prot'
#data %>% filter(controlling == "DDX5", controlled == "DDX17") %>% group_by(controlling, controlled) %>% summarise(cor = cor.test(as.numeric(as.character(cnv_controlling)), proteomics_controlled_)$estimate)
plot13 <- ggplot(data = data %>% filter(controlling == "DDX5", controlled == "DDX17"), mapping = aes(x = cnv_controlling, y = proteomics_controlled_)) +
    geom_boxplot(outlier.shape=NA, aes(fill=cnv_controlling)) +
    geom_jitter(width=0.1, alpha=0.5, pch=21, aes(fill=cnv_controlling), color = "black") +
    annotate("text", label=paste("r = ", "-0.27, ", "p < ", "9.88e-08"), x = 1.5, y = 3) +
    theme_classic() +
    theme(axis.title.y=element_text(colour="black", size=14),
          axis.title.x=element_text(colour="black", size=14),
          axis.text.y=element_text(colour="black", size=12),
          axis.text.x=element_text(colour="black", size=12),
          aspect.ratio=1) +
    scale_x_discrete(name = "DDX5 CNV (gistic2)") +
    scale_y_continuous(name = "DDX17 proteomics residuals (log2FC)") +
    scale_fill_brewer(type = "qual", palette = "GnBu", guide=F)
ggsave(filename="ddx5_cnv_ddx17_prot.png", plot=plot13, width=5, height=5, path = "./plots/protein_pairs_examples/")
unlink("ddx5_cnv_ddx17_prot.png")



#DDX5 Prot -> DDX17 Prot'
plot14 <- ggplot(data = data %>% filter(controlling == "DDX5", controlled == "DDX17"), mapping = aes(x = proteomics_controlling, y = proteomics_controlled_)) +
    geom_point(aes(color=batch)) +
    geom_smooth(method=lm , color="black", fill = "grey", se=T, size=0.5) +
    stat_cor() +
    theme_classic() +
    theme(axis.title.y=element_text(colour="black", size=14),
          axis.title.x=element_text(colour="black", size=14),
          axis.text.y=element_text(colour="black", size=12),
          axis.text.x=element_text(colour="black", size=12),
          aspect.ratio=1,
          legend.position="bottom") +
    scale_x_continuous(name = "DDX5 proteomics (log2FC)") +
    scale_y_continuous(name = "DDX17 proteomics residuals (log2FC)")
ggsave(filename="ddx5_prot_ddx17_prot.png", plot=plot14, width=5, height=5, path = "./plots/protein_pairs_examples/")
unlink("ddx5_prot_ddx17_prot.png")



#CYLD CNV -> HSPA8 Prot'
#data %>% filter(controlling == "CYLD", controlled == "HSPA8") %>% group_by(controlling, controlled) %>% summarise(cor = cor.test(as.numeric(as.character(cnv_controlling)), proteomics_controlled_)$estimate)
plot15 <- ggplot(data = data %>% filter(controlling == "CYLD", controlled == "HSPA8"), mapping = aes(x = cnv_controlling, y = proteomics_controlled_)) +
    geom_boxplot(outlier.shape=NA, aes(fill=cnv_controlling)) +
    geom_jitter(width=0.1, alpha=0.5, pch=21, aes(fill=cnv_controlling), color = "black") +
    annotate("text", label=paste("r = ", "-0.29, ", "p < ", "8.8e-09"), x = 1.5, y = 2.5, size=3) +
    theme_classic() +
    theme(axis.title.y=element_text(colour="black", size=12),
          axis.title.x=element_text(colour="black", size=12),
          axis.text.y=element_text(colour="black", size=10),
          axis.text.x=element_text(colour="black", size=10),
          aspect.ratio = 1) +
    scale_x_discrete(name = "CYLD CNV") +
    scale_y_continuous(name = "HSPA8 protein residuals (log2FC)", limits = c(-2.5, 2.5)) +
    scale_fill_brewer(type = "seq", palette = "YlOrBr", guide=F)
ggsave(filename="cyld_cnv_hspa8_prot.png", plot=plot15, width=4, height=4, path = "./plots/protein_pairs_examples/")
unlink("cyld_cnv_hspa8_prot.png")



#CYLD Prot -> HSPA8 Prot'
plot16 <- ggplot(data = data %>% filter(controlling == "CYLD", controlled == "HSPA8"), mapping = aes(x = proteomics_controlling, y = proteomics_controlled_)) +
    geom_point(aes(color=batch)) +
    geom_smooth(method=lm , color="black", fill = "grey", se=T, size=0.5) +
    stat_cor(size = 3) +
    theme_classic() +
    theme(axis.title.y=element_text(colour="black", size=10),
          axis.title.x=element_text(colour="black", size=10),
          axis.text.y=element_text(colour="black", size=8),
          axis.text.x=element_text(colour="black", size=8),
          aspect.ratio=1,
          legend.position="bottom",
          legend.title=element_text(colour="black", size=10),
          legend.text=element_text(colour="black", size=8)) +
    scale_color_brewer(type = "div", palette = "RdYlGn", labels=c("Lawrence et al", "Lapek et al", "TCGA BRCA", "Roumeliotis et al", "TCGA COREAD", "TCGA HGSC"), name="Batch", guide=guide_legend(ncol=2), direction=-1) +
    scale_x_continuous(name = "CYLD protein (log2FC)") +
    scale_y_continuous(name = "HSPA8 protein residuals (log2FC)")
ggsave(filename="cyld_prot_hspa8_prot.png", plot=plot16, width=4, height=4, path = "./plots/protein_pairs_examples/")
unlink("cyld_prot_hspa8_prot.png")



#UFL1 CNV -> CDK5RAP3 Prot'
#data2 %>% filter(controlling == "UFL1", controlled == "CDK5RAP3") %>% group_by(controlling, controlled) %>% summarise(cor = cor.test(as.numeric(as.character(cnv_controlling)), proteomics_controlled_)$estimate)
plot17 <- ggplot(data = data2 %>% filter(controlling == "UFL1", controlled == "CDK5RAP3"), mapping = aes(x = cnv_controlling, y = proteomics_controlled_)) +
    geom_boxplot(outlier.shape=NA, aes(fill=cnv_controlling)) +
    geom_jitter(width=0.1, alpha=0.5, pch=21, aes(fill=cnv_controlling), color = "black") +
    annotate("text", label=paste("r = ", "0.41, ", "p < ", "1.58e-16"), x = 1.5, y = 3, size=3) +
    theme_classic() +
    theme(axis.title.y=element_text(colour="black", size=12),
          axis.title.x=element_text(colour="black", size=12),
          axis.text.y=element_text(colour="black", size=10),
          axis.text.x=element_text(colour="black", size=10),
          aspect.ratio = 1) +
    scale_x_discrete(name = "UFL1 CNV") +
    scale_y_continuous(name = "CDK5RAP3 protein residuals (log2FC)") +
    scale_fill_brewer(type = "seq", palette = "YlOrBr", guide=F)
ggsave(filename="ufl1_cnv_cdk5rap3_prot.png", plot=plot17, width=4, height=4, path = "./plots/protein_pairs_examples/")
unlink("ufl1_cnv_cdk5rap3_prot.png")



#UFL1 Prot -> CDK5RAP3 Prot'
plot18 <- ggplot(data = data2 %>% filter(controlling == "UFL1", controlled == "CDK5RAP3"), mapping = aes(x = proteomics_controlling, y = proteomics_controlled_)) +
    geom_point(aes(color=batch)) +
    geom_smooth(method=lm , color="black", fill = "grey", se=T, size=0.5) +
    stat_cor(size=3) +
    theme_classic() +
    theme(axis.title.y=element_text(colour="black", size=12),
          axis.title.x=element_text(colour="black", size=12),
          axis.text.y=element_text(colour="black", size=10),
          axis.text.x=element_text(colour="black", size=10),
          aspect.ratio=1) +
    scale_color_brewer(type = "div", palette = "RdYlGn", labels=c("Lawrence et al", "Lapek et al", "TCGA BRCA", "Roumeliotis et al", "TCGA COREAD", "TCGA HGSC"), name="Samples", guide=F, direction=-1) +
    scale_x_continuous(name = "UFL1 protein (log2FC)") +
    scale_y_continuous(name = "CDK5RAP3 protein residuals (log2FC)")
ggsave(filename="ufl1_prot_cdk5rap3_prot.png", plot=plot18, width=4, height=4, path = "./plots/protein_pairs_examples/")
unlink("ufl1_prot_cdk5rap3_prot.png")


#USP5 CNV -> PHB Prot'
#data2 %>% filter(controlling == "USP5", controlled == "PHB") %>% group_by(controlling, controlled) %>% summarise(cor = cor.test(as.numeric(as.character(cnv_controlling)), proteomics_controlled_)$estimate)
plot19 <- ggplot(data = data2 %>% filter(controlling == "USP5", controlled == "PHB"), mapping = aes(x = cnv_controlling, y = proteomics_controlled_)) +
    geom_boxplot(outlier.shape=NA, aes(fill=cnv_controlling)) +
    geom_jitter(width=0.1, alpha=0.5, pch=21, aes(fill=cnv_controlling), color = "black") +
    annotate("text", label=paste("r = ", "0.41, ", "p < ", "3.75e-16"), x = 1.5, y = 3, size=3) +
    theme_classic() +
    theme(axis.title.y=element_text(colour="black", size=12),
          axis.title.x=element_text(colour="black", size=12),
          axis.text.y=element_text(colour="black", size=10),
          axis.text.x=element_text(colour="black", size=10),
          aspect.ratio = 1) +
    scale_x_discrete(name = "USP5 CNV") +
    scale_y_continuous(name = "PHB protein residuals (log2FC)") +
    scale_fill_brewer(type = "seq", palette = "YlOrBr", guide=F)
ggsave(filename="usp5_cnv_phb_prot.png", plot=plot19, width=4, height=4, path = "./plots/protein_pairs_examples/")
unlink("usp5_cnv_phb_prot.png")



#USP5 Prot -> PHB Prot'
plot20 <- ggplot(data = data2 %>% filter(controlling == "USP5", controlled == "PHB"), mapping = aes(x = proteomics_controlling, y = proteomics_controlled_)) +
    geom_point(aes(color=batch)) +
    geom_smooth(method=lm , color="black", fill = "grey", se=T, size=0.5) +
    stat_cor(size=3) +
    theme_classic() +
    theme(axis.title.y=element_text(colour="black", size=12),
          axis.title.x=element_text(colour="black", size=12),
          axis.text.y=element_text(colour="black", size=10),
          axis.text.x=element_text(colour="black", size=10),
          aspect.ratio=1) +
    scale_color_brewer(type = "div", palette = "RdYlGn", labels=c("Lawrence et al", "Lapek et al", "TCGA BRCA", "Roumeliotis et al", "TCGA COREAD", "TCGA HGSC"), name="Samples", guide=F, direction=-1) +
    scale_x_continuous(name = "USP5 protein (log2FC)") +
    scale_y_continuous(name = "PHB protein residuals (log2FC)")
ggsave(filename="usp5_prot_phb_prot.png", plot=plot20, width=4, height=4, path = "./plots/protein_pairs_examples/")
unlink("usp5_prot_phb_prot.png")


#USP5 Prot -> PHB Prot'
plot21 <- ggplot(data = data2 %>% filter(controlling == "USP5", controlled == "PHB"), mapping = aes(x = proteomics_controlling, y = proteomics_controlled_)) +
    geom_point(aes(color=batch)) +
    theme_classic() +
    theme(legend.title=element_text(colour="black", size=14),
        legend.text=element_text(colour="black", size=12),
        legend.position="left") +
    scale_color_brewer(type = "div", palette = "RdYlGn", labels=c("Lawrence et al", "Lapek et al", "TCGA BRCA", "Roumeliotis et al", "TCGA COREAD", "TCGA HGSC"), name="Batch", direction=-1)
legend <- cowplot::get_legend(plot21)
ggsave(filename="samples_legend.png", plot=grid.draw(legend), width=2, height=2, path = "./plots/protein_pairs_examples/")
unlink("samples_legend.png")


#XRCC6 CNV -> XRCC5 Prot'
#data %>% filter(controlling == "XRCC6", controlled == "XRCC5") %>% group_by(controlling, controlled) %>% summarise(cor = cor.test(as.numeric(as.character(cnv_controlling)), proteomics_controlled_)$estimate)
plot22 <- ggplot(data = data %>% filter(controlling == "XRCC6", controlled == "XRCC5"), mapping = aes(x = cnv_controlling, y = proteomics_controlled_)) +
    geom_boxplot(outlier.shape=NA, aes(fill=cnv_controlling)) +
    geom_jitter(width=0.1, alpha=0.5, pch=21, aes(fill=cnv_controlling), color = "black", size=1) +
    annotate("text", label=paste("r = ", "0.36, ", "p < ", "1.27e-12"), x = 1.5, y = 3, size=2) +
    theme_classic() +
    theme(axis.title.y=element_text(colour="black", size=8),
          axis.title.x=element_text(colour="black", size=8),
          axis.text.y=element_text(colour="black", size=6),
          axis.text.x=element_text(colour="black", size=6),
          aspect.ratio = 1) +
    scale_x_discrete(name = "XRCC6 CNV") +
    scale_y_continuous(name = "XRCC5 protein residuals (log2FC)") +
    scale_fill_brewer(type = "seq", palette = "YlOrBr", guide=F)
ggsave(filename="xrcc6_cnv_xrcc5_prot.png", plot=plot22, width=2, height=2, path = "./plots/protein_pairs_examples/")
unlink("xrcc6_cnv_xrcc5_prot.png")



#XRCC6 CNV -> XRCC5 Prot'
plot23 <- ggplot(data = data %>% filter(controlling == "XRCC6", controlled == "XRCC5"), mapping = aes(x = proteomics_controlling, y = proteomics_controlled_)) +
    geom_point(aes(color=batch), size=1) +
    geom_smooth(method=lm , color="black", fill = "grey", se=T, size=0.5) +
    stat_cor(size=2) +
    theme_classic() +
    theme(axis.title.y=element_text(colour="black", size=8),
          axis.title.x=element_text(colour="black", size=8),
          axis.text.y=element_text(colour="black", size=6),
          axis.text.x=element_text(colour="black", size=6),
          aspect.ratio=1) +
    scale_color_brewer(type = "div", palette = "RdYlGn", labels=c("Lawrence et al", "Lapek et al", "TCGA BRCA", "Roumeliotis et al", "TCGA COREAD", "TCGA HGSC"), name="Samples", guide=F, direction=-1) +
    scale_x_continuous(name = "XRCC6 protein (log2FC)") +
    scale_y_continuous(name = "XRCC5 protein residuals (log2FC)")
ggsave(filename="xrcc6_prot_xrcc5_prot.png", plot=plot23, width=2, height=2, path = "./plots/protein_pairs_examples/")
unlink("xrcc6_prot_xrcc5_prot.png")



#HSPD1 CNV -> DRG1 Prot'
#data2 %>% filter(controlling == "HSPD1", controlled == "DRG1") %>% group_by(controlling, controlled) %>% summarise(cor = cor.test(as.numeric(as.character(cnv_controlling)), proteomics_controlled_)$estimate)
plot24 <- ggplot(data = data2 %>% filter(controlling == "HSPD1", controlled == "DRG1"), mapping = aes(x = cnv_controlling, y = proteomics_controlled_)) +
    geom_boxplot(outlier.shape=NA, aes(fill=cnv_controlling)) +
    geom_jitter(width=0.1, alpha=0.5, pch=21, aes(fill=cnv_controlling), color = "black") +
    annotate("text", label=paste("r = ", "0.34, ", "p < ", "1.83e-11"), x = 1.5, y = 3, size=3) +
    theme_classic() +
    theme(axis.title.y=element_text(colour="black", size=12),
          axis.title.x=element_text(colour="black", size=12),
          axis.text.y=element_text(colour="black", size=10),
          axis.text.x=element_text(colour="black", size=10),
          aspect.ratio = 1) +
    scale_x_discrete(name = "HSPD1 CNV") +
    scale_y_continuous(name = "DRG1 protein residuals (log2FC)") +
    scale_fill_brewer(type = "seq", palette = "YlOrBr", guide=F)
ggsave(filename="hspd1_cnv_drg1_prot.png", plot=plot24, width=4, height=4, path = "./plots/protein_pairs_examples/")
unlink("hspd1_cnv_drg1_prot.png")



#HSPD1 Prot -> DRG1 Prot'
plot25 <- ggplot(data = data2 %>% filter(controlling == "HSPD1", controlled == "DRG1"), mapping = aes(x = proteomics_controlling, y = proteomics_controlled_)) +
    geom_point(aes(color=batch)) +
    geom_smooth(method=lm , color="black", fill = "grey", se=T, size=0.5) +
    stat_cor(size=3) +
    theme_classic() +
    theme(axis.title.y=element_text(colour="black", size=12),
          axis.title.x=element_text(colour="black", size=12),
          axis.text.y=element_text(colour="black", size=10),
          axis.text.x=element_text(colour="black", size=10),
          aspect.ratio=1) +
    scale_color_brewer(type = "div", palette = "RdYlGn", labels=c("Lawrence et al", "Lapek et al", "TCGA BRCA", "Roumeliotis et al", "TCGA COREAD", "TCGA HGSC"), name="Samples", guide=F, direction=-1) +
    scale_x_continuous(name = "HSPD1 protein (log2FC)") +
    scale_y_continuous(name = "DRG1 protein residuals (log2FC)")
ggsave(filename="hspd1_prot_drg1_prot.png", plot=plot25, width=4, height=4, path = "./plots/protein_pairs_examples/")
unlink("hspd1_prot_drg1_prot.png")


#MARCKS CNV -> HDAC1 Prot'
#data2 %>% filter(controlling == "MARCKS", controlled == "HDAC1") %>% group_by(controlling, controlled) %>% summarise(cor = cor.test(as.numeric(as.character(cnv_controlling)), proteomics_controlled_)$estimate)
plot26 <- ggplot(data = data2 %>% filter(controlling == "MARCKS", controlled == "HDAC1"), mapping = aes(x = cnv_controlling, y = proteomics_controlled_)) +
    geom_boxplot(outlier.shape=NA, aes(fill=cnv_controlling)) +
    geom_jitter(width=0.1, alpha=0.5, pch=21, aes(fill=cnv_controlling), color = "black") +
    annotate("text", label=paste("r = ", " -0.25, ", "p < ", "1.32e-06"), x = 1.5, y = 3, size=3) +
    theme_classic() +
    theme(axis.title.y=element_text(colour="black", size=12),
          axis.title.x=element_text(colour="black", size=12),
          axis.text.y=element_text(colour="black", size=10),
          axis.text.x=element_text(colour="black", size=10),
          aspect.ratio = 1) +
    scale_x_discrete(name = "MARCKS CNV") +
    scale_y_continuous(name = "HDAC1 protein residuals (log2FC)") +
    scale_fill_brewer(type = "seq", palette = "YlOrBr", guide=F)
ggsave(filename="MARCKS_cnv_HDAC1_prot.png", plot=plot26, width=4, height=4, path = "./plots/protein_pairs_examples/")
unlink("MARCKS_cnv_HDAC1_prot.png")



#HSPD1 Prot -> DRG1 Prot'
plot27 <- ggplot(data = data2 %>% filter(controlling == "MARCKS", controlled == "HDAC1"), mapping = aes(x = proteomics_controlling, y = proteomics_controlled_)) +
    geom_point(aes(color=batch)) +
    geom_smooth(method=lm , color="black", fill = "grey", se=T, size=0.5) +
    stat_cor(size=3) +
    theme_classic() +
    theme(axis.title.y=element_text(colour="black", size=12),
          axis.title.x=element_text(colour="black", size=12),
          axis.text.y=element_text(colour="black", size=10),
          axis.text.x=element_text(colour="black", size=10),
          aspect.ratio=1) +
    scale_color_brewer(type = "div", palette = "RdYlGn", labels=c("Lawrence et al", "Lapek et al", "TCGA BRCA", "Roumeliotis et al", "TCGA COREAD", "TCGA HGSC"), name="Samples", guide=F, direction=-1) +
    scale_x_continuous(name = "MARCKS protein (log2FC)") +
    scale_y_continuous(name = "HDAC1 protein residuals (log2FC)")
ggsave(filename="MARCKS_prot_HDAC1_prot.png", plot=plot27, width=4, height=4, path = "./plots/protein_pairs_examples/")
unlink("MARCKS_prot_HDAC1_prot.png")




save(list=ls(), file="./r_workspaces/protein_pairs_examples.RData")
