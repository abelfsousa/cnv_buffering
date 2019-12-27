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



# significant associations with CNV and Phospho
signf_cnv_rna_phospho <- read_tsv("./files/signf_protein_associations_cnv_rna_phospho.txt")



# load common proteomics phosphoproteomics rna cnv metadata samples
common_samples <- read_tsv("./files/cptac_cellLines_samples_prot_phospho_rna_cnv.txt", col_names = FALSE)

cptac_samples <- read_tsv("./files/cptac_samples.txt")

# load metadata
metadata <- read_tsv("./files/metadata_cptac_cellLines.txt") %>%
  arrange(sample)


cnv <- read_tsv("./files/cnv_tcga_cellLines.txt") %>%
  gather(key="sample", value="cnv_gistic2", -gene) %>%
  arrange(sample) %>%
  filter(sample %in% common_samples$X1)


rna <- read_tsv("./files/rna_tcga_cellLines.txt") %>%
  gather(key="sample", value="rna_log2CPM", -gene) %>%
  arrange(sample) %>%
  filter(sample %in% common_samples$X1) %>%
  group_by(gene) %>%
  mutate(rna_log2CPM = scale(rna_log2CPM)[,1]) %>%
  ungroup()


proteomics <- read_tsv("./files/proteomicsQ_cptac_cellLines.txt") %>%
  gather(key="sample", value="prot_log2FC", -gene) %>%
  arrange(sample) %>%
  filter(sample %in% common_samples$X1) %>%
  group_by(gene) %>%
  filter( sum(!is.na(prot_log2FC)) > n()*0.5 ) %>%
  mutate(prot_log2FC = scale(prot_log2FC)[,1]) %>%
  ungroup()


phospho <- read_tsv("./files/phosphoproteomicsQ_cptac_cellLines.txt") %>%
  gather(key="sample", value="phospho_log2FC", -phospho_site) %>%
  arrange(sample) %>%
  filter(sample %in% common_samples$X1) %>%
  group_by(phospho_site) %>%
  filter( sum(!is.na(phospho_log2FC)) > n()*0.5 ) %>%
  mutate(phospho_log2FC = scale(phospho_log2FC)[,1]) %>%
  ungroup()



linear_model1 <- function(data){
  data <- as.data.frame(data)
  rownames(data) <- data$sample

  covars <- data[, -c(1)]

  # remove nas
  covars <- na.exclude(covars)


  # remove covariates with less than two levels
  covars <- covars[, apply(covars, 2, function(x) length(levels(as.factor(x))) >= 2), drop=F]


  # fit the model
  lreg_f <- as.formula( paste( "phospho_controlling", "~", paste(colnames(covars[, -match("phospho_controlling", colnames(covars)), drop=F]), collapse="+") ) )
  lreg <- lm(lreg_f, data = covars)


  # model residuals
  resd <- as.tibble(data.frame(sample = names(residuals(lreg)), phospho_controlling_ = unname(residuals(lreg)), stringsAsFactors = F))

  return(resd)

}



data1 <- signf_cnv_rna_phospho %>%
  dplyr::select(controlling, controlling_phx, controlled) %>%
  inner_join(phospho, by=c("controlling_phx" = "phospho_site")) %>%
  dplyr::rename(phospho_controlling = phospho_log2FC) %>%
  inner_join(proteomics, by=c("controlling" = "gene", "sample")) %>%
  dplyr::rename(proteomics_controlling = prot_log2FC) %>%
  inner_join(metadata, by = "sample") %>%
  #regress-out protein from phospho-site abundance
  nest(-c(controlling, controlling_phx, controlled)) %>%
  mutate(phospho_controlling_ = purrr::map(data, ~ linear_model1(.))) %>%
  dplyr::select(-data) %>%
  unnest(phospho_controlling_)



linear_model2 <- function(data){
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


data2 <- signf_cnv_rna_phospho %>%
  dplyr::select(controlling, controlling_phx, controlled) %>%
  inner_join(proteomics, by=c("controlled" = "gene")) %>%
  dplyr::rename(proteomics_controlled = prot_log2FC) %>%
  inner_join(rna, by=c("controlled" = "gene", "sample")) %>%
  dplyr::rename(rna_controlled = rna_log2CPM) %>%
  inner_join(metadata, by = "sample") %>%
  #regress-out rna from controlled protein abundance
  nest(-c(controlling, controlling_phx, controlled)) %>%
  mutate(proteomics_controlled_ = purrr::map(data, ~ linear_model2(.))) %>%
  dplyr::select(-data) %>%
  unnest(proteomics_controlled_) %>%
  inner_join(cnv, by=c("controlling" = "gene", "sample")) %>%
  dplyr::rename(cnv_controlling = cnv_gistic2) %>%
  inner_join(proteomics, by=c("controlling" = "gene", "sample")) %>%
  dplyr::rename(proteomics_controlling = prot_log2FC) %>%
  mutate(pair = paste(controlling_phx, controlled, sep="_")) %>%
  dplyr::select(pair, controlling, controlling_phx, controlled, sample, everything()) %>%
  mutate(cnv_controlling = as.factor(cnv_controlling)) %>%
  inner_join(metadata[, c("sample", "batch")], by = "sample") %>%
  mutate_at(vars(contains("batch")), as.factor) %>%
  mutate(batch = fct_relevel(batch, "TCGA_BRCA", "RMLT", "TCGA_OV"))


data <- inner_join(data1, data2, by=c("controlling", "controlling_phx", "controlled", "sample"))


#POLD3 CNV -> POLD2 Prot'
#data %>% filter(controlling_phx == "POLD3_s307", controlled == "POLD2") %>% group_by(controlling_phx, controlled) %>% summarise(cor = cor.test(as.numeric(as.character(cnv_controlling)), proteomics_controlled_)$estimate)
plot1 <- ggplot(data = data %>% filter(controlling_phx == "POLD3_s307", controlled == "POLD2"), mapping = aes(x = cnv_controlling, y = proteomics_controlled_)) +
  geom_boxplot(outlier.shape=NA, aes(fill=cnv_controlling)) +
  geom_jitter(width=0.1, alpha=0.5, pch=21, aes(fill=cnv_controlling), color = "black") +
  annotate("text", label=paste("r = 0.35", ", ", "p < ", "6.94e-06"), x = 2, y = 3, size=3) +
  theme_classic() +
  theme(axis.title.y=element_text(colour="black", size=10),
    axis.title.x=element_text(colour="black", size=10),
    axis.text.y=element_text(colour="black", size=8),
    axis.text.x=element_text(colour="black", size=8),
    legend.position = "none",
    aspect.ratio=1) +
  scale_x_discrete(name = "POLD3 CNV") +
  scale_y_continuous(name = "POLD2 protein residuals (log2FC)") +
  scale_fill_brewer(type = "seq", palette = "YlOrBr", guide=F)
ggsave(filename="pold3_cnv_pold2_prot.png", plot=plot1, height=3, width=3, path = "./plots/phospho_pairs_examples/")
unlink("pold3_cnv_pold2_prot.png")



#POLD3 prot -> POLD2 Prot'
plot2 <- ggplot(data = data %>% filter(controlling_phx == "POLD3_s307", controlled == "POLD2"), mapping = aes(x = proteomics_controlling, y = proteomics_controlled_)) +
    geom_point(aes(color=batch)) +
    geom_smooth(method=lm , color="black", fill = "grey", se=T, size=0.5) +
    stat_cor(size=3) +
    theme_classic() +
    theme(axis.title.y=element_text(colour="black", size=10),
          axis.title.x=element_text(colour="black", size=10),
          axis.text.y=element_text(colour="black", size=8),
          axis.text.x=element_text(colour="black", size=8),
          legend.position = "none",
          aspect.ratio=1) +
    scale_color_brewer(type = "div", palette = "RdYlGn", labels=c("TCGA BRCA", "Roumeliotis et al", "TCGA HGSC"), name="Samples", guide=F, direction=-1) +
    scale_x_continuous(name = "POLD3 protein (log2FC)") +
    scale_y_continuous(name = "POLD2 protein residuals (log2FC)")
ggsave(filename="pold3_prot_pold2_prot.png", plot=plot2, width=3, height=3, path = "./plots/phospho_pairs_examples/")
unlink("pold3_prot_pold2_prot.png")


#POLD3_s307 phospho -> POLD2 Prot'
plot3 <- ggplot(data = data %>% filter(controlling_phx == "POLD3_s307", controlled == "POLD2"), mapping = aes(x = phospho_controlling_, y = proteomics_controlled_)) +
    geom_point(aes(color=batch)) +
    geom_smooth(method=lm , color="black", fill = "grey", se=T, size=0.5) +
    stat_cor() +
    theme_classic() +
    theme(axis.title.y=element_text(colour="black", size=14),
          axis.title.x=element_text(colour="black", size=14),
          axis.text.y=element_text(colour="black", size=12),
          axis.text.x=element_text(colour="black", size=12),
          legend.position = "bottom",
          aspect.ratio=1) +
    scale_x_continuous(name = "POLD3_s307 phospho residuals (log2FC)") +
    scale_y_continuous(name = "POLD2 proteomics residuals (log2FC)")
ggsave(filename="pold3s307_phospho_pold2_prot.png", plot=plot3, width=5, height=5, path = "./plots/phospho_pairs_examples/")
unlink("pold3s307_phospho_pold2_prot.png")


#EIF3A CNV -> EIF3D Prot'
#data %>% filter(controlling_phx == "EIF3A_s492", controlled == "EIF3D") %>% group_by(controlling_phx, controlled) %>% summarise(cor = cor.test(as.numeric(as.character(cnv_controlling)), proteomics_controlled_)$estimate)
plot4 <- ggplot(data = data %>% filter(controlling_phx == "EIF3A_s492", controlled == "EIF3D"), mapping = aes(x = cnv_controlling, y = proteomics_controlled_)) +
  geom_boxplot(outlier.shape=NA, aes(fill=cnv_controlling)) +
  geom_jitter(width=0.1, alpha=0.5, pch=21, aes(fill=cnv_controlling), color = "black") +
  annotate("text", label=paste("R = 0.29", ", ", "p < ", "2.57e-3"), x = 2, y = 3, size=3) +
  theme_classic() +
  theme(
    axis.title=element_text(colour="black", size=12),
    axis.text=element_text(colour="black", size=10),
    aspect.ratio=1) +
  scale_x_discrete(name = "EIF3A CNV") +
  scale_y_continuous(name = "EIF3D protein residuals (log2FC)") +
  scale_fill_brewer(type = "seq", palette = "YlOrBr", guide=F)
ggsave(filename="eif3a_cnv_eif3d_prot.png", plot=plot4, height=3, width=3, path = "./plots/phospho_pairs_examples/")
ggsave(filename="eif3a_cnv_eif3d_prot.pdf", plot=plot4, height=3, width=3, path = "./plots/phospho_pairs_examples/")
unlink("eif3a_cnv_eif3d_prot.png")
unlink("eif3a_cnv_eif3d_prot.pdf")


#EIF3A prot -> EIF3D Prot'
plot5 <- ggplot(data = data %>% filter(controlling_phx == "EIF3A_s492", controlled == "EIF3D"), mapping = aes(x = proteomics_controlling, y = proteomics_controlled_)) +
    geom_point(aes(color=batch)) +
    geom_smooth(method=lm , color="black", fill = "grey", se=T, size=0.5) +
    stat_cor(size=3) +
    theme_classic() +
    theme(
      axis.title=element_text(colour="black", size=12),
      axis.text=element_text(colour="black", size=10),
      aspect.ratio=1) +
    scale_color_manual(values=c("#1a9850", "#fc8d59", "#d73027"), labels=c("TCGA BRCA", "Roumeliotis et al", "TCGA HGSC"), name="Batch", guide=F) +
    scale_x_continuous(name = "EIF3A protein (log2FC)") +
    scale_y_continuous(name = "EIF3D protein residuals (log2FC)")
ggsave(filename="eif3a_prot_eif3d_prot.png", plot=plot5, width=3, height=3, path = "./plots/phospho_pairs_examples/")
ggsave(filename="eif3a_prot_eif3d_prot.pdf", plot=plot5, width=3, height=3, path = "./plots/phospho_pairs_examples/")
unlink("eif3a_prot_eif3d_prot.png")
unlink("eif3a_prot_eif3d_prot.pdf")


#EIF3A prot -> EIF3D Prot'
#TCGA BRCA
plot5_tcga_brca <- ggplot(data = data %>% filter(controlling_phx == "EIF3A_s492", controlled == "EIF3D") %>% filter(batch == "TCGA_BRCA"), mapping = aes(x = proteomics_controlling, y = proteomics_controlled_)) +
    geom_point(color = "red") +
    geom_smooth(method=lm , color="black", fill = "grey", se=T, size=0.5) +
    stat_cor(size=3) +
    theme_classic() +
    theme(axis.title.y=element_text(colour="black", size=10),
          axis.title.x=element_text(colour="black", size=10),
          axis.text.y=element_text(colour="black", size=8),
          axis.text.x=element_text(colour="black", size=8),
          aspect.ratio=1) +
    scale_x_continuous(name = "EIF3A protein (log2FC)") +
    scale_y_continuous(name = "EIF3D protein residuals (log2FC)") +
    labs(title = "TCGA BRCA")
ggsave(filename="eif3a_prot_eif3d_prot_tcga_brca.png", plot=plot5_tcga_brca, width=3, height=3, path = "./plots/phospho_pairs_examples/")
ggsave(filename="eif3a_prot_eif3d_prot_tcga_brca.pdf", plot=plot5_tcga_brca, width=3, height=3, path = "./plots/phospho_pairs_examples/")
unlink("eif3a_prot_eif3d_prot_tcga_brca.png")
unlink("eif3a_prot_eif3d_prot_tcga_brca.pdf")


#EIF3A prot -> EIF3D Prot'
#RMLT
plot5_rmlt <- ggplot(data = data %>% filter(controlling_phx == "EIF3A_s492", controlled == "EIF3D") %>% filter(batch == "RMLT"), mapping = aes(x = proteomics_controlling, y = proteomics_controlled_)) +
    geom_point(color = "red") +
    geom_smooth(method=lm , color="black", fill = "grey", se=T, size=0.5) +
    stat_cor(size=3) +
    theme_classic() +
    theme(axis.title.y=element_text(colour="black", size=10),
          axis.title.x=element_text(colour="black", size=10),
          axis.text.y=element_text(colour="black", size=8),
          axis.text.x=element_text(colour="black", size=8),
          aspect.ratio=1) +
    scale_x_continuous(name = "EIF3A protein (log2FC)") +
    scale_y_continuous(name = "EIF3D protein residuals (log2FC)") +
    labs(title = "Roumeliotis et al")
ggsave(filename="eif3a_prot_eif3d_prot_rmlt.png", plot=plot5_rmlt, width=3, height=3, path = "./plots/phospho_pairs_examples/")
ggsave(filename="eif3a_prot_eif3d_prot_rmlt.pdf", plot=plot5_rmlt, width=3, height=3, path = "./plots/phospho_pairs_examples/")
unlink("eif3a_prot_eif3d_prot_rmlt.png")
unlink("eif3a_prot_eif3d_prot_rmlt.pdf")



#EIF3A_s492 phospho -> EIF3D Prot'
plot6 <- ggplot(data = data %>% filter(controlling_phx == "EIF3A_s492", controlled == "EIF3D"), mapping = aes(x = phospho_controlling_, y = proteomics_controlled_)) +
    geom_point(aes(color=batch)) +
    geom_smooth(method=lm , color="black", fill = "grey", se=T, size=0.5) +
    stat_cor(size=3) +
    theme_classic() +
    theme(
      axis.title=element_text(colour="black", size=12),
      axis.text=element_text(colour="black", size=10),
      aspect.ratio=1) +
    scale_color_manual(values=c("#1a9850", "#fc8d59", "#d73027"), labels=c("TCGA BRCA", "Roumeliotis et al", "TCGA HGSC"), name="Batch", guide=F) +
    scale_x_continuous(name = "EIF3A_s492 phospho residuals (log2FC)") +
    scale_y_continuous(name = "EIF3D protein residuals (log2FC)")
ggsave(filename="eif3as492_phospho_eif3d_prot.png", plot=plot6, width=3, height=3, path = "./plots/phospho_pairs_examples/")
ggsave(filename="eif3as492_phospho_eif3d_prot.pdf", plot=plot6, width=3, height=3, path = "./plots/phospho_pairs_examples/")
unlink("eif3as492_phospho_eif3d_prot.png")
unlink("eif3as492_phospho_eif3d_prot.pdf")


#EIF3A_s492 phospho -> EIF3D Prot'
#TCGA BRCA
plot6_tcga_brca <- ggplot(data = data %>% filter(controlling_phx == "EIF3A_s492", controlled == "EIF3D") %>% filter(batch == "TCGA_BRCA"), mapping = aes(x = phospho_controlling_, y = proteomics_controlled_)) +
    geom_point(color="red") +
    geom_smooth(method=lm , color="black", fill = "grey", se=T, size=0.5) +
    stat_cor(size=3) +
    theme_classic() +
    theme(axis.title.y=element_text(colour="black", size=10),
          axis.title.x=element_text(colour="black", size=10),
          axis.text.y=element_text(colour="black", size=8),
          axis.text.x=element_text(colour="black", size=8),
          aspect.ratio=1) +
    scale_x_continuous(name = "EIF3A_s492 phospho residuals (log2FC)") +
    scale_y_continuous(name = "EIF3D protein residuals (log2FC)") +
    labs(title = "TCGA BRCA")
ggsave(filename="eif3as492_phospho_eif3d_prot_tcga_brca.png", plot=plot6_tcga_brca, width=3, height=3, path = "./plots/phospho_pairs_examples/")
ggsave(filename="eif3as492_phospho_eif3d_prot_tcga_brca.pdf", plot=plot6_tcga_brca, width=3, height=3, path = "./plots/phospho_pairs_examples/")
unlink("eif3as492_phospho_eif3d_prot_tcga_brca.png")
unlink("eif3as492_phospho_eif3d_prot_tcga_brca.pdf")


#EIF3A_s492 phospho -> EIF3D Prot'
#RMLT
plot6_tcga_rmlt <- ggplot(data = data %>% filter(controlling_phx == "EIF3A_s492", controlled == "EIF3D") %>% filter(batch == "RMLT"), mapping = aes(x = phospho_controlling_, y = proteomics_controlled_)) +
    geom_point(color="red") +
    geom_smooth(method=lm , color="black", fill = "grey", se=T, size=0.5) +
    stat_cor(size=3) +
    theme_classic() +
    theme(axis.title.y=element_text(colour="black", size=10),
          axis.title.x=element_text(colour="black", size=10),
          axis.text.y=element_text(colour="black", size=8),
          axis.text.x=element_text(colour="black", size=8),
          aspect.ratio=1) +
    scale_x_continuous(name = "EIF3A_s492 phospho residuals (log2FC)") +
    scale_y_continuous(name = "EIF3D protein residuals (log2FC)") +
    labs(title = "Roumeliotis et al")
ggsave(filename="eif3as492_phospho_eif3d_prot_rmlt.png", plot=plot6_tcga_rmlt, width=3, height=3, path = "./plots/phospho_pairs_examples/")
ggsave(filename="eif3as492_phospho_eif3d_prot_rmlt.pdf", plot=plot6_tcga_rmlt, width=3, height=3, path = "./plots/phospho_pairs_examples/")
unlink("eif3as492_phospho_eif3d_prot_rmlt.png")
unlink("eif3as492_phospho_eif3d_prot_rmlt.pdf")






#EIF3A CNV -> EIF3E Prot'
#data %>% filter(controlling_phx == "EIF3A_s492", controlled == "EIF3E") %>% group_by(controlling_phx, controlled) %>% summarise(cor = cor.test(as.numeric(as.character(cnv_controlling)), proteomics_controlled_)$estimate)
plot7 <- ggplot(data = data %>% filter(controlling_phx == "EIF3A_s492", controlled == "EIF3E"), mapping = aes(x = cnv_controlling, y = proteomics_controlled_)) +
  geom_boxplot(outlier.shape=NA, aes(fill=cnv_controlling)) +
  geom_jitter(width=0.1, alpha=0.5, pch=21, aes(fill=cnv_controlling), color = "black") +
  annotate("text", label=paste("r = 0.23", ", ", "p < ", "1.44e-2"), x = 1.5, y = 3) +
  theme_classic() +
  theme(axis.title.y=element_text(colour="black", size=14),
    axis.title.x=element_text(colour="black", size=14),
    axis.text.y=element_text(colour="black", size=12),
    axis.text.x=element_text(colour="black", size=12),
    legend.position = "none",
    aspect.ratio=1) +
  scale_x_discrete(name = "EIF3A CNV (gistic2)") +
  scale_y_continuous(name = "EIF3E proteomics residuals (log2FC)") +
  scale_fill_brewer(type = "qual", palette = "GnBu")
ggsave(filename="eif3a_cnv_eif3e_prot.png", plot=plot7, height=5, width=5, path = "./plots/phospho_pairs_examples/")
unlink("eif3a_cnv_eif3e_prot.png")


#EIF3A prot -> EIF3E Prot'
plot8 <- ggplot(data = data %>% filter(controlling_phx == "EIF3A_s492", controlled == "EIF3E"), mapping = aes(x = proteomics_controlling, y = proteomics_controlled_)) +
    geom_point(aes(color=batch)) +
    geom_smooth(method=lm , color="black", fill = "grey", se=T, size=0.5) +
    stat_cor() +
    theme_classic() +
    theme(axis.title.y=element_text(colour="black", size=14),
          axis.title.x=element_text(colour="black", size=14),
          axis.text.y=element_text(colour="black", size=12),
          axis.text.x=element_text(colour="black", size=12),
          legend.position = "bottom",
          aspect.ratio=1) +
    scale_x_continuous(name = "EIF3A proteomics (log2FC)") +
    scale_y_continuous(name = "EIF3E proteomics residuals (log2FC)")
ggsave(filename="eif3a_prot_eif3e_prot.png", plot=plot8, width=5, height=5, path = "./plots/phospho_pairs_examples/")
unlink("eif3a_prot_eif3e_prot.png")


#EIF3A_s492 phospho -> EIF3E Prot'
plot9 <- ggplot(data = data %>% filter(controlling_phx == "EIF3A_s492", controlled == "EIF3E"), mapping = aes(x = phospho_controlling_, y = proteomics_controlled_)) +
    geom_point(aes(color=batch)) +
    geom_smooth(method=lm , color="black", fill = "grey", se=T, size=0.5) +
    stat_cor() +
    theme_classic() +
    theme(axis.title.y=element_text(colour="black", size=14),
          axis.title.x=element_text(colour="black", size=14),
          axis.text.y=element_text(colour="black", size=12),
          axis.text.x=element_text(colour="black", size=12),
          legend.position = "bottom",
          aspect.ratio=1) +
    scale_x_continuous(name = "EIF3A_s492 phospho residuals (log2FC)") +
    scale_y_continuous(name = "EIF3E proteomics residuals (log2FC)")
ggsave(filename="eif3as492_phospho_eif3e_prot.png", plot=plot9, width=5, height=5, path = "./plots/phospho_pairs_examples/")
unlink("eif3as492_phospho_eif3e_prot.png")


#POLD3_s458 phospho -> POLD2 Prot'
plot10 <- ggplot(data = data %>% filter(controlling_phx == "POLD3_s458", controlled == "POLD2"), mapping = aes(x = phospho_controlling_, y = proteomics_controlled_)) +
    geom_point(aes(color=batch)) +
    geom_smooth(method=lm , color="black", fill = "grey", se=T, size=0.5) +
    stat_cor(size=3) +
    theme_classic() +
    theme(axis.title.y=element_text(colour="black", size=10),
          axis.title.x=element_text(colour="black", size=10),
          axis.text.y=element_text(colour="black", size=8),
          axis.text.x=element_text(colour="black", size=8),
          aspect.ratio=1) +
    scale_color_brewer(type = "div", palette = "RdYlGn", labels=c("TCGA BRCA", "Roumeliotis et al", "TCGA HGSC"), name="Samples", guide=F, direction=-1) +
    scale_x_continuous(name = "POLD3_s458 phospho residuals (log2FC)") +
    scale_y_continuous(name = "POLD2 protein residuals (log2FC)")
ggsave(filename="pold3s458_phospho_pold2_prot.png", plot=plot10, width=3, height=3, path = "./plots/phospho_pairs_examples/")
unlink("pold3s458_phospho_pold2_prot.png")



#EIF3A prot -> EIF3D Prot'
#legend
plot11 <- ggplot(data = data %>% filter(controlling_phx == "EIF3A_s492", controlled == "EIF3D"), mapping = aes(x = proteomics_controlling, y = proteomics_controlled_)) +
    geom_point(aes(color=batch)) +
    theme_classic() +
    theme(legend.title=element_text(colour="black", size=12),
        legend.text=element_text(colour="black", size=10),
        legend.position="left") +
    scale_color_manual(values=c("#1a9850", "#fc8d59", "#d73027"), labels=c("TCGA BRCA", "Roumeliotis et al", "TCGA HGSC"), name="Batch")
legend <- cowplot::get_legend(plot11)
ggsave(filename="samples_legend.png", plot=grid.draw(legend), width=1.5, height=1, path = "./plots/phospho_pairs_examples/")
ggsave(filename="samples_legend.pdf", plot=grid.draw(legend), width=1.5, height=1, path = "./plots/phospho_pairs_examples/")
unlink("samples_legend.png")
unlink("samples_legend.pdf")


#MYH9 CNV -> MYL6 Prot'
#data %>% filter(controlling_phx == "MYH9_s1943", controlled == "MYL6") %>% group_by(controlling_phx, controlled) %>% summarise(cor = cor.test(as.numeric(as.character(cnv_controlling)), proteomics_controlled_)$estimate)
plot12 <- ggplot(data = data %>% filter(controlling_phx == "MYH9_s1943", controlled == "MYL6"), mapping = aes(x = cnv_controlling, y = proteomics_controlled_)) +
  geom_boxplot(outlier.shape=NA, aes(fill=cnv_controlling)) +
  geom_jitter(width=0.1, alpha=0.5, pch=21, aes(fill=cnv_controlling), color = "black") +
  annotate("text", label=paste("r = 0.12", ", ", "p < ", "0.12"), x = 2, y = 3, size=3) +
  theme_classic() +
  theme(axis.title.y=element_text(colour="black", size=10),
    axis.title.x=element_text(colour="black", size=10),
    axis.text.y=element_text(colour="black", size=8),
    axis.text.x=element_text(colour="black", size=8),
    aspect.ratio=1) +
  scale_x_discrete(name = "MYH9 CNV") +
  scale_y_continuous(name = "MYL6 protein residuals (log2FC)") +
  scale_fill_brewer(type = "seq", palette = "YlOrBr", guide=F)
ggsave(filename="MYH9_cnv_MYL6_prot.png", plot=plot12, height=3, width=3, path = "./plots/phospho_pairs_examples/")
unlink("MYH9_cnv_MYL6_prot.png")


#MYH9 prot -> MYL6 Prot'
plot13 <- ggplot(data = data %>% filter(controlling_phx == "MYH9_s1943", controlled == "MYL6"), mapping = aes(x = proteomics_controlling, y = proteomics_controlled_)) +
    geom_point(aes(color=batch)) +
    geom_smooth(method=lm , color="black", fill = "grey", se=T, size=0.5) +
    stat_cor(size=3) +
    theme_classic() +
    theme(axis.title.y=element_text(colour="black", size=10),
          axis.title.x=element_text(colour="black", size=10),
          axis.text.y=element_text(colour="black", size=8),
          axis.text.x=element_text(colour="black", size=8),
          aspect.ratio=1) +
    scale_color_manual(values=c("#1a9850", "#fc8d59", "#d73027"), labels=c("TCGA BRCA", "Roumeliotis et al", "TCGA HGSC"), name="Batch", guide=F) +
    scale_x_continuous(name = "MYH9 protein (log2FC)") +
    scale_y_continuous(name = "MYL6 protein residuals (log2FC)")
ggsave(filename="MYH9_prot_MYL6_prot.png", plot=plot13, width=3, height=3, path = "./plots/phospho_pairs_examples/")
unlink("MYH9_prot_MYL6_prot.png")


#MYH9_s1943 phospho -> MYL6 Prot'
plot14 <- ggplot(data = data %>% filter(controlling_phx == "MYH9_s1943", controlled == "MYL6") %>% filter(phospho_controlling_ > -1.5), mapping = aes(x = phospho_controlling_, y = proteomics_controlled_)) +
    geom_point(aes(color=batch)) +
    geom_smooth(method=lm , color="black", fill = "grey", se=T, size=0.5) +
    stat_cor(size=3) +
    theme_classic() +
    theme(axis.title.y=element_text(colour="black", size=10),
          axis.title.x=element_text(colour="black", size=10),
          axis.text.y=element_text(colour="black", size=8),
          axis.text.x=element_text(colour="black", size=8),
          aspect.ratio=1) +
    scale_color_manual(values=c("#1a9850", "#fc8d59", "#d73027"), labels=c("TCGA BRCA", "Roumeliotis et al", "TCGA HGSC"), name="Batch", guide=F) +
    scale_x_continuous(name = "MYH9_s1943 phospho residuals (log2FC)") +
    scale_y_continuous(name = "MYL6 protein residuals (log2FC)")
ggsave(filename="MYH9_s1943_phospho_MYL6_prot.png", plot=plot14, width=3, height=3, path = "./plots/phospho_pairs_examples/")
unlink("MYH9_s1943_phospho_MYL6_prot.png")


save(list=ls(), file="./r_workspaces/phospho_pairs_examples.RData")
