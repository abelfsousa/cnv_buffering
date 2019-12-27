# Abel Sousa 2018

# PhD project
# Study of Control of Protein Abundance using Multi-omics Data from Cancer Samples


# -- regress-out batch effects from protein and rna abundance

# -- PCA the results



suppressMessages(library(tidyverse))
suppressMessages(library(gplots))
suppressMessages(library(gridExtra))

options(bitmapType = "cairo")
set.seed(123)



# load datasets


# samples with CNV, RNA and Protein measurements
common_samples <- read_tsv(file = "./files/cptac_cellLines_samples_prot_rna_cnv.txt", col_names = FALSE)


proteomics <- read_tsv("./files/proteomicsQ_cptac_cellLines.txt") %>%
    gather(key="sample", value="prot_log2FC", -gene) %>%
    filter(sample %in% common_samples$X1) %>%
    group_by(gene) %>%
    filter( sum(!is.na(prot_log2FC)) > n()*0.25 ) %>%
    mutate(prot_log2FC = scale(prot_log2FC)[,1]) %>%
    #mutate(prot_log2FC = (prot_log2FC-min(prot_log2FC, na.rm=T))/(max(prot_log2FC, na.rm=T)-min(prot_log2FC, na.rm=T))) %>%
    ungroup() %>%
    spread(key="sample", value="prot_log2FC") %>%
    as.data.frame()
rownames(proteomics) <- proteomics$gene
proteomics <- proteomics[, -c(1)]


rna <- read_tsv("./files/rna_tcga_cellLines.txt") %>%
    gather(key="sample", value="rna_log2CPM", -gene) %>%
    filter(sample %in% common_samples$X1) %>%
    group_by(gene) %>%
    mutate(rna_log2CPM = scale(rna_log2CPM)[,1]) %>%
    #mutate(rna_log2CPM = (rna_log2CPM-min(rna_log2CPM))/(max(rna_log2CPM)-min(rna_log2CPM))) %>%
    ungroup() %>%
    spread(key="sample", value="rna_log2CPM") %>%
    as.data.frame()
rownames(rna) <- rna$gene
rna <- rna[, -c(1)]


metadata <- read_tsv("./files/metadata_cptac_cellLines.txt") %>%
    filter(sample %in% common_samples$X1) %>%
    mutate(age = scale(age)[,1]) %>%
    #mutate(age = (age-min(age))/(max(age)-min(age))) %>%
    as.data.frame()
rownames(metadata) <- metadata$sample
metadata <- metadata[, -c(1)]


# First PCA
pca_before_rna <- prcomp(t(rna), center = TRUE, scale. = TRUE)$x %>% as.data.frame
pca_before_rna <- cbind.data.frame(pca_before_rna, batch = as.character(metadata[rownames(pca_before_rna), "batch"]))

pca_before_rna_var <- prcomp(t(rna), center = TRUE, scale. = TRUE)
pca_before_rna_var <- tibble(var = pca_before_rna_var$sdev^2, perc_var = (var/sum(var))*100) %>%
  mutate(pr_comp = fct_inorder(as.factor(paste("PC", seq(1:nrow(.)), sep=""))))


pca_before_protein <- prcomp(t(na.omit(proteomics)), center = TRUE, scale. = TRUE)$x %>% as.data.frame
pca_before_protein <- cbind.data.frame(pca_before_protein, batch = as.character(metadata[rownames(pca_before_protein), "batch"]))

pca_before_protein_var <- prcomp(t(na.omit(proteomics)), center = TRUE, scale. = TRUE)
pca_before_protein_var <- tibble(var = pca_before_protein_var$sdev^2, perc_var = (var/sum(var))*100) %>%
  mutate(pr_comp = fct_inorder(as.factor(paste("PC", seq(1:nrow(.)), sep=""))))


#scatterplots of PC1 vs PC2
pca_before_rna_plot <- ggplot( data=pca_before_rna, mapping=aes(x=PC1, y=PC2, color=batch) ) +
  geom_point() +
  theme_classic() +
  theme(
    axis.title=element_text(colour="black", size=15),
    axis.text=element_text(colour="black", size=12),
    legend.text=element_text(colour="black", size=10),
    legend.title=element_text(colour="black", size=15),
    legend.position = "bottom",
    aspect.ratio=1) +
  scale_color_discrete(labels=c("Lawrence et al", "Lapek et al", "Roumeliotis et al", "TCGA BRCA", "TCGA COREAD", "TCGA HGSC"), name="Batch") +
  guides(color = guide_legend(nrow = 3))
ggsave(filename="pca_before_rna.png", plot=pca_before_rna_plot, height=5, width=5, path = "./plots/pre_processing/")
ggsave(filename="pca_before_rna.pdf", plot=pca_before_rna_plot, height=5, width=5, path = "./plots/pre_processing/")
unlink("pca_before_rna.png")
unlink("pca_before_rna.pdf")


pca_before_protein_plot <- ggplot( data=pca_before_protein, mapping=aes(x=PC1, y=PC2, color=batch) ) +
  geom_point() +
  theme_classic() +
  theme(
    axis.title=element_text(colour="black", size=15),
    axis.text=element_text(colour="black", size=12),
    legend.text=element_text(colour="black", size=10),
    legend.title=element_text(colour="black", size=15),
    legend.position="bottom",
    aspect.ratio=1) +
  scale_color_discrete(labels=c("Lawrence et al", "Lapek et al", "Roumeliotis et al", "TCGA BRCA", "TCGA COREAD", "TCGA HGSC"), name="Batch") +
  guides(color = guide_legend(nrow = 3))
ggsave(filename="pca_before_protein.png", plot=pca_before_protein_plot, height=5, width=5, path = "./plots/pre_processing/")
ggsave(filename="pca_before_protein.pdf", plot=pca_before_protein_plot, height=5, width=5, path = "./plots/pre_processing/")
unlink("pca_before_protein.png")
unlink("pca_before_protein.pdf")


# correlation of principal components (1:10) with covariates
rna_before_cor <- cbind(pca_before_rna[, 1:10], metadata[rownames(pca_before_rna[, 1:10]), ])
rna_before_cor$gender <- as.numeric(as.factor(rna_before_cor$gender))
rna_before_cor$cancer <- as.numeric(as.factor(rna_before_cor$cancer))
rna_before_cor$batch <- as.numeric(as.factor(rna_before_cor$batch))
rna_before_cor$proteomics <- as.numeric(as.factor(rna_before_cor$proteomics))
rna_before_cor <- cor(rna_before_cor)


protein_before_cor <- cbind(pca_before_protein[, 1:10], metadata[rownames(pca_before_protein[, 1:10]), ])
protein_before_cor$gender <- as.numeric(as.factor(protein_before_cor$gender))
protein_before_cor$cancer <- as.numeric(as.factor(protein_before_cor$cancer))
protein_before_cor$batch <- as.numeric(as.factor(protein_before_cor$batch))
protein_before_cor$proteomics <- as.numeric(as.factor(protein_before_cor$proteomics))
protein_before_cor <- cor(protein_before_cor)


pdf(file="./plots/pre_processing/pca_before_rna_heatmap.pdf", height=6, width=6)
par(oma=c(1, 0, 0, 1))
heatmap.2( x = rna_before_cor, trace="none", Rowv=TRUE, Colv=TRUE, dendrogram="both", cexRow=1.2, cexCol = 1.2, key=TRUE, keysize=1, symkey=TRUE, key.title="Pearson's r", key.xlab = NA, key.ylab = NA, density.info = "none", col=bluered(100) )
dev.off()

pdf(file="./plots/pre_processing/pca_before_protein_heatmap.pdf", height=6, width=6)
par(oma=c(1, 0, 0, 1))
heatmap.2( x = protein_before_cor, trace="none", Rowv=TRUE, Colv=TRUE, dendrogram="both", cexRow=1.2, cexCol=1.2, key=TRUE, keysize=1, symkey=TRUE, key.title="Pearson's r", key.xlab = NA, key.ylab = NA, density.info = "none", col=bluered(100) )
dev.off()


# barplot of %var explained by principal component
pca_before_rna_var_plot <- ggplot( data=head(pca_before_rna_var, 10), mapping=aes(x=pr_comp, y=perc_var) ) +
  geom_bar(stat="identity") +
  theme_classic() +
  theme(plot.title=element_text(colour="black", size=15),
    axis.title=element_text(colour="black", size=15),
    axis.text=element_text(colour="black", size=10),
    aspect.ratio=1) +
  scale_y_continuous(limits = c(NA, 20)) +
  labs(x = "Principal components", y = "Proportion of variance explained (%)", title = "")
ggsave(filename="pca_before_rna_var_expl.png", plot=pca_before_rna_var_plot, height=5, width=5, path = "./plots/pre_processing/")
ggsave(filename="pca_before_rna_var_expl.pdf", plot=pca_before_rna_var_plot, height=5, width=5, path = "./plots/pre_processing/")
unlink("pca_before_rna_var_expl.png")
unlink("pca_before_rna_var_expl.pdf")

pca_before_protein_var_plot <- ggplot( data=head(pca_before_protein_var, 10), mapping=aes(x=pr_comp, y=perc_var) ) +
  geom_bar(stat="identity") +
  theme_classic() +
  theme(plot.title=element_text(colour="black", size=15),
    axis.title=element_text(colour="black", size=15),
    axis.text=element_text(colour="black", size=10),
    aspect.ratio=1) +
  scale_y_continuous(limits = c(NA, 8)) +
  labs(x = "Principal components", y = "Proportion of variance explained (%)", title = "")
ggsave(filename="pca_before_protein_var_expl.png", plot=pca_before_protein_var_plot, height=5, width=5, path = "./plots/pre_processing/")
ggsave(filename="pca_before_protein_var_expl.pdf", plot=pca_before_protein_var_plot, height=5, width=5, path = "./plots/pre_processing/")
unlink("pca_before_protein_var_expl.png")
unlink("pca_before_protein_var_expl.pdf")






# regress-out covariates

linear_model <- function(data){
  data <- as.data.frame(data)
  rownames(data) <- data$sample

  covars <- data[, -c(1)]

  # remove nas
  covars <- na.exclude(covars)


  # remove covariates with less than two levels
  covars <- covars[, apply(covars, 2, function(x) length(levels(as.factor(x))) >= 2), drop=F]


  # fit the model
  lreg_f <- as.formula( paste( "rna_log2CPM", "~", paste(colnames(covars[, -match("rna_log2CPM", colnames(covars)), drop=F]), collapse="+") ) )
  lreg <- lm(lreg_f, data = covars)
  lreg_r <- residuals(lreg)


  # model residuals
  resd <- as.tibble(data.frame(sample = names(lreg_r), rna_log2CPM = unname(lreg_r), stringsAsFactors = F))

  return(resd)

}


linear_model2 <- function(data){
  data <- as.data.frame(data)
  rownames(data) <- data$sample

  covars <- data[, -c(1)]

  # remove nas
  covars <- na.exclude(covars)


  # remove covariates with less than two levels
  covars <- covars[, apply(covars, 2, function(x) length(levels(as.factor(x))) >= 2), drop=F]


  # fit the model
  lreg_f <- as.formula( paste( "prot_log2FC", "~", paste(colnames(covars[, -match("prot_log2FC", colnames(covars)), drop=F]), collapse="+") ) )
  lreg <- lm(lreg_f, data = covars)
  lreg_r <- residuals(lreg)


  # model residuals
  resd <- as.tibble(data.frame(sample = names(lreg_r), prot_log2FC = unname(lreg_r), stringsAsFactors = F))

  return(resd)

}


metadata2 <- cbind.data.frame(sample = rownames(metadata), metadata)
metadata2 <- as.tibble(metadata2)


rna2 <- cbind.data.frame(gene = rownames(rna), rna)
rna2 <- rna2 %>%
  as.tibble() %>%
  gather(key="sample", value="rna_log2CPM", -gene) %>%
  inner_join(metadata2, by = "sample") %>%
  nest(-c(gene)) %>%
  mutate(rna_log2CPM_ = purrr::map(data, ~ linear_model(.))) %>%
  dplyr::select(-data) %>%
  unnest(rna_log2CPM_)


rna2 <- rna2 %>%
  spread(key="sample", value="rna_log2CPM") %>%
  as.data.frame()
rownames(rna2) <- rna2$gene
rna2 <- rna2[, -c(1)]


proteomics2 <- cbind.data.frame(gene = rownames(proteomics), proteomics)
proteomics2 <- proteomics2 %>%
  as.tibble() %>%
  gather(key="sample", value="prot_log2FC", -gene) %>%
  inner_join(metadata2, by = "sample") %>%
  nest(-c(gene)) %>%
  mutate(prot_log2FC_ = purrr::map(data, ~ linear_model2(.))) %>%
  dplyr::select(-data) %>%
  unnest(prot_log2FC_)


proteomics2 <- proteomics2 %>%
  spread(key="sample", value="prot_log2FC") %>%
  as.data.frame()
rownames(proteomics2) <- proteomics2$gene
proteomics2 <- proteomics2[, -c(1)]


#second PCA
pca_after_rna <- prcomp(t(rna2), center = TRUE, scale. = TRUE)$x %>% as.data.frame
pca_after_rna <- cbind.data.frame(pca_after_rna, batch = as.character(metadata[rownames(pca_after_rna), "batch"]))

pca_after_rna_var <- prcomp(t(rna2), center = TRUE, scale. = TRUE)
pca_after_rna_var <- tibble(var = pca_after_rna_var$sdev^2, perc_var = (var/sum(var))*100) %>%
  mutate(pr_comp = fct_inorder(as.factor(paste("PC", seq(1:nrow(.)), sep=""))))


pca_after_protein <- prcomp(t(na.omit(proteomics2)), center = TRUE, scale. = TRUE)$x %>% as.data.frame
pca_after_protein <- cbind.data.frame(pca_after_protein, batch = as.character(metadata[rownames(pca_after_protein), "batch"]))

pca_after_protein_var <- prcomp(t(na.omit(proteomics2)), center = TRUE, scale. = TRUE)
pca_after_protein_var <- tibble(var = pca_after_protein_var$sdev^2, perc_var = (var/sum(var))*100) %>%
  mutate(pr_comp = fct_inorder(as.factor(paste("PC", seq(1:nrow(.)), sep=""))))


#scatterplots of PC1 vs PC2
pca_after_rna_plot <- ggplot( data=pca_after_rna, mapping=aes(x=PC1, y=PC2, color=batch) ) +
  geom_point() +
  theme_classic() +
  theme(
    axis.title=element_text(colour="black", size=15),
    axis.text=element_text(colour="black", size=12),
    legend.text=element_text(colour="black", size=10),
    legend.title=element_text(colour="black", size=15),
    legend.position="bottom",
    aspect.ratio=1) +
  scale_color_discrete(labels=c("Lawrence et al", "Lapek et al", "Roumeliotis et al", "TCGA BRCA", "TCGA COREAD", "TCGA HGSC"), name="Batch") +
  guides(color = guide_legend(nrow = 3))
ggsave(filename="pca_after_rna.png", plot=pca_after_rna_plot, height=5, width=5, path = "./plots/pre_processing/")
ggsave(filename="pca_after_rna.pdf", plot=pca_after_rna_plot, height=5, width=5, path = "./plots/pre_processing/")
unlink("pca_after_rna.png")
unlink("pca_after_rna.pdf")


pca_after_rna_plot2 <- ggplot( data=pca_after_rna, mapping=aes(x=PC2, y=PC3, color=batch) ) +
  geom_point() +
  theme_classic() +
  theme(
    axis.title=element_text(colour="black", size=15),
    axis.text=element_text(colour="black", size=12),
    legend.text=element_text(colour="black", size=10),
    legend.title=element_text(colour="black", size=15),
    legend.position="bottom",
    aspect.ratio=1) +
  scale_color_discrete(labels=c("Lawrence et al", "Lapek et al", "Roumeliotis et al", "TCGA BRCA", "TCGA COREAD", "TCGA HGSC"), name="Batch") +
  guides(color = guide_legend(nrow = 3))
ggsave(filename="pca_after_rna2.png", plot=pca_after_rna_plot2, height=5, width=5, path = "./plots/pre_processing/")
unlink("pca_after_rna2.png")

pca_after_protein_plot <- ggplot( data=pca_after_protein, mapping=aes(x=PC1, y=PC2, color=batch) ) +
  geom_point() +
  theme_classic() +
  theme(
    axis.title=element_text(colour="black", size=15),
    axis.text=element_text(colour="black", size=12),
    legend.text=element_text(colour="black", size=10),
    legend.title=element_text(colour="black", size=15),
    legend.position="bottom",
    aspect.ratio=1) +
  scale_color_discrete(labels=c("Lawrence et al", "Lapek et al", "Roumeliotis et al", "TCGA BRCA", "TCGA COREAD", "TCGA HGSC"), name="Batch") +
  guides(color = guide_legend(nrow = 3))
ggsave(filename="pca_after_protein.png", plot=pca_after_protein_plot, height=5, width=5, path = "./plots/pre_processing/")
ggsave(filename="pca_after_protein.pdf", plot=pca_after_protein_plot, height=5, width=5, path = "./plots/pre_processing/")
unlink("pca_after_protein.png")
unlink("pca_after_protein.pdf")


# correlation of principal components (1:10) with covariates
rna_after_cor <- cbind(pca_after_rna[, 1:10], metadata[rownames(pca_after_rna[, 1:10]), ])
rna_after_cor$gender <- as.numeric(as.factor(rna_after_cor$gender))
rna_after_cor$cancer <- as.numeric(as.factor(rna_after_cor$cancer))
rna_after_cor$batch <- as.numeric(as.factor(rna_after_cor$batch))
rna_after_cor$proteomics <- as.numeric(as.factor(rna_after_cor$proteomics))
rna_after_cor <- cor(rna_after_cor)

protein_after_cor <- cbind(pca_after_protein[, 1:10], metadata[rownames(pca_after_protein[, 1:10]), ])
protein_after_cor$gender <- as.numeric(as.factor(protein_after_cor$gender))
protein_after_cor$cancer <- as.numeric(as.factor(protein_after_cor$cancer))
protein_after_cor$batch <- as.numeric(as.factor(protein_after_cor$batch))
protein_after_cor$proteomics <- as.numeric(as.factor(protein_after_cor$proteomics))
protein_after_cor <- cor(protein_after_cor)

pdf(file="./plots/pre_processing/pca_after_rna_heatmap.pdf", height=6, width=6)
par(oma=c(1, 0, 0, 1))
heatmap.2( x = rna_after_cor, trace="none", Rowv=TRUE, Colv=TRUE, dendrogram="both", cexRow=1.2, cexCol = 1.2, key=TRUE, keysize=1, symkey=TRUE, key.title="Pearson's r", key.xlab = NA, key.ylab = NA, density.info = "none", col=bluered(100) )
dev.off()

pdf(file="./plots/pre_processing/pca_after_protein_heatmap.pdf", height=6, width=6)
par(oma=c(1, 0, 0, 1))
heatmap.2( x = protein_after_cor, trace="none", Rowv=TRUE, Colv=TRUE, dendrogram="both", cexRow=1.2, cexCol=1.2, key=TRUE, keysize=1, symkey=TRUE, key.title="Pearson's r", key.xlab = NA, key.ylab = NA, density.info = "none", col=bluered(100) )
dev.off()


# barplot of %var explained by principal component
pca_after_rna_var_plot <- ggplot( data=head(pca_after_rna_var, 10), mapping=aes(x=pr_comp, y=perc_var) ) +
  geom_bar(stat="identity") +
  theme_classic() +
  theme(plot.title=element_text(colour="black", size=15),
    axis.title=element_text(colour="black", size=15),
    axis.text=element_text(colour="black", size=10),
    aspect.ratio=1) +
  scale_y_continuous(limits = c(NA, 20)) +
  labs(x = "Principal components", y = "Proportion of variance explained (%)", title = "")
ggsave(filename="pca_after_rna_var_expl.png", plot=pca_after_rna_var_plot, height=5, width=5, path = "./plots/pre_processing/")
ggsave(filename="pca_after_rna_var_expl.pdf", plot=pca_after_rna_var_plot, height=5, width=5, path = "./plots/pre_processing/")
unlink("pca_after_rna_var_expl.png")
unlink("pca_after_rna_var_expl.pdf")


pca_after_protein_var_plot <- ggplot( data=head(pca_after_protein_var, 10), mapping=aes(x=pr_comp, y=perc_var) ) +
  geom_bar(stat="identity") +
  theme_classic() +
  theme(plot.title=element_text(colour="black", size=15),
    axis.title=element_text(colour="black", size=15),
    axis.text=element_text(colour="black", size=10),
    aspect.ratio=1) +
  scale_y_continuous(limits = c(NA, 8)) +
  labs(x = "Principal components", y = "Proportion of variance explained (%)", title = "")
ggsave(filename="pca_after_protein_var_expl.png", plot=pca_after_protein_var_plot, height=5, width=5, path = "./plots/pre_processing/")
ggsave(filename="pca_after_protein_var_expl.pdf", plot=pca_after_protein_var_plot, height=5, width=5, path = "./plots/pre_processing/")
unlink("pca_after_protein_var_expl.png")
unlink("pca_after_protein_var_expl.pdf")




# export data

write.table(cbind.data.frame(gene = rownames(rna2), rna2), "./files/rna_tcga_cellLines_corrected.txt", sep="\t", quote=F, row.names=F)
write.table(cbind.data.frame(gene = rownames(proteomics2), proteomics2), "./files/proteomicsQ_cptac_cellLines_corrected.txt", sep="\t", quote=F, row.names=F)



save(list=ls(), file="./r_workspaces/regress_out_protein_rna.RData")
