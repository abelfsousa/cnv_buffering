# Abel Sousa 2018

# PhD project
# Study of Control of Protein Abundance using Multi-omics Data from Cancer Samples

# -- Analysis of protein associations using normal samples from mouse







suppressMessages(library(tidyverse))
suppressMessages(library(ggpubr))
suppressMessages(library(data.table))
suppressMessages(library(VennDiagram))

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



# load human mouse orthologous
# http://www.ensembl.info/2009/01/21/how-to-get-all-the-orthologous-genes-between-two-species/#comments
# https://www.ensembl.org/Help/View?id=135
ortho_human_mouse <-  read_tsv("./data/ensembl/mart_export-human-mouse-orthologous.txt") %>%
  filter(Mouse_homology_type == "ortholog_one2one")



# load mouse eQTLs/pQTLs data from paper nature18270
# select the mouse genes with human orthologous
mice_eQTLs <- fread("./data/mouse_data/nature18270-s1/DO_mice_eQTLs.txt") %>%
  as.tibble() %>%
  dplyr::select(Gene_Symbol, Gene_Chromosome, Gene_Start_Mbp, Gene_End_Mbp, eQTL_LocusID, eQTL_Chromosome, QTL_Type) %>%
  inner_join(ortho_human_mouse[, c("Gene_name", "Mouse_gene_name")], by=c("Gene_Symbol" = "Mouse_gene_name")) %>%
  dplyr::select(Gene_mouse = Gene_Symbol, Gene_human = Gene_name, everything())

mice_pQTLs <- fread("./data/mouse_data/nature18270-s1/DO_mice_pQTLs.txt") %>%
  as.tibble() %>%
  dplyr::select(Gene_Symbol, Gene_Chromosome, Gene_Start_Mbp, Gene_End_Mbp, pQTL_LocusID, pQTL_Chromosome, QTL_Type) %>%
  inner_join(ortho_human_mouse[, c("Gene_name", "Mouse_gene_name")], by=c("Gene_Symbol" = "Mouse_gene_name")) %>%
  dplyr::select(Gene_mouse = Gene_Symbol, Gene_human = Gene_name, everything())




# highly-significative protein associations
# FDR CNV < 1e-3
p_associations_hsignf <- protein_associations_rna_cnv %>%
  filter(fdr_cnv < 1e-3) %>%
  dplyr::select(controlling, controlled) %>%
  inner_join(mice_pQTLs[, c("Gene_human", "pQTL_LocusID", "QTL_Type")], by=c("controlling" = "Gene_human")) %>%
  dplyr::rename(pQTL_controlling = QTL_Type) %>%
  inner_join(mice_pQTLs[, c("Gene_human", "pQTL_LocusID", "QTL_Type")], by=c("controlled" = "Gene_human", "pQTL_LocusID")) %>%
  dplyr::rename(pQTL_controlled = QTL_Type) %>%
  mutate(signf = "highly-significant") %>%
  filter(pQTL_controlling == "LOCAL", pQTL_controlled == "DISTANT")
write.table(p_associations_hsignf, "./files/p_associations_hsignf_pQTLs_mice.txt", sep = "\t", quote = F, row.names = F)



# merge with significative protein associations
# 5e-2 >= FDR CNV >= 1e-3
p_associations_signf <- protein_associations_rna_cnv %>%
  filter(fdr_cnv >= 1e-3) %>%
  dplyr::select(controlling, controlled) %>%
  inner_join(mice_pQTLs[, c("Gene_human", "pQTL_LocusID", "QTL_Type")], by=c("controlling" = "Gene_human")) %>%
  dplyr::rename(pQTL_controlling = QTL_Type) %>%
  inner_join(mice_pQTLs[, c("Gene_human", "pQTL_LocusID", "QTL_Type")], by=c("controlled" = "Gene_human", "pQTL_LocusID")) %>%
  dplyr::rename(pQTL_controlled = QTL_Type) %>%
  mutate(signf = "significant") %>%
  filter(pQTL_controlling == "LOCAL", pQTL_controlled == "DISTANT")
write.table(p_associations_signf, "./files/p_associations_signf_pQTLs_mice.txt", sep = "\t", quote = F, row.names = F)



# merge with non-significative protein associations
# FDR CNV > 5e-2
# sample 516 associations
p_associations_nonsignf <- non_signf_protein_associations_cnv_rna %>%
  dplyr::select(controlling, controlled) %>%
  inner_join(mice_pQTLs[, c("Gene_human", "pQTL_LocusID", "QTL_Type")], by=c("controlling" = "Gene_human")) %>%
  dplyr::rename(pQTL_controlling = QTL_Type) %>%
  inner_join(mice_pQTLs[, c("Gene_human", "pQTL_LocusID", "QTL_Type")], by=c("controlled" = "Gene_human", "pQTL_LocusID")) %>%
  dplyr::rename(pQTL_controlled = QTL_Type) %>%
  mutate(signf = "non-significant") %>%
  filter(pQTL_controlling == "LOCAL", pQTL_controlled == "DISTANT")



#fisher test
contingency_table <- data.frame(QTL = c(2, 4), NQTL = c(514, 357544), row.names=c("signf", "nsignf"))
fisher.test(t(contingency_table))$p.value #5.88781e-06
chisq.test(t(contingency_table))$p.value #3.604731e-40



# Venn Diagram
venn <- draw.pairwise.venn(516, 26, 5,
    category = c("", ""),
    #lty = rep("blank", 2),
    fill = c("deepskyblue", "red"),
    alpha = rep(0.5, 2),
    cex = rep(2, 3),
    cat.cex = rep(1.5, 2),
    ext.text=TRUE,
    ext.pos = c(0, 0),
    ext.dist = c(-0.1, 0.03),
    ext.line.lty = "solid",
    ext.line.lwd = 1)
ggsave(filename="venn_diagram.png", plot=venn, path = "./plots/mouse_QTLs/", width=3, height=3)
unlink("venn_diagram.png")


save(list=ls(), file="./r_workspaces/protein_associations_QTLs_mouse_samples.RData")
