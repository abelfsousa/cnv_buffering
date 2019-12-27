# Abel Sousa 2018

# PhD project
# Study of Control of Protein Abundance using Multi-omics Data from Cancer Samples

# -- Analysis of protein associations using hipsci project pQTLs







suppressMessages(library(tidyverse))
suppressMessages(library(ggpubr))

options(bitmapType = "cairo")



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


# load significative protein associations with CNV and RNA (filtered and not filtered by genomic co-localization)
protein_associations_filtered <- read_tsv("./files/signf_protein_associations_cnv_mrna_filtered_close_controlling.txt")
protein_associations_nonfiltered <- read_tsv("./files/signf_protein_associations_cnv_mrna.txt")


# load protein associations with CNV
protein_associations_cnv <- read_tsv("./files/protein_associations_tcga_ccle_cnv_reg.txt") %>%
  dplyr::select(controlling, controlled, beta_cnv, pval_cnv = lrt_p, fdr_cnv = fdr)


signf_protein_associations_cnv <- protein_associations_cnv %>%
  dplyr::filter(fdr_cnv < 0.05)


non_signf_protein_associations_cnv <- protein_associations_cnv %>%
  dplyr::filter(fdr_cnv > 0.05)




# load genome-wide protein associations with CNV
protein_associations_cnv_gw <- readRDS("./files/protein_associations_tcga_ccle_cnv_reg_genome_wide.rds") %>%
  dplyr::select(controlling, controlled, beta_cnv, pval_cnv = lrt_p, fdr_cnv = fdr)



# load hipsci pQTLs data
hipsci_pQTLs <- read_tsv("./data/hipsci/hipsci_cis_trans_genetic_associations.txt") %>%
  na.exclude() %>%
  dplyr::select("snp_id", "cis_gene_name", "trans_gene_name")




# merge with highly-significative protein associations
# FDR CNV < 1e-3
p_associations_hsignf <- signf_protein_associations_cnv %>%
  filter(fdr_cnv < 1e-3) %>%
  dplyr::select(controlling, controlled) %>%
  inner_join(hipsci_pQTLs, by=c("controlling" = "cis_gene_name", "controlled" = "trans_gene_name")) %>%
  mutate(signf = "highly-significant")
write.table(p_associations_hsignf, "./files/p_associations_hsignf_pQTLs_hipsci.txt", sep = "\t", quote = F, row.names = F)



# merge with significative protein associations
# 5e-2 >= FDR CNV >= 1e-3
p_associations_signf <- signf_protein_associations_cnv %>%
  filter(fdr_cnv >= 1e-3) %>%
  dplyr::select(controlling, controlled) %>%
  inner_join(hipsci_pQTLs, by=c("controlling" = "cis_gene_name", "controlled" = "trans_gene_name")) %>%
  mutate(signf = "significant")
write.table(p_associations_signf, "./files/p_associations_signf_pQTLs_hipsci.txt", sep = "\t", quote = F, row.names = F)



# merge with non-significative protein associations
# FDR CNV > 5e-2
p_associations_nonsignf <- non_signf_protein_associations_cnv %>%
  dplyr::select(controlling, controlled) %>%
  inner_join(hipsci_pQTLs, by=c("controlling" = "cis_gene_name", "controlled" = "trans_gene_name")) %>%
  mutate(signf = "non-significant")



#fisher test
contingency_table <- data.frame(QTL = c(3, 17), NQTL = c(1215, 410346), row.names=c("signf", "nsignf"))
fisher.test(contingency_table)$p.value #2.838365e-05
chisq.test(contingency_table)$p.value #9.37922e-24




# all protein pairs
biogrid_corum_ER_pairs <- read_tsv("./files/biogrid_corum_ER_pairs2.txt")

hipsci_pQTLs_physical <- hipsci_pQTLs %>%
  inner_join(biogrid_corum_ER_pairs, by=c("cis_gene_name" = "interactor.A", "trans_gene_name" = "interactor.B"))

hipsci_pQTLs_tested <- hipsci_pQTLs_physical %>%
  inner_join(protein_associations_cnv, by=c("cis_gene_name" = "controlling", "trans_gene_name" = "controlled"))

hipsci_pQTLs_signf <- hipsci_pQTLs_tested %>%
  inner_join(signf_protein_associations_cnv, by=c("cis_gene_name" = "controlling", "trans_gene_name" = "controlled", "beta_cnv", "pval_cnv", "fdr_cnv"))


hipsci_pQTLs_tested_pval <- hipsci_pQTLs_tested %>%
  mutate(pval_cnv = -log10(pval_cnv), pair = paste(cis_gene_name, trans_gene_name, sep="_")) %>%
  mutate(pair = as.factor(pair)) %>%
  mutate(pair = fct_reorder(pair, pval_cnv)) %>%
  ggplot(mapping = aes(x = pair, y = pval_cnv)) +
  geom_bar(stat="identity", fill = "grey", color="black") +
  geom_hline(yintercept=-log10(0.05), linetype="dashed") +
  theme_classic() +
  theme(plot.margin = margin(1, 0, 0, 1, "cm")) +
  theme(axis.title.y=element_text(colour="black", size=15),
        axis.title.x=element_text(colour="black", size=15),
        axis.text.y=element_text(colour="black", size=12),
        axis.text.x=element_text(colour="black", size=12, angle = 45, hjust = 1)) +
  labs(x = "Protein association", y = "CNV model P-value (-log10)")
ggsave("hipsci_pQTLs_tested_pval.png", plot = hipsci_pQTLs_tested_pval, path="./plots/hipsci_QTLs/")
unlink("hipsci_pQTLs_tested_pval.png")





hipsci_pQTLs_tested_gw <- hipsci_pQTLs %>%
  anti_join(biogrid_corum_ER_pairs, by=c("cis_gene_name" = "interactor.A", "trans_gene_name" = "interactor.B")) %>%
  inner_join(protein_associations_cnv_gw, by=c("cis_gene_name" = "controlling", "trans_gene_name" = "controlled"))

hipsci_pQTLs_tested_gw_pval <- hipsci_pQTLs_tested_gw %>%
  mutate(pval_cnv = -log10(pval_cnv), pair = paste(cis_gene_name, trans_gene_name, sep="_")) %>%
  mutate(pair = as.factor(pair)) %>%
  mutate(pair = fct_reorder(pair, pval_cnv)) %>%
  ggplot(mapping = aes(x = pair, y = pval_cnv)) +
  geom_bar(stat="identity", fill = "grey", color="black") +
  geom_hline(yintercept=-log10(0.05), linetype="dashed") +
  theme_classic() +
  theme(plot.margin = margin(1, 0, 0, 1, "cm")) +
  theme(axis.title.y=element_text(colour="black", size=15),
    axis.title.x=element_text(colour="black", size=15),
    axis.text.y=element_text(colour="black", size=12),
    axis.text.x=element_text(colour="black", size=7, angle = 45, hjust = 1)) +
  labs(x = "Protein association", y = "CNV model P-value (-log10)")
ggsave("hipsci_pQTLs_tested_gw_pval.png", plot = hipsci_pQTLs_tested_gw_pval, path="./plots/hipsci_QTLs/")
unlink("hipsci_pQTLs_tested_gw_pval.png")



save(list=ls(), file="./r_workspaces/protein_associations_QTLs_hipsci_samples.RData")
