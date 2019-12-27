# Abel Sousa 2018

# PhD project
# Study of Control of Protein Abundance using Multi-omics Data from Cancer Samples



# -- Mapping between cis-eQTLs and GWAS hits across attenuation states


suppressMessages(library(tidyverse))
suppressMessages(library(ggpubr))
suppressMessages(library(RColorBrewer))
suppressMessages(library(data.table))

options(bitmapType = "cairo")



# import protein attenuation data
protein_attenuation <- read_tsv("./files/protein_attenuation.txt") %>%
  mutate_at(vars(contains("class")), as.factor) %>%
  mutate_if(is.factor, fct_infreq)



# import GTEx cis-eQTLs
gtex_eqtls <- fread("./data/gtex/all_tissues.v6p.signif_snpgene_pairs.txt")
gtex_eqtls <- gtex_eqtls %>%
  as.tibble() %>%
  dplyr::rename(tissue = Adipose_Subcutaneous)


#gtex_eqtls_sel <- gtex_eqtls %>%
#  dplyr::select(variant_id, gene_id, tissue)


# import GENCODE v19 annotation
gene_annot_v19 <- read_tsv("./data/gtex/geneAnnot.gencode.v19.txt")


gtex_eqtls <- gtex_eqtls %>%
  inner_join(gene_annot_v19[, c("geneID", "geneName", "geneType")], by=c("gene_id" = "geneID")) %>%
  dplyr::rename(gene_name=geneName, gene_type=geneType)


# import GTEx cis-eQTLs rsids
gtex_rsids <- fread("./data/gtex/GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_635Ind_PASS_AB02_GQ20_HETX_MISS15_PLINKQC.lookup_table.txt")
gtex_rsids <- gtex_rsids %>%
  as.tibble()


# join rsids
gtex_eqtls <- gtex_eqtls %>%
  inner_join(gtex_rsids, by="variant_id")


# join attenuation class
attenuation_gtex_cis_eqtls <- gtex_eqtls %>%
  inner_join(protein_attenuation[, c("gene", "attenuation", "class")], by=c("gene_name" = "gene")) %>%
  dplyr::select(variant_id, rs_id_dbSNP147_GRCh37p13, chr, variant_pos, ref, alt, num_alt_per_site, gene_id, gene_name, gene_type, tissue, attenuation, attenuation_class = class, everything()) %>%
  dplyr::select(-gene_type)
saveRDS(attenuation_gtex_cis_eqtls, file = "./files/attenuation_gtex_cis_eqtls.rds")


# import GWAS catalog
gwas_catalog <- fread("./data/gwas_catalog/gwas_catalog_v1.0-associations_e94_r2018-09-30.tsv")
gwas_catalog <- gwas_catalog %>%
  as.tibble() %>%
  dplyr::select(gene_name = `REPORTED GENE(S)`, SNP_ID_CURRENT) %>%
  na.exclude() %>%
  filter(gene_name != "" & SNP_ID_CURRENT != "") %>%
  group_by(gene_name) %>%
  mutate(gene_name2 = str_split(gene_name, ", ")) %>%
  unnest() %>%
  ungroup() %>%
  dplyr::select(-gene_name) %>%
  dplyr::rename(gene_name = gene_name2) %>%
  distinct() %>%
  mutate(SNP_ID_CURRENT = paste0("rs", SNP_ID_CURRENT))




# merge GTEx eQTLs and GWAS catalog
attenuation_eqtls <- attenuation_gtex_cis_eqtls %>%
  #dplyr::select(-tissue) %>%
  #group_by(gene_name) %>%
  #distinct() %>%
  #ungroup() %>%
  group_by(attenuation_class) %>%
  #summarise(snps = length(unique(rs_id_dbSNP147_GRCh37p13))) %>%
  summarise(cis_eqtls = n()) %>%
  ungroup()


attenuation_eqtls_gwas <- attenuation_gtex_cis_eqtls %>%
  #dplyr::select(-tissue) %>%
  #group_by(gene_name) %>%
  #distinct() %>%
  #ungroup() %>%
  inner_join(gwas_catalog, by=c("rs_id_dbSNP147_GRCh37p13" = "SNP_ID_CURRENT", "gene_name")) %>%
  group_by(attenuation_class) %>%
  #summarise(cis_eqtls_gwas = length(unique(rs_id_dbSNP147_GRCh37p13))) %>%
  summarise(cis_eqtls_gwas = n()) %>%
  ungroup()


#
attenuation_eqtls_gwas_prop <- inner_join(attenuation_eqtls, attenuation_eqtls_gwas, by="attenuation_class") %>% mutate(perc = (cis_eqtls_gwas/cis_eqtls)*100)

contingency_table <- as.data.frame(attenuation_eqtls_gwas_prop[, c(2,3)])
rownames(contingency_table) = c("non-attenuated", "low-attenuated-protein", "high-attenuated-protein")

chisq.test(contingency_table)$p.value


attenuation_eqtls2 <- attenuation_gtex_cis_eqtls %>%
  #dplyr::select(-tissue) %>%
  #group_by(gene_name) %>%
  #distinct() %>%
  #ungroup() %>%
  group_by(gene_name, attenuation_class, attenuation) %>%
  #summarise(snps = length(unique(rs_id_dbSNP147_GRCh37p13))) %>%
  summarise(cis_eqtls = n()) %>%
  ungroup()


attenuation_eqtls_gwas2 <- attenuation_gtex_cis_eqtls %>%
  #dplyr::select(-tissue) %>%
  #group_by(gene_name) %>%
  #distinct() %>%
  #ungroup() %>%
  inner_join(gwas_catalog, by=c("rs_id_dbSNP147_GRCh37p13" = "SNP_ID_CURRENT", "gene_name")) %>%
  group_by(gene_name, attenuation_class, attenuation) %>%
  #summarise(snps_gwas = length(unique(rs_id_dbSNP147_GRCh37p13))) %>%
  summarise(cis_eqtls_gwas = n()) %>%
  ungroup()


attenuation_eqtls_gwas_prop2 <- full_join(attenuation_eqtls2, attenuation_eqtls_gwas2, by=c("gene_name", "attenuation_class", "attenuation")) %>%
  mutate(cis_eqtls_gwas = replace_na(cis_eqtls_gwas, 0)) %>%
  mutate(perc = (cis_eqtls_gwas/cis_eqtls)*100)


attenuation_eqtls_gwas2_boxplot <- ggplot(data = attenuation_eqtls_gwas_prop2 %>% filter(perc != 0), mapping = aes(x=attenuation_class, y=perc, fill=attenuation_class)) +
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(width=0.1, alpha=0.1) +
  theme_classic() +
  stat_compare_means(comparisons = list( c(1, 2), c(1, 3))) +
  theme(axis.title = element_text(colour="black", size=12),
    axis.text.y = element_text(colour="black", size=10),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.text = element_text(colour="black", size=10)) +
  scale_y_continuous(trans="log10") +
  labs(x = "Attenuation state", y = "Percentage of SNPs cis-eQTLs that are GWAS hits", fill = "")
ggsave(filename="attenuation_perc_gwas_snps.png", plot=attenuation_eqtls_gwas2_boxplot, width=5, height=5, path = "./plots/protein_associations_normal")
unlink("attenuation_perc_gwas_snps.png")


attenuation_eqtls_gwas2_scatter <- ggplot(data = attenuation_eqtls_gwas_prop2, mapping = aes(x=attenuation, y=perc)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_classic() +
  #stat_compare_means(comparisons = list( c(1, 3))) +
  theme(axis.title = element_text(colour="black", size=15),
      axis.text.y = element_text(colour="black", size=12),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()) +
  scale_y_continuous(trans="log2") +
  labs(x = "Attenuation", y = "Percentage of SNPs cis-eQTLs that are GWAS hits (log2)")
ggsave(filename="attenuation_perc_gwas_snps_scatter.png", plot=attenuation_eqtls_gwas2_scatter, path = "./plots/protein_associations_normal")
unlink("attenuation_perc_gwas_snps_scatter.png")
