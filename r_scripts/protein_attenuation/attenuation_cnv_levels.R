# Abel Sousa 2018

# PhD project
# Study of Control of Protein Abundance using Multi-omics Data from Cancer Samples



# -- Exploration of CNV levels across attenuated genes



suppressMessages(library(tidyverse))
suppressMessages(library(ggpubr))
suppressMessages(library(RColorBrewer))
suppressMessages(library(data.table))

options(bitmapType = "cairo")

# import protein attenuation data
protein_attenuation <- read_tsv("./files/protein_attenuation.txt") %>%
  mutate_at(vars(contains("class")), as.factor) %>%
  mutate_if(is.factor, fct_infreq)


# import expression data
expression_data <- read_tsv("./files/protein_attenuation_cnv_rna_protein.txt")



# merge protein attenuation info with expression data
attenuation_expr <- protein_attenuation %>%
  dplyr::select(gene, attenuation, class) %>%
  inner_join(expression_data, by = "gene")


attenuation_expr %>%
  mutate(cnv_gistic2 = as.character(cnv_gistic2)) %>%
  filter(class == "high-attenuated-protein") %>%
  group_by(cnv_gistic2) %>%
  count() %>%
  ungroup() %>%
  mutate(freq = n/sum(n))

attenuation_expr %>%
  mutate(cnv_gistic2 = as.character(cnv_gistic2)) %>%
  filter(class == "high-attenuated-protein") %>%
  group_by(gene, cnv_gistic2) %>%
  summarise(counts = n()) %>%
  group_by(cnv_gistic2) %>%
  summarise(n_median_samples = median(counts))



foo <- attenuation_expr %>%
  mutate(cnv_gistic2 = as.character(cnv_gistic2)) %>%
  filter(class == "high-attenuated-protein") %>%
  group_by(gene, cnv_gistic2) %>%
  summarise(counts = n()) %>%
  arrange(desc(cnv_gistic2), desc(counts)) %>%
  ungroup() %>%
  mutate(class = "high-att")


foo2 <- attenuation_expr %>%
  mutate(cnv_gistic2 = as.character(cnv_gistic2)) %>%
  filter(class == "non-attenuated") %>%
  group_by(gene, cnv_gistic2) %>%
  summarise(counts = n()) %>%
  arrange(desc(cnv_gistic2), desc(counts)) %>%
  ungroup() %>%
  mutate(class = "non-att")


foo3 <- attenuation_expr %>%
  mutate(cnv_gistic2 = as.character(cnv_gistic2)) %>%
  filter(class == "low-attenuated-protein") %>%
  group_by(gene, cnv_gistic2) %>%
  summarise(counts = n()) %>%
  arrange(desc(cnv_gistic2), desc(counts)) %>%
  ungroup() %>%
  mutate(class = "low-att")

foo4 <- bind_rows(foo, foo2, foo3)


plot1 <- attenuation_expr %>%
  filter(gene %in% (foo %>% head(10) %>% pull(gene))) %>%
  dplyr::select(gene, sample, prot_log2FC, rna_log2CPM, cnv_gistic2) %>%
  gather(key = "measure", value = "expression", -c("gene", "sample", "cnv_gistic2")) %>%
  mutate_at(vars(contains("cnv_gistic2")), as.factor) %>%
  ggplot(mapping=aes(x=cnv_gistic2, y=expression, fill=measure)) +
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(mapping=aes(color=measure), alpha = 0.4, size=1, width = 0.1) +
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  facet_wrap( ~ gene, nrow=2, ncol=5) +
  theme_classic() +
  theme(
    axis.title=element_text(colour="black", size=15),
    axis.text=element_text(colour="black", size=12),
    strip.text=element_text(colour="black", size=13),
    strip.background=element_blank(),
    #legend.position="bottom",
    legend.text=element_text(size=12),
    legend.title=element_text(size=13)) +
  scale_fill_manual(values = c("black", "chocolate4"), labels = c("Protein", "RNA"), name = "Data") +
  scale_color_manual(values = c("black", "chocolate4"), labels = c("Protein", "RNA"), name = "Data") +
  labs(x = "CNV (gistic2)", y = "Expression (z-scores)")
ggsave("CNV_prot_rna_h_attenuated.png", plot=plot1, path="./plots/protein_attenuation_tcga_ccle/", height=10, width=15)
unlink("CNV_prot_rna_h_attenuated.png")

plot2 <- attenuation_expr %>%
  filter(gene %in% (foo2 %>% head(10) %>% pull(gene))) %>%
  dplyr::select(gene, sample, prot_log2FC, rna_log2CPM, cnv_gistic2) %>%
  gather(key = "measure", value = "expression", -c("gene", "sample", "cnv_gistic2")) %>%
  mutate_at(vars(contains("cnv_gistic2")), as.factor) %>%
  ggplot(mapping=aes(x=cnv_gistic2, y=expression, fill=measure)) +
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(mapping=aes(color=measure), alpha = 0.4, size=1, width = 0.1) +
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  facet_wrap( ~ gene, nrow=2, ncol=5) +
  theme_classic() +
  theme(
    axis.title=element_text(colour="black", size=15),
    axis.text=element_text(colour="black", size=12),
    strip.text=element_text(colour="black", size=13),
    strip.background=element_blank(),
    #legend.position="bottom",
    legend.text=element_text(size=12),
    legend.title=element_text(size=13)) +
  scale_fill_manual(values = c("black", "chocolate4"), labels = c("Protein", "RNA"), name = "Data") +
  scale_color_manual(values = c("black", "chocolate4"), labels = c("Protein", "RNA"), name = "Data") +
  labs(x = "CNV (gistic2)", y = "Expression (z-scores)")
ggsave("CNV_prot_rna_non_attenuated.png", plot=plot2, path="./plots/protein_attenuation_tcga_ccle/", height=10, width=15)
unlink("CNV_prot_rna_non_attenuated.png")


genes <- foo4 %>%
  filter(cnv_gistic2 == "-2" & (class == "high-att" | class == "non-att")) %>%
  group_by(class) %>%
  top_n(n=20, wt=counts) %>%
  pull(gene)


plot3 <- attenuation_expr %>%
  filter(gene %in% genes) %>%
  dplyr::select(gene, sample, class, prot_log2FC, rna_log2CPM, cnv_gistic2) %>%
  gather(key = "measure", value = "expression", -c("gene", "sample", "class", "cnv_gistic2")) %>%
  mutate_at(vars(contains("cnv_gistic2")), as.factor) %>%
  ggplot(mapping=aes(x=cnv_gistic2, y=expression)) +
  geom_boxplot(mapping=aes(fill=measure), outlier.shape=NA) +
  geom_jitter(mapping=aes(color=measure), alpha = 0.4, size=1, width = 0.1) +
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  facet_wrap( ~ class,
    labeller = labeller(class = c("high-attenuated-protein" = "Highly-attenuated", "non-attenuated" = "Non-attenuated"))) +
  theme_classic() +
  theme(
    axis.title=element_text(colour="black", size=15),
    axis.text=element_text(colour="black", size=13),
    strip.text=element_text(colour="black", size=15),
    strip.background=element_blank(),
    #legend.position="bottom",
    legend.text=element_text(size=12),
    legend.title=element_text(size=13)) +
  scale_fill_manual(values = c("black", "chocolate4"), labels = c("Protein", "RNA"), name = "Data") +
  scale_color_manual(values = c("black", "chocolate4"), labels = c("Protein", "RNA"), name = "Data") +
  labs(x = "CNV (gistic2)", y = "Protein/mRNA Expression (z-score standardized)")
ggsave("CNV_prot_rna_high_non_attenuated.png", plot=plot3, path="./plots/protein_attenuation_tcga_ccle/", height=8, width=12)
unlink("CNV_prot_rna_high_non_attenuated.png")




# load raw rna and proteomics
raw_rna <- fread("./files/tcga_cell_lines_rna_raw.txt") %>% as.tibble()
raw_prot <- fread("./files/cptac_cell_lines_proteomics_raw.txt") %>% as.tibble()



high_att <- attenuation_expr %>%
  dplyr::select(gene, attenuation, class, sample, cnv_gistic2) %>%
  filter(gene %in% (foo %>% filter(cnv_gistic2 == "-2" & counts > 10) %>% pull(gene))) %>%
  filter(cnv_gistic2 == "-2") %>%
  inner_join(raw_rna %>% dplyr::rename(rna = counts), by=c("gene", "sample")) %>%
  inner_join(raw_prot %>% dplyr::rename(protein = measure), by=c("gene", "sample", "cancer")) %>%
  arrange(rna) %>%
  dplyr::select(gene, sample, cancer, cnv_gistic2, rna, protein, proteomics) %>%
  as.data.frame()
write.table(high_att, "./files/high_att_21_genes_raw_rna_protein.txt", sep = "\t", quote = F, row.names=F)



low_att <- attenuation_expr %>%
  dplyr::select(gene, attenuation, class, sample, cnv_gistic2) %>%
  filter(gene %in% (foo3 %>% filter(cnv_gistic2 == "-2" & counts > 10) %>% pull(gene))) %>%
  filter(cnv_gistic2 == "-2") %>%
  inner_join(raw_rna %>% dplyr::rename(rna = counts), by=c("gene", "sample")) %>%
  inner_join(raw_prot %>% dplyr::rename(protein = measure), by=c("gene", "sample", "cancer")) %>%
  arrange(rna) %>%
  dplyr::select(gene, sample, cancer, cnv_gistic2, rna, protein, proteomics) %>%
  as.data.frame()
write.table(low_att, "./files/low_att_104_genes_raw_rna_protein.txt", sep = "\t", quote = F, row.names=F)



non_att <- attenuation_expr %>%
  dplyr::select(gene, attenuation, class, sample, cnv_gistic2) %>%
  filter(gene %in% (foo2 %>% filter(cnv_gistic2 == "-2" & counts > 10) %>% pull(gene))) %>%
  filter(cnv_gistic2 == "-2") %>%
  inner_join(raw_rna %>% dplyr::rename(rna = counts), by=c("gene", "sample")) %>%
  inner_join(raw_prot %>% dplyr::rename(protein = measure), by=c("gene", "sample", "cancer")) %>%
  arrange(rna) %>%
  dplyr::select(gene, sample, cancer, cnv_gistic2, rna, protein, proteomics) %>%
  as.data.frame()
write.table(non_att, "./files/non_att_142_genes_raw_rna_protein.txt", sep = "\t", quote = F, row.names=F)




all_genes <- attenuation_expr %>%
  dplyr::select(gene, attenuation, class, sample, cnv_gistic2) %>%
  #filter(gene %in% (foo %>% filter(cnv_gistic2 == "-2" & counts > 10) %>% pull(gene))) %>%
  #filter(cnv_gistic2 == "-2") %>%
  inner_join(raw_rna %>% dplyr::rename(rna = counts), by=c("gene", "sample")) %>%
  inner_join(raw_prot %>% dplyr::rename(protein = measure), by=c("gene", "sample", "cancer"))
  #arrange(rna) %>%
  #dplyr::select(gene, sample, cancer, cnv_gistic2, rna, protein, proteomics) %>%
  #as.data.frame()
