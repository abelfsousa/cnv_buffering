# Abel Sousa 2018

# PhD project
# Study of Control of Protein Abundance using Multi-omics Data from Cancer Samples



# -- Ubiquitination sites fold-change upon proteasome inhibition by attenuation state


suppressMessages(library(tidyverse))
suppressMessages(library(ggpubr))
suppressMessages(library(viridis))
suppressMessages(library(RColorBrewer))

options(bitmapType = "cairo")




# import protein attenuation data
protein_attenuation <- read_tsv("./files/protein_attenuation.txt")




# import proteasome inhibition data
# Kim et al., 2011 (emanuel data)
protesome_inhibition <- read_csv("./data/emanuel/proteasome_inhibition.csv") %>%
  dplyr::select(gene = gene_symbol, site_position, bortezomib_2hrs_median, bortezomib_4hrs_median, bortezomib_8hrs_median) %>%
  gather(key="time", value="log2FC", -c(gene, site_position)) %>%
  filter(!is.na(log2FC))
  #group_by(gene, time) %>%
  #summarise(log2FC = mean(log2FC)) %>%
  #ungroup()




protesome_inhibition_boxpl <- protein_attenuation %>%
  dplyr::select(gene, class) %>%
  inner_join(protesome_inhibition, by="gene") %>%
  mutate_at(vars(contains("class")), as.factor) %>%
  mutate_if(is.factor, fct_infreq) %>%
  ggplot(mapping=aes(x=class, y=log2FC, fill=class)) +
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(fill="black", width=0.1, alpha=0.05) +
  theme_classic() +
  stat_compare_means(method = "wilcox.test", comparisons = list( c(1, 3)), label.y=5 ) +
  facet_wrap( ~ time, strip.position = "bottom", labeller=labeller(time = c(bortezomib_2hrs_median = "2 hours", bortezomib_4hrs_median = "4 hours", bortezomib_8hrs_median = "8 hours"))) +
  theme(axis.title.y=element_text(colour="black", size=15),
    axis.title.x=element_blank(),
    axis.text.y=element_text(colour="black", size=12),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    strip.text=element_text(colour="black", size=15),
    strip.background=element_blank(),
    legend.position="bottom",
    legend.text=element_text(size=10),
    legend.title=element_text(size=12)) +
  scale_fill_brewer(type = "seq", palette = "Blues", labels = c("non-attenuated", "lowly-attenuated", "highly-attenuated"), name = "Attenuation") +
  scale_y_continuous(limits = c(-5, 5)) +
  labs(y = "Ubiquitination sites (log2FC)")
ggsave("prot_attenuation_proteasome_inhibition.png", plot=protesome_inhibition_boxpl, path="./plots/protein_attenuation_tcga_ccle/", height=5, width=6)
unlink("prot_attenuation_proteasome_inhibition.png")




# map of uniprot/gene symbol/protein length
uniprot2gene <- read_tsv("./data/uniprot/uniprot_gene_plength_19_06_2018.tab") %>%
    dplyr::rename(UNIPROT=Entry, SYMBOL=`Gene names  (primary )`, LENGTH=Length) %>%
    dplyr::select(-Status, -Organism)



# import proteasome inhibition data
# Inigo compilation
proteasome_inhibition2 <- readRDS("./data/inigo/Sites_from_proteasomal_inhibition.rds") %>%
  as.tibble() %>%
  dplyr::select(Ratio, Entry, IDpos, log2) %>%
  filter(Entry != "") %>%
  filter(!grepl("NA", Ratio, fixed=T)) %>%
  mutate(log2 = as.numeric(log2)) %>%
  mutate(Ratio = replace(Ratio, grepl("Bortezomib 1um_vs_DMSO_time(min)_480", Ratio, fixed=T), "Bortezomib_1um_vs_DMSO time(min)_480")) %>%
  group_by(Ratio, Entry, IDpos) %>%
  summarise(log2 = median(log2)) %>%
  ungroup() %>%
  mutate(Ratio = replace(Ratio, grepl("Bortezomib 1um_vs_Control_time(min)_240_exp_372", Ratio, fixed=T), "Bortezomib_1um_vs_Control time(min)_240")) %>%
  mutate(Ratio = replace(Ratio, grepl("Bortezomib 1um_vs_DMSO_time(min)_120_exp_425", Ratio, fixed=T), "Bortezomib_1um_vs_DMSO time(min)_120")) %>%
  mutate(Ratio = replace(Ratio, grepl("Bortezomib 1um_vs_DMSO_time(min)_240_exp_425", Ratio, fixed=T), "Bortezomib_1um_vs_DMSO time(min)_240")) %>%
  mutate(Ratio = replace(Ratio, grepl("Epoxomicin 1um_vs_DMSO_time(min)_480_exp_423", Ratio, fixed=T), "Epoxomicin_1um_vs_DMSO time(min)_480")) %>%
  mutate(Ratio = replace(Ratio, grepl("Epoxomicin 20um_vs_Control_time(min)_480_exp_436", Ratio, fixed=T), "Epoxomicin_20um_vs_Control time(min)_480")) %>%
  mutate(Ratio = replace(Ratio, grepl("MG-132 10um _vs_Control_time(min)_240_exp_354", Ratio, fixed=T), "MG-132_10um_vs_Control time(min)_240")) %>%
  mutate(Ratio = replace(Ratio, grepl("MG-132 5um_vs_DMSO_time(min)_240_exp_364", Ratio, fixed=T), "MG-132_5um_vs_DMSO time(min)_240")) %>%
  mutate(Ratio = replace(Ratio, grepl("PR-619 5um_vs_DMSO_time(min)_240_exp_364", Ratio, fixed=T), "PR-619_5um_vs_DMSO time(min)_240")) %>%
  filter(!(Ratio %in% c("Bortezomib_1um_vs_Control time(min)_240", "PR-619_5um_vs_DMSO time(min)_240"))) %>%
  inner_join(uniprot2gene[, c("UNIPROT", "SYMBOL")], by=c("Entry" = "UNIPROT")) %>%
  inner_join(protein_attenuation[, c("gene", "class")], by=c("SYMBOL" = "gene")) %>%
  separate(Ratio, c("inhibitor", "time"), sep="[ ]") %>%
  #mutate(Ratio = sapply(Ratio, function(x) paste(strsplit(x, split=" ", fixed=T)[[1]], collapse="\n"))) %>%
  mutate_if(is.character, as.factor) %>%
  mutate(class = fct_infreq(class)) %>%
  dplyr::select(inhibitor, time, SYMBOL, class, Entry, IDpos, log2)




proteasome_inhibition2_distribution <- proteasome_inhibition2 %>%
  ggplot(mapping=aes(x=inhibitor, y=log2, fill=inhibitor)) +
  geom_boxplot() +
  geom_jitter(proteasome_inhibition2 %>% filter(SYMBOL == "PHB"), mapping=aes(x=inhibitor, y=log2), colour="red", width=0) +
  facet_wrap( ~ time, scales = "free_x", strip.position = "bottom", labeller=labeller(time = c(`time(min)_120` = "2 hours", `time(min)_240` = "4 hours", `time(min)_480` = "8 hours"))) +
  theme_classic() +
  theme(axis.title.y=element_text(colour="black", size=15),
    axis.title.x=element_blank(),
    axis.text.y=element_text(colour="black", size=12),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    strip.text=element_text(colour="black", size=15),
    strip.background=element_blank(),
    legend.position="right",
    legend.text=element_text(size=8),
    legend.title=element_text(size=10)) +
  scale_fill_viridis(discrete=T) +
  labs(y = "Ubiquitination sites (log2FC)")
ggsave("prot_attenuation_proteasome_inhibition_inigo_distribution.png", plot=proteasome_inhibition2_distribution, width=10, height=5, path="./plots/protein_attenuation_tcga_ccle/")
unlink("prot_attenuation_proteasome_inhibition_inigo_distribution.png")



proteasome_inhibition2_distribution_att <- proteasome_inhibition2 %>%
  mutate(inhibitor = fct_recode(inhibitor, `Bortezomib 1um vs DMSO` = "Bortezomib_1um_vs_DMSO",
  `Epoxomicin 1um vs DMSO` = "Epoxomicin_1um_vs_DMSO", `Epoxomicin 20um vs Control` = "Epoxomicin_20um_vs_Control",
  `MG-132 10um vs Control` = "MG-132_10um_vs_Control", `MG-132 5um vs DMSO` = "MG-132_5um_vs_DMSO")) %>%
  ggplot(mapping=aes(x=inhibitor, y=log2, fill=inhibitor, color=class)) +
  geom_boxplot(outlier.size=1) +
  #geom_jitter(proteasome_inhibition2 %>% filter(SYMBOL == "PHB"), mapping=aes(x=inhibitor, y=log2), colour="red", width=0, size=1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap( ~ time, scales = "free_x", strip.position = "bottom", labeller=labeller(time = c(`time(min)_120` = "2 hours", `time(min)_240` = "4 hours", `time(min)_480` = "8 hours"))) +
  #stat_compare_means(method = "wilcox.test", comparisons = list( c(1, 3)), label.y=5 ) +
  theme_classic() +
  theme(axis.title.y=element_text(colour="black", size=15),
    axis.title.x=element_blank(),
    axis.text.y=element_text(colour="black", size=12),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    strip.text=element_text(colour="black", size=15),
    strip.background=element_blank(),
    legend.position="bottom",
    legend.text=element_text(size=10),
    legend.title=element_text(size=12)) +
  scale_fill_viridis(discrete=T) +
  #scale_color_brewer(type = "seq", palette = "Blues", labels = c("non-attenuated", "lowly-attenuated", "highly-attenuated")) +
  scale_color_manual(values=c("#bdbdbd", "#969696", "#636363"), labels = c("non-attenuated", "lowly-attenuated", "highly-attenuated")) +
  guides(fill = guide_legend(ncol=2), color = guide_legend(ncol=1) ) +
  labs(y = "Ubiquitination sites (log2FC)", fill = "Condition", color="Attenuation")
ggsave("prot_attenuation_proteasome_inhibition_inigo_distribution_att.png", plot=proteasome_inhibition2_distribution_att, height=6, width=8.5, path="./plots/protein_attenuation_tcga_ccle/")
ggsave("prot_attenuation_proteasome_inhibition_inigo_distribution_att.pdf", plot=proteasome_inhibition2_distribution_att, height=6, width=8.5, path="./plots/protein_attenuation_tcga_ccle/")
unlink("prot_attenuation_proteasome_inhibition_inigo_distribution_att.png")
unlink("prot_attenuation_proteasome_inhibition_inigo_distribution_att.pdf")
