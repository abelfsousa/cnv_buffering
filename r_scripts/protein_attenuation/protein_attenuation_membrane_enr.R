# Abel Sousa 2018

# PhD project
# Study of Control of Protein Abundance using Multi-omics Data from Cancer Samples



# -- Enrichment of attenuated proteins on membrane proteins


suppressMessages(library(tidyverse))
suppressMessages(library(clusterProfiler))
suppressMessages(library(org.Hs.eg.db))
suppressMessages(library(igraph))
suppressMessages(library(viridis))

options(bitmapType = "cairo")



# import protein attenuation data
protein_attenuation <- read_tsv("./files/protein_attenuation.txt")


non_attenuated <- protein_attenuation %>% filter(class == "non-attenuated") %>% pull(gene)
low_attenuated <- protein_attenuation %>% filter(class == "low-attenuated-protein") %>% pull(gene)
high_attenuated <- protein_attenuation %>% filter(class == "high-attenuated-protein") %>% pull(gene)
attenuated <- protein_attenuation %>% filter(class == "low-attenuated-protein" | class == "high-attenuated-protein") %>% pull(gene)


# import GO CC terms

go_cc <- read.gmt("./data/go/c5.cc.v6.2.symbols.gmt") %>%
  as_tibble() %>%
  mutate(ont = str_replace(ont, "GO_", "")) %>%
  filter(ont == "INTRINSIC_COMPONENT_OF_PLASMA_MEMBRANE")


# fisher-test

contingency_table <- matrix(c(68, 303, 3367, 4386), nrow=2, ncol=2, dimnames = list(c("att", "no_att"), c("in_membrane", "not_in_membrane")))
contingency_table <- t(contingency_table)

fisher.test(contingency_table, alternative = "two.sided")$p.value
fisher.test(contingency_table, alternative = "greater")$p.value
fisher.test(contingency_table, alternative = "less")$p.value
