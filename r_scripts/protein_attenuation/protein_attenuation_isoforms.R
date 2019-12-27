# Abel Sousa 2018

# PhD project
# Study of Control of Protein Abundance using Multi-omics Data from Cancer Samples



# -- Number of transcripts per attenuation class


suppressMessages(library(tidyverse))
suppressMessages(library(ggpubr))

options(bitmapType = "cairo")



# import protein attenuation data
protein_attenuation <- read_tsv("./files/protein_attenuation.txt")


# import gencode v28 annotation
gene_annotation <- read_tsv("./data/gencode/gencode.v28.annotation_gene_transcript.txt") %>%
  filter(gene_type == "protein_coding") %>%
  group_by(gene_name) %>%
  summarise(isoforms = n()) %>%
  ungroup()


attenuation_iso <- protein_attenuation %>%
  dplyr::select(gene, class) %>%
  inner_join(gene_annotation, by = c("gene" = "gene_name"))%>%
  ggplot(mapping = aes(x = class, y = isoforms)) +
  geom_boxplot(aes(fill=class)) +
  theme_classic() +
  stat_compare_means(comparisons = list( c(1,2), c(1,3), c(1,4), c(2,3) )) +
  theme(axis.title.y=element_text(colour="black", size=14),
    axis.title.x=element_text(colour="black", size=14),
    axis.text.y=element_text(colour="black", size=12),
    axis.text.x=element_blank()) +
  scale_x_discrete(name = "Attenuation class") +
  scale_y_continuous(name = "Number of isoforms per gene")
ggsave(filename="attenuation_isoforms.png", plot=attenuation_iso, width=5, height=5, path = "./plots/protein_attenuation_tcga_ccle/")
unlink("attenuation_isoforms.png")




gene_size <- read_tsv("./data/gencode/gencode.v28.annotation_gene_coordinates.txt") %>%
  filter(gene_type == "protein_coding") %>%
  mutate(width = end-(start-1))


attenuation_gene_size <- protein_attenuation %>%
  dplyr::select(gene, class) %>%
  inner_join(gene_size[, c("gene_name", "width")], by = c("gene" = "gene_name")) %>%
  ggplot(mapping = aes(x = class, y = width)) +
  geom_boxplot(aes(fill=class)) +
  theme_classic() +
  stat_compare_means(comparisons = list( c(1,2), c(1,3), c(1,4), c(2,3) )) +
  theme(axis.title.y=element_text(colour="black", size=14),
    axis.title.x=element_text(colour="black", size=14),
    axis.text.y=element_text(colour="black", size=12),
    axis.text.x=element_blank()) +
  scale_x_discrete(name = "Attenuation class") +
  scale_y_continuous(name = "Gene size")
ggsave(filename="attenuation_gene_size.png", plot=attenuation_gene_size, width=5, height=5, path = "./plots/protein_attenuation_tcga_ccle/")
unlink("attenuation_gene_size.png")


transcript_size <- read_tsv("./data/gencode/gencode.v28.annotation_transcript_coordinates.txt") %>%
  filter(gene_type == "protein_coding") %>%
  mutate(width = end-(start-1))


attenuation_transcript_size <- protein_attenuation %>%
  dplyr::select(gene, class) %>%
  inner_join(transcript_size[, c("gene_name", "width")], by = c("gene" = "gene_name")) %>%
  ggplot(mapping = aes(x = class, y = width)) +
  geom_boxplot(aes(fill=class)) +
  theme_classic() +
  #stat_compare_means(comparisons = list( c(1,2), c(1,3), c(1,4), c(2,3) )) +
  theme(axis.title.y=element_text(colour="black", size=14),
    axis.title.x=element_text(colour="black", size=14),
    axis.text.y=element_text(colour="black", size=12),
    axis.text.x=element_blank()) +
  scale_x_discrete(name = "Attenuation class") +
  scale_y_continuous(name = "Transcript size")
ggsave(filename="attenuation_transcript_size.png", plot=attenuation_transcript_size, width=5, height=5, path = "./plots/protein_attenuation_tcga_ccle/")
unlink("attenuation_transcript_size.png")
