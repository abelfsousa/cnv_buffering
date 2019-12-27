# Abel Sousa 2018

# PhD project
# Study of Control of Protein Abundance using Multi-omics Data from Cancer Samples


# -- Pre-processing of protein half-lifes data


suppressMessages(library(tidyverse))
options(bitmapType = "cairo")


# load dataset
half_lifes <- read_tsv("./data/protein_turnover/protein_half_life_high_quality.txt", na = c("", "NA", "#N/A")) %>%
	dplyr::select("gene_name", matches("half_life")) %>%
	transmute(gene = gene_name,
		Bcells = rowMeans(data.frame(`Bcells replicate 1 half_life`, `Bcells replicate 2 half_life`), na.rm=TRUE),
		NKcells = rowMeans(data.frame(`NK cells replicate 1 half_life`, `NK cells replicate 2 half_life`), na.rm=TRUE),
		Hepatocytes = rowMeans(data.frame(`Hepatocytes replicate 1 half_life`, `Hepatocytes replicate 2 half_life`), na.rm=TRUE),
		Monocytes = rowMeans(data.frame(`Monocytes replicate 1 half_life`, `Monocytes replicate 2 half_life`), na.rm=TRUE),
		MouseNeurons = rowMeans(data.frame(`Mouse Neurons, replicate 3 half_life`, `Mouse Neurons, replicate 4 half_life`), na.rm=TRUE)) %>%
	mutate(Bcells = log10(Bcells), NKcells = log10(NKcells), Hepatocytes = log10(Hepatocytes), Monocytes = log10(Monocytes), MouseNeurons = log10(MouseNeurons))
write.table(half_lifes, "./files/protein_half_lifes.txt", sep="\t", quote=F, row.names=F)




# plot boxplot of complex type versus attenuation potential
half_lifes_boxplot <- ggplot(data=gather(half_lifes, key="cell", value="half_life", -gene)) +
    geom_boxplot(mapping=aes(x=cell, y=half_life, fill=cell), outlier.shape = NA) +
    geom_jitter(mapping=aes(x=cell, y=half_life), width=0.1, colour="black", alpha=0.1) +
    theme_classic() +
    theme(axis.title.y=element_text(colour="black", size=12),
          axis.title.x=element_text(colour="black", size=12),
          axis.text.y=element_text(colour="black", size=10),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) +
    labs(x = "Cell", y = "Protein Half-life (log10)", title = "", fill="Cell")
ggsave(filename="protein_half_lifes.png", plot=half_lifes_boxplot, path = "./plots/pre_processing/")
unlink("protein_half_lifes.png")





save(list=ls(), file = "./r_workspaces/preprocessing_protein_half_lifes.RData")



