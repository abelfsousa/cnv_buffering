# Abel Sousa 2018

# PhD project
# Study of Control of Protein Abundance using Multi-omics Data from Cancer Samples

# -- plots about number of samples by cancer type, protein associations, etc



suppressMessages(library(tidyverse))
suppressMessages(library(RColorBrewer))
suppressMessages(library(scales))
options(bitmapType = "cairo")



# load samples metadata
samples_metadata <- read_tsv("./files/metadata_cptac_cellLines.txt") %>%
    mutate_if(is.character, as.factor)


samples_metadata_summary <- samples_metadata %>%
    group_by(cancer, batch) %>%
    summarise(counts = n()) %>%
    ungroup()


sample_type_barplot <- ggplot(samples_metadata) +
    geom_bar(mapping = aes(fill=batch, y=..count.., x=cancer), position="stack") +
    theme(axis.title.y=element_text(colour="black", size=15), axis.title.x=element_text(colour="black", size=15)) +
    theme(axis.text.y=element_text(colour="black", size=12), axis.text.x=element_text(colour="black", size=12)) +
    scale_y_continuous(limits = c(0, 200)) +
    scale_fill_brewer(type = "qual", palette = "Set3") +
    labs(x = "Cancer type", y = "Number of samples", title = "")
ggsave(filename="tcga_ccle_total_samples_cnv.png", plot=sample_type_barplot, path = "./plots/summary_plots/")
unlink("tcga_ccle_total_samples_cnv.png")


# load common samples with protein, rna and cnv
samples_prot_rna_cnv <- read_tsv("./files/cptac_cellLines_samples_prot_rna_cnv.txt", col_names = FALSE) %>%
    rename(sample=X1) %>%
    inner_join(samples_metadata, by="sample") %>%
    mutate_if(is.character, as.factor) %>%
    mutate(batch = fct_relevel(batch, "LAW", "LPK", "TCGA_BRCA", "RMLT", "TCGA_COREAD", "TCGA_OV"))


sample_type_barplot2 <- ggplot(samples_prot_rna_cnv %>% group_by(cancer, batch) %>% summarise(counts = n())) +
    geom_bar(mapping = aes(fill=batch, y=counts, x=cancer), position="dodge", stat="identity") +
    theme_classic() +
    theme(
        axis.title=element_text(colour="black", size=15),
        axis.text=element_text(colour="black", size=12),
        plot.title = element_blank(),
        legend.text=element_text(size=10),
        legend.title=element_text(size=12)) +
    scale_y_continuous(limits = c(0, 150), name="Number of samples") +
    scale_x_discrete(labels=c("Breast", "Colorectal", "Ovarian"), name="Cancer") +
    scale_fill_brewer(type = "div", palette = "RdYlGn", labels=c("Lawrence et al", "Lapek et al", "TCGA BRCA", "Roumeliotis et al", "TCGA COREAD", "TCGA HGSC"), direction=-1, name="Batch")
    #scale_fill_brewer(type = "qual", palette = "Set3", labels=c("Lawrence et al", "Lapek et al", "TCGA BRCA", "Roumeliotis et al", "TCGA COREAD", "TCGA HGSC"), name="Samples")
ggsave(filename="tcga_ccle_common_samples_prot_rna_cnv.png", plot=sample_type_barplot2, height=5, width=5, path = "./plots/summary_plots/")
ggsave(filename="tcga_ccle_common_samples_prot_rna_cnv.pdf", plot=sample_type_barplot2, height=5, width=5, path = "./plots/summary_plots/")
unlink("tcga_ccle_common_samples_prot_rna_cnv.png")
unlink("tcga_ccle_common_samples_prot_rna_cnv.pdf")


# load common samples with protein, phospho, rna and cnv
samples_prot_phospho_rna_cnv <- read_tsv("./files/cptac_cellLines_samples_prot_phospho_rna_cnv.txt", col_names = FALSE) %>%
    rename(sample=X1) %>%
    inner_join(samples_metadata, by="sample") %>%
    mutate_if(is.character, as.factor) %>%
    mutate(batch = fct_drop(batch)) %>%
    mutate(batch = fct_relevel(batch, "TCGA_BRCA", "RMLT", "TCGA_OV"))

sample_type_barplot3 <- ggplot(samples_prot_phospho_rna_cnv) +
    geom_bar(mapping = aes(fill=batch, y=..count.., x=cancer), position="dodge") +
    theme_classic() +
    coord_flip() +
    theme(axis.title=element_text(colour="black", size=20),
          axis.text.x=element_text(colour="black", size=16),
          axis.text.y=element_text(colour="black", size=16),
          legend.text=element_text(colour="black", size=16),
          legend.title=element_text(colour="black", size=20),
          plot.title=element_blank(),
          legend.position = "bottom") +
    scale_y_continuous(limits = c(0, 80), name="Number of samples") +
    scale_x_discrete(labels=c("Breast", "Colorectal", "Ovarian"), name="Cancer") +
    scale_fill_manual(values=c("#1a9850", "#fc8d59", "#d73027"), labels=c("TCGA BRCA", "Roumeliotis et al", "TCGA HGSC"), name="Batch") +
    labs(x = "", y = "Number of samples")
ggsave(filename="tcga_ccle_common_samples_prot_phospho_rna_cnv.png", plot=sample_type_barplot3, height=5, width=10, path = "./plots/summary_plots/")
unlink("tcga_ccle_common_samples_prot_phospho_rna_cnv.png")


sample_type_barplot3 <- ggplot(samples_prot_phospho_rna_cnv) +
    geom_bar(mapping = aes(fill=batch, y=..count.., x=cancer), position="dodge") +
    theme_classic() +
    theme(axis.title=element_text(colour="black", size=17),
          axis.text=element_text(colour="black", size=14),
          legend.text=element_text(colour="black", size=14),
          legend.title=element_text(colour="black", size=17),
          plot.title=element_blank(),
          legend.position = "right") +
    scale_y_continuous(limits = c(0, 80), name="Number of samples") +
    scale_x_discrete(labels=c("Breast", "Colorectal", "Ovarian"), name="Cancer") +
    scale_fill_manual(values=c("#1a9850", "#fc8d59", "#d73027"), labels=c("TCGA BRCA", "Roumeliotis et al", "TCGA HGSC"), name="Batch") +
    labs(x = "", y = "Number of samples")
ggsave(filename="tcga_ccle_common_samples_prot_phospho_rna_cnv.png", plot=sample_type_barplot3, height=7, width=6, path = "./plots/summary_plots/")
ggsave(filename="tcga_ccle_common_samples_prot_phospho_rna_cnv.pdf", plot=sample_type_barplot3, height=7, width=6, path = "./plots/summary_plots/")
unlink("tcga_ccle_common_samples_prot_phospho_rna_cnv.png")
unlink("tcga_ccle_common_samples_prot_phospho_rna_cnv.pdf")




# load BioGRID + CORUM + ER membrane protein interactions
protein_interactions <- read_tsv("./files/biogrid_corum_ER_pairs.txt") %>%
  mutate(experiment.type = fct_infreq(as.factor(experiment.type)))



# load protein associations with CNV
protein_associations_cnv <- read_tsv("./files/protein_associations_tcga_ccle_cnv_reg.txt") %>%
    mutate(pair = paste(controlling, controlled, sep="_")) %>%
    dplyr::select(pair, everything())


# load protein associations with RNA
protein_associations_rna <- read_tsv("./files/protein_associations_tcga_ccle_mRNA_reg.txt") %>%
      mutate(pair = paste(controlling, controlled, sep="_")) %>%
      dplyr::select(pair, everything())


# load protein associations with Phospho
protein_associations_phospho <- read_tsv("./files/protein_associations_tcga_ccle_phospho_reg.txt") %>%
    mutate(pair = paste(controlling, controlled, sep="_")) %>%
    dplyr::select(pair, everything())



protein_interactions_system <- ggplot(protein_interactions %>% filter(experiment.type != "curated_membrane_complexes")) +
    geom_bar(mapping = aes(x=experiment.type, y=(..count..)/sum(..count..)*100)) +
    theme_classic() +
    theme(axis.title.y=element_text(colour="black", size=15),
      axis.title.x=element_text(colour="black", size=15),
      axis.text.y=element_text(colour="black", size=13),
      axis.text.x=element_text(colour="black", size=10, angle = 45, hjust = 1)) +
    scale_y_continuous(limits = c(0, 50)) +
    theme(plot.margin = margin(1, 1, 0, 2, "cm")) +
    labs(x = "Experimental System", y = "Percentage of interactions")
ggsave(filename="biogrid_corum_curated_complexes_interactions.png", plot=protein_interactions_system, path = "./plots/protein_interactions/")
unlink("biogrid_corum_curated_complexes_interactions.png")



protein_interactions_counts <- protein_interactions %>%
    mutate(experiment.type = str_replace(experiment.type, "^((?!corum|curated_membrane_complexes).)*$", "BioGRID")) %>%
    group_by(experiment.type) %>%
    distinct() %>%
    summarise(counts = n()) %>%
    add_row(experiment.type = "total_interactions", counts = nrow(unique(protein_interactions[,c(1,2)]))) %>%
    add_row(experiment.type = "tested_cnv", counts = nrow(protein_associations_cnv)) %>%
    add_row(experiment.type = "tested_phospho", counts = nrow(protein_associations_phospho)) %>%
    filter(experiment.type != "curated_membrane_complexes") %>%
    mutate_if(is.character, as.factor) %>%
    mutate(experiment.type = fct_relevel(experiment.type, "BioGRID", "corum", "total_interactions")) %>%
    mutate(experiment.type = fct_rev(experiment.type))


protein_interactions_counts_barplot <- ggplot(data=protein_interactions_counts, mapping=aes(x=experiment.type, y=counts)) +
    geom_bar(stat="identity", position="dodge", fill="darkcyan") +
    theme_classic() +
    theme(
        axis.title.x=element_text(colour="black", size=20),
        axis.text.y=element_text(colour="black", size=20),
        axis.text.x=element_text(colour="black", size=17),
        plot.title = element_text(size=25, hjust = 0.5)) +
    #scale_fill_brewer(type="qual", palette="Set2", guide=F) +
    scale_y_continuous(name="Number of protein interactions", labels = scales::comma) +
    scale_x_discrete(labels = rev(c("BioGRID", "CORUM", "Total interactions", "Tested interactions\nwith CNV model", "Tested interactions\nwith Phospho model")), name="") +
    labs(title = "Protein interactions") +
    theme(plot.margin = margin(0.5, 2, 0.5, 0, "cm")) +
    coord_flip()
ggsave(filename="protein_interactions_total_number.png", plot=protein_interactions_counts_barplot,  width = 12, path = "./plots/protein_interactions/")
unlink("protein_interactions_total_number.png")






save(list=ls(), file="./r_workspaces/summary_plots.RData")
