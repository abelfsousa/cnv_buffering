# Abel Sousa 2018

# PhD project
# Study of Control of Protein Abundance using Multi-omics Data from Cancer Samples


# -- Protein attenuation by control status


suppressMessages(library(tidyverse))
suppressMessages(library(ggpubr))
suppressMessages(library(RColorBrewer))
suppressMessages(library(gridExtra))
options(bitmapType = "cairo")




# load datasets




# corr(CNV, RNA) and corr(CNV, Protein)
cnv_rna_protein <- read_tsv("./files/protein_attenuation.txt")



# import protein associations with CNV and RNA
# filtered by genomic co-localization
protein_associations_cnv_mrna <- read_tsv("./files/signf_protein_associations_cnv_mrna_filtered_close_controlling.txt") %>%
    mutate(pair = paste(controlling, controlled, sep="_")) %>%
    dplyr::select(pair, everything())



# significative protein associations by fdr and cnv beta value
signf_protein_associations_cnv_mrna <- cbind(protein_associations_cnv_mrna[, c("controlling", "controlled")], data.frame(signf = "fdr < 0.05", stringsAsFactors=F))
signf_protein_associations_cnv_mrna <- signf_protein_associations_cnv_mrna %>%
    bind_rows( cbind(protein_associations_cnv_mrna[protein_associations_cnv_mrna$fdr_cnv < 0.01, c("controlling", "controlled")], data.frame(signf = "fdr < 0.01", stringsAsFactors=F)) ) %>%
    bind_rows( cbind(protein_associations_cnv_mrna[protein_associations_cnv_mrna$fdr_cnv < 0.001, c("controlling", "controlled")], data.frame(signf = "fdr < 0.001", stringsAsFactors=F)) ) %>%
    bind_rows( cbind(protein_associations_cnv_mrna[protein_associations_cnv_mrna$fdr_cnv < 0.05 & abs(protein_associations_cnv_mrna$beta_cnv) > 0.5, c("controlling", "controlled")], data.frame(signf = "fdr < 0.05 & cnv |beta| > 0.5", stringsAsFactors=F)) ) %>%
    bind_rows( cbind(protein_associations_cnv_mrna[protein_associations_cnv_mrna$fdr_cnv < 0.01 & abs(protein_associations_cnv_mrna$beta_cnv) > 0.5, c("controlling", "controlled")], data.frame(signf = "fdr < 0.01 & cnv |beta| > 0.5", stringsAsFactors=F)) ) %>%
    bind_rows( cbind(protein_associations_cnv_mrna[protein_associations_cnv_mrna$fdr_cnv < 0.001 & abs(protein_associations_cnv_mrna$beta_cnv) > 0.5, c("controlling", "controlled")], data.frame(signf = "fdr < 0.001 & cnv |beta| > 0.5", stringsAsFactors=F)) ) %>%
    as.tibble()



# control status
protein_list_control <- signf_protein_associations_cnv_mrna %>%
    group_by(signf) %>%
    mutate(gene = list(c(controlling, controlled))) %>%
    unnest() %>%
    mutate(control_status = if_else(gene %in% setdiff(controlling, controlled), "controlling", if_else(gene %in% setdiff(controlled, controlling), "controlled", "both"))) %>%
    dplyr::select(-c(controlling, controlled)) %>%
    group_by(signf, gene) %>%
    distinct() %>%
    ungroup()



# number of times each protein is found
protein_list_interactions <- signf_protein_associations_cnv_mrna %>%
    group_by(signf) %>%
    mutate(gene = list(c(controlling, controlled))) %>%
    dplyr::select(signf, gene) %>%
    unique() %>%
    unnest() %>%
    group_by(signf, gene) %>%
    summarise(interactions = n()) %>%
    ungroup()



# join control status and number of interactions
# add attenuation status
protein_list <- inner_join(protein_list_control, protein_list_interactions, by=c("signf", "gene"))
protein_list <- inner_join(protein_list, cnv_rna_protein, by="gene") %>%
    dplyr::select(-c(r_cnv_rna, r_cnv_prot, p_cnv_rna, p_cnv_prot)) %>%
    mutate(signf = factor(signf, levels = c("fdr < 0.05", "fdr < 0.01", "fdr < 0.001", "fdr < 0.05 & cnv |beta| > 0.5", "fdr < 0.01 & cnv |beta| > 0.5", "fdr < 0.001 & cnv |beta| > 0.5")))


# boxplots of control status by attenuation potential
protein_attenuation_control_status1 <- ggplot(data = protein_list[protein_list$signf == "fdr < 0.05", ], mapping = aes(x = control_status, y = attenuation, fill = control_status)) +
    geom_boxplot(outlier.shape=NA) +
    geom_jitter(width=0.1, alpha=0.1) +
    stat_compare_means(comparisons = list( c("controlling", "controlled") ) ) +
    #stat_compare_means(comparisons = list( c("both", "controlled"), c("controlled", "controlling"), c("both", "controlling") )) +
    theme_classic() +
    coord_flip() +
    #facet_grid(. ~ signf, scales = "free") +
    theme(axis.title=element_text(colour="black", size=12),
          axis.text.x=element_text(colour="black", size=10),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.line.y = element_blank(),
          legend.text=element_text(colour="black", size=10),
          legend.title=element_text(colour="black", size=12),
          plot.title=element_blank(),
          legend.position = "none") +
    ggpubr::fill_palette("jco") +
    labs(x = "", y = "Attenuation potential", fill = "Control status")
ggsave(filename="protein_attenuation_control_status_allbetas.png", plot=protein_attenuation_control_status1, height=2, width=5, path = "./plots/protein_attenuation_tcga_ccle/")
unlink("protein_attenuation_control_status_allbetas.png")



protein_attenuation_control_status1 <- ggplot(data = protein_list[protein_list$signf == "fdr < 0.05", ], mapping = aes(x = control_status, y = attenuation, fill = control_status)) +
    geom_boxplot(outlier.shape=NA) +
    geom_jitter(width=0.1, alpha=0.1) +
    stat_compare_means(comparisons = list( c("controlling", "controlled") ) ) +
    #stat_compare_means(comparisons = list( c("both", "controlled"), c("controlled", "controlling"), c("both", "controlling") )) +
    theme_classic() +
    #facet_grid(. ~ signf, scales = "free") +
    theme(axis.title=element_text(colour="black", size=15),
      axis.text.y=element_text(colour="black", size=12),
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      legend.position="none") +
    ggpubr::fill_palette("jco") +
    labs(x = "Control status", y = "Attenuation potential", fill = "")
ggsave(filename="protein_attenuation_control_status_allbetas.png", plot=protein_attenuation_control_status1, height=4, width=3, path = "./plots/protein_attenuation_tcga_ccle/")
ggsave(filename="protein_attenuation_control_status_allbetas.pdf", plot=protein_attenuation_control_status1, height=4, width=3, path = "./plots/protein_attenuation_tcga_ccle/")
unlink("protein_attenuation_control_status_allbetas.png")
unlink("protein_attenuation_control_status_allbetas.pdf")


protein_attenuation_control_status2 <- ggplot(data = protein_list[protein_list$signf == "fdr < 0.05 & cnv |beta| > 0.5", ], mapping = aes(x = control_status, y = attenuation, fill = control_status)) +
    geom_boxplot(outlier.shape=NA) +
    geom_jitter(width=0.1, alpha=0.5) +
    stat_compare_means(comparisons = list(c("controlled", "controlling"))) +
    theme_classic() +
    #facet_grid(. ~ signf, scales = "free") +
    theme(axis.title.y=element_text(colour="black", size=14),
          axis.title.x=element_text(colour="black", size=14),
          axis.text.y=element_text(colour="black", size=13),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          legend.text=element_text(size=13),
          plot.title=element_text(size=13)) +
    labs(x = "Control status", y = "Attenuation potential\ncorr(CNV, RNA) - corr(CNV, Protein)", fill="", title="Control status and attenuation potential\n(FDR < 5% and CNV |beta| > 0.5)")
ggsave(filename="protein_attenuation_control_status_betas_0_5.png", plot=protein_attenuation_control_status2, height=6, width=4, path = "./plots/protein_attenuation_tcga_ccle/")
unlink("protein_attenuation_control_status_betas_0_5.png")




# endoplasmic reticulum membrane proteins
ER <- read_tsv("./data/protein_complexes/subunits_of_multisubunit_complexes_ER_INM.txt") %>%
    mutate(gene_name = toupper(gene_name)) %>%
    dplyr::rename(gene = gene_name)

ER_attenuation <- ER %>%
    dplyr::select(complex_name, gene) %>%
    inner_join(cnv_rna_protein[, -c(2:5)], by = "gene") %>%
    left_join(protein_list_control[protein_list_control$signf == "fdr < 0.05", -c(1)], by = "gene") %>%
    arrange(complex_name) %>%
    mutate_at(vars(contains("control_status")), as.factor)
write.table(ER_attenuation, "./files/ER_membrane_proteins_attenuation.txt", sep="\t", quote=F, row.names=F)



# plot boxplot of complex type versus attenuation potential
ER_attenuation_boxplot <- ggplot(data=ER_attenuation) +
    geom_boxplot(mapping=aes(x=complex_name, y=attenuation, fill=complex_name), outlier.shape = NA) +
    geom_jitter(mapping=aes(x=complex_name, y=attenuation, colour=class, shape=control_status), width=0.1) +
    theme_classic() +
    scale_colour_manual(values = c(rev(brewer.pal(n = 6, "Greys")[4:6]), "red")) +
    scale_shape_manual(values = c(18, 17, 15), labels=c("both", "controlled", "controlling", "unknown"), name="Controlling status", na.value = 16) +
    theme(axis.title.y=element_text(colour="black", size=12),
          axis.title.x=element_text(colour="black", size=12),
          axis.text.y=element_text(colour="black", size=10),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) +
    labs(x = "Complex", y = "Attenuation potential: [CNV ~ RNA] - [CNV ~ Protein]", title = "", fill="Complex name", colour="Attenuation level")
ggsave(filename="protein_attenuation_ER.png", plot=ER_attenuation_boxplot, path = "./plots/protein_attenuation_tcga_ccle/")
unlink("protein_attenuation_ER.png")








save(list=ls(), file = "./r_workspaces/protein_attenuation_control_status.RData")
