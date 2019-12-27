# Abel Sousa 2018

# PhD project
# Study of Control of Protein Abundance using Multi-omics Data from Cancer Samples




# -- Analysis of Protein-Protein Interfaces




suppressMessages(library(tidyverse))
suppressMessages(library(ggpubr))
options(bitmapType = "cairo")




# import protein associations with CNV and RNA
# filtered by genomic co-localization
protein_associations_cnv_mrna <- read_tsv("./files/signf_protein_associations_cnv_mrna_filtered_close_controlling.txt") %>%
    mutate(pair = paste(controlling, controlled, sep="_")) %>%
    dplyr::select(pair, everything())




# number of protein interactions and control status
protein_list <- read_tsv("./files/protein_list_control_status.txt")




# import protein attenuation state
protein_attenuation <- read_tsv("./files/protein_attenuation.txt")




# map of uniprot/gene symbol/protein length
uniprot2gene <- read_tsv("./data/uniprot/uniprot_gene_plength_19_06_2018.tab") %>%
    dplyr::rename(UNIPROT=Entry, SYMBOL=`Gene names  (primary )`, LENGTH=Length) %>%
    dplyr::select(-Status, -Organism)



# load amino acids hydropathy
hydropathy_scale <- read_tsv("./data/aa_hydrophobicity/aa_hydropathy_scale.txt")




# list of all residues that change accessibility when in/out complex
# proxy for interface size
# bind gene symbol and protein length to interfaces
interfaces <- read_tsv("../int3dInterfaces/results/interfaces.tab") %>%
    #filter(STRUCT == "Structure") %>%
    inner_join(uniprot2gene[, c(1,2)], by=c("PROTEIN"="UNIPROT")) %>% dplyr::rename(PROTEIN_GENE = SYMBOL) %>%
    inner_join(uniprot2gene[, c(1,2)], by=c("PARTNER"="UNIPROT")) %>% dplyr::rename(PARTNER_GENE = SYMBOL) %>%
    inner_join(uniprot2gene[, c(1,3)], by=c("PROTEIN"="UNIPROT")) %>% dplyr::rename(PROTEIN_LENGTH = LENGTH) %>%
    inner_join(uniprot2gene[, c(1,3)], by=c("PARTNER"="UNIPROT")) %>% dplyr::rename(PARTNER_LENGTH = LENGTH) %>%
    left_join(hydropathy_scale[, c(2,3)], by=c("AA"="symbol")) %>%
    mutate(GENE_PAIR = paste(PROTEIN_GENE, PARTNER_GENE, sep="_")) %>%
    dplyr::select(GENE_PAIR, PROTEIN, PROTEIN_GENE, PROTEIN_LENGTH, PARTNER, PARTNER_GENE, PARTNER_LENGTH, everything()) %>%
    #remove cases where one gene (symbol) maps to different proteins (uniprot IDs)
    group_by(PROTEIN_GENE) %>%
    filter(length(unique(PROTEIN)) == 1) %>%
    ungroup() %>%
    group_by(PARTNER_GENE) %>%
    filter(length(unique(PARTNER)) == 1) %>%
    ungroup() %>%
    #remove cases where PROTEIN and PARTNER are the same
    filter(PROTEIN_GENE != PARTNER_GENE) %>%
    group_by(INTERFACE, PROTEIN, PARTNER) %>%
    #remove cases where chain length is bigger than uniprot protein length
    filter(unique(PROTEIN_CHAIN_LENGTH) <= unique(PROTEIN_LENGTH) & unique(PARTNER_CHAIN_LENGTH) <= unique(PARTNER_LENGTH)) %>%
    #keep only the pairs where the PROTEIN and PARTNER have at least 100 amino acids in the stucture
    filter(unique(PROTEIN_CHAIN_LENGTH) >= 100 & unique(PARTNER_CHAIN_LENGTH) >= 100) %>%
    #keep only the pairs where the PROTEIN and PARTNER have different chain lengths
    filter(unique(PROTEIN_CHAIN_LENGTH) != unique(PARTNER_CHAIN_LENGTH)) %>%
    ungroup()



# get number/percentages of unique residues on interfaces by PROTEIN
interfaces_size_info <- interfaces %>%
    group_by(PROTEIN_GENE) %>%
    summarise(PROTEIN = unique(PROTEIN), PROTEIN_LENGTH = unique(PROTEIN_LENGTH), UNIQUE_POS = length(unique(POS_UNIPROT)), MEDIAN_HYDROPATHY = median(unique(hydropathy_index)), NUMBER_HYDROPHOBIC_RESIDUES = sum(unique(hydropathy_index) > 0)) %>%
    mutate(AA_PERC_INTERFACE = UNIQUE_POS/PROTEIN_LENGTH, HYDROPHOBIC_RESIDUES_PERC = NUMBER_HYDROPHOBIC_RESIDUES/UNIQUE_POS, HYDROPHOBIC_RESIDUES_PERC_PROTEIN = NUMBER_HYDROPHOBIC_RESIDUES/PROTEIN_LENGTH)











# correlate interface size with effect size (CNV coefficient) and FDR of protein associations
# by protein pair
# interface size measured on the controlling protein

interface_size_controlling <- interfaces %>%
    group_by(PROTEIN_GENE, PARTNER_GENE) %>%
    summarise(INTERFACE_RESIDUES = length(POS_UNIPROT), INTERFACE_RESIDUES_PERC = (length(POS_UNIPROT)/unique(PROTEIN_LENGTH))*100) %>%
    ungroup() %>%
    inner_join(dplyr::select(protein_associations_cnv_mrna, controlling, controlled, beta_cnv, fdr_cnv), by=c("PROTEIN_GENE" = "controlling", "PARTNER_GENE" = "controlled")) %>%
    mutate(fdr_cnv = -log10(fdr_cnv)) %>%
    mutate(beta_cnv = abs(beta_cnv)) %>%
    gather(key="Var", value="value", -c("PROTEIN_GENE", "PARTNER_GENE", "INTERFACE_RESIDUES", "INTERFACE_RESIDUES_PERC")) %>%
    arrange(PROTEIN_GENE)

interface_size_controlling_cor <- interface_size_controlling %>%
    arrange(Var) %>%
    group_by(Var) %>%
    do(broom::tidy(cor.test(.$INTERFACE_RESIDUES, .$value, method = "pearson"))) %>%
    ungroup()

interface_size_controlling_cor2 <- interface_size_controlling %>%
    arrange(Var) %>%
    group_by(Var) %>%
    do(broom::tidy(cor.test(.$INTERFACE_RESIDUES_PERC, .$value, method = "pearson"))) %>%
    ungroup()



interface_size_controlling_plot1 <- ggplot(data = interface_size_controlling, mapping = aes(x = value, y = INTERFACE_RESIDUES)) +
    geom_point() +
    geom_smooth(method=lm , color="black", fill = "grey", se=T) +
    theme_classic() +
    stat_cor() +
    #facet_grid(. ~ Var, scales = "free", labeller=labeller(Var = c(beta_cnv = paste("|beta| CNV", paste("r =", round(interface_size_controlling_cor$estimate[1],2), "p.value =", round(interface_size_controlling_cor$p.value[1], 5)), sep="\n"), fdr_cnv = paste("FDR (-log10)", paste("r =", round(interface_size_controlling_cor$estimate[2], 2), "p.value =", round(interface_size_controlling_cor$p.value[2], 5)), sep="\n")))) +
    facet_wrap( ~ Var, scales = "free", labeller=labeller(Var = c(beta_cnv = "|beta|", fdr_cnv = "FDR (-log10)")) ) +
    theme(axis.title=element_text(colour="black", size=15),
          axis.text=element_text(colour="black", size=12),
          strip.text.x = element_text(size=12),
          strip.background = element_blank(),
          legend.position="none") +
    labs(x = "Protein association statistics", y = "Protein interface size\n Nº of residues in the controlling protein")
ggsave(filename="interface_size_controlling_effect_size_fdr_correlation.png", plot = interface_size_controlling_plot1, path = "./plots/interfaces_protein_associations_cnv/", width=5, height=5)
unlink("interface_size_controlling_effect_size_fdr_correlation.png")



interface_size_controlling_plot2 <- ggplot(data = interface_size_controlling, mapping = aes(x = value, y = INTERFACE_RESIDUES_PERC)) +
    geom_point() +
    geom_smooth(method=lm , color="black", fill = "grey", se=T) +
    theme_classic() +
    stat_cor() +
    #facet_grid(. ~ Var, scales = "free", labeller=labeller(Var = c(beta_cnv = paste("|beta| CNV", paste("r =", round(interface_size_controlling_cor2$estimate[1],2), "p.value =", round(interface_size_controlling_cor2$p.value[1], 5)), sep="\n"), fdr_cnv = paste("FDR (-log10)", paste("r =", round(interface_size_controlling_cor2$estimate[2], 2), "p.value =", round(interface_size_controlling_cor2$p.value[2], 5)), sep="\n")))) +
    facet_wrap( ~ Var, scales = "free", labeller=labeller(Var = c(beta_cnv = "|beta|", fdr_cnv = "FDR (-log10)")) ) +
    theme(axis.title=element_text(colour="black", size=15),
          axis.text=element_text(colour="black", size=12),
          strip.text.x = element_text(size=12),
          strip.background = element_blank(),
          legend.position="none") +
    labs(x = "Protein association statistics", y = "Protein interface size\n % of residues in the controlling protein")
ggsave(filename="interface_size_perc_controlling_effect_size_fdr_correlation.png", plot=interface_size_controlling_plot2, path = "./plots/interfaces_protein_associations_cnv/", width=5, height=5)
unlink("interface_size_perc_controlling_effect_size_fdr_correlation.png")






# correlate interface size with effect size (CNV coefficient) and FDR of protein associations
# by protein pair
# interface size measured on the controlled protein

interface_size_controlled <- interfaces %>%
    group_by(PROTEIN_GENE, PARTNER_GENE) %>%
    summarise(INTERFACE_RESIDUES = length(POS_UNIPROT), INTERFACE_RESIDUES_PERC = (length(POS_UNIPROT)/unique(PROTEIN_LENGTH))*100) %>%
    ungroup() %>%
    inner_join(dplyr::select(protein_associations_cnv_mrna, controlling, controlled, beta_cnv, fdr_cnv), by=c("PROTEIN_GENE" = "controlled", "PARTNER_GENE" = "controlling")) %>%
    mutate(fdr_cnv = -log10(fdr_cnv)) %>%
    mutate(beta_cnv = abs(beta_cnv)) %>%
    gather(key="Var", value="value", -c("PROTEIN_GENE", "PARTNER_GENE", "INTERFACE_RESIDUES", "INTERFACE_RESIDUES_PERC")) %>%
    arrange(PROTEIN_GENE)

interface_size_controlled_cor <- interface_size_controlled %>%
    arrange(Var) %>%
    group_by(Var) %>%
    do(broom::tidy(cor.test(.$INTERFACE_RESIDUES, .$value, method = "pearson"))) %>%
    ungroup()

interface_size_controlled_cor2 <- interface_size_controlled %>%
    arrange(Var) %>%
    group_by(Var) %>%
    do(broom::tidy(cor.test(.$INTERFACE_RESIDUES_PERC, .$value, method = "pearson"))) %>%
    ungroup()



interface_size_controlled_plot1 <- ggplot(data = interface_size_controlled, mapping = aes(x = value, y = INTERFACE_RESIDUES)) +
    geom_point() +
    geom_point(interface_size_controlled %>% filter(PROTEIN_GENE == "CSNK2B" & PARTNER_GENE == "CSNK2A1"), mapping=aes(x=value, y=INTERFACE_RESIDUES), colour="red") +
    geom_point(interface_size_controlled %>% filter(PROTEIN_GENE == "TRMT61A" & PARTNER_GENE == "TRMT6"), mapping=aes(x=value, y=INTERFACE_RESIDUES), colour="red") +
    geom_smooth(method=lm , color="black", fill = "grey", se=T) +
    theme_classic() +
    stat_cor() +
    #facet_grid(. ~ Var, scales = "free", labeller=labeller(Var = c(beta_cnv = paste("|beta| CNV", paste("r =", round(interface_size_controlled_cor$estimate[1],2), "p.value =", round(interface_size_controlled_cor$p.value[1], 5)), sep="\n"), fdr_cnv = paste("FDR (-log10)", paste("r =", round(interface_size_controlled_cor$estimate[2], 2), "p.value =", round(interface_size_controlled_cor$p.value[2], 5)), sep="\n")))) +
    facet_wrap( ~ Var, scales = "free", labeller=labeller(Var = c(beta_cnv = "|beta|", fdr_cnv = "FDR (-log10)")) ) +
    theme(axis.title=element_text(colour="black", size=15),
          axis.text=element_text(colour="black", size=12),
          strip.text.x = element_text(size=12),
          strip.background = element_blank(),
          legend.position="none") +
    scale_y_continuous(limits = c(0, 225)) +
    labs(x = "Protein association statistics", y = "Protein interface size (nº of residues)")
ggsave(filename="interface_size_controlled_effect_size_fdr_correlation.png", plot = interface_size_controlled_plot1, path = "./plots/interfaces_protein_associations_cnv/", width=5, height=5)
ggsave(filename="interface_size_controlled_effect_size_fdr_correlation.pdf", plot = interface_size_controlled_plot1, path = "./plots/interfaces_protein_associations_cnv/", width=5, height=5)
unlink("interface_size_controlled_effect_size_fdr_correlation.png")
unlink("interface_size_controlled_effect_size_fdr_correlation.pdf")



interface_size_controlled_plot2 <- ggplot(data = interface_size_controlled, mapping = aes(x = value, y = INTERFACE_RESIDUES_PERC)) +
    geom_point() +
    geom_smooth(method=lm , color="black", fill = "grey", se=T) +
    theme_classic() +
    stat_cor() +
    #facet_grid(. ~ Var, scales = "free", labeller=labeller(Var = c(beta_cnv = paste("|beta| CNV", paste("r =", round(interface_size_controlled_cor2$estimate[1],2), "p.value =", round(interface_size_controlled_cor2$p.value[1], 5)), sep="\n"), fdr_cnv = paste("FDR (-log10)", paste("r =", round(interface_size_controlled_cor2$estimate[2], 2), "p.value =", round(interface_size_controlled_cor2$p.value[2], 5)), sep="\n")))) +
    facet_wrap( ~ Var, scales = "free", labeller=labeller(Var = c(beta_cnv = "|beta|", fdr_cnv = "FDR (-log10)")) ) +
    theme(axis.title=element_text(colour="black", size=15),
          axis.text=element_text(colour="black", size=12),
          strip.text.x = element_text(size=12),
          strip.background = element_blank(),
          legend.position="none") +
    labs(x = "Protein association statistics", y = "Protein interface size\n% of residues in the controlled protein")
ggsave(filename="interface_size_perc_controlled_effect_size_fdr_correlation.png", plot=interface_size_controlled_plot2, path = "./plots/interfaces_protein_associations_cnv/", width=5, height=5)
unlink("interface_size_perc_controlled_effect_size_fdr_correlation.png")















# correlate number of hydrophobic residues on the interface with effect size (CNV coefficient) and FDR of protein associations
# by protein pair
# number of hydrophobic residues on the interface calculated on the controlling protein
interface_hydrophobicity_controlling <- interfaces %>%
    group_by(PROTEIN_GENE, PARTNER_GENE) %>%
    summarise(NUMBER_HYDROPHOBIC_RESIDUES = sum(hydropathy_index > 0), PERC_HYDROPHOBIC_RESIDUES_INTERFACE = sum(hydropathy_index > 0)/length(POS_UNIPROT), PERC_HYDROPHOBIC_RESIDUES_PROTEIN = sum(hydropathy_index > 0)/unique(PROTEIN_LENGTH)) %>%
    ungroup() %>%
    inner_join(dplyr::select(protein_associations_cnv_mrna, controlling, controlled, beta_cnv, fdr_cnv), by=c("PROTEIN_GENE" = "controlling", "PARTNER_GENE" = "controlled")) %>%
    mutate(fdr_cnv = -log10(fdr_cnv)) %>%
    mutate(beta_cnv = abs(beta_cnv)) %>%
    gather(key="Var", value="value", -c("PROTEIN_GENE", "PARTNER_GENE", "NUMBER_HYDROPHOBIC_RESIDUES", "PERC_HYDROPHOBIC_RESIDUES_INTERFACE", "PERC_HYDROPHOBIC_RESIDUES_PROTEIN")) %>%
    arrange(PROTEIN_GENE)


interface_hydrophobicity_controlling_plot1 <- ggplot(data = interface_hydrophobicity_controlling, mapping = aes(x = value, y = NUMBER_HYDROPHOBIC_RESIDUES)) +
    geom_point() +
    geom_smooth(method=lm , color="blue", se=T) +
    stat_cor() +
    theme_classic() +
    facet_grid(. ~ Var, scales = "free", labeller=labeller(Var = c(beta_cnv = "|beta| CNV", fdr_cnv = "FDR (-log10)"))) +
    theme(axis.title=element_text(colour="black", size=15),
        axis.text=element_text(colour="black", size=12),
        strip.text.x = element_text(size=15),
        strip.background = element_blank(),
        plot.title=element_text(size=15)) +
    labs(x = "Protein association statistics", y = "Number of hydrophobic residues in the controlling protein interface", title="Correlation between the number of hydrophobic residues in the interface\nand CNV beta / -log10(FDR)")
ggsave(filename="interface_hydrophobicity_controlling_effect_size_fdr_correlation.png", plot = interface_hydrophobicity_controlling_plot1, path = "./plots/interfaces_protein_associations_cnv/")
unlink("interface_hydrophobicity_controlling_effect_size_fdr_correlation.png")


interface_hydrophobicity_controlling_plot2 <- ggplot(data = interface_hydrophobicity_controlling, mapping = aes(x = value, y = PERC_HYDROPHOBIC_RESIDUES_INTERFACE)) +
    geom_point() +
    geom_smooth(method=lm , color="blue", se=T) +
    stat_cor() +
    theme_classic() +
    facet_grid(. ~ Var, scales = "free", labeller=labeller(Var = c(beta_cnv = "|beta| CNV", fdr_cnv = "FDR (-log10)"))) +
    theme(axis.title=element_text(colour="black", size=15),
        axis.text=element_text(colour="black", size=12),
        strip.text.x = element_text(size=15),
        strip.background = element_blank(),
        plot.title=element_text(size=15)) +
    labs(x = "Protein association statistics", y = "Percentage of hydrophobic residues in the controlling protein interface\nDivided by the interface size", title="Correlation between the % of hydrophobic residues in the interface\nand CNV beta / -log10(FDR)")
ggsave(filename="interface_hydrophobicity_controlling_effect_size_fdr_correlation_2.png", plot = interface_hydrophobicity_controlling_plot2, path = "./plots/interfaces_protein_associations_cnv/")
unlink("interface_hydrophobicity_controlling_effect_size_fdr_correlation_2.png")


interface_hydrophobicity_controlling_plot3 <- ggplot(data = interface_hydrophobicity_controlling, mapping = aes(x = value, y = PERC_HYDROPHOBIC_RESIDUES_PROTEIN)) +
    geom_point() +
    geom_smooth(method=lm , color="blue", se=T) +
    stat_cor() +
    theme_classic() +
    facet_grid(. ~ Var, scales = "free", labeller=labeller(Var = c(beta_cnv = "|beta| CNV", fdr_cnv = "FDR (-log10)"))) +
    theme(axis.title=element_text(colour="black", size=15),
        axis.text=element_text(colour="black", size=12),
        strip.text.x = element_text(size=15),
        strip.background = element_blank(),
        plot.title=element_text(size=15)) +
    labs(x = "Protein association statistics", y = "Percentage of hydrophobic residues in the controlling protein interface\nDivided by the controlling protein size", title="Correlation between the % of hydrophobic residues in the interface\nand CNV beta / -log10(FDR)")
ggsave(filename="interface_hydrophobicity_controlling_effect_size_fdr_correlation_3.png", plot = interface_hydrophobicity_controlling_plot3, path = "./plots/interfaces_protein_associations_cnv/")
unlink("interface_hydrophobicity_controlling_effect_size_fdr_correlation_3.png")













# correlate number of hydrophobic residues on the interface with effect size (CNV coefficient) and FDR of protein associations
# by protein pair
# number of hydrophobic residues on the interface calculated on the controlled protein
interface_hydrophobicity_controlled <- interfaces %>%
    group_by(PROTEIN_GENE, PARTNER_GENE) %>%
    summarise(NUMBER_HYDROPHOBIC_RESIDUES = sum(hydropathy_index > 0), PERC_HYDROPHOBIC_RESIDUES_INTERFACE = sum(hydropathy_index > 0)/length(POS_UNIPROT), PERC_HYDROPHOBIC_RESIDUES_PROTEIN = sum(hydropathy_index > 0)/unique(PROTEIN_LENGTH)) %>%
    ungroup() %>%
    inner_join(dplyr::select(protein_associations_cnv_mrna, controlling, controlled, beta_cnv, fdr_cnv), by=c("PROTEIN_GENE" = "controlled", "PARTNER_GENE" = "controlling")) %>%
    mutate(fdr_cnv = -log10(fdr_cnv)) %>%
    mutate(beta_cnv = abs(beta_cnv)) %>%
    gather(key="Var", value="value", -c("PROTEIN_GENE", "PARTNER_GENE", "NUMBER_HYDROPHOBIC_RESIDUES", "PERC_HYDROPHOBIC_RESIDUES_INTERFACE", "PERC_HYDROPHOBIC_RESIDUES_PROTEIN")) %>%
    arrange(PROTEIN_GENE)


interface_hydrophobicity_controlled_plot1 <- ggplot(data = interface_hydrophobicity_controlled, mapping = aes(x = value, y = NUMBER_HYDROPHOBIC_RESIDUES)) +
    geom_point() +
    geom_smooth(method=lm , color="blue", se=T) +
    stat_cor() +
    theme_classic() +
    facet_grid(. ~ Var, scales = "free", labeller=labeller(Var = c(beta_cnv = "|beta| CNV", fdr_cnv = "FDR (-log10)"))) +
    theme(axis.title=element_text(colour="black", size=15),
        axis.text=element_text(colour="black", size=12),
        strip.text.x = element_text(size=15),
        strip.background = element_blank(),
        plot.title=element_text(size=15)) +
    labs(x = "Protein association statistics", y = "Number of hydrophobic residues in the controlled protein interface", title="Correlation between the number of hydrophobic residues in the interface\nand CNV beta / -log10(FDR)")
ggsave(filename="interface_hydrophobicity_controlled_effect_size_fdr_correlation.png", plot = interface_hydrophobicity_controlled_plot1, path = "./plots/interfaces_protein_associations_cnv/")
unlink("interface_hydrophobicity_controlled_effect_size_fdr_correlation.png")


interface_hydrophobicity_controlled_plot2 <- ggplot(data = interface_hydrophobicity_controlled, mapping = aes(x = value, y = PERC_HYDROPHOBIC_RESIDUES_INTERFACE)) +
    geom_point() +
    geom_smooth(method=lm , color="blue", se=T) +
    stat_cor() +
    theme_classic() +
    facet_grid(. ~ Var, scales = "free", labeller=labeller(Var = c(beta_cnv = "|beta| CNV", fdr_cnv = "FDR (-log10)"))) +
    theme(axis.title=element_text(colour="black", size=15),
        axis.text=element_text(colour="black", size=12),
        strip.text.x = element_text(size=15),
        strip.background = element_blank(),
        plot.title=element_text(size=15)) +
    labs(x = "Protein association statistics", y = "Percentage of hydrophobic residues in the controlled protein interface\nDivided by the interface size", title="Correlation between the % of hydrophobic residues in the interface\nand CNV beta / -log10(FDR)")
ggsave(filename="interface_hydrophobicity_controlled_effect_size_fdr_correlation_2.png", plot = interface_hydrophobicity_controlled_plot2, path = "./plots/interfaces_protein_associations_cnv/")
unlink("interface_hydrophobicity_controlled_effect_size_fdr_correlation_2.png")


interface_hydrophobicity_controlled_plot3 <- ggplot(data = interface_hydrophobicity_controlled, mapping = aes(x = value, y = PERC_HYDROPHOBIC_RESIDUES_PROTEIN)) +
    geom_point() +
    geom_smooth(method=lm , color="blue", se=T) +
    stat_cor() +
    theme_classic() +
    facet_grid(. ~ Var, scales = "free", labeller=labeller(Var = c(beta_cnv = "|beta| CNV", fdr_cnv = "FDR (-log10)"))) +
    theme(axis.title=element_text(colour="black", size=15),
        axis.text=element_text(colour="black", size=12),
        strip.text.x = element_text(size=15),
        strip.background = element_blank(),
        plot.title=element_text(size=15)) +
    labs(x = "Protein association statistics", y = "Percentage of hydrophobic residues in the controlled protein interface\nDivided by the controlled protein size", title="Correlation between the % of hydrophobic residues in the interface\nand CNV beta / -log10(FDR)")
ggsave(filename="interface_hydrophobicity_controlled_effect_size_fdr_correlation_3.png", plot = interface_hydrophobicity_controlled_plot3, path = "./plots/interfaces_protein_associations_cnv/")
unlink("interface_hydrophobicity_controlled_effect_size_fdr_correlation_3.png")









# select proteins involved in the significative associations
# calculate medians of FDR and effect size
# select the proteins with information about interface size
# correlate interface size with median of effect size (CNV coefficient) and FDR by control status

# protein_list_a <- protein_list %>%
#     inner_join(signf_protein_associations_cnv[, c("controlling", "beta_cnv", "fdr")], by=c("gene" = "controlling"))

# protein_list_b <- protein_list %>%
#     inner_join(signf_protein_associations_cnv[, c("controlled", "beta_cnv", "fdr")], by=c("gene" = "controlled"))

# protein_list_c <- bind_rows(protein_list_a, protein_list_b) %>%
#     group_by(gene) %>%
#     summarise(median_CNV_beta = median(abs(beta_cnv)), median_log10_FDR = median(-log10(fdr)), control_status = unique(control_status)) %>%
#     inner_join(interfaces_size_info, by=c("gene" = "PROTEIN_GENE")) %>%
#     dplyr::select(-PROTEIN_LENGTH, -PROTEIN) %>%
#     gather(key = "var", value = "value", -c(gene, control_status, UNIQUE_POS, AA_PERC_INTERFACE)) %>%
#     arrange(gene)


# interface_size_plot3 <- ggplot(data = protein_list_c, mapping = aes(x = value, y = AA_PERC_INTERFACE)) +
#     geom_point() +
#     geom_smooth(method=lm, color="red", se=T) +
#     facet_grid(control_status ~ var, scales = "free", labeller=labeller(var = c(median_CNV_beta = "median value of |Beta| CNV", median_log10_FDR = "median value of -log10(FDR)"))) +
#     theme(axis.title.y=element_text(colour="black", size=15),
#           axis.title.x=element_text(colour="black", size=15),
#           axis.text.y=element_text(colour="black", size=12),
#           axis.text.x=element_text(colour="black", size=12)) +
#     theme(strip.text = element_text(size=15)) +
#     labs(x = "", y = "Percentage of Unique Residues on Protein Interface")
# ggsave(filename="interface_size_effect_size_fdr_correlation_normalized_plen_all_proteins.png", plot=interface_size_plot3, path = "./plots/interfaces_protein_associations_cnv/")
# unlink("interface_size_effect_size_fdr_correlation_normalized_plen_all_proteins.png")



# interface_size_plot4 <- ggplot(data = protein_list_c, mapping = aes(x = value, y = AA_PERC_INTERFACE)) +
#     geom_point() +
#     geom_smooth(method=lm, color="red", se=T) +
#     facet_grid(. ~ var, scales = "free", labeller=labeller(var = c(median_CNV_beta = "median value of |Beta| CNV", median_log10_FDR = "median value of -log10(FDR)"))) +
#     theme(axis.title.y=element_text(colour="black", size=12),
#           axis.title.x=element_text(colour="black", size=12),
#           axis.text.y=element_text(colour="black", size=10),
#           axis.text.x=element_text(colour="black", size=10)) +
#     theme(strip.text.x = element_text(size=10)) +
#     labs(x = "", y = "Percentage of Unique Residues on Protein Interface")
# ggsave(filename="interface_size_effect_size_fdr_correlation_normalized_plen_all_proteins2.png", plot=interface_size_plot4, path = "./plots/interfaces_protein_associations_cnv/")
# unlink("interface_size_effect_size_fdr_correlation_normalized_plen_all_proteins2.png")


# interface_size_plot5 <- ggplot(data = protein_list_c, mapping = aes(x = value, y = UNIQUE_POS)) +
#     geom_point() +
#     geom_smooth(method=lm, color="red", se=T) +
#     facet_grid(control_status ~ var, scales = "free", labeller=labeller(var = c(median_CNV_beta = "median value of |Beta| CNV", median_log10_FDR = "median value of -log10(FDR)"))) +
#     theme(axis.title.y=element_text(colour="black", size=15),
#           axis.title.x=element_text(colour="black", size=15),
#           axis.text.y=element_text(colour="black", size=12),
#           axis.text.x=element_text(colour="black", size=12)) +
#     theme(strip.text = element_text(size=15)) +
#     labs(x = "", y = "Number of Unique Residues on Protein Interface")
# ggsave(filename="interface_size_effect_size_fdr_correlation_all_proteins.png", plot=interface_size_plot5, path = "./plots/interfaces_protein_associations_cnv/")
# unlink("interface_size_effect_size_fdr_correlation_all_proteins.png")



# interface_size_plot6 <- ggplot(data = protein_list_c, mapping = aes(x = value, y = UNIQUE_POS)) +
#     geom_point() +
#     geom_smooth(method=lm, color="red", se=T) +
#     facet_grid(. ~ var, scales = "free", labeller=labeller(var = c(median_CNV_beta = "median value of |Beta| CNV", median_log10_FDR = "median value of -log10(FDR)"))) +
#     theme(axis.title.y=element_text(colour="black", size=12),
#           axis.title.x=element_text(colour="black", size=12),
#           axis.text.y=element_text(colour="black", size=10),
#           axis.text.x=element_text(colour="black", size=10)) +
#     theme(strip.text.x = element_text(size=10)) +
#     labs(x = "", y = "Number of Unique Residues on Protein Interface")
# ggsave(filename="interface_size_effect_size_fdr_correlation_all_proteins2.png", plot=interface_size_plot6, path = "./plots/interfaces_protein_associations_cnv/")
# unlink("interface_size_effect_size_fdr_correlation_all_proteins2.png")









# comparison between controlled/controlling/both proteins
# percentage of unique interface residues
# percentage of hydrophobic residues on the interface
controlStatusInt <- protein_list %>%
    inner_join(interfaces_size_info, by=c("gene" = "PROTEIN_GENE"))

controlStatusInt_plot1 <- ggplot(data=controlStatusInt, aes(x=control_status, y=AA_PERC_INTERFACE, fill=control_status)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.2) +
    stat_compare_means(comparisons = list( c("both", "controlled"), c("controlled", "controlling"), c("both", "controlling") )) +
    theme_classic() +
    theme(axis.title.y=element_text(colour="black", size=18),
          axis.title.x=element_text(colour="black", size=18),
          axis.text.y=element_text(colour="black", size=15),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          legend.text=element_text(colour="black", size=15),
          plot.title=element_text(colour="black", size=18)) +
    labs(x = "Control status", y = "Percentage of residues annotated\non protein interfaces", title = "Residues annotated on protein interfaces\nand control status", fill="")
ggsave(filename="interfaceSize_distribution_controlStatus.png", plot=controlStatusInt_plot1, path = "./plots/interfaces_protein_associations_cnv/")
unlink("interfaceSize_distribution_controlStatus.png")

controlStatusInt_plot2 <- ggplot(data=controlStatusInt, aes(x=control_status, y=UNIQUE_POS, fill=control_status)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.2) +
    stat_compare_means(comparisons = list( c("both", "controlled"), c("controlled", "controlling"), c("both", "controlling") )) +
    theme_classic() +
    theme(axis.title.y=element_text(colour="black", size=18),
          axis.title.x=element_text(colour="black", size=18),
          axis.text.y=element_text(colour="black", size=15),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          legend.text=element_text(colour="black", size=15),
          plot.title=element_text(colour="black", size=18)) +
    labs(x = "Control status", y = "Number of residues annotated\non protein interfaces", title = "Residues annotated on protein interfaces\nand control status", fill="")
ggsave(filename="interfaceSize_distribution_controlStatus_2.png", plot=controlStatusInt_plot2, path = "./plots/interfaces_protein_associations_cnv/")
unlink("interfaceSize_distribution_controlStatus_2.png")

controlStatusInt_plot3 <- ggplot(data=controlStatusInt, aes(x=control_status, y=HYDROPHOBIC_RESIDUES_PERC, fill=control_status)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.2) +
    stat_compare_means(comparisons = list( c("both", "controlled"), c("controlled", "controlling"), c("both", "controlling") )) +
    theme_classic() +
    theme(axis.title.y=element_text(colour="black", size=18),
          axis.title.x=element_text(colour="black", size=18),
          axis.text.y=element_text(colour="black", size=15),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          legend.text=element_text(colour="black", size=15),
          plot.title=element_text(colour="black", size=18)) +
    labs(x = "Control status", y = "Percentage of hydrophobic residues on protein interfaces\nDivided by number of residues on protein interfaces", title = "Hydrophobic residues annotated on protein\ninterfaces and control status", fill="")
ggsave(filename="hydrophobicResidues_distribution_controlStatus.png", plot=controlStatusInt_plot3, path = "./plots/interfaces_protein_associations_cnv/")
unlink("hydrophobicResidues_distribution_controlStatus.png")

controlStatusInt_plot4 <- ggplot(data=controlStatusInt, aes(x=control_status, y=NUMBER_HYDROPHOBIC_RESIDUES, fill=control_status)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.2) +
    stat_compare_means(comparisons = list( c("both", "controlled"), c("controlled", "controlling"), c("both", "controlling") )) +
    theme_classic() +
    theme(axis.title.y=element_text(colour="black", size=18),
          axis.title.x=element_text(colour="black", size=18),
          axis.text.y=element_text(colour="black", size=15),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          legend.text=element_text(colour="black", size=15),
          plot.title=element_text(colour="black", size=18)) +
    labs(x = "Control status", y = "Number of hydrophobic residues\non protein interfaces", title = "Hydrophobic residues annotated on protein\ninterfaces and control status", fill="")
ggsave(filename="hydrophobicResidues_distribution_controlStatus_2.png", plot=controlStatusInt_plot4, path = "./plots/interfaces_protein_associations_cnv/")
unlink("hydrophobicResidues_distribution_controlStatus_2.png")

controlStatusInt_plot5 <- ggplot(data=controlStatusInt, aes(x=control_status, y=HYDROPHOBIC_RESIDUES_PERC_PROTEIN, fill=control_status)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.2) +
    stat_compare_means(comparisons = list( c("both", "controlled"), c("controlled", "controlling"), c("both", "controlling") )) +
    theme_classic() +
    theme(axis.title.y=element_text(colour="black", size=18),
          axis.title.x=element_text(colour="black", size=18),
          axis.text.y=element_text(colour="black", size=15),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          legend.text=element_text(colour="black", size=15),
          plot.title=element_text(colour="black", size=18)) +
    labs(x = "Control status", y = "Percentage of hydrophobic residues on protein interfaces\nDivided by protein size", title = "Hydrophobic residues annotated on protein\ninterfaces and control status", fill="")
ggsave(filename="hydrophobicResidues_distribution_controlStatus_3.png", plot=controlStatusInt_plot5, path = "./plots/interfaces_protein_associations_cnv/")
unlink("hydrophobicResidues_distribution_controlStatus_3.png")






# analysis of interface size difference across protein pairs
ppairsInterfaceDiff <- protein_associations_cnv_mrna[,c("pair", "controlling", "controlled")] %>%
    inner_join(interfaces_size_info[, c("PROTEIN_GENE", "PROTEIN", "PROTEIN_LENGTH", "UNIQUE_POS", "AA_PERC_INTERFACE")], by=c("controlling" = "PROTEIN_GENE")) %>%
    dplyr::rename(controlling_p = PROTEIN, controlling_plength = PROTEIN_LENGTH, controlling_unique_res = UNIQUE_POS, controlling_perc_resInt = AA_PERC_INTERFACE) %>%
    inner_join(interfaces_size_info[, c("PROTEIN_GENE", "PROTEIN", "PROTEIN_LENGTH", "UNIQUE_POS", "AA_PERC_INTERFACE")], by=c("controlled" = "PROTEIN_GENE")) %>%
    dplyr::rename(controlled_p = PROTEIN, controlled_plength = PROTEIN_LENGTH, controlled_unique_res = UNIQUE_POS, controlled_perc_resInt = AA_PERC_INTERFACE) %>%
    mutate(diff_unique_res = controlling_unique_res - controlled_unique_res, diff_perc_resInt = controlling_perc_resInt-controlled_perc_resInt)


ppairsInterfaceDiff_plot1 <- ggplot(data=ppairsInterfaceDiff) +
    geom_boxplot(mapping = aes(x=factor(0), y=diff_perc_resInt), width=0.3, fill="lightblue") +
    theme_classic() +
    theme(axis.title.y=element_text(colour="black", size=15),
          axis.title.x=element_text(colour="black", size=15),
          axis.text.y=element_text(colour="black", size=13),
          axis.text.x=element_blank(),
          axis.ticks.x = element_blank()) +
    labs(x = "", y = "Difference Between the Percentage of Residues Involved in Interfaces\nby Protein-Protein Pair", title = "")
ggsave(filename="ppairsInterfaceDiff_plot1.png", plot=ppairsInterfaceDiff_plot1, path = "./plots/interfaces_protein_associations_cnv/")
unlink("ppairsInterfaceDiff_plot1.png")


ppairsInterfaceDiff_plot2 <- ggplot(data=ppairsInterfaceDiff) +
    geom_point(mapping = aes(x=controlling_perc_resInt, y=controlled_perc_resInt)) +
    theme_classic() +
    theme(axis.title.y=element_text(colour="black", size=15),
          axis.title.x=element_text(colour="black", size=15),
          axis.text.y=element_text(colour="black", size=13),
          axis.text.x=element_text(colour="black", size=13)) +
    labs(x = "Controlling proteins (% of residues involved in protein interfaces)", y = "Controlled (% of residues involved in protein interfaces)", title = "")
ggsave(filename="ppairsInterfaceDiff_plot2.png", plot=ppairsInterfaceDiff_plot2, path = "./plots/interfaces_protein_associations_cnv/")
unlink("ppairsInterfaceDiff_plot2.png")










# comparison between attenuation states
# percentage of unique interface residues
# percentage of hydrophobic residues on the interface
proteinAttenuationInt <- interfaces_size_info %>%
    inner_join(protein_attenuation, by = c("PROTEIN_GENE" = "gene")) %>%
    mutate_at(vars(contains("class")), as.factor) %>%
    mutate_if(is.factor, fct_infreq)
write.table(proteinAttenuationInt, file = "./files/protein_attenuation_interface.txt", sep="\t", quote=F, row.names=F)


proteinAttenuationInt_plot1 <- ggplot(data=proteinAttenuationInt, mapping = aes(x=class, y=AA_PERC_INTERFACE, fill=class)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.1, alpha = 0.1) +
    stat_compare_means( method = "wilcox.test", label.y = 0.6, comparisons = list( c(1, 3)) ) +
    theme_classic() +
    coord_flip() +
    theme(axis.title=element_text(colour="black", size=10),
      axis.text.x=element_text(colour="black", size=8),
      axis.text.y=element_blank(),
      axis.ticks.y=element_blank(),
      axis.line.y = element_blank(),
      legend.text=element_blank(),
      legend.title=element_blank(),
      plot.title=element_blank(),
      legend.position = "none") +
    scale_y_continuous(limits = c(0, 0.6)) +
    scale_fill_brewer(type = "seq", palette = "Blues", labels = c("non-attenuated", "lowly-attenuated", "highly-attenuated"), name = "", guide=F) +
labs(x = "", y = "Percentage of residues in protein interfaces")
ggsave(filename="interfaceSize_distribution_attenuationState.png", plot=proteinAttenuationInt_plot1, path = "./plots/interfaces_protein_associations_cnv/", width=5, height=1.5)
unlink("interfaceSize_distribution_attenuationState.png")



proteinAttenuationInt_plot1 <- ggplot(data=proteinAttenuationInt, mapping = aes(x=class, y=AA_PERC_INTERFACE, fill=class)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.1) +
  stat_compare_means( method = "wilcox.test", label.y = 0.6, comparisons = list( c(1, 3)) ) +
  theme_classic() +
  theme(axis.title=element_text(colour="black", size=15),
    axis.text.y=element_text(colour="black", size=12),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
  plot.title=element_blank()) +
  scale_y_continuous(limits = c(0, 0.6)) +
  scale_fill_brewer(type = "seq", palette = "Blues", labels = c("non-attenuated", "lowly-attenuated", "highly-attenuated"), name = "", guide = F) +
  labs(x = "Attenuation", y = "Percentage of residues at protein interfaces")
ggsave(filename="interfaceSize_distribution_attenuationState.png", plot=proteinAttenuationInt_plot1, path = "./plots/interfaces_protein_associations_cnv/", width=2, height=5)
ggsave(filename="interfaceSize_distribution_attenuationState.pdf", plot=proteinAttenuationInt_plot1, path = "./plots/interfaces_protein_associations_cnv/", width=2, height=5)
unlink("interfaceSize_distribution_attenuationState.png")
unlink("interfaceSize_distribution_attenuationState.pdf")



proteinAttenuationInt_plot2 <- ggplot(data=proteinAttenuationInt, mapping = aes(x=class, y=UNIQUE_POS, fill=class)) +
    geom_boxplot(outlier.shape = NA, show.legend=T) +
    geom_jitter(width = 0.2, alpha = 0.1) +
    stat_compare_means(comparisons = list( c(1, 2), c(2, 3), c(1, 3) )) +
    theme_classic() +
    theme(axis.title.y=element_text(colour="black", size=18),
          axis.title.x=element_text(colour="black", size=18),
          axis.text.y=element_text(colour="black", size=15),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          legend.text=element_text(colour="black", size=15),
          plot.title=element_text(colour="black", size=18)) +
    labs(x = "Attenuation state", y = "Number of residues annotated\non protein interfaces", title = "Residues annotated on protein interfaces\nand attenuation state", fill="")
ggsave(filename="interfaceSize_distribution_attenuationState_2.png", plot=proteinAttenuationInt_plot2, path = "./plots/interfaces_protein_associations_cnv/")
unlink("interfaceSize_distribution_attenuationState_2.png")


proteinAttenuationInt_plot3 <- ggplot(data = proteinAttenuationInt, mapping = aes(x = attenuation, y = AA_PERC_INTERFACE)) +
    geom_point() +
    geom_smooth(method=lm, color="red", se=T) +
    theme(axis.title.y=element_text(colour="black", size=12),
          axis.title.x=element_text(colour="black", size=12),
          axis.text.y=element_text(colour="black", size=10),
          axis.text.x=element_text(colour="black", size=10)) +
    labs(x = "Attenuation potential", y = "Percentage of Unique Residues on Protein Interface")
ggsave(filename="interfaceSize_correlation_attenuationPotential.png", plot=proteinAttenuationInt_plot3, path = "./plots/interfaces_protein_associations_cnv/")
unlink("interfaceSize_correlation_attenuationPotential.png")

proteinAttenuationInt_plot4 <- ggplot(data=proteinAttenuationInt, mapping = aes(x=class, y=NUMBER_HYDROPHOBIC_RESIDUES, fill=class)) +
    geom_boxplot(outlier.shape = NA, show.legend=T) +
    geom_jitter(width = 0.2, alpha = 0.1) +
    stat_compare_means(comparisons = list( c(1, 2), c(2, 3), c(1, 3) )) +
    theme_classic() +
    theme(axis.title.y=element_text(colour="black", size=18),
          axis.title.x=element_text(colour="black", size=18),
          axis.text.y=element_text(colour="black", size=15),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          legend.text=element_text(colour="black", size=15),
          plot.title=element_text(colour="black", size=18)) +
    labs(x = "Attenuation state", y = "Number of hydrophobic residues\non protein interfaces", title = "Hydrophobic residues annotated on protein\ninterfaces and attenuation state", fill="")
ggsave(filename="hydrophobicResidues_distribution_attenuationState.png", plot=proteinAttenuationInt_plot4, path = "./plots/interfaces_protein_associations_cnv/")
unlink("hydrophobicResidues_distribution_attenuationState.png")

proteinAttenuationInt_plot5 <- ggplot(data=proteinAttenuationInt, mapping = aes(x=class, y=HYDROPHOBIC_RESIDUES_PERC, fill=class)) +
    geom_boxplot(outlier.shape = NA, show.legend=T) +
    geom_jitter(width = 0.2, alpha = 0.1) +
    stat_compare_means(comparisons = list( c(1, 2), c(2, 3), c(1, 3) )) +
    theme_classic() +
    theme(axis.title.y=element_text(colour="black", size=18),
          axis.title.x=element_text(colour="black", size=18),
          axis.text.y=element_text(colour="black", size=15),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          legend.text=element_text(colour="black", size=15),
          plot.title=element_text(colour="black", size=18)) +
    labs(x = "Attenuation state", y = "Percentage of hydrophobic residues on protein interfaces\nDivided by number of residues on protein interfaces", title = "Hydrophobic residues annotated on protein\ninterfaces and attenuation state", fill="")
ggsave(filename="hydrophobicResidues_distribution_attenuationState_2.png", plot=proteinAttenuationInt_plot5, path = "./plots/interfaces_protein_associations_cnv/")
unlink("hydrophobicResidues_distribution_attenuationState_2.png")

proteinAttenuationInt_plot6 <- ggplot(data=proteinAttenuationInt, mapping = aes(x=class, y=HYDROPHOBIC_RESIDUES_PERC_PROTEIN, fill=class)) +
    geom_boxplot(outlier.shape = NA, show.legend=T) +
    geom_jitter(width = 0.2, alpha = 0.1) +
    stat_compare_means(comparisons = list( c(1, 2), c(2, 3), c(1, 3) )) +
    theme_classic() +
    theme(axis.title.y=element_text(colour="black", size=18),
          axis.title.x=element_text(colour="black", size=18),
          axis.text.y=element_text(colour="black", size=15),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          legend.text=element_text(colour="black", size=15),
          plot.title=element_text(colour="black", size=18)) +
    labs(x = "Attenuation state", y = "Percentage of hydrophobic residues on protein interfaces\nDivided by protein size", title = "Hydrophobic residues annotated on protein\ninterfaces and attenuation state", fill="")
ggsave(filename="hydrophobicResidues_distribution_attenuationState_3.png", plot=proteinAttenuationInt_plot6, path = "./plots/interfaces_protein_associations_cnv/")
unlink("hydrophobicResidues_distribution_attenuationState_3.png")







# comparison of protein size between controlled/controlling/both proteins
controlStatusProtSize <- protein_list %>%
    inner_join(uniprot2gene, by=c("gene" = "SYMBOL")) %>%
    group_by(gene) %>%
    filter(n() == 1) %>%
    ungroup()


controlStatusProtSize_plot1 <- ggplot(data=controlStatusProtSize, mapping = aes(x=control_status, y=LENGTH, fill=control_status)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.1, alpha = 0.1) +
    stat_compare_means(comparisons = list( c("controlling", "controlled") ), label.y=log10(3000)) +
    theme_classic() +
    theme(axis.title=element_text(colour="black", size=15),
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      axis.text.y=element_text(colour="black", size=12),
      plot.title=element_blank(),
      legend.position = "none") +
    scale_y_continuous(trans="log10", limits = c(NA, 3000), breaks = c(100, 500, 3000)) +
    ggpubr::fill_palette("jco") +
    labs(x = "Control status", y = "Protein length (log10)")
ggsave(filename="controlStatus_proteinSize_distribution.png", plot=controlStatusProtSize_plot1, path = "./plots/interfaces_protein_associations_cnv/", width=2, height=5)
ggsave(filename="controlStatus_proteinSize_distribution.pdf", plot=controlStatusProtSize_plot1, path = "./plots/interfaces_protein_associations_cnv/", width=2, height=5)
unlink("controlStatus_proteinSize_distribution.png")
unlink("controlStatus_proteinSize_distribution.pdf")









# analysis of protein size difference across protein pairs
pairs_size_diff <- protein_associations_cnv_mrna[,c("pair", "controlling", "controlled")] %>%
    inner_join(uniprot2gene, by=c("controlling"="SYMBOL")) %>%
    dplyr::rename(length_controlling=LENGTH) %>%
    dplyr::select(-UNIPROT) %>%
    group_by(pair) %>%
    summarise(controlling=unique(controlling), controlled=unique(controlled), length_controlling=max(length_controlling)) %>%
    ungroup() %>%
    inner_join(uniprot2gene, by=c("controlled"="SYMBOL")) %>%
    dplyr::rename(length_controlled=LENGTH) %>%
    dplyr::select(-UNIPROT) %>%
    group_by(pair) %>%
    summarise(controlling=unique(controlling), controlled=unique(controlled), length_controlling=max(length_controlling), length_controlled=max(length_controlled)) %>%
    ungroup() %>%
    mutate(length_diff = length_controlling-length_controlled, length_diff_3root = kader:::cuberoot(length_controlling-length_controlled), stat = "size_difference")


pairs_size_diff_plot1 <- ggplot(data=pairs_size_diff) +
    geom_boxplot(mapping = aes(x=stat, y=length_diff_3root), width=0.3, fill="lightblue") +
    theme_classic() +
    theme(axis.title.y=element_text(colour="black", size=15),
          axis.title.x=element_text(colour="black", size=15),
          axis.text.y=element_text(colour="black", size=13),
          axis.text.x=element_blank(),
          axis.ticks.x = element_blank()) +
    labs(x = "", y = "Cube Root of the Protein Size Difference by Protein-Protein Pair", title = "")
ggsave(filename="pairs_size_diff_plot1.png", plot=pairs_size_diff_plot1, path = "./plots/interfaces_protein_associations_cnv/")
unlink("pairs_size_diff_plot1.png")





# comparison of protein size between attenuation state
attenuationStatusProtSize <- protein_attenuation[,c(1,6,7)] %>%
    inner_join(uniprot2gene[, c(2,3)], by=c("gene" = "SYMBOL")) %>%
    group_by(gene, attenuation, class) %>%
    filter(n() == 1) %>%
    ungroup() %>%
    mutate_at(vars(contains("class")), as.factor) %>%
    mutate_if(is.factor, fct_infreq)



attenuationStatusProtSize_plot1 <- ggplot(data=attenuationStatusProtSize, mapping = aes(x=class, y=LENGTH, fill=class)) +
    geom_boxplot(outlier.shape=NA) +
    #geom_hline(yintercept = median(attenuationStatusProtSize$LENGTH), linetype = 2) +
    geom_jitter(width = 0.1, alpha = 0.01) +
    stat_compare_means( method = "wilcox.test", comparisons = list( c(1, 3)), label.y=log10(4000) ) +
    #stat_compare_means(label = "p.signif", method = "wilcox.test", ref.group = ".all.") +
    theme_classic() +
    theme(axis.title.y=element_text(colour="black", size=15),
        axis.title.x=element_text(colour="black", size=15),
        axis.text.y=element_text(colour="black", size=12),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        plot.title = element_blank()) +
    scale_y_continuous(trans="log10", limits = c(NA, 4000), breaks = c(100, 500, 4000)) +
    #scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x), labels = scales::trans_format("log10", scales::math_format(10^.x)), limits=c(NA, 4000) ) +
    #annotation_logticks(sides = "l") +
    scale_fill_brewer(type = "seq", palette = "Blues", labels = c("non-attenuated", "lowly-attenuated", "highly-attenuated"), name = "", guide=F) +
    labs(x = "Attenuation", y = "Protein length (log10)")
ggsave(filename="attenuationStatus_proteinSize_distribution.png", plot=attenuationStatusProtSize_plot1, path = "./plots/interfaces_protein_associations_cnv/", width=2, height=5)
ggsave(filename="attenuationStatus_proteinSize_distribution.pdf", plot=attenuationStatusProtSize_plot1, path = "./plots/interfaces_protein_associations_cnv/", width=2, height=5)
unlink("attenuationStatus_proteinSize_distribution.png")
unlink("attenuationStatus_proteinSize_distribution.pdf")



attenuationPotentialProtSize_plot1 <- ggplot(data = attenuationStatusProtSize, mapping = aes(x = attenuation, y = LENGTH)) +
    geom_point() +
    geom_smooth(method=lm, color="red", se=T) +
    theme(axis.title.y=element_text(colour="black", size=12),
          axis.title.x=element_text(colour="black", size=12),
          axis.text.y=element_text(colour="black", size=10),
          axis.text.x=element_text(colour="black", size=10)) +
    labs(x = "Attenuation potential", y = "Protein length (log2)")
ggsave(filename="attenuationPotential_proteinSize_scatterplot.png", plot=attenuationPotentialProtSize_plot1, path = "./plots/interfaces_protein_associations_cnv/")
unlink("attenuationPotential_proteinSize_scatterplot.png")


attenuationPotentialProtSize_plot2 <- ggplot(data = attenuationStatusProtSize, mapping = aes(x = attenuation, y = LENGTH)) +
    geom_point() +
    geom_smooth(method=lm, color="red", se=T) +
    facet_grid(. ~ class, scales = "free") +
    theme(axis.title.y=element_text(colour="black", size=12),
          axis.title.x=element_text(colour="black", size=12),
          axis.text.y=element_text(colour="black", size=10),
          axis.text.x=element_text(colour="black", size=10)) +
    labs(x = "Attenuation potential", y = "Protein length (log2)")
ggsave(filename="attenuationPotential_proteinSize_scatterplot2.png", plot=attenuationPotentialProtSize_plot2, path = "./plots/interfaces_protein_associations_cnv/")
unlink("attenuationPotential_proteinSize_scatterplot2.png")







# check if there is a phosphosite in a protein interface

signf_cnv_rna_phospho <- read_tsv("./files/signf_protein_associations_cnv_rna_phospho.txt") %>%
    dplyr::select(controlling, controlling_phx, controlled) %>%
    separate(col = controlling_phx, into = c("x", "phosphosite"), sep = "_") %>%
    dplyr::select(-x) %>%
    mutate(phosphosite = gsub("[a-z]", "", phosphosite)) %>%
    filter(controlling != "SRSF9", controlled != "SKIV2L2") %>%
    mutate(phosphosite = as.numeric(phosphosite))


phosphosite_interface_controlling <- interfaces %>%
    dplyr::select(PROTEIN_GENE, PARTNER_GENE, POS_UNIPROT) %>%
    inner_join(signf_cnv_rna_phospho, by=c("PROTEIN_GENE" = "controlling", "POS_UNIPROT" = "phosphosite"))




save(list=ls(), file="./r_workspaces/protein_interfaces.RData")
