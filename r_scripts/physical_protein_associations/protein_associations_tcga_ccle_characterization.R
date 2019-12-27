# Abel Sousa 2018

# PhD project
# Study of Control of Protein Abundance using Multi-omics Data from Cancer Samples



# -- Exploration of protein associations found with CNV



suppressMessages(library(tidyverse))
suppressMessages(library(ggrepel))
suppressMessages(library(RColorBrewer))
suppressMessages(library(clusterProfiler))
suppressMessages(library(org.Hs.eg.db))
suppressMessages(library(igraph))
suppressMessages(library(GenomicRanges))

options(bitmapType = "cairo")





# import protein associations with CNV and RNA
# filtered by genomic co-localization
protein_associations_cnv_mrna <- read_tsv("./files/signf_protein_associations_cnv_mrna_filtered_close_controlling.txt") %>%
    mutate(pair = paste(controlling, controlled, sep="_")) %>%
    dplyr::select(pair, everything())



# import protein associations with CNV
protein_associations_cnv <- read_tsv("./files/protein_associations_tcga_ccle_cnv_reg.txt") %>%
    mutate(pair = paste(controlling, controlled, sep="_")) %>%
    dplyr::select(pair, controlling, controlled, beta_cnv, fdr) %>%
    mutate(fdr_label = if_else(fdr > 0.05, "fdr > 0.05", "fdr <= 0.05")) %>%
    mutate(in_RNA_model = if_else(fdr > 0.05 & !(pair %in% protein_associations_cnv_mrna$pair), "1", if_else(fdr < 0.05 & !(pair %in% protein_associations_cnv_mrna$pair), "2", "3")))




# volcano plot
# plot CNV beta by FDR
# color by FDR and mark as red the associations found with the RNA
protein_associations_cnv_plot <- ggplot( data=protein_associations_cnv, mapping=aes(x=beta_cnv, y=-log10(fdr), color=in_RNA_model, fill = fdr_label) ) +
    geom_point(pch = 21) +
    theme_classic() +
    geom_text_repel(
      data = bind_rows(
        protein_associations_cnv %>% filter(in_RNA_model == "3" & beta_cnv > 0) %>% arrange(fdr) %>% head(2),
        protein_associations_cnv %>% filter(in_RNA_model == "3" & beta_cnv > 0) %>% arrange(fdr) %>% filter(grepl("TCP1", controlling, fixed=T)) %>% head(4),
        protein_associations_cnv %>% filter(in_RNA_model == "3" & beta_cnv > 0) %>% arrange(fdr) %>% filter(grepl("COPS3", controlling, fixed=T)) %>% head(2),
        protein_associations_cnv %>% filter(in_RNA_model == "3" & beta_cnv > 0) %>% arrange(fdr) %>% filter(grepl("IDH3A", controlling, fixed=T)) %>% head(2),
        protein_associations_cnv %>% filter(in_RNA_model == "3" & beta_cnv < 0) %>% arrange(fdr) %>% head(5)),
      aes(x=beta_cnv, y=-log10(fdr), label = paste(controlling, controlled, sep="~")),
      size = 3, box.padding = unit(0.3, "lines"), point.padding = unit(0.3, "lines"), colour="black") +
    scale_fill_manual(values=c("#6baed6", "#bdd7e7"), labels = c("FDR <= 0.05", "FDR > 0.05"), guide=F) +
    scale_colour_manual(values=c("#bdd7e7", "#6baed6", "#2171b5"), guide=F) +
    scale_x_continuous(limits = c(-0.6, 0.7), breaks=round(seq(-0.6, 0.7, 0.2), 1)) +
    geom_line(aes(x=0),colour="black", linetype=2, size = 0.3) +
    geom_line(aes(y=-log10(0.05)),colour="black", linetype=2, size = 0.3) +
    theme(axis.title=element_text(colour="black", size=15),
          axis.text=element_text(colour="black", size=12),
          plot.title = element_blank(),
          legend.text=element_text(colour="black", size=12),
          legend.title=element_blank()) +
    labs(x = "CNV Beta", y = "CNV model FDR (-log10)")
ggsave(filename="protein_associations_tcga_ccle_cnv_rna_reg.png", plot=protein_associations_cnv_plot, width=6, height=6, path="./plots/protein_associations_tcga_ccle/")
ggsave(filename="protein_associations_tcga_ccle_cnv_rna_reg.pdf", plot=protein_associations_cnv_plot, width=6, height=6, path="./plots/protein_associations_tcga_ccle/")
unlink("protein_associations_tcga_ccle_cnv_rna_reg.png")
unlink("protein_associations_tcga_ccle_cnv_rna_reg.pdf")





# import complete data from BioGRID
biogrid <- read_tsv("./data/biogrid/BIOGRID-ORGANISM-3/BIOGRID-ORGANISM-Homo_sapiens-3.4.157.tab2.txt") %>%
    filter(`Organism Interactor A` == 9606, `Organism Interactor B` == 9606)



# import BioGRID experiments
# select only experiments defined as stable or transient
biogrid_type <- read_csv("./data/biogrid/physical_experimental_systems.csv") %>%
    filter(Group == "stable" | Group == "transient")



# load protein interactions
all_interactions <- read_tsv("./files/biogrid_corum_ER_pairs.txt") %>%
    mutate(pair = paste(interactor.A, interactor.B, sep="_"))


all_interactions2 <- read_tsv("./files/biogrid_corum_ER_pairs2.txt") %>%
    mutate(pair = paste(interactor.A, interactor.B, sep="_"))




# load genomic coordinates
load("/nfs/research1/beltrao/memon/databases/protein/ensembl/ensembl_v73/ensembl_gtf_v73_gene_anno_grange.Rdata")
ensembl_gtf_v73_gene_anno_grange <- as.tibble(as.data.frame(ensembl_gtf_v73_gene_anno_grange))








# significant protein associations
# add column with Pubmed IDs that support the interaction
protein_associations_cnv_mrna2 <- protein_associations_cnv_mrna %>%
    left_join(all_interactions2[, c("pair", "experiment.type")], by="pair")

protein_associations_cnv_mrna2 <- protein_associations_cnv_mrna2 %>%
    mutate( pmids = apply( as.matrix(dplyr::select(protein_associations_cnv_mrna2, c(controlling, controlled) )), 1, function(x) filter(biogrid, `Official Symbol Interactor A` == x[1] & `Official Symbol Interactor B` == x[2] | `Official Symbol Interactor A` == x[2] & `Official Symbol Interactor B` == x[1]) %>% pull(`Pubmed ID`) %>% unique() %>% paste(collapse="|") ) )
write.table(protein_associations_cnv_mrna2, "./files/signf_protein_associations_cnv_mrna_filtered_close_controlling2.txt", sep="\t", quote=F, row.names=F)








# proportion of experiment types between tested/significant protein pairs
exp_type_proportion_tested <- all_interactions %>%
    filter(pair %in% protein_associations_cnv$pair) %>%
    group_by(experiment.type) %>%
    summarise(freq=n()) %>%
    mutate(perc=(freq/sum(freq))*100) %>%
    mutate(label = "Tested pairs")
exp_type_proportion_signf <- all_interactions %>%
    filter(pair %in% protein_associations_cnv_mrna$pair) %>%
    group_by(experiment.type) %>%
    summarise(freq=n()) %>%
    mutate(perc=(freq/sum(freq))*100) %>%
    mutate(label = "Significative pairs")
exp_type_proportion_ppairs <- bind_rows(exp_type_proportion_tested, exp_type_proportion_signf) %>%
    left_join(dplyr::select(biogrid_type, c(Experiment, Group)), by = c("experiment.type" = "Experiment")) %>%
    replace(., is.na(.), c("corum", "curated_membrane_complexes"))


exp_type_proportion_ppairs_plot <- ggplot(data=exp_type_proportion_ppairs) +
    geom_bar(mapping=aes(x=experiment.type, y=perc, fill=label, color=Group), stat="identity", position="dodge") +
    theme_classic() +
    scale_color_manual(values = c("stable"="black", "transient"="grey", "corum"="blue", "curated_membrane_complexes" = "red")) +
    theme(legend.title=element_blank()) +
    theme(plot.margin = margin(1, 0, 0, 2, "cm")) +
    theme(axis.title.y=element_text(colour="black", size=12),
          axis.text.y=element_text(colour="black", size=10),
          axis.title.x=element_text(colour="black", size=12),
          axis.text.x=element_text(colour="black", size=10, angle = 45, hjust = 1)) +
    labs(x = "Experimental types", y = "Percentage of Protein Pairs Identified")
ggsave("perc_experiment_types_tested_sign_ppairs.png", plot = exp_type_proportion_ppairs_plot, path="./plots/protein_associations_tcga_ccle/")
unlink("perc_experiment_types_tested_sign_ppairs.png")



# proportion of pubmed IDs number
biogrid_pubmed <- biogrid %>%
    filter(`Experimental System Type` == "physical") %>%
    filter(`Experimental System` %in% as.character(biogrid_type$Experiment)) %>%
    dplyr::select(`Official Symbol Interactor A`, `Official Symbol Interactor B`, `Pubmed ID`) %>%
    mutate_all(as.character) %>%
    dplyr::rename(interactor.A=`Official Symbol Interactor A`, interactor.B=`Official Symbol Interactor B`, pubmed_ID=`Pubmed ID`)

biogrid_pubmed <- biogrid_pubmed %>%
    bind_rows(setNames(dplyr::select(biogrid_pubmed, c(interactor.B, interactor.A, pubmed_ID)), c("interactor.A", "interactor.B", "pubmed_ID"))) %>%
    distinct() %>%
    filter(interactor.A != interactor.B)

biogrid_pubmed <- biogrid_pubmed %>%
    mutate(pair = paste(interactor.A, interactor.B, sep="_")) %>%
    dplyr::select(pair, everything())


pubmed_number_tested <- protein_associations_cnv %>%
    dplyr::select(pair) %>%
    left_join(biogrid_pubmed, by="pair") %>%
    filter(!is.na(pubmed_ID)) %>%
    group_by(pair) %>%
    summarise(n_diff_papers=length(unique(pubmed_ID))) %>%
    group_by(n_diff_papers) %>%
    summarise(freq = n())

pubmed_number_tested_10 <- pubmed_number_tested %>%
    filter(n_diff_papers > 10) %>%
    mutate(n_diff_papers = ">10", freq = sum(freq)) %>%
    distinct()

pubmed_number_tested <- pubmed_number_tested %>%
    filter(n_diff_papers <= 10) %>%
    mutate(n_diff_papers = as.character(n_diff_papers)) %>%
    bind_rows(pubmed_number_tested_10) %>%
    mutate(label = "Tested pairs")

pubmed_number_signf <- protein_associations_cnv_mrna %>%
    dplyr::select(pair) %>%
    left_join(biogrid_pubmed, "pair") %>%
    filter(!is.na(pubmed_ID)) %>%
    group_by(pair) %>%
    summarise(n_diff_papers=length(unique(pubmed_ID))) %>%
    group_by(n_diff_papers) %>%
    summarise(freq = n())

pubmed_number_signf_10 <- pubmed_number_signf %>%
    filter(n_diff_papers > 10) %>%
    mutate(n_diff_papers = ">10", freq = sum(freq)) %>%
    distinct()

pubmed_number_signf <- pubmed_number_signf %>%
    filter(n_diff_papers <= 10) %>%
    mutate(n_diff_papers = as.character(n_diff_papers)) %>%
    bind_rows(pubmed_number_signf_10) %>%
    mutate(label = "Significative pairs")

pubmed_number_tested_signf <- bind_rows(pubmed_number_tested, pubmed_number_signf) %>%
    group_by(label) %>%
    mutate(perc = (freq/sum(freq))*100) %>%
    ungroup() %>%
    mutate(n_diff_papers = factor(n_diff_papers, levels=c(1:10, ">10")))



pubmed_number_tested_signf_plot <- ggplot(data=pubmed_number_tested_signf) +
    geom_bar(mapping=aes(x=n_diff_papers, y=perc, fill=label), stat="identity", position="dodge") +
    theme_classic() +
    theme(axis.title.y=element_text(colour="black", size=12),
          axis.title.x=element_text(colour="black", size=12),
          axis.text.y=element_text(colour="black", size=10),
          axis.text.x=element_text(colour="black", size=10)) +
    labs(x = "Number of different articles supporting the interaction", y = "Percentage of protein interactions")
ggsave("perc_ppairs_tested_signf_articles.png", plot = pubmed_number_tested_signf_plot, path="./plots/protein_associations_tcga_ccle/")
unlink("perc_ppairs_tested_signf_articles.png")






# number of times each protein is found and control status
protein_list <- protein_associations_cnv_mrna %>%
    dplyr::select(controlling, controlled) %>%
    dplyr::combine() %>%
    as_tibble() %>%
    dplyr::rename(gene = value) %>%
    group_by(gene) %>%
    dplyr::summarise(interactions = n()) %>%
    arrange(desc(interactions)) %>%
    dplyr::mutate(control_status = if_else(gene %in% setdiff(protein_associations_cnv_mrna$controlling, protein_associations_cnv_mrna$controlled), "controlling", if_else(gene %in% setdiff(protein_associations_cnv_mrna$controlled, protein_associations_cnv_mrna$controlling), "controlled", "both")))
nrow(protein_list)
write.table(protein_list, "./files/protein_list_control_status.txt", sep="\t", quote=F, row.names=F)
#836


control_status_barplot <- ggplot(data=protein_list) +
    geom_bar(mapping=aes(x=control_status, y=..count.., fill=control_status)) +
    theme_classic() +
    coord_flip() +
    theme(axis.title=element_text(colour="black", size=12),
          axis.text.x=element_text(colour="black", size=10),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          legend.text=element_text(colour="black", size=10),
          legend.title=element_text(colour="black", size=12),
          plot.title=element_blank(),
          legend.position = "bottom") +
    ggpubr::fill_palette("jco") +
    labs(x = "", y = "Number of proteins", fill = "Control status")
ggsave("number_controlling_controlled_both_proteins.png", plot = control_status_barplot, height=2, width=5, path="./plots/protein_associations_tcga_ccle/")
unlink("number_controlling_controlled_both_proteins.png")


control_status_barplot <- ggplot(data=protein_list) +
    geom_bar(mapping=aes(x=control_status, y=..count.., fill=control_status)) +
    theme_classic() +
    theme(axis.title=element_text(colour="black", size=15),
          axis.text.y=element_text(colour="black", size=12),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          legend.text=element_text(colour="black", size=12),
          legend.title=element_text(colour="black", size=15),
          plot.title=element_blank(),
          legend.position = "right") +
    ggpubr::fill_palette("jco") +
    labs(x = "Control status", y = "Number of proteins", fill = "")
ggsave("number_controlling_controlled_both_proteins.png", plot = control_status_barplot, height=4, width=4, path="./plots/protein_associations_tcga_ccle/")
ggsave("number_controlling_controlled_both_proteins.pdf", plot = control_status_barplot, height=4, width=4, path="./plots/protein_associations_tcga_ccle/")
unlink("number_controlling_controlled_both_proteins.png")
unlink("number_controlling_controlled_both_proteins.pdf")


control_status_interactions_dist <- ggplot(data=protein_list) +
    geom_boxplot(mapping=aes(x=control_status, y=interactions, fill=control_status)) +
    theme_classic() +
    theme(axis.title.y=element_text(colour="black", size=12),
          axis.title.x=element_text(colour="black", size=12),
          axis.text.y=element_text(colour="black", size=10),
          axis.text.x=element_text(colour="black", size=10)) +
    scale_fill_brewer(palette = "Set2", type="qual") +
    labs(x = "Control status", y = "Number of interactions")
ggsave("distribution_interactions_controlling_controlled_both_proteins.png", plot=control_status_interactions_dist, path="./plots/protein_associations_tcga_ccle/")
unlink("distribution_interactions_controlling_controlled_both_proteins.png")





# characterization of proteins that act as controlling and controlled (both)
both_proteins <- protein_list %>%
    filter(control_status == "both") %>%
    inner_join(protein_associations_cnv_mrna[, c(2,4,5)], by=c("gene" = "controlling")) %>%
    mutate(fdr_cnv = -log10(fdr_cnv)) %>%
    group_by(gene) %>%
    summarise(interactions = unique(interactions), control_status = unique(control_status), controlling_times = n(), beta_cnv_controlling = median(beta_cnv), fdr_controlling = median(fdr_cnv)) %>%
    inner_join(protein_associations_cnv_mrna[, c(3,4,5)], by=c("gene" = "controlled")) %>%
    mutate(fdr_cnv = -log10(fdr_cnv)) %>%
    group_by(gene) %>%
    summarise(interactions = unique(interactions), control_status = unique(control_status), controlling_times = unique(controlling_times), beta_cnv_controlling = unique(beta_cnv_controlling), fdr_controlling = unique(fdr_controlling), controlled_times = n(), beta_cnv_controlled = median(beta_cnv), fdr_controlled = median(fdr_cnv)) %>%
    arrange(desc(interactions)) %>%
    mutate(beta_cnv_diff = beta_cnv_controlling - beta_cnv_controlled) %>%
    mutate(log10_fdr_diff = fdr_controlling - fdr_controlled)


both_proteins_interactions <- both_proteins[, c("interactions", "controlling_times", "controlled_times")] %>%
    gather(key="stat", value="value") %>%
    mutate_if(is.character, as.factor) %>%
    mutate(stat = fct_relevel(stat, "interactions", "controlled_times", "controlling_times"))


both_proteins_interactions_barplot <- ggplot(data=both_proteins_interactions) +
    geom_bar(mapping=aes(x=stat, y=value, fill=stat), stat="identity") +
    theme(axis.title.y=element_text(colour="black", size=15),
          axis.title.x=element_text(colour="black", size=15),
          axis.text.y=element_text(colour="black", size=13),
          axis.text.x=element_text(colour="black", size=13),
          legend.position = "none") +
    #scale_fill_brewer(palette = "Accent", type="qual") +
    scale_x_discrete(labels = c("Total interactions", "Interactions as controlled", "Interactions as controlling"), name="") +
    scale_y_continuous(limits = c(0,700), name="Number of interactions") +
    coord_flip()
ggsave("number_interactions_bothProteins.png", plot = both_proteins_interactions_barplot, path="./plots/protein_associations_tcga_ccle/")
unlink("number_interactions_bothProteins.png")


both_proteins_interactions_beta_diff <- ggplot(data=both_proteins) +
    geom_boxplot(mapping=aes(x=factor(0), y=beta_cnv_diff), width=0.5, fill = "lightblue") +
    theme(axis.title=element_text(colour="black", size=15),
          axis.text.y=element_text(colour="black", size=13),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) +
    scale_x_discrete(labels = c("Both proteins")) +
    labs(x = "Both proteins", y = "Difference between median beta values")
ggsave("beta_cnv_diff_bothProteins.png", plot=both_proteins_interactions_beta_diff, path="./plots/protein_associations_tcga_ccle/")
unlink("beta_cnv_diff_bothProteins.png")


both_proteins_interactions_fdr_diff <- ggplot(data=both_proteins) +
    geom_boxplot(mapping=aes(x=factor(0), y=log10_fdr_diff), width=0.5, fill = "lightblue") +
    theme(axis.title=element_text(colour="black", size=15),
          axis.text.y=element_text(colour="black", size=13),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) +
    scale_x_discrete(labels = c("Both proteins")) +
    labs(x = "Both proteins", y = "Difference between median FDR (-log10) values")
ggsave("fdr_diff_bothProteins.png", plot=both_proteins_interactions_fdr_diff, path="./plots/protein_associations_tcga_ccle/")
unlink("fdr_diff_bothProteins.png")





# distribution of effect sizes and FDR by control status
protein_list_a <- protein_list %>%
    left_join(protein_associations_cnv_mrna[, c("controlling", "beta_cnv", "fdr_cnv")], by=c("gene" = "controlling")) %>%
    na.exclude()

protein_list_b <- protein_list %>%
    left_join(protein_associations_cnv_mrna[, c("controlled", "beta_cnv", "fdr_cnv")], by=c("gene" = "controlled")) %>%
    na.exclude()

protein_list_c <- bind_rows(protein_list_a, protein_list_b) %>%
    arrange(gene) %>%
    group_by(gene) %>%
    summarise(median_CNV_beta = median(beta_cnv), median_FDR = median(fdr_cnv), control_status = unique(control_status)) %>%
    gather(key = "var", value = "value", -c(gene, control_status))


effect_size_fdr_distribution <- ggplot(data = protein_list_c, mapping = aes(x = control_status, y = value, fill=control_status)) +
    geom_boxplot() +
    facet_grid(var ~ ., scales = "free") +
    theme(axis.title.y=element_text(colour="black", size=13), axis.title.x=element_text(colour="black", size=13)) +
    theme(axis.text.y=element_text(colour="black", size=10), axis.text.x=element_text(colour="black", size=10)) +
    theme(strip.text.x = element_text(size=10)) +
    labs(x = "", y = "")
ggsave(filename="effect_size_fdr_distribution_control_status.png", plot=effect_size_fdr_distribution, path = "./plots/protein_associations_tcga_ccle/")
unlink("effect_size_fdr_distribution_control_status.png")






# number of interactions by control status
protein_list_Ninteractions <- protein_list %>%
    group_by(control_status, interactions) %>%
    summarise(counts = n()) %>%
    ungroup()

protein_list_Ninteractions_5 <- protein_list_Ninteractions %>%
    group_by(control_status) %>%
    filter(interactions >= 5) %>%
    summarise(interactions = ">=5", counts=sum(counts))

protein_list_Ninteractions <- protein_list_Ninteractions %>%
    group_by(control_status) %>%
    filter(interactions < 5) %>%
    mutate(interactions = as.character(interactions)) %>%
    bind_rows(protein_list_Ninteractions_5) %>%
    arrange(control_status) %>%
    mutate(perc = counts/sum(counts)*100) %>%
    ungroup()


control_status_number_interactions <- ggplot(data=protein_list_Ninteractions) +
    geom_bar(mapping=aes(x=control_status, y=perc, fill=interactions), stat = "identity", position = "dodge") +
    theme(axis.title.y=element_text(colour="black", size=15),
          axis.title.x=element_text(colour="black", size=15),
          axis.text.y=element_text(colour="black", size=13),
          axis.text.x=element_text(colour="black", size=13)) +
    labs(x = "Control status", y = "Percentage of interactions")
ggsave("number_interactions_controlling_controlled_both_proteins.png", plot=control_status_number_interactions, path="./plots/protein_associations_tcga_ccle/")
unlink("number_interactions_controlling_controlled_both_proteins.png")





# positive and negative interactions
positive_interactions <- protein_associations_cnv_mrna %>%
    filter(beta_cnv > 0) %>%
    dplyr::select(pair, controlling, controlled)

negative_interactions <- protein_associations_cnv_mrna %>%
    filter(beta_cnv < 0) %>%
    dplyr::select(pair, controlling, controlled)



# mutual interactions X <-> Y
mutual_interactions <- protein_associations_cnv_mrna %>%
    dplyr::select(controlling, controlled)

mutual_interactions <- mutual_interactions %>%
    bind_rows(setNames(dplyr::select(mutual_interactions, controlled, controlling), c("controlling", "controlled"))) %>%
    mutate(pair=paste(controlling, controlled, sep="_")) %>%
    group_by(pair) %>%
    summarise(counts=n()) %>%
    filter(counts>1) %>%
    dplyr::select(pair) %>%
    arrange(pair)
nrow(mutual_interactions)
#20


# number of times a single protein acts as controlling
interactions_controlling <- protein_associations_cnv_mrna %>%
    dplyr::select(controlling, controlled) %>%
    group_by(controlling) %>%
    summarise(counts=n()) %>%
    group_by(counts) %>%
    summarise(freq=n())


# number of times a single protein acts as controlled
interactions_controlled <- protein_associations_cnv_mrna %>%
    dplyr::select(controlling, controlled) %>%
    group_by(controlled) %>%
    summarise(counts=n()) %>%
    group_by(counts) %>%
    summarise(freq=n())




ppairs_distribution_by_controlling <- ggplot(data=interactions_controlling, mapping=aes(x=as.factor(counts), y=freq)) +
    geom_bar(fill="#386cb0", stat="identity") +
    labs(x = "Number of Times a Single Protein Controls other Proteins", y = "Number of Interactions") +
    theme(axis.title.y=element_text(colour="black", size=15),
          axis.title.x=element_text(colour="black", size=15),
          axis.text.y=element_text(colour="black", size=13),
          axis.text.x=element_text(colour="black", size=13),
          legend.position="none")
ggsave(filename="ppairs_distribution_by_controlling.png", plot=ppairs_distribution_by_controlling, path="./plots/protein_associations_tcga_ccle/")
unlink("ppairs_distribution_by_controlling.png")


ppairs_distribution_by_controlled <- ggplot(data=interactions_controlled, mapping=aes(x=as.factor(counts), y=freq)) +
    geom_bar(fill="#386cb0", stat="identity") +
    labs(x = "Number of Times a Single Protein is Controlled by other Proteins", y = "Number of Interactions") +
    theme(axis.title.y=element_text(colour="black", size=15),
          axis.title.x=element_text(colour="black", size=15),
          axis.text.y=element_text(colour="black", size=13),
          axis.text.x=element_text(colour="black", size=13),
          legend.position="none")
ggsave(filename="ppairs_distribution_by_controlled.png", plot=ppairs_distribution_by_controlled, path="./plots/protein_associations_tcga_ccle/")
unlink("ppairs_distribution_by_controlled.png")









# proteins controlling at least 5 proteins
interactions_controlling_5 <- protein_associations_cnv_mrna %>%
    dplyr::select(controlling, controlled) %>%
    group_by(controlling) %>%
    summarise(counts=n()) %>%
    filter(counts >= 5)
write.table(as.data.frame(interactions_controlling_5), "./files/proteins_controlling_5_more_cnv.txt", row.names=F, quote=F, sep="\t")




# nested interactions
# X -> Y & Z -> X: Z -> X -> Y
nested_interactions_number <- protein_associations_cnv_mrna %>%
    dplyr::select(controlling, controlled) %>%
    mutate(Var=if_else(controlling %in% controlled, 1, 0)) %>%
    group_by(Var) %>%
    summarise(counts = n())

nested_interactions <- protein_associations_cnv_mrna %>%
    dplyr::select(controlling, controlled) %>%
    mutate(control_controlling=sapply(controlling, function(x) paste(controlling[which(controlled == x)], collapse="|"))) %>%
    filter(control_controlling != "")
write.table(as.data.frame(nested_interactions), "./files/nested_interactions_cnv.txt", row.names=F, quote=F, sep="\t")



# protein interactions graph
p_pairs <- protein_associations_cnv_mrna %>%
    dplyr::select(controlling, controlled) %>%
    group_by(controlling) %>%
    filter(n() >= 2) %>%
    arrange(controlling)
p_pairs_graph <- graph_from_data_frame(d = as.data.frame(p_pairs), directed = TRUE)

pdf("./plots/protein_associations_tcga_ccle/p_pairs_graph_controlling_2more.pdf")
par(mar=c(0,0,0,0))
plot(p_pairs_graph,
     layout=layout_with_fr,
     vertex.color="lightblue",
     vertex.frame.color= "black",
     vertex.label.color="black",
     vertex.size=3,
     vertex.label = V(p_pairs_graph)$name,
     vertex.label.cex = 0.3,
     edge.color="grey",
     edge.width=0.2,
     edge.arrow.size=0.3,
     edge.arrow.width=0.6)
dev.off()



# enrichment analysis

protein_set <- protein_list$gene
protein_set_controlling <- protein_list %>% filter(control_status == "controlling") %>% pull(gene)
protein_set_controlled <- protein_list %>% filter(control_status == "controlled") %>% pull(gene)
protein_set_both <- protein_list %>% filter(control_status == "both") %>% pull(gene)
universe_set <- unique(c(protein_associations_cnv$controlling, protein_associations_cnv$controlled))



# GO over-representation test
go_test <- function(gene_set, universe, ontology, p_adj, q_value){
  enr <- enrichGO(
    gene = gene_set,
    universe = universe,
    ont = ontology,
    OrgDb = org.Hs.eg.db,
    keyType = "SYMBOL",
    pvalueCutoff = p_adj,
    qvalueCutoff = q_value,
    pAdjustMethod = "BH",
    minGSSize = 5,
    maxGSSize = 500)
}


protein_set_controlling_overRepGO <- go_test(protein_set_controlling, universe_set, "ALL", 0.05, 1)
protein_set_controlled_overRepGO <- go_test(protein_set_controlled, universe_set, "ALL", 0.05, 1)
protein_set_both_overRepGO <- go_test(protein_set_both, universe_set, "ALL", 0.05, 1)


protein_set_controlling_overRepGO_dot <- dotplot(protein_set_controlling_overRepGO, x = "geneRatio", colorBy = "p.adjust", showCategory = 5, split = "ONTOLOGY", font.size = 20, title = "Controlling proteins\nTop 5 GO categories for BP, MF and CC")
ggsave(filename="protein_set_controlling_over_representation_GO.png", plot=protein_set_controlling_overRepGO_dot, path="./plots/protein_associations_tcga_ccle/", width = 15, height = 12)
unlink("protein_set_controlling_over_representation_GO.png")

protein_set_controlled_overRepGO_dot <- dotplot(protein_set_controlled_overRepGO, x = "geneRatio", colorBy = "p.adjust", showCategory = 5, split = "ONTOLOGY", font.size = 20, title = "Controlled proteins\nTop 5 GO categories for BP, MF and CC")
ggsave(filename="protein_set_controlled_over_representation_GO.png", plot=protein_set_controlled_overRepGO_dot, path="./plots/protein_associations_tcga_ccle/", width = 15, height = 12)
unlink("protein_set_controlled_over_representation_GO.png")

protein_set_both_overRepGO_dot <- dotplot(protein_set_both_overRepGO, x = "geneRatio", colorBy = "p.adjust", showCategory = 5, split = "ONTOLOGY", font.size = 20, title = "Controlling and controlled proteins\nTop 5 GO categories for BP, MF and CC")
ggsave(filename="protein_set_both_over_representation_GO.png", plot=protein_set_both_overRepGO_dot, path="./plots/protein_associations_tcga_ccle/", width = 15, height = 12)
unlink("protein_set_both_over_representation_GO.png")





# load corum data
corum <- read_tsv("./data/corum/coreComplexes_29_05_2018.txt")


overlap_complexes <- function(complexes){
  x <- rep(1, length(complexes))

  for(i in 1:(length(complexes)-1)){
    for(j in (i+1):length(complexes)){
      if(x[i] == 1 & x[j] == 1){
        a <- complexes[[i]]
        b <- complexes[[j]]
        jac <- (length(intersect(a, b)))/(length(union(a, b)))
        if(jac >= 0.9){
          x[j] = 0
        }
      }
    }
  }
  return(x)
}


# protein subunits by corum complex
corum_complexes <- corum %>%
	filter(Organism == "Human") %>%
	dplyr::select(ComplexName, `subunits(Gene name)`) %>%
	dplyr::rename(GeneName = `subunits(Gene name)`) %>%
	group_by(ComplexName) %>%
	summarise(GeneName = paste(GeneName, collapse=";")) %>%
	ungroup() %>%
	mutate(GeneName = str_split(GeneName, ";\\s?")) %>%
	unnest() %>%
	filter(GeneName != "") %>%
	distinct() %>%
  group_by(ComplexName) %>%
  summarise(GeneName = list(GeneName)) %>%
  ungroup() %>%
  mutate(keep = overlap_complexes(GeneName)) %>%
  filter(keep == 1) %>%
  dplyr::select(-keep) %>%
  unnest(GeneName)


corum_test <- function(gene_set, universe, terms, p_adj, q_value){
  enr <- enricher(
    gene = gene_set,
    universe = universe,
    TERM2GENE = terms,
    pvalueCutoff = p_adj,
    qvalueCutoff = q_value,
    pAdjustMethod = "BH",
    minGSSize = 5,
    maxGSSize = 500)
}


both_CORUMenr <- corum_test(protein_set_both, universe_set, corum_complexes, 0.05, 1)
controlling_CORUMenr <- corum_test(protein_set_controlling, universe_set, corum_complexes, 0.05, 1)
controlled_CORUMenr <- corum_test(protein_set_controlled, universe_set, corum_complexes, 0.05, 1)





# analysis of protein complexes

#OST
ost_proteins <- c("DDOST", "KRTCAP2", "MLEC", "OSTC", "RPN1", "DAD1", "TUSC3", "MAGT1", "OST4", "RPN2", "STT3A", "STT3B")
#EMC
emc_proteins <- c("MMGT1", "EMC5", "EMC4", "EMC1", "EMC3", "EMC7", "EMC10", "EMC9", "EMC8", "EMC6", "EMC2")

# with synonyms
#ost_proteins <- c("DDOST", "OST", "GATD6", "KIAA0115", "OST48", "WBP1", "KRTCAP2", "KCP2", "MLEC", "KIAA0152", "OSTC", "DC2", "RPN1", "OST1", "DAD1", "OST2", "TUSC3", "MAGT2", "MGC13453", "MRT7", "N33", "OST3A", "MAGT1", "DKFZP564K142", "IAP", "MRX95", "OST3B", "OST4", "RPN2", "RIBIIR", "RPN-II", "RPNII", "SWP1", "STT3A", "ITM1", "STT3", "MGC9042", "STT3-A", "TMC", "STT3B", "FLJ90106", "SIMP", "STT3-B")
#emc_proteins <- c("MMGT1", "EMC5", "EMC4", "FLJ90746", "MGC24415", "PIG17", "EMC1", "EMC3", "EMC7", "C11ORF3", "EMC10", "HSM1", "HSS1", "INM02", "EMC9", "CGI-112", "EMC8", "FAM158B", "EMC6", "MGC2963", "RAB5IFL", "EMC2")


# Spliceosome
spl_proteins <- read_delim("./data/corum/spliceosome_complex.txt", delim="\t") %>%
    filter(ComplexName == "Spliceosome" & Organism == "Human") %>%
    pull(`subunits(Gene name)`)

spl_proteins <- spl_proteins %>%
    strsplit(split=";") %>%
    unlist()


ost_interactions <- protein_associations_cnv_mrna %>%
    dplyr::select(controlling, controlled) %>%
    mutate(OST = controlling %in% ost_proteins | controlled %in% ost_proteins) %>%
    filter(OST)

emc_interactions <- protein_associations_cnv_mrna %>%
    dplyr::select(controlling, controlled) %>%
    mutate(EMC = controlling %in% emc_proteins | controlled %in% emc_proteins) %>%
    filter(EMC)

spl_interactions <- protein_associations_cnv_mrna %>%
    dplyr::select(controlling, controlled) %>%
    mutate(SPL = controlling %in% spl_proteins | controlled %in% spl_proteins) %>%
    filter(SPL)


# SMARCC1 and SMARCC2
negative_interactions %>%
    filter(controlling %in% c("SMARCC1", "SMARCC2") | controlled %in% c("SMARCC1", "SMARCC2") )

positive_interactions %>%
    filter(controlling %in% c("SMARCC1", "SMARCC2") | controlled %in% c("SMARCC1", "SMARCC2") )


# endoplasmic reticulum membrane proteins
ER <- read_tsv("./data/protein_complexes/subunits_of_multisubunit_complexes_ER_INM.txt") %>%
    mutate(gene_name = toupper(gene_name))

# all pairs within protein complexes
ER_pairs <- ER %>%
    group_by(complex_name) %>%
    mutate(pair = list(tidyr::crossing(gene_name, gene_name))) %>%
    ungroup() %>%
    unnest() %>%
    dplyr::select("gene_name1", "gene_name2") %>%
    distinct() %>%
    filter(gene_name1 != gene_name2) %>%
    dplyr::rename(p1=gene_name1, p2=gene_name2) %>%
    mutate(pair = paste(p1, p2, sep="_"))

# intersection with significative protein associations
ER_pairs %>%
    pull(pair) %>%
    intersect(protein_associations_cnv_mrna$pair)








save(list=ls(), file="./r_workspaces/protein_associations_tcga_ccle_characterization.RData")
