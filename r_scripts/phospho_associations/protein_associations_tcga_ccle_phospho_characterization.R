# Abel Sousa 2018

# PhD project
# Study of Control of Protein Abundance using Multi-omics Data from Cancer Samples



# -- Exploration of protein associations using phosphoproteomics



suppressMessages(library(tidyverse))
suppressMessages(library(ggrepel))
suppressMessages(library(RColorBrewer))
suppressMessages(library(clusterProfiler))
suppressMessages(library(org.Hs.eg.db))
suppressMessages(library(igraph))
suppressMessages(library(ggpubr))
options(bitmapType = "cairo")




# import protein associations with CNV and RNA
# filtered by genomic co-localization
signf_protein_associations_cnv_rna <- read_tsv("./files/signf_protein_associations_cnv_mrna_filtered_close_controlling.txt") %>%
    mutate(pair = paste(controlling, controlled, sep="_")) %>%
    dplyr::select(pair, everything())



protein_associations_phospho <- read_tsv("./files/protein_associations_tcga_ccle_phospho_reg.txt") %>%
    mutate(pair = paste(controlling, controlled, sep="_")) %>%
    dplyr::select(pair, controlling, controlling_phx, controlled, beta_phospho, fdr_phospho = fdr) %>%
    mutate(fdr_label = if_else(fdr_phospho > 0.05, "fdr > 0.05", "fdr <= 0.05"))



# significant phospho associations at FDR < 0.05
signf_phospho <- protein_associations_phospho %>%
    dplyr::select(pair, controlling, controlling_phx, controlled, beta_phospho, fdr_phospho) %>%
    filter(fdr_phospho<0.05)





# significant associations with CNV/RNA and Phospho
signf_cnv_rna_phospho <- inner_join(signf_protein_associations_cnv_rna, signf_phospho, by=c("pair", "controlling", "controlled")) %>%
    dplyr::select(pair, controlling, controlling_phx, controlled, everything()) %>%
    dplyr::arrange(fdr_phospho)
write.table(signf_cnv_rna_phospho, "./files/signf_protein_associations_cnv_rna_phospho.txt", row.names=F, col.names=T, sep="\t", quote=F)






# volcano plots
# plot phospho beta by fdr
# just common associations with CNV/RNA model


# associations with phospho and present in the significative associations with CNV/RNA
phospho_associations <- protein_associations_phospho %>%
    inner_join(signf_protein_associations_cnv_rna[, c(1,2,3)], by=c("pair", "controlling", "controlled"))


protein_associations_cnv_rna_phospho_plot <- ggplot( data=phospho_associations, mapping=aes(x=beta_phospho, y=-log10(fdr_phospho), colour=fdr_label) ) +
    geom_point() +
    theme_classic() +
    geom_text_repel(
      data=bind_rows(
        phospho_associations %>% filter(fdr_phospho < 0.05, beta_phospho > 0) %>% arrange(fdr_phospho) %>% head(5),
        phospho_associations %>% filter(fdr_phospho < 0.05, beta_phospho > 0) %>% filter(controlling_phx == "POLD3_s458"),
        phospho_associations %>% filter(fdr_phospho < 0.05, beta_phospho < 0) %>% arrange(fdr_phospho) %>% head(1)),
      mapping=aes(x=beta_phospho, y=-log10(fdr_phospho),
      label = paste(controlling_phx, controlled, sep="~")), size = 3.2, box.padding = unit(0.3, "lines"), point.padding = unit(0.3, "lines"), colour="black" ) +
    scale_colour_manual(values=c("#2171B5", "#BDD7E7"), name = "", labels = c("FDR <= 0.05", "FDR > 0.05"), guide=F) +
    scale_x_continuous(limits = c(-0.35, 0.5), breaks=c(-0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5)) +
    geom_line(aes(x=0),colour="black", linetype=2, size = 0.3) +
    geom_line(aes(y=-log10(0.05)),colour="black", linetype=2, size = 0.3) +
    theme(axis.title=element_text(colour="black", size=15),
          axis.text=element_text(colour="black", size=12),
          plot.title = element_blank(),
          legend.text=element_text(colour="black", size=12),
          legend.title=element_blank()) +
    labs(x = "Phospho Beta", y = "Phospho model FDR (-log10)")
ggsave(filename="protein_associations_tcga_ccle_cnv_phospho.png", plot=protein_associations_cnv_rna_phospho_plot, width=5, height=6, path="./plots/protein_associations_tcga_ccle_phospho/")
ggsave(filename="protein_associations_tcga_ccle_cnv_phospho.pdf", plot=protein_associations_cnv_rna_phospho_plot, width=5, height=6, path="./plots/protein_associations_tcga_ccle_phospho/")
unlink("protein_associations_tcga_ccle_cnv_phospho.png")
unlink("protein_associations_tcga_ccle_cnv_phospho.pdf")










# number of times a protein is controlled
# number of times a phosphosite/gene is controlling
signfPhospho_Proteincontrolled <- signf_phospho %>%
    group_by(controlled) %>%
    summarise(controlled_times = n())


signfPhospho_Genecontrolling <- signf_phospho %>%
    group_by(controlling) %>%
    summarise(controlling_times_gene = n())


signfPhospho_Phosphocontrolling <- signf_phospho %>%
    group_by(controlling_phx) %>%
    summarise(controlling_times_phospho = n())


signf_phospho2 <- full_join(signfPhospho_Proteincontrolled, signfPhospho_Genecontrolling, by=c("controlled" = "controlling")) %>%
    full_join(signfPhospho_Phosphocontrolling, by=c("controlled" = "controlling_phx")) %>%
    dplyr::rename(feature=controlled) %>%
    gather(key="stat", value="freq", -feature)


signf_phospho3 <- signf_phospho2 %>%
    filter(!is.na(freq)) %>%
    group_by(stat) %>%
    summarise(number = n(), mean_freq = mean(freq))



signf_phospho2_barplot<- ggplot(data=signf_phospho3, mapping=aes(x = stat, y=mean_freq, fill=stat)) +
    geom_bar(position = "dodge", stat = "identity") +
    theme(axis.title=element_text(colour="black", size=15),
          axis.text.y=element_text(colour="black", size=13),
          axis.text.x=element_text(colour="black", size=13),
          legend.position="none") +
    scale_x_discrete(name = "", labels=c("Controlled proteins", "Controlling proteins", "Controlling phosphosites")) +
    scale_y_continuous(name = "Mean of the number of interactions") +
    scale_fill_discrete(name = "") +
    coord_flip()
ggsave("signf_phospho_controlling_controlled_times.png", plot=signf_phospho2_barplot, path="./plots/protein_associations_tcga_ccle_phospho/")
unlink("signf_phospho_controlling_controlled_times.png")



signf_phospho3_barplot <- ggplot(data=signf_phospho3, mapping=aes(x=stat, y=number, fill=stat)) +
    geom_bar(position = "dodge", stat = "identity") +
    theme(axis.title=element_text(colour="black", size=15),
          axis.text.y=element_text(colour="black", size=13),
          axis.text.x=element_text(colour="black", size=13),
          legend.position="none") +
    scale_x_discrete(name = "", labels=c("Controlled proteins", "Controlling proteins", "Controlling phosphosites")) +
    scale_y_continuous(name = "Number of proteins/phosphosites") +
    scale_fill_discrete(name = "") +
    coord_flip()
ggsave("signf_phospho_controlling_controlled_proteins.png", plot=signf_phospho3_barplot, path="./plots/protein_associations_tcga_ccle_phospho/")
unlink("signf_phospho_controlling_controlled_proteins.png")








# load proteomics and phosphoproteomics data
common_samples <- read_tsv("./files/cptac_cellLines_samples_prot_phospho_rna_cnv.txt", col_names = FALSE)

proteomics <- read_tsv("./files/proteomicsQ_cptac_cellLines.txt") %>%
    gather(key="sample", value="prot_log2FC", -gene) %>%
    arrange(sample) %>%
    filter(sample %in% common_samples$X1) %>%
    group_by(gene) %>%
    filter( sum(!is.na(prot_log2FC)) > n()*0.5 ) %>%
    mutate(prot_log2FC = scale(prot_log2FC)[,1]) %>%
    ungroup()


phospho <- read_tsv("./files/phosphoproteomicsQ_cptac_cellLines.txt") %>%
    gather(key="sample", value="phospho_log2FC", -phospho_site) %>%
    arrange(sample) %>%
    filter(sample %in% common_samples$X1) %>%
    group_by(phospho_site) %>%
    filter( sum(!is.na(phospho_log2FC)) > n()*0.5 ) %>%
    mutate(phospho_log2FC = scale(phospho_log2FC)[,1]) %>%
    ungroup()




# correlate phosphoproteomics and proteomics for the significative phospho associations
var <- signf_phospho %>%
    dplyr::select(controlling, controlling_phx) %>%
    inner_join(proteomics, by=c("controlling" = "gene")) %>%
    inner_join(phospho, by=c("controlling_phx" = "phospho_site", "sample")) %>%
    group_by(controlling, controlling_phx) %>%
    do(broom::tidy(cor.test(.$prot_log2FC, .$phospho_log2FC, method = "pearson"))) %>%
    ungroup() %>%
    dplyr::select(protein=controlling, phosphosite=controlling_phx, cor=estimate, p.value) %>%
    mutate(var = "controlling phospho vs controlling protein FDR < 5%")



# correlate phosphoproteomics and proteomics for all the other possible proteomics/phosphoproteomics pairs
gene_phosphosites <- phospho %>%
    dplyr::select(phospho_site) %>%
    distinct() %>%
    mutate(gene = str_split_fixed(string=phospho_site, pattern="_", n=2)[,1]) %>%
    dplyr::select(gene, phospho_site)



var2 <- gene_phosphosites %>%
    inner_join(proteomics, by="gene") %>%
    inner_join(phospho, by=c("phospho_site", "sample")) %>%
    filter(!phospho_site %in% var$phosphosite) %>%
    group_by(gene, phospho_site) %>%
    do(broom::tidy(cor.test(.$prot_log2FC, .$phospho_log2FC, method = "pearson"))) %>%
    ungroup() %>%
    dplyr::select(protein=gene, phosphosite=phospho_site, cor=estimate, p.value) %>%
    mutate(var = "other possible pairs phosphosite / protein")


var2_5 <- rbind(var, var2)
var2_5_sum <- var2_5 %>%
    group_by(var) %>%
    summarise(counts=n()) %>%
    mutate(y = -0.75)


var2_5_boxplot <- ggplot(data=var2_5, mapping=aes(x = var, y=cor, fill=var)) +
    geom_boxplot() +
    geom_text(data = var2_5_sum, mapping=aes(x = var, y = y, label=counts), color = "black") +
    theme(axis.title=element_text(colour="black", size=15),
          axis.text.y=element_text(colour="black", size=13),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) +
    scale_x_discrete(name = "") +
    scale_y_continuous(name = "Pearson correlation coefficient", limits = c(-1,1)) +
    scale_fill_discrete(name = "")
ggsave("phospho_associations_controlling_phosphosite_protein_cor.png", plot=var2_5_boxplot, path="./plots/protein_associations_tcga_ccle_phospho/")
unlink("phospho_associations_controlling_phosphosite_protein_cor.png")





# correlate phosphoproteomics and proteomics between controlling and controlled proteins, for the significative phospho associations
var3 <- signf_phospho %>%
    dplyr::select(controlling_phx, controlled) %>%
    inner_join(phospho, by=c("controlling_phx" = "phospho_site")) %>%
    inner_join(proteomics, by=c("controlled" = "gene", "sample")) %>%
    group_by(controlling_phx, controlled) %>%
    do(broom::tidy(cor.test(.$phospho_log2FC, .$prot_log2FC, method = "pearson"))) %>%
    ungroup() %>%
    dplyr::select(controlling = controlling_phx, controlled, cor=estimate, p.value) %>%
    mutate(var = "Controlling Phospho vs Controlled Protein \n FDR < 5%")




# correlate proteomics between controlling and controlled proteins, for the significative phospho associations
var4 <- signf_phospho %>%
    dplyr::select(controlling, controlled) %>%
    inner_join(proteomics, by=c("controlling" = "gene")) %>%
    inner_join(proteomics, by=c("controlled" = "gene", "sample")) %>%
    group_by(controlling, controlled) %>%
    do(broom::tidy(cor.test(.$prot_log2FC.x, .$prot_log2FC.y, method = "pearson"))) %>%
    ungroup() %>%
    dplyr::select(controlling, controlled, cor=estimate, p.value) %>%
    mutate(var = "Controlling Protein vs Controlled Protein \n FDR < 5%")




# correlate phosphoproteomics and proteomics between controlling and controlled proteins, for the non-significative phospho associations
var5 <- protein_associations_phospho %>%
    filter(fdr>0.05) %>%
    dplyr::select(controlling_phx, controlled) %>%
    inner_join(phospho, by=c("controlling_phx" = "phospho_site")) %>%
    inner_join(proteomics, by=c("controlled" = "gene", "sample")) %>%
    group_by(controlling_phx, controlled) %>%
    do(broom::tidy(cor.test(.$phospho_log2FC, .$prot_log2FC, method = "pearson"))) %>%
    ungroup() %>%
    dplyr::select(controlling = controlling_phx, controlled, cor=estimate, p.value) %>%
    mutate(var = "Controlling Phospho vs Controlled Protein \n FDR > 5%")




# correlate proteomics between controlling and controlled proteins, for the non-significative phospho associations
var6 <- protein_associations_phospho %>%
    filter(fdr>0.05) %>%
    dplyr::select(controlling, controlled) %>%
    inner_join(proteomics, by=c("controlling" = "gene")) %>%
    inner_join(proteomics, by=c("controlled" = "gene", "sample")) %>%
    group_by(controlling, controlled) %>%
    do(broom::tidy(cor.test(.$prot_log2FC.x, .$prot_log2FC.y, method = "pearson"))) %>%
    ungroup() %>%
    dplyr::select(controlling, controlled, cor=estimate, p.value) %>%
    mutate(var = "Controlling Protein vs Controlled Protein \n FDR > 5%")


var7 <- rbind(var3, var4, var5, var6) %>%
    mutate(var = as.factor(var)) %>%
    mutate(var = fct_relevel(var, "Controlling Phospho vs Controlled Protein \n FDR < 5%", "Controlling Protein vs Controlled Protein \n FDR < 5%"))


var7_sum <- var7 %>%
    group_by(var) %>%
    summarise(counts=n()) %>%
    mutate(y = -0.75)


var7_boxplot <- ggplot(data=var7, mapping=aes(x = var, y=cor, fill=var)) +
    geom_boxplot() +
    geom_text(data = var7_sum, mapping=aes(x = var, y = y, label=counts), color = "black") +
    stat_compare_means(comparisons = list( c(1, 2), c(3, 4), c(1, 3), c(1, 4)) ) +
    theme(axis.title=element_text(colour="black", size=15),
          axis.text.y=element_text(colour="black", size=13),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) +
    scale_x_discrete(name = "") +
    scale_y_continuous(name = "Pearson correlation coefficient") +
    scale_fill_discrete(name = "")
ggsave("phospho_associations_controlling_phosphosite_controlled_protein_cor.png", plot=var7_boxplot, path="./plots/protein_associations_tcga_ccle_phospho/")
unlink("phospho_associations_controlling_phosphosite_controlled_protein_cor.png")







# emanuel phospho associations
emanuel_phospho <- read_csv("./data/emanuel/regressions_protein_phosphoproteomics.csv") %>%
    mutate(fdr2 = p.adjust(pval, method="BH"))

emanuel_phospho_signf <- emanuel_phospho %>%
    mutate(pair = paste(PHx, Py, sep="_")) %>%
    dplyr::select(pair, everything()) %>%
    filter(fdr2 <= 0.05)
#36,893


# comparison
intersect(emanuel_phospho_signf$pair, toupper(paste(signf_phospho$controlling_phx, signf_phospho$controlled, sep="_"))) %>% length
#3519
#30%




save(list=ls(), file="./r_workspaces/protein_associations_tcga_ccle_phospho_characterization.RData")
