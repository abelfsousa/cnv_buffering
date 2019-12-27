# Abel Sousa 2018

# PhD project
# Study of Control of Protein Abundance using Multi-omics Data from Cancer Samples

# -- Number of proteins controlled by at least two proteins on the same chromosome

# -- Correlation of protein pairs on the same chromosome


suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
options(bitmapType = "cairo")




# load datasets


# load genomic coordinates
load("/nfs/research1/beltrao/memon/databases/protein/ensembl/ensembl_v73/ensembl_gtf_v73_gene_anno_grange.Rdata")
ensembl_gtf_v73_gene_anno_grange <- as.tibble(as.data.frame(ensembl_gtf_v73_gene_anno_grange)) %>%
    group_by(gene_name) %>%
    filter(n() == 1) %>%
    ungroup()


# cnv correlation
cnv_correlation <- fread("./files/correlation_cnv.txt") %>% as.tibble()



# protein associations with CNV
protein_associations_cnv <- read_tsv("./files/protein_associations_tcga_ccle_cnv_reg.txt")


# significant protein associations with CNV
signf_protein_associations_cnv <- protein_associations_cnv %>%
    filter(fdr <= 0.05)


# protein associations with mRNA
protein_associations_mRNA <- read_tsv("./files/protein_associations_tcga_ccle_mRNA_reg.txt")


# significant protein associations with mRNA
signf_protein_associations_mRNA <- protein_associations_mRNA %>%
    filter(fdr <= 0.05)




# significant protein associations with CNV and mRNA
signf_protein_associations_cnv_mrna <- signf_protein_associations_cnv[, c("controlling", "controlled", "beta_cnv", "fdr")] %>%
    inner_join(protein_associations_mRNA[, c("controlling", "controlled", "beta_rna", "lrt_p")], by=c("controlling", "controlled")) %>%
    dplyr::rename(fdr_cnv = fdr) %>%
    mutate(fdr_rna = p.adjust(p=lrt_p, method="BH")) %>%
    filter(fdr_rna < 0.05) %>%
    dplyr::select(-lrt_p)
write.table(signf_protein_associations_cnv_mrna, file = "./files/signf_protein_associations_cnv_mrna.txt", sep = "\t", quote = F, row.names=F)




# number of times each protein is found and control status
protein_list <- signf_protein_associations_cnv_mrna %>%
    dplyr::select(controlling, controlled) %>%
    dplyr::combine() %>%
    as_tibble() %>%
    dplyr::rename(gene = value) %>%
    group_by(gene) %>%
    dplyr::summarise(interactions = n()) %>%
    arrange(desc(interactions)) %>%
    dplyr::mutate(control_status = if_else(
        gene %in% setdiff(signf_protein_associations_cnv_mrna$controlling, signf_protein_associations_cnv_mrna$controlled), "controlling",
        if_else(gene %in% setdiff(signf_protein_associations_cnv_mrna$controlled, signf_protein_associations_cnv_mrna$controlling), "controlled", "both")))





# proteins being controlled by at least 2 proteins on the same chromosome
proteins_controlled_2more <- signf_protein_associations_cnv_mrna %>%
    dplyr::select(controlling, controlled) %>%
    group_by(controlled) %>%
    summarise(controlled_times=n()) %>%
    filter(controlled_times > 1) %>%
    inner_join(signf_protein_associations_cnv_mrna[, c(1,2)], by="controlled") %>%
    left_join(ensembl_gtf_v73_gene_anno_grange[, c(6,1,2,3,5)], by=c("controlling" = "gene_name")) %>%
    dplyr::rename(chr=seqnames) %>%
    dplyr::select(controlled, controlled_times, controlling, chr) %>%
    group_by(controlled, chr) %>%
    #summarise(counts = n()) %>%
    #filter(counts > 1)
    filter(n()>1) %>%
    mutate(controlled_times=n()) %>%
    ungroup() %>%
    arrange(desc(controlled_times), chr, controlled, controlling)
write.table(as.data.frame(proteins_controlled_2more), "./files/proteins_controlled_2more_sameChr.txt", row.names=F, quote=F, sep="\t")




# proteins being controlled by at least 2 proteins on the same chromosome
# CNV correlation by protein pairs on the same chromosome
proteins_controlled_2more_pairsCor <- signf_protein_associations_cnv_mrna %>%
    dplyr::select(controlling, controlled) %>%
    group_by(controlled) %>%
    summarise(controlled_times=n()) %>%
    filter(controlled_times > 1) %>%
    inner_join(signf_protein_associations_cnv_mrna[, c(1,2)], by="controlled") %>%
    left_join(ensembl_gtf_v73_gene_anno_grange[, c(6,1,2,3,5)], by=c("controlling" = "gene_name")) %>%
    dplyr::rename(chr=seqnames) %>%
    dplyr::select(controlled, controlled_times, controlling, chr) %>%
    arrange(desc(controlled_times), chr, controlled, controlling) %>%
    group_by(controlled, chr) %>%
    filter(n()>1) %>%
    mutate(controlled_times=n()) %>%
	mutate(pair = list(as.tibble(t(combn(controlling, 2))))) %>%
	ungroup() %>%
	unnest() %>%
	dplyr::rename(pair1 = V1, pair2 = V2) %>%
	dplyr::select(-controlling) %>%
	distinct() %>%
	arrange(desc(controlled_times), chr, controlled) %>%
	inner_join(cnv_correlation, by=c("chr", "pair1"="gene1", "pair2"="gene2")) %>%
	mutate(controlled_chr = as.factor(paste(controlled, chr, sep="_"))) %>%
	mutate(controlled_chr = fct_infreq(controlled_chr)) %>%
	dplyr::select(controlled_chr, everything())
write.table(as.data.frame(proteins_controlled_2more_pairsCor), "./files/proteins_controlled_2more_pairsCor.txt", row.names=F, quote=F, sep="\t")





# scatterplot
proteins_controlled_2more_pairsCor_plot <- ggplot( data=proteins_controlled_2more_pairsCor, mapping=aes(x=controlled_chr, y=pearson) ) +
    geom_point() +
    theme_classic() +
    geom_hline(yintercept=0.5, linetype="dashed", color = "red") +
    theme(axis.title.y=element_text(colour="black", size=15),
          axis.title.x=element_text(colour="black", size=15),
          axis.text.y=element_text(colour="black", size=13),
          axis.text.x=element_text(colour="black", size=5, angle=90, hjust=1),
          plot.title = element_text(size=13)) +
    labs(x = "Controlled proteins", y = "CNV profile correlation of controlling protein pairs (r)",
          title=paste(paste0(length(unique(proteins_controlled_2more_pairsCor$controlled)), " proteins controlled by at least 2 proteins on the same chromosome"), paste0(length(unique(proteins_controlled_2more_pairsCor$controlled_chr)), " pairs controlled protein - controlling proteins chromosome"), sep="\n"))
ggsave(filename="proteins_controlled_2more_sameChr_pairsCor.png", plot=proteins_controlled_2more_pairsCor_plot, path="./plots/protein_associations_tcga_ccle/")
unlink("proteins_controlled_2more_sameChr_pairsCor.png")





save(list=ls(), file="./r_workspaces/protein_associations_tcga_ccle_controlling_genomic_proximity.RData")
