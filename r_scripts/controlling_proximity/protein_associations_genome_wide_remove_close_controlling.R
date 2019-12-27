# Abel Sousa 2018

# PhD project
# Study of Control of Protein Abundance using Multi-omics Data from Cancer Samples


# -- Filter cases where multiple controlling proteins are close in the genome


suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(TopKLists))
options(bitmapType = "cairo")




# load datasets


# significative protein associations with CNV and mRNA
signf_protein_associations_cnv_mrna <- read_tsv("./files/signf_protein_associations_cnv_mrna_genomeWide_fdr10.txt")



# proteins being controlled by at least 2 proteins on the same chromosome
proteins_controlled_2more <- read_tsv("./files/proteins_controlled_2more_sameChr_genomeWide_fdr10.txt")



# proteins being controlled by at least 2 proteins on the same chromosome
# CNV correlation by protein pairs on the same chromosome
proteins_controlled_2more_pairsCor <- read_tsv("./files/proteins_controlled_2more_pairsCor_genomeWide_fdr10.txt")



#CNV correlation
cnv_correlation <- fread("./files/correlation_cnv.txt") %>% as.tibble()



# load genomic coordinates
load("/nfs/research1/beltrao/memon/databases/protein/ensembl/ensembl_v73/ensembl_gtf_v73_gene_anno_grange.Rdata")
ensembl_gtf_v73_gene_anno_grange <- as.tibble(as.data.frame(ensembl_gtf_v73_gene_anno_grange)) %>%
	group_by(gene_name) %>%
	filter(n() == 1) %>%
	ungroup()




keep_controlling <- function(chr, controlling_proteins, cnv_correlation){
	#

	if(length(controlling_proteins) == 1){
		out <- 1
	}
	else{

		l <- list()
		for(i in 2:length(controlling_proteins)){
			cor <- c()
			for(j in 1:(i-1)){
				cor <- c(cor, cnv_correlation %>% filter(chr == chr, gene1 == controlling_proteins[i], gene2 == controlling_proteins[j]) %>% pull(pearson) %>% as.numeric)
			}
			l <- c(l, list(cor))
		}

		# use 0.5 as correlation cutoff
		out <- sapply(l, function(x) sum(sum(x > 0.5) == 0))
		out <- c(1, out)
	}

	return(out)

}



# combine CNV and mRNA pvalues for the same protein associations
# rank the associations based on the two p-values
cnv_mRNA_pvalue <- proteins_controlled_2more %>%
	inner_join(signf_protein_associations_cnv_mrna[, c("controlling", "controlled", "fdr_cnv", "fdr_rna")], by=c("controlling", "controlled")) %>%
	mutate(fdr_cnv = -log10(fdr_cnv), fdr_rna = -log10(fdr_rna)) %>%
	group_by(controlled, chr) %>%
	mutate(borda = Borda(list(order(fdr_cnv, decreasing=T), order(fdr_rna, decreasing=T)))$TopK[, "mean"]) %>%
	mutate(order_borda=1:n()) %>%
	mutate(order_borda = match(order_borda, borda)) %>%
	ungroup() %>%
	dplyr::arrange(desc(controlled_times), chr, controlled, order_borda) %>%
	dplyr::select(-c(controlled_times, borda, order_borda, fdr_cnv, fdr_rna)) %>%
	group_by(controlled, chr) %>%
	mutate(keep = keep_controlling(unique(chr), as.character(controlling), cnv_correlation))





# for the controlled proteins with more than one controlling protein on the same chromosome
# keep the non-correlated controlling proteins (pearson < 0.5)
signf_protein_associations_cnv_mrna <- signf_protein_associations_cnv_mrna %>%
	left_join(cnv_mRNA_pvalue[, c("controlling", "controlled", "keep")], by=c("controlling", "controlled")) %>%
	mutate(keep = replace_na(keep, 1)) %>%
	filter(keep == 1) %>%
	dplyr::select(-keep)
write.table(signf_protein_associations_cnv_mrna, file = "./files/signf_protein_associations_cnv_mrna_filtered_close_controlling_genomeWide_fdr10.txt", sep = "\t", quote = F, row.names=F)













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





# scatterplot
proteins_controlled_2more_pairsCor_plot <- ggplot( data=proteins_controlled_2more_pairsCor, mapping=aes(x=controlled_chr, y=pearson) ) +
    geom_point() +
    theme_classic() +
    geom_hline(yintercept=0.5, linetype="dashed", color = "red") +
    theme(axis.title.y=element_text(colour="black", size=15),
          axis.title.x=element_text(colour="black", size=15),
          axis.text.y=element_text(colour="black", size=13),
          axis.text.x=element_text(colour="black", size=13, angle=90, hjust=1),
          plot.title = element_text(size=13)) +
		labs(x = "Controlled proteins", y = "CNV profile correlation of controlling protein pairs (r)",
			title=paste(paste0(length(unique(proteins_controlled_2more_pairsCor$controlled)), " proteins controlled by at least 2 proteins on the same chromosome"), paste0(length(unique(proteins_controlled_2more_pairsCor$controlled_chr)), " pairs controlled protein - controlling proteins chromosome"), sep="\n"))
ggsave(filename="proteins_controlled_2more_sameChr_pairsCor2_genomeWide_fdr10.png", plot=proteins_controlled_2more_pairsCor_plot, path="./plots/protein_associations_genome_wide/")
unlink("proteins_controlled_2more_sameChr_pairsCor2_genomeWide_fdr10.png")








#save(list=ls(), file="./r_workspaces/protein_associations_genome_wide_remove_close_controlling.RData")
