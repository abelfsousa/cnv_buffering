# Abel Sousa 2018

# PhD project
# Study of Control of Protein Abundance using Multi-omics Data from Cancer Samples

# -- CNV correlation across TCGA and CCLE samples

# -- Genes by chromosome


suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
options(bitmapType = "cairo")



# load datasets

common_samples <- read_tsv(file = "./files/cptac_cellLines_samples_prot_rna_cnv.txt", col_names = FALSE)


load("/nfs/research1/beltrao/memon/databases/protein/ensembl/ensembl_v73/ensembl_gtf_v73_gene_anno_grange.Rdata")
ensembl_gtf_v73_gene_anno_grange <- as.tibble(as.data.frame(ensembl_gtf_v73_gene_anno_grange)) %>%
	group_by(gene_name) %>%
	filter(n() == 1) %>%
	ungroup()


cnv <- read_tsv("./files/cnv_tcga_cellLines.txt") %>%
    gather(key="sample", value="cnv_gistic2", -gene) %>%
    filter(sample %in% common_samples$X1) %>%
    inner_join(ensembl_gtf_v73_gene_anno_grange[, c("gene_name", "seqnames")], by=c("gene" = "gene_name")) %>%
    dplyr::rename(chr = seqnames) %>%
    spread(key="sample", value="cnv_gistic2") %>%
    as.data.frame()

rownames(cnv) <- cnv$gene



correlation_list <- list()

for (i in as.character(unique(cnv$chr)) ){
	temp <- as.data.frame(cor(t(cnv[cnv$chr == i, c(3:370)])))
	temp <- setNames(cbind(i, rownames(temp), temp), c("chr", "gene2", colnames(temp)))
	temp <- gather(temp, key="gene1", value="pearson", -gene2, -chr)
	temp <- temp[, c(1,3,2,4)]
	correlation_list[[match(i, as.character(unique(cnv$chr)))]] <- temp
}


correlation_cnv <- do.call(rbind, correlation_list) %>% as.tibble()
#write.table(as.data.frame(correlation_cnv), "./files/correlation_cnv.txt", row.names=F, quote=F, sep="\t")
fwrite(as.data.frame(correlation_cnv), "./files/correlation_cnv.txt", sep="\t")




save(list=ls(), file="./r_workspaces/correlation_cnv.RData")




