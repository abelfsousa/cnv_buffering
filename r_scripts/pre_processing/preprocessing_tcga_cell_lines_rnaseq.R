# Abel Sousa 2018

# PhD project
# Study of Control of Protein Abundance using Multi-omics Data from Cancer Samples


# -- Pre-processing of TCGA and cancer cell lines RNA-seq data


# -- TCGA RNA-seq data downloaded from GEO with the accession number GSE62944


# -- Useful cell lines database: https://web.expasy.org/cellosaurus/


suppressMessages(library(tidyverse))
suppressMessages(library(limma))
suppressMessages(library(edgeR))
options(bitmapType = "cairo")






# -- TCGA RNA-seq data


# load TCGA samples with proteomics
cptac_samples <- read.delim("./files/cptac_samples.txt", stringsAsFactors=F)


# load read counts
tcga_rna <- read.delim("./data/tcga/GSM1536837_06_01_15_TCGA_24.tumor_Rsubread_FeatureCounts.txt", row.names=c(1))


# load dataframe with cancer type
tcga_cancer_types <- read.delim("./data/tcga/GSE62944_06_01_15_TCGA_24_CancerType_Samples.txt", header=F, stringsAsFactors=F, col.names=c("sample", "cancer"))
tcga_cancer_types$sample <- gsub("-", ".", tcga_cancer_types$sample)
tcga_cancer_types$cancer[tcga_cancer_types$cancer == "COAD" | tcga_cancer_types$cancer == "READ"] <- "COREAD"
rownames(tcga_cancer_types) <- tcga_cancer_types$sample


# select primary tumour samples
tcga_rna <- tcga_rna[, substr(colnames(tcga_rna), 14, 15) == "01"]


# select samples with proteome
tcga_rna <- tcga_rna[, substr(colnames(tcga_rna), 1, 12) %in% cptac_samples$sample]


tcga_cancer_types <- tcga_cancer_types[colnames(tcga_rna), ]
all.equal(rownames(tcga_cancer_types), colnames(tcga_rna))
#TRUE



# remove low expressed genes
# genes with average cpm across samples lower than 1 are removed
tcga_rna <- tcga_rna[ rowMeans( cpm(tcga_rna) ) >= 1, ]



# normalise raw counts by TMM method
tcga_rna <- DGEList(counts=tcga_rna)
tcga_rna <- calcNormFactors(tcga_rna, method="TMM")



# derive log-CPM values from voom
tcga_rna <- voom(tcga_rna, design = model.matrix(~cancer, data=tcga_cancer_types))$E
#tcga_rna <- voom(tcga_rna)$E


# average replicated samples
colnames(tcga_rna) <- substr(colnames(tcga_rna), 1, 12)

tcga_rna_rep <- unique(colnames(tcga_rna)[duplicated(colnames(tcga_rna))])
for( i in tcga_rna_rep ){
	col_index <- which( i == colnames(tcga_rna) )
	new_column <- rowMeans( tcga_rna[, col_index] )
	tcga_rna <- tcga_rna[,-col_index]
	tcga_rna <- cbind(tcga_rna, new_column)
	colnames(tcga_rna)[length(colnames(tcga_rna))] <- i
}
rm(i, col_index, new_column)



# define gene column
tcga_rna <- cbind(data.frame(gene=rownames(tcga_rna), stringsAsFactors=F), tcga_rna)
rownames(tcga_rna) <- NULL


# export rna data
write.table(tcga_rna, file = "./files/tcga_rna_logCPM.txt", quote = FALSE, sep = "\t", row.names=FALSE)



# gather data
tcga_rna_gather <- tcga_rna %>%
    gather(key="sample", value="log2CPM", -gene) %>%
    as.tibble() %>%
    inner_join(cptac_samples, by="sample")


# boxplot by sample
tcga_rna_gather_boxp <- ggplot(tcga_rna_gather, aes(x=sample, y=log2CPM, fill=cancer)) +
    geom_boxplot(outlier.size = 0.1) +
    theme(axis.title.y=element_text(colour="black", size=13), axis.title.x=element_text(colour="black", size=13)) +
    theme(axis.text.y=element_text(colour="black", size=10), axis.text.x=element_blank()) +
    theme(axis.ticks.x=element_blank()) +
    labs(x = "Sample", y = "TCGA RNA log2-counts per million (logCPM)", title = "")
ggsave(filename="tcga_rna_log2CPM.png", plot=tcga_rna_gather_boxp, path = "./plots/", width=20)
unlink("tcga_rna_log2CPM.png")








# -- Cell lines RNA-seq data


# load read counts
ccle_counts <- read.table("./data/ccle/CCLE_RNAseq_081117.reads.gct.txt", h=T, sep="\t", stringsAsFactors=F, skip=c(2), check.names=F)
dim(ccle_counts)
#56318  1021



# cell and tissue map
cells <- colnames(ccle_counts)[3:1021]
tissues <- sub("_", "|", cells, fixed=T)
tissues <- sapply(strsplit(tissues, "|", fixed=T), function(x) x[2])
ccle_cell_lines <- cbind(cells, toupper(str_extract(cells, "^[^_]+")), tissues)
ccle_cell_lines[is.na(ccle_cell_lines[, "tissues"]), "tissues"] <- "OTHER"
ccle_cell_lines <- setNames(as.data.frame(ccle_cell_lines, stringsAsFactors=F), c("cell_id", "cell", "tissue"))



# load cell lines with proteomics
cell_lines_proteomics <- read.table("./files/cell_lines_proteomics_metadata.txt", h=T, sep="\t", stringsAsFactors=F)
nrow(cell_lines_proteomics)
#96



#select cell lines that have proteomics measurements
ccle_cell_lines <- ccle_cell_lines[(ccle_cell_lines$cell %in% cell_lines_proteomics$cell_line), ]
nrow(ccle_cell_lines)
#74



sel_ccle_cell_lines <- c("Name", "Description", as.character(ccle_cell_lines$cell_id))
ccle_counts <- ccle_counts[, sel_ccle_cell_lines]
dim(ccle_counts)
#56318    76



# repeated genes: select the gene that has the highest median across samples
ccle_counts_gene_rep <- ccle_counts[ccle_counts$Description %in% ccle_counts$Description[duplicated(ccle_counts$Description)], ]
dim(ccle_counts_gene_rep)
#2170   76


rmv.genes <- function(mat, index){

    genes <- unique(mat$Description)
    rmv <- c()

    for(gene in genes){
        sub_mat <- mat[mat$Description == gene, c(index:ncol(mat))]
        medians <- apply(sub_mat, 1, function(x) median(x, na.rm=T))
        sub_mat <- sub_mat[-match(max(medians), medians), ]
        rmv <- c(rmv, as.numeric(rownames(sub_mat)))
    }

    return(rmv)
}

rmv <- rmv.genes(ccle_counts_gene_rep, 3)
ccle_counts <- ccle_counts[-rmv, ]
dim(ccle_counts)
#54356    76



# define rownames
rownames(ccle_counts) <- ccle_counts$Description
ccle_counts <- ccle_counts[,-c(1:2)]
dim(ccle_counts)
#54356    74




# simplify colnames
all.equal(colnames(ccle_counts), ccle_cell_lines$cell_id)
#TRUE
colnames(ccle_counts) <- ccle_cell_lines$cell





# remove low expressed genes
# genes with average cpm across samples lower than 1 are removed
ccle_counts <- ccle_counts[ rowMeans( cpm(ccle_counts) ) >= 1, ]
dim(ccle_counts)
#16993    74



# normalise raw counts by TMM method
ccle_counts <- DGEList(counts=ccle_counts)
ccle_counts <- calcNormFactors(ccle_counts, method="TMM")



# derive log-CPM values from voom
ccle_counts <- voom(ccle_counts, design=model.matrix(~tissue, data=ccle_cell_lines))$E



# define gene column
ccle_counts <- cbind(data.frame(gene=rownames(ccle_counts), stringsAsFactors=F), ccle_counts)
rownames(ccle_counts) <- NULL


# export rna data
write.table(ccle_counts, file = "./files/ccle_rna_logCPM.txt", quote = FALSE, sep = "\t", row.names=FALSE)



# gather data
ccle_counts_gather <- ccle_counts %>%
    gather(key="cell_line", value="log2CPM", -gene) %>%
    as.tibble() %>%
    inner_join(cell_lines_proteomics, by="cell_line")


# boxplot by sample
ccle_counts_gather_boxp <- ggplot(ccle_counts_gather, aes(x=cell_line, y=log2CPM, fill=batch)) +
    geom_boxplot(outlier.size = 0.1) +
    theme(axis.title.y=element_text(colour="black", size=13), axis.title.x=element_text(colour="black", size=13)) +
    theme(axis.text.y=element_text(colour="black", size=10), axis.text.x=element_blank()) +
    theme(axis.ticks.x=element_blank()) +
    labs(x = "Cell line", y = "CCLE RNA log2-counts per million (logCPM)", title = "")
ggsave(filename="ccle_rna_log2CPM.png", plot=ccle_counts_gather_boxp, path = "./plots/", width=20)
unlink("ccle_rna_log2CPM.png")





# -- Merge TCGA and cell lines data


# join datasets keeping the common genes
tcga_cellLines_rna <- inner_join(tcga_rna, ccle_counts, by="gene")



# export transcriptomics data
write.table(tcga_cellLines_rna, "./files/rna_tcga_cellLines.txt", sep = "\t", quote = F, row.names=F)


# make dataset with sample batch
sample_batch <- cptac_samples %>%
	mutate(cancer = paste("TCGA", cancer, sep="_")) %>%
	rename(batch=cancer) %>%
	bind_rows(setNames(cell_lines_proteomics[, c("cell_line", "batch")], c("sample", "batch")))



# gather dataset
tcga_cellLines_rna_gather <- tcga_cellLines_rna %>%
	gather(key="sample", value="log2CPM", -gene) %>%
	as.tibble() %>%
	inner_join(sample_batch, by = "sample")


# boxplot by sample
tcga_cellLines_rna_boxp <- ggplot(tcga_cellLines_rna_gather, aes(x=sample, y=log2CPM, fill=batch)) +
    geom_boxplot(outlier.size = 0.1) +
    theme(axis.title.y=element_text(colour="black", size=13), axis.title.x=element_text(colour="black", size=13)) +
    theme(axis.text.y=element_text(colour="black", size=10), axis.text.x=element_blank()) +
    theme(axis.ticks.x=element_blank()) +
    labs(x = "Sample", y = "RNA log2-counts per million (logCPM)", title = "")
ggsave(filename="tcga_cellLines_rna_log2CPM.png", plot=tcga_cellLines_rna_boxp, path = "./plots/", width=20)
unlink("tcga_cellLines_rna_log2CPM.png")











save(list=ls(), file="./r_workspaces/preprocessing_tcga_cell_lines_rnaseq.RData")
