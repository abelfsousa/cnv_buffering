# Abel Sousa 2018

# PhD project
# Study of Control of Protein Abundance using Multi-omics Data from Cancer Samples


# -- Pre-processing of cell lines data


# useful cell lines database:
# https://web.expasy.org/cellosaurus/



library(tidyverse)
library(limma)
library(Biostrings)
options(bitmapType = "cairo")

source("./r_scripts/pre_processing/utils.R")


# -- Proteomics



# Roumeliotis et al
# Genomic Determinants of Protein Abundance Variation in Colorectal Cancer Cells
# https://www.sciencedirect.com/science/article/pii/S2211124717311002?via%3Dihub#app3

# Protein quantification values for 50 colorectal cancer cell lines.

coread_prot_rmlt <- read.csv("./data/cell_lines_proteomics/Roumeliotis_etal/roumeliotis_prot_50_coread_cell_lines.csv", stringsAsFactors=F, check.names=T, na.strings=c(""))
dim(coread_prot_rmlt)
#9410   52



# remove NA genes
coread_prot_rmlt <- coread_prot_rmlt[!is.na(coread_prot_rmlt$Gene.name), ]
rownames(coread_prot_rmlt) <- NULL

dim(coread_prot_rmlt)
#9409   52




# select genes with multiple protein isoforms
coread_prot_rmlt_dup <- coread_prot_rmlt[(coread_prot_rmlt$Gene.name %in% coread_prot_rmlt$Gene.name[duplicated(coread_prot_rmlt$Gene.name)]), ]



# select protein isoform for the same gene with highest median across samples
rmv.proteins <- function(mat, index){

    genes <- unique(mat$Gene.name)
    rmv <- c()

    for(gene in genes){
        sub_mat <- mat[mat$Gene.name == gene, c(index:ncol(mat))]
        medians <- apply(sub_mat, 1, function(x) median(x, na.rm=T))
        sub_mat <- sub_mat[-match(max(medians), medians), ]
        rmv <- c(rmv, as.numeric(rownames(sub_mat)))
    }

    return(rmv)
}

rmv <- rmv.proteins(coread_prot_rmlt_dup, 3)
coread_prot_rmlt <- coread_prot_rmlt[-rmv, ]



# define row names
rownames(coread_prot_rmlt) <- coread_prot_rmlt$Gene.name



# remove gene symbol and uniprot IDs columns
coread_prot_rmlt <- coread_prot_rmlt[, c(3:52)]
dim(coread_prot_rmlt)
#9395   50



# format cell line IDs
colnames(coread_prot_rmlt) <- gsub(".", "", toupper(colnames(coread_prot_rmlt)), fixed=T)



# divide by 100 and transform to log2 (values were multiplied by 100 (see methods))
coread_prot_rmlt <- coread_prot_rmlt/100
coread_prot_rmlt <- log2(coread_prot_rmlt)






# Lapek et al
# Detection of dysregulated protein-association networks by high-throughput proteomics predicts cancer vulnerabilities
# https://www.nature.com/articles/nbt.3955#main

# Protein quantification values for 41 breast cancer cell lines.


brca_prot_lpk <- read.csv("./data/cell_lines_proteomics/Lapek_etal/lapek_prot_41_brca_cell_lines.csv", stringsAsFactors=F, check.names=T, na.strings=c("n/a"))
dim(brca_prot_lpk)
#10535    43


# simplify protein names
brca_prot_lpk$Protein.name <- sapply(strsplit(brca_prot_lpk$Protein.name, "|", fixed=T), function(x) x[2])


# uniprot protein to gene map
uniprot <- uniprot2genename("./data/uniprot/uniprot_sprot_2017_11_23.fasta")
uniprot <- do.call(rbind, uniprot)


# keep only proteins that map to a gene name
brca_prot_lpk <- brca_prot_lpk[brca_prot_lpk$Protein.name %in% rownames(uniprot), ]
dim(brca_prot_lpk)
#9870   43


# reset row names
rownames(brca_prot_lpk) <- NULL


# obtain gene names
gene_names <- unname(uniprot[brca_prot_lpk$Protein.name, c(1)])


# bind gene names
brca_prot_lpk <- setNames( cbind(gene_names, brca_prot_lpk), c("Gene.name", colnames(brca_prot_lpk)) )


# select genes with multiple protein isoforms
brca_prot_lpk_dup <- brca_prot_lpk[(brca_prot_lpk$Gene.name %in% brca_prot_lpk$Gene.name[duplicated(brca_prot_lpk$Gene.name)]), ]


# select protein isoform for the same gene with highest median across samples
rmv <- rmv.proteins(brca_prot_lpk_dup, 4)
brca_prot_lpk <- brca_prot_lpk[-rmv, ]


# define row names
rownames(brca_prot_lpk) <- brca_prot_lpk$Gene.name


# remove gene symbol and uniprot IDs columns
brca_prot_lpk <- brca_prot_lpk[, c(4:44)]
dim(brca_prot_lpk)
#9839   41



# calculate fold changes (ratio of protein intensities over the median intensities across samples)
brca_prot_lpk <- apply( brca_prot_lpk, 1, function(x) x/median(x, na.rm=T) )
brca_prot_lpk <- as.data.frame( t(brca_prot_lpk), stringsAsFactors=F )



# transform to log2
brca_prot_lpk <- log2(brca_prot_lpk)





# Lawrence et al
# The Proteomic Landscape of Triple-Negative Breast Cancer
# http://www.cell.com/cell-reports/abstract/S2211-1247(15)00341-1

# Label-free deep proteome analysis of 24 human breast specimens.

# Values represent the summed peptide intensity normalized by total proteome intensity for that sample and by the number of theoretically observable peptides for the protein.
# For each sample (including all fractions), summed peak areas for all peptides matching to the same protein were divided by the maximum number of observable peptides, and normalized to the total sum intensity of the observed proteome.


brca_prot_law <- read.csv("./data/cell_lines_proteomics/Lawrence_etal/lawrence_prot_20_brca_cell_lines.csv", stringsAsFactors=F, check.names=T, na.strings=c(""))
dim(brca_prot_law)
#15524    45


# protein names
protein.names <- sapply(strsplit(brca_prot_law$Uniprot.Entry, "|", fixed=T), function(x) x[2])

dim(uniprot[protein.names[protein.names %in% rownames(uniprot)], ])
#7964    2



# gene names
gene.names <- sapply(brca_prot_law$Uniprot.Entry, function(x) {
    v <- strsplit(x, " ", fixed=T)[[1]];
    if(sum(grepl("GN=", v, fixed=T)) == 0){ NA }
    else{
        i <- grep("GN=", v, fixed=T);
        strsplit(v[i], "=", fixed=T)[[1]][2];
        }
    }, USE.NAMES=F)


brca_prot_law <- setNames(cbind(gene.names, protein.names, brca_prot_law), c("Gene.name", "Protein.name", colnames(brca_prot_law)))

dim(brca_prot_law)
#15524    47


# remove NA genes
brca_prot_law <- brca_prot_law[!is.na(brca_prot_law$Gene.name), ]
rownames(brca_prot_law) <- NULL


# remove tumour samples
brca_prot_law <- brca_prot_law[,-c(44, 45, 46, 47)]


# remove Uniprot.Entry column
brca_prot_law <- brca_prot_law[,-c(3)]


dim(brca_prot_law)
#15404    42



# average sample replicates
colnames(brca_prot_law)[3:42] <- sapply(strsplit(colnames(brca_prot_law)[3:42], ".", fixed=T), function(x) x[1])
brca_prot_law <- avg.duplicatedSamp(brca_prot_law)




# select genes with multiple protein isoforms
brca_prot_law_dup <- brca_prot_law[(brca_prot_law$Gene.name %in% brca_prot_law$Gene.name[duplicated(brca_prot_law$Gene.name)]), ]




# select protein isoform for the same gene with highest median across samples
rmv <- rmv.proteins(brca_prot_law_dup, 3)
brca_prot_law <- brca_prot_law[-rmv, ]

dim(brca_prot_law)
#11362    22


# define row names
rownames(brca_prot_law) <- brca_prot_law$Gene.name


# remove gene symbol and protein IDs columns
brca_prot_law <- brca_prot_law[, c(3:22)]
dim(brca_prot_law)
#11362    20



# calculate fold changes (ratio of protein intensities over the median intensities across samples)
brca_prot_law <- apply( brca_prot_law, 1, function(x) x/median(x, na.rm=T) )
brca_prot_law <- as.data.frame( t(brca_prot_law), stringsAsFactors=F )



# transform to log2
brca_prot_law <- log2(brca_prot_law)








# cell lines map
cell_map <- c( colnames(coread_prot_rmlt), colnames(brca_prot_lpk), colnames(brca_prot_law) )
cell_map <- cbind( cell_map, c(rep("coread", ncol(coread_prot_rmlt)), rep("brca", ncol(brca_prot_lpk)+ncol(brca_prot_law))) )
cell_map <- cbind( cell_map, c(rep("rmlt", ncol(coread_prot_rmlt)), rep("lpk", ncol(brca_prot_lpk)), rep("law", ncol(brca_prot_law))) )
cell_map <- cbind( cell_map, c(rep("tmt", ncol(coread_prot_rmlt)), rep("tmt", ncol(brca_prot_lpk)), rep("lf", ncol(brca_prot_law))) )

cell_map <- as.data.frame(cell_map)
colnames(cell_map) <- c("cell_line", "cancer_type", "batch", "proteomics")



# remove repeated cell lines between experiments
# for the cell lines duplicated between datasets I will keep those on tmt datasets

# common cell lines between Lapek et al and Lawrence et al
common_lpk_law <- intersect(cell_map[cell_map$batch == "lpk", "cell_line"], cell_map[cell_map$batch == "law", "cell_line"])


# cell lines to remove from Lawrence et al dataset
rmv <- rownames(cell_map[(cell_map$batch == "law" & cell_map$cell_line %in% common_lpk_law), ])



cell_map <- cell_map[-as.numeric(rmv),]
rownames(cell_map) <- NULL



# remove repetead cell lines and add column with gene name
coread_prot_rmlt <- cbind(data.frame(gene = rownames(coread_prot_rmlt), stringsAsFactors=F), coread_prot_rmlt)
brca_prot_lpk <- cbind(data.frame(gene = rownames(brca_prot_lpk), stringsAsFactors=F), brca_prot_lpk)

brca_prot_law <- brca_prot_law[, as.character(cell_map[cell_map$batch == "law", "cell_line"])]
brca_prot_law <- cbind(data.frame(gene = rownames(brca_prot_law), stringsAsFactors=F), brca_prot_law)




# join all datasets

cell_lines_proteomics <- full_join(coread_prot_rmlt, brca_prot_lpk, by="gene")
cell_lines_proteomics <- cell_lines_proteomics %>%
    full_join(brca_prot_law, by="gene")




# load ccle cell lines
ccle_cells <- read.table("./data/ccle/ccle_cell_lines_ids.txt", h=F, stringsAsFactors=F)

length(intersect(colnames(cell_lines_proteomics), ccle_cells$V1))
#75

length(intersect(colnames(coread_prot_rmlt), ccle_cells$V1))
#41

length(intersect(colnames(brca_prot_lpk), ccle_cells$V1))
#31

length(intersect(colnames(brca_prot_law), ccle_cells$V1))
#3





sample_code <- cell_map$batch
names(sample_code) <- cell_map$cell_line


cell_lines_proteomics_gather <- gather(cell_lines_proteomics, key="sample", value="fc", -gene)
cell_lines_proteomics_gather <- cell_lines_proteomics_gather %>%
	mutate(batch = sample_code[cell_lines_proteomics_gather$sample])

cell_lines_proteomics_distribution_1 <- ggplot(cell_lines_proteomics_gather, aes(sample, fc, fill=batch)) +
	theme_classic() +
	geom_boxplot(outlier.size = 0.1) +
	theme(axis.title=element_text(colour="black", size=12),
		  axis.text.y=element_text(colour="black", size=10),
		  axis.text.x=element_blank(),
		  axis.ticks.x=element_blank()) +
	labs(x = "Sample", y = "log2 FC")
ggsave(filename="cell_lines_proteomics_distribution_1.png", plot = cell_lines_proteomics_distribution_1, path = "./plots/", width = 20)
unlink("cell_lines_proteomics_distribution_1.png")





# normalize samples using quantile normalization
genes <- cell_lines_proteomics$gene
cell_lines_proteomicsQ <- normalizeQuantiles( as.matrix(cell_lines_proteomics[, -c(1)]) )


# define gene column
cell_lines_proteomicsQ <- cbind(data.frame(gene=genes, stringsAsFactors=F), cell_lines_proteomicsQ)
rownames(cell_lines_proteomicsQ) <- NULL



cell_lines_proteomicsQ_gather <- gather(cell_lines_proteomicsQ, key="sample", value="fc", -gene)
cell_lines_proteomicsQ_gather <- cell_lines_proteomicsQ_gather %>%
	mutate(batch = sample_code[cell_lines_proteomicsQ_gather$sample])

cell_lines_proteomics_distribution_2 <- ggplot(cell_lines_proteomicsQ_gather, aes(sample, fc, fill=batch)) +
	theme_classic() +
	geom_boxplot(outlier.size = 0.1) +
	theme(axis.title=element_text(colour="black", size=12),
		  axis.text.y=element_text(colour="black", size=10),
		  axis.text.x=element_blank(),
		  axis.ticks.x=element_blank()) +
	labs(x = "Sample", y = "log2 FC")
ggsave(filename="cell_lines_proteomics_distribution_2.png", plot = cell_lines_proteomics_distribution_2, path = "./plots/", width = 20)
unlink("cell_lines_proteomics_distribution_2.png")







# export datasets
write.table(cell_map, "./files/cell_lines_proteomics_metadata.txt", sep="\t", row.names=F, quote=F)
write.table(cell_lines_proteomics, "./files/cell_lines_proteomics.txt", sep="\t", row.names=F, quote=F)
write.table(cell_lines_proteomicsQ, "./files/cell_lines_proteomicsQ.txt", sep="\t", row.names=F, quote=F)




save(list=ls(), file="./r_workspaces/preprocessing_cell_lines_proteomics.RData")
