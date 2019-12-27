# Abel Sousa 2018

# PhD project
# Study of Control of Protein Abundance using Multi-omics Data from Cancer Samples


library(tidyverse)
library(Biostrings)
options(bitmapType = "cairo")
source("./r_scripts/pre_processing/utils.R")


# -- CPTAC



# -- COREAD

# read mass-spectrometry data
coread <- read.delim("./data/cptac/TCGA_Colon_VU_Proteome_CDAP.r2.precursor_area.tsv", row.names=c(1), na.strings=0)
coread <- coread[, -c(191:194)]


# select only unambigous mapping peptides
coread <- coread[, grep("Unshared.Area", colnames(coread), fixed = TRUE)]
colnames(coread) <- substr(colnames(coread), 1, 14)


# normalize counts by total proteome for each sample
#coread <- apply( coread, 2, function(x) x/sum(x, na.rm=T) )


# calculate fold changes (ratio of protein counts over the median counts across samples)
#coread <- apply( coread, 1, function(x) x/median(x, na.rm=T) )
#coread <- as.data.frame( t(coread), stringsAsFactors=F )


# transform to log2
#coread <- log2(coread)


# match sample ids
coread_clinical <- read.csv("./data/cptac/CPTAC_TCGA_ColorectalCancer_select_clinical_data_release1_090413.csv")
colnames(coread) <- gsub("-", ".", as.character(coread_clinical$TCGA.barcode[match( colnames(coread),  substr(gsub("-", ".", coread_clinical$TCGA.barcode, fixed=T),6,19) )]), fixed=T)


# define gene column names
coread <- cbind(data.frame(gene=rownames(coread), stringsAsFactors=F), coread)
rownames(coread) <- NULL





# -- BRCA

# read mass-spectrometry data
brca <- read.delim("./data/cptac/TCGA_Breast_BI_Proteome_CDAP.r2.itraq.tsv", row.names=c(1))
brca <- brca[-c(1:3), -c(217:220)]


# select only unambigous mapping peptides
brca <- brca[, grep("Unshared.Log.Ratio", colnames(brca), fixed = TRUE)]
colnames(brca) <- sapply( strsplit(colnames(brca), ".Unshared.Log.Ratio", fixed=TRUE), function(x) x[1] )


# average duplicated samples
colnames(brca) <- substr(colnames(brca), 1, 11)
brca <- avg.duplicatedSamp(brca)


# match sample ids
brca_clinical <- read.csv("./data/cptac/CPTAC_TCGA_BreastCancer_select_clinical_data_r1.csv")
brca <- brca[, colnames(brca) %in% substr(gsub("-", ".", brca_clinical$TCGA.barcode, fixed=T),6,16)]
colnames(brca) <- gsub("-", ".", as.character(brca_clinical$TCGA.barcode[match( colnames(brca),  substr(gsub("-", ".", brca_clinical$TCGA.barcode, fixed=T),6,16) )]), fixed=T)


# QC status samples
brca_samples_qc <- read.csv("./data/cptac/CPTAC_TCGA_BreastCancer_select_clinical_data_r1_QC_status.csv")
brca <- brca[, colnames(brca) %in% gsub("-", ".", as.character(brca_samples_qc[(brca_samples_qc$QC.Status == "pass"), "Biospecimen.Barcode.Side"]), fixed=T)]


# define gene column
brca <- cbind(data.frame(gene=rownames(brca), stringsAsFactors=F), brca)
rownames(brca) <- NULL




# -- HGSC

# JHU
hgsc_jhu <- read.delim("./data/cptac/TCGA_Ovarian_JHU_Proteome_CDAP.r2.itraq.tsv", row.names=c(1))
hgsc_jhu <- hgsc_jhu[-c(1:3), -c(265:268)]
colnames(hgsc_jhu) <- sub("X", "", colnames(hgsc_jhu), fixed = TRUE)


# select only unambigous mapping peptides
hgsc_jhu <- hgsc_jhu[, grep("Unshared.Log.Ratio", colnames(hgsc_jhu), fixed = TRUE)]
colnames(hgsc_jhu) <- sapply( strsplit(colnames(hgsc_jhu), ".Unshared.Log.Ratio", fixed=TRUE), function(x) x[1] )


# exclude control samples
hgsc_jhu <- hgsc_jhu[, !grepl("CONTROL", colnames(hgsc_jhu), fixed=T)]



# define gene column
hgsc_jhu <- cbind(data.frame(gene=rownames(hgsc_jhu), stringsAsFactors=F), hgsc_jhu)
rownames(hgsc_jhu) <- NULL



# PNNL
hgsc_pnnl <- read.delim("./data/cptac/TCGA_Ovarian_PNNL_Proteome_CDAP.r2.itraq.tsv", row.names=c(1))
hgsc_pnnl <- hgsc_pnnl[-c(1:3), -c(169:172)]
colnames(hgsc_pnnl) <- sub("X", "", colnames(hgsc_pnnl), fixed = TRUE)


# select only unambigous mapping peptides
hgsc_pnnl <- hgsc_pnnl[, grep("Unshared.Log.Ratio", colnames(hgsc_pnnl), fixed = TRUE)]
colnames(hgsc_pnnl) <- sapply( strsplit(colnames(hgsc_pnnl), ".Unshared.Log.Ratio", fixed=TRUE), function(x) x[1] )


# define gene column
hgsc_pnnl <- cbind(data.frame(gene=rownames(hgsc_pnnl), stringsAsFactors=F), hgsc_pnnl)
rownames(hgsc_pnnl) <- NULL



# merge JHU and PNNL datasets
hgsc <- full_join(hgsc_jhu, hgsc_pnnl, by="gene")


rownames(hgsc) <- hgsc$gene
hgsc <- hgsc[,-c(1)]
colnames(hgsc) <- substr(colnames(hgsc), 1, 11)


# average duplicated samples
hgsc <- avg.duplicatedSamp(hgsc)


# match sample ids
hgsc_clinical <- read.csv("./data/cptac/OV_All_clinical_features_TCGAbiotab_CPTAC_S020.csv")
hgsc <- hgsc[, colnames(hgsc) %in% substr(gsub("-", ".", hgsc_clinical$TCGA.barcode, fixed=T),6,16)]
colnames(hgsc) <- gsub("-", ".", as.character(hgsc_clinical$TCGA.barcode[match( colnames(hgsc),  substr(gsub("-", ".", hgsc_clinical$TCGA.barcode, fixed=T),6,16) )]), fixed=T)


# define gene column
hgsc <- cbind(data.frame(gene=rownames(hgsc), stringsAsFactors=F), hgsc)
rownames(hgsc) <- NULL







# -- Assemble proteomics data-sets into one data set


# join by all genes
proteomics <- full_join(coread, brca, by="gene")
proteomics <- proteomics %>%
	full_join(hgsc, by="gene")


rownames(proteomics) <- proteomics$gene
proteomics <- proteomics[, -c(1)]


# select only primary tumour samples
proteomics <- proteomics[, substr(colnames(proteomics), 14, 15) == "01" ]
colnames(proteomics) <- substr(colnames(proteomics), 1, 12)


# average duplicated (coread) samples
proteomics <- avg.duplicatedSamp(proteomics)




# cptac samples
brca_samples <- data.frame(colnames(proteomics)[colnames(proteomics) %in% substr(colnames(brca), 1, 12)], rep("BRCA", sum(colnames(proteomics) %in% substr(colnames(brca), 1, 12))))
colnames(brca_samples) <- c("sample", "cancer")

hgsc_samples <- data.frame(colnames(proteomics)[colnames(proteomics) %in% substr(colnames(hgsc), 1, 12)], rep("HGSC", sum(colnames(proteomics) %in% substr(colnames(hgsc), 1, 12))))
colnames(hgsc_samples) <- c("sample", "cancer")

coread_samples <- data.frame(colnames(proteomics)[colnames(proteomics) %in% substr(colnames(coread), 1, 12)], rep("COREAD", sum(colnames(proteomics) %in% substr(colnames(coread), 1, 12))))
colnames(coread_samples) <- c("sample", "cancer")


cptac_samples <- rbind(brca_samples, coread_samples, hgsc_samples)
#write.table(cptac_samples, file="./files/cptac_samples.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=T)




# define gene column
proteomics <- cbind(data.frame(gene=rownames(proteomics), stringsAsFactors=F), proteomics)
rownames(proteomics) <- NULL
#write.table(proteomics, file = "./files/cptac_proteomics.txt", quote = FALSE, sep = "\t", row.names=FALSE, col.names=T)




# define sample - cancer code
sample_code <- cptac_samples$cancer
names(sample_code) <- cptac_samples$sample



# gather proteomics dataset
proteomics_gather <- gather(proteomics, key="sample", value="measure", -gene)
proteomics_gather <- proteomics_gather %>%
	mutate(cancer = sample_code[proteomics_gather$sample]) %>%
  mutate(proteomics = if_else(cancer == "COREAD", "LF", "TMT"), data = "CPTAC") %>%
  as.tibble()


cptac_meta <- proteomics_gather %>%
  dplyr::select(-c(gene, measure)) %>%
  distinct()






# -- CCLE cell lines




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
#coread_prot_rmlt <- coread_prot_rmlt/100
#coread_prot_rmlt <- log2(coread_prot_rmlt)



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
#brca_prot_lpk <- apply( brca_prot_lpk, 1, function(x) x/median(x, na.rm=T) )
#brca_prot_lpk <- as.data.frame( t(brca_prot_lpk), stringsAsFactors=F )



# transform to log2
#brca_prot_lpk <- log2(brca_prot_lpk)





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
#brca_prot_law <- apply( brca_prot_law, 1, function(x) x/median(x, na.rm=T) )
#brca_prot_law <- as.data.frame( t(brca_prot_law), stringsAsFactors=F )



# transform to log2
#brca_prot_law <- log2(brca_prot_law)








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


cell_lines_proteomics_gather <- gather(cell_lines_proteomics, key="sample", value="measure", -gene)
cell_lines_proteomics_gather <- cell_lines_proteomics_gather %>%
	mutate(cancer = sample_code[cell_lines_proteomics_gather$sample]) %>%
  mutate(proteomics = if_else(cancer == "law", "LF", "TMT"), data = "CCLE") %>%
  as.tibble()

cell_lines_meta <- cell_lines_proteomics_gather %>%
  dplyr::select(-c(gene, measure)) %>%
  distinct()




# merge CPTAC and cell lines data
#cptac_cell_lines_prot <- bind_rows(proteomics_gather, cell_lines_proteomics_gather)

cptac_cell_lines_meta <- bind_rows(cptac_meta, cell_lines_meta)

cptac_cell_lines_prot <- full_join(
  proteomics_gather %>% dplyr::select(gene, sample, measure) %>% spread(key = "sample", value = "measure"),
  cell_lines_proteomics_gather %>% dplyr::select(gene, sample, measure) %>% spread(key = "sample", value = "measure"),
  by = "gene") %>%
  gather(key="sample", value="measure", -gene) %>%
  inner_join(cptac_cell_lines_meta, by = "sample")



write.table(cptac_cell_lines_prot, file = "./files/cptac_cell_lines_proteomics_raw.txt", quote = FALSE, sep = "\t", row.names=FALSE, col.names=T)
