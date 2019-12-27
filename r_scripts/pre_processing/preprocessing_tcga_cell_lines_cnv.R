# Abel Sousa 2018

# PhD project
# Study of Control of Protein Abundance using Multi-omics Data from Cancer Samples


# -- Pre-processing of TCGA and cancer cell lines CNV data


# -- TCGA CNV data downloaded from firebrowse


# -- Cell lines CNV data download from...


# -- Useful cell lines database: https://web.expasy.org/cellosaurus/


suppressMessages(library(tidyverse))
options(bitmapType = "cairo")






# -- TCGA CNV data


# load TCGA samples with proteomics
cptac_samples <- read_tsv("./files/cptac_samples.txt")


# import CNV datasets
brca_cnv <- read_tsv("./data/tcga/cnv_firebrowse/BRCA/all_thresholded.by_genes.txt") %>%
    rename(gene=`Gene Symbol`) %>%
    dplyr::select(-c(`Locus ID`, Cytoband))


coread_cnv <- read_tsv("./data/tcga/cnv_firebrowse/COREAD/all_thresholded.by_genes.txt") %>%
    rename(gene=`Gene Symbol`) %>%
    dplyr::select(-c(`Locus ID`, Cytoband))


hgsc_cnv <- read_tsv("./data/tcga/cnv_firebrowse/HGSC/all_thresholded.by_genes.txt") %>%
    rename(gene=`Gene Symbol`) %>%
    dplyr::select(-c(`Locus ID`, Cytoband))



# merge data

tcga_cnv <- full_join(brca_cnv, coread_cnv, by="gene") %>%
    full_join(hgsc_cnv, by="gene")


# select only samples with proteomics measurements
tcga_cnv <- tcga_cnv %>%
    setNames(substr(gsub("-", ".", colnames(tcga_cnv)), 1, 12))

tcga_cnv <- tcga_cnv %>%
    dplyr::select(c("gene", intersect(colnames(tcga_cnv), cptac_samples$sample)))



#export dataset
write.table(tcga_cnv, "./files/tcga_cnv.txt", sep = "\t", quote = F, row.names=F)






# -- cell lines CNV data


# load cell lines with proteomics measurements
cell_lines_proteomics <- read_tsv("./files/cell_lines_proteomics_metadata.txt")
nrow(cell_lines_proteomics)
#96


# load CNV data
ccle_cnv <- read_tsv("./data/ccle/data_CNA.txt")
dim(ccle_cnv)
#23082   998


# define gene column
ccle_cnv <- ccle_cnv %>%
    rename(gene = Hugo_Symbol) %>%
    dplyr::select(-c(Entrez_Gene_Id, Cytoband))




# cell / tissue map
cells <- colnames(ccle_cnv)[2:996]
tissues <- sub("_", "|", cells, fixed=T)
tissues <- sapply(strsplit(tissues, "|", fixed=T), function(x) x[2])
ccle_cell_lines <- cbind(cells, toupper(str_extract(cells, "^[^_]+")), tissues)
ccle_cell_lines <- setNames(as.data.frame(ccle_cell_lines, stringsAsFactors=F), c("cell_id", "cell", "tissue"))
nrow(ccle_cell_lines)
#995




#select cell lines that overlap with proteomics dataset
ccle_cell_lines <- ccle_cell_lines[(ccle_cell_lines$cell %in% cell_lines_proteomics$cell_line), ]
nrow(ccle_cell_lines)
#75

ccle_cnv <- ccle_cnv %>%
    dplyr::select(c("gene", ccle_cell_lines$cell_id))



# simplify colnames
all.equal(colnames(ccle_cnv)[2:76], ccle_cell_lines$cell_id)
#TRUE

ccle_cnv <- ccle_cnv %>%
    setNames( c("gene", ccle_cell_lines$cell) )


dim(ccle_cnv)
#23082    76


# export dataset
write.table(ccle_cnv, "./files/ccle_cnv.txt", sep = "\t", quote = F, row.names=F)





# -- Merge TCGA and cell lines data


# join datasets keeping the common genes
tcga_cellLines_cnv <- inner_join(tcga_cnv, ccle_cnv, by="gene")



# export proteomics data
write.table(tcga_cellLines_cnv, "./files/cnv_tcga_cellLines.txt", sep = "\t", quote = F, row.names=F)










save(list=ls(), file="./r_workspaces/preprocessing_tcga_cell_lines_cnv.RData")
