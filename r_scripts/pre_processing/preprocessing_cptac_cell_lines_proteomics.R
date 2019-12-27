# Abel Sousa 2018

# PhD project
# Study of Control of Protein Abundance using Multi-omics Data from Cancer Samples



# -- Pre-processing of proteomics from cptac and cell lines


library(tidyverse)
library(limma)
options(bitmapType = "cairo")



# load proteomics datasets
cptac_proteomics <- read_tsv("./files/cptac_proteomics.txt")
cptac_samples <- read_tsv("./files/cptac_samples.txt") %>%
	mutate(cancer = paste("CPTAC", cancer, sep="_"))

cell_lines_proteomics <- read_tsv("./files/cell_lines_proteomics.txt")
cell_lines_samples <- read_tsv("./files/cell_lines_proteomics_metadata.txt") %>%
	mutate(batch = paste("cell_lines", batch, sep="_"))



# batch code
cptac_code <- cptac_samples$cancer %>% setNames(cptac_samples$sample)
cell_lines_code <- cell_lines_samples$batch %>% setNames(cell_lines_samples$cell_line)
samples_code <- c(cptac_code, cell_lines_code)

batch <- full_join(cptac_samples, cell_lines_samples, by=c("sample" = "cell_line", "cancer" = "batch")) %>%
	dplyr::select(-c(cancer_type, proteomics))




# join datasets keeping all genes
cptac_cellLines_proteomics <- full_join(cptac_proteomics, cell_lines_proteomics, by="gene")



# export proteomics data
write.table(cptac_cellLines_proteomics, "./files/proteomics_cptac_cellLines.txt", sep = "\t", quote = F, row.names=F)


# gather dataset
cptac_cellLines_proteomics_gather <- cptac_cellLines_proteomics %>%
	gather(key="sample", value="log2FC", -gene) %>%
	inner_join(batch, by = "sample") %>%
    mutate_if(is.character, as.factor) %>%
    mutate(sample = fct_relevel(sample, batch$sample)) %>%
    mutate(cancer = fct_relevel(cancer, "CPTAC_BRCA", "CPTAC_COREAD", "CPTAC_HGSC", "cell_lines_rmlt", "cell_lines_lpk", "cell_lines_law"))


# boxplot by sample
cptac_cellLines_proteomics_boxp <- ggplot(cptac_cellLines_proteomics_gather, aes(x=sample, y=log2FC, fill=cancer)) +
    geom_boxplot(outlier.size = 0.1) +
    theme(axis.title.y=element_text(colour="black", size=13), axis.title.x=element_text(colour="black", size=13)) +
    theme(axis.text.y=element_text(colour="black", size=10), axis.text.x=element_blank()) +
    theme(axis.ticks.x=element_blank()) +
    labs(x = "Sample", y = "log2FC proteomics", title = "")
ggsave(filename="cptac_cellLines_proteomics.png", plot=cptac_cellLines_proteomics_boxp, path = "./plots/", width=20)
unlink("cptac_cellLines_proteomics.png")






# normalize samples using quantile normalization
genes <- cptac_cellLines_proteomics$gene
cptac_cellLines_proteomicsQ <- normalizeQuantiles( as.matrix(cptac_cellLines_proteomics[, -c(1)]) )


# define gene column
cptac_cellLines_proteomicsQ <- cbind(data.frame(gene=genes, stringsAsFactors=F), cptac_cellLines_proteomicsQ)
rownames(cptac_cellLines_proteomicsQ) <- NULL


# export proteomics data
write.table(cptac_cellLines_proteomicsQ, file = "./files/proteomicsQ_cptac_cellLines.txt", quote = FALSE, sep = "\t", row.names=FALSE)



# gather data
cptac_cellLines_proteomicsQ_gather <- cptac_cellLines_proteomicsQ %>%
    gather(key="sample", value="log2FC", -gene) %>%
    as.tibble() %>%
    inner_join(batch, by="sample") %>%
    mutate_if(is.character, as.factor) %>%
    mutate(sample = fct_relevel(sample, batch$sample)) %>%
    mutate(cancer = fct_relevel(cancer, "CPTAC_BRCA", "CPTAC_COREAD", "CPTAC_HGSC", "cell_lines_rmlt", "cell_lines_lpk", "cell_lines_law"))


# boxplot by sample
cptac_cellLines_proteomicsQ_boxp <- ggplot(cptac_cellLines_proteomicsQ_gather, aes(x=sample, y=log2FC, fill=cancer)) +
    geom_boxplot(outlier.size = 0.1) +
    theme(axis.title.y=element_text(colour="black", size=13), axis.title.x=element_text(colour="black", size=13)) +
    theme(axis.text.y=element_text(colour="black", size=10), axis.text.x=element_blank()) +
    theme(axis.ticks.x=element_blank()) +
    labs(x = "Sample", y = "log2FC proteomics", title = "")
ggsave(filename="cptac_cellLines_proteomicsQ.png", plot=cptac_cellLines_proteomicsQ_boxp, path = "./plots/", width=20)
unlink("cptac_cellLines_proteomicsQ.png")





save(list=ls(), file="./r_workspaces/preprocessing_cptac_cell_lines_proteomics.RData")









