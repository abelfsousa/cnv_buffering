# Abel Sousa 2018

# PhD project
# Study of Control of Protein Abundance using Multi-omics Data from Cancer Samples



# -- Pre-processing of phosphoproteomics from cptac and cell lines


library(tidyverse)
library(limma)
options(bitmapType = "cairo")



# -- CPTAC data




# import brca phosphoproteomics
brca_phospho <- read_tsv("./data/cptac/TCGA_Breast_BI_Phosphoproteome.phosphosite.itraq.tsv")
brca_qc <- read_csv("./data/cptac/CPTAC_TCGA_BreastCancer_select_clinical_data_r1_QC_status.csv")


#qc samples
brca_qc_samples <- brca_qc %>%
    dplyr::filter(`QC Status` == "pass") %>% pull(`Biospecimen Barcode Side`) %>% substr(6,16)


# change colnames
brca_phospho <- brca_phospho %>%
    setNames(str_split(colnames(brca_phospho), " Log Ratio", simplify=T)[,1])


# add phosphosite
brca_phospho <- brca_phospho %>%
    dplyr::mutate(phospho_site = paste(Gene, str_split(brca_phospho$Phosphosite, ":", simplify=T)[,2], sep="_"))

brca_phospho <- brca_phospho %>%
    dplyr::select(-c(Phosphosite, Peptide, Gene, Organism))

brca_phospho <- brca_phospho %>%
    dplyr::select(phospho_site, everything())


# select samples that passed QC
select_samples <- brca_phospho %>%
    dplyr::select(-c(1)) %>%
    colnames %>%
    substr(1,11) %>%
    intersect(brca_qc_samples)

brca_phospho <- brca_phospho[, c("phospho_site", substr(colnames(brca_phospho)[2:109],1,11)) %in% c("phospho_site", select_samples)]


# average duplicated samples
brca_phospho <- brca_phospho %>%
    dplyr::mutate(`C8-A131-01A` = rowMeans(brca_phospho[, c("C8-A131-01A.1", "C8-A131-01A.2")], na.rm=T)) %>%
    dplyr::mutate(`AO-A12D-01A` = rowMeans(brca_phospho[, c("AO-A12D-01A.1", "AO-A12D-01A.2")], na.rm=T)) %>%
    dplyr::mutate(`AO-A12B-01A` = rowMeans(brca_phospho[, c("AO-A12B-01A.1", "AO-A12B-01A.2")], na.rm=T))

brca_phospho <- brca_phospho %>%
    dplyr::select(-c(`C8-A131-01A.1`, `C8-A131-01A.2`, `AO-A12D-01A.1`, `AO-A12D-01A.2`, `AO-A12B-01A.1`, `AO-A12B-01A.2`))


# median of phosphosites by sample
brca_phospho <- brca_phospho %>%
    gather(key="sample", value="log2FC", -phospho_site) %>%
    group_by(phospho_site, sample) %>%
    summarise(log2FC = median(log2FC, na.rm=T)) %>%
    ungroup() %>%
    spread(key="sample", value="log2FC")


# format colnames
brca_phospho <- brca_phospho %>%
    setNames(c("phospho_site", paste("TCGA", colnames(brca_phospho)[2:77], sep="-")))

brca_phospho <- brca_phospho %>%
    setNames(gsub("-", ".", colnames(brca_phospho)))

brca_phospho <- brca_phospho %>%
    setNames(c("phospho_site", substr(colnames(brca_phospho)[2:77], 1, 12)))




# import hgsc phosphoproteomics
hgsc_phospho <- read_tsv("./data/cptac/TCGA_Ovarian_PNNL_Phosphoproteome.phosphosite.itraq.tsv")


# change colnames
hgsc_phospho <- hgsc_phospho %>%
    setNames(str_split(colnames(hgsc_phospho), " Log Ratio", simplify=T)[,1])


# add phosphosite
hgsc_phospho <- hgsc_phospho %>%
    dplyr::mutate(phospho_site = paste(Gene, str_split(hgsc_phospho$Phosphosite, ":", simplify=T)[,2], sep="_"))

hgsc_phospho <- hgsc_phospho %>%
    dplyr::select(-c(Phosphosite, Peptide, Gene, Organism))

hgsc_phospho <- hgsc_phospho %>%
    dplyr::select(phospho_site, everything())


# median of phosphosites by sample
hgsc_phospho <- hgsc_phospho %>%
    gather(key="sample", value="log2FC", -phospho_site) %>%
    group_by(phospho_site, sample) %>%
    summarise(log2FC = median(log2FC, na.rm=T)) %>%
    ungroup() %>%
    spread(key="sample", value="log2FC")


# format colnames
hgsc_phospho <- hgsc_phospho %>%
    setNames(c("phospho_site", paste("TCGA", colnames(hgsc_phospho)[2:70], sep="-")))

hgsc_phospho <- hgsc_phospho %>%
    setNames(gsub("-", ".", colnames(hgsc_phospho)))

hgsc_phospho <- hgsc_phospho %>%
    setNames(c("phospho_site", substr(colnames(hgsc_phospho)[2:70], 1, 12)))


# merge CPTAC phosphoproteomics
# join all phosphosites
cptac_phospho <- full_join(brca_phospho, hgsc_phospho, by="phospho_site")

# export dataset
write.table(cptac_phospho, "./files/cptac_phospho.txt", sep = "\t", quote = F, row.names=F)




# -- cell lines data
# Phosphoproteome from Roumeliotis et al
# Genomic Determinants of Protein Abundance Variation in Colorectal Cancer Cells

coread_phospho <- read_csv("./data/cell_lines_proteomics/Roumeliotis_etal/roumeliotis_phosphoprot_50_coread_cell_lines.csv")

# format protein site
coread_phospho <- coread_phospho %>%
    mutate(`Protein site` = tolower(gsub("_", "", `Protein site`)) )

# add phosphosite
coread_phospho <- coread_phospho %>%
    mutate(phospho_site = paste(`Gene name`, `Protein site`, sep="_"))

coread_phospho <- coread_phospho %>%
    select(-c(`Gene name`, `Protein site`))

coread_phospho <- coread_phospho %>%
    dplyr::select(phospho_site, everything())

# format colnames
coread_phospho <- coread_phospho %>%
    setNames(c("phospho_site", gsub("-", "", toupper(colnames(coread_phospho)[2:51]))))


# divide each value by 100 (they were multiplied by 100) and transform to log2
coread_phospho <- coread_phospho %>%
    gather(key="sample", value="log2FC", -phospho_site) %>%
    mutate(log2FC = log2(log2FC/100))


# median of phosphosites by sample
coread_phospho <- coread_phospho %>%
    group_by(phospho_site, sample) %>%
    summarise(log2FC = median(log2FC, na.rm=T)) %>%
    ungroup() %>%
    spread(key="sample", value="log2FC")


# export dataset
write.table(coread_phospho, "./files/cell_lines_phospho.txt", sep = "\t", quote = F, row.names=F)



# -- merge CPTAC and cell lines phosphoproteomics


# join all phosphosites
cptac_cellLines_phospho <- full_join(brca_phospho, hgsc_phospho, by="phospho_site") %>%
    full_join(coread_phospho, by="phospho_site")



# export phosphoproteomics data
write.table(cptac_cellLines_phospho, "./files/phosphoproteomics_cptac_cellLines.txt", sep = "\t", quote = F, row.names=F)



# gather data
batch <- rbind(cbind(colnames(brca_phospho)[2:77], "TCGA_BRCA"), cbind(colnames(hgsc_phospho)[2:70], "TCGA_HGSC"), cbind(colnames(coread_phospho)[2:51], "ROUMELIOTIS_COREAD")) %>%
    as.tibble() %>%
    setNames(c("sample", "batch"))

cptac_cellLines_phospho_gather <- cptac_cellLines_phospho %>%
    gather(key="sample", value="log2FC", -phospho_site) %>%
    inner_join(batch, by="sample")


# boxplot by sample
cptac_cellLines_phospho_boxp <- ggplot(cptac_cellLines_phospho_gather, aes(x=sample, y=log2FC, fill=batch)) +
    geom_boxplot() +
    theme(axis.title.y=element_text(colour="black", size=13), axis.title.x=element_text(colour="black", size=13)) +
    theme(axis.text.y=element_text(colour="black", size=10), axis.text.x=element_blank()) +
    theme(axis.ticks.x=element_blank()) +
    labs(x = "Sample", y = "log2FC phosphoproteomics", title = "")
ggsave(filename="cptac_cellLines_phosphoproteomics.png", plot=cptac_cellLines_phospho_boxp, path = "./plots/", width=20)
unlink("cptac_cellLines_phosphoproteomics.png")




# normalize samples using quantile normalization
phosphosites <- cptac_cellLines_phospho$phospho_site
cptac_cellLines_phosphoQ <- normalizeQuantiles( as.matrix(cptac_cellLines_phospho[, -c(1)]) )


# define phosphosite column
cptac_cellLines_phosphoQ <- cbind(data.frame(phospho_site=phosphosites, stringsAsFactors=F), cptac_cellLines_phosphoQ)
rownames(cptac_cellLines_phosphoQ) <- NULL
write.table(cptac_cellLines_phosphoQ, file = "./files/phosphoproteomicsQ_cptac_cellLines.txt", quote = FALSE, sep = "\t", row.names=FALSE)



# gather data
cptac_cellLines_phosphoQ_gather <- cptac_cellLines_phosphoQ %>%
    gather(key="sample", value="log2FC", -phospho_site) %>%
    as.tibble() %>%
    inner_join(batch, by="sample")


# boxplot by sample
cptac_cellLines_phosphoQ_boxp <- ggplot(cptac_cellLines_phosphoQ_gather, aes(x=sample, y=log2FC, fill=batch)) +
    geom_boxplot() +
    theme(axis.title.y=element_text(colour="black", size=13), axis.title.x=element_text(colour="black", size=13)) +
    theme(axis.text.y=element_text(colour="black", size=10), axis.text.x=element_blank()) +
    theme(axis.ticks.x=element_blank()) +
    labs(x = "Sample", y = "log2FC phosphoproteomics", title = "")
ggsave(filename="cptac_cellLines_phosphoproteomicsQ.png", plot=cptac_cellLines_phosphoQ_boxp, path = "./plots/", width=20)
unlink("cptac_cellLines_phosphoproteomicsQ.png")





save(list=ls(), file="./r_workspaces/preprocessing_cptac_cell_lines_phosphoproteomics.RData")
