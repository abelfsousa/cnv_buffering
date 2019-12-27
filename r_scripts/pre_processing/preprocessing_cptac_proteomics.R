# Abel Sousa 2018

# PhD project
# Study of Control of Protein Abundance using Multi-omics Data from Cancer Samples


# -- Pre-processing of CPTAC proteomics data

# -- Proteomic data downloaded from CPTAC Data Portal




library(tidyverse)
library(limma)
options(bitmapType = "cairo")

source("./r_scripts/pre_processing/utils.R")




# -- COREAD

# read mass-spectrometry data
coread <- read.delim("./data/cptac/TCGA_Colon_VU_Proteome_CDAP.r2.precursor_area.tsv", row.names=c(1), na.strings=0)
coread <- coread[, -c(191:194)]


# select only unambigous mapping peptides
coread <- coread[, grep("Unshared.Area", colnames(coread), fixed = TRUE)]
colnames(coread) <- substr(colnames(coread), 1, 14)


# normalize counts by total proteome for each sample
coread <- apply( coread, 2, function(x) x/sum(x, na.rm=T) )


# calculate fold changes (ratio of protein counts over the median counts across samples)
coread <- apply( coread, 1, function(x) x/median(x, na.rm=T) )
coread <- as.data.frame( t(coread), stringsAsFactors=F )


# transform to log2
coread <- log2(coread)


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
write.table(cptac_samples, file="./files/cptac_samples.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=T)




# define gene column
proteomics <- cbind(data.frame(gene=rownames(proteomics), stringsAsFactors=F), proteomics)
rownames(proteomics) <- NULL
write.table(proteomics, file = "./files/cptac_proteomics.txt", quote = FALSE, sep = "\t", row.names=FALSE, col.names=T)




# define sample - cancer code
sample_code <- cptac_samples$cancer
names(sample_code) <- cptac_samples$sample



# gather proteomics dataset
proteomics_gather <- gather(proteomics, key="sample", value="fc", -gene)
proteomics_gather <- proteomics_gather %>%
	mutate(cancer = sample_code[proteomics_gather$sample])

cptac_proteomics_distribution <- ggplot(proteomics_gather, aes(sample, fc, fill=cancer)) +
	theme_classic() +
	geom_boxplot(outlier.size = 0.1) +
	theme(axis.title=element_text(colour="black", size=12),
		  axis.text.y=element_text(colour="black", size=10),
		  axis.text.x=element_blank(),
		  axis.ticks.x=element_blank()) +
	labs(x = "Sample", y = "log2 FC")
ggsave(filename="cptac_proteomics_distribution.png", plot = cptac_proteomics_distribution, path = "./plots/", width = 20)
unlink("cptac_proteomics_distribution.png")





# normalize samples using quantile normalization
genes <- proteomics$gene
proteomicsQ <- normalizeQuantiles( as.matrix(proteomics[, -c(1)]) )


# define gene column
proteomicsQ <- cbind(data.frame(gene=genes, stringsAsFactors=F), proteomicsQ)
rownames(proteomicsQ) <- NULL
write.table(proteomicsQ, file = "./files/cptac_proteomicsQ.txt", quote = FALSE, sep = "\t", row.names=FALSE, col.names=T)



proteomicsQ_gather <- gather(proteomicsQ, key="sample", value="fc", -gene)
proteomicsQ_gather <- proteomicsQ_gather %>%
	mutate(cancer = sample_code[proteomicsQ_gather$sample])

cptac_proteomicsQ_distribution <- ggplot(proteomicsQ_gather, aes(sample, fc, fill=cancer)) +
	theme_classic() +
	geom_boxplot(outlier.size = 0.1) +
	theme(axis.title=element_text(colour="black", size=12),
		  axis.text.y=element_text(colour="black", size=10),
		  axis.text.x=element_blank(),
		  axis.ticks.x=element_blank()) +
	labs(x = "Sample", y = "log2 FC")
ggsave(filename="cptac_proteomicsQ_distribution.png", plot = cptac_proteomicsQ_distribution, path = "./plots/", width = 20)
unlink("cptac_proteomicsQ_distribution.png")







# remove outlier samples of COREAD cancer (samples whose median is 1.5x the standard deviation of all samples)
coread_outlier_samples <- proteomics_gather %>%
	filter(cancer == "COREAD") %>%
	group_by(sample) %>%
	summarise(median_fc = median(fc, na.rm=T)) %>%
	mutate(sd_all = sd(median_fc)) %>%
	mutate(remove = abs(median_fc) > sd_all*1.5) %>%
	filter(remove) %>%
	pull(sample)


proteomics2 <- proteomics %>%
	dplyr::select(-coread_outlier_samples)
write.table(proteomics2, file = "./files/cptac_proteomics2.txt", quote = FALSE, sep = "\t", row.names=FALSE, col.names=T)


# gather proteomics dataset
proteomics2_gather <- gather(proteomics2, key="sample", value="fc", -gene)
proteomics2_gather <- proteomics2_gather %>%
	mutate(cancer = sample_code[proteomics2_gather$sample])

cptac_proteomics2_distribution <- ggplot(proteomics2_gather, aes(sample, fc, fill=cancer)) +
	theme_classic() +
	geom_boxplot(outlier.size = 0.1) +
	theme(axis.title=element_text(colour="black", size=12),
		  axis.text.y=element_text(colour="black", size=10),
		  axis.text.x=element_blank(),
		  axis.ticks.x=element_blank()) +
	labs(x = "Sample", y = "log2 FC")
ggsave(filename="cptac_proteomics2_distribution.png", plot = cptac_proteomics2_distribution, path = "./plots/", width = 20)
unlink("cptac_proteomics2_distribution.png")




save(list=ls(), file="./r_workspaces/preprocessing_cptac_proteomics.RData")
