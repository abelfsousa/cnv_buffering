# Abel Sousa 2018

# PhD project
# Study of Control of Protein Abundance using Multi-omics Data from Cancer Samples



# -- Pre-processing of metadata from CPTAC/TCGA and cell lines


suppressMessages(library(tidyverse))
options(bitmapType = "cairo")




# load TCGA metadata

brca_metadata <- read_tsv("./data/cptac/CPTAC_BreastCancer_metadata.txt") %>%
	mutate(TCGA_ID = gsub("-", ".", TCGA_ID), Proteomics = "TMT", Batch = paste("TCGA", Disease_code, sep="_")) %>%
	rename(cancer = Disease_code, sample = TCGA_ID) %>%
	rename_all(tolower) %>%
	dplyr::select(sample, cancer, batch, proteomics, gender, age)


coread_metadata <- read_tsv("./data/cptac/CPTAC_ColoRectalCancer_metadata.txt") %>%
	mutate(TCGA_ID = substr(gsub("-", ".", TCGA_ID), 1, 12), Proteomics = "LF", Gender = toupper(Gender), Disease_code = "COREAD", Batch = "TCGA_COREAD") %>%
	rename(cancer = Disease_code, sample = TCGA_ID) %>%
	rename_all(tolower) %>%
	dplyr::select(sample, cancer, batch, proteomics, gender, age) %>%
	distinct()

ov_metadata <- read_tsv("./data/cptac/CPTAC_OvarianCancer_metadata.txt") %>%
	mutate(TCGA_ID = substr(gsub("-", ".", TCGA_ID), 1, 12), Proteomics = "TMT", Gender = toupper(Gender), Batch = paste("TCGA", Disease_code, sep="_")) %>%
	rename(cancer = Disease_code, sample = TCGA_ID) %>%
	rename_all(tolower) %>%
	dplyr::select(sample, cancer, batch, proteomics, gender, age) %>%
	distinct()




# load cell lines metadata

cell_lines_metadata <- read_tsv("./files/cell_lines_proteomics_metadata.txt") %>%
	mutate(cancer_type = toupper(cancer_type), batch = toupper(batch), proteomics = toupper(proteomics), gender = "NONE", age = 0) %>%
	rename(sample = cell_line, cancer = cancer_type)





# merge datasets

cptac_cell_lines_metadata <- bind_rows(brca_metadata, coread_metadata, ov_metadata, cell_lines_metadata)


# export

write.table(cptac_cell_lines_metadata, "./files/metadata_cptac_cellLines.txt", sep = "\t", quote = F, row.names=F)





save(list=ls(), file="./r_workspaces/preprocessing_tcga_cell_lines_metadata.RData")







