# Abel Sousa 2018

# PhD project
# Study of Control of Protein Abundance using Multi-omics Data from Cancer Samples

# -- CCLE and CPATC datasets

# -- Compute Linear Regressions Between CNVs and Protein Abundances

# -- All Covars are Used in the Regression Models

# -- Linear Models are compared using a Likelihood Ratio Test

# -- Genome-wide protein associations



suppressMessages(library(lmtest))
suppressMessages(library(parallel))
suppressMessages(library(tidyverse))

options(bitmapType = "cairo")



# import protein attenuation data
protein_attenuation <- read_tsv("./files/protein_attenuation.txt")



# make tibble with attenuated - non-attenuated protein pairs
protein_interactions <- tibble(
	attenuated = protein_attenuation %>% filter(class == "high-attenuated-protein" | class == "low-attenuated-protein") %>% pull(gene),
	non_attenuated = list(protein_attenuation %>% filter(class == "non-attenuated") %>% pull(gene))) %>%
	unnest()


protein_interactions <- protein_interactions %>%
  dplyr::select(controlling = non_attenuated, controlled = attenuated) %>%
	mutate_if(is.character, as.factor)


cat("Interactions to calculate:", nrow(protein_interactions), "\n")
#16,106,715



# import common samples with CNV, RNA and Protein abundances
common_samples <- read_tsv("./files/cptac_cellLines_samples_prot_rna_cnv.txt", col_names = FALSE)




# load datasets
proteomics <- read_tsv("./files/proteomicsQ_cptac_cellLines.txt") %>%
  gather(key="sample", value="prot_log2FC", -gene) %>%
  filter(sample %in% common_samples$X1) %>%
  group_by(gene) %>%
  filter( sum(!is.na(prot_log2FC)) > n()*0.25 ) %>%
  mutate(prot_log2FC = scale(prot_log2FC)[,1]) %>%
  ungroup() %>%
	mutate_if(is.character, as.factor) %>%
	arrange(gene, sample) %>%
	dplyr::select(-sample) %>%
	nest(-(gene)) %>%
	dplyr::rename(prot_log2FC = data)


rna <- read_tsv("./files/rna_tcga_cellLines.txt") %>%
  gather(key="sample", value="rna_log2CPM", -gene) %>%
  filter(sample %in% common_samples$X1) %>%
  group_by(gene) %>%
  mutate(rna_log2CPM = scale(rna_log2CPM)[,1]) %>%
  ungroup() %>%
	mutate_if(is.character, as.factor) %>%
	arrange(gene, sample) %>%
	dplyr::select(-sample) %>%
	nest(-(gene)) %>%
	dplyr::rename(rna_log2CPM = data)


cnv <- read_tsv("./files/cnv_tcga_cellLines.txt") %>%
  gather(key="sample", value="cnv_gistic2", -gene) %>%
  filter(sample %in% common_samples$X1) %>%
	mutate_if(is.character, as.factor) %>%
	arrange(gene, sample) %>%
	dplyr::select(-sample) %>%
	nest(-(gene)) %>%
	dplyr::rename(cnv_gistic2 = data)



metadata <- read_tsv("./files/metadata_cptac_cellLines.txt") %>%
  filter(sample %in% common_samples$X1) %>%
	arrange(sample) %>%
  mutate(age = scale(age)[,1]) %>%
	mutate_if(is.character, as.factor)





linear_model <- function(py_prot, py_rna, px_cnv, gend, canc, bat, prot, ag){

  ## define covars matrix
	covars <- bind_cols(py_prot, py_rna, px_cnv, gender = gend, cancer = canc, batch = bat, proteomics = prot, age = ag) %>%
		dplyr::rename(py_prot = prot_log2FC, py_rna = rna_log2CPM, px_cnv = cnv_gistic2)

	covars <- as.data.frame(covars)
  covars <- covars[, c("py_prot", "py_rna", "age", "cancer", "gender", "proteomics", "batch", "px_cnv")]

  ## remove nas
  covars <- na.exclude(covars)

	## remove covariates with less than two levels
  covars <- covars[, apply(covars, 2, function(x) length(levels(as.factor(x))) >= 2), drop=F]


	## fit the models
  #restricted/reduced model (without cnv)
  lreg1_f <- as.formula( paste( "py_prot", "~", paste(colnames(covars[, -match(c("py_prot", "px_cnv"), colnames(covars)), drop=F]), collapse="+") ) )
  lreg1 <- lm(lreg1_f, data = covars)


  #unrestricted model (with cnv)
  lreg2_f <- update.formula(formula(lreg1), ~ . + px_cnv)
  lreg2 <- lm(lreg2_f, data = covars)


  ## summary the models
  suml2 <- summary(lreg2)

  ## coefficients and p-values
  beta_cnv <- suml2$coefficients["px_cnv", "Estimate"]

  ## Likelihood Ratio Test
  ## H0: Restricted model (full or null) is statistically better than Unrestricted (alternative) model
  lrt <- lrtest(lreg1, lreg2)
  lrt_p <- lrt[2,5]


  return( list(c("beta_cnv" = beta_cnv, "lrt_p" = lrt_p)) )

}



# compute linear regressions
model_data <- protein_interactions %>%
	inner_join(proteomics, by=c("controlled" = "gene")) %>%
  dplyr::rename(controlled_protein = prot_log2FC) %>%
  inner_join(rna, by=c("controlled" = "gene")) %>%
  dplyr::rename(controlled_rna = rna_log2CPM) %>%
  inner_join(cnv, by=c("controlling" = "gene")) %>%
  dplyr::rename(controlling_cnv = cnv_gistic2) %>%
  mutate(gender = list(metadata$gender), cancer = list(metadata$cancer), batch = list(metadata$batch), proteomics = list(metadata$proteomics), age = list(metadata$age) ) %>%
	rowwise() %>%
	mutate(res = linear_model(controlled_protein, controlled_rna, controlling_cnv, gender, cancer, batch, proteomics, age)) %>%
	ungroup() %>%
	dplyr::select(controlling, controlled, res)
	#mutate(res = purrr::invoke_map(tibble, res)) %>%
	#mutate(res = map(res, ~ data.frame(t(.)))) %>%
	#unnest(res) %>%
  #mutate(fdr = p.adjust(p=lrt_p, method="BH")) %>%
  #dplyr::arrange(fdr)
saveRDS(model_data, file = "./files/protein_associations_tcga_ccle_cnv_reg_genome_wide.rds")




#cat("Total associations FDR < 0.05:", sum(model_data$fdr<0.05), "\n")





save(list=ls(), file="./r_workspaces/protein_associations_tcga_ccle_cnv_reg_genome_wide.RData")
