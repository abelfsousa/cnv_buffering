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



cat("Interactions to calculate:", nrow(protein_interactions), "\n")
#16,106,715



# import common samples with CNV, RNA and Protein abundances
common_samples <- read_tsv("./files/cptac_cellLines_samples_prot_rna_cnv.txt", col_names = FALSE)




# load datasets
proteomics <- read_tsv("./files/proteomicsQ_cptac_cellLines.txt") %>%
    gather(key="sample", value="prot_log2FC", -gene) %>%
		arrange(sample) %>%
    filter(sample %in% common_samples$X1) %>%
    group_by(gene) %>%
    filter( sum(!is.na(prot_log2FC)) > n()*0.25 ) %>%
    mutate(prot_log2FC = scale(prot_log2FC)[,1]) %>%
    ungroup() %>%
    spread(key="sample", value="prot_log2FC") %>%
    as.data.frame()
rownames(proteomics) <- proteomics$gene
proteomics <- proteomics[, -c(1)]


rna <- read_tsv("./files/rna_tcga_cellLines.txt") %>%
    gather(key="sample", value="rna_log2CPM", -gene) %>%
		arrange(sample) %>%
    filter(sample %in% common_samples$X1) %>%
    group_by(gene) %>%
    mutate(rna_log2CPM = scale(rna_log2CPM)[,1]) %>%
    ungroup() %>%
    spread(key="sample", value="rna_log2CPM") %>%
    as.data.frame()
rownames(rna) <- rna$gene
rna <- rna[, -c(1)]


cnv <- read_tsv("./files/cnv_tcga_cellLines.txt") %>%
    gather(key="sample", value="cnv_gistic2", -gene) %>%
		arrange(sample) %>%
    filter(sample %in% common_samples$X1) %>%
    spread(key="sample", value="cnv_gistic2") %>%
    as.data.frame()
rownames(cnv) <- cnv$gene
cnv <- cnv[, -c(1)]



metadata <- read_tsv("./files/metadata_cptac_cellLines.txt") %>%
		arrange(sample) %>%
    filter(sample %in% common_samples$X1) %>%
    mutate(age = scale(age)[,1]) %>%
    as.data.frame()
rownames(metadata) <- metadata$sample
metadata <- metadata[, -c(1)]







linear_model <- function(px, py, proteomics, rna, cnv, metadata){
    #print(paste(px, py, sep=" "))

	py_prot <- t(proteomics[py, ])
	py_rna <- t(rna[py, ])
	px_cnv <- t(cnv[px, ])


  ## define covars matrix
	## bind prot/rna of Py and cnv of Px
	covars <- setNames(cbind(metadata, py_prot, py_rna, px_cnv), c(colnames(metadata), "py_prot", "py_rna", "px_cnv"))
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
  suml1 <- summary(lreg1)
  suml2 <- summary(lreg2)

  ## coefficients and p-values
  beta_cnv <- suml2$coefficients["px_cnv", "Estimate"]
  beta_cnv_p <- suml2$coefficients["px_cnv", "Pr(>|t|)"]

  ## F-statistic
  if( "fstatistic" %in% names(suml1) ){
    f_val1 <- suml1$fstatistic[[1]]
    f_pval1 <- 1-pf(suml1$fstatistic[[1]], suml1$fstatistic[[2]], suml1$fstatistic[[3]])
  } else {
    f_val1 <- NA; f_pval1 <- NA;
  }

  if( "fstatistic" %in% names(suml2) ){
    f_val2 <- suml2$fstatistic[[1]]
    f_pval2 <- 1-pf(suml2$fstatistic[[1]], suml2$fstatistic[[2]], suml2$fstatistic[[3]])
  } else {
    f_val2 <- NA; f_pval2 <- NA;
  }

  ## adjusted R-squared
  r2_adj1 <- suml1$adj.r.squared
  r2_adj2 <- suml2$adj.r.squared


  ## Likelihood Ratio Test
  ## H0: Restricted model (full or null) is statistically better than Unrestricted (alternative) model
  lrt <- lrtest(lreg1, lreg2)
  logl1 <- lrt$LogLik[1]
  logl2 <- lrt$LogLik[2]
  lrt_s <- lrt[2,4]
  lrt_p <- lrt[2,5]


  return( c(px, py, beta_cnv, beta_cnv_p, f_val1, f_pval1, f_val2, f_pval2, r2_adj1, r2_adj2, logl1, logl2, lrt_s, lrt_p) )

}





# compute linear regressions
cl <- makeCluster(detectCores())
clusterExport(cl, c("proteomics", "rna", "cnv", "metadata", "linear_model"))
clusterEvalQ(cl, library(lmtest))

protein_associations <- parApply( cl, as.matrix(protein_interactions), 1, function(x) linear_model(x[[2]], x[[1]], proteomics, rna, cnv, metadata) )

stopCluster(cl)


protein_associations <- protein_associations %>%
    t() %>%
    as.tibble() %>%
    setNames(c("controlling", "controlled", "beta_cnv", "beta_cnv_p", "f_val1", "f_pval1", "f_val2", "f_pval2", "r2_adj1", "r2_adj2", "logl1", "logl2", "lrt_s", "lrt_p"))

protein_associations <- protein_associations %>%
    mutate_at(colnames(protein_associations)[3:14], as.numeric) %>%
    mutate(fdr = p.adjust(p=lrt_p, method="BH")) %>%
    dplyr::arrange(fdr)
saveRDS(protein_associations, file = "./files/protein_associations_tcga_ccle_cnv_reg_genome_wide.rds")


cat("Total associations FDR < 0.05:", sum(protein_associations$fdr<0.05), "\n")
#10,843






save(list=ls(), file="./r_workspaces/protein_associations_tcga_ccle_cnv_reg_genome_wide.RData")
