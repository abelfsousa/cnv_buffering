# Abel Sousa 2018

# PhD project
# Study of Control of Protein Abundance using Multi-omics Data from Cancer Samples

# -- CCLE and CPATC datasets

# -- Compute Linear Regressions Between Phosphosites and Protein Abundances

# -- All Covars are Used in the Regression Models

# -- Linear Models are compared using a Likelihood Ratio Test



suppressMessages(library(lmtest))
suppressMessages(library(parallel))
suppressMessages(library(tidyverse))
options(bitmapType = "cairo")



# load datasets
proteomics <- read_tsv("./files/proteomicsQ_cptac_cellLines.txt") %>%
    gather(key="sample", value="prot_log2FC", -gene) %>%
    arrange(sample)

phospho <- read_tsv("./files/phosphoproteomicsQ_cptac_cellLines.txt") %>%
    gather(key="sample", value="phospho_log2FC", -phospho_site) %>%
    arrange(sample)

rna <- read_tsv("./files/rna_tcga_cellLines.txt") %>%
    gather(key="sample", value="rna_log2CPM", -gene) %>%
    arrange(sample)

cnv <- read_tsv("./files/cnv_tcga_cellLines.txt") %>%
    gather(key="sample", value="cnv_gistic2", -gene) %>%
    arrange(sample)

metadata <- read_tsv("./files/metadata_cptac_cellLines.txt") %>%
    arrange(sample)



common_samples <- unique(proteomics$sample) %>%
    intersect(unique(phospho$sample)) %>%
    intersect(unique(rna$sample)) %>%
    intersect(unique(cnv$sample)) %>%
    intersect(unique(metadata$sample))
write.table(common_samples, file = "./files/cptac_cellLines_samples_prot_phospho_rna_cnv.txt", sep = "\t", quote = F, row.names=F, col.names=F)



proteomics <- proteomics %>%
    filter(sample %in% common_samples) %>%
    group_by(gene) %>%
    filter( sum(!is.na(prot_log2FC)) > n()*0.5 ) %>%
    mutate(prot_log2FC = scale(prot_log2FC)[,1]) %>%
    #mutate(prot_log2FC = (prot_log2FC-min(prot_log2FC, na.rm=T))/(max(prot_log2FC, na.rm=T)-min(prot_log2FC, na.rm=T))) %>%
    ungroup() %>%
    spread(key="sample", value="prot_log2FC") %>%
    as.data.frame()
rownames(proteomics) <- proteomics$gene
proteomics <- proteomics[, -c(1)]



phospho <- phospho %>%
    filter(sample %in% common_samples) %>%
    group_by(phospho_site) %>%
    filter( sum(!is.na(phospho_log2FC)) > n()*0.5 ) %>%
    mutate(phospho_log2FC = scale(phospho_log2FC)[,1]) %>%
    #mutate(phospho_log2FC = (phospho_log2FC-min(phospho_log2FC, na.rm=T))/(max(phospho_log2FC, na.rm=T)-min(phospho_log2FC, na.rm=T))) %>%
    ungroup() %>%
    spread(key="sample", value="phospho_log2FC") %>%
    as.data.frame()
rownames(phospho) <- phospho$phospho_site
phospho <- phospho[, -c(1)]



rna <- rna %>%
    filter(sample %in% common_samples) %>%
    group_by(gene) %>%
    mutate(rna_log2CPM = scale(rna_log2CPM)[,1]) %>%
    #mutate(rna_log2CPM = (rna_log2CPM-min(rna_log2CPM))/(max(rna_log2CPM)-min(rna_log2CPM))) %>%
    ungroup() %>%
    spread(key="sample", value="rna_log2CPM") %>%
    as.data.frame()
rownames(rna) <- rna$gene
rna <- rna[, -c(1)]


cnv <- cnv %>%
    filter(sample %in% common_samples) %>%
    spread(key="sample", value="cnv_gistic2") %>%
    as.data.frame()
rownames(cnv) <- cnv$gene
cnv <- cnv[, -c(1)]


metadata <- metadata %>%
    filter(sample %in% common_samples) %>%
    mutate(age = scale(age)[,1]) %>%
    #mutate(age = (age-min(age))/(max(age)-min(age))) %>%
    dplyr::select(-c(cancer, proteomics)) %>%
    as.data.frame()
rownames(metadata) <- metadata$sample
metadata <- metadata[, -c(1)]






# load protein interactions
all_interactions <- read_tsv("./files/biogrid_corum_ER_pairs.txt")

keep_interactions <- all_interactions %>%
    dplyr::select(1,2) %>%
    distinct() %>%
    as.data.frame()
nrow(keep_interactions)
#572,856


gene_phosphosites <- str_split(string=rownames(phospho), pattern="_", simplify=T) %>%
    as.data.frame(stringsAsFactors=F) %>%
    rename(gene=V1, phosphosite=V2)
nrow(gene_phosphosites)
#5733


# select only the proteins pairs for which interactor.A %in% phosphoproteomics & proteomics & CNV and interactor.B %in% proteomics & transcriptomics
keep <- apply( keep_interactions, 1, function(ppair) ppair[1] %in% gene_phosphosites$gene & ppair[1] %in% rownames(proteomics) & ppair[1] %in% rownames(cnv) & ppair[2] %in% rownames(proteomics) & ppair[2] %in% rownames(rna) )
keep_interactions <- keep_interactions[keep, ]


# convert protein-protein interactions in phosphosite-protein interactions
keep_interactions <- keep_interactions %>%
    inner_join(gene_phosphosites, by=c("interactor.A" = "gene")) %>%
    mutate(phosphosite = paste(interactor.A, phosphosite, sep="_")) %>%
    dplyr::select(interactor.A, phosphosite, interactor.B) %>%
    rename(phosphosite.A = phosphosite)


cat("Interactions to calculate:", nrow(keep_interactions), "\n")
#315,772




linear_model <- function(px, phospho_x, py, proteomics, phospho, rna, cnv, metadata){
    #print(paste(px, phospho_x, py, sep=" "))

	py_prot <- t(proteomics[py, ])
	py_rna <- t(rna[py, ])
    px_prot <- t(proteomics[px, ])
	px_phospho <- t(phospho[phospho_x, ])
    px_cnv <- t(cnv[px, ])


    ## define covars matrix
	## bind prot/rna of Py and prot/phospho of Px
	covars <- setNames( cbind(metadata, py_prot, py_rna, px_cnv, px_prot, px_phospho), c(colnames(metadata), "py_prot", "py_rna", "px_cnv", "px_prot", "px_phospho") )
    covars <- covars[, c("py_prot", "py_rna", "px_prot", "px_cnv", "age", "gender", "batch", "px_phospho")]

    ## remove nas
    covars <- na.exclude(covars)


    ## remove covariates with less than two levels
    covars <- covars[, apply(covars, 2, function(x) length(levels(as.factor(x))) >= 2), drop=F]


    ## fit the models
    #restricted/reduced model (without phospho)
    lreg1_f <- as.formula( paste( "py_prot", "~", paste(colnames(covars[, -match(c("py_prot", "px_phospho"), colnames(covars)), drop=F]), collapse="+") ) )
    lreg1 <- lm(lreg1_f, data = covars)


    #unrestricted model (with phospho)
    lreg2_f <- update.formula(formula(lreg1), ~ . + px_phospho)
    lreg2 <- lm(lreg2_f, data = covars)

    ## summary the models
    suml1 <- summary(lreg1)
    suml2 <- summary(lreg2)

    ## coefficients and p-values
    beta_phospho <- suml2$coefficients["px_phospho", "Estimate"]
    beta_phospho_p <- suml2$coefficients["px_phospho", "Pr(>|t|)"]

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


    return( c(px, phospho_x, py, beta_phospho, beta_phospho_p, f_val1, f_pval1, f_val2, f_pval2, r2_adj1, r2_adj2, logl1, logl2, lrt_s, lrt_p) )

}



# compute linear regressions
cl <- makeCluster(detectCores())
clusterExport(cl, c("proteomics", "phospho", "rna", "cnv", "metadata", "linear_model"))
clusterEvalQ(cl, library(lmtest))

protein_associations <- parApply( cl, keep_interactions, 1, function(x) linear_model(x[[1]], x[[2]], x[[3]], proteomics, phospho, rna, cnv, metadata) )

stopCluster(cl)


protein_associations <- protein_associations %>%
    t() %>%
    as.tibble() %>%
    setNames(c("controlling", "controlling_phx", "controlled", "beta_phospho", "beta_phospho_p", "f_val1", "f_pval1", "f_val2", "f_pval2", "r2_adj1", "r2_adj2", "logl1", "logl2", "lrt_s", "lrt_p"))

protein_associations <- protein_associations %>%
    mutate_at(colnames(protein_associations)[4:15], as.numeric) %>%
    mutate(fdr = p.adjust(p=lrt_p, method="BH")) %>%
    dplyr::arrange(fdr)
write.table(protein_associations, file = "./files/protein_associations_tcga_ccle_phospho_reg.txt", sep = "\t", quote = F, row.names=F)



cat("Total associations FDR < 0.05:", sum(protein_associations$fdr<0.05), "\n")
#11,672


cat("Total associations FDR < 0.01:", sum(protein_associations$fdr<0.01), "\n")
#4386


cat("Total associations FDR < 0.001:", sum(protein_associations$fdr<0.001), "\n")
#1571








save(list=ls(), file="./r_workspaces/protein_associations_tcga_ccle_phospho_reg.RData")
