# Abel Sousa 2018

# PhD project
# Study of Control of Protein Abundance using Multi-omics Data from Cancer Samples


# -- Evaluation of proteomics datasets using ROC curves




suppressMessages(library(tidyverse))
suppressMessages(library(ROCR))
options(bitmapType = "cairo")
set.seed(123)




# load proteomics datasets
proteomics <- read_tsv("./files/cptac_proteomics.txt") %>%
	gather(key="sample", value="log2FC", -gene) %>%
	group_by(gene) %>%
	filter(sum(!is.na(log2FC)) > n()*0.5)

proteomicsQ <- read_tsv("./files/cptac_proteomicsQ.txt") %>%
	gather(key="sample", value="log2FC", -gene) %>%
	group_by(gene) %>%
	filter(sum(!is.na(log2FC)) > n()*0.5)

proteomics2 <- read_tsv("./files/cptac_proteomics2.txt") %>%
	gather(key="sample", value="log2FC", -gene) %>%
	group_by(gene) %>%
	filter(sum(!is.na(log2FC)) > n()*0.5)



# intersect proteins across 3 datasets
proteins <- unique(proteomics$gene) %>%
	intersect(unique(proteomicsQ$gene)) %>%
	intersect(unique(proteomics2$gene))




# load corum protein pairs
# select corum pairs whose interactor.A and B are present in the proteomics datasets
corum_pairs <- read_tsv("./files/corum_protein_pairs.txt") %>%
	filter(interactor.A %in% proteins & interactor.B %in% proteins)



# select 5000 random corum pairs
corum_pairs_random <- corum_pairs[sample(nrow(corum_pairs), 5000), ]



# select 5000 random protein pairs
protein_pairs_random <- data.frame(
	interactor.A = proteins[ sample(length(proteins), length(proteins)) ],
	interactor.B = proteins[ sample(length(proteins), length(proteins)) ],
	stringsAsFactors = F) %>%
	as.tibble() %>%
	distinct() %>%
	filter(interactor.A != interactor.B) %>%
	mutate(experiment.type = "random") %>%
	slice(1:5000)



length(
	intersect(
		paste(corum_pairs_random$interactor.A, corum_pairs_random$interactor.B, sep="_"),
		paste(protein_pairs_random$interactor.A, protein_pairs_random$interactor.B, sep="_")
	)
)
#0




# bind random and corum protein pairs
protein_pairs_corum_random <- bind_rows(corum_pairs_random, protein_pairs_random)





# -- ROC curve using dataset 1 (proteomics)



# bind proteomics data and calculate protein-protein correlations
protein_pairs_cor_proteomics <- protein_pairs_corum_random %>%
	left_join(proteomics, by=c("interactor.A"="gene")) %>%
	rename(proteomics.interactor.A = log2FC) %>%
	left_join(proteomics, by=c("interactor.B"="gene", "sample")) %>%
	rename(proteomics.interactor.B = log2FC) %>%
	group_by(interactor.A, interactor.B) %>%
	do(broom::tidy(cor.test(.$proteomics.interactor.A, .$proteomics.interactor.B))) %>%
	ungroup() %>%
	inner_join(protein_pairs_corum_random, by=c("interactor.A", "interactor.B")) %>%
	na.exclude() %>%
	mutate(log10.p.value = -log10(p.value)) %>%
	mutate(experiment.type = if_else(experiment.type == "corum", 1, 0))


proteomics_pred_r <- prediction(protein_pairs_cor_proteomics$estimate, protein_pairs_cor_proteomics$experiment.type)
proteomics_prf_r <- performance(proteomics_pred_r, measure = "tpr", x.measure = "fpr")
proteomics_auc_r <- performance(proteomics_pred_r, measure = "auc")@y.values[[1]]

proteomics_pred_p <- prediction(protein_pairs_cor_proteomics$log10.p.value, protein_pairs_cor_proteomics$experiment.type)
proteomics_prf_p <- performance(proteomics_pred_p, measure = "tpr", x.measure = "fpr")
proteomics_auc_p <- performance(proteomics_pred_p, measure = "auc")@y.values[[1]]


pdf("./plots/proteomics_roc_curve.pdf")
plot(proteomics_prf_r, main = paste("CPTAC proteomics (without alterations)\nPearson r", paste("AUC: ", proteomics_auc_r), sep="\n"))
plot(proteomics_prf_p, main = paste("CPTAC proteomics (without alterations)\nlog10 p-value", paste("AUC: ", proteomics_auc_p), sep="\n"))
dev.off()





# -- ROC curve using dataset 2 (proteomicsQ)


# bind proteomics data and calculate protein-protein correlations
protein_pairs_cor_proteomicsQ <- protein_pairs_corum_random %>%
	left_join(proteomicsQ, by=c("interactor.A"="gene")) %>%
	rename(proteomics.interactor.A = log2FC) %>%
	left_join(proteomicsQ, by=c("interactor.B"="gene", "sample")) %>%
	rename(proteomics.interactor.B = log2FC) %>%
	group_by(interactor.A, interactor.B) %>%
	do(broom::tidy(cor.test(.$proteomics.interactor.A, .$proteomics.interactor.B))) %>%
	ungroup() %>%
	inner_join(protein_pairs_corum_random, by=c("interactor.A", "interactor.B")) %>%
	na.exclude() %>%
	mutate(log10.p.value = -log10(p.value)) %>%
	mutate(experiment.type = if_else(experiment.type == "corum", 1, 0))


proteomicsQ_pred_r <- prediction(protein_pairs_cor_proteomicsQ$estimate, protein_pairs_cor_proteomicsQ$experiment.type)
proteomicsQ_prf_r <- performance(proteomicsQ_pred_r, measure = "tpr", x.measure = "fpr")
proteomicsQ_auc_r <- performance(proteomicsQ_pred_r, measure = "auc")@y.values[[1]]

proteomicsQ_pred_p <- prediction(protein_pairs_cor_proteomicsQ$log10.p.value, protein_pairs_cor_proteomicsQ$experiment.type)
proteomicsQ_prf_p <- performance(proteomicsQ_pred_p, measure = "tpr", x.measure = "fpr")
proteomicsQ_auc_p <- performance(proteomicsQ_pred_p, measure = "auc")@y.values[[1]]


pdf("./plots/proteomicsQ_roc_curve.pdf")
plot(proteomicsQ_prf_r, main = paste("CPTAC proteomics (quantile normalization)\nPearson r", paste("AUC: ", proteomicsQ_auc_r), sep="\n"))
plot(proteomicsQ_prf_p, main = paste("CPTAC proteomics (quantile normalization)\nlog10 p-value", paste("AUC: ", proteomicsQ_auc_p), sep="\n"))
dev.off()




# -- ROC curve using dataset 3 (proteomics without 8 outlier samples)


# bind proteomics data and calculate protein-protein correlations
protein_pairs_cor_proteomics2 <- protein_pairs_corum_random %>%
	left_join(proteomics2, by=c("interactor.A"="gene")) %>%
	rename(proteomics.interactor.A = log2FC) %>%
	left_join(proteomics2, by=c("interactor.B"="gene", "sample")) %>%
	rename(proteomics.interactor.B = log2FC) %>%
	group_by(interactor.A, interactor.B) %>%
	do(broom::tidy(cor.test(.$proteomics.interactor.A, .$proteomics.interactor.B))) %>%
	ungroup() %>%
	inner_join(protein_pairs_corum_random, by=c("interactor.A", "interactor.B")) %>%
	na.exclude() %>%
	mutate(log10.p.value = -log10(p.value)) %>%
	mutate(experiment.type = if_else(experiment.type == "corum", 1, 0))


proteomics2_pred_r <- prediction(protein_pairs_cor_proteomics2$estimate, protein_pairs_cor_proteomics2$experiment.type)
proteomics2_prf_r <- performance(proteomics2_pred_r, measure = "tpr", x.measure = "fpr")
proteomics2_auc_r <- performance(proteomics2_pred_r, measure = "auc")@y.values[[1]]

proteomics2_pred_p <- prediction(protein_pairs_cor_proteomics2$log10.p.value, protein_pairs_cor_proteomics2$experiment.type)
proteomics2_prf_p <- performance(proteomics2_pred_p, measure = "tpr", x.measure = "fpr")
proteomics2_auc_p <- performance(proteomics2_pred_p, measure = "auc")@y.values[[1]]


pdf("./plots/proteomics2_roc_curve.pdf")
plot(proteomics2_prf_r, main = paste("CPTAC proteomics (without 8 outlier samples)\nPearson r", paste("AUC: ", proteomics2_auc_r), sep="\n"))
plot(proteomics2_prf_p, main = paste("CPTAC proteomics (without 8 outlier samples)\nlog10 p-value", paste("AUC: ", proteomics2_auc_p), sep="\n"))
dev.off()




save(list=ls(), file="./r_workspaces/evaluation_proteomics_cptac_roc_curve.RData")



