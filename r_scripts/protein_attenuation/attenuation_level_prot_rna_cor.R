# Abel Sousa 2018

# PhD project
# Study of Control of Protein Abundance using Multi-omics Data from Cancer Samples



# -- Correlation of Protein/RNA abundance between gene pairs by attenuation group
# -- Distinction between interacting/non-interacting proteins



suppressMessages(library(tidyverse))
suppressMessages(library(ggpubr))
suppressMessages(library(RColorBrewer))
suppressMessages(library(Hmisc))
suppressMessages(library(ROCR))

options(bitmapType = "cairo")



# import protein attenuation data
protein_attenuation <- read_tsv("./files/protein_attenuation.txt") %>%
  mutate_at(vars(contains("class")), as.factor) %>%
  mutate_if(is.factor, fct_infreq)


# import expression data
expression_data <- read_tsv("./files/protein_attenuation_cnv_rna_protein.txt")



# merge protein attenuation info with expression data
attenuation_expr <- protein_attenuation %>%
  dplyr::select(gene, attenuation, class) %>%
  inner_join(expression_data, by = "gene")




# compute correlation of protein/rna abundance between gene pairs by attenuation group

# highly-attenuated proteins
h_attenuated_prot <- attenuation_expr %>%
  filter(class == "high-attenuated-protein") %>%
  dplyr::select(gene, sample, prot_log2FC) %>%
  spread(key = "sample", value="prot_log2FC") %>%
  as.data.frame() %>%
  column_to_rownames(var = "gene") %>%
  t() %>%
  rcorr(x = ., type="pearson")
  #cor(x = ., use="pairwise.complete.obs")
h_attenuated_prot <- h_attenuated_prot$r
h_attenuated_prot[upper.tri(h_attenuated_prot)] <- NA
h_attenuated_prot <- h_attenuated_prot %>%
  as.data.frame() %>%
  rownames_to_column(var="g1") %>%
  as.tibble() %>%
  gather(key="g2", value="r", -g1) %>%
  na.exclude() %>%
  mutate(class = "highly-attenuated", data = "protein")


h_attenuated_rna <- attenuation_expr %>%
  filter(class == "high-attenuated-protein") %>%
  dplyr::select(gene, sample, rna_log2CPM) %>%
  spread(key = "sample", value="rna_log2CPM") %>%
  as.data.frame() %>%
  column_to_rownames(var = "gene") %>%
  t() %>%
  rcorr(x = ., type="pearson")
  #cor(x = ., use="pairwise.complete.obs")
h_attenuated_rna <- h_attenuated_rna$r
h_attenuated_rna[upper.tri(h_attenuated_rna)] <- NA
h_attenuated_rna <- h_attenuated_rna %>%
  as.data.frame() %>%
  rownames_to_column(var="g1") %>%
  as.tibble() %>%
  gather(key="g2", value="r", -g1) %>%
  na.exclude() %>%
  mutate(class = "highly-attenuated", data = "rna")



# lowly-attenuated proteins
l_attenuated_prot <- attenuation_expr %>%
  filter(class == "low-attenuated-protein") %>%
  dplyr::select(gene, sample, prot_log2FC) %>%
  spread(key = "sample", value="prot_log2FC") %>%
  as.data.frame() %>%
  column_to_rownames(var = "gene") %>%
  t() %>%
  rcorr(x = ., type="pearson")
  #cor(x = ., use="pairwise.complete.obs")
l_attenuated_prot <- l_attenuated_prot$r
l_attenuated_prot[upper.tri(l_attenuated_prot)] <- NA
l_attenuated_prot <- l_attenuated_prot %>%
  as.data.frame() %>%
  rownames_to_column(var="g1") %>%
  as.tibble() %>%
  gather(key="g2", value="r", -g1) %>%
  na.exclude() %>%
  mutate(class = "lowly-attenuated", data = "protein")


l_attenuated_rna <- attenuation_expr %>%
  filter(class == "low-attenuated-protein") %>%
  dplyr::select(gene, sample, rna_log2CPM) %>%
  spread(key = "sample", value="rna_log2CPM") %>%
  as.data.frame() %>%
  column_to_rownames(var = "gene") %>%
  t() %>%
  rcorr(x = ., type="pearson")
  #cor(x = ., use="pairwise.complete.obs")
l_attenuated_rna <- l_attenuated_rna$r
l_attenuated_rna[upper.tri(l_attenuated_rna)] <- NA
l_attenuated_rna <- l_attenuated_rna %>%
  as.data.frame() %>%
  rownames_to_column(var="g1") %>%
  as.tibble() %>%
  gather(key="g2", value="r", -g1) %>%
  na.exclude() %>%
  mutate(class = "lowly-attenuated", data = "rna")



# non-attenuated proteins
non_attenuated_prot <- attenuation_expr %>%
  filter(class == "non-attenuated") %>%
  dplyr::select(gene, sample, prot_log2FC) %>%
  spread(key = "sample", value="prot_log2FC") %>%
  as.data.frame() %>%
  column_to_rownames(var = "gene") %>%
  t() %>%
  rcorr(x = ., type="pearson")
  #cor(x = ., use="pairwise.complete.obs")
non_attenuated_prot <- non_attenuated_prot$r
non_attenuated_prot[upper.tri(non_attenuated_prot)] <- NA
non_attenuated_prot <- non_attenuated_prot %>%
  as.data.frame() %>%
  rownames_to_column(var="g1") %>%
  as.tibble() %>%
  gather(key="g2", value="r", -g1) %>%
  na.exclude() %>%
  mutate(class = "non-attenuated", data = "protein")


non_attenuated_rna <- attenuation_expr %>%
  filter(class == "non-attenuated") %>%
  dplyr::select(gene, sample, rna_log2CPM) %>%
  spread(key = "sample", value="rna_log2CPM") %>%
  as.data.frame() %>%
  column_to_rownames(var = "gene") %>%
  t() %>%
  rcorr(x = ., type="pearson")
  #cor(x = ., use="pairwise.complete.obs")
non_attenuated_rna <- non_attenuated_rna$r
non_attenuated_rna[upper.tri(non_attenuated_rna)] <- NA
non_attenuated_rna <- non_attenuated_rna %>%
  as.data.frame() %>%
  rownames_to_column(var="g1") %>%
  as.tibble() %>%
  gather(key="g2", value="r", -g1) %>%
  na.exclude() %>%
  mutate(class = "non-attenuated", data = "rna")



# import protein interactions
all_interactions <- read_tsv("./files/biogrid_corum_ER_pairs.txt") %>%
  mutate(pair = paste(interactor.A, interactor.B, sep="_"))

corum_interactions <- all_interactions %>%
  filter(experiment.type == "corum")

biogrid_interactions <- all_interactions %>%
  filter(!(experiment.type == "corum" | experiment.type == "curated_membrane_complexes"))


# merge all protein pairs
protein_pairs <- bind_rows(h_attenuated_prot, h_attenuated_rna, l_attenuated_prot, l_attenuated_rna, non_attenuated_prot, non_attenuated_rna) %>%
  mutate(type = "all")


# add protein interactions reported on corum and biogrid
protein_pairs2 <- protein_pairs %>%
  bind_rows(protein_pairs %>% filter(paste(g1,g2,sep="_") %in% biogrid_interactions$pair | paste(g2,g1,sep="_") %in% biogrid_interactions$pair) %>% mutate(type = "biogrid")) %>%
  bind_rows(protein_pairs %>% filter(paste(g1,g2,sep="_") %in% corum_interactions$pair | paste(g2,g1,sep="_") %in% corum_interactions$pair) %>% mutate(type = "corum"))



protein_pairs_boxpl <- protein_pairs2 %>%
  mutate_at(vars(contains("class")), as.factor) %>%
  mutate_if(is.factor, fct_infreq) %>%
  ggplot(mapping=aes(x=class, y=r, fill=class, color = data)) +
  geom_boxplot(outlier.shape=NA) +
  theme_classic() +
  #stat_compare_means(method = "wilcox.test", comparisons = list( c(1, 3)), label.y=5 ) +
  facet_wrap( ~ type, labeller=labeller(type = c("all" = "All\nprotein pairs", "biogrid" = "BioGRID\ninteractions", "corum" = "CORUM\ninteractions"))) +
  theme(
    axis.title=element_text(colour="black", size=15),
    axis.text.x=element_blank(),
    axis.text.y=element_text(colour="black", size=12),
    axis.ticks.x=element_blank(),
    strip.text=element_text(colour="black", size=13),
    strip.background=element_blank(),
    #legend.position="bottom",
    legend.text=element_text(size=12),
    legend.title=element_text(size=13)) +
  scale_fill_brewer(type = "seq", palette = "Blues", labels = c("Non-attenuated", "Lowly-attenuated", "Highly-attenuated"), name = "Attenuation") +
  scale_color_manual(values = c("black", "chocolate4"), labels = c("Protein", "RNA"), name = "Data") +
  labs(x = "Attenuation", y = "Pearson's r")
ggsave("protein_pairs_prot_rna_cor_boxpl.png", plot=protein_pairs_boxpl, path="./plots/protein_attenuation_tcga_ccle/", height=5, width=6)
unlink("protein_pairs_prot_rna_cor_boxpl.png")


# -- ROC curve

# corum complexes
ppairs_roc <- protein_pairs %>%
  mutate(corum = if_else( (paste(g1,g2,sep="_") %in% corum_interactions$pair | paste(g2,g1,sep="_") %in% corum_interactions$pair), 1, 0 ))


# highly-attenuated genes
h_att_prot <- ppairs_roc %>%
  filter(class == "highly-attenuated", data == "protein")

h_att_prot_pred <- prediction(h_att_prot$r, h_att_prot$corum)
h_att_prot_perf <- performance(h_att_prot_pred, measure = "tpr", x.measure = "fpr")
h_att_prot_auc <- performance(h_att_prot_pred, measure = "auc")@y.values[[1]]

h_att_rna <- ppairs_roc %>%
  filter(class == "highly-attenuated", data == "rna")

h_att_rna_pred <- prediction(h_att_rna$r, h_att_rna$corum)
h_att_rna_perf <- performance(h_att_rna_pred, measure = "tpr", x.measure = "fpr")
h_att_rna_auc <- performance(h_att_rna_pred, measure = "auc")@y.values[[1]]

#pdf("./plots/protein_attenuation_tcga_ccle/roc_curve.pdf")
#plot(h_att_prot_perf)
#dev.off()

# lowly-attenuated genes
l_att_prot <- ppairs_roc %>%
  filter(class == "lowly-attenuated", data == "protein")

l_att_prot_pred <- prediction(l_att_prot$r, l_att_prot$corum)
l_att_prot_perf <- performance(l_att_prot_pred, measure = "tpr", x.measure = "fpr")
l_att_prot_auc <- performance(l_att_prot_pred, measure = "auc")@y.values[[1]]

l_att_rna <- ppairs_roc %>%
  filter(class == "lowly-attenuated", data == "rna")

l_att_rna_pred <- prediction(l_att_rna$r, l_att_rna$corum)
l_att_rna_perf <- performance(l_att_rna_pred, measure = "tpr", x.measure = "fpr")
l_att_rna_auc <- performance(l_att_rna_pred, measure = "auc")@y.values[[1]]


# non-attenuated genes
non_att_prot <- ppairs_roc %>%
  filter(class == "non-attenuated", data == "protein")

non_att_prot_pred <- prediction(non_att_prot$r, non_att_prot$corum)
non_att_prot_perf <- performance(non_att_prot_pred, measure = "tpr", x.measure = "fpr")
non_att_prot_auc <- performance(non_att_prot_pred, measure = "auc")@y.values[[1]]

non_att_rna <- ppairs_roc %>%
  filter(class == "non-attenuated", data == "rna")

non_att_rna_pred <- prediction(non_att_rna$r, non_att_rna$corum)
non_att_rna_perf <- performance(non_att_rna_pred, measure = "tpr", x.measure = "fpr")
non_att_rna_auc <- performance(non_att_rna_pred, measure = "auc")@y.values[[1]]

roc_data <- bind_rows(
  tibble( fpr=h_att_prot_perf@x.values[[1]], tpr=h_att_prot_perf@y.values[[1]], class = "high-att", data = "protein"),
  tibble( fpr=h_att_rna_perf@x.values[[1]], tpr=h_att_rna_perf@y.values[[1]], class = "high-att", data = "rna"),
  tibble( fpr=l_att_prot_perf@x.values[[1]], tpr=l_att_prot_perf@y.values[[1]], class = "low-att", data = "protein"),
  tibble( fpr=l_att_rna_perf@x.values[[1]], tpr=l_att_rna_perf@y.values[[1]], class = "low-att", data = "rna"),
  tibble( fpr=non_att_prot_perf@x.values[[1]], tpr=non_att_prot_perf@y.values[[1]], class = "non-att", data = "protein"),
  tibble( fpr=non_att_rna_perf@x.values[[1]], tpr=non_att_rna_perf@y.values[[1]], class = "non-att", data = "rna")
)


roc_auc <- bind_rows(
  tibble( auc=h_att_prot_auc, class = "high-att", data = "protein"),
  tibble( auc=h_att_rna_auc, class = "high-att", data = "rna"),
  tibble( auc=l_att_prot_auc, class = "low-att", data = "protein"),
  tibble( auc=l_att_rna_auc, class = "low-att", data = "rna"),
  tibble( auc=non_att_prot_auc, class = "non-att", data = "protein"),
  tibble( auc=non_att_rna_auc, class = "non-att", data = "rna")) %>%
  mutate(auc = round(auc, 2))



roc_curve1 <- roc_data %>%
  #mutate_at(vars(contains("class")), as.factor) %>%
  #mutate_if(is.factor, fct_infreq) %>%
  filter(class == "high-att") %>%
  ggplot(mapping=aes(x = fpr, y = tpr, color = data)) +
  geom_line() +
  geom_text(data = roc_auc %>% filter(class == "high-att") %>% mutate(x = c(0.05,0.05), y = c(1,0.95)), mapping = aes(x = x, y = y, label = auc, color = data), size = 4,  inherit.aes = FALSE) +
  theme_classic() +
  theme(
    axis.title=element_text(colour="black", size=15),
    axis.text=element_text(colour="black", size=12),
    strip.text=element_text(colour="black", size=13),
    strip.background=element_blank(),
    legend.text=element_text(size=12),
    legend.title=element_text(size=13)) +
  scale_color_manual(values = c("black", "chocolate4"), labels = c("Protein", "RNA"), name = "Data") +
  labs(x = "FPR", y = "TPR", title = "Highly-attenuated")
ggsave("protein_pairs_roc_curve_h_att.png", plot=roc_curve1, path="./plots/protein_attenuation_tcga_ccle/", height=4, width=5)
unlink("protein_pairs_roc_curve_h_att.png")
ggsave("protein_pairs_roc_curve_h_att.pdf", plot=roc_curve1, path="./plots/protein_attenuation_tcga_ccle/", height=4, width=5)
unlink("protein_pairs_roc_curve_h_att.pdf")

roc_curve2 <- roc_data %>%
  #mutate_at(vars(contains("class")), as.factor) %>%
  #mutate_if(is.factor, fct_infreq) %>%
  filter(class == "low-att") %>%
  ggplot(mapping=aes(x = fpr, y = tpr, color = data)) +
  geom_line() +
  geom_text(data = roc_auc %>% filter(class == "low-att") %>% mutate(x = c(0.05,0.05), y = c(1,0.95)), mapping = aes(x = x, y = y, label = auc, color = data), size = 4,  inherit.aes = FALSE) +
  theme_classic() +
  theme(
    axis.title=element_text(colour="black", size=15),
    axis.text=element_text(colour="black", size=12),
    strip.text=element_text(colour="black", size=13),
    strip.background=element_blank(),
    legend.text=element_text(size=12),
    legend.title=element_text(size=13)) +
  scale_color_manual(values = c("black", "chocolate4"), labels = c("Protein", "RNA"), name = "Data") +
  labs(x = "FPR", y = "TPR", title = "Lowly-attenuated")
ggsave("protein_pairs_roc_curve_l_att.png", plot=roc_curve2, path="./plots/protein_attenuation_tcga_ccle/", height=4, width=5)
unlink("protein_pairs_roc_curve_l_att.png")
ggsave("protein_pairs_roc_curve_l_att.pdf", plot=roc_curve2, path="./plots/protein_attenuation_tcga_ccle/", height=4, width=5)
unlink("protein_pairs_roc_curve_l_att.pdf")

roc_curve3 <- roc_data %>%
  #mutate_at(vars(contains("class")), as.factor) %>%
  #mutate_if(is.factor, fct_infreq) %>%
  filter(class == "non-att") %>%
  ggplot(mapping=aes(x = fpr, y = tpr, color = data)) +
  geom_line() +
  geom_text(data = roc_auc %>% filter(class == "non-att") %>% mutate(x = c(0.05,0.05), y = c(1,0.95)), mapping = aes(x = x, y = y, label = auc, color = data), size = 4,  inherit.aes = FALSE) +
  theme_classic() +
  theme(
    axis.title=element_text(colour="black", size=15),
    axis.text=element_text(colour="black", size=12),
    strip.text=element_text(colour="black", size=13),
    strip.background=element_blank(),
    legend.text=element_text(size=12),
    legend.title=element_text(size=13)) +
  scale_color_manual(values = c("black", "chocolate4"), labels = c("Protein", "RNA"), name = "Data") +
  labs(x = "FPR", y = "TPR", title = "Non-attenuated")
ggsave("protein_pairs_roc_curve_non_att.png", plot=roc_curve3, path="./plots/protein_attenuation_tcga_ccle/", height=4, width=5)
unlink("protein_pairs_roc_curve_non_att.png")
ggsave("protein_pairs_roc_curve_non_att.pdf", plot=roc_curve3, path="./plots/protein_attenuation_tcga_ccle/", height=4, width=5)
unlink("protein_pairs_roc_curve_non_att.pdf")





save(list=ls(), file = "./r_workspaces/attenuation_level_prot_rna_cor.RData")
