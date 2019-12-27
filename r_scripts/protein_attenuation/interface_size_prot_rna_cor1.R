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
protein_attenuation <- read_tsv("./files/protein_attenuation_interface.txt") %>%
  mutate_at(vars(contains("class")), as.factor) %>%
  mutate_if(is.factor, fct_infreq) %>%
  mutate(interface1 = if_else(AA_PERC_INTERFACE > 0.1, "big" , "small"), interface2 = if_else(AA_PERC_INTERFACE < 0.05, "small" , if_else(AA_PERC_INTERFACE > 0.18, "big", "intermediate")))


# import expression data
expression_data <- read_tsv("./files/protein_attenuation_cnv_rna_protein.txt")


# merge protein attenuation info with expression data
attenuation_expr <- protein_attenuation %>%
  dplyr::select(gene = PROTEIN_GENE, attenuation, class, interface1, interface2) %>%
  inner_join(expression_data, by = "gene")


# import protein interactions
all_interactions <- read_tsv("./files/biogrid_corum_ER_pairs.txt") %>%
  mutate(pair = paste(interactor.A, interactor.B, sep="_"))

corum_interactions <- all_interactions %>%
  filter(experiment.type == "corum")

biogrid_interactions <- all_interactions %>%
  filter(!(experiment.type == "corum" | experiment.type == "curated_membrane_complexes"))





# compute correlation of protein/rna abundance between gene pairs by interface group

# small interface
small_interface_prot <- attenuation_expr %>%
  filter(interface1 == "small") %>%
  dplyr::select(gene, sample, prot_log2FC) %>%
  spread(key = "sample", value="prot_log2FC") %>%
  as.data.frame() %>%
  column_to_rownames(var = "gene") %>%
  t() %>%
  rcorr(x = ., type="pearson")
  #cor(x = ., use="pairwise.complete.obs")
small_interface_prot <- small_interface_prot$r
small_interface_prot[upper.tri(small_interface_prot)] <- NA
small_interface_prot <- small_interface_prot %>%
  as.data.frame() %>%
  rownames_to_column(var="g1") %>%
  as.tibble() %>%
  gather(key="g2", value="r", -g1) %>%
  na.exclude() %>%
  mutate(interface = "small", data = "protein") %>%
  mutate(corum = if_else( (paste(g1,g2,sep="_") %in% corum_interactions$pair | paste(g2,g1,sep="_") %in% corum_interactions$pair), 1, 0 ))


small_interface_rna <- attenuation_expr %>%
  filter(interface1 == "small") %>%
  dplyr::select(gene, sample, rna_log2CPM) %>%
  spread(key = "sample", value="rna_log2CPM") %>%
  as.data.frame() %>%
  column_to_rownames(var = "gene") %>%
  t() %>%
  rcorr(x = ., type="pearson")
  #cor(x = ., use="pairwise.complete.obs")
small_interface_rna <- small_interface_rna$r
small_interface_rna[upper.tri(small_interface_rna)] <- NA
small_interface_rna <- small_interface_rna %>%
  as.data.frame() %>%
  rownames_to_column(var="g1") %>%
  as.tibble() %>%
  gather(key="g2", value="r", -g1) %>%
  na.exclude() %>%
  mutate(interface = "small", data = "rna") %>%
  mutate(corum = if_else( (paste(g1,g2,sep="_") %in% corum_interactions$pair | paste(g2,g1,sep="_") %in% corum_interactions$pair), 1, 0 ))


# big interface
big_interface_prot <- attenuation_expr %>%
  filter(interface1 == "big") %>%
  dplyr::select(gene, sample, prot_log2FC) %>%
  spread(key = "sample", value="prot_log2FC") %>%
  as.data.frame() %>%
  column_to_rownames(var = "gene") %>%
  t() %>%
  rcorr(x = ., type="pearson")
  #cor(x = ., use="pairwise.complete.obs")
big_interface_prot <- big_interface_prot$r
big_interface_prot[upper.tri(big_interface_prot)] <- NA
big_interface_prot <- big_interface_prot %>%
  as.data.frame() %>%
  rownames_to_column(var="g1") %>%
  as.tibble() %>%
  gather(key="g2", value="r", -g1) %>%
  na.exclude() %>%
  mutate(interface = "big", data = "protein") %>%
  mutate(corum = if_else( (paste(g1,g2,sep="_") %in% corum_interactions$pair | paste(g2,g1,sep="_") %in% corum_interactions$pair), 1, 0 ))


big_interface_rna <- attenuation_expr %>%
  filter(interface1 == "big") %>%
  dplyr::select(gene, sample, rna_log2CPM) %>%
  spread(key = "sample", value="rna_log2CPM") %>%
  as.data.frame() %>%
  column_to_rownames(var = "gene") %>%
  t() %>%
  rcorr(x = ., type="pearson")
  #cor(x = ., use="pairwise.complete.obs")
big_interface_rna <- big_interface_rna$r
big_interface_rna[upper.tri(big_interface_rna)] <- NA
big_interface_rna <- big_interface_rna %>%
  as.data.frame() %>%
  rownames_to_column(var="g1") %>%
  as.tibble() %>%
  gather(key="g2", value="r", -g1) %>%
  na.exclude() %>%
  mutate(interface = "big", data = "rna") %>%
  mutate(corum = if_else( (paste(g1,g2,sep="_") %in% corum_interactions$pair | paste(g2,g1,sep="_") %in% corum_interactions$pair), 1, 0 ))




# -- ROC curve

# small interfaces
small_interface_prot_pred <- prediction(small_interface_prot$r, small_interface_prot$corum)
small_interface_prot_perf <- performance(small_interface_prot_pred, measure = "tpr", x.measure = "fpr")
small_interface_prot_auc <- performance(small_interface_prot_pred, measure = "auc")@y.values[[1]]

small_interface_rna_pred <- prediction(small_interface_rna$r, small_interface_rna$corum)
small_interface_rna_perf <- performance(small_interface_rna_pred, measure = "tpr", x.measure = "fpr")
small_interface_rna_auc <- performance(small_interface_rna_pred, measure = "auc")@y.values[[1]]



# big interfaces
big_interface_prot_pred <- prediction(big_interface_prot$r, big_interface_prot$corum)
big_interface_prot_perf <- performance(big_interface_prot_pred, measure = "tpr", x.measure = "fpr")
big_interface_prot_auc <- performance(big_interface_prot_pred, measure = "auc")@y.values[[1]]

big_interface_rna_pred <- prediction(big_interface_rna$r, big_interface_rna$corum)
big_interface_rna_perf <- performance(big_interface_rna_pred, measure = "tpr", x.measure = "fpr")
big_interface_rna_auc <- performance(big_interface_rna_pred, measure = "auc")@y.values[[1]]


#pdf("./plots/protein_attenuation_tcga_ccle/roc_curve.pdf")
#plot(small_interface_rna_perf)
#dev.off()



roc_data <- bind_rows(
  tibble( fpr=small_interface_prot_perf@x.values[[1]], tpr=small_interface_prot_perf@y.values[[1]], interface = "small", data = "protein"),
  tibble( fpr=small_interface_rna_perf@x.values[[1]], tpr=small_interface_rna_perf@y.values[[1]], interface = "small", data = "rna"),
  tibble( fpr=big_interface_prot_perf@x.values[[1]], tpr=big_interface_prot_perf@y.values[[1]], interface = "big", data = "protein"),
  tibble( fpr=big_interface_rna_perf@x.values[[1]], tpr=big_interface_rna_perf@y.values[[1]], interface = "big", data = "rna"),
)


roc_curve1 <- roc_data %>%
  filter(interface == "small") %>%
  ggplot(mapping=aes(x = fpr, y = tpr, color = data)) +
  geom_line() +
  theme_classic() +
  theme(
    axis.title=element_text(colour="black", size=15),
    axis.text=element_text(colour="black", size=12),
    strip.text=element_text(colour="black", size=13),
    strip.background=element_blank(),
    legend.text=element_text(size=12),
    legend.title=element_text(size=13)) +
  scale_color_manual(values = c("black", "chocolate4"), labels = c("Protein", "RNA"), name = "Data") +
  labs(x = "FPR", y = "TPR", title = "Small interfaces")
ggsave("protein_pairs_roc_curve_small_interfaces.png", plot=roc_curve1, path="./plots/protein_attenuation_tcga_ccle/", height=4, width=5)
unlink("protein_pairs_roc_curve_small_interfaces.png")

roc_curve2 <- roc_data %>%
  filter(interface == "big") %>%
  ggplot(mapping=aes(x = fpr, y = tpr, color = data)) +
  geom_line() +
  theme_classic() +
  theme(
    axis.title=element_text(colour="black", size=15),
    axis.text=element_text(colour="black", size=12),
    strip.text=element_text(colour="black", size=13),
    strip.background=element_blank(),
    legend.text=element_text(size=12),
    legend.title=element_text(size=13)) +
  scale_color_manual(values = c("black", "chocolate4"), labels = c("Protein", "RNA"), name = "Data") +
  labs(x = "FPR", y = "TPR", title = "Big interfaces")
ggsave("protein_pairs_roc_curve_big_interfaces.png", plot=roc_curve2, path="./plots/protein_attenuation_tcga_ccle/", height=4, width=5)
unlink("protein_pairs_roc_curve_big_interfaces.png")
