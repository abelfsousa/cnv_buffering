# Abel Sousa 2018

# PhD project
# Study of Control of Protein Abundance using Multi-omics Data from Cancer Samples


# -- Protein attenuation
# -- Spearman correlation


suppressMessages(library(mclust))
suppressMessages(library(tidyverse))
suppressMessages(library(RColorBrewer))
suppressMessages(library(gridExtra))
suppressMessages(library(ggpubr))
suppressMessages(library(cowplot))
options(bitmapType = "cairo")
set.seed(123)





# load datasets

merged_data <- read_tsv("./files/protein_attenuation_cnv_rna_protein.txt")




# compute correlation between CNV and RNA by gene across samples
cnv_rna <- merged_data %>%
    group_by(gene) %>%
    do(broom::tidy(cor.test(.$cnv_gistic2, .$rna_log2CPM, method = "spearman"))) %>%
    ungroup()



# compute correlation between CNV and PROTEIN by gene across samples
cnv_protein <- merged_data %>%
    group_by(gene) %>%
    do(broom::tidy(cor.test(.$cnv_gistic2, .$prot_log2FC, method = "spearman"))) %>%
    ungroup()



# merge correlations and p-values
# compute attenuation potential: cor(CNV, RNA) - cor(CNV, Protein)
# gaussian mixture model to select proteins with high attenuation levels
cnv_rna_protein <- inner_join(cnv_rna, cnv_protein, by="gene") %>%
    dplyr::select(gene, estimate.x, p.value.x, estimate.y, p.value.y) %>%
    dplyr::rename(s_cnv_rna=estimate.x, p_cnv_rna=p.value.x, s_cnv_prot=estimate.y, p_cnv_prot=p.value.y) %>%
    mutate(attenuation = s_cnv_rna - s_cnv_prot)


gmm <- mclust::Mclust(cnv_rna_protein$attenuation, G=4)


cnv_rna_protein <- cnv_rna_protein %>%
    mutate(class = gmm$classification) %>%
    mutate(class = if_else(class==1 | class==2, "non-attenuated", if_else(class==3, "low-attenuated-protein", "high-attenuated-protein" ))) %>%
    arrange(desc(attenuation)) %>%
    mutate_at(vars(contains("class")), as.factor) %>%
    mutate_if(is.factor, fct_infreq)
write.table(cnv_rna_protein, "./files/protein_attenuation_spearman.txt", sep="\t", quote=F, row.names=F)



# some stats
attenuation_stats <- cnv_rna_protein %>%
    group_by(class) %>%
    summarise(counts = n(), perc = n()/nrow(cnv_rna_protein), mean = mean(attenuation))



# comparison with previous analysis
cnv_rna_protein1 <- read_tsv("./files/protein_attenuation.txt")

comparison_att <- cnv_rna_protein1 %>%
  dplyr::select(gene, attenuation1 = attenuation, class1 = class) %>%
  inner_join(cnv_rna_protein %>% dplyr::select(gene, attenuation2 = attenuation, class2 = class), by="gene") %>%
  mutate(same_class = if_else(class1 != class2, "0", "1"))


comparison_att_plot <- ggplot(comparison_att) +
  geom_point(mapping = aes(x = attenuation1, y = attenuation2, color = same_class), size=1) +
  theme_classic() +
  coord_fixed() +
  stat_cor(mapping = aes(x = attenuation1, y = attenuation2)) +
  theme(axis.title=element_text(colour="black", size=15),
    axis.text=element_text(colour="black", size=12),
    plot.title = element_blank(),
    legend.text=element_text(size=12),
    legend.title=element_text(size=14),
    #legend.position="none",
    panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  labs(x = "Attenuation (Pearson corr)", y = "Attenuation (Spearman corr)")
ggsave(filename="comparison_pearson_spearman_cor.png", plot = comparison_att_plot, path = "./plots/protein_attenuation_tcga_ccle/", width=5, height=5)
unlink("comparison_pearson_spearman_cor.png")







# plots







correlation_plot1 <- ggplot(cnv_rna_protein %>% filter(s_cnv_rna > -0.2 & s_cnv_prot > -0.2) %>% inner_join(comparison_att[, c("gene", "same_class")], by = "gene")) +
  geom_point(aes(x = s_cnv_rna, y = s_cnv_prot, color = same_class), size= 1) +
  geom_abline(slope=1, intercept=0, linetype=2, color="black") +
  geom_vline(xintercept=0, color="black", size = 0.2, alpha = 0.5) +
  geom_hline(yintercept=0, color="black", size = 0.2, alpha = 0.5) +
  theme_classic() +
  coord_fixed() +
  theme(axis.title=element_text(colour="black", size=15),
    axis.text=element_text(colour="black", size=12),
    plot.title = element_blank(),
    legend.text=element_text(size=12),
    legend.title=element_text(size=14),
    legend.position="bottom",
    panel.border = element_rect(colour = "black", fill=NA, size=1)) +
    scale_x_continuous(limits = c(-0.2, 0.8), breaks=seq(-0.2, 0.8, by = 0.2)) +
    scale_y_continuous(limits = c(-0.2, 0.8), breaks=seq(-0.2, 0.8, by = 0.2)) +
  labs(x = "corr(CNV, RNA)", y = "corr(CNV, Protein)")
ggsave(filename="protein_attenuation_scatterplot_spearman1.png", plot = correlation_plot1, path = "./plots/protein_attenuation_tcga_ccle/", width=5, height=5)
unlink("protein_attenuation_scatterplot_spearman1.png")



correlation_plot2 <- ggplot(cnv_rna_protein %>% filter(s_cnv_rna > -0.2 & s_cnv_prot > -0.2)) +
  geom_point(aes(x = s_cnv_rna, y = s_cnv_prot, color = class), size= 1) +
  geom_density2d(mapping=aes(x=s_cnv_rna, y=s_cnv_prot, group=class), color="grey", size = 0.2) +
  geom_abline(slope=1, intercept=0, linetype=2, color="black") +
  geom_vline(xintercept=0, color="grey", size = 0.2, alpha = 0.5) +
  geom_hline(yintercept=0, color="grey", size = 0.2, alpha = 0.5) +
  theme_classic() +
  coord_fixed() +
  theme(axis.title=element_text(colour="black", size=15),
    axis.text=element_text(colour="black", size=12),
    plot.title = element_blank(),
    legend.text=element_text(size=12),
    legend.title=element_text(size=14),
    legend.position="bottom") +
  scale_x_continuous(limits = c(-0.2, 0.8), breaks=seq(-0.2, 0.8, by = 0.2)) +
  scale_y_continuous(limits = c(-0.2, 0.8), breaks=seq(-0.2, 0.8, by = 0.2)) +
  scale_color_brewer(type = "seq", palette = "Blues", labels = c("non-attenuated", "lowly-attenuated", "highly-attenuated"), name = "Attenuation") +
  guides(color = guide_legend(title.position = "top", override.aes = list(size=2.5, shape=15))) +
  labs(x = "corr(CNV, RNA)", y = "corr(CNV, Protein)")
  # Marginal densities along x axis
xdens <- axis_canvas(correlation_plot2, axis = "x") +
  geom_density(data = cnv_rna_protein %>% filter(s_cnv_rna > -0.2 & s_cnv_prot > -0.2), aes(x = s_cnv_rna, fill = class),
  alpha = 0.7, size = 0.2) +
  ggpubr::fill_palette(palette = brewer.pal(3, "Blues"))
  # Marginal densities along y axis
  # Need to set coord_flip = TRUE, if you plan to use coord_flip()
ydens <- axis_canvas(correlation_plot2, axis = "y", coord_flip = TRUE)+
  geom_density(data = cnv_rna_protein %>% filter(s_cnv_rna > -0.2 & s_cnv_prot > -0.2), aes(x = s_cnv_prot, fill = class),
  alpha = 0.7, size = 0.2) +
  coord_flip() +
  ggpubr::fill_palette(palette = brewer.pal(3, "Blues"))
p1 <- insert_xaxis_grob(correlation_plot2, xdens, grid::unit(0.2, "null"), position = "top")
p2 <- insert_yaxis_grob(p1, ydens, grid::unit(0.2, "null"), position = "right")
ggsave(filename="protein_attenuation_scatterplot_spearman2.png", plot = p2, path = "./plots/protein_attenuation_tcga_ccle/", width=5, height=5)
unlink("protein_attenuation_scatterplot_spearman2.png")




save(list=ls(), file = "./r_workspaces/protein_attenuation_spearman.RData")
