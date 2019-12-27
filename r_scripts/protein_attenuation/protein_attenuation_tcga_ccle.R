# Abel Sousa 2018

# PhD project
# Study of Control of Protein Abundance using Multi-omics Data from Cancer Samples


# -- Protein attenuation



suppressMessages(library(mclust))
suppressMessages(library(tidyverse))
suppressMessages(library(RColorBrewer))
suppressMessages(library(gridExtra))
suppressMessages(library(ggpubr))
suppressMessages(library(cowplot))
options(bitmapType = "cairo")
set.seed(123)




# load datasets


# samples with CNV, RNA and Protein measurements
common_samples <- read_tsv(file = "./files/cptac_cellLines_samples_prot_rna_cnv.txt", col_names = FALSE)


proteomics <- read_tsv("./files/proteomicsQ_cptac_cellLines_corrected.txt") %>%
    gather(key="sample", value="prot_log2FC", -gene)

rna <- read_tsv("./files/rna_tcga_cellLines_corrected.txt") %>%
    gather(key="sample", value="rna_log2CPM", -gene)


cnv <- read_tsv("./files/cnv_tcga_cellLines.txt") %>%
    gather(key="sample", value="cnv_gistic2", -gene) %>%
    filter(sample %in% common_samples$X1)





# merge datasets by gene and sample
merged_data <- inner_join(proteomics, rna, by=c("gene", "sample")) %>%
    inner_join(cnv, by=c("gene", "sample"))
write.table(merged_data, "./files/protein_attenuation_cnv_rna_protein.txt", sep="\t", quote=F, row.names=F)



# compute correlation between CNV and RNA by gene across samples
cnv_rna <- merged_data %>%
    group_by(gene) %>%
    do(broom::tidy(cor.test(.$cnv_gistic2, .$rna_log2CPM, method = "pearson"))) %>%
    ungroup()



# compute correlation between CNV and PROTEIN by gene across samples
cnv_protein <- merged_data %>%
    group_by(gene) %>%
    do(broom::tidy(cor.test(.$cnv_gistic2, .$prot_log2FC, method = "pearson"))) %>%
    ungroup()


# compute correlation between PROTEIN and RNA by gene across samples
protein_rna <- merged_data %>%
    group_by(gene) %>%
    do(broom::tidy(cor.test(.$prot_log2FC, .$rna_log2CPM, method = "pearson"))) %>%
    ungroup()


# merge correlations and p-values
# compute attenuation potential: cor(CNV, RNA) - cor(CNV, Protein)
# gaussian mixture model to select proteins with high attenuation levels
cnv_rna_protein <- inner_join(cnv_rna, cnv_protein, by="gene") %>%
    dplyr::select(gene, estimate.x, p.value.x, estimate.y, p.value.y) %>%
    dplyr::rename(r_cnv_rna=estimate.x, p_cnv_rna=p.value.x, r_cnv_prot=estimate.y, p_cnv_prot=p.value.y) %>%
    mutate(attenuation = r_cnv_rna - r_cnv_prot)


gmm <- mclust::Mclust(cnv_rna_protein$attenuation, G=4)


cnv_rna_protein <- cnv_rna_protein %>%
    mutate(class = gmm$classification) %>%
    mutate(class = if_else(class==1 | class==2, "non-attenuated", if_else(class==3, "low-attenuated-protein", "high-attenuated-protein" ))) %>%
    arrange(desc(attenuation)) %>%
    mutate_at(vars(contains("class")), as.factor) %>%
    mutate_if(is.factor, fct_infreq)
write.table(cnv_rna_protein, "./files/protein_attenuation.txt", sep="\t", quote=F, row.names=F)



# some stats
attenuation_stats <- cnv_rna_protein %>%
    group_by(class) %>%
    summarise(counts = n(), perc = n()/nrow(cnv_rna_protein), mean = mean(attenuation))




correlation_plot1 <- ggplot(cnv_rna_protein %>% filter(r_cnv_rna > -0.2 & r_cnv_prot > -0.2)) +
  geom_point(aes(x = r_cnv_rna, y = r_cnv_prot), size= 1, color="grey") +
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
    legend.position="none",
    panel.border = element_rect(colour = "black", fill=NA, size=1)) +
    scale_x_continuous(limits = c(-0.2, 0.8), breaks=seq(-0.2, 0.8, by = 0.2)) +
    scale_y_continuous(limits = c(-0.2, 0.8), breaks=seq(-0.2, 0.8, by = 0.2)) +
  labs(x = "corr(CNV, RNA)", y = "corr(CNV, Protein)")
ggsave(filename="protein_attenuation_tcga_ccle_scatterplot.png", plot = correlation_plot1, path = "./plots/protein_attenuation_tcga_ccle/", width=5, height=5)
unlink("protein_attenuation_tcga_ccle_scatterplot.png")



correlation_plot2 <- ggplot(cnv_rna_protein %>% filter(r_cnv_rna > -0.2 & r_cnv_prot > -0.2)) +
  geom_point(aes(x = r_cnv_rna, y = r_cnv_prot, color = class), size= 1) +
  geom_density2d(mapping=aes(x=r_cnv_rna, y=r_cnv_prot, group=class), color="grey", size = 0.2) +
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
  guides(color = guide_legend(title.position = "top", override.aes = list(size=3, shape=19))) +
  labs(x = "corr(CNV, RNA)", y = "corr(CNV, Protein)")
# Marginal densities along x axis
xdens <- axis_canvas(correlation_plot2, axis = "x") +
  geom_density(data = cnv_rna_protein %>% filter(r_cnv_rna > -0.2 & r_cnv_prot > -0.2), aes(x = r_cnv_rna, fill = class),
  alpha = 0.95, size = 0.6) +
  ggpubr::fill_palette(palette = brewer.pal(3, "Blues"))
# Marginal densities along y axis
# Need to set coord_flip = TRUE, if you plan to use coord_flip()
ydens <- axis_canvas(correlation_plot2, axis = "y", coord_flip = TRUE)+
  geom_density(data = cnv_rna_protein %>% filter(r_cnv_rna > -0.2 & r_cnv_prot > -0.2), aes(x = r_cnv_prot, fill = class),
  alpha = 0.95, size = 0.6) +
  coord_flip() +
  ggpubr::fill_palette(palette = brewer.pal(3, "Blues"))
p1 <- insert_xaxis_grob(correlation_plot2, xdens, grid::unit(0.2, "null"), position = "top")
p2 <- insert_yaxis_grob(p1, ydens, grid::unit(0.2, "null"), position = "right")
ggsave(filename="protein_attenuation_tcga_ccle_scatterplot2.png", plot = p2, path = "./plots/protein_attenuation_tcga_ccle/", width=5, height=5)
ggsave(filename="protein_attenuation_tcga_ccle_scatterplot2.pdf", plot = p2, path = "./plots/protein_attenuation_tcga_ccle/", width=5, height=5)
unlink("protein_attenuation_tcga_ccle_scatterplot2.png")
unlink("protein_attenuation_tcga_ccle_scatterplot2.pdf")



# plot scatterplot with CNV ~ RNA against CNV ~ Protein
correlation_plot3 <- ggplot(data=cnv_rna_protein %>% mutate(sign = if_else(attenuation > 0, "positive", "negative"))) +
  geom_point(mapping=aes(x=r_cnv_rna, y=r_cnv_prot, colour=sign)) +
  geom_density2d(mapping=aes(x=r_cnv_rna, y=r_cnv_prot, group=class), color="grey", size = 0.5) +
  geom_abline(slope=1, intercept=0, linetype=2, color="black") +
  theme_classic() +
  coord_fixed() +
  theme(axis.title.y=element_text(colour="black", size=12),
    axis.title.x=element_text(colour="black", size=12),
    axis.text.y=element_text(colour="black", size=10),
    axis.text.x=element_text(colour="black", size=10),
    plot.title = element_blank(),
    legend.text=element_text(size=10),
    legend.title=element_text(size=12),
    legend.position="right",
    panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
  scale_x_continuous(limits = c(-0.4, 0.8), breaks=seq(-0.4, 0.8, by = 0.2)) +
  scale_y_continuous(limits = c(-0.4, 0.8), breaks=seq(-0.4, 0.8, by = 0.2)) +
  labs(x = "corr(CNV, RNA)", y = "corr(CNV, Protein)", title = "")
ggsave(filename="protein_attenuation_tcga_ccle_scatterplot3.png", plot = correlation_plot3, path = "./plots/protein_attenuation_tcga_ccle/", width=5, height=5)
unlink("protein_attenuation_tcga_ccle_scatterplot3.png")



# plot histogram/density lines
histogram1 <- ggplot(data=cnv_rna_protein) +
    #geom_histogram(mapping=aes(x=attenuation, y=..density.., fill=class), color="black", alpha=0.6, bins=20) +
    geom_density(mapping=aes(x=attenuation, fill=class, colour=class), size=0.8, alpha=0.6) +
    #geom_vline(data=attenuation_stats, mapping=aes(xintercept=mean, group=class), color="black", linetype="dashed", size=0.5) +
    theme_classic() +
    theme(axis.title.y=element_text(colour="black", size=15),
      axis.title.x=element_text(colour="black", size=15),
      axis.text.y=element_text(colour="black", size=12),
      axis.text.x=element_text(colour="black", size=12),
      plot.title = element_blank(),
      legend.text=element_text(colour="black", size=10),
      legend.title=element_text(colour="black", size=12)) +
    scale_fill_brewer(type = "seq", palette = "Blues", labels = c("non-attenuated", "lowly-attenuated", "highly-attenuated"), name = "Attenuation class") +
    scale_color_brewer(type = "seq", palette = "Blues", labels = c("non-attenuated", "lowly-attenuated", "highly-attenuated"), name = "Attenuation class") +
    labs(x = "Attenuation potential", y = "Density")
ggsave(filename="protein_attenuation_tcga_ccle_histogram.png", plot = histogram1, height=3, width=5, path = "./plots/protein_attenuation_tcga_ccle/")
unlink("protein_attenuation_tcga_ccle_histogram.png")


histogram2 <- ggplot(data=cnv_rna_protein %>% filter(attenuation < 0)) +
    geom_density(mapping=aes(x=attenuation, fill=class, colour=class), size=0.8, alpha=0.6) +
    theme_classic() +
    theme_classic() +
    theme(axis.title.y=element_text(colour="black", size=15),
      axis.title.x=element_text(colour="black", size=15),
      axis.text.y=element_text(colour="black", size=12),
      axis.text.x=element_text(colour="black", size=12),
      plot.title = element_blank(),
      legend.text=element_text(colour="black", size=10),
      legend.title=element_text(colour="black", size=12)) +
    scale_fill_brewer(type = "seq", palette = "Blues", labels = c("non-attenuated", "lowly-attenuated", "highly-attenuated"), name = "Attenuation class") +
    scale_color_brewer(type = "seq", palette = "Blues", labels = c("non-attenuated", "lowly-attenuated", "highly-attenuated"), name = "Attenuation class") +
    labs(x = "Attenuation potential", y = "Density")
ggsave(filename="protein_attenuation_tcga_ccle_histogram2.png", plot = histogram2, path = "./plots/protein_attenuation_tcga_ccle/")
unlink("protein_attenuation_tcga_ccle_histogram2.png")





# plot boxplot
boxplot1 <- ggplot(data=cnv_rna_protein) +
    geom_boxplot(mapping=aes(x=class, y=attenuation, fill=class)) +
    theme_classic() +
    theme(axis.title.y=element_text(colour="black", size=13), axis.title.x=element_text(colour="black", size=13)) +
    theme(axis.text.y=element_text(colour="black", size=10), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
    scale_y_continuous(limits = c(-1, 1)) +
    labs(x = "Classification", y = "Attenuation potential: [CNV ~ RNA] - [CNV ~ Protein]", title = "")
ggsave(filename="protein_attenuation_tcga_ccle_boxplot.png", plot = boxplot1, path = "./plots/protein_attenuation_tcga_ccle/")
unlink("protein_attenuation_tcga_ccle_boxplot.png")




# plot barplot
barplot1 <- ggplot(data=attenuation_stats) +
    geom_bar(mapping=aes(x=class, y=counts, fill=class), stat="identity") +
    theme_classic() +
    coord_flip() +
    theme(axis.title=element_text(colour="black", size=12),
          axis.text.x=element_text(colour="black", size=10),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          legend.text=element_blank(),
          legend.title=element_blank(),
          plot.title=element_blank(),
          legend.position = "none") +
    scale_fill_brewer(type = "seq", palette = "Blues", labels = c("non-attenuated", "lowly-attenuated", "highly-attenuated"), name = "", guide=F) +
    scale_y_continuous(limits = c(0, 5000), breaks=seq(0, 5000, by = 1000), name = "Number of proteins") +
    scale_x_discrete(name = "")
ggsave(filename="protein_attenuation_tcga_ccle_barplot.png", height=2, width=4, plot = barplot1, path = "./plots/protein_attenuation_tcga_ccle/")
unlink("protein_attenuation_tcga_ccle_barplot.png")


# plot barplot
barplot1 <- ggplot(data=attenuation_stats) +
    geom_bar(mapping=aes(x=class, y=counts, fill=class), stat="identity") +
    theme_classic() +
    theme(
      axis.title=element_text(colour="black", size=15),
      axis.text.x=element_blank(),
      axis.text.y=element_text(colour="black", size=12),
      axis.ticks.x=element_blank(),
      plot.title=element_blank()) +
    scale_fill_brewer(type = "seq", palette = "Blues", labels = c("non-attenuated", "lowly-attenuated", "highly-attenuated"), name = "", guide=F) +
    scale_y_continuous(limits = c(0, 5000), breaks=seq(0, 5000, by = 1000), name = "Number of proteins") +
    scale_x_discrete(name = "Attenuation")
ggsave(filename="protein_attenuation_tcga_ccle_barplot.png", height=5, width=3, plot = barplot1, path = "./plots/protein_attenuation_tcga_ccle/")
ggsave(filename="protein_attenuation_tcga_ccle_barplot.pdf", height=5, width=3, plot = barplot1, path = "./plots/protein_attenuation_tcga_ccle/")
unlink("protein_attenuation_tcga_ccle_barplot.png")
unlink("protein_attenuation_tcga_ccle_barplot.pdf")



# plot barplot of cnv distribution by attenuation state
barplot2 <- merged_data %>%
  dplyr::select(gene, sample, cnv_gistic2) %>%
  inner_join(cnv_rna_protein[, c("gene", "class")], by="gene") %>%
  group_by(gene, class, cnv_gistic2) %>%
  summarise(count = n()) %>%
  group_by(class, cnv_gistic2) %>%
  summarise(count = sum(count))

barplot2 <- barplot2 %>%
  group_by(class) %>%
  mutate(count2 = sum(count)) %>%
  ungroup() %>%
  mutate(count3 = count/count2) %>%
  mutate(cnv_gistic2 = as.factor(as.character(cnv_gistic2))) %>%
  mutate(cnv_gistic2 = fct_relevel(cnv_gistic2, c("-2", "-1", "0", "1", "2"))) %>%
  ggplot() +
  geom_bar(mapping=aes(x=class, y=count3, fill=cnv_gistic2), stat="identity", position="dodge") +
  theme_classic() +
  theme(
    axis.title=element_text(colour="black", size=15),
    axis.text=element_text(colour="black", size=12),
    plot.title=element_blank()) +
  scale_fill_brewer(type = "seq", palette = "YlOrBr") +
  #scale_y_continuous(limits = c(0, 5000), breaks=seq(0, 5000, by = 1000), name = "Number of proteins") +
  scale_x_discrete(labels = c("non-attenuated", "lowly-attenuated", "highly-attenuated"), name = "Attenuation") +
  scale_y_continuous(name = "Percentage of measurements")
ggsave(filename="protein_attenuation_tcga_ccle_barplot_cnv_dist.png", plot = barplot2, path = "./plots/protein_attenuation_tcga_ccle/")
unlink("protein_attenuation_tcga_ccle_barplot_cnv_dist.png")


save(list=ls(), file = "./r_workspaces/protein_attenuation_tcga_ccle.RData")
