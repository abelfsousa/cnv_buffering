# Abel Sousa 2018

# PhD project
# Study of Control of Protein Abundance using Multi-omics Data from Cancer Samples



# -- Manhattan pot



suppressMessages(library(tidyverse))
suppressMessages(library(ggrepel))
suppressMessages(library(RColorBrewer))

options(bitmapType = "cairo")





# import physical protein associations with CNV and RNA
# filtered by genomic co-localization
protein_associations_cnv_rna_physical <- read_tsv("./files/signf_protein_associations_cnv_mrna_filtered_close_controlling.txt") %>%
  mutate(pair = paste(controlling, controlled, sep="_")) %>%
  dplyr::select(pair, everything())




# import genome-wide protein associations with CNV and RNA
# filtered by genomic co-localization
protein_associations_cnv_mrna <- read_tsv("./files/signf_protein_associations_cnv_mrna_filtered_close_controlling_genomeWide.txt") %>%
  mutate(pair = paste(controlling, controlled, sep="_")) %>%
  dplyr::select(pair, everything()) %>%
  filter(!pair %in% protein_associations_cnv_rna_physical$pair)



# import protein associations with CNV
protein_associations_cnv <- readRDS("./files/protein_associations_tcga_ccle_cnv_reg_genome_wide.rds") %>%
  as.tibble() %>%
  mutate(pair = paste(controlling, controlled, sep="_")) %>%
  dplyr::select(pair, everything()) %>%
  mutate(in_RNA_model = pair %in% protein_associations_cnv_mrna$pair)


# import go enr categories
go_enr <- read_tsv("./files/protein_set_controlling_over_representation_GO.txt") %>%
  dplyr::select(Description, geneID) %>%
  mutate(geneID = str_split(geneID, "/")) %>%
  unnest()



# import all go terms
all_go_terms <- readRDS("./data/go/all_go_terms_symbols.rds")
all_go_terms <- tibble(term = names(all_go_terms), gene = all_go_terms) %>%
  unnest(gene)


degradation_go_terms <- all_go_terms  %>%
  filter(str_detect(term, "UBIQUITIN|LIGASE"))



# # volcano plot
# # plot CNV beta by FDR
# # color by FDR and mark as red the associations found with the RNA
# protein_associations_cnv_plot <- ggplot( data=protein_associations_cnv, mapping=aes(x=beta_cnv, y=-log10(fdr), color=in_RNA_model, fill = fdr_label) ) +
#   geom_point(pch = 21) +
#   theme_classic() +
#   #geom_text_repel(data = head(protein_associations_cnv[protein_associations_cnv$in_RNA_model == "3",], 20), aes(x=beta_cnv, y=-log10(fdr), label = paste(controlling, controlled, sep="~")), size = 2.5, box.padding = unit(0.3, "lines"), point.padding = unit(0.3, "lines"), colour="black" ) +
#   #scale_fill_manual(values=rev(brewer.pal(4, "Blues")[3:4])) +
#   scale_fill_manual(values=c("#2171B5", "#BDD7E7"), name = "", labels = c("FDR <= 0.05", "FDR > 0.05")) +
#   scale_colour_manual(values=c("#BDD7E7", "#2171B5", "red"), guide=F) +
#   scale_x_continuous(limits = c(-0.8, 0.8)) +
#   geom_line(aes(x=0),colour="black", linetype=2, size = 0.3) +
#   geom_line(aes(y=-log10(0.05)),colour="black", linetype=2, size = 0.3) +
#   theme(axis.title.y=element_text(colour="black", size=15),
#     axis.title.x=element_text(colour="black", size=15),
#     axis.text.y=element_text(colour="black", size=13),
#     axis.text.x=element_text(colour="black", size=13),
#     plot.title = element_text(size=15, hjust = 0.5),
#     legend.text=element_text(size=12),
#     legend.title=element_text(size=15)) +
#   labs(x = "CNV Beta", y = "FDR (-log10)", title="Protein associations")
# ggsave(filename="protein_associations_tcga_ccle_cnv_rna_reg.png", plot=protein_associations_cnv_plot, path="./plots/protein_associations_genome_wide/")
# unlink("protein_associations_tcga_ccle_cnv_rna_reg.png")



# manhattan plot
# plot controlling protein chromosome by protein-association FDR
# color by FDR and mark as red the associations found with the RNA

# load genomic coordinates
load("/nfs/research1/beltrao/memon/databases/protein/ensembl/ensembl_v73/ensembl_gtf_v73_gene_anno_grange.Rdata")
ensembl_gtf_v73_gene_anno_grange <- ensembl_gtf_v73_gene_anno_grange %>%
  as.data.frame() %>%
  as.tibble() %>%
  dplyr::select(chr = seqnames, gene_name, gene_biotype) %>%
  mutate(chr = fct_drop(chr)) %>%
  filter(gene_biotype == "protein_coding") %>%
  filter(!(gene_name == "CKS1B" & chr == "5")) %>%
  filter(!(gene_name == "LSP1" & chr == "13")) %>%
  distinct()

protein_associations_manht <- protein_associations_cnv %>%
  dplyr::select(controlling, controlled, beta_cnv, fdr, in_RNA_model) %>%
  inner_join(ensembl_gtf_v73_gene_anno_grange[, c("chr", "gene_name")], by = c("controlling" = "gene_name")) %>%
  mutate(chr = as.factor(chr)) %>%
  mutate(chr = fct_drop(chr)) %>%
  mutate(chr = fct_relevel(chr, as.character(1:22))) %>%
  arrange(chr, controlling) %>%
  mutate(xpos = 1:nrow(.))

protein_associations_manht2 <- protein_associations_manht %>%
  group_by(chr) %>%
  summarize(center=(max(xpos) + min(xpos)) / 2)

protein_associations_manht_plot <-  ggplot(protein_associations_manht, mapping=aes(x=xpos, y=-log10(fdr))) +
    geom_point(mapping=aes(colour=chr), alpha=0.8, size=1.3) +
    #scale_color_manual(values = rep(c("#deebf7", "#9ecae1"), 22 ), guide = F) +
    scale_color_manual(values = rep(c("#e5f5f9", "#99d8c9"), 22 ), guide = F) +
    geom_point(protein_associations_manht %>% filter(controlling %in% intersect(go_enr[go_enr$Description == "catalytic complex", "geneID"]$geneID, degradation_go_terms$gene) & in_RNA_model) %>% arrange(fdr), mapping=aes(x=xpos, y=-log10(fdr)), colour="red", alpha=0.5, size=1.3) +
    geom_point(protein_associations_manht %>% filter(in_RNA_model) %>% arrange(fdr) %>% head(5), mapping=aes(x=xpos, y=-log10(fdr)), colour="black", alpha=0.5, size=1.3) +
    theme_classic() +
    geom_text_repel(
      data = bind_rows(
        protein_associations_manht %>% filter(controlling %in% intersect(go_enr[go_enr$Description == "catalytic complex", "geneID"]$geneID, degradation_go_terms$gene) & in_RNA_model) %>% arrange(fdr),
        protein_associations_manht %>% filter(in_RNA_model) %>% arrange(fdr) %>% head(5)),
      aes(x=xpos, y=-log10(fdr), label = paste(controlling, controlled, sep="~")), size = 3, box.padding = unit(0.3, "lines"), point.padding = unit(0.3, "lines"), colour="black" ) +
    scale_x_continuous(label = protein_associations_manht2$chr, breaks= protein_associations_manht2$center) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 15) ) +
    geom_line(aes(y=-log10(0.05)),colour="black", linetype=2, size = 0.5) +
    theme(axis.title=element_text(colour="black", size=15),
      axis.text.x=element_text(colour="black", size=10),
      axis.text.y=element_text(colour="black", size=12),
      plot.title = element_blank()) +
    labs(x = "Chromosome", y = "CNV model FDR (-log10)")
ggsave(filename="genome_wide_associations_manht.png", plot=protein_associations_manht_plot, path="./plots/protein_associations_genome_wide/")
unlink("genome_wide_associations_manht.png")
