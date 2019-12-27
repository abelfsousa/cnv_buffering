# Abel Sousa 2018

# PhD project
# Study of Control of Protein Abundance using Multi-omics Data from Cancer Samples

# -- Comparison of protein associations using normal samples

# -- Correlations of protein/RNA expression across normal tissues




suppressMessages(library(tidyverse))
suppressMessages(library(ggpubr))
suppressMessages(library(gplots))

options(bitmapType = "cairo")



# load R workspace
load(file="./r_workspaces/normal_samples_gtex_hpm.RData")



# top controlling-controlled pairs with highest protein correlation
top_pairs <- gtex_hpm_control_status_cor$data %>%
	filter(signf == "highly-significant", correlation == "cor_Pcontrolling_Pcontrolled") %>%
	arrange(desc(estimate)) %>%
	head(50) %>%
	dplyr::select(controlling, controlled)

#gtex_hpm_control_status_cor$data %>%
#spread(key="correlation", value="estimate") %>%
#group_by(controlling, controlled) %>%
#filter(cor_Pcontrolling_Pcontrolled > 0.4, cor_RNAcontrolled_Pcontrolled < 0.2) %>%
#filter(signf == "highly-significant") %>%
#arrange(desc(cor_Pcontrolling_Pcontrolled)) %>%
#as.data.frame



pdf(file="./plots/protein_associations_normal/examples_controlling_controlled_heatmap.pdf", w=10, h=4)
par(oma=c(2,0,0,3))


for(i in 1:nrow(top_pairs)){

example <- p_associations %>%
	filter(controlling == as.character(top_pairs[i,1]), controlled == as.character(top_pairs[i,2])) %>%
	dplyr::arrange(Protein_controlling) %>%
	dplyr::select(-c(controlling, controlled, signf)) %>%
	as.data.frame()

rownames(example) <- example$Tissue
example <- example[, -c(1)]


heatmap.2( x=t(example),
	Rowv=T,
	Colv=T,
	dendrogram="both",
	scale = "none",
	trace = "none",
	col = greenred(100),
	main = paste (as.character(top_pairs[i,1]), as.character(top_pairs[i,2]), sep = " vs "),
	xlab = "",
	ylab = "",
	key.title=NA,
	key.xlab=NA,
	key.ylab=NA,
	cexCol = 1,
	cexRow = 1,
	density.info="none")

}

dev.off()




#pdf(file="./plots/protein_associations_normal/TRAPPC8_controlling_TRAPPC11_controlled_example.pdf", w=12, h=4)
#par(oma=c(0,0,0,11))
png(file="./plots/protein_associations_normal/TRAPPC8_controlling_TRAPPC11_controlled_example.png", width = 1000, height = 350)
par(oma=c(0,0,0,11))

example2 <- p_associations %>%
	filter(controlling == "TRAPPC8", controlled == "TRAPPC11") %>%
	dplyr::arrange(Protein_controlling) %>%
	dplyr::select(-c(controlling, controlled, signf)) %>%
	as.data.frame()

rownames(example2) <- example2$Tissue
example2 <- example2[, -c(1)]

heatmap.2( x=t(example2),
	Rowv=T,
	Colv=T,
	dendrogram="both",
	scale = "none",
	trace = "none",
	col = greenred(100),
	main = "TRAPPC8 and TRAPPC11",
	xlab = "",
	ylab = "",
	key=T,
	key.title="RNA/Protein expression",
	keysize=1.5,
	key.xlab=NA,
	key.ylab=NA,
	key.par=list(cex = 1, cex.main = 1),
	cexRow = 2,
	density.info="none",
	labCol = NA,
	labRow = c("Protein TRAPPC8", "RNA TRAPPC8", "Protein TRAPPC11", "RNA TRAPPC11"))

dev.off()




#pdf(file="./plots/protein_associations_normal/ARCN1_controlling_COPA_controlled_example.pdf", w=12, h=4)
#par(oma=c(0,0,0,11))
png(file="./plots/protein_associations_normal/ARCN1_controlling_COPA_controlled_example.png", width = 1000, height = 350)
par(oma=c(0,0,0,11))

example2 <- p_associations %>%
	filter(controlling == "ARCN1", controlled == "COPA") %>%
	dplyr::arrange(Protein_controlling) %>%
	dplyr::select(-c(controlling, controlled, signf)) %>%
	as.data.frame()

rownames(example2) <- example2$Tissue
example2 <- example2[, -c(1)]

heatmap.2( x=t(example2),
	Rowv=T,
	Colv=T,
	dendrogram="both",
	scale = "none",
	trace = "none",
	col = greenred(100),
	main = "ARCN1 and COPA",
	xlab = "",
	ylab = "",
	key=T,
	key.title="RNA/Protein expression",
	keysize=1.5,
	key.xlab=NA,
	key.ylab=NA,
	key.par=list(cex = 1, cex.main = 1),
	cexRow = 2,
	density.info="none",
	labCol = NA,
	labRow = c("Protein ARCN1", "RNA ARCN1", "Protein COPA", "RNA COPA"))

dev.off()










example2 <- p_associations %>%
	filter(controlling == "TRAPPC8", controlled == "TRAPPC11") %>%
	dplyr::select(-c(controlling, controlled, signf)) %>%
	gather(key="data", value="expression", -Tissue) %>%
	ggplot(aes(x = Tissue, y = data, fill = expression)) +
	geom_tile() +
	theme_classic() +
	scale_x_discrete(expand = c(0, 0)) +
	scale_y_discrete(expand = c(0, 0)) +
	scale_fill_gradient(low = "black", high = "blue") +
	theme(axis.title = element_text(colour="black", size=15),
		axis.text.y = element_text(colour="black", size=13),
		axis.text.x = element_text(colour="black", size=13, angle=45, hjust=1),
		axis.ticks = element_blank(),
		axis.line = element_blank()) +
	labs(x = "", y = "")
ggsave(filename="TRAPPC8_controlling_TRAPPC11_controlled_example2.png", plot=example2, width = 14, height = 4, path = "./plots/protein_associations_normal/")
unlink("TRAPPC8_controlling_TRAPPC11_controlled_example2.png")
