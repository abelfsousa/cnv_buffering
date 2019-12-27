# Abel Sousa 2018

# PhD project
# Study of Control of Protein Abundance using Multi-omics Data from Cancer Samples


# -- Compendium of protein-protein interactions from CORUM and BioGRID databases



suppressMessages(library(tidyverse))
options(bitmapType = "cairo")




# load biogrid data
biogrid <- read_tsv("./data/biogrid/BIOGRID-ORGANISM-3/BIOGRID-ORGANISM-Homo_sapiens-3.4.157.tab2.txt")
biogrid_type <- read_csv("./data/biogrid/physical_experimental_systems.csv")


# select only physical interactions and experiments defined as stable or transient
# remove duplicated interactor.A/interactor.B/experiment.type
# remove protein homodimers
biogrid_pairs <- biogrid %>%
	filter(`Experimental System` %in% biogrid_type[biogrid_type$Group %in% c("stable", "transient"), "Experiment"]$Experiment) %>%
	filter(`Organism Interactor A` == 9606, `Organism Interactor B` == 9606) %>%
	dplyr::select(`Official Symbol Interactor A`, `Official Symbol Interactor B`, `Experimental System`) %>%
	rename(interactor.A = `Official Symbol Interactor A`, interactor.B = `Official Symbol Interactor B`, experiment.type = `Experimental System`) %>%
	#mutate(interactor.A = toupper(interactor.A), interactor.B = toupper(interactor.B)) %>%
	distinct() %>%
	filter(interactor.A != interactor.B)
nrow(biogrid_pairs)
#318,716


# duplicate matrix with interactions inverted
biogrid_pairs <- biogrid_pairs %>%
	rbind(setNames(biogrid_pairs[, c(2,1,3)], colnames(biogrid_pairs))) %>%
	distinct() %>%
	filter(interactor.A != interactor.B)
write.table(biogrid_pairs, file = "./files/biogrid_pairs.txt", sep = "\t", quote = F, row.names=F)


cat("biogrid interactions:", nrow(unique(biogrid_pairs[,c(1:2)])), "\n")
#524,148





# load corum data
corum <- read_tsv("./data/corum/coreComplexes_29_05_2018.txt")
corum_pairs <- corum %>%
	filter(Organism == "Human") %>%
	dplyr::select(ComplexName, `subunits(Gene name)`) %>%
	rename(GeneName = `subunits(Gene name)`) %>%
	#mutate(GeneName = toupper(GeneName)) %>%
	group_by(ComplexName) %>%
	summarise(GeneName = paste(GeneName, collapse=";")) %>%
	ungroup() %>%
	mutate(GeneName = str_split(GeneName, ";\\s?")) %>%
	unnest() %>%
	filter(GeneName != "") %>%
	distinct() %>%
	group_by(ComplexName) %>%
	mutate(pair = list(tidyr::crossing(GeneName, GeneName))) %>%
	ungroup() %>%
	unnest() %>%
	dplyr::select("GeneName1", "GeneName2") %>%
	distinct() %>%
	filter(GeneName1 != GeneName2) %>%
	rename(interactor.A = GeneName1, interactor.B = GeneName2) %>%
	mutate(experiment.type="corum")
write.table(corum_pairs, file = "./files/corum_protein_pairs.txt", sep = "\t", quote = F, row.names=F)


cat("corum interactions:", nrow(corum_pairs), "\n")
#74712





## import ER membrane proteins
ER <- read_tsv("./data/protein_complexes/subunits_of_multisubunit_complexes_ER_INM.txt") %>%
	mutate(gene_name = toupper(gene_name))

# all interactions within protein complexes
ER_pairs <- ER %>%
	group_by(complex_name) %>%
	mutate(pair = list(tidyr::crossing(gene_name, gene_name))) %>%
	ungroup() %>%
	unnest() %>%
	dplyr::select("gene_name1", "gene_name2") %>%
	distinct() %>%
	filter(gene_name1 != gene_name2) %>%
	dplyr::rename(interactor.A = gene_name1, interactor.B = gene_name2) %>%
	mutate(experiment.type="curated_membrane_complexes")
write.table(ER_pairs, file = "./files/ER_protein_pairs.txt", sep = "\t", quote = F, row.names=F)


cat("ER membrane proteins interactions:", nrow(ER_pairs), "\n")
#1196







# bind corum, biogrid and ER interactions
biogrid_corum_ER_pairs <- bind_rows(biogrid_pairs, corum_pairs, ER_pairs)
write.table(biogrid_corum_ER_pairs, file = "./files/biogrid_corum_ER_pairs.txt", sep = "\t", quote = F, row.names=F)



# collapse experiment types
biogrid_corum_ER_pairs2 <- biogrid_corum_ER_pairs %>%
	group_by(interactor.A, interactor.B) %>%
	summarise(experiment.type = paste(experiment.type, collapse = "|"))
write.table(biogrid_corum_ER_pairs2, file = "./files/biogrid_corum_ER_pairs2.txt", sep = "\t", quote = F, row.names=F)



cat("BioGRID + CORUM + ER membrane interactions:", nrow(biogrid_corum_ER_pairs2), "\n")
#572,856






save(list=ls(), file="./r_workspaces/protein_pairs_biogrid_corum.RData")
