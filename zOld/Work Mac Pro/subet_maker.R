# setwd("/Users/madhavmakkena/Downloads/RNAseq/Subsets/Raw_Goterm_lists")

setwd("/Users/madhavmakkena/Desktop/R/Sets")


# skn_list <- read.csv("tissue_category_rna_skin_Tissue.tsv", header = T, sep = '\t')
# skn_list <- subset(skn_list, select = c(Gene))
# skn_list <- distinct(skn_list)
# write.csv(skn_list, "skn_list.csv", row.names = F)
# rm(skn_list)
skn_list <- read.csv("skn_list.csv", header = F)$V1

# cncr_list <- read.csv("protein_class_COSMIC_Cancer_related.tsv", header = T, sep = '\t')
# cncr_list <- subset(cncr_list, select = c(Gene))
# cncr_list <- distinct(cncr_list)
# write.csv(cncr_list, "cncr_list.csv", row.names = F)
# rm(cncr_list)
cncr_list <- read.csv("cncr_list.csv", header = F)$V1

# brain_list <- read.csv("tissue_category_rna_brain_Tissue.tsv", header = T, sep = '\t')
# brain_list <- subset(brain_list, select = c(Gene))
# brain_list <- distinct(brain_list)
# write.csv(brain_list, "brain_list.csv", row.names = F)
# rm(brain_list)
brain_list <- read.csv("brain_list.csv", header = F)$V1

# RNApol_list <- read.csv("protein_class_RNA_RNApoly_All_tisues.tsv", header = T, sep = '\t')
# RNApol_list <- subset(RNApol_list, select = c(Gene))
# RNApol_list <- distinct(RNApol_list)
# write.csv(RNApol_list, "RNApol_list.csv", row.names = F)
# rm(RNApol_list)
RNApol_list <- read.csv("RNApol_list.csv", header = F)$V1

# global_list <- read.csv("tissue_category_rna_Any_Detected.tsv", header = T, sep = '\t')
# global_list <- subset(global_list, select = c(Gene))
# global_list <- distinct(global_list)
# write.csv(global_list, "global_list.csv", row.names = F)
# rm(global_list)
global_list <- read.csv("global_list.csv", header = F)$V1

chol_homeo_list <- read.csv("Cholesterol_Homeostasis.csv", header = F)$V1
chol_homeo_list_gmt <- list(GeneSet1=(chol_homeo_list))

chol_transport_list <- read.csv("Cholesterol_Transport.csv", header = F)$V1
chol_transport_list_gmt <- list(GeneSet1=(chol_transport_list))

chol_list <- c(chol_homeo_list, chol_transport_list)
chol_list_gmt <- list(GeneSet1=(chol_homeo_list), GeneSet2=(chol_transport_list))
write.csv(chol_list, "chol_list.csv", row.names = F)


Endosome_list <- read.csv("Endosome.csv", header = F)$V1
Endosome_list_gmt <- list(GeneSet1=(Endosome_list))

Aaron_list <- c(chol_homeo_list, chol_transport_list, Endosome_list)
Aaron_list_gmt <- list(GeneSet1=(chol_homeo_list), GeneSet2=(chol_transport_list), GeneSet3=(Endosome_list))

# multivesicular_list <- read.csv("QuickGO-annotations-multivesicular.tsv", header = T, sep = '\t')
# multivesicular_list <- subset(multivesicular_list, select = c(SYMBOL))
# multivesicular_list <- distinct(multivesicular_list)
# write.csv(multivesicular_list, "multivesicular_list.csv", row.names = F)
# rm(multivesicular_list)
multivesicular_list <- read.csv("multivesicular_list.csv", header = F)$V1

# mito_list <- read.csv("QuickGO-annotations-mitochondrion.tsv", header = T, sep = '\t')
# mito_list <- subset(mito_list, select = c(SYMBOL))
# mito_list <- distinct(mito_list)
# write.csv(mito_list, "mito_list.csv", row.names = F)
# rm(mito_list)
mito_list <- read.csv("mito_list.csv", header = F)$V1

# mitotic_cell_cycle_list <- read.csv("QuickGO-annotations-mitotic cell cycle.tsv", header = T, sep = '\t')
# mitotic_cell_cycle_list <- subset(mitotic_cell_cycle_list, select = c(SYMBOL))
# mitotic_cell_cycle_list <- distinct(mitotic_cell_cycle_list)
# write.csv(mitotic_cell_cycle_list, "mitotic_cell_cycle_list.csv", row.names = F)
# rm(mitotic_cell_cycle_list)
mitotic_cell_cycle_list <- read.csv("mitotic_cell_cycle_list.csv", header = F)$V1

# cell_pop_maintain_list <- read.csv("QuickGO-annotations-stem cell population maintenance.tsv", header = T, sep = '\t')
# cell_pop_maintain_list <- subset(cell_pop_maintain_list, select = c(SYMBOL))
# cell_pop_maintain_list <- distinct(cell_pop_maintain_list)
# write.csv(cell_pop_maintain_list, "cell_pop_maintain_list.csv", row.names = F)
# rm(cell_pop_maintain_list)
cell_pop_maintain_list <- read.csv("cell_pop_maintain_list.csv", header = F)$V1

# growth_fctr_list <- read.csv("QuickGO-annotations-growth factor activity.tsv", header = T, sep = '\t')
# growth_fctr_list <- subset(growth_fctr_list, select = c(SYMBOL))
# growth_fctr_list <- distinct(growth_fctr_list)
# write.csv(growth_fctr_list, "growth_fctr_list.csv", row.names = F)
# rm(growth_fctr_list)
growth_fctr_list <- read.csv("growth_fctr_list.csv", header = F)$V1

# chol_list <- read.csv("1_cholesterol_gene_sets_combined.csv", header = F)$V2
# chol_list_gmt <- list(GeneSet1=c(chol_list))
# 
# chol_list_c <- read.csv("1_cholesterol_gene_sets_combined copy.csv", header = F)$V2
# chol_list_gmt_c <- list(GeneSet1=c(chol_list))
# 
# volcano_gene_list <- c("NPC2", "NPC1", "ABCA1", "STARD3", "STARD4", "HMGCS1", "INHBA","SMO","SCARB2","SULT2B1","SOAT2")
# volcano_gene_lista <- c("NPC2", "NPC1", "ABCA1", "STARD3", "STARD4")
# volcano_gene_listb <- c("HMGCS1", "INHBA","SMO","SCARB2","SULT2B1","SOAT2")
# 
# 
# 
# summary(volcano_gene_list)
# geneSetz<-list(GeneSet1=c("NPC2","NPC1","ABCA1"))
# 
# geneSet<-list(GeneSet1=c("NPC2","NPC1","ABCA1"), GeneSet2=c("STARD3","STARD4","AACS","HMGCS1"), GeneSet3=c("INHBA","SMO","SCARB2","SULT2B1","SOAT2"))
# geneSets<- c("NPC2","NPC1","ABCA1", "STARD3","STARD4","AACS","HMGCS1", "INHBA","SMO","SCARB2","SULT2B1","SOAT2")
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# library("dplyr")
# library("stringr")
# #creating a list of dataframes and all the data from the file names
# comb_go_lists <- lapply(go_lists,function(x) {
#   read.table(file = x, 
#              sep = '\t', 
#              header = TRUE)
# })
# 
# 
# #combining the dataframe lists into a dataframe
# comb_go_lists <- bind_rows(comb_go_lists)
# 
# #Ordering based on the Go.Names and subsetting to have only go.names and symbols
# comb_go_lists <- comb_go_lists[order(comb_go_lists$GO.NAME),]
# comb_go_lists <- subset(comb_go_lists, select = c(SYMBOL,GO.NAME) )
# 
# # removing duplicates (only removed when common in both go.name and symbol)
# comb_go_lists <- distinct(comb_go_lists)
# 
# # creating a list of all the unique gonames and getting the length of this list
# go_name_list <- c(unique(comb_go_lists$GO.NAME))
# go_name_list_len <- length(go_name_list)
# 
# #creating a new vector used in the for loop to temporarily store names
# go_name <- vector(,1)
# 
# comb_go_lists_subset <- comb_go_lists
# 
# #creating a directory
# setwd("/Users/madhavmakkena/Downloads/RNAseq/Subsets")
# dir.create("all_subsets")
# 
# #subsetting and creating multiple csv files based on the go-names with the symbols as the data
# for (i in 1:go_name_list_len) {
#   go_name <- go_name_list[i]
#   comb_go_lists_subset <- subset(comb_go_lists, GO.NAME==go_name)
#   comb_go_lists_subset <- subset(comb_go_lists_subset, select = -c(GO.NAME) )
#   write.csv(comb_go_lists_subset, str_interp("all_subsets/${go_name}.csv"), row.names = FALSE)
# }
