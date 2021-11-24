setwd("/Users/madhavmakkena/Desktop/R/Sets/atlas_GO_lists")
skn_list <- read.csv("skn_list.csv", header = F)$V1
cncr_list <- read.csv("cncr_list.csv", header = F)$V1
brain_list <- read.csv("brain_list.csv", header = F)$V1
RNApol_list <- read.csv("RNApol_list.csv", header = F)$V1
global_list <- read.csv("global_list.csv", header = F)$V1
multivesicular_list <- read.csv("multivesicular_list.csv", header = F)$V1
mito_list <- read.csv("mito_list.csv", header = F)$V1
mitotic_cell_cycle_list <- read.csv("mitotic_cell_cycle_list.csv", header = F)$V1
cell_pop_maintain_list <- read.csv("cell_pop_maintain_list.csv", header = F)$V1
growth_fctr_list <- read.csv("growth_fctr_list.csv", header = F)$V1

#
#
#
setwd("/Users/madhavmakkena/Desktop/R/Sets")
chol_homeo_list <- read.csv("Cholesterol_Homeostasis.csv", header = F)$V1
chol_homeo_list_gmt <- list(GeneSet1=(chol_homeo_list))

chol_transport_list <- read.csv("Cholesterol_Transport.csv", header = F)$V1
chol_transport_list_gmt <- list(GeneSet1=(chol_transport_list))

chol_list <- read.csv("chol_list.csv", header = F)$V1
chol_list_gmt <- list(GeneSet1=(chol_list))

Endosome_list <- read.csv("Endo_list.csv", header = F)$V1
Endosome_list_gmt <- list(GeneSet1=(Endosome_list))

# Aaron_list <- c(chol_homeo_list, chol_transport_list, Endosome_list)
# Aaron_list_gmt <- list(GeneSet1=(chol_homeo_list), GeneSet2=(chol_transport_list), GeneSet3=(Endosome_list))













