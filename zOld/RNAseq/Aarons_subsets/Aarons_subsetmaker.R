setwd("/Users/madhavmakkena/Downloads/RNAseq/Aarons_subsets")

#loading the data
Endosome <- read.csv("QuickGO-annotations-1632334000633-20210922.tsv", sep = "\t", header = TRUE)
Endosome <- subset(Endosome, select = c(SYMBOL))
Endosome <- distinct(Endosome)
write.csv(Endosome, "Endosome.csv", row.names = FALSE)

Cholesterol_Transport <- read.csv("QuickGO-annotations-1632417583680-20210923.tsv", sep = "\t", header = TRUE)
Cholesterol_Transport <- subset(Cholesterol_Transport, select = c(SYMBOL))
Cholesterol_Transport <- distinct(Cholesterol_Transport)
write.csv(Cholesterol_Transport, "Cholesterol_Transport.csv", row.names = FALSE)

Cholesterol_Homeostasis <- read.csv("QuickGO-annotations-1632417508178-20210923.tsv", sep = "\t", header = TRUE)
Cholesterol_Homeostasis <- subset(Cholesterol_Homeostasis, select = c(SYMBOL))
Cholesterol_Homeostasis <- distinct(Cholesterol_Homeostasis)
write.csv(Cholesterol_Homeostasis, "Cholesterol_Homeostasis.csv", row.names = FALSE)

#making the lists for graphs

#combined GMT list