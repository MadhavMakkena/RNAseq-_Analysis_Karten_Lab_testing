setwd("/Users/madhavmakkena/Desktop/RNAseq_Analysis/Amigo_lists/Raw_lists")


Mitochondrion <- read.csv(file="GO_0005739_mitochondrion.csv", header=TRUE)
Mitochondrion <- distinct(Mitochondrion)
write.csv(Mitochondrion, "Mitochondrion.csv", row.names=F)
rm(Mitochondrion)
