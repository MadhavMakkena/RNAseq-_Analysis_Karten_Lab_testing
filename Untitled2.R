setwd("/Users/madhavmakkena/Desktop/RNAseq_Analysis/complete_data")

Up_Data <- read.csv("Up_Data.csv", row.names=1, header = T)
Down_Data <- read.csv("Down_Data.csv", row.names=1, header = T)
Complete_Data <- read.csv("Complete_Data.csv", row.names=1, header = T)
Nochange_Data <- read.csv("Nochange_Data.csv", row.names=1, header = T)


test <- read.csv(file="Cholesterol_Homeostasis.csv", header = FALSE)$V1

test3 <- c("NPC1", "NPC2")

setwd("/Users/madhavmakkena/Desktop/RNAseq_Analysis/Amigo_lists")

test2 <- Up_Data[row.names(Up_Data) %in% test, ]

test4 <- Complete_Data[row.names(Complete_Data) %in% test3, ]

Mitotic_Checkpoint.csv

test


Cholesterol_Homeostasis.csv