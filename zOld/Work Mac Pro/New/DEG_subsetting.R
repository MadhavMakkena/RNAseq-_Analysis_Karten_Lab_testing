setwd("/Users/madhavmakkena/Desktop/R/New")

setwd("/Users/madhavmakkena/Desktop/R/New/Gene_Sets")
# Endosome <- read.table(file = 'QuickGO-annotations-Endosome.tsv', sep = '\t', header = TRUE)
# Endosome <- Endosome["SYMBOL"]
# Endosome <- distinct(Endosome)
# write.csv(Endosome, "Endosome.csv", row.names = F)
Endosome <- read.csv('Endosome.csv', header = F)$V1

# Cholesterol_Homeostasis <- read.table(file = 'QuickGO-annotations-Cholesterol_Homeostasis.tsv', sep = '\t', header = TRUE)
# Cholesterol_Homeostasis <- Cholesterol_Homeostasis["SYMBOL"]
# Cholesterol_Homeostasis <- distinct(Cholesterol_Homeostasis)
# write.csv(Cholesterol_Homeostasis, "Cholesterol_Homeostasis.csv", row.names = F)
Cholesterol_Homeostasis <- read.csv('Cholesterol_Homeostasis.csv', header = F)$V1

# Cholesterol_Transport <- read.table(file = 'QuickGO-annotations-Cholesterol_Transport.tsv', sep = '\t', header = TRUE)
# Cholesterol_Transport <- Cholesterol_Transport["SYMBOL"]
# Cholesterol_Transport <- distinct(Cholesterol_Transport)
# write.csv(Cholesterol_Transport, "Cholesterol_Transport.csv", row.names = F)
Cholesterol_Transport <- read.csv('Cholesterol_Transport.csv', header = F)$V1

# Cholesterol_Total <- c(Cholesterol_Homeostasis, Cholesterol_Transport)
# Cholesterol_Total <- as.data.frame(Cholesterol_Total)
# Cholesterol_Total <- Cholesterol_Total[Cholesterol_Total != "SYMBOL",]
# Cholesterol_Total <- as.data.frame(Cholesterol_Total)
# colnames(Cholesterol_Total) <- c("SYMBOL")
# Cholesterol_Total <- distinct(Cholesterol_Total)
# write.csv(Cholesterol_Total, "Cholesterol_Total.csv", row.names = F)
Cholesterol_Total <- read.csv('Cholesterol_Total.csv', header = F)$V1

# Complete_Data <- merge(count, result, by = 0)
# rownames(Complete_Data)=Complete_Data$Row.names
# Complete_Data <- subset(Complete_Data, select = -c(Row.names) )
# write.csv(Complete_Data, "Complete_Data.csv", row.names=T)
Complete_Data <- read.csv("Complete_Data.csv", row.names=1, header = T)
# all (Complete_Data$Fibroblast_1 %in% count$Fibroblast_1)
# all (Complete_Data$Fibroblast_2 %in% count$Fibroblast_2)
# all (Complete_Data$SHSY5Y_1 %in% count$SHSY5Y_1) 
# all (Complete_Data$SHSY5Y_2 %in% count$SHSY5Y_2)

# Up_Data <- merge(count, result_up, by = 0)
# rownames(Up_Data)=Up_Data$Row.names
# Up_Data <- subset(Up_Data, select = -c(Row.names) )
# write.csv(Up_Data, "Up_Data.csv", row.names=T)
Up_Data <- read.csv("Up_Data.csv", row.names=1, header = T)

# Down_Data <- merge(count, result_down, by = 0)
# rownames(Down_Data)=Down_Data$Row.names
# Down_Data <- subset(Down_Data, select = -c(Row.names) )
# write.csv(Down_Data, "Down_Data.csv", row.names=T)
Down_Data <- read.csv("Down_Data.csv", row.names=1, header = T)

# Nochange_Data <- merge(count, result_nochange, by = 0)
# rownames(Nochange_Data)=Nochange_Data$Row.names
# Nochange_Data <- subset(Nochange_Data, select = -c(Row.names) )
# write.csv(Nochange_Data, "Nochange_Data.csv", row.names=T)
Nochange_Data <- read.csv("Nochange_Data.csv", row.names=1, header = T)

# Endosome_Data <- Complete_Data[row.names(Complete_Data) %in% Endosome, ]
# write.csv(Endosome_Data, "Endosome_Data.csv", row.names = T)
Endosome_Data <- read.csv("Endosome_Data.csv", row.names=1, header = T)

# Cholesterol_Total_Data <- Complete_Data[row.names(Complete_Data) %in% Cholesterol_Total, ]
# write.csv(Cholesterol_Total_Data, "Cholesterol_Total_Data.csv", row.names = T)
Cholesterol_Total_Data <- read.csv("Cholesterol_Total_Data.csv", row.names=1, header = T)

# Endosome_Up_Data <- Up_Data[row.names(Up_Data) %in% Endosome, ]
# write.csv(Endosome_Up_Data, "Endosome_Up_Data.csv", row.names = T)
Endosome_Up_Data <- read.csv("Endosome_Up_Data.csv", row.names=1, header = T)

# Endosome_Down_Data <- Down_Data[row.names(Down_Data) %in% Endosome, ]
# write.csv(Endosome_Down_Data, "Endosome_Down_Data.csv", row.names = T)
Endosome_Down_Data <- read.csv("Endosome_Down_Data.csv", row.names=1, header = T)

# Endosome_Nochange_Data <- Nochange_Data[row.names(Nochange_Data) %in% Endosome, ]
# write.csv(Endosome_Nochange_Data, "Endosome_Nochange_Data.csv", row.names = T)
Endosome_Nochange_Data <- read.csv("Endosome_Nochange_Data.csv", row.names=1, header = T)

# Cholesterol_Total_Up_Data <- Up_Data[row.names(Up_Data) %in% Cholesterol_Total, ]
# write.csv(Cholesterol_Total_Up_Data, "Cholesterol_Total_Up_Data.csv", row.names = T)
Cholesterol_Total_Up_Data <- read.csv("Cholesterol_Total_Up_Data.csv", row.names=1, header = T)

# Cholesterol_Total_Down_Data <- Down_Data[row.names(Down_Data) %in% Cholesterol_Total, ]
# write.csv(Cholesterol_Total_Down_Data, "Cholesterol_Total_Down_Data.csv", row.names = T)
Cholesterol_Total_Down_Data <- read.csv("Cholesterol_Total_Down_Data.csv", row.names=1, header = T)

# Cholesterol_Total_Nochange_Data <- Nochange_Data[row.names(Nochange_Data) %in% Cholesterol_Total, ]
# write.csv(Cholesterol_Total_Nochange_Data, "Cholesterol_Total_Nochange_Data.csv", row.names = T)
Cholesterol_Total_Nochange_Data <- read.csv("Cholesterol_Total_Nochange_Data.csv", row.names=1, header = T)

















