setwd("/Users/madhavmakkena/Desktop/R/Sets/atlas_GO_lists")
# MultiVesicular <- read.table(file = 'multivesicular_list.csv', sep = '\t', header = TRUE)
# MultiVesicular <- MultiVesicular["SYMBOL"]
# MultiVesicular <- distinct(MultiVesicular)
# write.csv(MultiVesicular, "MultiVesicular.csv", row.names = F)
MultiVesicular <- read.csv('multivesicular_list.csv', header = F)$V1

setwd("/Users/madhavmakkena/Desktop/R/New/Gene_Sets")
MultiVesicular_Data <- Complete_Data[row.names(Complete_Data) %in% MultiVesicular, ]
write.csv(MultiVesicular_Data, "MultiVesicular_Data.csv", row.names = T)
MultiVesicular_Data <- read.csv("MultiVesicular_Data.csv", row.names=1, header = T)

MultiVesicular_Up_Data <- Up_Data[row.names(Up_Data) %in% MultiVesicular, ]
write.csv(MultiVesicular_Up_Data, "MultiVesicular_Up_Data.csv", row.names = T)
MultiVesicular_Up_Data <- read.csv("MultiVesicular_Up_Data.csv", row.names=1, header = T)

MultiVesicular_Down_Data <- Down_Data[row.names(Down_Data) %in% MultiVesicular, ]
write.csv(MultiVesicular_Down_Data, "MultiVesicular_Down_Data.csv", row.names = T)
MultiVesicular_Down_Data <- read.csv("MultiVesicular_Down_Data.csv", row.names=1, header = T)

MultiVesicular_Nochange_Data <- Nochange_Data[row.names(Nochange_Data) %in% MultiVesicular, ]
write.csv(MultiVesicular_Nochange_Data, "MultiVesicular_Nochange_Data.csv", row.names = T)
MultiVesicular_Nochange_Data <- read.csv("MultiVesicular_Nochange_Data.csv", row.names=1, header = T)






