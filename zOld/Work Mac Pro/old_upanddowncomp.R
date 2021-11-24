setwd("/Users/madhavmakkena/Desktop/R")

UP <- read.csv("up_1.5FC.csv", header = F)$V1
up_list <- Complete_Data[row.names(Complete_Data) %in% UP, ]
write.csv(up_list, "up_data_1.5FC.csv", row.names = T)

DOWN <- read.csv("down_1.5FC.csv", header = F)$V1
down_list <- Complete_Data[row.names(Complete_Data) %in% DOWN, ]
write.csv(down_list, "down_data_1.5FC.csv", row.names = T)

const_list <- Complete_Data[row.names(Complete_Data) %in% CONST, ]
write.csv(const_list, "const_data_1.5FC.csv", row.names = T)



#
#
#
#
#
#
#
# Complete_Data <- merge(count, result, by = 0)
# rownames(Complete_Data)=Complete_Data$Row.names
# Complete_Data <- subset(Complete_Data, select = -c(Row.names) )
# colnames(Complete_Data)[1] <- "Fibroblast_1"
# colnames(Complete_Data)[2] <- "Fibroblast_2"
# colnames(Complete_Data)[3] <- "SH-SY5Y_1"
# colnames(Complete_Data)[4] <- "SH-SY5Y_2"
# 
# write.csv(Complete_Data, "Complete_Data.csv", row.names = T)


#
#
#
Endosome_list_data <- Complete_Data[row.names(Complete_Data) %in% Endosome_list, ]
chol_list_data <- Complete_Data[row.names(Complete_Data) %in% chol_list, ]

UP_data <- Complete_Data[row.names(Complete_Data) %in% UP, ]
DOWN_data <- Complete_Data[row.names(Complete_Data) %in% DOWN, ]


Endo_Up_List <- intersect(Endosome_list, UP)
write.csv(Endo_Up_List, "Endo_Up_List.csv", row.names = F)
Endo_Up_List_Data <- intersect(Endosome_list_data, UP_data)
write.csv(Endo_Up_List_Data, "Endo_Up_List_Data.csv", row.names = T)

Endo_Down_List <- intersect(Endosome_list, DOWN)
write.csv(Endo_Down_List, "Endo_Down_List.csv", row.names = F)
Endo_Down_List_Data <- intersect(Endosome_list_data, DOWN_data)
write.csv(Endo_Down_List_Data, "Endo_Down_List_Data.csv", row.names = T)

Chol_Up_List <- intersect(chol_list, UP)
write.csv(Chol_Up_List, "Chol_Up_List.csv", row.names = F)
Chol_Up_List_Data <- intersect(chol_list_data, UP_data)
write.csv(Chol_Up_List_Data, "Chol_Up_List_Data.csv", row.names = T)

Chol_Down_List <- intersect(chol_list, DOWN)
write.csv(Chol_Down_List, "Chol_Down_List.csv", row.names = F)
Chol_Down_List_Data <- intersect(chol_list_data, DOWN_data)
write.csv(Chol_Down_List_Data, "Chol_Down_List_Data.csv", row.names = T)

install.packages('prob')

Endo_const_List_Data <- Complete_Data[row.names(Complete_Data) %in% Endo_const_list, ]
write.csv(Endo_const_List_Data, "Endo_const_List_Data.csv", row.names = T)
Endo_const_list <- c(row.names(Endo_const_List_Data))
write.csv(Endo_const_list, "Endo_const_list.csv", row.names = F)

Chol_const_List_Data <- Complete_Data[row.names(Complete_Data) %in% Chol_const_list, ]
write.csv(Chol_const_List_Data, "Chol_const_List_Data.csv", row.names = T)
Chol_const_list <- c(row.names(Chol_const_List_Data))
write.csv(Chol_const_list, "Chol_const_list.csv", row.names = F)


#
#
#