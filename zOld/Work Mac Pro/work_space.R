setwd("/Users/madhavmakkena/Desktop/R/Sets")

Complete_Data <- read.csv("Complete_Data.csv", header = T, row.names = 1)

UP <- read.csv("up_1.5FC.csv", header = F)$V1
up_data <- read.csv("up_data_1.5FC.csv", header = T, row.names = 1)

DOWN <- read.csv("down_1.5FC.csv", header = F)$V1
down_data <- read.csv("down_data_1.5FC.csv", header = T, row.names = 1)

CONST <- read.csv("const_1.5FC.csv", header = F)$V1
const_data <- read.csv("const_data_1.5FC.csv", header = T, row.names = 1)


Endosome_list <- read.csv("Endo_list.csv", header = F)$V1
# Endosome_list_data <- Complete_Data[row.names(Complete_Data) %in% Endosome_list, ]
# write.csv(Endosome_list_data, "Endosome_list_data.csv", row.names = T)
Endosome_list_data <- read.csv("Endosome_list_data.csv", header = T, row.names = 1)

chol_list <- read.csv("chol_list.csv", header = F)$V1
# chol_list_data <- Complete_Data[row.names(Complete_Data) %in% chol_list, ]
# write.csv(chol_list_data, "chol_list_data.csv", row.names = T)
chol_list_data <- read.csv("chol_list_data.csv", header = T, row.names = 1)

# Chol_up_data <- chol_list_data[row.names(chol_list_data) %in% row.names(up_data),]
# Chol_down_data <- chol_list_data[row.names(chol_list_data) %in% row.names(down_data),]
# Chol_const_data <- chol_list_data[row.names(chol_list_data) %in% row.names(const_data),]
# write.csv(Chol_up_data, "Chol_up_data.csv", row.names = T)
# write.csv(Chol_down_data, "Chol_down_data.csv", row.names = T)
# write.csv(Chol_const_data, "Chol_const_data.csv", row.names = T)
Chol_up_data <- read.csv("Chol_up_data.csv", header = T, row.names = 1)
Chol_down_data <- read.csv("Chol_down_data.csv", header = T, row.names = 1)
Chol_const_data <- read.csv("Chol_const_data.csv", header = T, row.names = 1)

# Endo_up_data <- Endosome_list_data[row.names(Endosome_list_data) %in% row.names(up_data),]
# Endo_down_data <- Endosome_list_data[row.names(Endosome_list_data) %in% row.names(down_data),]
# Endo_const_data <- Endosome_list_data[row.names(Endosome_list_data) %in% row.names(const_data),]
# write.csv(Endo_up_data, "Endo_up_data.csv", row.names = T)
# write.csv(Endo_down_data, "Endo_down_data.csv", row.names = T)
# write.csv(Endo_const_data, "Endo_const_data.csv", row.names = T)
Endo_up_data <- read.csv("Endo_up_data.csv", header = T, row.names = 1)
Endo_down_data <- read.csv("Endo_down_data.csv", header = T, row.names = 1)
Endo_const_data <- read.csv("Endo_const_data.csv", header = T, row.names = 1)


test1 <- c(row.names(Endo_up_data))
test2 <- c(row.names(Endo_down_data))
test3 <- c(row.names(Endo_const_data))




any(row.names(Chol_CONST_data)=="APOA5")
test




statVolcano<-function(res){
  up<-subset(res, padj<0.005 & log2FoldChange > 1)
  down<-subset(res, padj<0.005 & log2FoldChange < -1)
  plot(abs(res$stat)~res$log2FoldChange, pch=16, cex=0.6, col="black", xlab="Log2FoldChange", ylab="Absolute Wald Statistic", main="DEGs")
  points(abs(up$stat)~up$log2FoldChange, pch=16, cex=0.6, col="red")
  points(abs(down$stat)~down$log2FoldChange, pch=16, cex=0.6, col="blue")
  # test <- c("APOA5")
  test <- subset(res, row.names = "APOA5")
  points(abs(test$stat)~test$log2FoldChange, pch=16, cex=0.6, col="green")
  abline(h=min(abs(subset(res, padj<0.05)$stat)), v=c(1,-1), lty=2)
}

statVolcano(result)

