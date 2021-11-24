setwd("/Users/madhavmakkena/Desktop/R/New")

# Expression <- function(result){
#   # result_subset<-subset(result[c("ANXA6", "NPC1", "ALK", "COL1A1"),])
#   # result_subset<-subset(result[c(row.names(Cholesterol_Total_Data)),])
#   My_Up<-subset(result, pvalue<=0.05 & log2FoldChange >= 1.5)
#   My_Down<-subset(result, pvalue<=0.05 & log2FoldChange <= -1.5)
#   My_NoChange<-subset(result, (pvalue<=0.05 & (log2FoldChange >-1.5 &  log2FoldChange < 1.5)) | (pvalue>0.05))
#   plot(abs(result$stat)~result$log2FoldChange, pch=16, cex=0.6, col="black", ylim = c(-10,(ceiling(max(abs(result$stat))))), xlim = c(-(ceiling(max(abs(result$log2FoldChange)))),(ceiling(max(abs(result$log2FoldChange))))), xlab="Log2FoldChange", ylab="Absolute Wald Statistic", main="Differentially Expressed Genes [Fibroblast vs. SH_SY5Y]", cex.main=.85)
#   text(-22,130, "Down", col = "Blue", pos = 4, cex = .85)
#   text(-22,120, "regulated", col = "Blue", pos = 4, cex = .85)
#   text(-22,140, (count(result_up)), col = "Blue", pos = 4, cex = .85)
#   text(22,130, "Up", col = "Red", pos = 2, cex = .85)
#   text(22,120, "regulated", col = "Red", pos = 2, cex = .85)
#   text(22,140, (count(result_down)), col = "Red", pos = 2, cex = .85)
#   text(-22,-8, (count(result_nochange)), col = "Black", pos = 4, cex = .85)
#   points(abs(result_up$stat)~result_up$log2FoldChange, pch=16, cex=0.6, col="pink")
#   points(abs(result_down$stat)~result_down$log2FoldChange, pch=16, cex=0.6, col="lightblue")
#   # points(abs(result_subset$stat)~result_subset$log2FoldChange, pch=16, cex=0.6, col="darkgreen")
#   abline(h=min(abs(subset(result, pvalue<0.05)$stat)), v=c(1.5,-1.5), lty=2)
# }
# 
# Expression(result)
devtools::update_packages("BinfTools")


test <- volcanoPlota(res=result, #Results object
                         title="l",
                         pval=0.05, #adjusted p-value threshold for DEGs
                         FC=log2(1.5), #log2FoldChange threshold for DEGs (can be 0)
                         lab=NULL, #list of genes to label (NULL to not label any)
                         col=NULL, # in red list of genes to colour (NULL to not colour any)
                         fclim=NULL, #x-axis (log2FoldChange) limits, genes passing this limit will be represented as triangles on the edge of the plot - good if you have some extreme outliers
                         showNum=T, #Show the numbers of genes on the plot?
                         returnDEG=T, #Return list of DEGs (Down, Up) - this is good for running GO later on
                         expScale=F, #Scale point size to mean expression?
                         upcol="lightblue3", #Colour value for upregulated genes, NULL will be red
                         dncol="rosybrown3") #Colour value for downregulated genes, NULL will be blue)

Old_up_l
Old_down_l
Old_nochange_l


Old_Up <- result[row.names(result) %in% Old_up_l,]
Old_Down <- result[row.names(result) %in% Old_down_l,]
Old_NoChange <- result[row.names(result) %in% Old_nochange_l,]

setwd("/Users/madhavmakkena/Desktop/R/Compare_method")


write.csv(My_Down, "My_Down.csv", row.names = T)
write.csv(My_NoChange, "My_NoChange.csv", row.names = T)
write.csv(My_Up, "My_Up.csv", row.names = T)
write.csv(Old_Down, "Old_Down.csv", row.names = T)
write.csv(Old_NoChange, "Old_NoChange.csv", row.names = T)
write.csv(Old_Up, "Old_Up.csv", row.names = T)
#




































Old_up_l <- c((testlist_up=test$Up))
Old_down_l <- c((testlist_down=test$Down))

summary(Old_up_l)
summary(Old_down_l)

result_l<-c((row.names(result)))
summary(result_l)

Old_nochange_l <- setdiff(result_l, Old_up_l)
Old_nochange_l <- setdiff(Old_nochange_l, Old_down_l)

summary(Old_nochange_l)

#
rm(temp)





















rm(count)




setwd("/Users/madhavmakkena/Desktop/R/Sets")
W_NC <- read.csv("const_data_1.5FC.csv", header = T, row.names = 1)
W_D <- read.csv("down_data_1.5FC.csv", header = T, row.names = 1)
W_U <- read.csv("up_data_1.5FC.csv", header = T, row.names = 1)



