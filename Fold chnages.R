setwd("/Users/madhavmakkena/Desktop/R/New")
list.files()
result <- read.csv("result.csv", header = T, row.names = 1)
count <- read.csv("count.csv", header = T, row.names = 1)

Expression <- function(result){
  result_subset<-subset(result[c("HMGCR","NPC2"),])
  # result_subset<-subset(result[c(row.names(Cholesterol_Total_Data)),])
  result_up<-subset(result, pvalue<=0.05 & log2FoldChange >= 1.5)
  result_down<-subset(result, pvalue<=0.05 & log2FoldChange <= -1.5)
  result_nochange<-subset(result, (pvalue<=0.05 & (log2FoldChange >-1.5 &  log2FoldChange < 1.5)) | (pvalue>0.05))
  plot(abs(result$stat)~result$log2FoldChange, pch=16, cex=0.6, col="black", ylim = c(-10,(ceiling(max(abs(result$stat))))), xlim = c(-(ceiling(max(abs(result$log2FoldChange)))),(ceiling(max(abs(result$log2FoldChange))))), xlab="Log2FoldChange", ylab="Absolute Wald Statistic", main="Differentially Expressed Genes [Fibroblast vs. SH_SY5Y]", cex.main=.85)
  text(-22,130, "Down", col = "Blue", pos = 4, cex = .85)
  text(-22,120, "regulated", col = "Blue", pos = 4, cex = .85)
  text(-22,140, (nrow(result_down)), col = "Blue", pos = 4, cex = .85)
  text(22,130, "Up", col = "Red", pos = 2, cex = .85)
  text(22,120, "regulated", col = "Red", pos = 2, cex = .85)
  # text(22,140, (nrow(result_up)), col = "Red", pos = 2, cex = .85)
  # text(-22,-8, (nrow(result_nochange)), col = "Black", pos = 4, cex = .85)
  # print((nrow(result_up)))
  points(abs(result_up$stat)~result_up$log2FoldChange, pch=16, cex=0.6, col="grey")
  points(abs(result_down$stat)~result_down$log2FoldChange, pch=16, cex=0.6, col="grey")
  points(abs(result_subset$stat)~result_subset$log2FoldChange, pch=16, cex=0.6, col="darkgreen")
  abline(h=min(abs(subset(result, pvalue<0.05)$stat)), v=c(1.5,-1.5), lty=2)
}

Expression(result)


result_up


install.packages("plyr")

test5<-c("a")

# write.csv(result_up, "result_up.csv", row.names=TRUE)
# write.csv(result_down, "result_down.csv", row.names=TRUE)
# write.csv(result_nochange, "result_nochange.csv", row.names=TRUE)

Expression(result)

subset(count[c("ANXA6", "NPC1", "ALK"),])

result[c("C8A"),]
count[c("C8A"),]

statVolcano<-function(res){
  up<-subset(res, padj<0.005 & log2FoldChange > 1)
  down<-subset(res, padj<0.005 & log2FoldChange < -1)
  plot(abs(res$stat)~res$log2FoldChange, pch=16, cex=0.6, col="black", xlab="Log2FoldChange", ylab="Absolute Wald Statistic", main="DEGs")
  points(abs(up$stat)~up$log2FoldChange, pch=16, cex=0.6, col="red")
  points(abs(down$stat)~down$log2FoldChange, pch=16, cex=0.6, col="blue")
  abline(h=min(abs(subset(res, padj<0.05)$stat)), v=c(1,-1), lty=2)
}
