#histogram of counts 3x2 plot
par(mfrow=c(2,3))
plotCounts(DESeq, gene="PSMB2", intgroup="sample_type", transform=F, ylim=c(0, 4000))
plotCounts(DESeq, gene="NPC1", intgroup="sample_type", transform=F, ylim=c(0, 4000))
plotCounts(DESeq, gene="NPC2", intgroup="sample_type", transform=F, ylim=c(0, 4000))
plotCounts(DESeq, gene="STARD3", intgroup="sample_type", transform=F, ylim=c(0, 4000))
plotCounts(DESeq, gene="STARD4", intgroup="sample_type", transform=F, ylim=c(0, 4000))
plotCounts(DESeq, gene="ANXA6", intgroup="sample_type", transform=F)

#
library("BinfTools")

test <- c("NPC1", "ANXA6")
par(mfrow=c(1,1))
testlist <- volcanoPlot(res=result, #Results object
                         title="Fibroblast vs. SH-SY5Y [Endo_const; p=0.005; FC=1.5]",
                         pval=0.005, #adjusted p-value threshold for DEGs
                         FC=log2(1.5), #log2FoldChange threshold for DEGs (can be 0)
                         lab=test, #list of genes to label (NULL to not label any)
                         col=test, # in red list of genes to colour (NULL to not colour any)
                         fclim=NULL, #x-axis (log2FoldChange) limits, genes passing this limit will be represented as triangles on the edge of the plot - good if you have some extreme outliers
                         showNum=T, #Show the numbers of genes on the plot?
                         returnDEG=T, #Return list of DEGs (Down, Up) - this is good for running GO later on
                         expScale=F, #Scale point size to mean expression?
                         upcol="lightblue3", #Colour value for upregulated genes, NULL will be red
                         dncol="rosybrown3") #Colour value for downregulated genes, NULL will be blue)


testlist_up <- c(list(testlist_up=testlist$Up))
testlist_down <- c(list(testlist_down=testlist$Down))
write.csv(testlist_up, "up.csv", row.names = FALSE)
write.csv(testlist_down, "down.csv", row.names = FALSE)

statVolcano<-function(result){
  up<-subset(res, pvalue<0.005 & log2FoldChange > 1.5)
  down<-subset(res, pvalue<0.005 & log2FoldChange < -1.5)
  plot(abs(res$stat)~res$log2FoldChange, pch=16, cex=0.6, col="black", xlab="Log2FoldChange", ylab="Absolute Wald Statistic", main="DEGs")
  points(abs(up$stat)~up$log2FoldChange, pch=16, cex=0.6, col="red")
  points(abs(down$stat)~down$log2FoldChange, pch=16, cex=0.6, col="blue")
  abline(h=min(abs(subset(res, pvalue<0.05)$stat)), v=c(1.5,-1.5), lty=2)
}

statVolcano(result)


MA_Plot(res=result,
        title="KO vs WT",
        p=0.05,
        pval=NULL,
        FC=log2(1.5),
        lab=test,
        col=test,
        fclim=NULL, #Same as volcano plot, but will act on y-axis, not x
        showNum=T,
        returnDEG=F,
        sigScale=F, #Scale point size to significance?
        upcol=NULL,
        dncol=NULL)




