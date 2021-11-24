setwd("/Users/madhavmakkena/Desktop/RNAseq_Analysis")
# DESeq Analysis
library(DESeq2)
# DESeq <- DESeqDataSetFromMatrix(countData = count_raw, colData = condition, design = ~ sample_type)

# DESeq histogram
# DESeq_hist <- sort(rowMeans(counts(estimateSizeFactors(DESeq), normalized = T)), decreasing = T)
# summary(DESeq)
# summary(DESeq_hist)
# DESeq_hist_plot <- hist(DESeq_hist, breaks = 1000000, col = "grey", xlim=c(0,100), , ylim=c(0,100))

# performing DESeq2 on DESeq
# DESeq <-DESeq(DESeq)
# plotDispEsts(DESeq)

# count <- as.data.frame(counts(DESeq, normalized=TRUE))
# write.csv(count, "count.csv", row.names = T)
count <- read.csv("count.csv", header = TRUE, row.names = 1)

# extracting the results from DESeq
# result <- as.data.frame(results(DESeq))
# write.csv(result, "result.csv", row.names = T)
result <- read.csv("result.csv", header = TRUE, row.names = 1)

cond<-as.character(DESeq$sample_type)


#removing non necessary data
rm(count_raw)
rm(test)
rm(DESeq_hist)
rm(DESeq_hist_plot)
