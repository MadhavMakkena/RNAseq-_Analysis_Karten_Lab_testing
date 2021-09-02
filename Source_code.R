# setting the directory
setwd("/Users/madhavmakkena/Downloads/RNAseq")

# loading the count data (already in a single file and names replaced with HGNC)
raw_counts <- read.csv("raw_counts_HGNC.csv", header = TRUE, row.names = 1)

# loading the sample information
sample_cond <- read.csv('sample_conditions.csv', header = TRUE, row.names = "sample")


# checking if the sample names in sample data match the count data
all (rownames(sample_cond) %in% colnames(raw_counts))
all (rownames(sample_cond) == colnames(raw_counts))
#both should be TRUE


# cleaning up raw_counts to remove genes with 0 counts in ALL samples
raw_counts<- raw_counts[rowSums(raw_counts) != 0, ]


#DESeq prep
library(DESeq2)
DESeq <- DESeqDataSetFromMatrix(countData = raw_counts, colData = sample_cond, design = ~ condition)
DESeq_hist <- sort(rowMeans(counts(estimateSizeFactors(DESeq), normalized = T)), decreasing = T)
summary(DESeq)
summary(DESeq_hist)

#plotting expression count histograms
# hist(DESeq_hist, breaks = 1000000, col = "grey", xlim=c(0,100), , ylim=c(0,100))

#performing DESeq2 on DESeq
DESeq <-DESeq(DESeq)
# plotDispEsts(DESeq)

#normalising the counts exp
count <- as.data.frame(counts(DESeq, normalized=TRUE))
head(count,10)

#extracting the results from DESeq
result <- as.data.frame(results(DESeq))
head(result,10)

# extracting the conditions from DESeq 
# cond<-as.character(DESeq$condition)

#installing packages for BinfTools
# BiocManager::install("SAGx")
# BiocManager::install("GSVA")
# BiocManager::install("fgsea")
# BiocManager::install("gage")
# BiocManager::install("qusage")
# devtools::install_github("kevincjnixon/gpGeneSets")
# devtools::install_github("kevincjnixon/BinfTools", force = TRUE)
# devtools::update_packages("BinfTools")
# library("BinfTools")
lapply(c("SAGx", "GSVA", "fgsea", "gage", "qusage", "gpGeneSets", "BinfTools"), require, character.only = TRUE)

#
count_HGNC_cholesterol <- read.csv("count_HGNC_cholesterol.csv", row.names = 1, header = TRUE)
result_HGNC_cholesterol <- read.csv("result_HGNC_cholesterol.csv", row.names = 1, header = TRUE)

#
par(mfrow=c(2,3))
plotCounts(DESeq, gene="PSMB2", intgroup="condition", transform=F, ylim=c(50, 4000))
plotCounts(DESeq, gene="NPC1", intgroup="condition", transform=F, ylim=c(50, 4000))
plotCounts(DESeq, gene="NPC2", intgroup="condition", transform=F, ylim=c(50, 4000))
plotCounts(DESeq, gene="ABCA1", intgroup="condition", transform=F, ylim=c(50, 4000))
plotCounts(DESeq, gene="STARD3", intgroup="condition", transform=F, ylim=c(50, 4000))
plotCounts(DESeq, gene="STARD4", intgroup="condition", transform=F, ylim=c(50, 4000))
plotCounts(DESeq, gene="HMGCR", intgroup="condition", transform=F, ylim=c(50, 4000))





#log works test <- log(count)



# PSMB2 

plotCounts(DESeq, gene="PSMB2", intgroup="condition", transform=F, ylim=c(100, 4000))

install.packages("remotes")
remotes::install_github("acidgenomics/acidplots")
install.packages("acidgenomics/acidplots")

dds <- DESeq
res <- result
res <- res[order(res$padj),]
head(res)


