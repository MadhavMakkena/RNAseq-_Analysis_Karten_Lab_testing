#setting the directory
setwd("/Users/madhavmakkena/Downloads/RNAseq")

#loading the count data (only when its already in a single file)
count_data <- read.csv("combined_single.csv", header = TRUE, sep = ",", row.names = "ENSGene")
#head(count_data, 10)

#loading the sample information
sample_data <- read.csv('sample_conditions.csv', header = TRUE, row.names = "sample")
#sample_data

#checking if the sample names in sample data match the count data
all (rownames(sample_data) %in% colnames(count_data))
all (rownames(sample_data) == colnames(count_data))
#both should be TRUE

#cleaning up count_data to remove genes with 0 counts in all samples
count_data<- count_data[rowSums(count_data) != 0, ]
#viewing the first 10 records to check
#head(count_data,10)

#DESeq analysis
#storing the input values from count_data_clean
library(DESeq2)
DESeq <- DESeqDataSetFromMatrix(countData = count_data, colData = sample_data, design = ~ condition)

#Excluding rows with less than 20 combined total reads or 10 mean total reads
keep <- rowSums(counts(DESeq)) >= 20
DESeq <- DESeq[keep,]
rm(keep)
keep <- rowMeans(counts(DESeq)) >= 10
DESeq <- DESeq[keep,]
rm(keep)

DESeq_hist <- sort(rowMeans(counts(estimateSizeFactors(DESeq), normalized = T)), decreasing = T)
summary(DESeq)
summary(DESeq_hist)

#plotting expression count histograms
hist(DESeq_hist, breaks = 1000000, col = "grey", xlim=c(0,100), ylim=c(0,100))

#performing DESeq2 on exp_100_sort
DESeq <-DESeq(DESeq)
plotDispEsts(DESeq)

#normalising the counts exp
count <- as.data.frame(counts(DESeq, normalized=TRUE))
head(count,10)

#extracting the results from DESeq
result <- as.data.frame(results(DESeq, contrast=c("condition", "Fibroblast", "SH-SY5Y")))
head(result,10)

#extracting the conditions from DESeq 
cond<-as.character(DESeq$condition)

#installing packages for BinfTools
# BiocManager::install("SAGx")
# BiocManager::install("GSVA")
# BiocManager::install("fgsea")
# BiocManager::install("gage")
# BiocManager::install("qusage")
# devtools::install_github("kevincjnixon/gpGeneSets")
# devtools::install_github("kevincjnixon/BinfTools")
lapply(c("SAGx", "GSVA", "fgsea", "gage", "qusage", "gpGeneSets", "BinfTools"), require, character.only = TRUE)

# dds is DESeq
# res is result
# counts is counts

# getting HGNC Names from ENSG names for result_exp100
result_HGCN <- result
result_HGCN <- getSym(object=result_HGCN,
               obType="res",
               species="hsapiens",
               target="HGNC",
               addCol=F)
head(result_HGCN)

# getting HGNC Names from ENSG names for counts_100
counts_HGCN <- counts
counts_HGCN <- getSym(object=counts_HGCN,
                        obType="counts",
                        species="hsapiens",
                        target="HGNC",
                        addCol=F)
head(counts_HGCN)








#specific gene sets onloading $V1 takes the first row and reads it as a vector
# Gene_set_chol_biosyn <- read.csv('Gene_set_cholesterol_biosynthesis.csv', header = FALSE)$V1
# Gene_set_chol_homeo <- read.csv('Gene_set_cholesterol_homeostasis.csv', header = FALSE)$V1
# Gene_set_KEGG_endocytosis <- read.csv('Gene_set_KEGG_endocytosis.csv', header = FALSE)$V1
#
# #gene sets exp
# Gene_set_chol_biosyn <- subset(result_HGCN, rownames(result_HGCN) %in% Gene_set_chol_biosyn)
# Gene_set_chol_homeo <- subset(result_HGCN, rownames(result_HGCN) %in% Gene_set_chol_homeo)
# Gene_set_KEGG_endocytosis <- subset(result_HGCN, rownames(result_HGCN) %in% Gene_set_KEGG_endocytosis)