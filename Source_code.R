#setting the directory
setwd("/Users/madhavmakkena/Downloads/RNAseq")

#loading the count data (only when its already in a single file)
count_data <- read.csv("combined_single.csv", header = TRUE, sep = ",", row.names = "ENSGene")
head(count_data, 10)

#loading the sample information
sample_data <- read.csv('sample_conditions.csv', header = TRUE, row.names = "sample")
sample_data

#checking if the sample names in sample data match the count data
all (rownames(sample_data) %in% colnames(count_data))
all (rownames(sample_data) == colnames(count_data))
#both should be TRUE

#cleaning up count_data to remove genes with 0 counts in all samples
count_data_clean <- count_data
count_data_clean <- count_data_clean[rowSums(count_data_clean) !=0, ]
#viewing the first 10 records to check
head(count_data,10)
head(count_data_clean,10)

#DESeq analysis
library(DESeq2)
#storing the input values from count_data_clean
DESeq_input <- DESeqDataSetFromMatrix(countData = count_data_clean, colData = sample_data, design = ~ condition)

#Keeping all the genes
exp_100 <-estimateSizeFactors(DESeq_input)
exp_100_sort <- sort(rowMeans(counts(exp_100, normalized = T)), decreasing = T)
exp_100_sort <- exp_100_sort[1:(round(length(exp_100_sort)*1))]
exp_100 <- exp_100[rownames(counts(exp_100)) %in% names(exp_100_sort),]
summary(exp_100_sort)
summary(exp_100)

#Keeping the top 60% expressing the genes
exp_60 <-estimateSizeFactors(DESeq_input)
exp_60_sort <- sort(rowMeans(counts(exp_60, normalized = T)), decreasing = T)
exp_60_sort <- exp_60_sort[1:(round(length(exp_60_sort)*.60))]
exp_60 <- exp_60[rownames(counts(exp_60)) %in% names(exp_60_sort),]
summary(exp_60_sort)
summary(exp_60)

#plotting expression count histograms
hist(exp_100_sort, breaks = 1000000, col = "grey", xlim=c(0,100), ylim=c(0,100))
hist(exp_60_sort, breaks = 1000000, col = "grey", xlim=c(0,100), ylim=c(0,100))

#performing DESeq2 on exp_100_sort
DESeq_run_exp_100 <-DESeq(exp_100)
plotDispEsts(DESeq_run_exp_100)

#performing DESeq2 on exp_60_sort
DESeq_run_exp_60 <-DESeq(exp_60)
plotDispEsts(DESeq_run_exp_60)

#normalising the counts exp_100
# normal_count_exp_100 <-counts(DESeq_run_exp_100, normalized=TRUE)
# normal_count_exp_100

#normalising the counts exp_60
# normal_count_exp_60 <-counts(DESeq_run_exp_60, normalized=TRUE)
# normal_count_exp_60

#extracting the results from DESeq_run_exp_100
result_exp100 <- results(DESeq_run_exp_100)
result_exp100

#extracting the results from DESeq_run_exp_60
result_exp60 <- results(DESeq_run_exp_60)
result_exp60















#adding gene names
result$ensembl <- sapply( strsplit( rownames(result), split="\\+" ), "[", 1 )
library( "biomaRt" )
ensembl = useMart( "ensembl", dataset = "hsapiens_gene_ensembl" )
genemap <- getBM( attributes = c("ensembl_gene_id", "hgnc_symbol"),
                  filters = "ensembl_gene_id",
                  values = result$ensembl,
                  mart = ensembl )
idx <- match( result$ensembl, genemap$ensembl_gene_id )
result$entrez <- genemap$entrezgene[ idx ]
result$hgnc_symbol <- genemap$hgnc_symbol[ idx ]
head(result,4)
write.csv( as.data.frame(result), file="results_all_dds_expressed_GeneName.csv")


