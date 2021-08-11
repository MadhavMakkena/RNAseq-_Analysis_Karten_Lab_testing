#setting the directory
setwd("/Users/madhavmakkena/Downloads/RNAseq")

#loading the count data (already in a single file)
count_data <- read.csv("combined_single.csv", header = TRUE, sep = ",", row.names = "ENSGene")



#loading the sample information
#sample_files <- grep(".txt", list.files("/Users/madhavmakkena/Downloads/RNAseq/Counts"), value = TRUE)
sample_data <- read.csv('metadata.csv', header = TRUE, row.names = 1)



#checking if the sample names in sample data match the count data
#should be TRUE
all (rownames(sample_data) %in% colnames(count_data))
all (rownames(sample_data) == colnames(count_data))



#cleaning up count data to remove genes with 0 coutns in all conditions
count_data_clean <- count_data
count_data_clean <- count_data_clean[rowSums(count_data_clean) !=0, ]
#viewing the first 10 records to check
head(count_data,10)
head(count_data_clean,10)



#writing the clean count files
write.csv(count_data_clean,"count_data_clean.csv", row.names = TRUE)




#using the coutn matrix csv as the input
dds <- DESeqDataSetFromMatrix(countData = count_data_clean, colData = sample_data, design = ~ dex)




#removing the lowest 1/3rd expressing genes
dds_high <-estimateSizeFactors(dds)
dds_high_sort <- sort(rowMeans(counts(dds_high, normalized = T)), decreasing = T)
dds_high_sort <- dds_high_sort[1:(round(length(dds_high_sort)*.67))]
dds_high <- dds_high[rownames(counts(dds_high)) %in% names(dds_high_sort),]
summary(dds_high_sort)
summary(dds_high)



dds_all <-estimateSizeFactors(dds)
dds_all_sort <- sort(rowMeans(counts(dds_all, normalized = T)), decreasing = T)
dds_all_sort <- dds_all_sort[1:(round(length(dds_all_sort)*1))]
summary(dds_all_sort)
summary(dds_all)




#plotting dds histogram
detach(mtcars)
attach(mtcars)
par(mfrow=c(1,2))
hist(dds_all_sort, breaks = 1000000, col = "grey", xlim=c(0,100), ylim=c(0,100))
hist(dds_high_sort, breaks = 1000000, col = "grey", xlim=c(0,100), ylim=c(0,100))





#performing DESeq2 and plotting the dispersion
dds_run <-DESeq(dds)
detach(mtcars)
attach(mtcars)
plotDispEsts(dds_run)



dds_run_high <-DESeq(dds)
detach(mtcars)
attach(mtcars)
plotDispEsts(dds_run_high)






#normalising the counts
normal_counts <-counts(dds_run, normalized=TRUE)
result <- results(dds_run)
result
write.csv(as.data.frame(normal_counts), file = "normal_counts_all_dds_expressed.csv")



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



