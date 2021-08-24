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
count_data<- count_data[rowSums(count_data) !=0, ]
#viewing the first 10 records to check
head(count_data,10)



#DESeq analysis
library(DESeq2)
#storing the input values from count_data_clean
DESeq_input <- DESeqDataSetFromMatrix(countData = count_data, colData = sample_data, design = ~ condition)




#Keeping all the genes
exp_100 <-estimateSizeFactors(DESeq_input)
exp_100_sort <- sort(rowMeans(counts(exp_100, normalized = T)), decreasing = T)
exp_100_sort <- exp_100_sort[1:(round(length(exp_100_sort)*1))]
exp_100 <- exp_100[rownames(counts(exp_100)) %in% names(exp_100_sort),]
summary(exp_100_sort)
summary(exp_100)
#
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
DESeq_100 <-DESeq(exp_100)
plotDispEsts(DESeq_100)
#
#performing DESeq2 on exp_60_sort
DESeq_60 <-DESeq(exp_60)
plotDispEsts(DESeq_60)
#getting the conditions of DESeq
cond<-as.character(DESeq_100$condition)


#normalising the counts exp_100
normal_count_exp_100 <-counts(DESeq_100, normalized=TRUE)
counts_100 <- as.data.frame(normal_count_exp_100)
normal_count_exp_100
counts_100
#
#normalising the counts exp_60
normal_count_exp_60 <-counts(DESeq_60, normalized=TRUE)
counts_60 <- as.data.frame(normal_count_exp_60)
normal_count_exp_60
counts_60



#extracting the results from DESeq_100
result_100 <- as.data.frame(results(DESeq_100))
result_100
#
#extracting the results from DESeq_60
result_60 <- results(DESeq_60)
result_60



#installing packages for BinfTools
BiocManager::install("SAGx")
BiocManager::install("GSVA")
BiocManager::install("fgsea")
BiocManager::install("gage")
BiocManager::install("qusage")
devtools::install_github("kevincjnixon/gpGeneSets")
devtools::install_github("kevincjnixon/BinfTools")
library(BinfTools)



# dds is DESeq_100 or DESeq_60
# res is result_100 or result_60
# counts is counts_100 or counts_60
#
# getting HGNC Names from ENSG names for result_exp100
sym_result_100 <- result_100
sym_result_100 <- getSym(object=sym_result_100,
               obType="res",
               species="hsapiens",
               target="HGNC",
               addCol=F)
head(sym_result_100)
#
# getting HGNC Names from ENSG names for result_exp60
sym_result_60 <- result_60
sym_result_60 <- getSym(object=sym_result_60,
                         obType="res",
                         species="hsapiens",
                         target="HGNC",
                         addCol=F)
head(sym_result_60)
#
# getting HGNC Names from ENSG names for counts_100
sym_counts_100 <- counts_100
sym_counts_100 <- getSym(object=sym_counts_100,
                        obType="counts",
                        species="hsapiens",
                        target="HGNC",
                        addCol=F)
head(sym_counts_100)
#
# getting HGNC Names from ENSG names for counts_60
sym_counts_60 <- counts_60
sym_counts_60 <- getSym(object=sym_counts_60,
                        obType="counts",
                        species="hsapiens",
                        target="HGNC",
                        addCol=F)
head(sym_counts_60)



#specific gene sets onloading $V1 takes the first row and reads it as a vector
Gene_set_chol_biosyn <- read.csv('Gene_set_cholesterol_biosynthesis.csv', header = FALSE)$V1
Gene_set_chol_homeo <- read.csv('Gene_set_cholesterol_homeostasis.csv', header = FALSE)$V1
Gene_set_KEGG_endocytosis <- read.csv('Gene_set_KEGG_endocytosis.csv', header = FALSE)$V1
#
# #gene sets exp100
# exp100_Gene_set_chol_biosyn <- subset(sym_result_exp100, rownames(sym_result_exp100) %in% Gene_set_chol_biosyn)
# exp100_Gene_set_chol_homeo <- subset(sym_result_exp100, rownames(sym_result_exp100) %in% Gene_set_chol_homeo)
# exp100_Gene_set_KEGG_endocytosis <- subset(sym_result_exp100, rownames(sym_result_exp100) %in% Gene_set_KEGG_endocytosis)
# #
# #gene sets exp100
# exp60_Gene_set_chol_biosyn <- subset(sym_result_exp60, rownames(sym_result_exp60) %in% Gene_set_chol_biosyn)
# exp60_Gene_set_chol_homeo <- subset(sym_result_exp60, rownames(sym_result_exp60) %in% Gene_set_chol_homeo)
# exp60_Gene_set_KEGG_endocytosis <- subset(sym_result_exp60, rownames(sym_result_exp60) %in% Gene_set_KEGG_endocytosis)



# volcano plots and MA plots use either sym_result_100 or sym_result_60
# barGene(), count_plot(), and gsva_plot() use either sym_counts_100 or sym_counts_100








dim(sym_result_exp100)
dim(cond)
view(cond)
cond<-as.character(DESeq_run_exp_100$condition)
barGene(genes=Gene_set_chol_biosyn,
        counts=sym_result_exp100,
        conditions=cond,
        title="sym_result_exp100",
        norm="WT",
        eb="se",
        col=c("blue","yellow"),
        ord=NULL)

barGene(genes=Gene_set_chol_biosyn,
        counts=log10(1+counts_exp100),
        conditions=cond,
        title="sym_result_exp100",
        norm="WT",
        eb="se",
        col=c("blue","yellow"),
        ord=NULL)
