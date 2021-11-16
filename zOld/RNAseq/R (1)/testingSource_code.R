# setting the directory
setwd("/Users/madhavmakkena/Desktop/R")

# loading the count data (only when its already in a single file)
# raw_counts <- read.csv("combined_single.csv", header = TRUE, sep = ",", row.names = "ENSGene")
# barplot(colSums(raw_counts), col= c("C_1" = "orange","C_2" = "orange","S_1" = "lightblue", "S_2" = "lightblue"))
# head(raw_counts, 10)
# raw_counts <- getSym(object=raw_counts,
#                         obType="counts",
#                         species="hsapiens",
#                         target="HGNC",
#                         addCol=F)
# head(raw_counts)
# write.csv(raw_counts,"raw_counts_HGNC.csv", row.names = TRUE)
# if already performed ^ then
raw_counts <- read.csv("raw_counts_HGNC.csv", header = TRUE, row.names = 1)
head(raw_counts)
# loading the sample information
sample_cond <- read.csv('sample_conditions.csv', header = TRUE, row.names = "sample")
sample_cond

# checking if the sample names in sample data match the count data
all (rownames(sample_cond) %in% colnames(raw_counts))
all (rownames(sample_cond) == colnames(raw_counts))
#both should be TRUE


# cleaning up raw_counts to remove genes with 0 counts in ALL samples
raw_counts<- raw_counts[rowSums(raw_counts) != 0, ]
# cleaning up raw_counts to remove genes with 0 counts in ANY samples
#raw_counts <- raw_counts[apply(raw_counts!=0, 1, all),]
# viewing the first 10 records to check
# head(raw_counts,10)

#DESeq analysis
#storing the input values from raw_counts
library(DESeq2)
DESeq <- DESeqDataSetFromMatrix(countData = raw_counts, colData = sample_cond, design = ~ condition)

#Excluding rows with less than 20 combined total reads or 10 mean total reads
# keep <- rowSums(counts(DESeq)) >= 20
# DESeq <- DESeq[keep,]
# rm(keep)
# keep <- rowMeans(counts(DESeq)) >= 10
# DESeq <- DESeq[keep,]
# rm(keep)

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
# write.csv(count, "count.csv", row.names = TRUE)



#extracting the results from DESeq
# result <- as.data.frame(results(DESeq, contrast=c("condition", "Fibroblast", "SH-SY5Y")))
result <- as.data.frame(results(DESeq))
head(result,10)
# write.csv(result, "result.csv", row.names = TRUE)


#extracting the conditions from DESeq 
cond<-as.character(DESeq$condition)

#exploreData(res=result, counts=count, cond=cond)

#installing packages for BinfTools
# BiocManager::install("SAGx")
# BiocManager::install("GSVA")
# BiocManager::install("fgsea")
# BiocManager::install("gage")
# BiocManager::install("qusage")
# install.packages("devtools")
# devtools::install_github("kevincjnixon/gpGeneSets")
# devtools::install_github("kevincjnixon/BinfTools", force = TRUE)
# devtools::update_packages("BinfTools")
# library("BinfTools")
lapply(c("SAGx", "GSVA", "fgsea", "gage", "qusage", "gpGeneSets", "BinfTools"), require, character.only = TRUE)

# dds is DESeq
# res is result
# counts is counts

# getting HGNC Names from ENSG names for result
# result_HGNC <- result
# result_HGNC <- getSym(object=result_HGNC,
#                obType="res",
#                species="hsapiens",
#                target="HGNC",
#                addCol=F)
# head(result_HGNC)
# write.csv(result_HGNC,"result_HGNC.csv", row.names = TRUE)
#if already performed
# result_HGNC <- read.csv("result_HGNC.csv", row.names = 1, header = TRUE)

# getting HGNC Names from ENSG names for count
# count_HGNC <- count
# count_HGNC <- getSym(object=count_HGNC,
#                         obType="counts",
#                         species="hsapiens",
#                         target="HGNC",
#                         addCol=F)
# head(count_HGNC)
# write.csv(count_HGNC,"count_HGNC.csv", row.names = TRUE)
#if already performed
# count_HGNC <- read.csv("count_HGNC.csv", row.names = 1, header = TRUE)


#cholesterol subset
# gene_subset <- read.csv("cholesterol_subsets/1_cholesterol_gene_sets_combined.csv", header = TRUE, row.names = 1)
# 
# count_HGNC_cholesterol <-count_HGNC
# result_HGNC_cholesterol <-result_HGNC
# 
# library(dplyr)
# count_HGNC_cholesterol <- subset(count_HGNC_cholesterol, row.names(count_HGNC_cholesterol) %in% gene_subset$SYMBOL)
# result_HGNC_cholesterol <- subset(result_HGNC_cholesterol, row.names(result_HGNC_cholesterol) %in% gene_subset$SYMBOL)
#
# write.csv(count_HGNC_cholesterol, "count_HGNC_cholesterol.csv")
# write.csv(result_HGNC_cholesterol, "result_HGNC_cholesterol.csv")

# count_HGNC_cholesterol <- read.csv("count_HGNC_cholesterol.csv", row.names = 1, header = TRUE)
# result_HGNC_cholesterol <- read.csv("result_HGNC_cholesterol.csv", row.names = 1, header = TRUE)


#        , ylim=c(0, 7000)
#
par(mfrow=c(2,2))
plotCounts(DESeq, gene="PSMB2", intgroup="condition", transform=F, ylim=c(0, 16000))
plotCounts(DESeq, gene="NPC1", intgroup="condition", transform=F, ylim=c(0, 16000))
plotCounts(DESeq, gene="NPC2", intgroup="condition", transform=F, ylim=c(0, 16000))
plotCounts(DESeq, gene="C7", intgroup="condition", transform=F, ylim=c(0, 16000))
plotCounts(DESeq, gene="DBH", intgroup="condition", transform=F, ylim=c(0, 16000))
plotCounts(DESeq, gene="ALK", intgroup="condition", transform=F)
plotCounts(DESeq, gene="PHOX2B", intgroup="condition", transform=F)
plotCounts(DESeq, gene="PRRX1", intgroup="condition", transform=F)
plotCounts(DESeq, gene="COL5A1", intgroup="condition", transform=F)



# plotCounts(DESeq, gene="STARD3NL", intgroup="condition", transform=F)
# plotCounts(DESeq, gene="STAR", intgroup="condition", transform=F)

chol_list <- read.csv("1_cholesterol_gene_sets_combined.csv", header = F)$V2
chol_list_gmt <- list(GeneSet1=c(chol_list))

chol_list_c <- read.csv("1_cholesterol_gene_sets_combined copy.csv", header = F)$V2
chol_list_gmt_c <- list(GeneSet1=c(chol_list))

volcano_gene_list <- c("NPC2", "NPC1", "ABCA1", "STARD3", "STARD4", "HMGCS1", "INHBA","SMO","SCARB2","SULT2B1","SOAT2")
volcano_gene_lista <- c("NPC2", "NPC1", "ABCA1", "STARD3", "STARD4")
volcano_gene_listb <- c("HMGCS1", "INHBA","SMO","SCARB2","SULT2B1","SOAT2")



summary(volcano_gene_list)
geneSetz<-list(GeneSet1=c("NPC2","NPC1","ABCA1"))

geneSet<-list(GeneSet1=c("NPC2","NPC1","ABCA1"), GeneSet2=c("STARD3","STARD4","AACS","HMGCS1"), GeneSet3=c("INHBA","SMO","SCARB2","SULT2B1","SOAT2"))
geneSets<- c("NPC2","NPC1","ABCA1", "STARD3","STARD4","AACS","HMGCS1", "INHBA","SMO","SCARB2","SULT2B1","SOAT2")
par(mfrow=c(1,1))
testlist <- volcanoPlota(res=result, #Results object
                        title="F",
                        p=0.005, #adjusted p-value threshold for DEGs
                        FC=log2(1.5), #log2FoldChange threshold for DEGs (can be 0)
                        lab=NULL, #list of genes to label (NULL to not label any)
                        col=chol_list, # in red list of genes to colour (NULL to not colour any)
                        cola=brain_list, # in green list of genes to colour (NULL to not colour any)
                        fclim=NULL, #x-axis (log2FoldChange) limits, genes passing this limit will be represented as triangles on the edge of the plot - good if you have some extreme outliers
                        showNum=T, #Show the numbers of genes on the plot?
                        returnDEG=T, #Return list of DEGs (Down, Up) - this is good for running GO later on
                        expScale=F, #Scale point size to mean expression?
                        upcol="lightblue3", #Colour value for upregulated genes, NULL will be red
                        dncol="rosybrown3") #Colour value for downregulated genes, NULL will be blue)
#figure out exactly what the axis mean


testlist <- volcanoPlot(res=result, #Results object
                        title="Fibroblast vs SH_SY5Y",
                        p=0.05, #adjusted p-value threshold for DEGs
                        FC=log2(1.5), #log2FoldChange threshold for DEGs (can be 0)
                        lab=volcano_gene_list, #list of genes to label (NULL to not label any)
                        col=volcano_gene_list, #list of genes to colour (NULL to not colour any)
                        fclim=NULL, #x-axis (log2FoldChange) limits, genes passing this limit will be represented as triangles on the edge of the plot - good if you have some extreme outliers
                        showNum=T, #Show the numbers of genes on the plot?
                        returnDEG=T, #Return list of DEGs (Down, Up) - this is good for running GO later on
                        expScale=F, #Scale point size to mean expression?
                        upcol=NULL, #Colour value for upregulated genes, NULL will be red
                        dncol=NULL) #Colour value for downregulated genes, NULL will be blue)
#figure out exactly what the axis mean




testlist_up <- c(list(testlist_up=testlist$Up))
testlist_down <- c(list(testlist_down=testlist$Down))
# write.csv(testlist_up, "testlist_up.csv", row.names = TRUE)
# write.csv(testlist_down, "testlist_down.csv", row.names = FALSE)

library("BinfTools")
gsva_plot(counts=as.matrix(count), #counts object (as matrix), make sure rownames are the same nomenclature as the gene symbols in geneset
          geneset=tot_list_gmt,
          method="gsva", #Method for gsva plot - see documentation for options
          condition=cond,
          con="Fibroblast", #Indicate the control condition
          title="Rhodospin ssGSEA", 
          compare=NULL, #for pairwise t-tests, leave NULL to do all possible comparisons, or provide a list of vectors, length 2 indicating the conditions to compare
          col="Dark2", #Colour scheme, can be RColourBrewer palette name, or vector of rgb(), hexadecimal, or colour names
          style="violin") #If not 'violin' it will be a box plot


?gsva_plot()

count_plot(
  counts=count,
  scaling = "zscore",
  genes=chol_list,
  condition=cond,
  con ="Fibroblast",
  title = "Volcano Gene Set",
  compare = NULL,
  col = "Dark2",
  method = "ind",
  pair = F,
  pc = 1,
  yax = NULL,
  showStat = T,
  style = "violin"
)

zheat(
  genes = chol_list,
  counts=count,
  conditions=cond,
  con = "Fibroblast",
  title = "Volcano Gene Set",
  labgenes = NULL, #leave null to label all genes, set to “” to label no genes, or type a character vector of select genes to label.
  zscore = T, #Plot z-score normalized expression or set to FALSE to plot raw expression, you may want to use counts=log10(1+count) if you do this though
  rclus = F, #Do you want to cluster the rows using hierarchical clustering?
  hmcol = NULL,
  retClus = F
)



