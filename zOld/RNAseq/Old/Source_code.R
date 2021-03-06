# setting the directory
setwd("/Users/madhavmakkena/Downloads/RNAseq")

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

# loading the sample information
sample_cond <- read.csv('sample_conditions.csv', header = TRUE, row.names = "sample")
# sample_cond

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

count_HGNC_cholesterol <- read.csv("count_HGNC_cholesterol.csv", row.names = 1, header = TRUE)
result_HGNC_cholesterol <- read.csv("result_HGNC_cholesterol.csv", row.names = 1, header = TRUE)



dds <- DESeq


vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)
head(assay(vsd), 3)
ntd <- normTransform(dds)
library("vsn")
meanSdPlot(assay(ntd))
meanSdPlot(assay(vsd))
meanSdPlot(assay(rld))

library("pheatmap")
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("condition")])
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)
sampleDists <- dist(t(assay(vsd)))

library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$type, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

plotPCA(vsd, intgroup=c("condition")) #when you have types of reads which doesnt apply to us


dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)
dds <- nbinomWaldTest(dds)


results(dds, contrast=c("condition", "C1_Fibroblast", "C2_Fibroblast"))

summary(res)
par(mar=c(8,5,2,2))
boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2)


























#creating a GO directory
dir.create("GO")


# volcano plots and MA plots use result_HGNC
# barGene(), count_plot(), and gsva_plot() use count_HGNC
#creating a volcano plot
DEG<-volcanoPlot(res=result_HGNC, #Results object
                 title="Fibroblast vs SH-SY5",
                 p=0.05, #adjusted p-value threshold for DEGs
                 pval=NULL, #unadjusted p-value threshold for DEGs (in case you don't want to use adjusted)
                 FC=log2(1.5), #log2FoldChange threshold for DEGs (can be 0)
                 lab=NULL, #list of genes to label (NULL to not label any)
                 col=NULL, #list of genes to colour (NULL to not colour any)
                 fclim=NULL, #x-axis (log2FoldChange) limits, genes passing this limit will be represented as triangles on the edge of the plot - good if you have some extreme outliers
                 showNum=T, #Show the numbers of genes on the plot?
                 returnDEG=T, #Return list of DEGs (Down, Up) - this is good for running GO later on
                 expScale=F, #Scale point size to mean expression?
                 upcol=NULL, #Colour value for upregulated genes, NULL will be red
                 dncol=NULL) #Colour value for downregulated genes, NULL will be blue)

GO_res<-GO_GEM(geneList=DEG,
               species="hsapiens",
               bg=rownames(result_HGNC), #A character vector of genes indicating the background for the GO analysis. Leave NULL to use all genes (if you don't have one)
               source="GO:BP", #A character indicating the source - see documentation for all of them
               corr="fdr", #How to correct the p-values
               iea=F, #Remove genes in terms 'inferred by electronic analysis' ?
               prefix="GO/", #Character for output prefix. If named list is provided as geneList, names of the list will be added to the prefix
               ts=c(10,500), #numeric of length 2 indicating minimum/maximum term size for results plots
               pdf=T, #print figures to pdf?
               fig=T, #Show figures in plots in RStudio?
               figCols=c("blue","orange"), #colours for enrichment/significance in plots
               returnGost=T, #Return gprofiler2 gost object
               writeRes=T, #Write results to ".GO.txt" file
               writeGem=T, #Write gem file?
               writeGene=F, #Write genes in query to file?
               returnRes=T) #Return the results table (only one of returnRes or returnGost can be T, not both)

GO_gost<-GO_GEM(geneList=DEG,
                species="hsapiens",
                bg=rownames(result_HGNC),
                source="GO:BP",
                prefix="GO/",
                pdf=F,
                fig=F,
                returnGost=T,
                writeRes=F)

BinfTools:::combGO_plot(GOresList=GO_res,
                        title="Biological process - significant",
                        ts=c(10,500),
                        sig=T,
                        numTerm=10,
                        upcols=c("lightpink","red"),
                        downcols=c("lightblue","blue"))

BinfTools:::combGO_plot(GOresList=GO_res,
                        title="Biological process - enriched",
                        ts=c(10,500),
                        sig=F,
                        numTerm=10,
                        upcols=c("lightpink","red"),
                        downcols=c("lightblue","blue"))

rnk<-GenerateGSEA(res=result_HGNC, #Results object
                  filename="GSEA.rnk", #Output rnk file name for GSEA preranked outside of R
                  bystat=T, #Rank genes by stats? will use Wald statistic or if not nere, -log10(p-value) with the direction from the log2FoldChange
                  byFC=F, #Rank genes by log2FoldChange? I like to use this with the shrunken log fold change from DESEq2
                  retRNK=T) #Return the RNK object? yes, to run GSEA in R


library(gpGeneSets)
gsea_res<-GSEA(rnk=rnk, #Rnk object
               gmt=gp_dm, #either .gmt filename or a loaded gene set
               pval=1, #adjusted p-value threshold for terms to return, set to 1 to return all terms and filter later
               ts=c(10,600), #min/max term sizes to filter terms BEFORE analysis
               nperm=10000, #number of permuations for p-value generation
               parseBader=F) #Set to TRUE if using Gary's genesets - it will parse the term names so it looks neater. I'm not using gary's genesets here, so we will set to FALSE

rows<-c(1, nrow(gsea_res))
enPlot(gseaRes=gsea_res[rows,], #GSEA results table subset into the rows of interest to make a plot
       rnk=rnk, #Original rnk object used to make gsea_res
       gmt=gp_dm, #Original gmt object/filename used to make gsea_res
       title=NULL) #character vector of length (nrow(gseaRes)) for custom titles, or leave NULL for automatic titles - works better for Gary's genesets

cholesterol<-c(customGMT(gost=GO_gost$Down, #gost object from GO_GEM and returnGost=T
                         key="cholesterol", #keyword to pull gene sets - this is a grep, so anything with this key will be pulled - the resulting geneset may require some manual curation so check the names
                         gmt=gp_dm), #The gpGeneSets object containing the complete gene sets 'gp_hs' for human, 'gp_mm' for mouse and 'gp_dm' for drosophila
               customGMT(gost=GO_gost$Up,
                         key="cholesterol",
                         gmt=gp_dm))


write.gmt(geneSet=cholesterol,
          filename="cholesterol.gmt")


forVenn<-list(DE_Up=DEG$Up,
              DE_Down=DEG$Down,
              cholesterol=unique(unlist(cholesterol)))

plotVenn(x=forVenn, #name list for plotting Venn diagram
         title="cholesterol related genes",
         cols="Dark2", #Colour scheme for plot
         lty="blank", #outlines for circles
         scale=F, #Scale to list sizes?
         retVals=F) #Return list of values in overlaps?

count_HGNC_matrix <- data.matrix(count_HGNC)

gsva_plot(counts=as.matrix(count_HGNC_matrix), #counts object (as matrix), make sure rownames are the same nomenclature as the gene symbols in geneset
          geneset=forVenn,
          method="gsva", #Method for gsva plot - see documentation for options
          condition=cond,
          con="Fibroblast", #Indicate the control condition
          title="cholesterol ssGSEA", 
          compare=NULL, #for pairwise t-tests, leave NULL to do all possible comparisons, or provide a list of vectors, length 2 indicating the conditions to compare
          col="Dark2", #Colour scheme, can be RColourBrewer palette name, or vector of rgb(), hexadecimal, or colour names
          style="violin")

count_plot(counts=count_HGNC_matrix,
           scaling="none", #Can be "zscore" to emphasize differences, or 'log10', or "none"
           genes=unique(unlist(forVenn)), #Character vector of gene names - need to unlist the geneset for this
           condition=cond,
           con="Fibroblast",
           title="Cholesterol Genes Expression",
           compare=NULL,
           col="Dark2",
           method="perMean", #What method to plot? "mean", "median", "perMean", "ind", "geoMean"
           pair=F, #Paired t-tests?
           pc=1, #pseudocount if scaling="log10"
           yax="Percent Mean Expression", #y-axis label if default isn't descriptive enough
           showStat=T, #Show statistics on plot?
           style="box") #Default is violin


htree<-zheat(genes=unique(unlist(forVenn)), #Character vector of genes to plot in heatmap, NULL will plot all genes
             counts=count_HGNC_matrix,
             conditions=cond,
             con="Fibroblast",
             title="cholesterol genes",
             labgenes="",#Character vecotr of gene names to label in heatmap, NULL labels all, "" will label none
             zscore=T, #Plot row-zscore? if FALSE, probably want to log transform counts
             rclus=T, #TRUE=hierarchical clustering, FALSE=order in decreasing expression in control condition, can also give it a dataframe with rownames=gene names and the first column with an identifier to cluster genes
             hmcol=NULL, #colorRampPalette of length 100 for custom heatmap colours (NULL=default colours)
             retClus=T) #return clustered objects if rclus=T - will be used to pull clustered genes later

annodf<-BinfTools:::heatClus(out=htree,
                             level=3)
head(annodf)

zheat(genes=unique(unlist(forVenn)),
      counts=count_HGNC_matrix,
      conditions=cond,
      con="Fibroblast",
      title="cholesterol - cut tree clusters",
      labgenes=NULL,
      zscore=T,
      rclus=annodf)

app<-exploreData(res=result_HGNC, #Results object, or named list of results objects
                 counts=count_HGNC_matrix, #Normalized counts or name list of normalized counts - must be same names as res
                 cond=cond)

app



#specific gene sets onloading $V1 takes the first row and reads it as a vector
# Gene_set_chol_biosyn <- read.csv('Gene_set_cholesterol_biosynthesis.csv', header = FALSE)$V1
# Gene_set_chol_homeo <- read.csv('Gene_set_cholesterol_homeostasis.csv', header = FALSE)$V1
# Gene_set_KEGG_endocytosis <- read.csv('Gene_set_KEGG_endocytosis.csv', header = FALSE)$V1
#
# #gene sets exp
# Gene_set_chol_biosyn <- subset(result_HGNC, rownames(result_HGNC) %in% Gene_set_chol_biosyn)
# Gene_set_chol_homeo <- subset(result_HGNC, rownames(result_HGNC) %in% Gene_set_chol_homeo)
# Gene_set_KEGG_endocytosis <- subset(result_HGNC, rownames(result_HGNC) %in% Gene_set_KEGG_endocytosis)