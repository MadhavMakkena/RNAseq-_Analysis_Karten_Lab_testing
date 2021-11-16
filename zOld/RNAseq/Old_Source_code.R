#setting the directory
setwd("/Users/madhavmakkena/Downloads/RNAseq")


#loading the count data (only when its already in a single file)
count_data_whole <- read.csv("combined_single.csv", header = TRUE, sep = ",", row.names = "ENSGene")


#loading the sample information
sample_data <- read.csv('sample_conditions.csv', header = TRUE, row.names = "sample")


#checking if the sample names in sample data match the count data
all (rownames(sample_data) %in% colnames(count_data_whole))
all (rownames(sample_data) == colnames(count_data_whole))
#both should be TRUE



#cleaning up count_data to remove genes with 0 counts in all samples
count_data<- count_data_whole[rowSums(count_data_whole) !=0, ]
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
# exp_60 <-estimateSizeFactors(DESeq_input)
# exp_60_sort <- sort(rowMeans(counts(exp_60, normalized = T)), decreasing = T)
# exp_60_sort <- exp_60_sort[1:(round(length(exp_60_sort)*.60))]
# exp_60 <- exp_60[rownames(counts(exp_60)) %in% names(exp_60_sort),]
# summary(exp_60_sort)
# summary(exp_60)



#plotting expression count histograms
hist(exp_100_sort, breaks = 1000000, col = "grey", xlim=c(0,100), ylim=c(0,100))
# hist(exp_60_sort, breaks = 1000000, col = "grey", xlim=c(0,100), ylim=c(0,100))


#performing DESeq2 on exp_100_sort
DESeq_100 <-DESeq(exp_100)
plotDispEsts(DESeq_100)
#
#performing DESeq2 on exp_60_sort
# DESeq_60 <-DESeq(exp_60)
# plotDispEsts(DESeq_60)

#getting the conditions of DESeq
cond<-as.character(DESeq_100$condition)


#normalising the counts exp_100
normal_count_exp_100 <-counts(DESeq_100, normalized=TRUE)
counts_100 <- as.data.frame(normal_count_exp_100)
normal_count_exp_100
counts_100
#
#normalising the counts exp_60
# normal_count_exp_60 <-counts(DESeq_60, normalized=TRUE)
# counts_60 <- as.data.frame(normal_count_exp_60)
# normal_count_exp_60
# counts_60


#extracting the results from DESeq_100
result_100 <- results(DESeq_100)
result_100
result_100<-as.data.frame(result_100)
result_100
#
#extracting the results from DESeq_60
# result_60 <- results(DESeq_60)
# result_60
# result_60<-as.data.frame(result_60)
# result_60


#installing packages for BinfTools
BiocManager::install("SAGx")
BiocManager::install("GSVA")
BiocManager::install("fgsea")
BiocManager::install("gage")
BiocManager::install("qusage")
devtools::install_github("kevincjnixon/gpGeneSets")
devtools::install_github("kevincjnixon/BinfTools")
library(SAGx)
library(GSVA)
library(fgsea)
library(gage)
library(qusage)
library(gpGeneSets)
library(BinfTools)

library('biomaRt')
G <- result_100
mart = useDataset("hsapiens_gene_ensembl", useEnsembl(biomart="ensembl"))
genes <- rownames(G)
G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=genes,mart= mart)

duplicated(G_list)
G_list[duplicated(G_list)]














summary(G_list)
# dds is DESeq_100 or DESeq_60
# res is result_100 or result_60
# counts is counts_100 or counts_60
#
# getting HGNC Names from ENSG names for result_exp100
sym_result_100 <- result_100
sym_result_100 <- getSym(object=sym_result_100, obType="res", species="hsapiens", target="HGNC", addCol=F)
head(sym_result_100)
head(result_100)


sym_result_100 <- gConvert(result_100, in.pre = NULL, in.var = NULL, out.type = c("inc", "rad"),
         out.pre = out.type)

gconvert(query = c("REAC:R-HSA-3928664", "rs17396340", "NLRP1"), organism = "hsapiens",
         target="ENSG", mthreshold = Inf, filter_na = TRUE)




#
# getting HGNC Names from ENSG names for result_exp60
sym_result_60 <- result_60
sym_result_60 <- getSym(object=sym_result_60,
                         obType="res",
                         species="hsapiens",
                         target="HGNC",
                         addCol=F)
head(sym_result_60)
head(result_60)
#
# # getting HGNC Names from ENSG names for counts_100
sym_counts_100 <- counts_100
sym_counts_100 <- getSym(object=sym_counts_100,
                        obType="counts",
                        species="hsapiens",
                        target="HGNC",
                        addCol=F)
head(sym_counts_100)
# 
# #
# # getting HGNC Names from ENSG names for counts_60
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


dir.create("GO")

DEG<-volcanoPlot(res=sym_result_100, #Results object
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
               bg=rownames(sym_result_100), #A character vector of genes indicating the background for the GO analysis. Leave NULL to use all genes (if you don't have one)
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
                bg=rownames(sym_result_100),
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

rnk<-GenerateGSEA(res=sym_result_100, #Results object
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

sym_counts_100_1 <- data.matrix(sym_counts_100)

gsva_plot(counts=as.matrix(sym_counts_100_1), #counts object (as matrix), make sure rownames are the same nomenclature as the gene symbols in geneset
          geneset=forVenn,
          method="gsva", #Method for gsva plot - see documentation for options
          condition=cond,
          con="Fibroblast", #Indicate the control condition
          title="cholesterol ssGSEA", 
          compare=NULL, #for pairwise t-tests, leave NULL to do all possible comparisons, or provide a list of vectors, length 2 indicating the conditions to compare
          col="Dark2", #Colour scheme, can be RColourBrewer palette name, or vector of rgb(), hexadecimal, or colour names
          style="violin")

count_plot(counts=sym_counts_100_1,
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
             counts=sym_counts_100_1,
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
      counts=sym_counts_100_1,
      conditions=cond,
      con="Fibroblast",
      title="cholesterol - cut tree clusters",
      labgenes=NULL,
      zscore=T,
      rclus=annodf)

app<-exploreData(res=sym_result_100, #Results object, or named list of results objects
                 counts=sym_counts_100, #Normalized counts or name list of normalized counts - must be same names as res
                 cond=cond)

app











barGene(genes=Gene_set_chol_biosyn,
        counts=log10(1+sym_counts_100),
        conditions=cond,
        title="Cholesterol Biosynthesis",
        norm="Fibroblast",
        eb="se",
        returnDat=T,
        ord=c("Fibroblast","SH-SY5Y"),
        col="Dark2")

barGene(genes=Gene_set_chol_homeo,
        counts=log10(1+sym_counts_100),
        conditions=cond,
        title="Cholesterol Homeostasis",
        norm="Fibroblast",
        eb="se",
        returnDat=T,
        ord=c("Fibroblast","SH-SY5Y"),
        col=c("blue","yellow"))

DEG<-volcanoPlot(res=sym_result_100, #Results object
                 title="Fibroblast vs SH-SY5",
                 p=0.05, #adjusted p-value threshold for DEGs
                 pval=NULL, #unadjusted p-value threshold for DEGs (in case you don't want to use adjusted)
                 FC=log2(1.5), #log2FoldChange threshold for DEGs (can be 0)
                 lab=NULL, #list of genes to label (NULL to not label any)
                 col=Gene_set_chol_biosyn, #list of genes to colour (NULL to not colour any)
                 fclim=NULL, #x-axis (log2FoldChange) limits, genes passing this limit will be represented as triangles on the edge of the plot - good if you have some extreme outliers
                 showNum=T, #Show the numbers of genes on the plot?
                 returnDEG=T, #Return list of DEGs (Down, Up) - this is good for running GO later on
                 expScale=F, #Scale point size to mean expression?
                 upcol=NULL, #Colour value for upregulated genes, NULL will be red
                 dncol=NULL) #Colour value for downregulated genes, NULL will be blue)

MA_Plot(res=sym_result_100,
        title="Fibroblast vs SH-SY5",
        p=0.05,
        pval=NULL,
        FC=log2(1.5),
        lab=Gene_set_chol_biosyn,
        col=Gene_set_chol_biosyn,
        fclim=NULL, #Same as volcano plot, but will act on y-axis, not x
        showNum=T,
        returnDEG=F,
        sigScale=F, #Scale point size to significance?
        upcol=NULL,
        dncol=NULL)
