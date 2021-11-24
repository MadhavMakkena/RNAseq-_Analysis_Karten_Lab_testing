# RNASeqAnalysis


# Only 2 conditions (e.g. fibroblast vs. SH-SY5Y) supported; number of replicates in each conditions have to be the same
# Only .csv files accepted
# Return the combined counts table to count_combined_ENSG if ENSG
# Return the combined counts table to count_combined_HNGC if HNGC
# Using single count files 
count_combined_ENSG <- single_files(folder_location="/Users/madhavmakkena/Desktop/RNAseq_Analysis/Counts/Seperate_files")
# Using already combined count files
# Include file extensions in the file name
# count_combined_ENSG <- combined_file(file_location="/Users/madhavmakkena/Desktop/RNAseq_Analysis/Counts", file_name= "count_combined.csv")


# Using above count_combined_ENSG and converting the ENSG names to gene symbols (HGNC)
# This uses kevin nixon's BinfTools package
# Return the combined counts with HGNC names table to count_combined_HGNC
count_combined_HGNC <- ENSG_to_HGNC(count_combined_ENSG)
# If already have a count_combined_HGNC file
# Set working directory using pathname
# setwd("/pathname")
# count_combined_HGNC <- read.csv("filename", header = TRUE, row.names = 1)


# Removing genes with a mean count of less than equal to (min_count_avg) in the count_combined_HGNC table
# i.e. only keeping genes with more than min_count_avg
# Return the cleaned combined counts table to count
count <- clean_up_min_counts(count_combined_HGNC, min_count_avg=0)
# THIS IS UNPROCESSED count FILE
# If already have a count file
# Set working directory using pathname
# setwd("/pathname")
# count <- read.csv("filename", header = TRUE, row.names = 1)


# setting up the conditions of the experiment using the column names from Count
# Return the conditions to Condition
Condition <- cond(Count)
# If error or already have a condition file
# Set working directory using pathname
# setwd("/pathname")
# condition <- read.csv("filename", header = TRUE, row.names = 1)


# Using the cleaned up UNPROCESSED count files and the extracted Condition file to run DESeq
# Make sure DESeq2 is installed
# Return the run_DESeq to DESeq_Result
# Both DESeq and the extracted PROCESSED Results are exported as lists
DESeq_Result <- run_DESeq(count,Condition)
# to seperate the lists run the following
# Return the DESeq to DESeq
DESeq <- (DESeq=DESeq_Result$DESeq)
# Return the result to Result
Result <-(Result=DESeq_Result$Result)
# Return the PROCESSED Count to Count
# NOTE: count is UNPROCESSED and Count is PROCESSED
Count <- (Count=DESeq_Result$Count)


# DESeq_hist <- sort(rowMeans(counts(estimateSizeFactors(DESeq), normalized = T)), decreasing = T)
# summary(DESeq)
# summary(DESeq_hist)
#plotting expression count histograms
# hist(DESeq_hist, breaks = 1000000, col = "grey", xlim=c(0,100), , ylim=c(0,100))
#performing DESeq2 on DESeq
# DESeq <-DESeq(DESeq)
# plotDispEsts(DESeq)


# plotting the histograms of PROCESSED Counts
# 2 rows x 3 columns
hist_plot_2x3(file_location="/Users/madhavmakkena/Desktop/RNAseq_Analysis", file_name="hist_plot_genes.csv")
# more than that the figures get too small to be useful


# Plotting the volcano plot
# result required Result
# Subset file and file location are mandatory for this function
# return the values to DEGs
DEGs <- volcano_plot_subset(result=Result, 
                    pval=0.05, 
                    log2FC=1.5, 
                    subset_file_location="/Users/madhavmakkena/Desktop/RNAseq_Analysis", 
                    subset_file_name="hist_plot_genes.csv")
# If you don't want any subsets use this instead
# DEGs <- volcano_plot(result=Result, pval=0.05, log2FC=1.5)
# to separate the DEG lists run the following
# Return the result_up to result_up
result_up <- (result_up=DEGs$result_up)
# Return the result_down to result_down
result_down <-(result_down=DEGs$result_down)
# Return the result_nochange to result_nochange
result_nochange <- (result_nochange=DEGs$result_nochange)
# Checking if the totals add up between the DEGs and Result
all (nrow(Result) == (nrow(result_up)+nrow(result_down)+nrow(result_nochange)))


# Combines the data from Result and Count
# Return DEG_data to DEG_Data
# Returns a list of Complete_Data, Up_Data, Down_Data, and Nochange_Data
DEG_Data <- DEG_data(result=Result, 
                     count=Count, 
                     result_up=result_up, 
                     result_down=result_down, 
                     result_nochange=result_nochange)
# To separate the DEG_Data lists run the following
# Return the Complete_Data to Complete_Data
Complete_Data <- (Complete_Data=DEG_Data$Complete_Data)
# Return the Up_Data to Up_Data
Up_Data <-(Up_Data=DEG_Data$Up_Data)
# Return the Down_Data to Down_Data
Down_Data <- (Down_Data=DEG_Data$Down_Data)
# Return the Nochange_Data to Nochange_Data
Nochange_Data <- (Nochange_Data=DEG_Data$Nochange_Data)


# Combines the data from Result and Count
# Return GO_geneset_DEG to GO_geneset_DEG
# Returns a list of Complete_Data, Up_Data, Down_Data, and Nochange_Data
GO_geneset_DEG <- GO_geneset_DEG(geneset_file_location="/Users/madhavmakkena/Desktop/RNAseq_Analysis/Amigo_lists/Raw_lists", 
                                 geneset_file_name="GO_0005739_mitochondrion.csv", 
                                 Complete_Data=Complete_Data, 
                                 Up_Data=Up_Data, 
                                 Down_Data=Down_Data, 
                                 Nochange_Data=Nochange_Data)
# To separate the DEG_Data lists run the following
# Return the subset_genes to subset_genes
# Replace the Name and you can save the data under the following pathname
setwd("pathname")
Name_genes <- (subset_genes=GO_geneset_DEG$subset_genes)
write.csv(Name_genes, "Name_genes.csv", row.names=TRUE)
# Return the GO_geneset_Data to GO_geneset_Data
Name_Data <- (GO_geneset_Data=GO_geneset_DEG$GO_geneset_Data)
write.csv(Name_Data, "Name_Data.csv", row.names=TRUE)
# Return the GO_geneset_Up_Data to GO_geneset_Up_Data
Name_Up_Data <-(GO_geneset_Up_Data=GO_geneset_DEG$GO_geneset_Up_Data)
write.csv(Name_Up_Data, "Name_Up_Data.csv", row.names=TRUE)
# Return the GO_geneset_Down_Data to GO_geneset_Down_Data
Name_Down_Data <- (GO_geneset_Down_Data=GO_geneset_DEG$GO_geneset_Down_Data)
write.csv(Name_Down_Data, "Name_Down_Data.csv", row.names=TRUE)
# Return the GO_geneset_Nochange_Data to GO_geneset_Nochange_Data
Name_Nochange_Data <- (GO_geneset_Nochange_Data=GO_geneset_DEG$GO_geneset_Nochange_Data)
write.csv(Name_Nochange_Data, "Name_Nochange_Data.csv", row.names=TRUE)
# Checking if the totals add up
all (nrow(GO_geneset_DEG$GO_geneset_Data) == (nrow(GO_geneset_DEG$GO_geneset_Up_Data)+nrow(GO_geneset_DEG$GO_geneset_Down_Data)+nrow(GO_geneset_DEG$GO_geneset_Nochange_Data)))



