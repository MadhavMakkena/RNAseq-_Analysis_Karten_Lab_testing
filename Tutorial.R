# Only 2 conditions supported; number of replicates in each conditions have to be the same
# Only .csv files accepted
# Return the combined counts table to count_combined_ENSG if ENSG
# Return the combined counts table to count_combined_HNGC if HNGC
# Using single count files 
count_combined_ENSG <- single_files(folder_location="/Users/madhavmakkena/Desktop/RNAseq_Analysis/Counts/Seperate_files")
# Using already combined count files
# Include file extensions in the file name
count_combined_ENSG <- combined_file(file_location="/Users/madhavmakkena/Desktop/RNAseq_Analysis/Counts",
                                     file_name= "count_combined.csv")


# Using above count_combined_ENSG and converting the ENSG names to gene symbols (HGNC)
# This uses kevin nixon's BinfTools package
# Return the combined counts with HGNC names table to count_combined_HGNC
count_combined_HGNC <- ENSG_to_HGNC(count_combined_ENSG)
# If already have a count_combined_HGNC file
# Set working directory using pathname
# setwd("/pathname")
# count_combined_HGNC <- read.csv("filename", header = TRUE, row.names = 1)


# Removing genes with a mean count of 0 in the count_combined_HGNC table
# Return the cleaned combined counts table to count
Count <- Clean_up_min_counts(count_combined_HGNC)
# If already have a Count file
# Set working directory using pathname
# setwd("/pathname")
# Count <- read.csv("filename", header = TRUE, row.names = 1)




