setwd("/Users/madhavmakkena/Downloads/RNAseq/Subsets/Raw_Goterm_lists")

#creating a list of all file names
go_lists <- list.files(pattern="*.tsv")


library("dplyr")
library("stringr")
#creating a list of dataframes and all the data from the file names
comb_go_lists <- lapply(go_lists,function(x) {
  read.table(file = x, 
             sep = '\t', 
             header = TRUE)
})


#combining the dataframe lists into a dataframe
comb_go_lists <- bind_rows(comb_go_lists)

#Ordering based on the Go.Names and subsetting to have only go.names and symbols
comb_go_lists <- comb_go_lists[order(comb_go_lists$GO.NAME),]
comb_go_lists <- subset(comb_go_lists, select = c(SYMBOL,GO.NAME) )

# removing duplicates (only removed when common in both go.name and symbol)
comb_go_lists <- distinct(comb_go_lists)

# creating a list of all the unique gonames and getting the length of this list
go_name_list <- c(unique(comb_go_lists$GO.NAME))
go_name_list_len <- length(go_name_list)

#creating a new vector used in the for loop to temporarily store names
go_name <- vector(,1)

comb_go_lists_subset <- comb_go_lists

total_go_list <- distinct(subset(comb_go_lists, select = -c(GO.NAME)))
# write.csv(total_go_list, "total_go_list.csv", row.names = FALSE, Header = FALSE)

#creating a directory
setwd("/Users/madhavmakkena/Downloads/RNAseq/Subsets")
dir.create("all_subsets")

#subsetting and creating multiple csv files based on the go-names with the symbols as the data
for (i in 1:go_name_list_len) {
  go_name <- go_name_list[i]
  comb_go_lists_subset <- subset(comb_go_lists, GO.NAME==go_name)
  comb_go_lists_subset <- subset(comb_go_lists_subset, select = -c(GO.NAME) )
  write.csv(comb_go_lists_subset, str_interp("all_subsets/${go_name}.csv"), row.names = FALSE)
}

