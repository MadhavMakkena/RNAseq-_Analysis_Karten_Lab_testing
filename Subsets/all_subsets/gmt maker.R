setwd("/Users/madhavmakkena/Downloads/RNAseq/Aarons_subsets")
subset_lists <- list.files(pattern="*.csv")
subset_lists <- lapply(subset_lists,function(x) {
  read.table(file = x, 
             sep = '\t', 
             header = TRUE)
})

subset_lists <- lapply(subset_list,function(x) {
  read.csv(file =x, header = FALSE)$V1
  chol_gene_list <- chol_gene_list[-1];
  
})


