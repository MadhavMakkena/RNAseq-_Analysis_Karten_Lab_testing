#setting the directory
setwd("/Users/madhavmakkena/Downloads/RNAseq")


#loading the data
GO_annotations_cholesterol_combined <- read.csv("GO_annotations_cholesterol_combined.csv", sep = ",")

library(dplyr)
order <- GO_annotations_cholesterol_combined[order(GO_annotations_cholesterol_combined$GO.NAME),]
order <- subset(order, select = c(SYMBOL,GO.NAME) )
order <- distinct(order)


combined <- subset(order, select = -c(GO.NAME) )
combined <- distinct(combined)
write.csv(combined, "cholesterol_subsets/1_cholesterol_gene_sets_combined.csv")
list <- c(unique(order$GO.NAME))
len <- length(list)

goname <- vector(,1)

subset <-order

library(stringr)

dir.create("cholesterol_subsets")

for (i in 1:len) {
  goname <- list[i]
  subset <- subset(order, GO.NAME==goname)
  subset <- subset(subset, select = -c(GO.NAME) )
  # write.csv(subset, str_interp("cholesterol_subsets/${goname}.csv"), row.names = FALSE)
}
