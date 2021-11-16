# combining the counts into one single file
setwd("/Users/madhavmakkena/Desktop/RNAseq_Analysis/Counts")
# extract all the individual count files (any number) and combine them into a single data frame
# library("plyr")
# count_files <- list.files()
# count_lists <- lapply(count_files, read.csv, sep=",", header=TRUE, row.names="ENSG")
# library("dplyr")
# library("purrr")
# count_combined <- map(count_lists, ~ .x %>%
#                         rownames_to_column("ENSG")) %>%
#   reduce(full_join, by = "ENSG") %>%
#   mutate(across(everything(), replace_na, 0))

# setting row names from the ENSG column
# library("tibble")
# count_combined <- column_to_rownames(count_combined, var = "ENSG")

# removing all the non ENSG rows
# count_raw_dirty <- count_combined %>% filter(str_detect(row.names(count_combined) , "ENSG*"))
# test <- subset(count_combined, row.names(count_combined) == "ENSG*")

# extracting all non ENSG rows into non_ENSG_count
# count_rownames_ENSG <- row.names(count_raw_dirty)
# non_ENSG_count <- subset(count_combined, rownames(count_combined) != count_rownames_ENSG)

# changing the ENSG names to Gene Symbol
# library("BinfTools")
# count_raw_dirty <- getSym(object=count_raw_dirty,
#                         obType="counts",
#                         species="hsapiens",
#                         target="HGNC",
#                         addCol=F)
# 
# setwd("/Users/madhavmakkena/Desktop/RNAseq_Analysis")
# write.csv(count_raw_dirty, "count_raw_dirty.csv", row.names = T)
# count_raw_dirty <- read.csv("count_raw_dirty.csv", header = TRUE, row.names = 1)

setwd("/Users/madhavmakkena/Desktop/RNAseq_Analysis")
# cleaning count_raw_dirty to exclude genes with 0 counts in ALL samples
# count_raw <- count_raw_dirty[rowSums(count_raw_dirty) != 0, ]
# write.csv(count_raw, "count_raw.csv", row.names = T)
count_raw <- read.csv("count_raw.csv", header = TRUE, row.names = 1)

# condition setup using the individual count_raw files
# setwd("/Users/madhavmakkena/Desktop/RNAseq_Analysis/Counts")
# sample_type <- trimws(count_files, whitespace = "_.*")
# sample_name <- trimws(count_files, whitespace = ".csv")
# condition <- data.frame(sample_name, sample_type)
# condition <- column_to_rownames(condition, var = "sample_name")
# setwd("/Users/madhavmakkena/Desktop/RNAseq_Analysis")
# write.csv(condition, "condition.csv", row.names = T)
condition <- read.csv("condition.csv", header = TRUE, row.names = 1)

# checking if the samples names match between count_raw and condition (including order)
rownames(condition) == colnames(count_raw)

# removing non necessary data
# rm(count_combined)
# rm(count_rownames_ENSG)
# rm(count_files)
# rm(sample_name)
# rm(sample_type)
# rm(test)
# rm(non_ENSG_count)
# rm(count_lists)
# rm(count_raw_dirty)