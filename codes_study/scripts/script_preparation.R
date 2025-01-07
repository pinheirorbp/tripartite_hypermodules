library(stringr)
library(readxl)
library(bipartite)
library(doSNOW)
# Functions
file.sources = list.files(path = "base_functions/",pattern = ".R",full.names = T)
sapply(file.sources,source)
# Dataset info
database_info <- read_excel("data/FINAL_DATABASE.xlsx")
x=str_split_fixed(database_info$Partitions,pattern = "-",n = 3)
database_info$level1=x[,1]
database_info$level2=x[,2]
database_info$level3=x[,3]
dataset_info=database_info[database_info$`Data set`==datasetID,]
#
rm(file.sources,x,database_info)
