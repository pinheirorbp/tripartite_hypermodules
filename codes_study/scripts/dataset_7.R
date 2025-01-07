#####
datasetID=7
# Functions, dataset info, packages
source("scripts/script_preparation.R")
#######
# NET files
# NETa
NETa=read.table("data/07_AB.txt", sep="\t", h=T, row.names = 1)
NETa[is.na(NETa)]=0
NETa=as.matrix(NETa)
NETa=NETa[rowSums(NETa)>0,colSums(NETa)>0] # removing empty rows / columns
# common group in rows

# NETb
NETb=read.table("data/07_NORWOOD_plant_pollinator.csv", sep=",", h=T, row.names = 1)
NETb[is.na(NETb)]=0
NETb=as.matrix(NETb)
NETb=NETb[rowSums(NETb)>0,colSums(NETb)>0]
# common group in rows

# Log Transformation
NETa=log(NETa+1, base=2)
NETb=log(NETb+1, base=2)
# Whole numbers -> next higher whole number
NETa=round(NETa, digits = 0)
NETb=round(NETb, digits = 0)

# Parameters analysis
NCores=5
N_null_topology=1000
N_null_congruence=1000
N_null_congruence2=1000

source("scripts/script_analysis.R",echo = T)