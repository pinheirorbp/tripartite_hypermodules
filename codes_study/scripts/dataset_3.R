datasetID=3
# Functions, dataset info, packages
source("scripts/script_preparation.R")
#######
# NET files
# NETa
NETa=read.table("data/03_AZORES_plant_herbivore.csv", sep=",", h=T, row.names = 1)
NETa[is.na(NETa)]=0
NETa=as.matrix(NETa)
NETa=NETa[rowSums(NETa)>0,colSums(NETa)>0] # removing empty rows / columns
NETa=t(NETa)
# common group in rows

# NETb
NETb=read.table("data/03_AZORES_herbivore_parasitoid.csv", sep=",", h=T, row.names = 1)
NETb[is.na(NETb)]=0
NETb=as.matrix(NETb)
NETb=NETb[rowSums(NETb)>0,colSums(NETb)>0]
# removing empty rows / columns
# common group in rows
## Problem with name of species that include (???)
rownames(NETb)[6:9]=c("Herbivore1","Herbivore2","Herbivore3","Herbivore4")
rownames(NETa)[21:24]=c("Herbivore1","Herbivore2","Herbivore3","Herbivore4")

# Parameters analysis
NCores=5
N_null_topology=1000
N_null_congruence=1000
N_null_congruence2=1000

source("scripts/script_analysis.R",echo = T)