#####
datasetID=14
# Functions, dataset info, packages
source("scripts/script_preparation.R")
#######

# NETa : Plant-Herbivore
NETa=read.table("data/14_BRISTOL_plant_herbivore.csv", sep=",", h=T, row.names = 1)
NETa[is.na(NETa)]=0
NETa=as.matrix(NETa)
NETa=NETa[rowSums(NETa)>0,colSums(NETa)>0] # removing empty rows / columns
NETa=t(NETa)
# common group (herbivores) in rows

# NETb : Herbivore-Parasitoid
NETb=read.table("data/14_BRISTOL_herbivore_parasitoid.csv", sep=",", h=T, row.names = 1)
NETb[is.na(NETb)]=0
NETb=as.matrix(NETb)
NETb=NETb[rowSums(NETb)>0,colSums(NETb)>0]

# Parameters analysis
NCores=5
N_null_topology=1000
N_null_congruence=1000
N_null_congruence2=1000

source("scripts/script_analysis.R",echo = T)