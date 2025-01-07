#####
datasetID=10
# Functions, dataset info, packages
source("scripts/script_preparation.R")
#######
# NET files
# NETa
NETa=read.table("data/10_COIMBRA_plant_disperser.csv", sep=",", h=T, row.names = 1)
NETa[is.na(NETa)]=0
NETa=as.matrix(NETa)
NETa=NETa[rowSums(NETa)>0,colSums(NETa)>0] # removing empty rows / columns
NETa=t(NETa)
# common group in rows

# NETb
NETb=read.table("data/10_COIMBRA_disperser_parasite.csv", sep=",", h=T, row.names = 1)
NETb[is.na(NETb)]=0
NETb=as.matrix(NETb)
NETb=NETb[rowSums(NETb)>0,colSums(NETb)>0]
# removing empty rows / columns
# common group in rows

# Parameters analysis
NCores=5
N_null_topology=1000
N_null_congruence=1000
N_null_congruence2=1000

source("scripts/script_analysis.R",echo = T)