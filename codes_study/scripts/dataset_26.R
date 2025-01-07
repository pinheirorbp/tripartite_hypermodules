#####
datasetID=26
# Functions, dataset info, packages
source("scripts/script_preparation.R")
#######

#######
# NET files
# NETa
NETa <- read_excel("data/26_MULLER_plant_herbivore.xlsx")
NETa=as.matrix(NETa)
row.names(NETa)<-NETa[,1]
NETa=NETa[,-1]
mode(NETa)<-"numeric"
NETa=t(NETa) # common group in rows
NETa[NETa>1]=1 # there was cell value=2 and =3 (the matrix should be binary)

# NETb
NETb <- read_excel("data/26_MULLER_herbivore_parasitoid.xlsx")
NETb=as.matrix(NETb)
row.names(NETb)<-NETb[,1]
NETb=NETb[,-1]
mode(NETb)<-"numeric"
# common group in rows

# Parameters analysis
NCores=5
N_null_topology=1000
N_null_congruence=1000
N_null_congruence2=1000

source("scripts/script_analysis.R",echo = T)