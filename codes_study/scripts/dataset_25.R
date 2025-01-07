#####
datasetID=25
# Functions, dataset info, packages
source("scripts/script_preparation.R")
#######

#######
# NET files
# NETa
NETa <- read_excel("data/25_COSTA_RICA_plant_herbivore.xlsx")
NETa=as.matrix(NETa)
row.names(NETa)<-NETa[,1]
NETa=NETa[,-1]
mode(NETa)<-"numeric"
NETa=NETa[rowSums(NETa)>0,colSums(NETa)>0] # removing empty rows / columns
NETa=t(NETa) # common group in rows
NETa[NETa>1]=1 # there was one cell value=2 (the matrix should be binary)

# NETb
NETb <- read_excel("data/25_COSTA_RICA_herbivore_parasitoid.xlsx")
NETb=as.matrix(NETb)
row.names(NETb)<-NETb[,1]
NETb=NETb[,-1]
mode(NETb)<-"numeric"
NETb=NETb[rowSums(NETb)>0,colSums(NETb)>0] # removing empty rows / columns
# common group in rows

# Parameters analysis
NCores=5
N_null_topology=1000
N_null_congruence=1000
N_null_congruence2=1000

source("scripts/script_analysis.R",echo = T)