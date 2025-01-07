#####
datasetID=38
# Functions, dataset info, packages
source("scripts/script_preparation.R")
#######

#######
# NET files
# NETa
NETa <- read_excel("data/38_BASSET_KHC_plant_pulp.feeder.xlsx")
NETa=as.matrix(NETa)
row.names(NETa)<-NETa[,1]
NETa=NETa[,-1]
mode(NETa)<-"numeric"
NETa=NETa[rowSums(NETa)>0,colSums(NETa)>0]
NETa=t(NETa)

# NETb
NETb <- read_excel("data/38_BASSET_KHC_pulp.feeder_parasitoid.xlsx")
NETb=as.matrix(NETb)
row.names(NETb)<-NETb[,1]
NETb=NETb[,-1]
mode(NETb)<-"numeric"
NETb=NETb[rowSums(NETb)>0,colSums(NETb)>0]
# common group in rows

# Parameters analysis
NCores=5
N_null_topology=1000
N_null_congruence=1000
N_null_congruence2=1000

source("scripts/script_analysis.R",echo = T)