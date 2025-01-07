#####
datasetID=32
# Functions, dataset info, packages
source("scripts/script_preparation.R")
#######

#######
# NET files
# NETa
NETa <- read_excel("data/32_SHINOHARA_plant_herbivore.xlsx")
NETa=as.matrix(NETa)
row.names(NETa)<-paste("Plant",NETa[,1], sep="")
colnames(NETa)<-paste("Herbivore",colnames(NETa),sep="")
NETa=NETa[,-1]
mode(NETa)<-"numeric"
NETa=NETa[rowSums(NETa)>0,colSums(NETa)>0]

# NETb
NETb <- read_excel("data/32_SHINOHARA_plant_flower_visitor.xlsx")
NETb=as.matrix(NETb)
row.names(NETb)<-paste("Plant",NETb[,1], sep="")
colnames(NETb)<-paste("Pollinator",colnames(NETb),sep="")
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