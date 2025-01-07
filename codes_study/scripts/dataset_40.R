#####
datasetID=40
# Functions, dataset info, packages
source("scripts/script_preparation.R")
#######

#######
# NET files
# NETa
NETa <- read_excel("data/40_BASSET_WAN_plant_seed.predator.xlsx")
NETa=as.matrix(NETa)
row.names(NETa)<-NETa[,1]
NETa=NETa[,-1]
mode(NETa)<-"numeric"
NETa=NETa[rowSums(NETa)>0,colSums(NETa)>0]
NETa=t(NETa)

# NETb
NETb <- read_excel("data/40_BASSET_WAN_seed.predator_parasitoid.xlsx")
NETb=as.matrix(NETb)
row.names(NETb)<-NETb[,1]
NETb=NETb[,-1]
mode(NETb)<-"numeric"
NETb=NETb[rowSums(NETb)>0,colSums(NETb)>0]
# common group in rows

## Checking info
sum(dim(NETa))==dataset_info$`Species in network #1`
sum(dim(NETb))==dataset_info$`Species in network #2`
length(unique(c(colnames(NETa),rownames(NETa),colnames(NETb),rownames(NETb))))==dataset_info$`total species`
any(is.element(rownames(NETa), rownames(NETb)))
!any(is.element(rownames(NETa), colnames(NETb)))
!any(is.element(colnames(NETa), rownames(NETb)))
!any(is.element(colnames(NETa), colnames(NETb)))
!(any(duplicated(rownames(NETa)))|any(duplicated(colnames(NETa))))
!(any(duplicated(rownames(NETb)))|any(duplicated(colnames(NETb))))
!any(NETa%%1!=0)&!any(NETb%%1!=0)
!any(c(rowSums(NETa),colSums(NETa),rowSums(NETb), colSums(NETb))==0)
hist(NETa[NETa!=0],breaks = 30)
hist(NETb[NETb!=0],breaks = 30)
dev.off()

# Parameters analysis
NCores=5
N_null_topology=1000
N_null_congruence=1000
N_null_congruence2=1000

source("scripts/script_analysis.R",echo = T)