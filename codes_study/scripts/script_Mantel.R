library(readxl)
library(stringr)
library(rdiversity)
library(vegan)
MantelResults=data.frame(TAXON_HM_r=numeric(0),TAXON_HM_p=numeric(0),DISTR_HM_r=numeric(0),DISTR_HM_p=numeric(0), nspecies=numeric(0), permutations=numeric(0))
Manteltests=list()
# # # # # # #
### DATASET 3 ####
datasetID=3
Manteltests$dataset3=list()
load(paste("files_results/dataset_",datasetID,".RData",sep=""))
# # # # # # #
## partition A ##
Manteltests$dataset3$partitionA=list()
species=colnames(NETa)[paste(MODNETa$Col_labels,"a",sep="")%in%rownames(CONGRUENCE$M1M2)]
# GEO
GEOdata= read_excel("spatial_taxonomic_data/03_AZORES_plantsxsite.xlsx")
GEOdata=as.matrix(GEOdata)
GEOdata[is.na(GEOdata)]=0
rownames(GEOdata)=GEOdata[,1]
GEOdata=GEOdata[,-1]
mode(GEOdata)<-"numeric"
# 2 species missing in GEOdata-> removed from the analyses
species=species[species%in%rownames(GEOdata)]
#
GEOdata=GEOdata[match(species,rownames(GEOdata)),]
# TAXONOMY
TAXOdata= read.csv("spatial_taxonomic_data/03_AZORES_plant taxonomy.csv")
TAXOdata=as.matrix(TAXOdata)
TAXOdata=TAXOdata[match(species,TAXOdata[,1]),]
rownames(TAXOdata)=TAXOdata[,1]
# Hypermodules (0= same, 1= different)
MODa=paste(MODNETa$Col_labels[match(species,colnames(NETa))],"a", sep="")
HYPERMOD=data.frame(row.names = species)
for (M in 1:nrow(CONGRUENCE$M1M2)){
  HYPERMOD$NETa[MODa==rownames(CONGRUENCE$M1M2)[M]]=unique(CONGRUENCE$M1M2[M,CONGRUENCE$M1M2[M,]!=0])
}
# Distances
TAXODIST=tax2dist(TAXOdata,tax_distance = c(Species=0, Genus=1, Family=2, Order=3, Kingdom=4))
TAXODIST=TAXODIST@distance
GEODIST=vegdist(GEOdata,method = "jaccard",binary = T)
HYPERMOD_DIST=matrix(NA, nrow=length(species), ncol=length(species), dimnames = list(species,species))
for (SP1 in 1:nrow(HYPERMOD)){
  for(SP2 in 1: nrow(HYPERMOD)){
    if(HYPERMOD$NETa[SP1]==HYPERMOD$NETa[SP2]){
      HYPERMOD_DIST[SP1,SP2]=0
    }else{HYPERMOD_DIST[SP1,SP2]=1}
  }
}
# Mantel Tests
Manteltests$dataset3$partitionA$TAXON_HM=
  mantel.partial(xdis=TAXODIST,ydis=HYPERMOD_DIST,zdis=GEODIST,method = "spearman",permutations = 10)
Manteltests$dataset3$partitionA$DISTR_HM= mantel.partial(xdis=GEODIST,ydis=HYPERMOD_DIST,zdis=TAXODIST,method = "spearman",permutations = 10)
## Table of results
MantelResults=rbind(MantelResults, data.frame(row.names = "dataset3_partitionA",
  TAXON_HM_r=Manteltests$dataset3$partitionA$TAXON_HM$statistic, 
  TAXON_HM_p=Manteltests$dataset3$partitionA$TAXON_HM$signif, 
  DISTR_HM_r=Manteltests$dataset3$partitionA$DISTR_HM$statistic,
  DISTR_HM_p=Manteltests$dataset3$partitionA$DISTR_HM$signif, 
  nspecies=length(species),
  permutations=Manteltests$dataset3$partitionA$TAXON_HM$permutations))
# # # # # # #
## partition B ##
Manteltests$dataset3$partitionB=list()
species=rownames(NETa)[rownames(NETa)%in%rownames(NETb)]
# GEO
GEOdata= read_excel("spatial_taxonomic_data/03_AZORES_herbivoresxsite.xlsx")
GEOdata=as.matrix(GEOdata)
GEOdata[is.na(GEOdata)]=0
rownames(GEOdata)=GEOdata[,1]
GEOdata=GEOdata[,-1]
mode(GEOdata)<-"numeric"
# Species missing in GEOdata-> removed from the analyses
species=species[species%in%rownames(GEOdata)]
#
GEOdata=GEOdata[match(species,rownames(GEOdata)),]
# TAXONOMY
TAXOdata= read_excel("spatial_taxonomic_data/03_AZORES_taxonomy.xlsx")
TAXOdata=as.matrix(TAXOdata)
TAXOdata=TAXOdata[match(species,TAXOdata[,1]),]
rownames(TAXOdata)=TAXOdata[,1]
colnames(TAXOdata)[1]="Species"
# Hypermodules (0= same in both interactions, .5= different in one interaction, 1= different in both interactions)
MODa=paste(MODNETa$Row_labels[match(species,rownames(NETa))],"a", sep="")
MODb=paste(MODNETb$Row_labels[match(species,rownames(NETb))],"b", sep="")
HYPERMOD=data.frame(row.names = species)
for (M in 1:nrow(CONGRUENCE$M1M2)){
  HYPERMOD$NETa[MODa==rownames(CONGRUENCE$M1M2)[M]]=unique(CONGRUENCE$M1M2[M,CONGRUENCE$M1M2[M,]!=0])
}
for (M in 1:ncol(CONGRUENCE$M1M2)){
  HYPERMOD$NETb[MODb==colnames(CONGRUENCE$M1M2)[M]]=unique(CONGRUENCE$M1M2[CONGRUENCE$M1M2[,M]!=0,M])
}
# Distances
TAXODIST=tax2dist(TAXOdata,tax_distance = c(Species=0, Genus=1, Family=2, Order=3))
TAXODIST=TAXODIST@distance
GEODIST=vegdist(GEOdata,method = "jaccard",binary = T)
HYPERMOD_DIST=matrix(1, nrow=length(species), ncol=length(species), dimnames = list(species,species))
for (SP1 in 1:nrow(HYPERMOD)){
  for(SP2 in 1: nrow(HYPERMOD)){
    if(HYPERMOD$NETa[SP1]==HYPERMOD$NETa[SP2]){
      HYPERMOD_DIST[SP1,SP2]=HYPERMOD_DIST[SP1,SP2]-.5}
    if(HYPERMOD$NETb[SP1]==HYPERMOD$NETb[SP2]){
      HYPERMOD_DIST[SP1,SP2]=HYPERMOD_DIST[SP1,SP2]-.5}
  }
}
# Mantel Tests
Manteltests$dataset3$partitionB$TAXON_HM=
  mantel.partial(xdis=TAXODIST,ydis=HYPERMOD_DIST,zdis=GEODIST,method = "spearman",permutations = 10)
Manteltests$dataset3$partitionB$DISTR_HM= mantel.partial(xdis=GEODIST,ydis=HYPERMOD_DIST,zdis=TAXODIST,method = "spearman",permutations = 10)
## Table of results
MantelResults=rbind(MantelResults, data.frame(row.names = "dataset3_partitionB",
  TAXON_HM_r=Manteltests$dataset3$partitionB$TAXON_HM$statistic, 
  TAXON_HM_p=Manteltests$dataset3$partitionB$TAXON_HM$signif, 
  DISTR_HM_r=Manteltests$dataset3$partitionB$DISTR_HM$statistic,
  DISTR_HM_p=Manteltests$dataset3$partitionB$DISTR_HM$signif, 
  nspecies=length(species),
  permutations=Manteltests$dataset3$partitionB$TAXON_HM$permutations))
# # # # # # #
## partition C ##
Manteltests$dataset3$partitionC=list()
species=colnames(NETb)[paste(MODNETb$Col_labels,"b",sep="")%in%colnames(CONGRUENCE$M1M2)]
# GEO - Unavailable data
# TAXONOMY
TAXOdata= read.csv("spatial_taxonomic_data/03_AZORES_parasitoid taxonomy.csv")
TAXOdata=as.matrix(TAXOdata)
TAXOdata=TAXOdata[match(species,TAXOdata[,1]),]
rownames(TAXOdata)=TAXOdata[,1]
# missing data
# Missing family -> remove species
MISSING_FAMILY=is.na(TAXOdata[,3])|TAXOdata[,3]=="NA"|TAXOdata[,3]==""
species=species[!MISSING_FAMILY]
TAXOdata=TAXOdata[!MISSING_FAMILY,]
# Hypermodules (0= same, 1= different)
MODb=paste(MODNETb$Col_labels[match(species,colnames(NETb))],"b", sep="")
HYPERMOD=data.frame(row.names = species)
for (M in 1:ncol(CONGRUENCE$M1M2)){
  HYPERMOD$NETb[MODb==colnames(CONGRUENCE$M1M2)[M]]=unique(CONGRUENCE$M1M2[CONGRUENCE$M1M2[,M]!=0,M])
}
# Distances
TAXODIST=tax2dist(TAXOdata,tax_distance = c(Species=0, Genus=1, Family=2, Order=3, Kingdom=4))
TAXODIST=TAXODIST@distance
HYPERMOD_DIST=matrix(NA, nrow=length(species), ncol=length(species), dimnames = list(species,species))
for (SP1 in 1:nrow(HYPERMOD)){
  for(SP2 in 1: nrow(HYPERMOD)){
    if(HYPERMOD$NETb[SP1]==HYPERMOD$NETb[SP2]){
      HYPERMOD_DIST[SP1,SP2]=0
    }else{HYPERMOD_DIST[SP1,SP2]=1}
  }
}
# Mantel Tests
Manteltests$dataset3$partitionC$TAXON_HM=
  mantel(xdis=TAXODIST,ydis=HYPERMOD_DIST,method = "spearman",permutations = 10)
Manteltests$dataset3$partitionC$DISTR_HM= NULL
## Table of results
MantelResults=rbind(MantelResults, data.frame(row.names = "dataset3_partitionC",
  TAXON_HM_r=Manteltests$dataset3$partitionC$TAXON_HM$statistic, 
  TAXON_HM_p=Manteltests$dataset3$partitionC$TAXON_HM$signif, 
  DISTR_HM_r=NA,
  DISTR_HM_p=NA, 
  nspecies=length(species),
  permutations=Manteltests$dataset3$partitionC$TAXON_HM$permutations))
save(MantelResults,Manteltests, file="files_results/Mantel_results.RData")
# # # # # # #
### DATASET 14 ####
datasetID=14
Manteltests$dataset14=list()
load(paste("files_results/dataset_",datasetID,".RData",sep=""))
# # # # # # # #
## partition A ##
Manteltests$dataset14$partitionA=list()
species=colnames(NETa)[paste(MODNETa$Col_labels,"a",sep="")%in%rownames(CONGRUENCE$M1M2)]
# GEO
GEOdata= read_excel("spatial_taxonomic_data/14_BRISTOL_plantsxsites.xlsx")
GEOdata=as.matrix(GEOdata)
GEOdata[is.na(GEOdata)]=0
GEOdata[,1]=str_replace_all(GEOdata[,1]," ","_")
GEOdata[,1]=str_replace_all(GEOdata[,1],"\\.","")
rownames(GEOdata)=GEOdata[,1]
GEOdata=GEOdata[,-1]
mode(GEOdata)<-"numeric"
GEOdata=GEOdata[match(species,rownames(GEOdata)),]
# TAXONOMY
TAXOdata= read.csv("spatial_taxonomic_data/14_BRISTOL_plant taxonomy.csv")
TAXOdata=as.matrix(TAXOdata)
TAXOdata=TAXOdata[match(species,TAXOdata[,1]),]
rownames(TAXOdata)=TAXOdata[,1]
# missing data
# Missing genus -> consider each as a different genus
TAXOdata[(is.na(TAXOdata[,2])|TAXOdata[,2]=="NA"|TAXOdata[,2]==""),2]=paste("NA",1:sum(is.na(TAXOdata[,2])|TAXOdata[,2]=="NA"|TAXOdata[,2]==""),sep="")
# # Missing family -> remove species
MISSING_FAMILY=is.na(TAXOdata[,3])|TAXOdata[,3]=="NA"|TAXOdata[,3]==""
species=species[!MISSING_FAMILY]
TAXOdata=TAXOdata[!MISSING_FAMILY,]
GEOdata=GEOdata[!MISSING_FAMILY,]
# Hypermodules (0= same, 1= different)
MODa=paste(MODNETa$Col_labels[match(species,colnames(NETa))],"a", sep="")
HYPERMOD=data.frame(row.names = species)
for (M in 1:nrow(CONGRUENCE$M1M2)){
  HYPERMOD$NETa[MODa==rownames(CONGRUENCE$M1M2)[M]]=unique(CONGRUENCE$M1M2[M,CONGRUENCE$M1M2[M,]!=0])
}
# Distances
TAXODIST=tax2dist(TAXOdata,tax_distance = c(Species=0, Genus=1, Family=2, Order=3, Kingdom=4))
TAXODIST=TAXODIST@distance
GEODIST=vegdist(GEOdata,method = "jaccard",binary = T)
HYPERMOD_DIST=matrix(NA, nrow=length(species), ncol=length(species), dimnames = list(species,species))
for (SP1 in 1:nrow(HYPERMOD)){
  for(SP2 in 1: nrow(HYPERMOD)){
    if(HYPERMOD$NETa[SP1]==HYPERMOD$NETa[SP2]){
      HYPERMOD_DIST[SP1,SP2]=0
    }else{HYPERMOD_DIST[SP1,SP2]=1}
  }
}
# Mantel Tests
Manteltests$dataset14$partitionA$TAXON_HM=
  mantel.partial(xdis=TAXODIST,ydis=HYPERMOD_DIST,zdis=GEODIST,method = "spearman",permutations = 10)
Manteltests$dataset14$partitionA$DISTR_HM= mantel.partial(xdis=GEODIST,ydis=HYPERMOD_DIST,zdis=TAXODIST,method = "spearman",permutations = 10)
## Table of results
MantelResults=rbind(MantelResults, data.frame(row.names = "dataset14_partitionA",
  TAXON_HM_r=Manteltests$dataset14$partitionA$TAXON_HM$statistic, 
  TAXON_HM_p=Manteltests$dataset14$partitionA$TAXON_HM$signif, 
  DISTR_HM_r=Manteltests$dataset14$partitionA$DISTR_HM$statistic,
  DISTR_HM_p=Manteltests$dataset14$partitionA$DISTR_HM$signif, 
  nspecies=length(species),
  permutations=Manteltests$dataset14$partitionA$TAXON_HM$permutations))
# # # # # # #
## partition B ##
Manteltests$dataset14$partitionB=list()
species=rownames(NETa)[rownames(NETa)%in%rownames(NETb)]
# GEO
GEOdata= read_excel("spatial_taxonomic_data/14_BRISTOL_herbivorexsites.xlsx")
GEOdata=as.matrix(GEOdata)
GEOdata[is.na(GEOdata)]=0
GEOdata[,1]=str_replace_all(GEOdata[,1]," ","_")
GEOdata[,1]=str_replace_all(GEOdata[,1],"\\.","")
rownames(GEOdata)=GEOdata[,1]
GEOdata=GEOdata[,-1]
mode(GEOdata)<-"numeric"
GEOdata=GEOdata[match(species,rownames(GEOdata)),]
# TAXONOMY
TAXOdata= read_excel("spatial_taxonomic_data/14_BRISTOL_taxonomy.xlsx")
TAXOdata=as.matrix(TAXOdata)
TAXOdata[,1]=str_replace_all(TAXOdata[,1]," ","_")
TAXOdata[,1]=str_replace_all(TAXOdata[,1],"\\.","")
TAXOdata=TAXOdata[match(species,TAXOdata[,1]),]
rownames(TAXOdata)=TAXOdata[,1]
colnames(TAXOdata)[1]="Species"
TAXOdata=cbind(TAXOdata, Class=rep("Insecta",153))
# missing data
# Missing genus -> consider each as a different genus
TAXOdata[(is.na(TAXOdata[,2])|TAXOdata[,2]=="NA"),2]=paste("NA",1:sum(is.na(TAXOdata[,2])|TAXOdata[,2]=="NA"),sep="")
# Missing family -> remove species
MISSING_FAMILY=is.na(TAXOdata[,3])|TAXOdata[,3]=="NA"
species=species[!MISSING_FAMILY]
TAXOdata=TAXOdata[!MISSING_FAMILY,]
GEOdata=GEOdata[!MISSING_FAMILY,]
# Hypermodules (0= same in both interactions, .5= different in one interaction, 1= different in both interactions)
MODa=paste(MODNETa$Row_labels[match(species,rownames(NETa))],"a", sep="")
MODb=paste(MODNETb$Row_labels[match(species,rownames(NETb))],"b", sep="")
HYPERMOD=data.frame(row.names = species)
for (M in 1:nrow(CONGRUENCE$M1M2)){
  HYPERMOD$NETa[MODa==rownames(CONGRUENCE$M1M2)[M]]=unique(CONGRUENCE$M1M2[M,CONGRUENCE$M1M2[M,]!=0])
}
for (M in 1:ncol(CONGRUENCE$M1M2)){
  HYPERMOD$NETb[MODb==colnames(CONGRUENCE$M1M2)[M]]=unique(CONGRUENCE$M1M2[CONGRUENCE$M1M2[,M]!=0,M])
}
# Distances
TAXODIST=tax2dist(TAXOdata,tax_distance = c(Species=0, Genus=1, Family=2, Order=3, Class=4))
TAXODIST=TAXODIST@distance
GEODIST=vegdist(GEOdata,method = "jaccard",binary = T)
HYPERMOD_DIST=matrix(1, nrow=length(species), ncol=length(species), dimnames = list(species,species))
for (SP1 in 1:nrow(HYPERMOD)){
  for(SP2 in 1: nrow(HYPERMOD)){
    if(HYPERMOD$NETa[SP1]==HYPERMOD$NETa[SP2]){
      HYPERMOD_DIST[SP1,SP2]=HYPERMOD_DIST[SP1,SP2]-.5}
    if(HYPERMOD$NETb[SP1]==HYPERMOD$NETb[SP2]){
      HYPERMOD_DIST[SP1,SP2]=HYPERMOD_DIST[SP1,SP2]-.5}
  }
}
# Mantel Tests
Manteltests$dataset14$partitionB$TAXON_HM=
  mantel.partial(xdis=TAXODIST,ydis=HYPERMOD_DIST,zdis=GEODIST,method = "spearman",permutations = 10)
Manteltests$dataset14$partitionB$DISTR_HM= mantel.partial(xdis=GEODIST,ydis=HYPERMOD_DIST,zdis=TAXODIST,method = "spearman",permutations = 10)
## Table of results
MantelResults=rbind(MantelResults, data.frame(row.names = "dataset14_partitionB",
  TAXON_HM_r=Manteltests$dataset14$partitionB$TAXON_HM$statistic, 
  TAXON_HM_p=Manteltests$dataset14$partitionB$TAXON_HM$signif, 
  DISTR_HM_r=Manteltests$dataset14$partitionB$DISTR_HM$statistic,
  DISTR_HM_p=Manteltests$dataset14$partitionB$DISTR_HM$signif, 
  nspecies=length(species),
  permutations=Manteltests$dataset14$partitionB$TAXON_HM$permutations))
# # # # # # #
## partition C ##
Manteltests$dataset14$partitionC=list()
species=colnames(NETb)[paste(MODNETb$Col_labels,"b",sep="")%in%colnames(CONGRUENCE$M1M2)]
species=str_replace_all(species,"\\.","")
species[species=="Apodesmia_similis_"]="Apodesmia_similis"
species=unique(species)
# GEO
GEOdata= read_excel("spatial_taxonomic_data/14_BRISTOL_parasitoidsxsites.xlsx")
GEOdata=as.matrix(GEOdata)
GEOdata[is.na(GEOdata)]=0
# The species names have a lot of special characters.
# I have to remove them in several places.
GEOdata[,1]=str_replace_all(GEOdata[,1]," ","_")
GEOdata[,1]=str_replace_all(GEOdata[,1],"\\.","")
GEOdata[,1]=str_replace_all(GEOdata[,1],"\\(","")
GEOdata[,1]=str_replace_all(GEOdata[,1],"\\)","")
rownames(GEOdata)=GEOdata[,1]
GEOdata=GEOdata[,-1]
mode(GEOdata)<-"numeric"
GEOdata=GEOdata[match(species,rownames(GEOdata)),]
# TAXONOMY
TAXOdata= read.csv("spatial_taxonomic_data/14_BRISTOL_parasitoid taxonomy.csv")
TAXOdata=as.matrix(TAXOdata)
TAXOdata[,1]=str_replace_all(TAXOdata[,1]," ","_")
TAXOdata[,1]=str_replace_all(TAXOdata[,1],"\\.","")
TAXOdata[,1]=str_replace_all(TAXOdata[,1],"\\(","")
TAXOdata[,1]=str_replace_all(TAXOdata[,1],"\\)","")
TAXOdata[TAXOdata[,1]=="Apodesmia_similis_",1]="Apodesmia_similis"
TAXOdata=TAXOdata[match(species,TAXOdata[,1]),]
rownames(TAXOdata)=TAXOdata[,1]
# Hypermodules (0= same, 1= different)
speciesinNETb=colnames(NETb)
speciesinNETb=str_replace_all(speciesinNETb," ","_")
speciesinNETb=str_replace_all(speciesinNETb,"\\.","")
speciesinNETb=str_replace_all(speciesinNETb,"\\(","")
speciesinNETb=str_replace_all(speciesinNETb,"\\)","")
speciesinNETb[speciesinNETb=="Apodesmia_similis_"]="Apodesmia_similis"
MODb=paste(MODNETb$Col_labels[match(species,speciesinNETb)],"b", sep="")
HYPERMOD=data.frame(row.names = species)
for (M in 1:ncol(CONGRUENCE$M1M2)){
  HYPERMOD$NETb[MODb==colnames(CONGRUENCE$M1M2)[M]]=unique(CONGRUENCE$M1M2[CONGRUENCE$M1M2[,M]!=0,M])
}
# Distances
TAXODIST=tax2dist(TAXOdata,tax_distance = c(Species=0, Genus=1, Family=2, Order=3, Kingdom=4))
TAXODIST=TAXODIST@distance
GEODIST=vegdist(GEOdata,method = "jaccard",binary = T)
HYPERMOD_DIST=matrix(NA, nrow=length(species), ncol=length(species), dimnames = list(species,species))
for (SP1 in 1:nrow(HYPERMOD)){
  for(SP2 in 1: nrow(HYPERMOD)){
    if(HYPERMOD$NETb[SP1]==HYPERMOD$NETb[SP2]){
      HYPERMOD_DIST[SP1,SP2]=0
    }else{HYPERMOD_DIST[SP1,SP2]=1}
  }
}
# Mantel Tests
Manteltests$dataset14$partitionC$TAXON_HM=
  mantel.partial(xdis=TAXODIST,ydis=HYPERMOD_DIST,zdis=GEODIST,method = "spearman",permutations = 10)
Manteltests$dataset14$partitionC$DISTR_HM= mantel.partial(xdis=GEODIST,ydis=HYPERMOD_DIST,zdis=TAXODIST,method = "spearman",permutations = 10)
## Table of results
MantelResults=rbind(MantelResults, data.frame(row.names = "dataset14_partitionC",
  TAXON_HM_r=Manteltests$dataset14$partitionC$TAXON_HM$statistic, 
  TAXON_HM_p=Manteltests$dataset14$partitionC$TAXON_HM$signif, 
  DISTR_HM_r=Manteltests$dataset14$partitionC$DISTR_HM$statistic,
  DISTR_HM_p=Manteltests$dataset14$partitionC$DISTR_HM$signif, 
  nspecies=length(species),
  permutations=Manteltests$dataset14$partitionC$TAXON_HM$permutations))
save(MantelResults,Manteltests, file="files_results/Mantel_results.RData")
# # # # # # #
### DATASET 18 ####
datasetID=18
Manteltests$dataset18=list()
load(paste("files_results/dataset_",datasetID,".RData",sep=""))
# # # # # # # #
## partition A ##
Manteltests$dataset18$partitionA=list()
species=colnames(NETa)[paste(MODNETa$Col_labels,"a",sep="")%in%rownames(CONGRUENCE$M1M2)]
# GEO
GEOdata= read.csv("spatial_taxonomic_data/18_plant_distribution_site.csv",row.names = 1)
GEOdata[GEOdata>1]=1
GEOdata=as.matrix(GEOdata)
species=species[species%in%rownames(GEOdata)]
GEOdata=GEOdata[match(species,rownames(GEOdata)),]
# TAXONOMY
TAXOdata= read.csv("spatial_taxonomic_data/18_plant_taxonomy.csv")
TAXOdata=as.matrix(TAXOdata)
TAXOdata=TAXOdata[match(species,TAXOdata[,1]),]
rownames(TAXOdata)=TAXOdata[,1]
# Hypermodules (0= same, 1= different)
MODa=paste(MODNETa$Col_labels[match(species,colnames(NETa))],"a", sep="")
HYPERMOD=data.frame(row.names = species)
for (M in 1:nrow(CONGRUENCE$M1M2)){
  HYPERMOD$NETa[MODa==rownames(CONGRUENCE$M1M2)[M]]=unique(CONGRUENCE$M1M2[M,CONGRUENCE$M1M2[M,]!=0])
}
# Distances
TAXODIST=tax2dist(TAXOdata,tax_distance = c(Species=0, Genus=1, Family=2, Order=3, Kingdom=4))
TAXODIST=TAXODIST@distance
GEODIST=vegdist(GEOdata,method = "jaccard",binary = T)
HYPERMOD_DIST=matrix(NA, nrow=length(species), ncol=length(species), dimnames = list(species,species))
for (SP1 in 1:nrow(HYPERMOD)){
  for(SP2 in 1: nrow(HYPERMOD)){
    if(HYPERMOD$NETa[SP1]==HYPERMOD$NETa[SP2]){
      HYPERMOD_DIST[SP1,SP2]=0
    }else{HYPERMOD_DIST[SP1,SP2]=1}
  }
}
# Mantel Tests
Manteltests$dataset18$partitionA$TAXON_HM=
  mantel.partial(xdis=TAXODIST,ydis=HYPERMOD_DIST,zdis=GEODIST,method = "spearman",permutations = 10)
Manteltests$dataset18$partitionA$DISTR_HM= mantel.partial(xdis=GEODIST,ydis=HYPERMOD_DIST,zdis=TAXODIST,method = "spearman",permutations = 10)
## Table of results
MantelResults=rbind(MantelResults, data.frame(row.names = "dataset18_partitionA",
  TAXON_HM_r=Manteltests$dataset18$partitionA$TAXON_HM$statistic, 
  TAXON_HM_p=Manteltests$dataset18$partitionA$TAXON_HM$signif, 
  DISTR_HM_r=Manteltests$dataset18$partitionA$DISTR_HM$statistic,
  DISTR_HM_p=Manteltests$dataset18$partitionA$DISTR_HM$signif, 
  nspecies=length(species),
  permutations=Manteltests$dataset18$partitionA$TAXON_HM$permutations))
# # # # # # #
## partition B ##
Manteltests$dataset14$partitionB=list()
species=rownames(NETa)[rownames(NETa)%in%rownames(NETb)]
# GEO
GEOdata= read.csv("spatial_taxonomic_data/18_herbivore_distribution_site.csv",row.names = 1)
GEOdata[GEOdata>1]=1
GEOdata=as.matrix(t(GEOdata))
GEOdata=GEOdata[match(species,rownames(GEOdata)),]
# TAXONOMY
TAXOdata= read.csv("spatial_taxonomic_data/18_herbivore_taxonomy.csv")
TAXOdata=as.matrix(TAXOdata)
TAXOdata=TAXOdata[match(species,TAXOdata[,1]),]
rownames(TAXOdata)=TAXOdata[,1]
# Hypermodules (0= same in both interactions, .5= different in one interaction, 1= different in both interactions)
MODa=paste(MODNETa$Row_labels[match(species,rownames(NETa))],"a", sep="")
MODb=paste(MODNETb$Row_labels[match(species,rownames(NETb))],"b", sep="")
HYPERMOD=data.frame(row.names = species)
for (M in 1:nrow(CONGRUENCE$M1M2)){
  HYPERMOD$NETa[MODa==rownames(CONGRUENCE$M1M2)[M]]=unique(CONGRUENCE$M1M2[M,CONGRUENCE$M1M2[M,]!=0])
}
for (M in 1:ncol(CONGRUENCE$M1M2)){
  HYPERMOD$NETb[MODb==colnames(CONGRUENCE$M1M2)[M]]=unique(CONGRUENCE$M1M2[CONGRUENCE$M1M2[,M]!=0,M])
}
# Distances
TAXODIST=tax2dist(TAXOdata,tax_distance = c(Species=0, Genus=1, Family=2, Order=3, Kingdom=4))
TAXODIST=TAXODIST@distance
GEODIST=vegdist(GEOdata,method = "jaccard",binary = T)
HYPERMOD_DIST=matrix(NA, nrow=length(species), ncol=length(species), dimnames = list(species,species))
for (SP1 in 1:nrow(HYPERMOD)){
  for(SP2 in 1: nrow(HYPERMOD)){
    if(HYPERMOD$NETa[SP1]==HYPERMOD$NETa[SP2]){
      HYPERMOD_DIST[SP1,SP2]=0
    }else{HYPERMOD_DIST[SP1,SP2]=1}
  }
}
# Mantel Tests
Manteltests$dataset18$partitionB$TAXON_HM=
  mantel.partial(xdis=TAXODIST,ydis=HYPERMOD_DIST,zdis=GEODIST,method = "spearman",permutations = 10)
Manteltests$dataset18$partitionB$DISTR_HM= mantel.partial(xdis=GEODIST,ydis=HYPERMOD_DIST,zdis=TAXODIST,method = "spearman",permutations = 10)
## Table of results
MantelResults=rbind(MantelResults, data.frame(row.names = "dataset18_partitionB",
  TAXON_HM_r=Manteltests$dataset18$partitionB$TAXON_HM$statistic, 
  TAXON_HM_p=Manteltests$dataset18$partitionB$TAXON_HM$signif, 
  DISTR_HM_r=Manteltests$dataset18$partitionB$DISTR_HM$statistic,
  DISTR_HM_p=Manteltests$dataset18$partitionB$DISTR_HM$signif, 
  nspecies=length(species),
  permutations=Manteltests$dataset18$partitionB$TAXON_HM$permutations))
# # # # # # #
## partition C ##
Manteltests$dataset18$partitionC=list()
species=colnames(NETb)[paste(MODNETb$Col_labels,"b",sep="")%in%colnames(CONGRUENCE$M1M2)]
# GEO
GEOdata= read.csv("spatial_taxonomic_data/18_parasitoid_distribution_site.csv",row.names = 1)
GEOdata[GEOdata>1]=1
GEOdata=as.matrix(t(GEOdata))
GEOdata=GEOdata[match(species,rownames(GEOdata)),]
# TAXONOMY
TAXOdata= read.csv("spatial_taxonomic_data/18_parasitoid_taxonomy.csv")
TAXOdata=as.matrix(TAXOdata)
TAXOdata=TAXOdata[match(species,TAXOdata[,1]),]
rownames(TAXOdata)=TAXOdata[,1]
# Hypermodules (0= same, 1= different)
MODb=paste(MODNETb$Col_labels[match(species,colnames(NETb))],"b", sep="")
HYPERMOD=data.frame(row.names = species)
for (M in 1:ncol(CONGRUENCE$M1M2)){
  HYPERMOD$NETb[MODb==colnames(CONGRUENCE$M1M2)[M]]=unique(CONGRUENCE$M1M2[CONGRUENCE$M1M2[,M]!=0,M])
}
# Distances
TAXODIST=tax2dist(TAXOdata,tax_distance = c(Species=0, Genus=1, Family=2, Order=3, Kingdom=4))
TAXODIST=TAXODIST@distance
GEODIST=vegdist(GEOdata,method = "jaccard",binary = T)
HYPERMOD_DIST=matrix(NA, nrow=length(species), ncol=length(species), dimnames = list(species,species))
for (SP1 in 1:nrow(HYPERMOD)){
  for(SP2 in 1: nrow(HYPERMOD)){
    if(HYPERMOD$NETb[SP1]==HYPERMOD$NETb[SP2]){
      HYPERMOD_DIST[SP1,SP2]=0
    }else{HYPERMOD_DIST[SP1,SP2]=1}
  }
}
# Mantel Tests
Manteltests$dataset18$partitionC$TAXON_HM=
  mantel.partial(xdis=TAXODIST,ydis=HYPERMOD_DIST,zdis=GEODIST,method = "spearman",permutations = 10)
Manteltests$dataset18$partitionC$DISTR_HM= mantel.partial(xdis=GEODIST,ydis=HYPERMOD_DIST,zdis=TAXODIST,method = "spearman",permutations = 10)
## Table of results
MantelResults=rbind(MantelResults, data.frame(row.names = "dataset18_partitionC",
  TAXON_HM_r=Manteltests$dataset18$partitionC$TAXON_HM$statistic, 
  TAXON_HM_p=Manteltests$dataset18$partitionC$TAXON_HM$signif, 
  DISTR_HM_r=Manteltests$dataset18$partitionC$DISTR_HM$statistic,
  DISTR_HM_p=Manteltests$dataset18$partitionC$DISTR_HM$signif, 
  nspecies=length(species),
  permutations=Manteltests$dataset18$partitionC$TAXON_HM$permutations))
### DATASET 19 ####
datasetID=19
Manteltests$dataset19=list()
load(paste("files_results/dataset_",datasetID,".RData",sep=""))
# # # # # # #
## partition A ##
Manteltests$dataset19$partitionA=list()
species=colnames(NETa)[paste(MODNETa$Col_labels,"a",sep="")%in%rownames(CONGRUENCE$M1M2)]
# GEO
GEOdata= read_excel("spatial_taxonomic_data/CORDOBA1_plantsxsites.xlsx")
GEOdata=as.matrix(GEOdata)
GEOdata[is.na(GEOdata)]=0
rownames(GEOdata)=GEOdata[,1]
GEOdata=GEOdata[,-1]
mode(GEOdata)<-"numeric"
GEOdata=GEOdata[match(species,rownames(GEOdata)),]
# TAXONOMY
TAXOdata= read_excel("spatial_taxonomic_data/CORDOBA1_plant_taxonomy.xlsx")
TAXOdata=as.matrix(TAXOdata)
TAXOdata=TAXOdata[match(species,TAXOdata[,1]),]
rownames(TAXOdata)=TAXOdata[,1]
# Hypermodules (0= same, 1= different)
MODa=paste(MODNETa$Col_labels[match(species,colnames(NETa))],"a", sep="")
HYPERMOD=data.frame(row.names = species)
for (M in 1:nrow(CONGRUENCE$M1M2)){
  HYPERMOD$NETa[MODa==rownames(CONGRUENCE$M1M2)[M]]=unique(CONGRUENCE$M1M2[M,CONGRUENCE$M1M2[M,]!=0])
}
# Distances
TAXODIST=tax2dist(TAXOdata,tax_distance = c(Species=0, Genus=1, Family=2, Order=3, Kingdom=4))
TAXODIST=TAXODIST@distance
GEODIST=vegdist(GEOdata,method = "jaccard",binary = T)
HYPERMOD_DIST=matrix(NA, nrow=length(species), ncol=length(species), dimnames = list(species,species))
for (SP1 in 1:nrow(HYPERMOD)){
  for(SP2 in 1: nrow(HYPERMOD)){
    if(HYPERMOD$NETa[SP1]==HYPERMOD$NETa[SP2]){
      HYPERMOD_DIST[SP1,SP2]=0
    }else{HYPERMOD_DIST[SP1,SP2]=1}
  }
}
# Mantel Tests
Manteltests$dataset19$partitionA$TAXON_HM=
  mantel.partial(xdis=TAXODIST,ydis=HYPERMOD_DIST,zdis=GEODIST,method = "spearman",permutations = 10)
Manteltests$dataset19$partitionA$DISTR_HM= mantel.partial(xdis=GEODIST,ydis=HYPERMOD_DIST,zdis=TAXODIST,method = "spearman",permutations = 10)
## Table of results
MantelResults=rbind(MantelResults, data.frame(row.names = "dataset19_partitionA",
  TAXON_HM_r=Manteltests$dataset19$partitionA$TAXON_HM$statistic, 
  TAXON_HM_p=Manteltests$dataset19$partitionA$TAXON_HM$signif, 
  DISTR_HM_r=Manteltests$dataset19$partitionA$DISTR_HM$statistic,
  DISTR_HM_p=Manteltests$dataset19$partitionA$DISTR_HM$signif, 
  nspecies=length(species),
  permutations=Manteltests$dataset19$partitionA$TAXON_HM$permutations))
# # # # # # # #
## partition B ##
Manteltests$dataset19$partitionB=list()
species=rownames(NETa)[rownames(NETa)%in%rownames(NETb)]
# GEO
GEOdata= read_excel("spatial_taxonomic_data/CORDOBA1_herbivoresxsites.xlsx")
GEOdata=as.matrix(GEOdata)
GEOdata[is.na(GEOdata)]=0
rownames(GEOdata)=GEOdata[,1]
GEOdata=GEOdata[,-1]
mode(GEOdata)<-"numeric"
GEOdata=GEOdata[match(species,rownames(GEOdata)),]
# TAXONOMY
TAXOdata= read_excel("spatial_taxonomic_data/CORDOBA1_taxonomy.xlsx")
TAXOdata=as.matrix(TAXOdata)
TAXOdata=TAXOdata[match(species,TAXOdata[,1]),]
rownames(TAXOdata)=TAXOdata[,1]
# missing data
# Missing genus -> consider each as a different genus
TAXOdata[(is.na(TAXOdata[,2])|TAXOdata[,2]=="NA"),2]=paste("NA",1:sum(is.na(TAXOdata[,2])|TAXOdata[,2]=="NA"),sep="")
# Missing family -> remove species
MISSING_FAMILY=is.na(TAXOdata[,3])|TAXOdata[,3]=="NA"
species=species[!MISSING_FAMILY]
TAXOdata=TAXOdata[!MISSING_FAMILY,]
GEOdata=GEOdata[!MISSING_FAMILY,]
TAXOdata=cbind(TAXOdata, Class=rep("Insecta",100))
# Hypermodules (0= same in both interactions, .5= different in one interaction, 1= different in both interactions)
MODa=paste(MODNETa$Row_labels[match(species,rownames(NETa))],"a", sep="")
MODb=paste(MODNETb$Row_labels[match(species,rownames(NETb))],"b", sep="")
HYPERMOD=data.frame(row.names = species)
for (M in 1:nrow(CONGRUENCE$M1M2)){
  HYPERMOD$NETa[MODa==rownames(CONGRUENCE$M1M2)[M]]=unique(CONGRUENCE$M1M2[M,CONGRUENCE$M1M2[M,]!=0])
}
for (M in 1:ncol(CONGRUENCE$M1M2)){
  HYPERMOD$NETb[MODb==colnames(CONGRUENCE$M1M2)[M]]=unique(CONGRUENCE$M1M2[CONGRUENCE$M1M2[,M]!=0,M])
}
# Distances
TAXODIST=tax2dist(TAXOdata,tax_distance = c(Species=0, Genus=1, Family=2, Order=3, Class=4))
TAXODIST=TAXODIST@distance
GEODIST=vegdist(GEOdata,method = "jaccard",binary = T)
HYPERMOD_DIST=matrix(1, nrow=length(species), ncol=length(species), dimnames = list(species,species))
for (SP1 in 1:nrow(HYPERMOD)){
  for(SP2 in 1: nrow(HYPERMOD)){
    if(HYPERMOD$NETa[SP1]==HYPERMOD$NETa[SP2]){
      HYPERMOD_DIST[SP1,SP2]=HYPERMOD_DIST[SP1,SP2]-.5}
    if(HYPERMOD$NETb[SP1]==HYPERMOD$NETb[SP2]){
      HYPERMOD_DIST[SP1,SP2]=HYPERMOD_DIST[SP1,SP2]-.5}
  }
}
# Mantel Tests
Manteltests$dataset19$partitionB$TAXON_HM=
  mantel.partial(xdis=TAXODIST,ydis=HYPERMOD_DIST,zdis=GEODIST,method = "spearman",permutations = 10)
Manteltests$dataset19$partitionB$DISTR_HM= mantel.partial(xdis=GEODIST,ydis=HYPERMOD_DIST,zdis=TAXODIST,method = "spearman",permutations = 10)
## Table of results
MantelResults=rbind(MantelResults, data.frame(row.names = "dataset19_partitionB",
  TAXON_HM_r=Manteltests$dataset19$partitionB$TAXON_HM$statistic, 
  TAXON_HM_p=Manteltests$dataset19$partitionB$TAXON_HM$signif, 
  DISTR_HM_r=Manteltests$dataset19$partitionB$DISTR_HM$statistic,
  DISTR_HM_p=Manteltests$dataset19$partitionB$DISTR_HM$signif, 
  nspecies=length(species),
  permutations=Manteltests$dataset19$partitionB$TAXON_HM$permutations))
# # # # # # #
## partition C ##
Manteltests$dataset19$partitionC=list()
species=colnames(NETb)[paste(MODNETb$Col_labels,"b",sep="")%in%colnames(CONGRUENCE$M1M2)]
# GEO
GEOdata= read_excel("spatial_taxonomic_data/CORDOBA1_parasitoidsxsites.xlsx")
GEOdata=as.matrix(GEOdata)
GEOdata[is.na(GEOdata)]=0
rownames(GEOdata)=GEOdata[,1]
GEOdata=GEOdata[,-1]
mode(GEOdata)<-"numeric"
GEOdata=GEOdata[match(species,rownames(GEOdata)),]
# TAXONOMY
TAXOdata= read.csv("spatial_taxonomic_data/CORDOBA1_parasitoid taxonomy.csv")
TAXOdata=as.matrix(TAXOdata)
TAXOdata=TAXOdata[match(species,TAXOdata[,1]),]
rownames(TAXOdata)=TAXOdata[,1]
# missing data
# Missing genus -> consider each as a different genus
TAXOdata[(is.na(TAXOdata[,2])|TAXOdata[,2]=="NA"|TAXOdata[,2]==""),2]=paste("NA",1:sum(is.na(TAXOdata[,2])|TAXOdata[,2]=="NA"|TAXOdata[,2]==""),sep="")
# Missing family -> remove species
MISSING_FAMILY=is.na(TAXOdata[,3])|TAXOdata[,3]=="NA"|TAXOdata[,3]==""
species=species[!MISSING_FAMILY]
TAXOdata=TAXOdata[!MISSING_FAMILY,]
GEOdata=GEOdata[!MISSING_FAMILY,]
# Hypermodules (0= same, 1= different)
MODb=paste(MODNETb$Col_labels[match(species,colnames(NETb))],"b", sep="")
HYPERMOD=data.frame(row.names = species)
for (M in 1:ncol(CONGRUENCE$M1M2)){
  HYPERMOD$NETb[MODb==colnames(CONGRUENCE$M1M2)[M]]=unique(CONGRUENCE$M1M2[CONGRUENCE$M1M2[,M]!=0,M])
}
# Distances
TAXODIST=tax2dist(TAXOdata,tax_distance = c(Species=0, Genus=1, Family=2, Order=3, Kingdom=4))
TAXODIST=TAXODIST@distance
GEODIST=vegdist(GEOdata,method = "jaccard",binary = T)
HYPERMOD_DIST=matrix(NA, nrow=length(species), ncol=length(species), dimnames = list(species,species))
for (SP1 in 1:nrow(HYPERMOD)){
  for(SP2 in 1: nrow(HYPERMOD)){
    if(HYPERMOD$NETb[SP1]==HYPERMOD$NETb[SP2]){
      HYPERMOD_DIST[SP1,SP2]=0
    }else{HYPERMOD_DIST[SP1,SP2]=1}
  }
}
# Mantel Tests
Manteltests$dataset19$partitionC$TAXON_HM=
  mantel.partial(xdis=TAXODIST,ydis=HYPERMOD_DIST,zdis=GEODIST,method = "spearman",permutations = 10)
Manteltests$dataset19$partitionC$DISTR_HM= mantel.partial(xdis=GEODIST,ydis=HYPERMOD_DIST,zdis=TAXODIST,method = "spearman",permutations = 10)
## Table of results
MantelResults=rbind(MantelResults, data.frame(row.names = "dataset19_partitionC",
  TAXON_HM_r=Manteltests$dataset19$partitionC$TAXON_HM$statistic, 
  TAXON_HM_p=Manteltests$dataset19$partitionC$TAXON_HM$signif, 
  DISTR_HM_r=Manteltests$dataset19$partitionC$DISTR_HM$statistic,
  DISTR_HM_p=Manteltests$dataset19$partitionC$DISTR_HM$signif, 
  nspecies=length(species),
  permutations=Manteltests$dataset19$partitionC$TAXON_HM$permutations))
save(MantelResults,Manteltests, file="files_results/Mantel_results.RData")
### DATASET 20 ####
datasetID=20
Manteltests$dataset20=list()
load(paste("files_results/dataset_",datasetID,".RData",sep=""))
# # # # # # #
## partition A ##
Manteltests$dataset20$partitionA=list()
species=colnames(NETa)[paste(MODNETa$Col_labels,"a",sep="")%in%rownames(CONGRUENCE$M1M2)]
species[species== "Celtis_ehrenbergiana_"]= "Celtis_ehrenbergiana"
species=unique(species)
# GEO
GEOdata= read_excel("spatial_taxonomic_data/CORDOBA2_plantsxsite.xlsx")
GEOdata=as.matrix(GEOdata)
GEOdata[is.na(GEOdata)]=0
rownames(GEOdata)=GEOdata[,1]
GEOdata=GEOdata[,-1]
mode(GEOdata)<-"numeric"
GEOdata=GEOdata[match(species,rownames(GEOdata)),]
# TAXONOMY
TAXOdata= read_excel("spatial_taxonomic_data/CORDOBA2_plants_taxonomy.xlsx")
TAXOdata=as.matrix(TAXOdata)
TAXOdata=TAXOdata[match(species,TAXOdata[,1]),]
rownames(TAXOdata)=TAXOdata[,1]
colnames(TAXOdata)[1]="Species"
# Hypermodules (0= same, 1= different)
MODa=paste(MODNETa$Col_labels[match(species,colnames(NETa))],"a", sep="")
HYPERMOD=data.frame(row.names = species)
for (M in 1:nrow(CONGRUENCE$M1M2)){
  HYPERMOD$NETa[MODa==rownames(CONGRUENCE$M1M2)[M]]=unique(CONGRUENCE$M1M2[M,CONGRUENCE$M1M2[M,]!=0])
}
# Distances
TAXODIST=tax2dist(TAXOdata,tax_distance = c(Species=0, Genus=1, Family=2, Order=3, Kingdom=4))
TAXODIST=TAXODIST@distance
GEODIST=vegdist(GEOdata,method = "jaccard",binary = T)
HYPERMOD_DIST=matrix(NA, nrow=length(species), ncol=length(species), dimnames = list(species,species))
for (SP1 in 1:nrow(HYPERMOD)){
  for(SP2 in 1: nrow(HYPERMOD)){
    if(HYPERMOD$NETa[SP1]==HYPERMOD$NETa[SP2]){
      HYPERMOD_DIST[SP1,SP2]=0
    }else{HYPERMOD_DIST[SP1,SP2]=1}
  }
}
# Mantel Tests
Manteltests$dataset20$partitionA$TAXON_HM=
  mantel.partial(xdis=TAXODIST,ydis=HYPERMOD_DIST,zdis=GEODIST,method = "spearman",permutations = 10)
Manteltests$dataset20$partitionA$DISTR_HM= mantel.partial(xdis=GEODIST,ydis=HYPERMOD_DIST,zdis=TAXODIST,method = "spearman",permutations = 10)
## Table of results
MantelResults=rbind(MantelResults, data.frame(row.names = "dataset20_partitionA",
  TAXON_HM_r=Manteltests$dataset20$partitionA$TAXON_HM$statistic, 
  TAXON_HM_p=Manteltests$dataset20$partitionA$TAXON_HM$signif, 
  DISTR_HM_r=Manteltests$dataset20$partitionA$DISTR_HM$statistic,
  DISTR_HM_p=Manteltests$dataset20$partitionA$DISTR_HM$signif, 
  nspecies=length(species),
  permutations=Manteltests$dataset20$partitionA$TAXON_HM$permutations))
# # # # # # # #
## partition B ##
Manteltests$dataset20$partitionB=list()
species=rownames(NETa)[rownames(NETa)%in%rownames(NETb)]
# GEO
GEOdata= read_excel("spatial_taxonomic_data/CORDOBA2_herbivoresxsite.xlsx")
GEOdata=as.matrix(GEOdata)
GEOdata[is.na(GEOdata)]=0
rownames(GEOdata)=GEOdata[,1]
GEOdata=GEOdata[,-1]
mode(GEOdata)<-"numeric"
GEOdata=GEOdata[match(species,rownames(GEOdata)),]
# TAXONOMY
TAXOdata= read_excel("spatial_taxonomic_data/CORDOBA2_herbivores_taxonomy.xlsx")
TAXOdata=as.matrix(TAXOdata)
TAXOdata=TAXOdata[match(species,TAXOdata[,1]),]
rownames(TAXOdata)=TAXOdata[,1]
# missing data
# Missing genus -> consider each as a different genus
TAXOdata[(is.na(TAXOdata[,2])|TAXOdata[,2]=="NA"),2]=paste("NA",1:sum(is.na(TAXOdata[,2])|TAXOdata[,2]=="NA"),sep="")
colnames(TAXOdata)=c("Species","Genus","Family","Superfamily","Suborder","Order")
# Hypermodules (0= same in both interactions, .5= different in one interaction, 1= different in both interactions)
MODa=paste(MODNETa$Row_labels[match(species,rownames(NETa))],"a", sep="")
MODb=paste(MODNETb$Row_labels[match(species,rownames(NETb))],"b", sep="")
HYPERMOD=data.frame(row.names = species)
for (M in 1:nrow(CONGRUENCE$M1M2)){
  HYPERMOD$NETa[MODa==rownames(CONGRUENCE$M1M2)[M]]=unique(CONGRUENCE$M1M2[M,CONGRUENCE$M1M2[M,]!=0])
}
for (M in 1:ncol(CONGRUENCE$M1M2)){
  HYPERMOD$NETb[MODb==colnames(CONGRUENCE$M1M2)[M]]=unique(CONGRUENCE$M1M2[CONGRUENCE$M1M2[,M]!=0,M])
}
# Distances
TAXODIST=tax2dist(TAXOdata,tax_distance = c(Species=0, Genus=1, Family=2, Superfamily=3, Suborder=4, Order=5))
TAXODIST=TAXODIST@distance
GEODIST=vegdist(GEOdata,method = "jaccard",binary = T)
HYPERMOD_DIST=matrix(1, nrow=length(species), ncol=length(species), dimnames = list(species,species))
for (SP1 in 1:nrow(HYPERMOD)){
  for(SP2 in 1: nrow(HYPERMOD)){
    if(HYPERMOD$NETa[SP1]==HYPERMOD$NETa[SP2]){
      HYPERMOD_DIST[SP1,SP2]=HYPERMOD_DIST[SP1,SP2]-.5}
    if(HYPERMOD$NETb[SP1]==HYPERMOD$NETb[SP2]){
      HYPERMOD_DIST[SP1,SP2]=HYPERMOD_DIST[SP1,SP2]-.5}
  }
}
# Mantel Tests
Manteltests$dataset20$partitionB$TAXON_HM=
  mantel.partial(xdis=TAXODIST,ydis=HYPERMOD_DIST,zdis=GEODIST,method = "spearman",permutations = 10)
Manteltests$dataset20$partitionB$DISTR_HM= mantel.partial(xdis=GEODIST,ydis=HYPERMOD_DIST,zdis=TAXODIST,method = "spearman",permutations = 10)
## Table of results
MantelResults=rbind(MantelResults, data.frame(row.names = "dataset20_partitionB",
  TAXON_HM_r=Manteltests$dataset20$partitionB$TAXON_HM$statistic, 
  TAXON_HM_p=Manteltests$dataset20$partitionB$TAXON_HM$signif, 
  DISTR_HM_r=Manteltests$dataset20$partitionB$DISTR_HM$statistic,
  DISTR_HM_p=Manteltests$dataset20$partitionB$DISTR_HM$signif, 
  nspecies=length(species),
  permutations=Manteltests$dataset20$partitionB$TAXON_HM$permutations))
# # # # # # #
## partition C ##
Manteltests$dataset20$partitionC=list()
species=colnames(NETb)[paste(MODNETb$Col_labels,"b",sep="")%in%colnames(CONGRUENCE$M1M2)]
# GEO
GEOdata= read_excel("spatial_taxonomic_data/CORDOBA2_parasitoidsxsite.xlsx")
GEOdata=as.matrix(GEOdata)
GEOdata[is.na(GEOdata)]=0
rownames(GEOdata)=GEOdata[,1]
GEOdata=GEOdata[,-1]
mode(GEOdata)<-"numeric"
# names in GEOdata (complete) in different format from networks (reduced)
X=str_split(rownames(GEOdata),pattern = "_")
rownames(GEOdata)=sapply(X, function (X){paste (substr(X[1],1,3),"_",substr(X[2],1,3), sep="")})
GEOdata=GEOdata[match(species,rownames(GEOdata)),]
# TAXONOMY
TAXOdata= read.csv("spatial_taxonomic_data/CORDOBA2_parasitoids taxonomy.csv")
TAXOdata=as.matrix(TAXOdata)
# names in TAXOdata (complete) in different format from networks (reduced)
X=str_split(TAXOdata[,1],pattern = "_")
TAXOdata[,1]=sapply(X, function (X){paste (substr(X[1],1,3),"_",substr(X[2],1,3), sep="")})
TAXOdata=TAXOdata[match(species,TAXOdata[,1]),]
rownames(TAXOdata)=TAXOdata[,1]
# missing data
# Missing genus -> consider each as a different genus
TAXOdata[(is.na(TAXOdata[,2])|TAXOdata[,2]=="NA"|TAXOdata[,2]==""),2]=paste("NA",1:sum(is.na(TAXOdata[,2])|TAXOdata[,2]=="NA"|TAXOdata[,2]==""),sep="")
# Missing family -> remove species
MISSING_FAMILY=is.na(TAXOdata[,3])|TAXOdata[,3]=="NA"|TAXOdata[,3]==""
species=species[!MISSING_FAMILY]
TAXOdata=TAXOdata[!MISSING_FAMILY,]
GEOdata=GEOdata[!MISSING_FAMILY,]
# Hypermodules (0= same, 1= different)
MODb=paste(MODNETb$Col_labels[match(species,colnames(NETb))],"b", sep="")
HYPERMOD=data.frame(row.names = species)
for (M in 1:ncol(CONGRUENCE$M1M2)){
  HYPERMOD$NETb[MODb==colnames(CONGRUENCE$M1M2)[M]]=unique(CONGRUENCE$M1M2[CONGRUENCE$M1M2[,M]!=0,M])
}
# Distances
TAXODIST=tax2dist(TAXOdata,tax_distance = c(Species=0, Genus=1, Family=2, Order=3, Kingdom=4))
TAXODIST=TAXODIST@distance
GEODIST=vegdist(GEOdata,method = "jaccard",binary = T)
HYPERMOD_DIST=matrix(NA, nrow=length(species), ncol=length(species), dimnames = list(species,species))
for (SP1 in 1:nrow(HYPERMOD)){
  for(SP2 in 1: nrow(HYPERMOD)){
    if(HYPERMOD$NETb[SP1]==HYPERMOD$NETb[SP2]){
      HYPERMOD_DIST[SP1,SP2]=0
    }else{HYPERMOD_DIST[SP1,SP2]=1}
  }
}
# Mantel Tests
Manteltests$dataset20$partitionC$TAXON_HM=
  mantel.partial(xdis=TAXODIST,ydis=HYPERMOD_DIST,zdis=GEODIST,method = "spearman",permutations = 10)
Manteltests$dataset20$partitionC$DISTR_HM= mantel.partial(xdis=GEODIST,ydis=HYPERMOD_DIST,zdis=TAXODIST,method = "spearman",permutations = 10)
## Table of results
MantelResults=rbind(MantelResults, data.frame(row.names = "dataset20_partitionC",
  TAXON_HM_r=Manteltests$dataset20$partitionC$TAXON_HM$statistic, 
  TAXON_HM_p=Manteltests$dataset20$partitionC$TAXON_HM$signif, 
  DISTR_HM_r=Manteltests$dataset20$partitionC$DISTR_HM$statistic,
  DISTR_HM_p=Manteltests$dataset20$partitionC$DISTR_HM$signif, 
  nspecies=length(species),
  permutations=Manteltests$dataset20$partitionC$TAXON_HM$permutations))
save(MantelResults,Manteltests, file="files_results/Mantel_results.RData")
# # # # # # #
### DATASET 21 ####
datasetID=21
Manteltests$dataset21=list()
load(paste("files_results/dataset_",datasetID,".RData",sep=""))
# # # # # # #
## partition A ##
Manteltests$dataset21$partitionA=list()
species=colnames(NETa)[paste(MODNETa$Col_labels,"a",sep="")%in%rownames(CONGRUENCE$M1M2)]
# GEO
GEOdata= read_excel("spatial_taxonomic_data/CORDOBA2_plant.defenderxsite.xlsx")
GEOdata=as.matrix(GEOdata)
GEOdata[is.na(GEOdata)]=0
rownames(GEOdata)=GEOdata[,1]
GEOdata=GEOdata[,-1]
mode(GEOdata)<-"numeric"
GEOdata=GEOdata[match(species,rownames(GEOdata)),]
# TAXONOMY
TAXOdata= read_excel("spatial_taxonomic_data/CORDOBA2_herbivore.defender_taxonomy.xlsx")
# Missing data that can be infered
TAXOdata=rbind(TAXOdata, data.frame(
  Species=c("Cephalotes_sp1","Cephalotes_sp2","Cephalotes_sp5","Cephalotes_sp6","Dorymyrmex_sp1","Dorymyrmex_sp2","Pheidole_sp1","Solenopsis_parva"),
  Genus=c("Cephalotes","Cephalotes","Cephalotes","Cephalotes","Dorymyrmex","Dorymyrmex","Pheidole","Solenopsis"),
  Family="Formicidae", Order="Hymenoptera", Kingdom="Animalia"
))
TAXOdata=as.matrix(TAXOdata)
TAXOdata=TAXOdata[match(species,TAXOdata[,1]),]
rownames(TAXOdata)=TAXOdata[,1]
# Hypermodules (0= same, 1= different)
MODa=paste(MODNETa$Col_labels[match(species,colnames(NETa))],"a", sep="")
HYPERMOD=data.frame(row.names = species)
for (M in 1:nrow(CONGRUENCE$M1M2)){
  HYPERMOD$NETa[MODa==rownames(CONGRUENCE$M1M2)[M]]=unique(CONGRUENCE$M1M2[M,CONGRUENCE$M1M2[M,]!=0])
}
# Distances
TAXODIST=tax2dist(TAXOdata,tax_distance = c(Species=0, Genus=1, Family=2, Order=3, Kingdom=4))
TAXODIST=TAXODIST@distance
GEODIST=vegdist(GEOdata,method = "jaccard",binary = T)
HYPERMOD_DIST=matrix(NA, nrow=length(species), ncol=length(species), dimnames = list(species,species))
for (SP1 in 1:nrow(HYPERMOD)){
  for(SP2 in 1: nrow(HYPERMOD)){
    if(HYPERMOD$NETa[SP1]==HYPERMOD$NETa[SP2]){
      HYPERMOD_DIST[SP1,SP2]=0
    }else{HYPERMOD_DIST[SP1,SP2]=1}
  }
}
# Mantel Tests
Manteltests$dataset21$partitionA$TAXON_HM=
  mantel.partial(xdis=TAXODIST,ydis=HYPERMOD_DIST,zdis=GEODIST,method = "spearman",permutations = 10)
Manteltests$dataset21$partitionA$DISTR_HM= mantel.partial(xdis=GEODIST,ydis=HYPERMOD_DIST,zdis=TAXODIST,method = "spearman",permutations = 10)
## Table of results
MantelResults=rbind(MantelResults, data.frame(row.names = "dataset21_partitionA",
  TAXON_HM_r=Manteltests$dataset21$partitionA$TAXON_HM$statistic, 
  TAXON_HM_p=Manteltests$dataset21$partitionA$TAXON_HM$signif, 
  DISTR_HM_r=Manteltests$dataset21$partitionA$DISTR_HM$statistic,
  DISTR_HM_p=Manteltests$dataset21$partitionA$DISTR_HM$signif, 
  nspecies=length(species),
  permutations=Manteltests$dataset21$partitionA$TAXON_HM$permutations))
# # # # # # # #
## partition B ##
Manteltests$dataset21$partitionB=list()
species=rownames(NETa)[rownames(NETa)%in%rownames(NETb)]
# GEO
GEOdata= read_excel("spatial_taxonomic_data/CORDOBA2_plantsxsite.xlsx")
GEOdata=as.matrix(GEOdata)
GEOdata[is.na(GEOdata)]=0
rownames(GEOdata)=GEOdata[,1]
GEOdata=GEOdata[,-1]
mode(GEOdata)<-"numeric"
GEOdata=GEOdata[match(species,rownames(GEOdata)),]
# TAXONOMY
TAXOdata= read_excel("spatial_taxonomic_data/CORDOBA2_plants_taxonomy.xlsx")
TAXOdata=as.matrix(TAXOdata)
TAXOdata=TAXOdata[match(species,TAXOdata[,1]),]
rownames(TAXOdata)=TAXOdata[,1]
colnames(TAXOdata)[1]="Species"
# Hypermodules (0= same in both interactions, .5= different in one interaction, 1= different in both interactions)
MODa=paste(MODNETa$Row_labels[match(species,rownames(NETa))],"a", sep="")
MODb=paste(MODNETb$Row_labels[match(species,rownames(NETb))],"b", sep="")
HYPERMOD=data.frame(row.names = species)
for (M in 1:nrow(CONGRUENCE$M1M2)){
  HYPERMOD$NETa[MODa==rownames(CONGRUENCE$M1M2)[M]]=unique(CONGRUENCE$M1M2[M,CONGRUENCE$M1M2[M,]!=0])
}
for (M in 1:ncol(CONGRUENCE$M1M2)){
  HYPERMOD$NETb[MODb==colnames(CONGRUENCE$M1M2)[M]]=unique(CONGRUENCE$M1M2[CONGRUENCE$M1M2[,M]!=0,M])
}
# Distances
TAXODIST=tax2dist(TAXOdata,tax_distance = c(Species=0, Genus=1, Family=2, Order=3, Kingdom=4))
TAXODIST=TAXODIST@distance
GEODIST=vegdist(GEOdata,method = "jaccard",binary = T)
HYPERMOD_DIST=matrix(1, nrow=length(species), ncol=length(species), dimnames = list(species,species))
for (SP1 in 1:nrow(HYPERMOD)){
  for(SP2 in 1: nrow(HYPERMOD)){
    if(HYPERMOD$NETa[SP1]==HYPERMOD$NETa[SP2]){
      HYPERMOD_DIST[SP1,SP2]=HYPERMOD_DIST[SP1,SP2]-.5}
    if(HYPERMOD$NETb[SP1]==HYPERMOD$NETb[SP2]){
      HYPERMOD_DIST[SP1,SP2]=HYPERMOD_DIST[SP1,SP2]-.5}
  }
}
# Mantel Tests
Manteltests$dataset21$partitionB$TAXON_HM=
  mantel.partial(xdis=TAXODIST,ydis=HYPERMOD_DIST,zdis=GEODIST,method = "spearman",permutations = 10)
Manteltests$dataset21$partitionB$DISTR_HM= mantel.partial(xdis=GEODIST,ydis=HYPERMOD_DIST,zdis=TAXODIST,method = "spearman",permutations = 10)
## Table of results
MantelResults=rbind(MantelResults, data.frame(row.names = "dataset21_partitionB",
  TAXON_HM_r=Manteltests$dataset21$partitionB$TAXON_HM$statistic, 
  TAXON_HM_p=Manteltests$dataset21$partitionB$TAXON_HM$signif, 
  DISTR_HM_r=Manteltests$dataset21$partitionB$DISTR_HM$statistic,
  DISTR_HM_p=Manteltests$dataset21$partitionB$DISTR_HM$signif, 
  nspecies=length(species),
  permutations=Manteltests$dataset21$partitionB$TAXON_HM$permutations))
# # # # # # #
## partition C ##
Manteltests$dataset21$partitionC=list()
species=colnames(NETb)[paste(MODNETb$Col_labels,"b",sep="")%in%colnames(CONGRUENCE$M1M2)]
# GEO
GEOdata= read_excel("spatial_taxonomic_data/CORDOBA2_herbivoresxsite.xlsx")
GEOdata=as.matrix(GEOdata)
GEOdata[is.na(GEOdata)]=0
rownames(GEOdata)=GEOdata[,1]
GEOdata=GEOdata[,-1]
mode(GEOdata)<-"numeric"
GEOdata=GEOdata[match(species,rownames(GEOdata)),]
# TAXONOMY
TAXOdata= read_excel("spatial_taxonomic_data/CORDOBA2_herbivores_taxonomy.xlsx")
TAXOdata=as.matrix(TAXOdata)
TAXOdata=TAXOdata[match(species,TAXOdata[,1]),]
rownames(TAXOdata)=TAXOdata[,1]
# missing data
# Missing genus -> consider each as a different genus
TAXOdata[(is.na(TAXOdata[,2])|TAXOdata[,2]=="NA"),2]=paste("NA",1:sum(is.na(TAXOdata[,2])|TAXOdata[,2]=="NA"),sep="")
colnames(TAXOdata)=c("Species","Genus","Family","Superfamily","Suborder","Order")
# Hypermodules (0= same, 1= different)
MODb=paste(MODNETb$Col_labels[match(species,colnames(NETb))],"b", sep="")
HYPERMOD=data.frame(row.names = species)
for (M in 1:ncol(CONGRUENCE$M1M2)){
  HYPERMOD$NETb[MODb==colnames(CONGRUENCE$M1M2)[M]]=unique(CONGRUENCE$M1M2[CONGRUENCE$M1M2[,M]!=0,M])
}
# Distances
TAXODIST=tax2dist(TAXOdata,tax_distance = c(Species=0, Genus=1, Family=2, Superfamily=3, Suborder=4, Order=5))
TAXODIST=TAXODIST@distance
GEODIST=vegdist(GEOdata,method = "jaccard",binary = T)
HYPERMOD_DIST=matrix(NA, nrow=length(species), ncol=length(species), dimnames = list(species,species))
for (SP1 in 1:nrow(HYPERMOD)){
  for(SP2 in 1: nrow(HYPERMOD)){
    if(HYPERMOD$NETb[SP1]==HYPERMOD$NETb[SP2]){
      HYPERMOD_DIST[SP1,SP2]=0
    }else{HYPERMOD_DIST[SP1,SP2]=1}
  }
}
# Mantel Tests
Manteltests$dataset21$partitionC$TAXON_HM=
  mantel.partial(xdis=TAXODIST,ydis=HYPERMOD_DIST,zdis=GEODIST,method = "spearman",permutations = 10)
Manteltests$dataset21$partitionC$DISTR_HM= mantel.partial(xdis=GEODIST,ydis=HYPERMOD_DIST,zdis=TAXODIST,method = "spearman",permutations = 10)
## Table of results
MantelResults=rbind(MantelResults, data.frame(row.names = "dataset21_partitionC",
  TAXON_HM_r=Manteltests$dataset21$partitionC$TAXON_HM$statistic, 
  TAXON_HM_p=Manteltests$dataset21$partitionC$TAXON_HM$signif, 
  DISTR_HM_r=Manteltests$dataset21$partitionC$DISTR_HM$statistic,
  DISTR_HM_p=Manteltests$dataset21$partitionC$DISTR_HM$signif, 
  nspecies=length(species),
  permutations=Manteltests$dataset21$partitionC$TAXON_HM$permutations))
save(MantelResults,Manteltests, file="files_results/Mantel_results.RData")
# # # # # # #

### DATASET 22 ####
datasetID=22
Manteltests$dataset22=list()
load(paste("files_results/dataset_",datasetID,".RData",sep=""))
# # # # # # #
## partition A ##
Manteltests$dataset22$partitionA=list()
species=colnames(NETa)[paste(MODNETa$Col_labels,"a",sep="")%in%rownames(CONGRUENCE$M1M2)]
species[species== "Celtis_ehrenbergiana_"]= "Celtis_ehrenbergiana"
species=unique(species)
# GEO
GEOdata= read_excel("spatial_taxonomic_data/CORDOBA2_plantsxsite.xlsx")
GEOdata=as.matrix(GEOdata)
GEOdata[is.na(GEOdata)]=0
rownames(GEOdata)=GEOdata[,1]
GEOdata=GEOdata[,-1]
mode(GEOdata)<-"numeric"
GEOdata=GEOdata[match(species,rownames(GEOdata)),]
# TAXONOMY
TAXOdata= read_excel("spatial_taxonomic_data/CORDOBA2_plants_taxonomy.xlsx")
TAXOdata=as.matrix(TAXOdata)
TAXOdata=TAXOdata[match(species,TAXOdata[,1]),]
rownames(TAXOdata)=TAXOdata[,1]
colnames(TAXOdata)[1]="Species"
# missing data
# Missing genus -> consider each as a different genus
TAXOdata[(is.na(TAXOdata[,2])|TAXOdata[,2]=="NA"),2]=paste("NA",1:sum(is.na(TAXOdata[,2])|TAXOdata[,2]=="NA"),sep="")
# Missing family -> remove species
MISSING_FAMILY=is.na(TAXOdata[,3])|TAXOdata[,3]=="NA"
species=species[!MISSING_FAMILY]
TAXOdata=TAXOdata[!MISSING_FAMILY,]
GEOdata=GEOdata[!MISSING_FAMILY,]
# Hypermodules (0= same, 1= different)
MODa=paste(MODNETa$Col_labels[match(species,colnames(NETa))],"a", sep="")
HYPERMOD=data.frame(row.names = species)
for (M in 1:nrow(CONGRUENCE$M1M2)){
  HYPERMOD$NETa[MODa==rownames(CONGRUENCE$M1M2)[M]]=unique(CONGRUENCE$M1M2[M,CONGRUENCE$M1M2[M,]!=0])
}
# Distances
TAXODIST=tax2dist(TAXOdata,tax_distance = c(Species=0, Genus=1, Family=2, Order=3, Kingdom=4))
TAXODIST=TAXODIST@distance
GEODIST=vegdist(GEOdata,method = "jaccard",binary = T)
HYPERMOD_DIST=matrix(NA, nrow=length(species), ncol=length(species), dimnames = list(species,species))
for (SP1 in 1:nrow(HYPERMOD)){
  for(SP2 in 1: nrow(HYPERMOD)){
    if(HYPERMOD$NETa[SP1]==HYPERMOD$NETa[SP2]){
      HYPERMOD_DIST[SP1,SP2]=0
    }else{HYPERMOD_DIST[SP1,SP2]=1}
  }
}
# Mantel Tests
Manteltests$dataset22$partitionA$TAXON_HM=
  mantel.partial(xdis=TAXODIST,ydis=HYPERMOD_DIST,zdis=GEODIST,method = "spearman",permutations = 10)
Manteltests$dataset22$partitionA$DISTR_HM= mantel.partial(xdis=GEODIST,ydis=HYPERMOD_DIST,zdis=TAXODIST,method = "spearman",permutations = 10)
## Table of results
MantelResults=rbind(MantelResults, data.frame(row.names = "dataset22_partitionA",
  TAXON_HM_r=Manteltests$dataset22$partitionA$TAXON_HM$statistic, 
  TAXON_HM_p=Manteltests$dataset22$partitionA$TAXON_HM$signif, 
  DISTR_HM_r=Manteltests$dataset22$partitionA$DISTR_HM$statistic,
  DISTR_HM_p=Manteltests$dataset22$partitionA$DISTR_HM$signif, 
  nspecies=length(species),
  permutations=Manteltests$dataset22$partitionA$TAXON_HM$permutations))
# # # # # # # #
## partition B ##
Manteltests$dataset22$partitionB=list()
species=rownames(NETa)[rownames(NETa)%in%rownames(NETb)]
# GEO
GEOdata= read_excel("spatial_taxonomic_data/CORDOBA2_herbivoresxsite.xlsx")
GEOdata=as.matrix(GEOdata)
GEOdata[is.na(GEOdata)]=0
rownames(GEOdata)=GEOdata[,1]
GEOdata=GEOdata[,-1]
mode(GEOdata)<-"numeric"
GEOdata=GEOdata[match(species,rownames(GEOdata)),]
# TAXONOMY
TAXOdata= read_excel("spatial_taxonomic_data/CORDOBA2_herbivores_taxonomy.xlsx")
TAXOdata=as.matrix(TAXOdata)
TAXOdata=TAXOdata[match(species,TAXOdata[,1]),]
rownames(TAXOdata)=TAXOdata[,1]
# missing data
# Missing genus -> consider each as a different genus
TAXOdata[(is.na(TAXOdata[,2])|TAXOdata[,2]=="NA"),2]=paste("NA",1:sum(is.na(TAXOdata[,2])|TAXOdata[,2]=="NA"),sep="")
colnames(TAXOdata)=c("Species","Genus","Family","Superfamily","Suborder","Order")
# Hypermodules (0= same in both interactions, .5= different in one interaction, 1= different in both interactions)
MODa=paste(MODNETa$Row_labels[match(species,rownames(NETa))],"a", sep="")
MODb=paste(MODNETb$Row_labels[match(species,rownames(NETb))],"b", sep="")
HYPERMOD=data.frame(row.names = species)
for (M in 1:nrow(CONGRUENCE$M1M2)){
  HYPERMOD$NETa[MODa==rownames(CONGRUENCE$M1M2)[M]]=unique(CONGRUENCE$M1M2[M,CONGRUENCE$M1M2[M,]!=0])
}
for (M in 1:ncol(CONGRUENCE$M1M2)){
  HYPERMOD$NETb[MODb==colnames(CONGRUENCE$M1M2)[M]]=unique(CONGRUENCE$M1M2[CONGRUENCE$M1M2[,M]!=0,M])
}
# Distances
TAXODIST=tax2dist(TAXOdata,tax_distance = c(Species=0, Genus=1, Family=2, Superfamily=3, Suborder=4, Order=5))
TAXODIST=TAXODIST@distance
GEODIST=vegdist(GEOdata,method = "jaccard",binary = T)
HYPERMOD_DIST=matrix(1, nrow=length(species), ncol=length(species), dimnames = list(species,species))
for (SP1 in 1:nrow(HYPERMOD)){
  for(SP2 in 1: nrow(HYPERMOD)){
    if(HYPERMOD$NETa[SP1]==HYPERMOD$NETa[SP2]){
      HYPERMOD_DIST[SP1,SP2]=HYPERMOD_DIST[SP1,SP2]-.5}
    if(HYPERMOD$NETb[SP1]==HYPERMOD$NETb[SP2]){
      HYPERMOD_DIST[SP1,SP2]=HYPERMOD_DIST[SP1,SP2]-.5}
  }
}
# Mantel Tests
Manteltests$dataset22$partitionB$TAXON_HM=
  mantel.partial(xdis=TAXODIST,ydis=HYPERMOD_DIST,zdis=GEODIST,method = "spearman",permutations = 10)
Manteltests$dataset22$partitionB$DISTR_HM= mantel.partial(xdis=GEODIST,ydis=HYPERMOD_DIST,zdis=TAXODIST,method = "spearman",permutations = 10)
## Table of results
MantelResults=rbind(MantelResults, data.frame(row.names = "dataset22_partitionB",
  TAXON_HM_r=Manteltests$dataset22$partitionB$TAXON_HM$statistic, 
  TAXON_HM_p=Manteltests$dataset22$partitionB$TAXON_HM$signif, 
  DISTR_HM_r=Manteltests$dataset22$partitionB$DISTR_HM$statistic,
  DISTR_HM_p=Manteltests$dataset22$partitionB$DISTR_HM$signif, 
  nspecies=length(species),
  permutations=Manteltests$dataset22$partitionB$TAXON_HM$permutations))
# # # # # # #
## partition C ##
Manteltests$dataset22$partitionC=list()
species=colnames(NETb)[paste(MODNETb$Col_labels,"b",sep="")%in%colnames(CONGRUENCE$M1M2)]
# GEO
GEOdata= read_excel("spatial_taxonomic_data/CORDOBA2_predatorsxsite.xlsx")
GEOdata=as.matrix(GEOdata)
GEOdata[is.na(GEOdata)]=0
rownames(GEOdata)=GEOdata[,1]
GEOdata=GEOdata[,-1]
mode(GEOdata)<-"numeric"
GEOdata=GEOdata[match(species,rownames(GEOdata)),]
# TAXONOMY
TAXOdata= read_excel("spatial_taxonomic_data/CORDOBA2_predators_taxonomy.xlsx")
TAXOdata=as.matrix(TAXOdata)
TAXOdata=TAXOdata[match(species,TAXOdata[,1]),]
rownames(TAXOdata)=TAXOdata[,1]
# missing data
# Missing genus -> consider each as a different genus
TAXOdata[(is.na(TAXOdata[,2])|TAXOdata[,2]=="NA"|TAXOdata[,2]==""),2]=paste("NA",1:sum(is.na(TAXOdata[,2])|TAXOdata[,2]=="NA"|TAXOdata[,2]==""),sep="")
# Missing family -> remove species
MISSING_FAMILY=is.na(TAXOdata[,3])|TAXOdata[,3]=="NA"|TAXOdata[,3]==""
species=species[!MISSING_FAMILY]
TAXOdata=TAXOdata[!MISSING_FAMILY,]
GEOdata=GEOdata[!MISSING_FAMILY,]
# Hypermodules (0= same, 1= different)
MODb=paste(MODNETb$Col_labels[match(species,colnames(NETb))],"b", sep="")
HYPERMOD=data.frame(row.names = species)
for (M in 1:ncol(CONGRUENCE$M1M2)){
  HYPERMOD$NETb[MODb==colnames(CONGRUENCE$M1M2)[M]]=unique(CONGRUENCE$M1M2[CONGRUENCE$M1M2[,M]!=0,M])
}
# Distances
TAXODIST=tax2dist(TAXOdata,tax_distance = c(Species=0, Genus=1, Family=2, Order=3, Kingdom=4))
TAXODIST=TAXODIST@distance
GEODIST=vegdist(GEOdata,method = "jaccard",binary = T)
HYPERMOD_DIST=matrix(NA, nrow=length(species), ncol=length(species), dimnames = list(species,species))
for (SP1 in 1:nrow(HYPERMOD)){
  for(SP2 in 1: nrow(HYPERMOD)){
    if(HYPERMOD$NETb[SP1]==HYPERMOD$NETb[SP2]){
      HYPERMOD_DIST[SP1,SP2]=0
    }else{HYPERMOD_DIST[SP1,SP2]=1}
  }
}
# Mantel Tests
Manteltests$dataset22$partitionC$TAXON_HM=
  mantel.partial(xdis=TAXODIST,ydis=HYPERMOD_DIST,zdis=GEODIST,method = "spearman",permutations = 10)
Manteltests$dataset22$partitionC$DISTR_HM= mantel.partial(xdis=GEODIST,ydis=HYPERMOD_DIST,zdis=TAXODIST,method = "spearman",permutations = 10)
## Table of results
MantelResults=rbind(MantelResults, data.frame(row.names = "dataset22_partitionC",
  TAXON_HM_r=Manteltests$dataset22$partitionC$TAXON_HM$statistic, 
  TAXON_HM_p=Manteltests$dataset22$partitionC$TAXON_HM$signif, 
  DISTR_HM_r=Manteltests$dataset22$partitionC$DISTR_HM$statistic,
  DISTR_HM_p=Manteltests$dataset22$partitionC$DISTR_HM$signif, 
  nspecies=length(species),
  permutations=Manteltests$dataset22$partitionC$TAXON_HM$permutations))
save(MantelResults,Manteltests, file="files_results/Mantel_results.RData")
# # # # # # #

### DATASET 23 ####
datasetID=23
Manteltests$dataset23=list()
load(paste("files_results/dataset_",datasetID,".RData",sep=""))
# # # # # # #
## partition A ##
Manteltests$dataset23$partitionA=list()
species=colnames(NETa)[paste(MODNETa$Col_labels,"a",sep="")%in%rownames(CONGRUENCE$M1M2)]
species[species== "Celtis_ehrenbergiana_"]= "Celtis_ehrenbergiana"
species=unique(species)
# GEO
GEOdata= read_excel("spatial_taxonomic_data/CORDOBA2_plantsxsite.xlsx")
GEOdata=as.matrix(GEOdata)
GEOdata[is.na(GEOdata)]=0
rownames(GEOdata)=GEOdata[,1]
GEOdata=GEOdata[,-1]
mode(GEOdata)<-"numeric"
GEOdata=GEOdata[match(species,rownames(GEOdata)),]
# TAXONOMY
TAXOdata= read_excel("spatial_taxonomic_data/CORDOBA2_plants_taxonomy.xlsx")
TAXOdata=as.matrix(TAXOdata)
TAXOdata=TAXOdata[match(species,TAXOdata[,1]),]
rownames(TAXOdata)=TAXOdata[,1]
colnames(TAXOdata)[1]="Species"
# missing data
# Missing genus -> consider each as a different genus
TAXOdata[(is.na(TAXOdata[,2])|TAXOdata[,2]=="NA"),2]=paste("NA",1:sum(is.na(TAXOdata[,2])|TAXOdata[,2]=="NA"),sep="")
# Missing family -> remove species
MISSING_FAMILY=is.na(TAXOdata[,3])|TAXOdata[,3]=="NA"
species=species[!MISSING_FAMILY]
TAXOdata=TAXOdata[!MISSING_FAMILY,]
GEOdata=GEOdata[!MISSING_FAMILY,]
# Hypermodules (0= same, 1= different)
MODa=paste(MODNETa$Col_labels[match(species,colnames(NETa))],"a", sep="")
HYPERMOD=data.frame(row.names = species)
for (M in 1:nrow(CONGRUENCE$M1M2)){
  HYPERMOD$NETa[MODa==rownames(CONGRUENCE$M1M2)[M]]=unique(CONGRUENCE$M1M2[M,CONGRUENCE$M1M2[M,]!=0])
}
# Distances
TAXODIST=tax2dist(TAXOdata,tax_distance = c(Species=0, Genus=1, Family=2, Order=3, Kingdom=4))
TAXODIST=TAXODIST@distance
GEODIST=vegdist(GEOdata,method = "jaccard",binary = T)
HYPERMOD_DIST=matrix(NA, nrow=length(species), ncol=length(species), dimnames = list(species,species))
for (SP1 in 1:nrow(HYPERMOD)){
  for(SP2 in 1: nrow(HYPERMOD)){
    if(HYPERMOD$NETa[SP1]==HYPERMOD$NETa[SP2]){
      HYPERMOD_DIST[SP1,SP2]=0
    }else{HYPERMOD_DIST[SP1,SP2]=1}
  }
}
# Mantel Tests
Manteltests$dataset23$partitionA$TAXON_HM=
  mantel.partial(xdis=TAXODIST,ydis=HYPERMOD_DIST,zdis=GEODIST,method = "spearman",permutations = 10)
Manteltests$dataset23$partitionA$DISTR_HM= mantel.partial(xdis=GEODIST,ydis=HYPERMOD_DIST,zdis=TAXODIST,method = "spearman",permutations = 10)
## Table of results
MantelResults=rbind(MantelResults, data.frame(row.names = "dataset23_partitionA",
  TAXON_HM_r=Manteltests$dataset23$partitionA$TAXON_HM$statistic, 
  TAXON_HM_p=Manteltests$dataset23$partitionA$TAXON_HM$signif, 
  DISTR_HM_r=Manteltests$dataset23$partitionA$DISTR_HM$statistic,
  DISTR_HM_p=Manteltests$dataset23$partitionA$DISTR_HM$signif, 
  nspecies=length(species),
  permutations=Manteltests$dataset23$partitionA$TAXON_HM$permutations))
# # # # # # # #
## partition B ##
Manteltests$dataset23$partitionB=list()
species=rownames(NETa)[rownames(NETa)%in%rownames(NETb)]
# GEO
GEOdata= read_excel("spatial_taxonomic_data/CORDOBA2_herbivoresxsite.xlsx")
GEOdata=as.matrix(GEOdata)
GEOdata[is.na(GEOdata)]=0
rownames(GEOdata)=GEOdata[,1]
GEOdata=GEOdata[,-1]
mode(GEOdata)<-"numeric"
GEOdata=GEOdata[match(species,rownames(GEOdata)),]
# TAXONOMY
TAXOdata= read_excel("spatial_taxonomic_data/CORDOBA2_herbivores_taxonomy.xlsx")
TAXOdata=as.matrix(TAXOdata)
TAXOdata=TAXOdata[match(species,TAXOdata[,1]),]
rownames(TAXOdata)=TAXOdata[,1]
# missing data
# Missing genus -> consider each as a different genus
TAXOdata[(is.na(TAXOdata[,2])|TAXOdata[,2]=="NA"),2]=paste("NA",1:sum(is.na(TAXOdata[,2])|TAXOdata[,2]=="NA"),sep="")
colnames(TAXOdata)=c("Species","Genus","Family","Superfamily","Suborder","Order")
# Hypermodules (0= same in both interactions, .5= different in one interaction, 1= different in both interactions)
MODa=paste(MODNETa$Row_labels[match(species,rownames(NETa))],"a", sep="")
MODb=paste(MODNETb$Row_labels[match(species,rownames(NETb))],"b", sep="")
HYPERMOD=data.frame(row.names = species)
for (M in 1:nrow(CONGRUENCE$M1M2)){
  HYPERMOD$NETa[MODa==rownames(CONGRUENCE$M1M2)[M]]=unique(CONGRUENCE$M1M2[M,CONGRUENCE$M1M2[M,]!=0])
}
for (M in 1:ncol(CONGRUENCE$M1M2)){
  HYPERMOD$NETb[MODb==colnames(CONGRUENCE$M1M2)[M]]=unique(CONGRUENCE$M1M2[CONGRUENCE$M1M2[,M]!=0,M])
}
# Distances
TAXODIST=tax2dist(TAXOdata,tax_distance = c(Species=0, Genus=1, Family=2, Superfamily=3, Suborder=4, Order=5))
TAXODIST=TAXODIST@distance
GEODIST=vegdist(GEOdata,method = "jaccard",binary = T)
HYPERMOD_DIST=matrix(1, nrow=length(species), ncol=length(species), dimnames = list(species,species))
for (SP1 in 1:nrow(HYPERMOD)){
  for(SP2 in 1: nrow(HYPERMOD)){
    if(HYPERMOD$NETa[SP1]==HYPERMOD$NETa[SP2]){
      HYPERMOD_DIST[SP1,SP2]=HYPERMOD_DIST[SP1,SP2]-.5}
    if(HYPERMOD$NETb[SP1]==HYPERMOD$NETb[SP2]){
      HYPERMOD_DIST[SP1,SP2]=HYPERMOD_DIST[SP1,SP2]-.5}
  }
}
# Mantel Tests
Manteltests$dataset23$partitionB$TAXON_HM=
  mantel.partial(xdis=TAXODIST,ydis=HYPERMOD_DIST,zdis=GEODIST,method = "spearman",permutations = 10)
Manteltests$dataset23$partitionB$DISTR_HM= mantel.partial(xdis=GEODIST,ydis=HYPERMOD_DIST,zdis=TAXODIST,method = "spearman",permutations = 10)
## Table of results
MantelResults=rbind(MantelResults, data.frame(row.names = "dataset23_partitionB",
  TAXON_HM_r=Manteltests$dataset23$partitionB$TAXON_HM$statistic, 
  TAXON_HM_p=Manteltests$dataset23$partitionB$TAXON_HM$signif, 
  DISTR_HM_r=Manteltests$dataset23$partitionB$DISTR_HM$statistic,
  DISTR_HM_p=Manteltests$dataset23$partitionB$DISTR_HM$signif, 
  nspecies=length(species),
  permutations=Manteltests$dataset23$partitionB$TAXON_HM$permutations))
# # # # # # #
## partition C ##
Manteltests$dataset23$partitionC=list()
species=colnames(NETb)[paste(MODNETb$Col_labels,"b",sep="")%in%colnames(CONGRUENCE$M1M2)]
species=species[species!="Myrmelachista_sp"]
# GEO
GEOdata= read_excel("spatial_taxonomic_data/CORDOBA2_plant.defenderxsite.xlsx")
GEOdata=as.matrix(GEOdata)
GEOdata[is.na(GEOdata)]=0
rownames(GEOdata)=GEOdata[,1]
GEOdata=GEOdata[,-1]
mode(GEOdata)<-"numeric"
GEOdata=GEOdata[match(species,rownames(GEOdata)),]
# TAXONOMY
TAXOdata= read_excel("spatial_taxonomic_data/CORDOBA2_plant.defender_taxonomy.xlsx")
TAXOdata=as.matrix(TAXOdata)
TAXOdata=TAXOdata[match(species,TAXOdata[,1]),]
rownames(TAXOdata)=TAXOdata[,1]
# missing data
# Missing genus -> consider each as a different genus
TAXOdata[(is.na(TAXOdata[,2])|TAXOdata[,2]=="NA"|TAXOdata[,2]==""),2]=paste("NA",1:sum(is.na(TAXOdata[,2])|TAXOdata[,2]=="NA"|TAXOdata[,2]==""),sep="")
# Missing family -> remove species
MISSING_FAMILY=is.na(TAXOdata[,3])|TAXOdata[,3]=="NA"|TAXOdata[,3]==""
species=species[!MISSING_FAMILY]
TAXOdata=TAXOdata[!MISSING_FAMILY,]
GEOdata=GEOdata[!MISSING_FAMILY,]
# Hypermodules (0= same, 1= different)
MODb=paste(MODNETb$Col_labels[match(species,colnames(NETb))],"b", sep="")
HYPERMOD=data.frame(row.names = species)
for (M in 1:ncol(CONGRUENCE$M1M2)){
  HYPERMOD$NETb[MODb==colnames(CONGRUENCE$M1M2)[M]]=unique(CONGRUENCE$M1M2[CONGRUENCE$M1M2[,M]!=0,M])
}
# Distances
TAXODIST=tax2dist(TAXOdata,tax_distance = c(Species=0, Genus=1, Family=2, Order=3, Kingdom=4))
TAXODIST=TAXODIST@distance
GEODIST=vegdist(GEOdata,method = "jaccard",binary = T)
HYPERMOD_DIST=matrix(NA, nrow=length(species), ncol=length(species), dimnames = list(species,species))
for (SP1 in 1:nrow(HYPERMOD)){
  for(SP2 in 1: nrow(HYPERMOD)){
    if(HYPERMOD$NETb[SP1]==HYPERMOD$NETb[SP2]){
      HYPERMOD_DIST[SP1,SP2]=0
    }else{HYPERMOD_DIST[SP1,SP2]=1}
  }
}
# Mantel Tests
Manteltests$dataset23$partitionC$TAXON_HM=
  mantel.partial(xdis=TAXODIST,ydis=HYPERMOD_DIST,zdis=GEODIST,method = "spearman",permutations = 10)
Manteltests$dataset23$partitionC$DISTR_HM= mantel.partial(xdis=GEODIST,ydis=HYPERMOD_DIST,zdis=TAXODIST,method = "spearman",permutations = 10)
## Table of results
MantelResults=rbind(MantelResults, data.frame(row.names = "dataset23_partitionC",
  TAXON_HM_r=Manteltests$dataset23$partitionC$TAXON_HM$statistic, 
  TAXON_HM_p=Manteltests$dataset23$partitionC$TAXON_HM$signif, 
  DISTR_HM_r=Manteltests$dataset23$partitionC$DISTR_HM$statistic,
  DISTR_HM_p=Manteltests$dataset23$partitionC$DISTR_HM$signif, 
  nspecies=length(species),
  permutations=Manteltests$dataset23$partitionC$TAXON_HM$permutations))
save(MantelResults,Manteltests, file="files_results/Mantel_results.RData")
# # # # # # #
### DATASET 24 ####
datasetID=24
Manteltests$dataset24=list()
load(paste("files_results/dataset_",datasetID,".RData",sep=""))
# # # # # # # #
## partition A ##
Manteltests$dataset24$partitionA=list()
species=colnames(NETa)[paste(MODNETa$Col_labels,"a",sep="")%in%rownames(CONGRUENCE$M1M2)]
# GEO
GEOdata= read_excel("spatial_taxonomic_data/CORDOBA2_herbivore.defenderxsite.xlsx")
GEOdata=as.matrix(GEOdata)
GEOdata[is.na(GEOdata)]=0
rownames(GEOdata)=GEOdata[,1]
GEOdata=GEOdata[,-1]
mode(GEOdata)<-"numeric"
GEOdata=GEOdata[match(species,rownames(GEOdata)),]
# TAXONOMY
TAXOdata= read_excel("spatial_taxonomic_data/CORDOBA2_herbivore.defender_taxonomy.xlsx")
TAXOdata=as.matrix(TAXOdata)
TAXOdata=TAXOdata[match(species,TAXOdata[,1]),]
rownames(TAXOdata)=TAXOdata[,1]
# missing data
# Missing genus -> consider each as a different genus
TAXOdata[(is.na(TAXOdata[,2])|TAXOdata[,2]=="NA"),2]=paste("NA",1:sum(is.na(TAXOdata[,2])|TAXOdata[,2]=="NA"),sep="")
# Missing family -> remove species
MISSING_FAMILY=is.na(TAXOdata[,3])|TAXOdata[,3]=="NA"
species=species[!MISSING_FAMILY]
TAXOdata=TAXOdata[!MISSING_FAMILY,]
GEOdata=GEOdata[!MISSING_FAMILY,]
# Hypermodules (0= same, 1= different)
MODa=paste(MODNETa$Col_labels[match(species,colnames(NETa))],"a", sep="")
HYPERMOD=data.frame(row.names = species)
for (M in 1:nrow(CONGRUENCE$M1M2)){
  HYPERMOD$NETa[MODa==rownames(CONGRUENCE$M1M2)[M]]=unique(CONGRUENCE$M1M2[M,CONGRUENCE$M1M2[M,]!=0])
}
# Distances
TAXODIST=tax2dist(TAXOdata,tax_distance = c(Species=0, Genus=1, Family=2, Order=3, Kingdom=4))
TAXODIST=TAXODIST@distance
GEODIST=vegdist(GEOdata,method = "jaccard",binary = T)
HYPERMOD_DIST=matrix(NA, nrow=length(species), ncol=length(species), dimnames = list(species,species))
for (SP1 in 1:nrow(HYPERMOD)){
  for(SP2 in 1: nrow(HYPERMOD)){
    if(HYPERMOD$NETa[SP1]==HYPERMOD$NETa[SP2]){
      HYPERMOD_DIST[SP1,SP2]=0
    }else{HYPERMOD_DIST[SP1,SP2]=1}
  }
}
# Mantel Tests
Manteltests$dataset24$partitionA$TAXON_HM=
  mantel.partial(xdis=TAXODIST,ydis=HYPERMOD_DIST,zdis=GEODIST,method = "spearman",permutations = 10)
Manteltests$dataset24$partitionA$DISTR_HM= mantel.partial(xdis=GEODIST,ydis=HYPERMOD_DIST,zdis=TAXODIST,method = "spearman",permutations = 10)
## Table of results
MantelResults=rbind(MantelResults, data.frame(row.names = "dataset24_partitionA",
  TAXON_HM_r=Manteltests$dataset24$partitionA$TAXON_HM$statistic, 
  TAXON_HM_p=Manteltests$dataset24$partitionA$TAXON_HM$signif, 
  DISTR_HM_r=Manteltests$dataset24$partitionA$DISTR_HM$statistic,
  DISTR_HM_p=Manteltests$dataset24$partitionA$DISTR_HM$signif, 
  nspecies=length(species),
  permutations=Manteltests$dataset24$partitionA$TAXON_HM$permutations))
# # # # # # # #
## partition B ##
Manteltests$dataset24$partitionB=list()
species=rownames(NETa)[rownames(NETa)%in%rownames(NETb)]
# GEO
GEOdata= read_excel("spatial_taxonomic_data/CORDOBA2_herbivoresxsite.xlsx")
GEOdata=as.matrix(GEOdata)
GEOdata[is.na(GEOdata)]=0
rownames(GEOdata)=GEOdata[,1]
GEOdata=GEOdata[,-1]
mode(GEOdata)<-"numeric"
GEOdata=GEOdata[match(species,rownames(GEOdata)),]
# TAXONOMY
TAXOdata= read_excel("spatial_taxonomic_data/CORDOBA2_herbivores_taxonomy.xlsx")
TAXOdata=as.matrix(TAXOdata)
TAXOdata=TAXOdata[match(species,TAXOdata[,1]),]
rownames(TAXOdata)=TAXOdata[,1]
# missing data
# Missing genus -> consider each as a different genus
TAXOdata[(is.na(TAXOdata[,2])|TAXOdata[,2]=="NA"),2]=paste("NA",1:sum(is.na(TAXOdata[,2])|TAXOdata[,2]=="NA"),sep="")
colnames(TAXOdata)=c("Species","Genus","Family","Superfamily","Suborder","Order")
# Hypermodules (0= same in both interactions, .5= different in one interaction, 1= different in both interactions)
MODa=paste(MODNETa$Row_labels[match(species,rownames(NETa))],"a", sep="")
MODb=paste(MODNETb$Row_labels[match(species,rownames(NETb))],"b", sep="")
HYPERMOD=data.frame(row.names = species)
for (M in 1:nrow(CONGRUENCE$M1M2)){
  HYPERMOD$NETa[MODa==rownames(CONGRUENCE$M1M2)[M]]=unique(CONGRUENCE$M1M2[M,CONGRUENCE$M1M2[M,]!=0])
}
for (M in 1:ncol(CONGRUENCE$M1M2)){
  HYPERMOD$NETb[MODb==colnames(CONGRUENCE$M1M2)[M]]=unique(CONGRUENCE$M1M2[CONGRUENCE$M1M2[,M]!=0,M])
}
# Distances
TAXODIST=tax2dist(TAXOdata,tax_distance = c(Species=0, Genus=1, Family=2, Superfamily=3, Suborder=4, Order=5))
TAXODIST=TAXODIST@distance
GEODIST=vegdist(GEOdata,method = "jaccard",binary = T)
HYPERMOD_DIST=matrix(1, nrow=length(species), ncol=length(species), dimnames = list(species,species))
for (SP1 in 1:nrow(HYPERMOD)){
  for(SP2 in 1: nrow(HYPERMOD)){
    if(HYPERMOD$NETa[SP1]==HYPERMOD$NETa[SP2]){
      HYPERMOD_DIST[SP1,SP2]=HYPERMOD_DIST[SP1,SP2]-.5}
    if(HYPERMOD$NETb[SP1]==HYPERMOD$NETb[SP2]){
      HYPERMOD_DIST[SP1,SP2]=HYPERMOD_DIST[SP1,SP2]-.5}
  }
}
# Mantel Tests
Manteltests$dataset24$partitionB$TAXON_HM=
  mantel.partial(xdis=TAXODIST,ydis=HYPERMOD_DIST,zdis=GEODIST,method = "spearman",permutations = 10)
Manteltests$dataset24$partitionB$DISTR_HM= mantel.partial(xdis=GEODIST,ydis=HYPERMOD_DIST,zdis=TAXODIST,method = "spearman",permutations = 10)
## Table of results
MantelResults=rbind(MantelResults, data.frame(row.names = "dataset24_partitionB",
  TAXON_HM_r=Manteltests$dataset24$partitionB$TAXON_HM$statistic, 
  TAXON_HM_p=Manteltests$dataset24$partitionB$TAXON_HM$signif, 
  DISTR_HM_r=Manteltests$dataset24$partitionB$DISTR_HM$statistic,
  DISTR_HM_p=Manteltests$dataset24$partitionB$DISTR_HM$signif, 
  nspecies=length(species),
  permutations=Manteltests$dataset24$partitionB$TAXON_HM$permutations))
# # # # # # #
## partition C ##
Manteltests$dataset24$partitionC=list()
species=colnames(NETb)[paste(MODNETb$Col_labels,"b",sep="")%in%colnames(CONGRUENCE$M1M2)]
# GEO
GEOdata= read_excel("spatial_taxonomic_data/CORDOBA2_enemiesxsite.xlsx")
GEOdata=as.matrix(GEOdata)
GEOdata[is.na(GEOdata)]=0
rownames(GEOdata)=GEOdata[,1]
GEOdata=GEOdata[,-1]
mode(GEOdata)<-"numeric"
GEOdata=GEOdata[match(species,rownames(GEOdata)),]
# TAXONOMY
TAXOdata= read_excel("spatial_taxonomic_data/CORDOBA2_enemies_taxonomy.xlsx")
TAXOdata=as.matrix(TAXOdata)
TAXOdata=TAXOdata[match(species,TAXOdata[,1]),]
rownames(TAXOdata)=TAXOdata[,1]
# missing data
# Missing genus -> consider each as a different genus
TAXOdata[(is.na(TAXOdata[,2])|TAXOdata[,2]=="NA"|TAXOdata[,2]==""),2]=paste("NA",1:sum(is.na(TAXOdata[,2])|TAXOdata[,2]=="NA"|TAXOdata[,2]==""),sep="")
# Missing family -> remove species
MISSING_FAMILY=is.na(TAXOdata[,3])|TAXOdata[,3]=="NA"|TAXOdata[,3]==""
species=species[!MISSING_FAMILY]
TAXOdata=TAXOdata[!MISSING_FAMILY,]
GEOdata=GEOdata[!MISSING_FAMILY,]
# Hypermodules (0= same, 1= different)
MODb=paste(MODNETb$Col_labels[match(species,colnames(NETb))],"b", sep="")
HYPERMOD=data.frame(row.names = species)
for (M in 1:ncol(CONGRUENCE$M1M2)){
  HYPERMOD$NETb[MODb==colnames(CONGRUENCE$M1M2)[M]]=unique(CONGRUENCE$M1M2[CONGRUENCE$M1M2[,M]!=0,M])
}
# Distances
TAXODIST=tax2dist(TAXOdata,tax_distance = c(Species=0, Genus=1, Family=2, Order=3, Kingdom=4))
TAXODIST=TAXODIST@distance
GEODIST=vegdist(GEOdata,method = "jaccard",binary = T)
HYPERMOD_DIST=matrix(NA, nrow=length(species), ncol=length(species), dimnames = list(species,species))
for (SP1 in 1:nrow(HYPERMOD)){
  for(SP2 in 1: nrow(HYPERMOD)){
    if(HYPERMOD$NETb[SP1]==HYPERMOD$NETb[SP2]){
      HYPERMOD_DIST[SP1,SP2]=0
    }else{HYPERMOD_DIST[SP1,SP2]=1}
  }
}
# Mantel Tests
Manteltests$dataset24$partitionC$TAXON_HM=
  mantel.partial(xdis=TAXODIST,ydis=HYPERMOD_DIST,zdis=GEODIST,method = "spearman",permutations = 10)
Manteltests$dataset24$partitionC$DISTR_HM= mantel.partial(xdis=GEODIST,ydis=HYPERMOD_DIST,zdis=TAXODIST,method = "spearman",permutations = 10)
## Table of results
MantelResults=rbind(MantelResults, data.frame(row.names = "dataset24_partitionC",
  TAXON_HM_r=Manteltests$dataset24$partitionC$TAXON_HM$statistic, 
  TAXON_HM_p=Manteltests$dataset24$partitionC$TAXON_HM$signif, 
  DISTR_HM_r=Manteltests$dataset24$partitionC$DISTR_HM$statistic,
  DISTR_HM_p=Manteltests$dataset24$partitionC$DISTR_HM$signif, 
  nspecies=length(species),
  permutations=Manteltests$dataset24$partitionC$TAXON_HM$permutations))
save(MantelResults,Manteltests, file="files_results/Mantel_results.RData")
# # # # # # #

### DATASET 29 ####
datasetID=29
Manteltests$dataset29=list()
load(paste("files_results/dataset_",datasetID,".RData",sep=""))
# # # # # # #
## partition A ##
Manteltests$dataset29$partitionA=list()
species=colnames(NETa)[paste(MODNETa$Col_labels,"a",sep="")%in%rownames(CONGRUENCE$M1M2)]
# GEO
GEOdata= read_excel("spatial_taxonomic_data/HACKETT_HH_herbivoresxsites.xlsx")
GEOdata=as.matrix(GEOdata)
GEOdata[is.na(GEOdata)]=0
rownames(GEOdata)=GEOdata[,1]
GEOdata=GEOdata[,-1]
mode(GEOdata)<-"numeric"
GEOdata=GEOdata[match(species,rownames(GEOdata)),]
# TAXONOMY
TAXOdata= read_excel("spatial_taxonomic_data/HACKETT_HH_herbivores_taxonomy.xlsx")
TAXOdata=as.matrix(TAXOdata)
TAXOdata=TAXOdata[match(species,TAXOdata[,1]),]
rownames(TAXOdata)=TAXOdata[,1]
TAXOdata=TAXOdata[,-5]
TAXOdata=cbind(TAXOdata,data.frame(Class=rep("Insecta",21)))
# Hypermodules (0= same, 1= different)
MODa=paste(MODNETa$Col_labels[match(species,colnames(NETa))],"a", sep="")
HYPERMOD=data.frame(row.names = species)
for (M in 1:nrow(CONGRUENCE$M1M2)){
  HYPERMOD$NETa[MODa==rownames(CONGRUENCE$M1M2)[M]]=unique(CONGRUENCE$M1M2[M,CONGRUENCE$M1M2[M,]!=0])
}
# Distances
TAXODIST=tax2dist(TAXOdata,tax_distance = c(Species=0, Genus=1, Family=2, Order=3, Class=4))
TAXODIST=TAXODIST@distance
GEODIST=vegdist(GEOdata,method = "jaccard",binary = T)
HYPERMOD_DIST=matrix(NA, nrow=length(species), ncol=length(species), dimnames = list(species,species))
for (SP1 in 1:nrow(HYPERMOD)){
  for(SP2 in 1: nrow(HYPERMOD)){
    if(HYPERMOD$NETa[SP1]==HYPERMOD$NETa[SP2]){
      HYPERMOD_DIST[SP1,SP2]=0
    }else{HYPERMOD_DIST[SP1,SP2]=1}
  }
}
# Mantel Tests
Manteltests$dataset29$partitionA$TAXON_HM=
  mantel.partial(xdis=TAXODIST,ydis=HYPERMOD_DIST,zdis=GEODIST,method = "spearman",permutations = 10)
Manteltests$dataset29$partitionA$DISTR_HM= mantel.partial(xdis=GEODIST,ydis=HYPERMOD_DIST,zdis=TAXODIST,method = "spearman",permutations = 10)
## Table of results
MantelResults=rbind(MantelResults, data.frame(row.names = "dataset29_partitionA",
  TAXON_HM_r=Manteltests$dataset29$partitionA$TAXON_HM$statistic, 
  TAXON_HM_p=Manteltests$dataset29$partitionA$TAXON_HM$signif, 
  DISTR_HM_r=Manteltests$dataset29$partitionA$DISTR_HM$statistic,
  DISTR_HM_p=Manteltests$dataset29$partitionA$DISTR_HM$signif, 
  nspecies=length(species),
  permutations=Manteltests$dataset29$partitionA$TAXON_HM$permutations))
# # # # # # # #
## partition B ##
Manteltests$dataset29$partitionB=list()
species=rownames(NETa)[rownames(NETa)%in%rownames(NETb)]
# GEO
GEOdata= read_excel("spatial_taxonomic_data/HACKETT_HH_plantsxsites.xlsx")
GEOdata=as.matrix(GEOdata)
GEOdata[is.na(GEOdata)]=0
rownames(GEOdata)=GEOdata[,1]
GEOdata=GEOdata[,-1]
mode(GEOdata)<-"numeric"
GEOdata=GEOdata[match(species,rownames(GEOdata)),]
# TAXONOMY
TAXOdata= read_excel("spatial_taxonomic_data/HACKETT_HH_plants_taxonomy.xlsx")
TAXOdata=as.matrix(TAXOdata)
TAXOdata=TAXOdata[match(species,TAXOdata[,1]),]
rownames(TAXOdata)=TAXOdata[,1]
# missing data
# Missing genus -> consider each as a different genus
TAXOdata[(is.na(TAXOdata[,2])|TAXOdata[,2]=="NA"),2]=paste("NA",1:sum(is.na(TAXOdata[,2])|TAXOdata[,2]=="NA"),sep="")
# Missing family -> remove species
MISSING_FAMILY=is.na(TAXOdata[,3])|TAXOdata[,3]=="NA"
species=species[!MISSING_FAMILY]
TAXOdata=TAXOdata[!MISSING_FAMILY,]
GEOdata=GEOdata[!MISSING_FAMILY,]
TAXOdata=cbind(TAXOdata, Kingdom=rep("Plantae",7))
# Hypermodules (0= same in both interactions, .5= different in one interaction, 1= different in both interactions)
MODa=paste(MODNETa$Row_labels[match(species,rownames(NETa))],"a", sep="")
MODb=paste(MODNETb$Row_labels[match(species,rownames(NETb))],"b", sep="")
HYPERMOD=data.frame(row.names = species)
for (M in 1:nrow(CONGRUENCE$M1M2)){
  HYPERMOD$NETa[MODa==rownames(CONGRUENCE$M1M2)[M]]=unique(CONGRUENCE$M1M2[M,CONGRUENCE$M1M2[M,]!=0])
}
for (M in 1:ncol(CONGRUENCE$M1M2)){
  HYPERMOD$NETb[MODb==colnames(CONGRUENCE$M1M2)[M]]=unique(CONGRUENCE$M1M2[CONGRUENCE$M1M2[,M]!=0,M])
}
# Distances
TAXODIST=tax2dist(TAXOdata,tax_distance = c(Species=0, Genus=1, Family=2, Order=3, Kingdom=4))
TAXODIST=TAXODIST@distance
GEODIST=vegdist(GEOdata,method = "jaccard",binary = T)
HYPERMOD_DIST=matrix(1, nrow=length(species), ncol=length(species), dimnames = list(species,species))
for (SP1 in 1:nrow(HYPERMOD)){
  for(SP2 in 1: nrow(HYPERMOD)){
    if(HYPERMOD$NETa[SP1]==HYPERMOD$NETa[SP2]){
      HYPERMOD_DIST[SP1,SP2]=HYPERMOD_DIST[SP1,SP2]-.5}
    if(HYPERMOD$NETb[SP1]==HYPERMOD$NETb[SP2]){
      HYPERMOD_DIST[SP1,SP2]=HYPERMOD_DIST[SP1,SP2]-.5}
  }
}
# Mantel Tests
Manteltests$dataset29$partitionB$TAXON_HM=
  mantel.partial(xdis=TAXODIST,ydis=HYPERMOD_DIST,zdis=GEODIST,method = "spearman",permutations = 10)
Manteltests$dataset29$partitionB$DISTR_HM= mantel.partial(xdis=GEODIST,ydis=HYPERMOD_DIST,zdis=TAXODIST,method = "spearman",permutations = 10)
## Table of results
MantelResults=rbind(MantelResults, data.frame(row.names = "dataset29_partitionB",
  TAXON_HM_r=Manteltests$dataset29$partitionB$TAXON_HM$statistic, 
  TAXON_HM_p=Manteltests$dataset29$partitionB$TAXON_HM$signif, 
  DISTR_HM_r=Manteltests$dataset29$partitionB$DISTR_HM$statistic,
  DISTR_HM_p=Manteltests$dataset29$partitionB$DISTR_HM$signif, 
  nspecies=length(species),
  permutations=Manteltests$dataset29$partitionB$TAXON_HM$permutations))
# # # # # # #
## partition C ##
Manteltests$dataset29$partitionC=list()
species=colnames(NETb)[paste(MODNETb$Col_labels,"b",sep="")%in%colnames(CONGRUENCE$M1M2)]
# GEO
GEOdata= read_excel("spatial_taxonomic_data/HACKETT_HH_flowervisitorxsites.xlsx")
GEOdata=as.matrix(GEOdata)
GEOdata[is.na(GEOdata)]=0
rownames(GEOdata)=GEOdata[,1]
GEOdata=GEOdata[,-1]
mode(GEOdata)<-"numeric"
GEOdata=GEOdata[match(species,rownames(GEOdata)),]
# TAXONOMY
TAXOdata= read.csv("spatial_taxonomic_data/HACKETT_HH_flowervisitor taxonomy.csv")
TAXOdata=as.matrix(TAXOdata)
TAXOdata=TAXOdata[match(species,TAXOdata[,1]),]
rownames(TAXOdata)=TAXOdata[,1]
# missing data
# Missing genus -> consider each as a different genus
TAXOdata[(is.na(TAXOdata[,2])|TAXOdata[,2]=="NA"|TAXOdata[,2]==""),2]=paste("NA",1:sum(is.na(TAXOdata[,2])|TAXOdata[,2]=="NA"|TAXOdata[,2]==""),sep="")
# Missing family -> remove species
MISSING_FAMILY=is.na(TAXOdata[,3])|TAXOdata[,3]=="NA"|TAXOdata[,3]==""
species=species[!MISSING_FAMILY]
TAXOdata=TAXOdata[!MISSING_FAMILY,]
GEOdata=GEOdata[!MISSING_FAMILY,]
# Hypermodules (0= same, 1= different)
MODb=paste(MODNETb$Col_labels[match(species,colnames(NETb))],"b", sep="")
HYPERMOD=data.frame(row.names = species)
for (M in 1:ncol(CONGRUENCE$M1M2)){
  HYPERMOD$NETb[MODb==colnames(CONGRUENCE$M1M2)[M]]=unique(CONGRUENCE$M1M2[CONGRUENCE$M1M2[,M]!=0,M])
}
# Distances
TAXODIST=tax2dist(TAXOdata,tax_distance = c(Species=0, Genus=1, Family=2, Order=3, Kingdom=4))
TAXODIST=TAXODIST@distance
GEODIST=vegdist(GEOdata,method = "jaccard",binary = T)
HYPERMOD_DIST=matrix(NA, nrow=length(species), ncol=length(species), dimnames = list(species,species))
for (SP1 in 1:nrow(HYPERMOD)){
  for(SP2 in 1: nrow(HYPERMOD)){
    if(HYPERMOD$NETb[SP1]==HYPERMOD$NETb[SP2]){
      HYPERMOD_DIST[SP1,SP2]=0
    }else{HYPERMOD_DIST[SP1,SP2]=1}
  }
}
# Mantel Tests
Manteltests$dataset29$partitionC$TAXON_HM=
  mantel.partial(xdis=TAXODIST,ydis=HYPERMOD_DIST,zdis=GEODIST,method = "spearman",permutations = 10)
Manteltests$dataset29$partitionC$DISTR_HM= mantel.partial(xdis=GEODIST,ydis=HYPERMOD_DIST,zdis=TAXODIST,method = "spearman",permutations = 10)
## Table of results
MantelResults=rbind(MantelResults, data.frame(row.names = "dataset29_partitionC",
  TAXON_HM_r=Manteltests$dataset29$partitionC$TAXON_HM$statistic, 
  TAXON_HM_p=Manteltests$dataset29$partitionC$TAXON_HM$signif, 
  DISTR_HM_r=Manteltests$dataset29$partitionC$DISTR_HM$statistic,
  DISTR_HM_p=Manteltests$dataset29$partitionC$DISTR_HM$signif, 
  nspecies=length(species),
  permutations=Manteltests$dataset29$partitionC$TAXON_HM$permutations))
save(MantelResults,Manteltests, file="files_results/Mantel_results.RData")
# # # # # # #

### DATASET 30 ####
datasetID=30
Manteltests$dataset30=list()
load(paste("files_results/dataset_",datasetID,".RData",sep=""))
# # # # # # #
## partition A ##
Manteltests$dataset30$partitionA=list()
species=colnames(NETa)[paste(MODNETa$Col_labels,"a",sep="")%in%rownames(CONGRUENCE$M1M2)]
# GEO
GEOdata= read_excel("spatial_taxonomic_data/HACKETT_HH_plantsxsites.xlsx")
GEOdata=as.matrix(GEOdata)
GEOdata[is.na(GEOdata)]=0
rownames(GEOdata)=GEOdata[,1]
GEOdata=GEOdata[,-1]
mode(GEOdata)<-"numeric"
GEOdata=GEOdata[match(species,rownames(GEOdata)),]
# TAXONOMY
TAXOdata= read_excel("spatial_taxonomic_data/HACKETT_HH_plants_taxonomy.xlsx")
TAXOdata=as.matrix(TAXOdata)
TAXOdata=TAXOdata[match(species,TAXOdata[,1]),]
rownames(TAXOdata)=TAXOdata[,1]
TAXOdata=cbind(TAXOdata,data.frame(Kingdom=rep("Plantae",8)))
# Hypermodules (0= same, 1= different)
MODa=paste(MODNETa$Col_labels[match(species,colnames(NETa))],"a", sep="")
HYPERMOD=data.frame(row.names = species)
for (M in 1:nrow(CONGRUENCE$M1M2)){
  HYPERMOD$NETa[MODa==rownames(CONGRUENCE$M1M2)[M]]=unique(CONGRUENCE$M1M2[M,CONGRUENCE$M1M2[M,]!=0])
}
# Distances
TAXODIST=tax2dist(TAXOdata,tax_distance = c(Species=0, Genus=1, Family=2, Order=3, Kingdom=4))
TAXODIST=TAXODIST@distance
GEODIST=vegdist(GEOdata,method = "jaccard",binary = T)
HYPERMOD_DIST=matrix(NA, nrow=length(species), ncol=length(species), dimnames = list(species,species))
for (SP1 in 1:nrow(HYPERMOD)){
  for(SP2 in 1: nrow(HYPERMOD)){
    if(HYPERMOD$NETa[SP1]==HYPERMOD$NETa[SP2]){
      HYPERMOD_DIST[SP1,SP2]=0
    }else{HYPERMOD_DIST[SP1,SP2]=1}
  }
}
# Mantel Tests
Manteltests$dataset30$partitionA$TAXON_HM=
  mantel.partial(xdis=TAXODIST,ydis=HYPERMOD_DIST,zdis=GEODIST,method = "spearman",permutations = 10)
Manteltests$dataset30$partitionA$DISTR_HM= mantel.partial(xdis=GEODIST,ydis=HYPERMOD_DIST,zdis=TAXODIST,method = "spearman",permutations = 10)
## Table of results
MantelResults=rbind(MantelResults, data.frame(row.names = "dataset30_partitionA",
  TAXON_HM_r=Manteltests$dataset30$partitionA$TAXON_HM$statistic, 
  TAXON_HM_p=Manteltests$dataset30$partitionA$TAXON_HM$signif, 
  DISTR_HM_r=Manteltests$dataset30$partitionA$DISTR_HM$statistic,
  DISTR_HM_p=Manteltests$dataset30$partitionA$DISTR_HM$signif, 
  nspecies=length(species),
  permutations=Manteltests$dataset30$partitionA$TAXON_HM$permutations))
# # # # # # # #
## partition B ##
Manteltests$dataset30$partitionB=list()
species=rownames(NETa)[rownames(NETa)%in%rownames(NETb)]
species=species[1:12]# species unknown_lm1 and unknown_lm2 without taxonomy
# GEO
GEOdata= read_excel("spatial_taxonomic_data/HACKETT_HH_herbivoresxsites.xlsx")
GEOdata=as.matrix(GEOdata)
GEOdata[is.na(GEOdata)]=0
rownames(GEOdata)=GEOdata[,1]
GEOdata=GEOdata[,-1]
mode(GEOdata)<-"numeric"
GEOdata=GEOdata[match(species,rownames(GEOdata)),]
# TAXONOMY
TAXOdata= read_excel("spatial_taxonomic_data/HACKETT_HH_herbivores_taxonomy.xlsx")
TAXOdata=as.matrix(TAXOdata)
TAXOdata=TAXOdata[match(species,TAXOdata[,1]),]
rownames(TAXOdata)=TAXOdata[,1]
TAXOdata=TAXOdata[,-5]
TAXOdata=cbind(TAXOdata, Class=rep("Insecta",12))
# Hypermodules (0= same in both interactions, .5= different in one interaction, 1= different in both interactions)
MODa=paste(MODNETa$Row_labels[match(species,rownames(NETa))],"a", sep="")
MODb=paste(MODNETb$Row_labels[match(species,rownames(NETb))],"b", sep="")
HYPERMOD=data.frame(row.names = species)
for (M in 1:nrow(CONGRUENCE$M1M2)){
  HYPERMOD$NETa[MODa==rownames(CONGRUENCE$M1M2)[M]]=unique(CONGRUENCE$M1M2[M,CONGRUENCE$M1M2[M,]!=0])
}
for (M in 1:ncol(CONGRUENCE$M1M2)){
  HYPERMOD$NETb[MODb==colnames(CONGRUENCE$M1M2)[M]]=unique(CONGRUENCE$M1M2[CONGRUENCE$M1M2[,M]!=0,M])
}
# Distances
TAXODIST=tax2dist(TAXOdata,tax_distance = c(Species=0, Genus=1, Family=2, Order=3, Class=4))
TAXODIST=TAXODIST@distance
GEODIST=vegdist(GEOdata,method = "jaccard",binary = T)
HYPERMOD_DIST=matrix(1, nrow=length(species), ncol=length(species), dimnames = list(species,species))
for (SP1 in 1:nrow(HYPERMOD)){
  for(SP2 in 1: nrow(HYPERMOD)){
    if(HYPERMOD$NETa[SP1]==HYPERMOD$NETa[SP2]){
      HYPERMOD_DIST[SP1,SP2]=HYPERMOD_DIST[SP1,SP2]-.5}
    if(HYPERMOD$NETb[SP1]==HYPERMOD$NETb[SP2]){
      HYPERMOD_DIST[SP1,SP2]=HYPERMOD_DIST[SP1,SP2]-.5}
  }
}
# Mantel Tests
Manteltests$dataset30$partitionB$TAXON_HM=
  mantel.partial(xdis=TAXODIST,ydis=HYPERMOD_DIST,zdis=GEODIST,method = "spearman",permutations = 10)
Manteltests$dataset30$partitionB$DISTR_HM= mantel.partial(xdis=GEODIST,ydis=HYPERMOD_DIST,zdis=TAXODIST,method = "spearman",permutations = 10)
## Table of results
MantelResults=rbind(MantelResults, data.frame(row.names = "dataset30_partitionB",
  TAXON_HM_r=Manteltests$dataset30$partitionB$TAXON_HM$statistic, 
  TAXON_HM_p=Manteltests$dataset30$partitionB$TAXON_HM$signif, 
  DISTR_HM_r=Manteltests$dataset30$partitionB$DISTR_HM$statistic,
  DISTR_HM_p=Manteltests$dataset30$partitionB$DISTR_HM$signif, 
  nspecies=length(species),
  permutations=Manteltests$dataset30$partitionB$TAXON_HM$permutations))
# # # # # # #
## partition C ##
Manteltests$dataset30$partitionC=list()
species=colnames(NETb)[paste(MODNETb$Col_labels,"b",sep="")%in%colnames(CONGRUENCE$M1M2)]
# GEO
GEOdata= read_excel("spatial_taxonomic_data/HACKETT_HH_herbivore.parasitoidsxsites.xlsx")
GEOdata=as.matrix(GEOdata)
GEOdata[is.na(GEOdata)]=0
rownames(GEOdata)=GEOdata[,1]
GEOdata=GEOdata[,-1]
mode(GEOdata)<-"numeric"
GEOdata=GEOdata[match(species,rownames(GEOdata)),]
# TAXONOMY
TAXOdata= read.csv("spatial_taxonomic_data/HACKETT_HH_herbivore.parasitoid taxonomy.csv")
TAXOdata=as.matrix(TAXOdata)
TAXOdata=TAXOdata[match(species,TAXOdata[,1]),]
rownames(TAXOdata)=TAXOdata[,1]
# missing data
# Missing genus -> consider each as a different genus
TAXOdata[(is.na(TAXOdata[,2])|TAXOdata[,2]=="NA"|TAXOdata[,2]==""),2]=paste("NA",1:sum(is.na(TAXOdata[,2])|TAXOdata[,2]=="NA"|TAXOdata[,2]==""),sep="")
# Missing family -> remove species
MISSING_FAMILY=is.na(TAXOdata[,3])|TAXOdata[,3]=="NA"|TAXOdata[,3]==""
species=species[!MISSING_FAMILY]
TAXOdata=TAXOdata[!MISSING_FAMILY,]
GEOdata=GEOdata[!MISSING_FAMILY,]
# Hypermodules (0= same, 1= different)
MODb=paste(MODNETb$Col_labels[match(species,colnames(NETb))],"b", sep="")
HYPERMOD=data.frame(row.names = species)
for (M in 1:ncol(CONGRUENCE$M1M2)){
  HYPERMOD$NETb[MODb==colnames(CONGRUENCE$M1M2)[M]]=unique(CONGRUENCE$M1M2[CONGRUENCE$M1M2[,M]!=0,M])
}
# Distances
TAXODIST=tax2dist(TAXOdata,tax_distance = c(Species=0, Genus=1, Family=2, Order=3, Kingdom=4))
TAXODIST=TAXODIST@distance
GEODIST=vegdist(GEOdata,method = "jaccard",binary = T)
HYPERMOD_DIST=matrix(NA, nrow=length(species), ncol=length(species), dimnames = list(species,species))
for (SP1 in 1:nrow(HYPERMOD)){
  for(SP2 in 1: nrow(HYPERMOD)){
    if(HYPERMOD$NETb[SP1]==HYPERMOD$NETb[SP2]){
      HYPERMOD_DIST[SP1,SP2]=0
    }else{HYPERMOD_DIST[SP1,SP2]=1}
  }
}
# Mantel Tests
Manteltests$dataset30$partitionC$TAXON_HM=
  mantel.partial(xdis=TAXODIST,ydis=HYPERMOD_DIST,zdis=GEODIST,method = "spearman",permutations = 10)
Manteltests$dataset30$partitionC$DISTR_HM= mantel.partial(xdis=GEODIST,ydis=HYPERMOD_DIST,zdis=TAXODIST,method = "spearman",permutations = 10)
## Table of results
MantelResults=rbind(MantelResults, data.frame(row.names = "dataset30_partitionC",
  TAXON_HM_r=Manteltests$dataset30$partitionC$TAXON_HM$statistic, 
  TAXON_HM_p=Manteltests$dataset30$partitionC$TAXON_HM$signif, 
  DISTR_HM_r=Manteltests$dataset30$partitionC$DISTR_HM$statistic,
  DISTR_HM_p=Manteltests$dataset30$partitionC$DISTR_HM$signif, 
  nspecies=length(species),
  permutations=Manteltests$dataset30$partitionC$TAXON_HM$permutations))
save(MantelResults,Manteltests, file="files_results/Mantel_results.RData")
# # # # # # #
### DATASET 32 ####
datasetID=32
Manteltests$dataset32=list()
load(paste("files_results/dataset_",datasetID,".RData",sep=""))
# # # # # # #
## partition A ##
Manteltests$dataset32$partitionA=list()
species=colnames(NETa)[paste(MODNETa$Col_labels,"a",sep="")%in%rownames(CONGRUENCE$M1M2)]
# GEO
GEOdata= read_excel("spatial_taxonomic_data/SHINOHARA_herbivoresxsites.xlsx",sheet = "Planilha1")
GEOdata=as.matrix(GEOdata)
GEOdata[is.na(GEOdata)]=0
rownames(GEOdata)=paste("Herbivore",GEOdata[,1],sep="") ##
GEOdata=GEOdata[,-1]
mode(GEOdata)<-"numeric"
GEOdata=GEOdata[match(species,rownames(GEOdata)),]
# Hypermodules (0= same, 1= different)
MODa=paste(MODNETa$Col_labels[match(species,colnames(NETa))],"a", sep="")
HYPERMOD=data.frame(row.names = species)
for (M in 1:nrow(CONGRUENCE$M1M2)){
  HYPERMOD$NETa[MODa==rownames(CONGRUENCE$M1M2)[M]]=unique(CONGRUENCE$M1M2[M,CONGRUENCE$M1M2[M,]!=0])
}
# Distances
GEODIST=vegdist(GEOdata,method = "jaccard",binary = T)
HYPERMOD_DIST=matrix(NA, nrow=length(species), ncol=length(species), dimnames = list(species,species))
for (SP1 in 1:nrow(HYPERMOD)){
  for(SP2 in 1: nrow(HYPERMOD)){
    if(HYPERMOD$NETa[SP1]==HYPERMOD$NETa[SP2]){
      HYPERMOD_DIST[SP1,SP2]=0
    }else{HYPERMOD_DIST[SP1,SP2]=1}
  }
}
# Mantel Tests
Manteltests$dataset32$partitionA$TAXON_HM=NULL
Manteltests$dataset32$partitionA$DISTR_HM= mantel(xdis=GEODIST,ydis=HYPERMOD_DIST,method = "spearman",permutations = 10)
## Table of results
MantelResults=rbind(MantelResults, data.frame(row.names = "dataset32_partitionA",
  TAXON_HM_r=NA, 
  TAXON_HM_p=NA, 
  DISTR_HM_r=Manteltests$dataset32$partitionA$DISTR_HM$statistic,
  DISTR_HM_p=Manteltests$dataset32$partitionA$DISTR_HM$signif, 
  nspecies=length(species),
  permutations=Manteltests$dataset32$partitionA$DISTR_HM$permutations))
# # # # # # # #
## partition B ##
Manteltests$dataset32$partitionB=list()
species=rownames(NETa)[rownames(NETa)%in%rownames(NETb)]
# GEO
GEOdata= read_excel("spatial_taxonomic_data/SHINOHARA_plantsxsites.xlsx")
GEOdata=as.matrix(GEOdata)
GEOdata[is.na(GEOdata)]=0
rownames(GEOdata)=paste("Plant",GEOdata[,1],sep="")
GEOdata=GEOdata[,-1]
mode(GEOdata)<-"numeric"
GEOdata=GEOdata[match(species,rownames(GEOdata)),]
# TAXONOMY
TAXOdata= read_excel("spatial_taxonomic_data/SHINOHARA_taxonomy.xlsx")
TAXOdata=as.matrix(TAXOdata)
TAXOdata[,1]=paste("Plant",as.numeric(TAXOdata[,1]),sep="")
TAXOdata=TAXOdata[match(species,TAXOdata[,1]),]
rownames(TAXOdata)=TAXOdata[,1]
colnames(TAXOdata)[1]="Species"
# missing data
# Missing genus -> consider each as a different genus
TAXOdata[(is.na(TAXOdata[,2])|TAXOdata[,2]=="NA"),2]=paste("NA",1:sum(is.na(TAXOdata[,2])|TAXOdata[,2]=="NA"),sep="")
# Missing family -> remove species
MISSING_FAMILY=is.na(TAXOdata[,3])|TAXOdata[,3]=="NA"
species=species[!MISSING_FAMILY]
TAXOdata=TAXOdata[!MISSING_FAMILY,]
GEOdata=GEOdata[!MISSING_FAMILY,]
TAXOdata=cbind(TAXOdata, data.frame(Kingdom="Plantae"))
# Hypermodules (0= same in both interactions, .5= different in one interaction, 1= different in both interactions)
MODa=paste(MODNETa$Row_labels[match(species,rownames(NETa))],"a", sep="")
MODb=paste(MODNETb$Row_labels[match(species,rownames(NETb))],"b", sep="")
HYPERMOD=data.frame(row.names = species)
for (M in 1:nrow(CONGRUENCE$M1M2)){
  HYPERMOD$NETa[MODa==rownames(CONGRUENCE$M1M2)[M]]=unique(CONGRUENCE$M1M2[M,CONGRUENCE$M1M2[M,]!=0])
}
for (M in 1:ncol(CONGRUENCE$M1M2)){
  HYPERMOD$NETb[MODb==colnames(CONGRUENCE$M1M2)[M]]=unique(CONGRUENCE$M1M2[CONGRUENCE$M1M2[,M]!=0,M])
}
# Distances
TAXODIST=tax2dist(TAXOdata,tax_distance = c(Species=0, Genus=1, Family=2, Order=3, Kingdom=4))
TAXODIST=TAXODIST@distance
GEODIST=vegdist(GEOdata,method = "jaccard",binary = T)
HYPERMOD_DIST=matrix(1, nrow=length(species), ncol=length(species), dimnames = list(species,species))
for (SP1 in 1:nrow(HYPERMOD)){
  for(SP2 in 1: nrow(HYPERMOD)){
    if(HYPERMOD$NETa[SP1]==HYPERMOD$NETa[SP2]){
      HYPERMOD_DIST[SP1,SP2]=HYPERMOD_DIST[SP1,SP2]-.5}
    if(HYPERMOD$NETb[SP1]==HYPERMOD$NETb[SP2]){
      HYPERMOD_DIST[SP1,SP2]=HYPERMOD_DIST[SP1,SP2]-.5}
  }
}
# Mantel Tests
Manteltests$dataset32$partitionB$TAXON_HM=
  mantel.partial(xdis=TAXODIST,ydis=HYPERMOD_DIST,zdis=GEODIST,method = "spearman",permutations = 10)
Manteltests$dataset32$partitionB$DISTR_HM= mantel.partial(xdis=GEODIST,ydis=HYPERMOD_DIST,zdis=TAXODIST,method = "spearman",permutations = 10)
## Table of results
MantelResults=rbind(MantelResults, data.frame(row.names = "dataset32_partitionB",
  TAXON_HM_r=Manteltests$dataset32$partitionB$TAXON_HM$statistic, 
  TAXON_HM_p=Manteltests$dataset32$partitionB$TAXON_HM$signif, 
  DISTR_HM_r=Manteltests$dataset32$partitionB$DISTR_HM$statistic,
  DISTR_HM_p=Manteltests$dataset32$partitionB$DISTR_HM$signif, 
  nspecies=length(species),
  permutations=Manteltests$dataset32$partitionB$TAXON_HM$permutations))
# # # # # # #
## partition C ##
Manteltests$dataset32$partitionC=list()
species=colnames(NETb)[paste(MODNETb$Col_labels,"b",sep="")%in%colnames(CONGRUENCE$M1M2)]
# GEO
GEOdata= read_excel("spatial_taxonomic_data/SHINOHARA_pollinatorsxsites.xlsx")
GEOdata=as.matrix(GEOdata)
GEOdata[is.na(GEOdata)]=0
rownames(GEOdata)=paste("Pollinator",GEOdata[,1], sep="")
GEOdata=GEOdata[,-1]
mode(GEOdata)<-"numeric"
GEOdata=GEOdata[match(species,rownames(GEOdata)),]
# Hypermodules (0= same, 1= different)
MODb=paste(MODNETb$Col_labels[match(species,colnames(NETb))],"b", sep="")
HYPERMOD=data.frame(row.names = species)
for (M in 1:ncol(CONGRUENCE$M1M2)){
  HYPERMOD$NETb[MODb==colnames(CONGRUENCE$M1M2)[M]]=unique(CONGRUENCE$M1M2[CONGRUENCE$M1M2[,M]!=0,M])
}
# Distances
GEODIST=vegdist(GEOdata,method = "jaccard",binary = T)
HYPERMOD_DIST=matrix(NA, nrow=length(species), ncol=length(species), dimnames = list(species,species))
for (SP1 in 1:nrow(HYPERMOD)){
  for(SP2 in 1: nrow(HYPERMOD)){
    if(HYPERMOD$NETb[SP1]==HYPERMOD$NETb[SP2]){
      HYPERMOD_DIST[SP1,SP2]=0
    }else{HYPERMOD_DIST[SP1,SP2]=1}
  }
}
# Mantel Tests
Manteltests$dataset32$partitionC$TAXON_HM=NULL
Manteltests$dataset32$partitionC$DISTR_HM= mantel(xdis=GEODIST,ydis=HYPERMOD_DIST,method = "spearman",permutations = 10)
## Table of results
MantelResults=rbind(MantelResults, data.frame(row.names = "dataset32_partitionC",
  TAXON_HM_r=NA, 
  TAXON_HM_p=NA, 
  DISTR_HM_r=Manteltests$dataset32$partitionC$DISTR_HM$statistic,
  DISTR_HM_p=Manteltests$dataset32$partitionC$DISTR_HM$signif, 
  nspecies=length(species),
  permutations=Manteltests$dataset32$partitionC$DISTR_HM$permutations))
save(MantelResults,Manteltests, file="files_results/Mantel_results.RData")
# # # # # # #