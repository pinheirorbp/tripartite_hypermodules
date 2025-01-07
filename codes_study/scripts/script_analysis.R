#
B=max(NETa)==1|max(NETb)==1

##NETa : Plant-Herbivore

# Modularity in NETa
MODNETa=DIRT_LPA_wb_plus(NETa)
# Null analysis
if(!B){
  # Proportional null model
  propNETa= propnull2(NETa, N = N_null_topology)
  # Cluster for parallel computing
  cl <- makeCluster(NCores)
  registerDoSNOW(cl)
  # Performing DIRT LPA wb+ for each null matrix of NETa
  modnull_NETa=foreach(j=1:N_null_topology,.verbose = T,.packages="bipartite",.combine="c")%dopar%{
    DIRT_LPA_wb_plus(propNETa[[j]])$modularity
  }
  stopCluster(cl)
}
if(B){ 
  propNETa=NULL
  modnull_NETa=NULL
}

# Low-level nestedness 
if(!B){
  Nest_SMDM_NETa=nest.smdm(NETa,constraints = c(MODNETa$Row_labels,MODNETa$Col_labels), weighted = T,decreasing = "abund")
}
if(B){
  Nest_SMDM_NETa=nest.smdm(NETa,constraints = c(MODNETa$Row_labels,MODNETa$Col_labels), weighted = F)
}
if(!B){
  # Restricted null model
  rest_NETa= restrictednull(NETa, Prior.Pij = "equiprobable",
                            conditional.level = "areas",N=N_null_topology,R.partitions = MODNETa$Row_labels,
                            C.partitions = MODNETa$Col_labels, byarea = T, connectance = F,allow.degeneration = F)
  # Cluster for parallel computing
  cl <- makeCluster(NCores)
  registerDoSNOW(cl)
  # WNODA SM for each null matrix of NETa
  NestSMnull_NETa=foreach(j=1:N_null_topology,.verbose = T,.packages="bipartite",.combine="c")%dopar%{
    nest.smdm(rest_NETa[[j]],constraints = c(MODNETa$Row_labels,MODNETa$Col_labels), weighted = T,decreasing = "abund")$WNODA_SM_matrix
  }
  stopCluster(cl)
}
if(B){
  rest_NETa=NULL
  NestSMnull_NETa=NULL
}
##NETb : Plant-Herbivore

# Modularity in NETb
MODNETb=DIRT_LPA_wb_plus(NETb)
# Null analysis
if(!B){
  # Proportional null model
  propNETb= propnull2(NETb, N = N_null_topology)
  # Cluster for parallel computing
  cl <- makeCluster(NCores)
  registerDoSNOW(cl)
  # Performing DIRT LPA wb+ for each null matrix of NETb
  modnull_NETb=foreach(j=1:N_null_topology,.verbose = T,.packages="bipartite",.combine="c")%dopar%{
    DIRT_LPA_wb_plus(propNETb[[j]])$modularity
  }
  stopCluster(cl)
}
if(B){ 
  propNETb=NULL
  modnull_NETb=NULL
}

# Low-level nestedness 
if(!B){
  Nest_SMDM_NETb=nest.smdm(NETb,constraints = c(MODNETb$Row_labels,MODNETb$Col_labels), weighted = T,decreasing = "abund")
}
if(B){
  Nest_SMDM_NETb=nest.smdm(NETb,constraints = c(MODNETb$Row_labels,MODNETb$Col_labels), weighted = F)
}
if(!B){
  # Restricted null model
  rest_NETb= restrictednull(NETb, Prior.Pij = "equiprobable",
                            conditional.level = "areas",N=N_null_topology,R.partitions = MODNETb$Row_labels,
                            C.partitions = MODNETb$Col_labels, byarea = T, connectance = F,allow.degeneration = F)
  # Cluster for parallel computing
  cl <- makeCluster(NCores)
  registerDoSNOW(cl)
  # WNODA SM for each null matrix of NETb
  NestSMnull_NETb=foreach(j=1:N_null_topology,.verbose = T,.packages="bipartite",.combine="c")%dopar%{
    nest.smdm(rest_NETb[[j]],constraints = c(MODNETb$Row_labels,MODNETb$Col_labels), weighted = T,decreasing = "abund")$WNODA_SM_matrix
  }
  stopCluster(cl)
}
if(B){
  rest_NETb=NULL
  NestSMnull_NETb=NULL
}
## Network-level metrics (Bipartite networks)
NETdata=list()
# Specialization in NETa (H2')
NETdata$NETa_H2=ifelse(!B,networklevel(NETa,index="H2"),NA)
# Connectance in NETa
NETdata$NETa_connectance=networklevel(NETa,index="connectance")
# Specialization in NETb (H2')
NETdata$NETb_H2=ifelse(!B,networklevel(NETb,index="H2"),NA)
# Connectance in NETb
NETdata$NETb_connectance=networklevel(NETb,index="connectance")
####

### module congruence ###
CONGRUENCE=module_congruence(NETa = NETa,NETb = NETb, MODNETa = MODNETa, MODNETb = MODNETb,plot_opt = F)

# Null model for hipermodule congruence
# Null1

nullCongruence= null_module_congruence(NETa = NETa, NETb = NETb, MODNETa = MODNETa, MODNETb = MODNETb, NCores = NCores, Nulls = N_null_congruence, M1M2 = CONGRUENCE$M1M2)

# Null2

nullCongruence2= null_module_congruence2(NETa = NETa, NETb = NETb, MODNETa = MODNETa, MODNETb = MODNETb, NCores = NCores, Nulls = N_null_congruence2)

# Hypermodularity

hypermod=hypermodularity(M1M2=CONGRUENCE$M1M2, MODNETa = MODNETa, MODNETb = MODNETb, NETa = NETa, NETb = NETb)

#

save(NETa, NETb,propNETa,modnull_NETa,MODNETa,Nest_SMDM_NETa,
  NestSMnull_NETa,rest_NETa, MODNETb,propNETb,
  modnull_NETb,Nest_SMDM_NETb,CONGRUENCE,rest_NETb,
  NestSMnull_NETb,nullCongruence,nullCongruence2, NETdata,dataset_info, hypermod, file=paste("files_results/dataset_",datasetID,".RData",sep=""))