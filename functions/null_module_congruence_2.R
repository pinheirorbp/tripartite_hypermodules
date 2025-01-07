null_module_congruence2= function(NETa=NULL,NETb=NULL,MODNETa=NULL, MODNETb= NULL, LAMBDA=0.9995,NSTEPS_MAX=50000,NSTEPS_STOP= 5000,TEMP0=4, NCores=1, Nulls=1){
  
  library(doSNOW)
  ## Filtering = species that are in both networks
  INT_SPEC=rownames(NETa)[is.element(rownames(NETa),rownames(NETb))]
  OSNETa=match(INT_SPEC, rownames(NETa))
  OSNETb=match(INT_SPEC, rownames(NETb))
  # modules per species
  M1_real=MODNETa$Row_labels[OSNETa]
  M1_real=paste(M1_real,"a", sep="")
  M2_real=MODNETb$Row_labels[OSNETb]
  M2_real=paste(M2_real,"b", sep="")
  # cleaning the environment
  rm(NETa,NETb,MODNETa,MODNETb)
  # setting cluster
  cl <- makeCluster(NCores, outfile="")
  registerDoSNOW(cl)
  
  NULL_CONGRUENCE=foreach(j=1:Nulls,.verbose = T,.export = "module_congruence")%dopar%{
    
    # randomizing
    M1=sample(M1_real)
    M2=sample(M2_real)
    CONGRUENCE= module_congruence(MODa = M1, MODb= M2, LAMBDA = LAMBDA, NSTEPS_MAX = NSTEPS_MAX, plot_opt = F, printR = F,printR_each = F,NSTEPS_STOP = NSTEPS_STOP, TEMP0 = TEMP0, NETa = NULL, NETb = NULL, MODNETa = NULL, MODNETb = NULL)
    list(CONGRUENCE$finalcongruence,CONGRUENCE$finalOPT)
  }
  stopCluster(cl)
  ## Reorganizing the results
  Nullcongruence=numeric()
  NullOPT=numeric()
  for (I in 1:Nulls){
    Nullcongruence[I]=NULL_CONGRUENCE[[I]][[1]]
    NullOPT[I]=NULL_CONGRUENCE[[I]][[2]]
  }
  list(Nullcongruence=Nullcongruence, NullOPT=NullOPT)
}