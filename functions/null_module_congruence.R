null_module_congruence= function(NETa=NULL,NETb=NULL,MODNETa=NULL, MODNETb= NULL,M1M2= NULL, NCores=1, Nulls=100){
  
  library(doSNOW)
  ## Filtering = species that are in both networks
  INT_SPEC=rownames(NETa)[is.element(rownames(NETa),rownames(NETb))]
  OSNETa=match(INT_SPEC, rownames(NETa))
  OSNETb=match(INT_SPEC, rownames(NETb))
  # module of each species
  M1_real=MODNETa$Row_labels[OSNETa]
  M1_real=paste(M1_real,"a", sep="")
  M2_real=MODNETb$Row_labels[OSNETb]
  M2_real=paste(M2_real,"b", sep="")
  # Hipermodule of each species
  HIPERM1_real=rep(NA,length(M1_real))
  for (M in 1:nrow(M1M2)){
    HIPERM1_real[M1_real==rownames(M1M2)[M]]=unique(M1M2[M,M1M2[M,]!=0])
  }
  HIPERM1_real=as.character(HIPERM1_real)
  HIPERM2_real=rep(NA,length(M2_real))
  for (M in 1:ncol(M1M2)){
    HIPERM2_real[M2_real==colnames(M1M2)[M]]=unique(M1M2[M1M2[,M]!=0,M])
  }
  HIPERM2_real=as.character(HIPERM2_real)
  # Calculating the expected congruence given sizes of hipermodules
  base_adj= function(M1M2=NULL, M1=NULL, M2=NULL){
    M1_comb2=numeric()
    for (M1n in 1:length(M1)){
      X1=unique(M1M2[M1[M1n],])
      M1_comb2[M1n]=X1[X1>0]
    }
    M2_comb2=numeric()
    for (M2n in 1:length(M2)){
      X2=unique(M1M2[,M2[M2n]])
      M2_comb2[M2n]=X2[X2>0]
    }
    X3=numeric()
    for (M12n in unique(M1_comb2)){
      X3[M12n]=sum(M1_comb2==M12n)/length(M1)*sum(M2_comb2==M12n)/length(M2)
    }
    base_adj=sum(X3,na.rm = T)
    base_adj
  }
  BASE=base_adj(M1M2 = M1M2, M1=M1_real, M2=M2_real)
  Nsp=length(HIPERM1_real)
  # cleaning the environment
  rm(NETa,NETb,MODNETa,MODNETb)
  # setting cluster
  cl <- makeCluster(NCores, outfile="")
  registerDoSNOW(cl)
  
  NULL_CONGRUENCE=foreach(j=1:Nulls,.verbose = T)%dopar%{
    
    # randomizing
    HIPERM1=sample(HIPERM1_real)
    HIPERM2=sample(HIPERM2_real)
    # congruence
    ADJ=rep(NA,Nsp)
    for (i in 1:Nsp){
      ADJ[i]=HIPERM1[i]==HIPERM2[i]
    }
    congruence=mean(ADJ)
    OPT=congruence-BASE
    list(congruence,OPT)
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