module_congruence=function (NETa=NULL,NETb=NULL,MODNETa=NULL, MODNETb= NULL, LAMBDA=0.9995,NSTEPS_MAX=50000,plot_opt=T, printR=T,printR_each=1000,
  NSTEPS_STOP= 5000,TEMP0=4, MODa=NULL, MODb=NULL){
  
  if(is.null(MODa)|is.null(MODb)){
    # error
    if(is.null(NETa)|is.null(NETb)|is.null(MODNETa)|is.null(MODNETb)){
    stop("must provide [NETa,NETb,MODNETa and MODNETb] or [MODa and MODb]")}
    ## Filtering = species that are in both networks
    INT_SPEC=rownames(NETa)[is.element(rownames(NETa),rownames(NETb))]
    OSNETa=match(INT_SPEC, rownames(NETa))
    OSNETb=match(INT_SPEC, rownames(NETb))
    # modules per species
    M1=MODNETa$Row_labels[OSNETa]
    M1=paste(M1,"a", sep="")
    M2=MODNETb$Row_labels[OSNETb]
    M2=paste(M2,"b", sep="")
  }else{
    M1=MODa
    M2=MODb
  }
  
  Nsp=length(M1)
  # LIST OF MODULES IN EACH NET
  M1L=unique(M1)
  M2L=unique(M2)
  M1L_N=length(M1L)
  M2L_N=length(M2L)
  # MATRIX MERGING MODULES IN HIPERMODULES
  M1M2_Mat0=matrix(0, nrow=M1L_N, ncol=M2L_N, dimnames = list(M1L, M2L))
  # Initial - Arbitrary
  if(M1L_N>M2L_N){
    for (i in 1:M1L_N){M1M2_Mat0[i,(i-1)%%(M2L_N)+1]=(i-1)%%(M2L_N)+1}
  } else{
    for (i in 1:M2L_N){M1M2_Mat0[(i-1)%%(M1L_N)+1,i]=(i-1)%%(M1L_N)+1}
  }
  # hipermodule for each separate module in M1M2_Mat0
  M1_comb=numeric()
  for (M1nn in 1:nrow(M1M2_Mat0)){
    X1=unique(M1M2_Mat0[M1nn,])
    if(sum(X1>0)!=1){stop()}
    M1_comb[M1nn]=X1[X1>0]
  }
  M2_comb=numeric()
  for (M2nn in 1:ncol(M1M2_Mat0)){
    X1=unique(M1M2_Mat0[,M2nn])
    if(sum(X1>0)!=1){stop()}
    M2_comb[M2nn]=X1[X1>0]
  }
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
  # Optimization
  ADJ=rep(NA,Nsp)
  for (i in 1:Nsp){
    ADJ[i]=M1M2_Mat0[M1[i],M2[i]]>0
  }
  congruence=numeric()
  congruence[1]=mean(ADJ)
  BASE=numeric()
  BASE[1]=base_adj(M1M2 = M1M2_Mat0, M1=M1, M2=M2)
  OPT=numeric()
  OPT[1]=congruence[1]-BASE[1]
  # SIMULATED ANNEALING ###
  ###########
  TEMP=numeric()
  TEMP[1]=TEMP0
  M1M2=numeric()
  M1M2=M1M2_Mat0
  R=2
  STOP=F
  Rlastchange=2
  while(R <=NSTEPS_MAX & !STOP){
    # for(R in 2:NSTEPS){
    newM1M2_Mat=M1M2
    
    ## Modification in hipermodules
    if(sample(c(F,T),1)){
    newM1M2_Mat[sample(1:M1L_N,1),]=0 # 50% disconnect a row
    }else{newM1M2_Mat[,sample(1:M2L_N,1)]=0} # 50% disconnect a column}
    
    # run until no row or column is disconnected
      while(any(c(colSums(newM1M2_Mat)==0,rowSums(newM1M2_Mat)==0))){
        # hipermodule of each module in updated newM1M2_Mat
        for (M1nn in 1:nrow(newM1M2_Mat)){
          X1=unique(newM1M2_Mat[M1nn,])
          M1_comb[M1nn]=ifelse(sum(X1)==0,0,X1[X1>0])
        }
        for (M2nn in 1:ncol(newM1M2_Mat)){
          X1=unique(newM1M2_Mat[,M2nn])
          M2_comb[M2nn]=ifelse(sum(X1)==0,0,X1[X1>0])
        }
        # connecting an empty row
        emptyrows=rowSums(newM1M2_Mat)==0
        emptycols=colSums(newM1M2_Mat)==0
        if(any(emptyrows)){ 
          RN=ifelse(sum(emptyrows)==1,which(emptyrows),sample(which(emptyrows)))
          CN1=sample(1:M2L_N,1)
          emptyrowsN=sum(emptyrows)
          probF=(M1L_N/(M1L_N+M2L_N))^(1/emptyrowsN)
          if(sample(c(T,F),1, prob = c(1-probF,probF))|emptycols[CN1]){
            newM1M2_Mat[,CN1]=0
            newM1M2_Mat[RN,CN1]=which(!is.element(1:max(M1L_N,M2L_N),newM1M2_Mat))[1]
          }else{
            newM1M2_Mat[RN,M2_comb==M2_comb[CN1]]=M2_comb[CN1]
          }
        }
        # hipermodule of each module in the updated newM1M2_Mat
        for (M1nn in 1:nrow(newM1M2_Mat)){
          X1=unique(newM1M2_Mat[M1nn,])
          M1_comb[M1nn]=ifelse(sum(X1)==0,0,X1[X1>0])
        }
        for (M2nn in 1:ncol(newM1M2_Mat)){
          X1=unique(newM1M2_Mat[,M2nn])
          M2_comb[M2nn]=ifelse(sum(X1)==0,0,X1[X1>0])
        }
        emptyrows=rowSums(newM1M2_Mat)==0
        emptycols=colSums(newM1M2_Mat)==0
        # connecting empty columns
        if(any(emptycols)){
          CN=ifelse(sum(emptycols)==1,which(emptycols),sample(which(emptycols)))
          RN1=sample(1:M1L_N,1)
          emptycolsN=sum(emptycols)
          probF=(M2L_N/(M1L_N+M2L_N))^(1/emptycolsN)
          if(sample(c(T,F),1, prob = c(1-probF,probF))|emptyrows[RN1]){
            newM1M2_Mat[RN1,]=0
            newM1M2_Mat[RN1,CN]=which(!is.element(1:max(M1L_N,M2L_N),newM1M2_Mat))[1]
          }else{
            newM1M2_Mat[M1_comb==M1_comb[RN1],CN]=M1_comb[RN1]}
        }
      }

    # Opitmization
    ADJ=rep(NA,Nsp)
    for (i in 1:Nsp){
      ADJ[i]=newM1M2_Mat[M1[i],M2[i]]>0
    }
    newcongruence=mean(ADJ)
    newbase_ADJ=base_adj(M1M2 = newM1M2_Mat, M1=M1, M2=M2)
    newOPT=newcongruence-newbase_ADJ
    # ACCEPTANCE / REJECTION #
    CHANGE= newOPT-OPT[R-1] 
    if(CHANGE>0){Rlastchange=R}
    if(CHANGE>=0){ 
      OPT[R]=newOPT
      congruence[R]=newcongruence
      BASE[R]=newbase_ADJ
      M1M2=newM1M2_Mat #if the change is positive, modify
    }else{
      PROB= exp(CHANGE/TEMP[R-1]) # Probability of acceptance, if CHANGE is negative
      if(PROB>1){PROB=1}
      if(rbinom(1,1,PROB)){
        OPT[R]=newOPT
        congruence[R]=newcongruence
        BASE[R]=newbase_ADJ
        M1M2=newM1M2_Mat
      }else{
        OPT[R]=OPT[R-1]
        congruence[R]=congruence[R-1]
        BASE[R]=BASE[R-1]
      }
    }
    TEMP[R]=TEMP[R-1]*LAMBDA
    R=R+1
    if(printR&R%%printR_each==0){print(R)}
    STOP=(R-Rlastchange)>=NSTEPS_STOP
  }
  # Plots
  if(plot_opt){
    plot(y=OPT[seq.int(from=500, to=R, by=100)],x=seq.int(from=500, to=R, by=100), pch=16, cex=.7, ylab="opt", xlab=NA)
    plot(y=congruence[seq.int(from=500, to=R, by=100)],x=seq.int(from=500, to=R, by=100), pch=16, cex=.7, ylab="congruence", xlab=NA)
    plot(y=BASE[seq.int(from=500, to=R, by=100)],x=seq.int(from=500, to=R, by=100), pch=16, cex=.7, ylab="base_adj", xlab=NA)
    }
  list(finalcongruence=congruence[length(congruence)],finalOPT=OPT[length(OPT)], M1M2=M1M2, OPT=OPT, congruence=congruence, BASE=BASE)
}