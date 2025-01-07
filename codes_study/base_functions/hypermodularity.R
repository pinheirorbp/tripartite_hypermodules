hypermodularity= function(M1M2=NULL,MODNETa=NULL, MODNETb=NULL, NETa=NULL, NETb=NULL){
  
  # Hypermodules
  # NETa
  HNETa=data.frame(MOD=rownames(M1M2), 
                   HIPERMOD=apply(M1M2,MARGIN = 1, 
                                  FUN = function (x) {
                                    unique(x)[unique(x)!=0]}))
  HIPNETa_row=HNETa$HIPERMOD[match(paste(MODNETa$Row_labels,"a",sep=""),HNETa$MOD)]
  HIPNETa_col=HNETa$HIPERMOD[match(paste(MODNETa$Col_labels,"a",sep=""),HNETa$MOD)]
  HIPNETa_row2=data.frame(HIPER=HIPNETa_row[is.element(rownames(NETa),rownames(NETb))], species=rownames(NETa)[is.element(rownames(NETa),rownames(NETb))])
  NETa2=NETa[is.element(rownames(NETa),rownames(NETb)),]
  NETa2=NETa2[,colSums(NETa2)>0]
  HIPNETa_col2=HIPNETa_col[match(colnames(NETa2),colnames(NETa))]
  HIPNETa_col2[is.na(HIPNETa_col2)]=max(HIPNETa_col2, na.rm = T)+10
  # NETb
  HNETb=data.frame(MOD=colnames(M1M2), 
                   HIPERMOD=apply(M1M2,MARGIN = 2, 
                                  FUN = function (x) {
                                    unique(x)[unique(x)!=0]}))
  HIPNETb_row=HNETb$HIPERMOD[match(paste(MODNETb$Row_labels,"b",sep=""),HNETb$MOD)]
  HIPNETb_col=HNETb$HIPERMOD[match(paste(MODNETb$Col_labels,"b",sep=""),HNETb$MOD)]
  HIPNETb_row2=data.frame(HIPER=HIPNETb_row[is.element(rownames(NETb),rownames(NETa))], species=rownames(NETb)[is.element(rownames(NETb),rownames(NETa))])
  NETb2=NETb[is.element(rownames(NETb),rownames(NETa)),]
  NETb2=NETb2[,colSums(NETb2)>0]
  HIPNETb_col2=HIPNETb_col[match(colnames(NETb2),colnames(NETb))]
  HIPNETb_col2[is.na(HIPNETb_col2)]=max(HIPNETb_col2, na.rm = T)+10
  
  #crossed hypermodularity of NETa
  crossedHIPNETa=HIPNETb_row2$HIPER[match(rownames(NETa2), HIPNETb_row2$species)]
  NETa3=NETa2/sum(NETa2)
  expNETa3=outer(rowSums(NETa3),colSums(NETa3), FUN = function(x,y){x*y})
  intraHIPNETa=outer(crossedHIPNETa, HIPNETa_col2, FUN=function(x,y){x==y})
  hipermod_NETa=sum((NETa3-expNETa3)[intraHIPNETa])

  
  #crossed hypermodularity of NETb
  crossedHIPNETb=HIPNETa_row2$HIPER[match(rownames(NETb2), HIPNETa_row2$species)]
  NETb3=NETb2/sum(NETb2)
  expNETb3=outer(rowSums(NETb3),colSums(NETb3), FUN = function(x,y){x*y})
  intraHIPNETb=outer(crossedHIPNETb, HIPNETb_col2, FUN=function(x,y){x==y})
  hipermod_NETb=sum((NETb3-expNETb3)[intraHIPNETb])
  
  # General hipermodularity
  hipermod=mean(c(hipermod_NETa,hipermod_NETb))
  list(overall= hipermod, NETa=hipermod_NETa, NETb=hipermod_NETb)
}