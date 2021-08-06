for(h in 1:38){
  i = Ind[h,1]; j = Ind[h,2]
  if(i<7){
    ind = Ind[Ind[,1]==i,2]
    BES = BE; LAS = LA; PPS = PP; ORS = OR # BES : current proposed beta's / BE : previous beta's
    # LACS = LAC; ORCS = ORC;
    BES[i,j,] = BE[i,j,]+rnorm(nbe,0,Tau[i,j])
    LAS[i,j,2:TT] = exp(XX%*%BES[i,j,])
    # LACS[i,j,2:(TT+tahead)] = exp(XXC%*%BES[i,j,])
    for(t in 2:TT) PPS[i,,t] = LAS[i,,t]/sum(LAS[i,,t])
    for(t in 2:TT) ORS[i,,t] = LAS[i,,t]/LAS[i,i,t]
    check = TRUE
    if(conbe) if(BES[i,j,3]>0) check=FALSE
    for(j1 in ind) if(j1!=i){
      if(!is.na(ORmin[i,j1])) if(any(ORS[i,j1,2:TT]<ORmin[i,j1])) check = FALSE
      if(!is.na(ORmax[i,j1])) if(any(ORS[i,j1,2:TT]>ORmax[i,j1])) check = FALSE
    }
    if(check){ 
      tmp = ldnorm1(BES[i,j,],BE[i,j,],si2)
  
      for(t in 2:TT){
        La = LA[,,t]; Las = LAS[,,t]
        tla = sum(La[i,ind]); tlas = sum(Las[i,ind])
        ntab = sum(TAB[i,ind,t])
        tmp = tmp + lgamma(tlas)-lgamma(ntab+tlas)+sum(lgamma(TAB[i,ind,t]+Las[i,ind])-lgamma(Las[i,ind]))-
          lgamma(tla)+lgamma(ntab+tla)-sum(lgamma(TAB[i,ind,t]+La[i,ind])-lgamma(La[i,ind]))
      }
      al = exp(tmp)
      if(is.nan(al)){
        print("NA in updating Beta")
        browser()
      }
      if(runif(1)<=al){
        BE = BES; LA = LAS; OR = ORS; PP = PPS
        accbe = accbe + 1/37
        Accbe[i,j] = Accbe[i,j]+1
      }
    }
  }
}
