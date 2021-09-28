library(pbmcapply)
library(splines)
library(extraDistr)


#'
#' MCMC for Dirichlet-multinomial AR model
#' TT total time
#' 
#' Categories must be ordered as
#' "susceptible", "recovered", "isolated",  "hospitalized", "threatened", "undetected", "extinct" 
#' 
#' @param Y matrix (TT x 7) of observed daily frequencies (TT = #days) for the 6 categories ordered from susceptibles to deceased
#' @param R number of iterations after burning in
#' @param burnin number of iteration of burning in
#' @param tin tinning step (store result each tinning iteration)
#' @param Tau matrix of standard deviations for each proposal of regression parameters (beta's)
#' @param mra maximum variation of each cell in updating the tables
#' @param si2BE variance of regression parameters (beta's)
#' @param degree maximum degree of time polynomial
#' @param tint times of intervention (dummies covariates)
#' @param disp display tables during estimation
#' @param ORmin matrix of lower limits of OR
#' @param ORmax matrix of upper limits of OR
#' 
#' 
#' 
#' @return a list with the following components:
#'         TTAB array of transition tables for each MCMC iteration
#'         BBE array of regression coefficients for each MCMC iteration
#'         PPE array of transition probabilities for each MCMC iteration
#'         acctab overall acceptance rate for tables
#'         accbe overal acceptance rate for regression parameters
#'         Accbe acceptance rate for each regression parameter vector
#'         
#'         

dirmultAR_mcmc <- function(Y, R=10000, burnin=1000, tin=10, Tau=matrix(0.1,7,7),
                           mra=10, si2BE = 100, degree=2, tint=NULL, disp=FALSE,
                           ORmin=matrix(NA,7,7), ORmax=matrix(NA,7,7)){
  

Y <- as.matrix(Y)
TT <- nrow(Y) # total number of time instants
categories <- colnames(Y)
ldnorm1 <- function(x1,x2,si2=1) -sum(x1^2-x2^2)/(2*si2)



# initial tables
TAB <- array(0,c(7,7,TT))
TAB[,,1] <- NA
for(t in 2:TT){ 
  mi <- round(pmin(Y[t,],Y[t-1,])*c(1,0.9,0.9,0.9,0.9,0.9,1))
  rtot <- Y[t-1,]-mi
  ctot <- Y[t,]-mi
  TAB[,,t] <- diag(mi)+r2dtable(1,rtot,ctot)[[1]] # generate 2-way tables with given margins
}




# design matrices
if(is.null(tint)){
  XX<- cbind(1,bs(1:(TT-1),degree=degree))
}else{
  XX <- cbind(1,bs(1:(TT-1),degree=degree,knots = tint-1))
}
XXC <- XX



# initial values of the parameters and multinomial probabilities
mTAB = apply(TAB[,,2:TT],c(1,2),mean) # matrix of the element-wise mean of the matrices
nbe = degree+1+length(tint)
BE = array(0,c(7,7,nbe)) # regression parameters
Ind = rbind(cbind(1,1:7),cbind(2,2:7),cbind(3,2:7),cbind(4,2:7),cbind(5,2:7),cbind(6,2:7),c(7,7))
BE[2:7,1,] = BE[7,2:6,] = NA
for(h in 1:38){ # iter over admissible couples
  i = Ind[h,1]; j = Ind[h,2]
  if(j!=i){
    # if(!is.null(mBE)) BE[i,j,] = mBE[i,j,]
    if(!is.na(ORmin[i,j]) & is.na(ORmax[i,j])) BE[i,j,1] = max(BE[i,j,1],log(ORmin[i,j]))
    if(is.na(ORmin[i,j]) & !is.na(ORmax[i,j])) BE[i,j,1] = min(BE[i,j,1],log(ORmax[i,j]))
    if(!is.na(ORmin[i,j]) & !is.na(ORmax[i,j])){
      if(BE[i,j,1]<log(ORmin[i,j]) | BE[i,j,1]>log(ORmax[i,j])) BE[i,j,1] = log((ORmin[i,j]+ORmax[i,j])/2)
    }
  }
}
LA = PP = array(0,c(7,7,TT))
LA[,,1] = PP[,,1] = NA
for(h in 1:38) LA[Ind[h,1],Ind[h,2],2:TT] = exp(XX%*%BE[Ind[h,1],Ind[h,2],]) # compute num transision probs
for(t in 2:TT) PP[,,t] = (1/rowSums(LA[,,t]))*LA[,,t] # compute transition probs

# LAC = array(0,c(7,7,TT))
# LAC[,,1] = NA
# for(h in 1:38) LAC[Ind[h,1],Ind[h,2],2:TT] = exp(XXC%*%BE[Ind[h,1],Ind[h,2],])
OR = LA

### SIMULATION ALGORITHM -------------------------------------------------------------------

  acctab = accbe = accnu = 0
  Accbe = matrix(0,7,7)
  Accbe[2:6,1] = NA
  TTAB = array(as.integer(0),c(7,7,TT,R/tin))
  TTAB[,,,1] = NA
  BBE = array(0,c(7,7,nbe,R/tin))
  BBE[2:7,1,,] = BBE[7,2:6,,] = NA
  PPP = array(0,c(7,7,TT,R/tin))
  PPP[,,,1] = NA
  it = 0
  t0 = proc.time()[3]; names(t0) = NULL
  seqra = c(-mra:-1,1:mra)
  
  pb <- progressBar(min = -(burnin-1), max = R, style = "ETA")
  
  for(r in -(burnin-1):R){
  it = it+1
  
  #### STEP 1: update tables
  for(t in 2:TT){
    La = LA[,,t]; P = PP[,,t]
    for(it1 in 1){
      tmp = sample(1:6,2); i1 = min(tmp); i2 = max(tmp)
      tmp = sample(2:7,2); j1 = min(tmp); j2 = max(tmp)
      ra = sample(seqra,1)
      Tab = Tabs = TAB[,,t] # Tabs : current proposed table / Tab : previous table
      Tabs[i1,j1] = Tab[i1,j1]+ra
      Tabs[i1,j2] = Tab[i1,j2]-ra
      Tabs[i2,j1] = Tab[i2,j1]-ra
      Tabs[i2,j2] = Tab[i2,j2]+ra
      # if(t==TT) browser()
      if(all(Tabs>=0)){
        lnum = sum(lgamma(Tabs[i1,c(j1,j2)]+La[i1,c(j1,j2)])-lgamma(Tabs[i1,c(j1,j2)]+1)+
                     lgamma(Tabs[i2,c(j1,j2)]+La[i2,c(j1,j2)])-lgamma(Tabs[i2,c(j1,j2)]+1))
        lden = sum(lgamma(Tab[i1,c(j1,j2)]+La[i1,c(j1,j2)])-lgamma(Tab[i1,c(j1,j2)]+1)+
                     lgamma(Tab[i2,c(j1,j2)]+La[i2,c(j1,j2)])-lgamma(Tab[i2,c(j1,j2)]+1))
        al = exp(lnum-lden)
        # al1 = exp(lnum-lden)
        # print(c(al,al1,al/al1-1))
        if(is.nan(al)){
          print("NA in updating tables")
          browser()
        }
        if(runif(1)<=al){
          TAB[,,t] = Tabs
          acctab = acctab + 1/(1*(TT-1))
        }
      }
    }
  }
  
  #### STEP 2: update beta's
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
      # if(conbe) if(BES[i,j,3]>0) check=FALSE
      for(j1 in ind) if(j1!=i){
        if(!is.na(ORmin[i,j1])) if(any(ORS[i,j1,2:TT]<ORmin[i,j1])) check = FALSE
        if(!is.na(ORmax[i,j1])) if(any(ORS[i,j1,2:TT]>ORmax[i,j1])) check = FALSE
      }
      if(check){ 
        tmp = ldnorm1(BES[i,j,],BE[i,j,],si2BE)
        
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
  
  
  #### STEP 4: store values
  if(r>0 & r%%tin==0){
    r1 = round(r/tin)
    TTAB[,,,r1] = TAB
    BBE[,,,r1] = BE
    PPP[,,,r1] = PP
  }
  
  #### STEP 2: display output
  if(r%%250==0){
    if(disp){
      tt = proc.time()[3]; names(tt) = NULL
      print(c(iteration=r,acc_tab=acctab/it,acc_be=accbe/it,time100=100*(tt-t0)/it))
      print(Accbe)
      cat("\n")
      print("last matrix of OR")
      print(round(OR[,,TT],6))
      cat("\n")
      print("last table")
      print(cbind(TAB[,,TT],NA,round(Y[TT-1,]*PP[,,TT])))
      cat("\n")
    }
  }
  
  # ProgressBar
  setTxtProgressBar(pb, r)
  }
  
  
### OUTPUT ----------------------------------------------------------------------------------
  out = list(TTAB = TTAB, BBE = BBE, PPP = PPP,
             acctab = acctab/(R+burnin),accbe = accbe/(R+burnin), Accbe = Accbe/(R+burnin),
             Y=Y,call = match.call())
  
  class(out) = "dirmultAR"
  return(out)
}
