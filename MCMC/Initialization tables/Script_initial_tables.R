library(stats)

setwd('C:/Users/edoar/Desktop/NAPDE/Project/Code/SUITHER-calibration/MCMC/Initialization tables')

remove(list = ls()[1:15])
Y <- Y[,-1]
col <- c(2,3,4,5,7)
perm <- gtools::permutations(n = length(col), r = length(col), v = col)

source('Functions.R') # load function

TT <- nrow(Y)

# initial tables
TAB <- array(0,c(7,7,TT))
TAB[,,1] <- NA
verify <- vector(mode = 'logical', length = TT)

for (t in 2:TT) {
  
  ### Set constrained values in contingency tables
  TAB[1,1,t] <- Y[t, 1]
  TAB[1,6,t] <- Y[t-1, 1] - TAB[1,1,t]
  TAB[6,6,t] <- Y[t, 6] - TAB[1,6,t]
  TAB[2,2,t] <- Y[t-1, 2]
  TAB[7,7,t] <- Y[t-1, 7]
  
  # new margins after subtraction already set values
  row_margin <- Y[t-1,] - apply(TAB[,,t], 1, sum) 
  col_margin <- Y[t,] - apply(TAB[,,t], 2, sum) 
  
  ### Generate new tables
  mi <- round(pmin(Y[t,],Y[t-1,])*c(0,0,0.9,0.9,0.9,0,0))
  rtot <- row_margin - mi
  ctot <- col_margin - mi
  TAB[,,t] <- diag(mi) + TAB[,,t]
  
  
  good_permutation <- FALSE
  totperm <- nrow(perm)
  
  for (rperm in 1:totperm) {
    
    if(!good_permutation) { # look for a permutation if not already found one
      
      tab <- check_permutation(TAB = TAB[,,t],
                               ctot = ctot,
                               rtot = rtot,
                               perm = perm[rperm,])
      
      if(!is.logical(tab)){
        
        TAB[,,t] <- tab
        # print(paste('Time',t,'ok'))
        good_permutation <- TRUE
        if(all(apply(TAB[,,t],1,sum) == Y[t-1,]) & all(apply(TAB[,,t],2,sum) == Y[t,])) verify[t] <- T;
        
      }
      
      if(rperm == totperm & !good_permutation) print(paste('Table at time', t, 'not initialized'));
      
    }
    
  }
  
}

save(TAB, file = 'initial_tables.RData')


