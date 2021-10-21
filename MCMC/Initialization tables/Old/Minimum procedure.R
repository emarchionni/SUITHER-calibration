TT <- nrow(Y)

# initial tables
TAB <- array(0,c(7,7,TT))
TAB[,,1] <- NA
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
  
  d <- diag(TAB[,,t])
  TAB[,,t] <- TAB[,,t] - diag(d)
  
  ind <- rbind(cbind(3,c(2,3,4,7)),cbind(4,c(2,4,5)),cbind(5,c(4,5,7)),cbind(6,c(2,3)))
  
  for (count in 1:12) {
    i <- ind[count,1]; j <- ind[count,2]
    TAB[i,j,t] <- min(rtot[i],ctot[j])
    ctot[j] <- ctot[j] - TAB[i,j,t]; rtot[i] <- rtot[i] - TAB[i,j,t]
  }
  
  TAB[,,t] <- TAB[,,t] + diag(d)
  
}