TAB <- array(0,c(7,7,TT))
TAB[,,1] <- NA
for(t in 2:TT){
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
  
  # fill coloumn 4 (rows: 3,4,5)
  d <- TAB[4,4,t]
  if(rtot[3] > ctot[4] & rtot[4] > ctot[4] & rtot[5] > ctot[4]){
    
    TAB[3,4,t] <- sample(1:(ctot[4]-1), size = 1); ctot[4] <- ctot[4] - TAB[3,4,t]
    TAB[5,4,t] <- sample(1:ctot[4], size = 1); ctot[4] <- ctot[4] - TAB[5,4,t]
    TAB[4,4,t] <- TAB[4,4,t] + ctot[4]; ctot[4] <- ctot[4] - TAB[4,4,t] + d
    
  } else if(rtot[3] < ctot[4] & rtot[4] > ctot[4] & rtot[5] > ctot[4]){
    
    TAB[3,4,t] <- rtot[3]; ctot[4] <- ctot[4] - TAB[3,4,t]
    TAB[5,4,t] <- sample(1:ctot[4], size = 1); ctot[4] <- ctot[4] - TAB[5,4,t]
    TAB[4,4,t] <- TAB[4,4,t] + ctot[4]; ctot[4] <- ctot[4] - TAB[4,4,t] - d
    
  } else if(rtot[3] > ctot[4] & rtot[4] < ctot[4] & rtot[5] > ctot[4]){
    
    TAB[4,4,t] <- TAB[4,4,t] + rtot[4]; ctot[4] <- ctot[4] - TAB[4,4,t] + d
    TAB[3,4,t] <- sample(1:(ctot[4]-1), size = 1); ctot[4] <- ctot[4] - TAB[3,4,t]
    TAB[5,4,t] <- ctot[4]; ctot[4] <- ctot[4] - TAB[5,4,t]
    
  } else if(rtot[3] > ctot[4] & rtot[4] > ctot[4] & rtot[5] < ctot[4]){
    
    TAB[5,4,t] <- rtot[5]; ctot[4] <- ctot[4] - TAB[5,4,t]
    TAB[3,4,t] <- sample(1:ctot[4], size = 1); ctot[4] <- ctot[4] - TAB[3,4,t]
    TAB[4,4,t] <- ctot[4]; ctot[4] <- ctot[4] - TAB[4,4,t] + d
    
  } else {
    if(rtot[3] + rtot[4] + rtot[5] >= ctot[4]){
      
      TAB[4,4,t] <- 0
      
      temp <- sort(c(rtot[3],rtot[4],rtot[5]))
      
      i <- which(rtot == temp[1])
      TAB[i,4,t] <- min(rtot[i],ctot[4]); ctot[4] <- ctot[4] - TAB[i,4,t]
      i <- which(rtot == temp[2])
      TAB[i,4,t] <- min(rtot[i],ctot[4]); ctot[4] <- ctot[4] - TAB[i,4,t]
      i <- which(rtot == temp[3])
      TAB[i,4,t] <- min(rtot[i],ctot[4]); ctot[4] <- ctot[4] - TAB[i,4,t]
      
      TAB[4,4,t] <- TAB[4,4,t] + d
    } else {
      stop('Error filling column 4')
    }
  }
  rtot[3] <- rtot[3] - TAB[3,4,t]
  rtot[4] <- rtot[4] - TAB[4,4,t] + d
  rtot[5] <- rtot[5] - TAB[5,4,t]
  
  # fill coloumn 7 (rows: 3,5)
  if(rtot[3] > ctot[7] & rtot[5] > ctot[7]){
    TAB[3,7,t] <- sample(1:ctot[7], size = 1); ctot[7] <- ctot[7] - TAB[3,7,t]
    TAB[5,7,t] <- ctot[7]
    ctot[7] <- ctot[7] - TAB[5,7,t]
  } else if(rtot[3] > ctot[7] & rtot[5] < ctot[7]){
    TAB[5,7,t] <- rtot[5]; ctot[7] <- ctot[7] - TAB[5,7,t]
    TAB[3,7,t] <- ctot[7]
    ctot[7] <- ctot[7] - TAB[3,7,t]
  } else if(rtot[3] < ctot[7] & rtot[5] > ctot[7]){
    TAB[3,7,t] <- rtot[3]; ctot[7] <- ctot[7] - TAB[3,7,t]
    TAB[5,7,t] <- ctot[7]
    ctot[7] <- ctot[7] - TAB[5,7,t]
  } else {
    if(rtot[3] + rtot[5] >= ctot[7]){
      if(rtot[3] > rtot[5]){
        TAB[3,7,t] <- rtot[3]
        TAB[5,7,t] <- ctot[7] - TAB[3,7,t]
      } else {
        TAB[5,7,t] <- rtot[5]
        TAB[3,7,t] <- ctot[7] - TAB[5,7,t]
      }
      ctot[7] <- ctot[7] - TAB[3,7,t] - TAB[3,7,t]
    } else {
      stop('Error filling column 7')
    }
  }
  rtot[3] <- rtot[3] - TAB[3,7,t]
  rtot[5] <- rtot[5] - TAB[5,7,t]
  
  # fill coloumn 5 (rows: 4,5)
  d <- TAB[5,5,t]
  if(rtot[4] > ctot[5] & rtot[5] > ctot[5]){
    TAB[4,5,t] <- sample(1:ctot[5], size = 1); ctot[5] <- ctot[5] - TAB[4,5,t]
    TAB[5,5,t] <- TAB[5,5,t] + ctot[5]
    ctot[5] <- ctot[5] - TAB[5,5,t] + d
    print(1)
  } else if(rtot[4] > ctot[5] & rtot[5] < ctot[5]){
    TAB[5,5,t] <- TAB[5,5,t] + rtot[5]; ctot[5] <- ctot[5] - TAB[5,5,t] + d 
    TAB[4,5,t] <- ctot[5]
    ctot[5] <- ctot[5] - TAB[4,5,t]
    print(2)
  } else if(rtot[4] < ctot[5] & rtot[5] > ctot[5]){
    TAB[4,5,t] <- rtot[4]; ctot[5] <- ctot[5] - TAB[4,5,t]
    TAB[5,5,t] <- TAB[5,5,t] + ctot[5]
    ctot[5] <- ctot[5] - TAB[5,5,t] + d
    print(3)
  } else {
    if(rtot[4] + rtot[5] >= ctot[5]){ 
      if(rtot[4] > rtot[5]){ 
        TAB[4,5,t] <- rtot[4]
        TAB[5,5,t] <- TAB[5,5,t] + ctot[5] - TAB[4,5,t]
        print(4)
      } else {
        TAB[5,5,t] <- TAB[5,5,t] + rtot[5]
        TAB[4,5,t] <- ctot[5] - TAB[5,5,t] + d
        print(5)
      }
      ctot[5] <- ctot[5] - TAB[5,5,t] - TAB[4,5,t] + d
    } else {
      stop('Error filling column 5')
    }
  }
  rtot[4] <- rtot[4] - TAB[4,5,t]
  rtot[5] <- rtot[5] - TAB[5,5,t] + d
  
  
  # fill coloumn 3 (rows: 3, 6)
  d <- TAB[3,3,t]
  if(rtot[6] > ctot[3] & rtot[3] > ctot[3]){
    TAB[6,3,t] <- sample(1:ctot[3], size = 1); ctot[3] <- ctot[5] - TAB[6,3,t]
    TAB[3,3,t] <- TAB[3,3,t] + ctot[3]; ctot[3] <- ctot[3] - TAB[3,3,t] + d
    print(1)
  } else if(rtot[6] > ctot[3] & rtot[3] < ctot[3]){
    TAB[3,3,t] <- TAB[3,3,t] + rtot[3]; ctot[3] <- ctot[3] - TAB[3,3,t] + d
    TAB[6,3,t] <- ctot[3]; ctot[3] <- ctot[3] - TAB[6,3,t]
    print(2)
  } else if(rtot[6] < ctot[3] & rtot[3] > ctot[3]){
    TAB[6,3,t] <- rtot[6]; ctot[3] <- ctot[3] - TAB[6,3,t]
    TAB[3,3,t] <- TAB[3,3,t] + ctot[3]; ctot[3] <- ctot[3] - TAB[3,3,t] + d
    print(3)
  } else {
    if(rtot[6] + rtot[3] >= ctot[3]){ 
      if(rtot[6] > rtot[3]){
        TAB[6,3,t] <- rtot[6]
        TAB[3,3,t] <- TAB[3,3,t] + ctot[3] - TAB[6,3,t]
        print(4)
      } else {
        TAB[3,3,t] <- TAB[3,3,t] + rtot[3]
        TAB[6,3,t] <- ctot[3] - TAB[3,3,t] + d
        print(5)
      }
      ctot[3] <- ctot[3] - TAB[3,3,t] - TAB[6,3,t] + d
    } else {
      stop('Error filling column 3')
    }
  }
  rtot[6] <- rtot[6] - TAB[6,3,t]
  rtot[3] <- rtot[3] - TAB[3,3,t] + d
  
  
  
  
  
  
  
  # # subtable rows: 3,6; coloumns: 2,3
  # ctotsub <- ctot[2:3]; rtotsub <- c(sum(ctotsub)-rtot[6],rtot[6])
  # TAB[c(3,6),c(2,3),t] <- TAB[c(3,6),c(2,3),t] + r2dtable(1,rtotsub,ctotsub)[[1]]
  # ctot[c(2,3)] <- ctot[c(2,3)] - ctotsub
  # rtot[c(3,6)] <- rtot[c(3,6)] - rtotsub
  # 
  # # subtable rows: 3,4; coloums: 2,4
  # ctotsub <- ctot[c(2,4)]; rtotsub <- c(sum(ctotsub)-rtot[6],rtot[6])
  # 
  # 
  # # row 3
  # maxel2 <- min(rtot[3],col[2])
  # TAB[3,2,t] <- sample(1:(rtot[3]-1), size = 1); rtot[3] <- rtot[3] - TAB[3,2,t] 
  # TAB[3,4,t] <- sample(1:(rtot[3]), size = 1); rtot[3] <- rtot[3] - TAB[3,4,t]   
  # TAB[3,7,t] <- rtot[3]
  # 
  # ctot[2] <- ctot[2] - TAB[3,2,t]
  # ctot[4] <- ctot[4] - TAB[3,4,t]
  # ctot[7] <- ctot[7] - TAB[3,7,t]
  # 
  # # row 4 (coerced)
  # TAB[4,5,t] <- ctot[5]
  # TAB[4,2,t] <- rtot[4] - TAB[4,5,t]
  # 
  # 
  # # row 5 (coerced)
  # TAB[5,4,t] <- ctot[4]
  # TAB[5,7,t] <- ctot[7]
  # 
  # #row 6 (coerced)
  # TAB[6,2,t] <- ctot[2]
  # TAB[6,3,t] <- ctot[3]
  
  
  
}
