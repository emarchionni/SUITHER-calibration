
check_permutation <- function(TAB, ctot, rtot, perm){
  
  
  out <- list(TAB = TAB,
              ctot = ctot,
              rtot = rtot)
  
  # check if each column can be filled in the order of permutation contained in perm
  for (j in perm) 
  {
    
    out <- col_filling(TAB = out$TAB,
                       ctot = as.array(unlist(out$ctot)),
                       rtot = as.array(unlist(out$rtot)),
                       j = j)
    
    if(is.logical(out)) return(FALSE);
    
    TAB[,j] <- out$TAB[,j]
  }
  
  return(TAB);
  
}




col_filling <- function(TAB, ctot, rtot, j){
  #'@param  j number of column
  
  ctot <- as.array(ctot)
  rtot <- as.array(rtot)
  TAB <- as.matrix(TAB)

  if(j == 2){
   
    # fill column 2 (rows: 3, 4, 6)
    if(rtot[3] > ctot[2] & rtot[4] > ctot[2] & rtot[6] > ctot[2]){
      
      TAB[3,2] <- sample(1:(ctot[2]-1), size = 1); ctot[2] <- ctot[2] - TAB[3,2]
      TAB[4,2] <- sample(1:ctot[2], size = 1); ctot[2] <- ctot[2] - TAB[4,2]
      TAB[6,2] <- ctot[2]; ctot[2] <- ctot[2] - TAB[6,2]
      
    } else if(rtot[3] < ctot[2] & rtot[4] > ctot[2] & rtot[6] > ctot[2]){
      
      TAB[3,2] <- rtot[3]; ctot[2] <- ctot[2] - TAB[3,2]
      TAB[4,2] <- sample(1:ctot[2], size = 1); ctot[2] <- ctot[2] - TAB[4,2]
      TAB[6,2] <- ctot[2]; ctot[2] <- ctot[2] - TAB[6,2]
      
    } else if(rtot[3] > ctot[2] & rtot[4] < ctot[2] & rtot[6] > ctot[2]){
      
      TAB[4,2] <- rtot[4]; ctot[2] <- ctot[2] - TAB[4,2]
      TAB[3,2] <- sample(1:ctot[2], size = 1); ctot[2] <- ctot[2] - TAB[3,2]
      TAB[6,2] <- ctot[2]; ctot[2] <- ctot[2] - TAB[6,2]
      
    } else if(rtot[3] > ctot[2] & rtot[4] > ctot[2] & rtot[6] < ctot[2]){
      
      TAB[6,2] <- rtot[6]; ctot[2] <- ctot[2] - TAB[6,2]
      TAB[4,2] <- sample(1:ctot[2], size = 1); ctot[2] <- ctot[2] - TAB[4,2]
      TAB[3,2] <- ctot[2]; ctot[2] <- ctot[2] - TAB[3,2]
      
    } else {
      
      if(rtot[3] + rtot[4] + rtot[6] >= ctot[2]){
        
        new_rtot <- c(0,0,rtot[3],rtot[4],0,rtot[6],0)
        temp <- sort(c(rtot[3],rtot[4],rtot[6]))
        
        i <- which(new_rtot == temp[1])
        TAB[i,2] <- min(rtot[i],ctot[2]); ctot[2] <- ctot[2] - TAB[i,2]
        i <- which(new_rtot == temp[2])
        TAB[i,2] <- min(rtot[i],ctot[2]); ctot[2] <- ctot[2] - TAB[i,2]
        i <- which(new_rtot == temp[3])
        TAB[i,2] <- min(rtot[i],ctot[2]); ctot[2] <- ctot[2] - TAB[i,2]
        
      } else {
        return(FALSE)
      }
    }
    
    rtot[3] <- rtot[3] - TAB[3,2]
    rtot[4] <- rtot[4] - TAB[4,2]
    rtot[6] <- rtot[6] - TAB[6,2]
    
  } else if (j == 3){
    
    # fill column 3 (rows: 3, 6)
    d <- TAB[3,3]
    if(rtot[6] > ctot[3] & rtot[3] > ctot[3]){
      
      TAB[6,3] <- sample(1:ctot[3], size = 1); ctot[3] <- ctot[5] - TAB[6,3]
      TAB[3,3] <- TAB[3,3] + ctot[3]; ctot[3] <- ctot[3] - TAB[3,3] + d
      
    } else if(rtot[6] > ctot[3] & rtot[3] < ctot[3]){
      
      TAB[3,3] <- TAB[3,3] + rtot[3]; ctot[3] <- ctot[3] - TAB[3,3] + d
      TAB[6,3] <- ctot[3]; ctot[3] <- ctot[3] - TAB[6,3]
      
    } else if(rtot[6] < ctot[3] & rtot[3] > ctot[3]){
      
      TAB[6,3] <- rtot[6]; ctot[3] <- ctot[3] - TAB[6,3]
      TAB[3,3] <- TAB[3,3] + ctot[3]; ctot[3] <- ctot[3] - TAB[3,3] + d

    } else {
      
      if(rtot[6] + rtot[3] >= ctot[3]){ 
        
        if(rtot[6] > rtot[3]){
          TAB[6,3] <- rtot[6]
          TAB[3,3] <- TAB[3,3] + ctot[3] - TAB[6,3]
        } else {
          TAB[3,3] <- TAB[3,3] + rtot[3]
          TAB[6,3] <- ctot[3] - TAB[3,3] + d
        }
        ctot[3] <- ctot[3] - TAB[3,3] - TAB[6,3] + d
        
      } else {
        return(FALSE)
      }
    }
    
    rtot[6] <- rtot[6] - TAB[6,3]
    rtot[3] <- rtot[3] - TAB[3,3] + d
    
  } else if (j == 4){
    
    # fill column 4 (rows: 3,4,5)
    d <- TAB[4,4]
    if(rtot[3] > ctot[4] & rtot[4] > ctot[4] & rtot[5] > ctot[4]){
      
      TAB[3,4] <- sample(1:(ctot[4]-1), size = 1); ctot[4] <- ctot[4] - TAB[3,4]
      TAB[5,4] <- sample(1:ctot[4], size = 1); ctot[4] <- ctot[4] - TAB[5,4]
      TAB[4,4] <- TAB[4,4] + ctot[4]; ctot[4] <- ctot[4] - TAB[4,4] + d
      
    } else if(rtot[3] < ctot[4] & rtot[4] > ctot[4] & rtot[5] > ctot[4]){
      
      TAB[3,4] <- rtot[3]; ctot[4] <- ctot[4] - TAB[3,4]
      TAB[5,4] <- sample(1:ctot[4], size = 1); ctot[4] <- ctot[4] - TAB[5,4]
      TAB[4,4] <- TAB[4,4] + ctot[4]; ctot[4] <- ctot[4] - TAB[4,4] - d
      
    } else if(rtot[3] > ctot[4] & rtot[4] < ctot[4] & rtot[5] > ctot[4]){
      
      TAB[4,4] <- TAB[4,4] + rtot[4]; ctot[4] <- ctot[4] - TAB[4,4] + d
      TAB[3,4] <- sample(1:ctot[4], size = 1); ctot[4] <- ctot[4] - TAB[3,4]
      TAB[5,4] <- ctot[4]; ctot[4] <- ctot[4] - TAB[5,4]
      
    } else if(rtot[3] > ctot[4] & rtot[4] > ctot[4] & rtot[5] < ctot[4]){
      
      TAB[5,4] <- rtot[5]; ctot[4] <- ctot[4] - TAB[5,4]
      TAB[3,4] <- sample(1:ctot[4], size = 1); ctot[4] <- ctot[4] - TAB[3,4]
      TAB[4,4] <- ctot[4]; ctot[4] <- ctot[4] - TAB[4,4] + d
      
    } else {
      if(rtot[3] + rtot[4] + rtot[5] >= ctot[4]){
        
        TAB[4,4] <- 0
        
        temp <- sort(c(rtot[3],rtot[4],rtot[5]))
        
        i <- which(rtot == temp[1])
        TAB[i,4] <- min(rtot[i],ctot[4]); ctot[4] <- ctot[4] - TAB[i,4]
        i <- which(rtot == temp[2])
        TAB[i,4] <- min(rtot[i],ctot[4]); ctot[4] <- ctot[4] - TAB[i,4]
        i <- which(rtot == temp[3])
        TAB[i,4] <- min(rtot[i],ctot[4]); ctot[4] <- ctot[4] - TAB[i,4]
    
        TAB[4,4] <- TAB[4,4] + d
        
      } else {
        return(FALSE)
      }
    }
    
    rtot[3] <- rtot[3] - TAB[3,4]
    rtot[4] <- rtot[4] - TAB[4,4] + d
    rtot[5] <- rtot[5] - TAB[5,4]
    
  } else if(j == 5){
    
    # fill column 5 (rows: 4,5)
    d <- TAB[5,5]
    if(rtot[4] > ctot[5] & rtot[5] > ctot[5]){
      
      TAB[4,5] <- sample(1:ctot[5], size = 1); ctot[5] <- ctot[5] - TAB[4,5]
      TAB[5,5] <- TAB[5,5] + ctot[5]
      ctot[5] <- ctot[5] - TAB[5,5] + d
      
    } else if(rtot[4] > ctot[5] & rtot[5] < ctot[5]){
      
      TAB[5,5] <- TAB[5,5] + rtot[5]; ctot[5] <- ctot[5] - TAB[5,5] + d 
      TAB[4,5] <- ctot[5]
      ctot[5] <- ctot[5] - TAB[4,5]
      
    } else if(rtot[4] < ctot[5] & rtot[5] > ctot[5]){
      
      TAB[4,5] <- rtot[4]; ctot[5] <- ctot[5] - TAB[4,5]
      TAB[5,5] <- TAB[5,5] + ctot[5]
      ctot[5] <- ctot[5] - TAB[5,5] + d
      
    } else {
      
      if(rtot[4] + rtot[5] >= ctot[5]){ 
        if(rtot[4] > rtot[5]){ 
          TAB[4,5] <- rtot[4]
          TAB[5,5] <- TAB[5,5] + ctot[5] - TAB[4,5]
          print(4)
        } else {
          TAB[5,5] <- TAB[5,5] + rtot[5]
          TAB[4,5] <- ctot[5] - TAB[5,5] + d
          print(5)
        }
        ctot[5] <- ctot[5] - TAB[5,5] - TAB[4,5] + d
        
      } else {
        return(FALSE)
      }
    }
    
    rtot[4] <- rtot[4] - TAB[4,5]
    rtot[5] <- rtot[5] - TAB[5,5] + d
    
  } else if (j == 7){
    
    # fill column 7 (rows: 3,5)
    if(rtot[3] > ctot[7] & rtot[5] > ctot[7]){
      
      TAB[3,7] <- sample(1:ctot[7], size = 1); ctot[7] <- ctot[7] - TAB[3,7]
      TAB[5,7] <- ctot[7]
      ctot[7] <- ctot[7] - TAB[5,7]
      
    } else if(rtot[3] > ctot[7] & rtot[5] < ctot[7]){
      
      TAB[5,7] <- rtot[5]; ctot[7] <- ctot[7] - TAB[5,7]
      TAB[3,7] <- ctot[7]
      ctot[7] <- ctot[7] - TAB[3,7]
      
    } else if(rtot[3] < ctot[7] & rtot[5] > ctot[7]){
      
      TAB[3,7] <- rtot[3]; ctot[7] <- ctot[7] - TAB[3,7]
      TAB[5,7] <- ctot[7]
      ctot[7] <- ctot[7] - TAB[5,7]
      
    } else {
      
      if(rtot[3] + rtot[5] >= ctot[7]){
        if(rtot[3] > rtot[5]){
          TAB[3,7] <- rtot[3]
          TAB[5,7] <- ctot[7] - TAB[3,7]
        } else {
          TAB[5,7] <- rtot[5]
          TAB[3,7] <- ctot[7] - TAB[5,7]
        }
        ctot[7] <- ctot[7] - TAB[3,7] - TAB[3,7]
        
      } else {
        return(FALSE)
      }
    }
    
    rtot[3] <- rtot[3] - TAB[3,7]
    rtot[5] <- rtot[5] - TAB[5,7]
  }
  
  out <- list(TAB = TAB,
              ctot = ctot,
              rtot = rtot)
  
  return(out)
}


