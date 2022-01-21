remove(list = ls())

library(coda)

setwd('C:/Users/edoar/Desktop/NAPDE/Project/Code/SUITHER-calibration/Posterior inference')

load('C:/Users/edoar/Desktop/NAPDE/Project/Code/model_7.RData')

Ind <- rbind(cbind(1,c(1,6)), cbind(2,2), cbind(3,c(2:4,7)),
            cbind(4,c(2,4,5)), cbind(5,c(4,5,7)),
            cbind(6,c(2,3,6)), c(7,7))

categories <- c('susceptible',
                'recovered',
                'isolated',
                'hospitalized',
                'threatened',
                'undetected',
                'extinct')

#### POSTERIOR ESTIMATE p OVER THE DIFFERENT PHASES ####

Ind_const <- rbind(cbind(3,c(2,7)), c(4,2), c(5,4), cbind(6,c(2,3)))
names_const <- c("p_3_2",
                 "p_3_7",
                 "p_4_2",
                 "p_5_4",
                 "p_6_2",
                 "p_6_3")
rownames(Ind_const) <- names_const

Ind_vary <- rbind(c(1,6), c(3,4), c(4,5), c(5,7))
names_vary <- c("p_1_6",
                "p_3_4",
                "p_4_5",
                "p_5_7")
rownames(Ind_vary) <- names_vary


dummies <- c(41, 54, 72, 83, 92, 97, 107, 114, 126)


attach(out_no_ORmax)

TT <- nrow(Y)


PP <- list()

### constant parameters ###

for (h in 1:nrow(Ind_const)) {
  
  
  i = Ind_const[h,1]; j = Ind_const[h,2]
  
  ## posterior mean (Bayes' estimate quadratic loss function)
  m <- apply(PPP[i,j,,], 1 , mean)
 

  ## mean
  PP[["const"]][[names_const[h]]][["mean"]] <- mean(m, na.rm = TRUE)
  
  ## variance
  PP[["const"]][[names_const[h]]][["variance"]] <- var(m, na.rm = TRUE)
  
  
  estimate <- cbind(PP[["const"]][[names_const[h]]][["mean"]], PP[["const"]][[names_const[h]]][["variance"]])
  
  write.table(estimate, paste0("p_phases/constant/",names_const[h],".txt"),row.names = FALSE, col.names = c("mean", "variance"))
  
  
  print(h)
  
}


#### varying parameters ###


for (h in 1:nrow(Ind_vary)) {
  
  
  i = Ind_vary[h,1]; j = Ind_vary[h,2]
  
  estimate <- c()
  tcnt <- 0
  tI <- 0
  for (t in dummies) {
    
    tcnt <- tcnt + 1
    
    ## posterior mean (Bayes' estimate quadratic loss function)
    m <- apply(PPP[i,j,tI:(t-1),], 1 , mean)
    
    
    
    ## mean
    PP[["varying"]][[names_vary[h]]][["mean"]][[tcnt]] <- mean(m, na.rm = TRUE)
    
    ## variance
    PP[["varying"]][[names_vary[h]]][["variance"]][[tcnt]] <- var(m, na.rm = TRUE)
    
    
    estimate <- rbind(estimate,
                      cbind(PP[["varying"]][[names_vary[h]]][["mean"]][[tcnt]], 
                            PP[["varying"]][[names_vary[h]]][["variance"]][[tcnt]]))
    
    
    tI <- t
    
    print(c(h,t))
    
  }
  
  write.table(estimate, paste0("p_phases/varying/",names_vary[h],".txt"),row.names = FALSE, col.names = c("mean", "variance"))
  
  
}


#### POSTERIOR ESTIMATE p OVER TIME INSTANTS ####

ntimes <-  dim(Y)[1]
PPP_mean <- array(0,c(7,7,134))

for(h in 1:nrow(Ind)){
  
  i = Ind[h,1]; j = Ind[h,2]
  
  ## posterior mean (Bayes' estimate quadratic loss function)
  PPP_mean[i,j,] <- apply(PPP[i,j,,], 1, mean)
  
}

for(t in 2:TT){
  write.table(PPP_mean[,,t], paste0("p_time_instants/",t,".txt"),row.names = categories, col.names = categories)
}


#### POSTERIOR ESTIMATE CONTINGENCY TABLES ####


mTAB <- apply(TTAB,1:3,mean)
lwTAB <- apply(TTAB,1:3,quantile,probs=0.025,na.rm=TRUE)
upTAB <- apply(TTAB,1:3,quantile,probs=0.975,na.rm=TRUE)

for(t in 2:TT){
  write.table(mTAB[,,t], paste0("tab/mean/",t,".txt"),row.names = categories, col.names = categories)
  write.table(lwTAB[,,t], paste0("tab/lower_bound/",t,".txt"),row.names = categories, col.names = categories)
  write.table(upTAB[,,t], paste0("tab/upper_bound/",t,".txt"),row.names = categories, col.names = categories)
}

View(mTAB[,,40])
matplot(mTAB[2,2,], type='l')
View(mTAB[2,2,])








###############################################

TAB_flux <- array(0,c(7,7,134))
Y_new <- Y

for(t in 2:TT){
  
  for(i in 1:7){
    
    for(j in 1:7){
      
      TAB_flux[i,j,t] <- floor(Y_new[t-1,i] * PPP_mean[i,j,t])
      
    }
    
    missing <- Y_new[t-1,i]  - sum(TAB_flux[i,,t])
    if(missing){
      
      if(missing<0){ warning('Value less than 0'); print(t,i,j)}
      
      if(i == 1){
        
        if(TAB_flux[i,1,t] == 0)
          col <- 1
        else if(TAB_flux[i,6,t] == 0)
          col <- 6
        else col <- sample(c(1,6), size = 1)
        
      } else if(i == 2){
        
        col <- 2
        
      } else if(i == 3){
        
        if(TAB_flux[i,2,t] == 0)
          col <- 2
        else if(TAB_flux[i,3,t] == 0)
          col <- 3
        else if(TAB_flux[i,4,t] == 0)
          col <- 4
        else if(TAB_flux[i,7,t] == 0)
          col <- 7
        else col <- sample(c(2,3,4,7), size = 1)
        
      } else if(i == 4){
        
        if(TAB_flux[i,2,t] == 0)
          col <- 2
        else if(TAB_flux[i,4,t] == 0)
          col <- 4
        else if(TAB_flux[i,5,t] == 0)
          col <- 5
        else col <- sample(c(2,4,5), size = 1)
        
      } else if(i == 5){
        
        if(TAB_flux[i,4,t] == 0)
          col <- 4
        else if(TAB_flux[i,5,t] == 0)
          col <- 5
        else if(TAB_flux[i,7,t] == 0)
          col <- 7
        else col <- sample(c(4,5,7), size = 1)
        
      } else if(i == 6){
        
        if(TAB_flux[i,2,t] == 0)
          col <- 2
        else if(TAB_flux[i,3,t] == 0)
          col <- 3
        else if(TAB_flux[i,6,t] == 0)
          col <- 6
        else col <- sample(c(2,3,6), size = 1)
        
      } else{
        
        col <- 7
        
      }
        
      
      TAB_flux[i,col,t] <- TAB_flux[i,col,t] + missing
    }
    
    
  }
  
  Y_new[t, ] <- colSums(TAB_flux[,,t])
  
}





# mm <- c()
# for (temp in 1:40) {
#   d <- mean(PPP[i,j,temp,])
#   mm <- rbind(mm, d)
# 
# }

# write.table(adj,paste0("UserLine/",nameline,"/",nameline,"_", as.character(giorno),".txt"),row.names = FALSE, col.names = TRUE)


# View(PP)




