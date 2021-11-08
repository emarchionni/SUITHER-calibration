remove(list = ls())

library(coda)

setwd('C:/Users/edoar/Desktop/NAPDE/Project/Code/SUITHER-calibration/Posterior inference')

load('C:/Users/edoar/Desktop/NAPDE/Project/Code/model_7.RData')

Ind <- rbind(cbind(1,c(1,6)), cbind(2,2), cbind(3,c(2:4,7)),
            cbind(4,c(2,4,5)), cbind(5,c(4,5,7)),
            cbind(6,c(2,3,6)), c(7,7))


Ind_const <- rbind(cbind(3,c(2,7)), c(4,2), c(5,4), cbind(6,c(2,3)))
names_const <- c("rho_I",
                  "gamma_I",
                  "rho_H",
                  "theta_T",
                  "rho_U",
                  "delta")
rownames(Ind_const) <- names_const

Ind_vary <- rbind(c(1,6), c(3,4), c(4,5), c(5,7))
names_vary <- c("beta_U",
                "omega_I",
                "omega_H",
                "gamma_T")
rownames(Ind_vary) <- names_vary


dummies <- c(41, 54, 72, 83, 92, 97, 107, 114, 126)


attach(out_no_ORmax)

TT <- nrow(Y)


PP <- list()

### constant parameters ####

for (h in 1:nrow(Ind_const)) {
  
  
  i = Ind_const[h,1]; j = Ind_const[h,2]
  
  ## posterior mean (Bayes' estimate quadratic loss function)
  m <- apply(PPP[i,j,,], 1 , mean)
 

  ## mean
  PP[["const"]][[names_const[h]]][["mean"]] <- mean(m, na.rm = TRUE)
  
  ## variance
  PP[["const"]][[names_const[h]]][["variance"]] <- var(m, na.rm = TRUE)
  
  
  estimate <- cbind(PP[["const"]][[names_const[h]]][["mean"]], PP[["const"]][[names_const[h]]][["variance"]])
  
  write.table(estimate, paste0("constant/",names_const[h],".txt"),row.names = FALSE, col.names = c("mean", "variance"))
  
  
  print(h)
  
}


#### varying parameters ####


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
  
  write.table(estimate, paste0("varying/",names_vary[h],".txt"),row.names = FALSE, col.names = c("mean", "variance"))
  
  
}

# mm <- c()
# for (temp in 1:40) {
#   d <- mean(PPP[i,j,temp,])
#   mm <- rbind(mm, d)
# 
# }

# write.table(adj,paste0("UserLine/",nameline,"/",nameline,"_", as.character(giorno),".txt"),row.names = FALSE, col.names = TRUE)


# View(PP)




