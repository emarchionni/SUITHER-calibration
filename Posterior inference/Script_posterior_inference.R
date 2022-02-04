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


attach(out_no_ORmax)
TT <- nrow(Y)


### saving plot of INFECTED
png(filename = 'Infected.png', height = 800, width = 800)
matplot(Y[,3]+Y[,4]+Y[,5]+Y[,6],
        type = 'l',
        lty = c(1),
        col = 'red',
        lwd = 2,
        xlab = 'time',
        ylab = 'frequencies',
        main = 'infected people')
for(time in dummies) abline(v = time, col='grey')
dev.off()


### saving plot of UNDETECTED
png(filename = 'Undetected.png', height = 800, width = 800)
matplot(Y[,6],
        type = 'l',
        lty = c(1),
        col = 'red',
        lwd = 2,
        xlab = 'time',
        ylab = 'frequencies',
        main = 'Undetected people')
for(time in dummies) abline(v = time, col='grey')
dev.off()

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



### plotting
for (h in 1:nrow(Ind)) {
  
  i = Ind[h,1]; j = Ind[h,2]
  
  x11()
  matplot(PPP_mean[i,j,],
          type = 'l',
          lty = c(1),
          col = 'orange',
          lwd = 2,
          xlab = 'time',
          ylab = 'transition probability',
          main = paste('Probability from', categories[i], 'to', categories[j]))
  for(time in dummies) abline(v = time, col='grey')
}

### saving plot form S to U
png(filename = 'Probability_fromS_toU.png', height = 800, width = 800)
matplot(PPP_mean[1,6,],
        type = 'l',
        lty = c(1),
        col = 'orange',
        lwd = 2,
        xlab = 'time',
        ylab = 'transition probability',
        main = 'Probability from S to U')
for(time in dummies) abline(v = time, col='grey')
dev.off()

### saving plot form U to I
png(filename = 'Probability_fromU_toI.png', height = 800, width = 800)
matplot(PPP_mean[6,3,],
        type = 'l',
        lty = c(1),
        col = 'orange',
        lwd = 2,
        xlab = 'time',
        ylab = 'transition probability',
        main = 'Probability from U to I')
for(time in dummies) abline(v = time, col='grey')
dev.off()



#### POSTERIOR ESTIMATE CONTINGENCY TABLES ####


mTAB <- apply(TTAB,1:3,mean)
lwTAB <- apply(TTAB,1:3,quantile,probs=0.025,na.rm=TRUE)
upTAB <- apply(TTAB,1:3,quantile,probs=0.975,na.rm=TRUE)

for(t in 2:TT){
  write.table(mTAB[,,t], paste0("tab/mean/",t,".txt"),row.names = categories, col.names = categories)
  write.table(lwTAB[,,t], paste0("tab/lower_bound/",t,".txt"),row.names = categories, col.names = categories)
  write.table(upTAB[,,t], paste0("tab/upper_bound/",t,".txt"),row.names = categories, col.names = categories)
}


for (h in 1:nrow(Ind)){
  
  i = Ind[h,1]; j = Ind[h,2]
  
  title <- paste("From", categories[i], "to", categories[j])
  x11()
  matplot(cbind(mTAB[i,j,], lwTAB[i,j,], upTAB[i,j,]) , 
          type = 'l',
          lty = c(1,2,2),
          col = c('dodgerblue2', 'olivedrab3', 'olivedrab3'),
          lwd = 1.5,
          xlab = 'time',
          ylab = 'counts',
          main = title)

  
}

View(Y)




### Recovered
Ind_red <- Ind[which(Ind[,2]==2),]
x11()
par(mfrow=c(3,2))
matplot(Y$recovered,
        type = 'l',
        lty = c(1),
        col = 'red',
        lwd = 1.5,
        xlab = 'time',
        ylab = 'counts',
        main = 'Recovered')

for(h in 1:nrow(Ind_red)){
  
  
  i = Ind_red[h,1]; j = Ind_red[h,2]
  
  title <- paste("From", categories[i], "to", categories[j]);
  matplot(cbind(mTAB[i,j,], lwTAB[i,j,], upTAB[i,j,]) , 
          type = 'l',
          lty = c(1,2,2),
          col = c('dodgerblue2', 'olivedrab3', 'olivedrab3'),
          lwd = 1.5,
          xlab = 'time',
          ylab = 'counts',
          main = title)
  
}



### Isolated
Ind_red <- Ind[which(Ind[,2]==3),]
x11()
par(mfrow=c(2,2))
matplot(Y$isolated,
        type = 'l',
        lty = c(1),
        col = 'red',
        lwd = 1.5,
        xlab = 'time',
        ylab = 'counts',
        main = 'Isolated')


for(h in 1:nrow(Ind_red)){
  
  
  i = Ind_red[h,1]; j = Ind_red[h,2]
  
  title <- paste("From", categories[i], "to", categories[j]);
  matplot(cbind(mTAB[i,j,], lwTAB[i,j,], upTAB[i,j,]) , 
          type = 'l',
          lty = c(1,2,2),
          col = c('dodgerblue2', 'olivedrab3', 'olivedrab3'),
          lwd = 1.5,
          xlab = 'time',
          ylab = 'counts',
          main = title)
  
}



### Hospitalized
Ind_red <- Ind[which(Ind[,2]==4),]
x11()
par(mfrow=c(2,2))
matplot(Y$hospitalized,
        type = 'l',
        lty = c(1),
        col = 'red',
        lwd = 1.5,
        xlab = 'time',
        ylab = 'counts',
        main = 'Hospitalized')


for(h in 1:nrow(Ind_red)){
  
  
  i = Ind_red[h,1]; j = Ind_red[h,2]
  
  title <- paste("From", categories[i], "to", categories[j]);
  matplot(cbind(mTAB[i,j,], lwTAB[i,j,], upTAB[i,j,]) , 
          type = 'l',
          lty = c(1,2,2),
          col = c('dodgerblue2', 'olivedrab3', 'olivedrab3'),
          lwd = 1.5,
          xlab = 'time',
          ylab = 'counts',
          main = title)
  
}


### Threatened
Ind_red <- Ind[which(Ind[,2]==5),]
x11()
par(mfrow=c(2,2))
matplot(Y$threatened,
        type = 'l',
        lty = c(1),
        col = 'red',
        lwd = 1.5,
        xlab = 'time',
        ylab = 'counts',
        main = 'Threatened')


for(h in 1:nrow(Ind_red)){
  
  
  i = Ind_red[h,1]; j = Ind_red[h,2]
  
  title <- paste("From", categories[i], "to", categories[j]);
  matplot(cbind(mTAB[i,j,], lwTAB[i,j,], upTAB[i,j,]) , 
          type = 'l',
          lty = c(1,2,2),
          col = c('dodgerblue2', 'olivedrab3', 'olivedrab3'),
          lwd = 1.5,
          xlab = 'time',
          ylab = 'counts',
          main = title)
  
}


### Undetected
Ind_red <- Ind[which(Ind[,2]==6),]
x11()
par(mfrow=c(2,2))
matplot(Y$undetected,
        type = 'l',
        lty = c(1),
        col = 'red',
        lwd = 1.5,
        xlab = 'time',
        ylab = 'counts',
        main = 'Undetected')


for(h in 1:nrow(Ind_red)){
  
  
  i = Ind_red[h,1]; j = Ind_red[h,2]
  
  title <- paste("From", categories[i], "to", categories[j]);
  matplot(cbind(mTAB[i,j,], lwTAB[i,j,], upTAB[i,j,]) , 
          type = 'l',
          lty = c(1,2,2),
          col = c('dodgerblue2', 'olivedrab3', 'olivedrab3'),
          lwd = 1.5,
          xlab = 'time',
          ylab = 'counts',
          main = title)
  
}


### Extinct
Ind_red <- Ind[which(Ind[,2]==7),]
x11()
par(mfrow=c(2,2))
matplot(Y$extinct,
        type = 'l',
        lty = c(1),
        col = 'red',
        lwd = 1.5,
        xlab = 'time',
        ylab = 'counts',
        main = 'Extinct')


for(h in 1:nrow(Ind_red)){
  
  
  i = Ind_red[h,1]; j = Ind_red[h,2]
  
  title <- paste("From", categories[i], "to", categories[j]);
  matplot(cbind(mTAB[i,j,], lwTAB[i,j,], upTAB[i,j,]) , 
          type = 'l',
          lty = c(1,2,2),
          col = c('dodgerblue2', 'olivedrab3', 'olivedrab3'),
          lwd = 1.5,
          xlab = 'time',
          ylab = 'counts',
          main = title)
  
}




### saving plot Isolated

Ind_red <- Ind[which(Ind[,2]==3),]
png(filename = 'Isolated.png', height = 800, width = 800)
#x11(width = 10, height = 10)
par(mfrow=c(2,2))
matplot(Y$isolated,
        type = 'l',
        lty = c(1),
        col = 'red',
        lwd = 2,
        xlab = 'time',
        ylab = 'counts',
        main = 'Isolated')


for(h in 1:nrow(Ind_red)){
  
  
  i = Ind_red[h,1]; j = Ind_red[h,2]
  
  title <- paste("From", categories[i], "to", categories[j]);
  matplot(mTAB[i,j,], 
          type = 'l',
          lty = c(1,2,2),
          col = c('dodgerblue2', 'olivedrab3', 'olivedrab3'),
          lwd = 2,
          xlab = 'time',
          ylab = 'counts',
          main = title)
  
}
dev.off()



### saving plot Undetected

Ind_red <- Ind[which(Ind[,2]==6),]
png(filename = 'Undetected.png', height = 800, width = 800)
par(mfrow=c(2,2))
matplot(Y$undetected,
        type = 'l',
        lty = c(1),
        col = 'red',
        lwd = 2,
        xlab = 'time',
        ylab = 'counts',
        main = 'Undetected')


for(h in 1:nrow(Ind_red)){
  
  
  i = Ind_red[h,1]; j = Ind_red[h,2]
  
  title <- paste("From", categories[i], "to", categories[j]);
  matplot(mTAB[i,j,], 
          type = 'l',
          lty = c(1,2,2),
          col = c('dodgerblue2', 'olivedrab3', 'olivedrab3'),
          lwd = 2,
          xlab = 'time',
          ylab = 'counts',
          main = title)
  
}
dev.off()

