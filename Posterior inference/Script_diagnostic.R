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



#### Betas coefficients ####

component <- 6
for (h in 1:nrow(Ind)) {
  
  i <- Ind[h,1]; j <- Ind[h,2]
  x11()
  traceplot(as.mcmc(BBE[i,j,component,]), 
            main = paste(paste0(component,"-th component"), "beta coefficient form", categories[i], "to", categories[j]))
  
}

### saving plots
png(filename = 'Traceplot_beta_4_5.png', height = 800, width = 800)
par(mfrow=c(2,2))
traceplot(as.mcmc(BBE[1,1,4,]), 
          main = paste(paste0(4,"th component"), "beta coefficient form S to S"))
traceplot(as.mcmc(BBE[1,6,4,]), 
          main = paste(paste0(4,"th component"), "beta coefficient form S to U"))
traceplot(as.mcmc(BBE[1,1,5,]), 
          main = paste(paste0(5,"th component"), "beta coefficient form S to S"))
traceplot(as.mcmc(BBE[1,6,5,]), 
          main = paste(paste0(5,"th component"), "beta coefficient form S to U"))
title('Beta 4 and beta 5', outer = T, line = -2)
dev.off()



#### Transition probabilities ####

time <- 70

for (h in 1:nrow(Ind)) {
  
  i <- Ind[h,1]; j <- Ind[h,2]
  x11()
  traceplot(as.mcmc(PPP[i,j,2,]), 
            main = paste("Probability form", categories[i], "to", categories[j]))
  
}


### saving plot

# time 2
png(filename = 'Traceplot_time2.png', height = 800, width = 800)
par(mfrow=c(2,2))
traceplot(as.mcmc(PPP[1,1,2,]), main = "Probability from S to S" )
traceplot(as.mcmc(PPP[1,6,2,]), main = "Probability from S to U" )
traceplot(as.mcmc(PPP[4,5,2,]), main = "Probability from H to T" )
traceplot(as.mcmc(PPP[5,4,2,]), main = "Probability from T to H")
title('Time 2', outer = T, line = -2)
dev.off()



