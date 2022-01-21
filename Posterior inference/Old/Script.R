remove(list = ls())

library(coda)

setwd('C:/Users/edoar/Desktop/NAPDE/Project/Code/SUITHER-calibration/Posterior inference')

load('C:/Users/edoar/Desktop/NAPDE/Project/Code/model_7.RData')

Ind <- rbind(cbind(1,c(1,6)), cbind(2,2), cbind(3,c(2:4,7)),
             cbind(4,c(2,4,5)), cbind(5,c(4,5,7)),
             cbind(6,c(2,3,6)), c(7,7))

names(out_no_ORmax)

attach(out_no_ORmax)


View(Y)

# figure data
categories <- c('R',
                'I',
                'H',
                'T',
                'U',
                'E')

col = c("royalblue2",
        "red",
        "olivedrab4",
        "gold2",
        "gray19",
        "wheat4")

ytick<-c(0,22000,45000,65500,90000, 108300)
xtick<-c(1,11,21,31,41,51,61,71)

x11()
matplot(Y[,c(-1,-2)],
        xlab=expression("Day ("*italic(t)*")"), 
        type="l",
        col = col,
        lwd=2, 
        lty=1,
        ylab = 'Observed frequencies')

legend("topleft",
       legend = categories,
       col = col,
       pch = 20, bty = "n", 
       x.intersp = 0.8,
       cex= 1,  pt.cex = 1,
       xpd = TRUE, text.width = 0.003)




### posterior TAB

niter <- dim(TTAB)[4]

PPP_mean <- apply(PPP[], margin, ...)

for(it in 1:niter){
  
}
