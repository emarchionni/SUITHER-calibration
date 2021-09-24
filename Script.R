remove(list = ls())


setwd('C:/Users/edoar/Desktop/NAPDE/Project/Code/SUITHER-calibration')

npop <- 60317000

italy <- read.csv('dpc-covid19-ita-andamento-nazionale.csv')
#View(italy)

beginning <- '2020-08-20T17:00:00'
end <- '2020-12-31T17:00:00'

beginning_row <- which(italy$data==beginning)
end_row <- which(italy$data==end)

italy_peak <- italy[beginning_row:end_row,]
italy_peak <- data.frame(italy_peak)


#### Building dataset ####
Y <- matrix(data = NA, nrow = dim(italy_peak)[1],  ncol = 8)
colnames(Y) <- c('date', 
                 'susceptible',
                 'recovered',
                 'isolated',
                 'hospitalized',
                 'threatened',
                 'undetected',
                 'extinct')
Y <- data.frame(Y)

Y$date <- italy_peak$data
Y$isolated <- italy_peak$isolamento_domiciliare
Y$hospitalized <- italy_peak$ricoverati_con_sintomi
Y$threatened <- italy_peak$terapia_intensiva
Y$extinct <- italy_peak$deceduti

#### recovered-time series
source('Build_recovered.R')
recovered <- build_recovered(IFR = 1.2/100,
                             time = Y$date,
                             extinct = Y$extinct)

Y$recovered <- recovered[,2]
# View(cbind(Y$recovered,italy_peak$dimessi_guariti))

#### undetected-time series

### CFR
dt <- 28

# CFR needed from 02/09 to 13/01 (confirmation to death delay d = 13)
beginning <- '2020-09-02T17:00:00'
end <- '2021-01-13T17:00:00'

beginning_row <- which(italy$data==beginning) - dt/2
end_row <- which(italy$data==end) + dt/2

italy_peak2 <- italy[beginning_row:end_row,]
italy_peak2 <- data.frame(italy_peak2)

source('Build_CFR.R')
CFR <- build_CFR(recovered.D = italy_peak2$dimessi_guariti,
                 extinct = italy_peak2$deceduti,
                 time = italy_peak2$data,
                 dt = 28)

### undetected 
source('Build_undetected.R')

undetected <- build_undetected(isolated = Y$isolated,
                               hospitalized = Y$hospitalized,
                               threatened = Y$threatened,
                               CFR = CFR[,2],
                               IFR = 1.2/100,
                               time = Y$date)
Y$undetected <- undetected[,2]


#### susceptible time series

Y$susceptible <- npop - rowSums(Y[,3:8])

View(Y)


#### Preliminary plots ####
x11()
par(mfrow=c(4,2))
titles <- colnames(Y)
for (i in 2:8) plot(Y[,i], xlab = 'Time', ylab = '', main = titles[i], type = 'p', col = 'grey', pch = 19)


#### MCMC ####
source('MCMC/MCMC_function.R')

