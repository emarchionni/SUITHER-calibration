#### Building CFR time series

#'@param recovered.D recovered from detected time series
#'@param extinct time series of extict people
#'@param time time instants
#'@param dt time-window
#'
#'
#'@return time series of CFR (case fatality ratio) for the central time instants (tails: dt/2)


build_CFR <- function(recovered.D, extinct, time, dt){
  
  if(!length(recovered.D)==length(extinct))
    stop('Dimensions do not agree')
  
  n <- length(recovered.D)  
  first_time <- dt/2 + 1 # each tail = dt/2
  last_time <- n - dt/2  # tot time instants = n - dt
  
  CFR <- c()
    
  for(i in first_time:last_time){
    delta_e <- extinct[i + dt/2] - extinct[i - dt/2]
    delta_r.D <- recovered.D[i + dt/2] - recovered.D[i - dt/2]
    CFR <- c(CFR, delta_e/(delta_r.D+delta_e))
  }
  
  time <- time[first_time:last_time]

  CFR <- data.frame(CFR)
  CFR <- cbind(time,CFR)
  colnames(CFR) <- c('date',
                    'CFR')
  
  
  return(CFR)
}