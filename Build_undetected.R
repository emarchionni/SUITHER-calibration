#### Building undetected time series

#'@param isolated time series of isolated individuals
#'@param hospitalized time series of extict individuals
#'@param threatened time series of threatened individuals
#'@param CFR time series of CFR in accordance with the confirmation to death delay
#'@param IFR infection fatality ratio
#'@param time time instants
#'
#'
#'@return time series of CFR (case fatality ratio) for the central time instants (tails: dt/2)


build_undetected <- function(isolated, hospitalized, threatened, CFR, IFR, time){
  
  if(!(length(isolated)==length(hospitalized) || length(hospitalized)==length(threatened) || 
       length(CFR) == length(isolated)))
    stop('Dimensions do not agree')
  
  undetected <- (CFR/IFR - 1) * (isolated + hospitalized + threatened)
  undetected <- floor(undetected)
  
  undetected <- data.frame(undetected)
  undetected <- cbind(time,undetected)
  colnames(undetected) <- c('date',
                           'undetected')
  return(undetected)
}