#### Building recovered time series

#'@param IFR infection fatality ratio
#'@param extinct time series of extict people
#'@param time time instants
#'
#'
#'@return time series of recovered people

build_recovered <- function(IFR, time, extinct){
  recovered <- (1/IFR - 1) * extinct
  recovered <- floor(recovered)
  
  recovered <- data.frame(recovered)
  recovered <- cbind(time,recovered)
  colnames(recovered) <- c('date',
                           'recovered')
  return(recovered)
}