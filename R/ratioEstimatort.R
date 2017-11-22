#Helper function to compute ratio Estimator total for bootstrapping

ratioEstimatort <- function(data, tau_x, indices){
  #data: 1st column:y, 2nd column:pis, rest: xsample_d
  d <- data[indices,]
  
  #y
  y <- d[,1]
  
  #pis 
  pis <- d[,2]
  
  #xsample
  xsample <- d[, 3]
  
  #Estimator
  tyHT <- horvitzThompson(y=y,pi=pis)$pop_total
  txHT <- horvitzThompson(y=xsample,pi=pis)$pop_total
  
  return(as.vector(tau_x/txHT*tyHT))
}

