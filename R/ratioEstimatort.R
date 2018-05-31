#Helper function to compute ratio Estimator total for bootstrapping

ratioEstimatort <- function(data, tau_x, indices){
  #data: 1st column:y, 2nd column:pis, rest: x_sample_d
  d <- data[indices,]
  
  #y
  y <- d[,1]
  
  #pis 
  pis <- d[,2]
  
  #x_sample
  x_sample <- d[, 3]
  
  #Estimator
  tyHT <- horvitzThompson(y=y,pi=pis)$pop_total
  txHT <- horvitzThompson(y=x_sample,pi=pis)$pop_total
  
  return(as.vector(tau_x/txHT*tyHT))
}

