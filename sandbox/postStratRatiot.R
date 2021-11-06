#Helper function to compute postStrat total for bootstrapping

postStratRatiot <- function(data, xpop, indices){

  
    #data: 4th column: y_den, 3rd column:y_num, 2nd column:pis, 1st column: xsample
  d <- data[indices,]
  
  #y
  y_num <- d[, 3]
  y_den <- d[, 4]
  
  #pis 
  pis <- d[, 2]
  
  #xsample
  xsample <- d[, 1]
  
  #Estimator
  y_num_ps <- postStrat(y = y_num, xsample = xsample,
                        xpop = xpop)
  y_den_ps <- postStrat(y = y_den, xsample = xsample,
                        xpop = xpop)
  
  #Ratio: Point Estimate
  out <- y_num_ps$pop_total/y_den_ps$pop_total

  return(out)
}

