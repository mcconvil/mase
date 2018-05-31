#Helper function to compute HT total for bootstrapping

htt <- function(data,  indices){
  #data: 1st column:y, 2nd column:pis, rest: x_sample_d
  d <- data[indices,]
  
  #y
  y <- d[,1]
  
  #pis 
  pis <- d[,2]
  
  return(pis^{-1} %*% as.vector(y))
}

