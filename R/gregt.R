#Helper function to compute linear GREG total for bootstrapping

gregt <- function(data, x_pop_d, indices){
  #data: 1st column:y, 2nd column:pis, rest: x_sample_d
  d <- data[indices,]
  
  #y
  y <- d[,1]
  
  #pis 
  pis <- d[,2]
  
  #Length of x_sample_d
  p <- dim(d)[2] - 2
  #x_sample_d
  x_sample_d <- d[, 3:(p+2)]
  
  #Survey weight
  w <- as.matrix(1 + t(as.matrix(x_pop_d) - t(x_sample_d) %*% pis^{-1} ) %*% solve(t(x_sample_d) %*% diag(pis^{-1}) %*% x_sample_d) %*% t(x_sample_d)) %*% diag(pis^{-1})
  
  return(w %*% as.vector(y))
}

