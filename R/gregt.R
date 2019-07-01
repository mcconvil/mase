#Helper function to compute linear GREG total for bootstrapping

gregt <- function(data, y,  x_pop_d, x_sample_d, indices){
  #data: 1st column:y, 2nd column:pis, rest: x_sample_d
  d <- data[indices,]
  #resids
  e <- d[,1]
  #pis 
  pis <- data[,2]
  #y star
  y_star <- y + e
  
  #Survey weight
  w <- as.matrix(1 + t(as.matrix(x_pop_d) - t(x_sample_d) %*% pis^{-1} ) %*% solve(t(x_sample_d) %*% diag(pis^{-1}) %*% x_sample_d) %*% t(x_sample_d)) %*% diag(pis^{-1})
  
  return(w %*% as.vector(y_star))
}

