#Helper function to compute linear GREG total for bootstrapping

gregt <- function(data, y_hat,  x_pop_d, x_sample_d, indices){
  #data: 1st column:y, 2nd column:pis, rest: x_sample_d
  #resids, wild bootstrap Mammen's dist.
  e <- data[,1]*ifelse(rbinom(length(y_hat), 1, 0.7236), -0.618, 1.618 )
  #pis 
  pis <- data[,2]
  #y star
  y_star <- y_hat + e
  
  #Survey weight
  w <- as.matrix(1 + t(as.matrix(x_pop_d) - t(x_sample_d) %*% pis^{-1} ) %*% solve(t(x_sample_d) %*% diag(pis^{-1}) %*% x_sample_d) %*% t(x_sample_d)) %*% diag(pis^{-1})
  
  return(w %*% as.vector(y_star))
}

