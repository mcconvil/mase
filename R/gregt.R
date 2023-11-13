#Helper function to compute linear GREG total for bootstrapping

gregt <- function(data, xpopd, domain_id, indices){
  
  #data: 1st column:y, 2nd column:pis, rest: xsample_d
  d <- data[indices, ]
  
  #y
  y <- d[,1]
  
  #pis 
  pis <- d[,2]
  weights <- as.vector(pis^(-1))
  
  #Length of xsample_d
  p <- dim(d)[2] - 2
  #xsample_d
  xsample_d <- d[, 3:(p+2)]
  
  one_mat <- matrix(rep(1, times = nrow(xsample_d)), nrow = 1)
  xpop_cpp <- as.matrix(xpopd)
  weight_mat <- diag(weights)
  
  w <- get_weights(xpop_cpp, xsample_d, weight_mat, one_mat)
  
  # Survey weight
  # w <- as.matrix(1 + t(as.matrix(xpopd) - t(xsample_d) %*% pis^{-1} ) %*% solve(t(xsample_d) %*% diag(pis^{-1}) %*% xsample_d) %*% t(xsample_d)) %*% diag(pis^{-1})
  t <- sum(as.numeric(w) * as.numeric(y))
  return(t)
}

