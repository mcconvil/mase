#Helper function to compute linear GREG total for bootstrapping

gregt <- function(data, xpopd, indices){
  #data: 1st column:y, 2nd column:pis, rest: xsample_d
  d <- data[indices,]
  
  #y
  y <- d[,1]
  
  #pis 
  pis <- d[,2]
  
  #Length of xsample_d
  p <- dim(d)[2] - 2
  #xsample_d
  xsample_d <- d[, 3:(p+2)]
  
  #Survey weight
  w <- as.matrix(1 + t(as.matrix(xpopd) - t(xsample_d) %*% pis^{-1} ) %*% solve(t(xsample_d) %*% diag(pis^{-1}) %*% xsample_d) %*% t(xsample_d)) %*% diag(pis^{-1})
  
  return(w %*% as.vector(y))
}

