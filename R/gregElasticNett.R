#Helper function to compute linear gregElasticNet total for bootstrapping
library(glmnet)

gregElasticNett <- function(data, xpopd, indices, alpha, lambda){
  #data: 1st column:y, 2nd column:pis, rest: xsample_d
  d <- data[indices,]
  
  #y
  y <- d[,1]
  
  #pis 
  pis <- d[,2]
  
  #Length of xsample_d
  p <- dim(d)[2] - 2
  #xsample_d
  xsample_d <- d[, 3:(p + 2)]
  
  #beta-hats
  pred.mod <- glmnet(x = as.matrix(xsample_d[,-1]), y = y, alpha = alpha, family = "gaussian", standardize = FALSE, weights = pis^{-1})
  beta_hat <- predict(pred.mod, s = lambda, type = "coefficients")[1:dim(xsample_d)[2],]

  return(beta_hat %*% (xpopd) + t(y - xsample_d %*% beta_hat) %*% pis^(-1))
  
}

