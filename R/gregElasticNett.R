#Helper function to compute linear gregElasticNet total for bootstrapping
library(glmnet)

gregElasticNett <- function(data, x_pop_d, indices, alpha, lambda){
  #data: 1st column:y, 2nd column:pis, rest: x_sample_d
  d <- data[indices,]
  
  #y
  y <- d[,1]
  
  #pis 
  pis <- d[,2]
  
  #Length of x_sample_d
  p <- dim(d)[2] - 2
  #x_sample_d
  x_sample_d <- d[, 3:(p + 2)]
  
  #beta-hats
  pred_mod <- glmnet(x = as.matrix(x_sample_d[,-1]), y = y, alpha = alpha, family = "gaussian", standardize = FALSE, weights = pis^{-1})
  beta_hat <- predict(pred_mod, s = lambda, type = "coefficients")[1:dim(x_sample_d)[2],]

  return(beta_hat %*% (x_pop_d) + t(y - x_sample_d %*% beta_hat) %*% pis^(-1))
}

