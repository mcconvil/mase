#Helper function to compute linear gregElasticNet total for bootstrapping
library(glmnet)

logisticGregElasticNett <- function(data, x_pop_d, indices, alpha, lambda){
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
  pred_mod <- glmnet(x = as.matrix(x_sample_d[,-1]), y = y, alpha = alpha, family = "binomial", standardize = FALSE, weights = pis^{-1})
  
  #Predictions over the universe
  y_hats_U <- predict(pred_mod,newx = x_pop_d[,-1], s = lambda, type = "response")
 
  #Predictions over the sample
  y_hats_s <- predict(pred_mod,newx = x_sample_d[,-1], s = lambda, type = "response")
  
  
  #Compute and return estimator
  return(sum(y_hats_U) + t(y - y_hats_s) %*% pis^(-1))
}

