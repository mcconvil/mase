#Helper function to compute linear gregElasticNet total for bootstrapping
library(glmnet)

logisticGregElasticNett <- function(data, xpopd, indices, alpha, lambda){
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
  pred.mod <- glmnet(x = as.matrix(xsample_d[,-1]), y = y, alpha = alpha, family = "binomial", standardize = FALSE, weights = pis^{-1})
  
  #Predictions over the universe
  y.hats.U <- predict(pred.mod,newx = xpopd[,-1], s = lambda, type = "response")
 
  #Predictions over the sample
  y.hats.s <- predict(pred.mod,newx = xsample_d[,-1], s = lambda, type = "response")
  
  
  #Compute and return estimator
  return(sum(y.hats.U) + t(y - y.hats.s) %*% pis^(-1))
}

