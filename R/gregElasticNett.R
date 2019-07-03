#Helper function to compute linear gregElasticNet total for bootstrapping
library(glmnet)

gregElasticNett <- function(data, x_pop_d, x_sample_d, y_hat, indices, alpha, lambda){
  #data: 1st column:resids, 2nd:pis
  #resids, wild bootstrap Mammen's dist.
  e <- data[,1]*ifelse(rbinom(length(y_hat), 1, 0.7236), -0.618, 1.618 )
  #pis 
  pis <- data[,2]
  #y star
  y_star <- y_hat + e
  
  #beta-hats
  pred_mod <- glmnet(x = as.matrix(x_sample_d[,-1]), y = y_star, alpha = alpha, family = "gaussian", standardize = FALSE, weights = pis^{-1})
  beta_hat <- predict(pred_mod, s = lambda, type = "coefficients")[1:dim(x_sample_d)[2],]

  return(beta_hat %*% (x_pop_d) + t(y_star - x_sample_d %*% beta_hat) %*% pis^(-1))
}

