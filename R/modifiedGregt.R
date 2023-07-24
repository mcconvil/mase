modifiedGregt <- function(data, xpopd, weight, domain, domain_labels, domain_col_name, indices) {
  
  # y, pi, xsample, domain
  d <- data[indices, ]
  
  y <- d[, 1]
  
  pis <- d[, 2]
  
  # last column contains domain labels
  xsamp <- d[, 3:(dim(d)[2])]
  xsample <- cbind(data.frame(xsamp[,-1, drop = FALSE]), domain_labels)
  names(xsample) <- c(colnames(xsamp[,-1, drop = FALSE]), domain_col_name)
  xsample_d <- as.matrix(xsamp)
  print(xsample_d)
  xsample_dt <- t(xsample_d)
  
  constant_component1 <- solve(xsample_dt %*% diag(weight) %*% xsample_d)
  constant_component2 <- t(weight * xsample_d)
  
  
  domain_indic_vec <- as.integer(xsample[ , ncol(xsample)] == domain)
  
  # xpop_domain <- xpop_d[xpop_d[domain_col_name] == domain_id, , drop = FALSE]
  xpop_d_domain <- xpopd
  xsample_domain <- xsample[xsample[ , ncol(xsample)] == domain, , drop = FALSE]
  xsample_d_domain <- model.matrix(~., data = data.frame(xsample_domain[ , 1:(ncol(xsample) - 1)]))
  xsample_dt_domain <- t(xsample_d_domain)
  weights_domain <- weight[which(domain_indic_vec == 1)]
  
  w <- as.matrix(
    weight*domain_indic_vec + (
      t(as.matrix(xpop_d_domain) - xsample_dt_domain %*% weights_domain) %*%
        constant_component1
    ) %*%
      constant_component2
  )
  
  t <- w %*% y
  
  return(t)
  
}
