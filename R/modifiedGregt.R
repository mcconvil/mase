modifiedGregt <- function(data,
                          xpopd,
                          ws,
                          domain,
                          domain_col_name,
                          indices) {
  
  # y, pi, xsample, domain
  d <- data[indices, ]
  
  y <- d[, 1]
  
  pis <- d[, 2]
  
  # last column contains domain labels
  xsamp <- d[, 3:(dim(d)[2])]
  xsample <- data.frame(xsamp[,-1, drop = FALSE])
  xsample_d <- as.matrix(xsamp[ , 1:(ncol(xsamp) - 1)])
  xsample_dt <- t(xsample_d)
  
  constant_component1 <- solve(xsample_dt %*% diag(ws) %*% xsample_d)
  constant_component2 <- t(ws * xsample_d)
  
  
  domain_indic_vec <- as.integer(xsample[domain_col_name] == domain)
  
  xpop_d_domain <- xpopd
  xsample_domain <- xsample[xsample[ , ncol(xsample)] == domain, , drop = FALSE]
  xsample_d_domain <- model.matrix(~., data = data.frame(xsample_domain[ , 1:(ncol(xsample) - 1)]))
  xsample_dt_domain <- t(xsample_d_domain)
  weights_domain <- ws[which(domain_indic_vec == 1)]
  
  w <- as.matrix(
    ws*domain_indic_vec + (
      t(as.matrix(xpop_d_domain) - xsample_dt_domain %*% weights_domain) %*%
        constant_component1
    ) %*%
      constant_component2
  )
  
  t <- w %*% y
  
  return(t)
  
}
