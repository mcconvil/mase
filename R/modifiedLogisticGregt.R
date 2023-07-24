modifiedLogisticGregt <- function(data,
                                 xpopd,
                                 xpop_sums,
                                 ws,
                                 domain,
                                 domain_col_name,
                                 lab,
                                 indices) {
  
  d <- data[indices, ]
  y <- d[, 1]
  pis <- d[, 2]
  xsamp <- d[, 3:(dim(d)[2])]
  xsample <- data.frame(xsamp[,-1, drop = FALSE])
  xsample_d <- as.matrix(xsamp[ , 1:(ncol(xsamp) - 1)])
  xsample_dt <- t(xsample_d)
  
  xsample_preds <- xsample[ , !(names(xsample) %in% domain_col_name)]
  dat <- data.frame(y, ws, xsample_preds)
  colnames(dat) <- c("y", "weight", names(xsample_preds))
  f <- paste(names(dat)[1], "~", paste(names(dat)[-c(1,2)], collapse = " + "))
  s_design <- survey::svydesign(ids = ~1, weights = ~weight, data = dat)
  mod <- survey::svyglm(f, design = s_design, family = quasibinomial())
  
  
  ## work
  
  domain_indic_vec <- as.integer(xsample[domain_col_name] == domain)
  
  xpop_domain <- xpopd
  xsample_domain <- xsample[xsample[domain_col_name] == domain, , drop = FALSE]
  domain_N <- xpop_sums[xpop_sums[domain_col_name] == domain, , drop = FALSE][lab]
  
  y_hats_U <- as.matrix(predict(mod, newdata = xpop_domain[ , -1], type = "response", family = quasibinomial()))
  y_hats_s <- as.matrix(predict(mod, newdata = xsample_domain, type = "response", family = quasibinomial()))
  
  y_domain <- y[which(domain_indic_vec == 1)]
  weights_domain <- ws[which(domain_indic_vec == 1)]
  
  t <- t(y_domain - y_hats_s) %*% weights_domain + sum(y_hats_U)
  
  return(t)
  
}
