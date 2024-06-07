gregt <- function(data,
                  xpopd,
                  domain_id,
                  indices){
  
  d <- data[indices, ]
  y <- d[,1]
  pis <- d[,2]
  weights <- as.vector(pis^(-1))
  p <- dim(d)[2] - 2
  xsample_d <- d[, 3:(p+2)]
  
  one_mat <- matrix(rep(1, times = nrow(xsample_d)), nrow = 1)
  xpop_cpp <- as.matrix(xpopd)
  weight_mat <- diag(weights)
  
  w <- get_weights_greg(xpop_cpp, xsample_d, weight_mat, one_mat)
  
  out <- sum(as.numeric(w) * as.numeric(y))
  
  return(out)
  
}

logisticGregt <- function(data,
                          xpopd,
                          indices){

  d <- data[indices,]
  
  f <- paste(names(d)[1], "~", paste(names(d)[-c(1,2)], collapse = " + "))
  s.design <- survey::svydesign(ids=~1,weights=~weight,data=d)
  
  mod <- survey::svyglm(f, design = s.design, family = quasibinomial())
  
  y.hats.U <- as.matrix(predict(mod,newdata = data.frame(xpopd[,-1]), type = "response", family = quasibinomial()))
  y.hats.s <- as.matrix(predict(mod, type = "response", family = quasibinomial()))
  
  out <- t(d$y - y.hats.s) %*% d$weight + sum(y.hats.U)

  return(out)
  
}


htt <- function(data,
                indices){

  d <- data[indices,]
  y <- d[,1]
  pis <- d[,2]
  
  out <- pis^(-1) %*% as.vector(y)
  
  return(out)
  
}


postStratt <- function(data,
                       xpoptab, 
                       indices){
  
  x <- N_h <- N_h_hats <- ps_h <- poptotal_h <- strat_pop_total <- NULL  
  
  d <- data[indices,]
  y <- d[,3]
  pis <- d[,2]
  xsample <- d[, 1]
  
  dat <- data.frame(x=xsample, pi=pis,y=y)
  tab <- dat |>
    group_by(x) |>
    summarize(poptotal_h = y%*%pi^(-1),N_h_hats = sum(pi^(-1)), var_h = var(y)) |>
    inner_join(xpoptab, by=c("x")) |>
    mutate(ps_h = N_h/N_h_hats, strat_pop_total = ps_h*poptotal_h, strat_pop_mean = strat_pop_total/N_h)
  
  out <- sum(tab$strat_pop_total)
  
  return(out)
  
}

gregElasticNett <- function(data,
                            xpopd,
                            indices,
                            alpha,
                            lambda){

  d <- data[indices,]
  y <- d[,1]
  pis <- d[,2]
  p <- dim(d)[2] - 2
  xsample_d <- d[, 3:(p + 2)]
  
  pred.mod <- glmnet(x = as.matrix(xsample_d[,-1]), y = y,
                     alpha = alpha, family = "gaussian",
                     standardize = FALSE, weights = pis^(-1))
  
  beta_hat <- predict(pred.mod, s = lambda, type = "coefficients")[1:dim(xsample_d)[2],]
  
  out <- beta_hat %*% (xpopd) + t(y - xsample_d %*% beta_hat) %*% pis^(-1)
  
  return(out)
  
}

logisticGregElasticNett <- function(data,
                                    xpopd,
                                    indices, 
                                    alpha,
                                    lambda) {

  d <- data[indices,]
  
  y <- d[,1]
  pis <- d[,2]
  p <- dim(d)[2] - 2
  xsample_d <- d[, 3:(p + 2)]
  
  pred.mod <- glmnet(x = as.matrix(xsample_d[,-1]),
                     y = y, alpha = alpha, family = "binomial",
                     standardize = FALSE, weights = pis^(-1))
  
  y.hats.U <- predict(pred.mod,newx = xpopd[,-1], s = lambda, type = "response")
  y.hats.s <- predict(pred.mod,newx = xsample_d[,-1], s = lambda, type = "response")
  
  out <- sum(y.hats.U) + t(y - y.hats.s) %*% pis^(-1)
  
  return(out)
  
}

gregTreet <- function(data, xpop, pval= pval, perm_reps = perm_reps, bin_size = bin_size, indices){

  d <- data[indices,]
  y <- d[,1]
  pis <- d[,2]
  
  p <- dim(d)[2] - 2
  xsample <- d[, 3:(p + 2)]

  f <- as.formula(paste("y ~ ", paste(names(xsample), collapse= "+")))
  treet <- rpms(rp_equ = f, data = d, weights = as.vector(pis^(-1)),
                pval= pval, perm_reps = perm_reps, bin_size = bin_size)
  
  xpop <- xpop[names(xsample)]
  xpop_treet <- box_ind(treet, xpop)
  xsample_treet <- box_ind(treet, xsample)
  
  w <- (1 + t(colSums(xpop_treet)- colSums(xsample_treet*pis^(-1))) %*%
          solve(t(xsample_treet)%*%diag(pis^(-1)) %*% 
                  as.matrix(xsample_treet))%*%t(xsample_treet)) %*% diag(pis^(-1)) 
  
  out <- w %*% y
  
  return(out)
  
}

modifiedGregt <- function(data,
                          xpopd,
                          ws,
                          domain,
                          domain_col_name,
                          indices) {
  
  d <- data[indices, ]
  y <- d[, 1]
  pis <- d[, 2]
  
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
  
  out <- w %*% y
  
  return(out)
  
}

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
  
  domain_indic_vec <- as.integer(xsample[domain_col_name] == domain)
  
  xpop_domain <- xpopd
  xsample_domain <- xsample[xsample[domain_col_name] == domain, , drop = FALSE]
  domain_N <- xpop_sums[xpop_sums[domain_col_name] == domain, , drop = FALSE][lab]
  
  y_hats_U <- as.matrix(predict(mod, newdata = xpop_domain[ , -1], type = "response", family = quasibinomial()))
  y_hats_s <- as.matrix(predict(mod, newdata = xsample_domain, type = "response", family = quasibinomial()))
  
  y_domain <- y[which(domain_indic_vec == 1)]
  weights_domain <- ws[which(domain_indic_vec == 1)]
  
  out <- t(y_domain - y_hats_s) %*% weights_domain + sum(y_hats_U)
  
  return(out)
  
}

ratioEstimatort <- function(data,
                            tau_x,
                            indices){

  d <- data[indices,]
  y <- d[,1]
  pis <- d[,2]
  xsample <- d[, 3]
  
  tyHT <- horvitzThompson(y = y, pi = pis)$pop_total
  txHT <- horvitzThompson(y = xsample, pi = pis)$pop_total
  
  out <- as.vector(tau_x/txHT*tyHT)
  
  return(out)
  
}
