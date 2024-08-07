#' Compute a ratio of two estimators
#' 
#' @inheritParams greg
#' @param y_num A vector containing the response value for each sampled unit in the numerator
#' @param y_den A vector containing the response value for each sampled unit in the denominator
#' @param estimator A string containing the name of the estimators of which you are taking a ratio of. The names follow the same format as the functions independently do in mase. Options are "horvitzThompson", "postStrat", and "greg".
#' @param datatype Default to "raw", takes values "raw", "totals" or "means" for whether the user is providing the raw population stratum memberships, the population totals of each stratum, or the population proportions of each stratum.
#' @param ... Any additional arguments that can be passed to mase::horvitzThompson, mase::greg, and mase::postStrat

#' @examples 
#' library(survey)
#' data(api) 
#' ratio(y_num = apisrs$api.stu, y_den = apisrs$enroll, xsample = apisrs$stype,
#' xpop = apipop$stype, pi = apisrs$pw^(-1), estimator = "postStrat",
#' var_est = TRUE, var_method = "LinHB", datatype = "raw")
#' 
#' @references 
#' \insertRef{coc77}{mase} 
#' \insertRef{sar92}{mase}
#'
#' @return 
#' A list of output containing:
#'
#' * ratio_est: Estimate of the ratio of the population totals/means of the two estimators.
#' 
#' * ratio_var_est: Estimate of the variance of the ratio of two estimators.
#' 
#' @importFrom rlang check_dots_used
#' @export ratio
#' @include bootHelpers.R
#' @include varMase.R

ratio <- function(y_num,
                  y_den,
                  xsample,
                  xpop,
                  pi = NULL,
                  pi2 = NULL,
                  N = NULL,
                  estimator = NULL,
                  var_est = F, 
                  var_method = "LinHTSRS",
                  datatype = "raw",
                  fpc = TRUE,
                  messages = TRUE,
                  ...){
  
  check_dots_used()
  
  if(is.null(estimator) || !(estimator %in% c("horvitzThompson", "postStrat", "greg"))) {
    stop("Must supply a valid value for argument 'estimator' \n  Support currently only exists for the following estimators: \n \t - horvitzThompson \n \t - postStrat \n \t - greg")
  }
  
  if (is.null(pi)) {
    if (messages) {
      message("Assuming simple random sampling") 
    }
    pi <- rep(length(y_num)/N, length(y_num))
  }
  
  weight <- as.vector(pi^(-1))
  
  if (estimator == "horvitzThompson") {
    
    ht_num <- horvitzThompson(
      y = y_num,
      pi = pi,
      N = N,
      pi2 = pi2,
      messages = messages,
      fpc = fpc,
      var_est = FALSE,
      ...
    )
    
    ht_den <- horvitzThompson(
      y = y_den,
      pi = pi,
      N = N,
      pi2 = pi2,
      messages = messages,
      fpc = fpc,
      var_est = FALSE,
      ...
    )
    
    rat <- ht_num$pop_total/ht_den$pop_total
    e_ratio <- y_num - as.vector(rat)*y_den
    est_den <- ht_den$pop_total
    
    
  } else if (estimator == "postStrat") {
    
    ps_num <- postStrat(
      y = y_num, 
      xsample = xsample,
      xpop = xpop,
      pi = pi,
      N = N,
      pi2 = pi2, 
      datatype = datatype,
      messages = messages,
      fpc = fpc,
      var_est = FALSE,
      ...
    )
    
    ps_den <- postStrat(
      y = y_den,
      xsample = xsample,
      xpop = xpop,
      pi = pi,
      N = N,
      pi2 = pi2,
      datatype = datatype,
      messages = messages,
      fpc = fpc,
      var_est = FALSE,
      ...
    )
    
    est_den <- ps_den$pop_total
    
    rat <- ps_num$pop_total/ps_den$pop_total
    
    m_hat <- function(y) {
      
      xy <- data.frame(y = y, xsample = xsample, pi = pi)
      strata_sums <- data.frame()
      
      for (i in unique(xsample)) {
        
        N_h <- length(xpop[xpop == i])
        
        strata_h <- xy[xy$xsample == i, ]
        strata_h$ratio_h <- strata_h$y/strata_h$pi
        sum_h <- (N_h^(-1))*sum(strata_h$ratio_h)
        row_h <- data.frame(beta_h = sum_h, index = i)
        strata_sums <- rbind(strata_sums, row_h)
        
      }
      
      m_hat_xi <- data.frame(index = xsample, id = 1:length(xsample))
      y_hat_unordered <- merge(m_hat_xi, strata_sums, by = "index")
      y_hat <- y_hat_unordered[order(y_hat_unordered$id), ]
      
      return(y_hat$beta_h)
      
    }
    
    y_hat_num <- m_hat(y_num)
    y_hat_den <- m_hat(y_den)
    
    e_num <- y_num - y_hat_num
    e_den <- y_den - y_hat_den
    
    e_ratio <- e_num - as.vector(rat)*e_den
    
  } else if (estimator == "greg") {
    
    greg_num <- greg(
      y = y_num,
      xsample = xsample,
      xpop = xpop,
      pi = pi,
      N = N,
      messages = messages,
      fpc = fpc,
      var_est = FALSE,
      ...
    )
    
    greg_den <- greg(
      y = y_den,
      xsample = xsample,
      xpop = xpop,
      pi = pi,
      N = N,
      messages = messages,
      fpc = fpc,
      var_est = FALSE,
      ...
    )
    
    est_den <- greg_den$pop_total
    
    rat <- greg_num$pop_total/greg_den$pop_total
    
    xsample_d_var <- model.matrix(~., data = data.frame(xsample))
    xsample_var <- data.frame(xsample_d_var[,-1, drop = FALSE])
    xsample_dt_var <- t(xsample_d_var) 
    
    y_hat_const <- xsample_d_var %*% solve(xsample_dt_var %*% diag(weight) %*% xsample_d_var) %*% 
      (xsample_dt_var) %*% diag(weight) 
    
    y_hat_num <- y_hat_const %*% y_num
    
    y_hat_den <- y_hat_const %*% y_den
    
    e_num <- y_num - as.vector(y_hat_num)
    e_den <- y_den - as.vector(y_hat_den)
    
    e_ratio <- e_num - as.vector(rat)*e_den
    
  }
  
  if (var_est == TRUE) {
    
    var_est <- varMase(
      y = as.numeric(e_ratio),
      pi = pi,
      pi2 = pi2,
      method = var_method,
      N = N,
      fpc = fpc
    )
    
    return(list(ratio_est = rat,
                ratio_var_est = (1/est_den^2)*var_est))
    
  } else {
    
    return(list(ratio_est = rat))
    
  }
  
}
