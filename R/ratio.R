#' Compute a ratio of two estimators
#' 
#' @inheritParams horvitzThompson
#' @param y_num A vector containing the response value for each sampled unit in the numerator
#' @param y_den A vector containing the response value for each sampled unit in the denominator
#' @param xsample A vector containing the appropriate form of xsample for the estimator of choice. For example see ?mase::greg() to see the appropriate input for xsample when computing a ratio of two greg estimators
#' @param xpop A vector containing the appropriate form of xpop for the estimator of choice.
#' @param pi A vector containing the first order inclusion probabilities for the sample.
#' @param pi2 Defaults to NULL. A vector containing the second order inclusion probabilities for the sample.
#' @param N The number of observations in the population.
#' @param estimator A string containing the name of the estimators of which you are taking a ratio of. The names follow the same format as the functions independently do in mase. Options are "horvitzThompson", "postStrat", and "greg".
#' @param datatype Default to "raw", takes values "raw", "totals" or "means" for whether the user is providing the raw population stratum memberships, the population totals of each stratum, or the population proportions of each stratum.
#' @param var_est Default to FALSE, logical for whether or not to compute estimate of variance
#' @param var_method The method to use when computing the variance estimator.  Options are a Taylor linearized technique: "LinHB"= Hajek-Berger estimator, "LinHH" = Hansen-Hurwitz estimator, "LinHTSRS" = Horvitz-Thompson estimator under simple random sampling without replacement, and "LinHT" = Horvitz-Thompson estimator.
#' @param ... Any additional arguments that can be passed to mase::horvitzThompson, mase::greg, and mase::postStrat

#' @examples 
#' library(survey)
#' data(api) 
#' ratio(y_num = apisrs$api.stu, y_den = apisrs$enroll, xsample = apisrs$stype,
#' xpop = apipop$stype, pi = apisrs$pw^(-1), estimator = "postStrat",
#' var_est = TRUE, var_method = "LinHB", datatype = "raw")
#' 
#'@references 
#'\insertRef{coc77}{mase} 
#' 
#'\insertRef{sar92}{mase}
#'
#' @return A list of output containing:
#' \itemize{
#' \item{ratio_est:}{Estimate of the ratio of the population totals/means of the two estimators}
#' \item{variance_est:}{Estimate of the variance of the ratio of two estimators}
#' }
#' @import dplyr
#' @import magrittr
#' @export ratio
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
                  var_method = "LinHB",
                  datatype = "raw",
                  messages = T,
                  ...){
  
  # I don't think we need to add any argument checks here, because the mase 
  # functions that the arguments are fed into will output the correct errors
  # for us. things like pi and N are needed for later variance calculations
  # but they will always be called as arguments to a mase:: function before
  # those steps so if they are in the wrong format those functions should
  # output the correct errors for us
  
  weight <- as.vector(pi^(-1))
  
  
  if (estimator == "horvitzThompson") {
    
    ht_num <- horvitzThompson(y_num, pi, N, pi2, messages, ...)
    ht_den <- horvitzThompson(y_den, pi, N, pi2, messages, ...)
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
      fpc = F,
      var_est = var_est,
      var_method = var_method,
      datatype = datatype,
      messages = messages,
      ...
    )
    
    ps_den <- postStrat(
      y = y_den,
      xsample = xsample,
      xpop = xpop,
      pi = pi,
      N = N,
      pi2 = pi2,
      var_est = var_est,
      var_method = var_method,
      datatype = datatype,
      messages = messages,
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
      ...
    )
    
    greg_den <- greg(
      y = y_den,
      xsample = xsample,
      xpop = xpop,
      pi = pi,
      N = N,
      messages = messages,
      ...
    )
    
    est_den <- greg_den$pop_total
    
    rat <- greg_num$pop_total/greg_den$pop_total
    
    xsample_d_var <- model.matrix(~., data = data.frame(xsample))
    xsample_var <- data.frame(xsample_d_var[,-1, drop = FALSE])
    xsample_dt_var <- t(xsample_d_var) 
    
    y_hat_num <- xsample_d_var %*% solve(xsample_dt_var %*% diag(weight) %*% xsample_d_var) %*% 
      (xsample_dt_var) %*% diag(weight)%*%y_num
    y_hat_den <- xsample_d_var %*% solve(xsample_dt_var %*% diag(weight) %*% xsample_d_var) %*%
      (xsample_dt_var) %*% diag(weight) %*% y_den
    
    e_num <- y_num - y_hat_num
    e_den <- y_den - y_hat_den
    
    e_ratio <- e_num - as.vector(rat)*e_den
    
  }
  
  if (var_est == T) {
    
    var_est <- varMase(
      y = e_ratio,
      pi = pi,
      pi2 = pi2,
      method = var_method,
      N = N
    )
    
    return(list(ratio_est = rat, variance_est = (1/est_den^2)*var_est))
    
  } else {
    
    return(list(ratio_est = e_ratio))
    
  }
  
}
