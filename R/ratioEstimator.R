#' Compute a ratio estimator
#'
#' Calculates a ratio estimator for a finite population mean/proportion or total based on sample data collected from a complex sampling design and auxiliary population data.  
#' 
#' @inheritParams horvitzThompson
#' @param xsample A numeric vector of the sampled auxiliary variable.
#' @param xpop A numeric vector of population level auxiliary information.  Must come in the form of raw data, population total or population mean.
#' @param datatype A string that specifies the form of population auxiliary data. The possible values are "raw", "total" or "mean".  If datatype = "raw", then xpop must contain a numeric vector of the auxiliary variable for each unit in the population. If datatype = "total" or "mean", then contains either the population total or population mean for the auxiliary variable.
#' 
#' @examples 
#' library(dplyr)
#' data(IdahoPop)
#' data(IdahoSamp)
#' 
#' xsample <- filter(IdahoSamp, COUNTYFIPS == "16055")
#' xpop <- filter(IdahoPop, COUNTYFIPS == "16055")
#' 
#' ratioEstimator(y = xsample$BA_TPA_ADJ,
#'                xsample = xsample$tcc,
#'                xpop = xpop$tcc,
#'                datatype = "means",
#'                N = xpop$npixels)
#' 
#' @references 
#' \insertRef{coc77}{mase} 
#' \insertRef{sar92}{mase}
#'
#' @returns 
#' List of output containing:
#' 
#' * pop_total: Estimate of population total.
#' 
#' * pop_mean: Estimate of population mean.
#' 
#' * pop_total_var: Estimated variance of population total estimate.
#' 
#' * pop_mean_var: Estimated variance of population mean estimate.
#' 
#' 
#' @export ratioEstimator
#' @import boot
#' @include bootHelpers.R
#' @include varMase.R


ratioEstimator <- function(y,
                           xsample,
                           xpop,
                           datatype = "raw",
                           pi = NULL, 
                           N = NULL,
                           pi2 = NULL,
                           var_est = FALSE,
                           var_method = "LinHB",
                           B = 1000,
                           fpc = TRUE,
                           messages = TRUE) {

  if (!(typeof(y) %in% c("numeric", "integer", "double"))) {
    stop("Must supply numeric y.  For binary variable, convert to 0/1's.")
  }
  
  #Make sure the var_method is valid
  if (!is.element(var_method, c("LinHB", "LinHH", "LinHTSRS", "LinHT", "bootstrapSRS"))) {
    stop("Variance method input incorrect. It has to be \"LinHB\", \"LinHH\", \"LinHT\", \"LinHTSRS\", or \"bootstrapSRS\".")
  }
  
  if (is.null(pi)) {
    if (messages) {
      message("Assuming simple random sampling") 
    }
  }  
  
  #Determine N if not provided
  if (is.null(N)) {
    if (datatype == "raw") {
      N <- dim(data.frame(xpop))[1]
    }else {
      N <- sum(pi^(-1))
      if (messages) {
        message("Assume N can be approximated by the sum of the inverse inclusion probabilities.") 
      }
    }
  }
  
  
  # create equal SRS weights if not provided
  if (is.null(pi)) {
    pi <- rep(length(y)/N, length(y))
  }
  
  #calculate Horvitz-Thompson estimator for y and x
  tyHT <- horvitzThompson(y=y,pi=pi)$pop_total
  txHT <- horvitzThompson(y=xsample,pi=pi)$pop_total
  r <- as.vector(tyHT/txHT)
  
  if (datatype == "raw") {
    tau_x <- sum(xpop)
  }
  if (datatype %in% c("total", "totals")) {
    tau_x <- xpop
  }
  if (datatype %in% c("mean", "means")) {
    tau_x <- xpop*N
  }
  
  # Mean of x
  mu_x <- tau_x/N
  
  # calculate estimates
  mu_r <- r * mu_x
  tau_r <- r * tau_x
  
  # Estimate the variance
  if (var_est == TRUE) {
    
    if (var_method != "bootstrapSRS") {
      y.hat <- r*as.vector(xsample)
      varEst <- varMase(y = (y-y.hat),
                        pi = pi,
                        pi2 = pi2,
                        method = var_method,
                        N = N,
                        fpc = fpc)
      varEstMu <- varEst*N^(-2)
    }
    
    if(var_method =="bootstrapSRS") {

      #Find bootstrap variance
      dat <- cbind(y,pi, xsample)
      #Bootstrap total estimates
      t_boot <- boot(data = dat, statistic = ratioEstimatort, R = B, tau_x = tau_x)
      
      #Adjust for bias and without replacement sampling
      n <- length(y)
      
      if (fpc == T) {
        varEst <- var(t_boot$t)*n/(n-1)*(N-n)/(N-1)
      }
      if (fpc == F) {
        varEst <- var(t_boot$t)*n/(n-1)
      }
      
      varEstMu <- varEst*N^(-2)
    }
    
    return(list(pop_total = tau_r,
                pop_mean = mu_r,
                pop_total_var = varEst,
                pop_mean_var = varEstMu))
    
  } else {
    return(list(pop_total = tau_r,
                pop_mean = mu_r))
  }
  
}
  
