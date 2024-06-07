#' Compute the Horvitz-Thompson Estimator
#' 
#' Calculate the Horvitz-Thompson Estimator for a finite population mean/proportion or total based on sample data collected from a complex sampling design.  
#'  
#' @param y  A numeric vector of the sampled response variable.
#' @param pi A numeric vector of inclusion probabilities for each sampled unit in y.  If NULL, then simple random sampling without replacement is assumed.
#' @param N A numeric value of the population size. If NULL, it is estimated with the sum of the inverse of the pis.
#' @param var_est A logical indicating whether or not to compute a variance estimator.  Default is FALSE.
#' @param var_method The method to use when computing the variance estimator.  Options are a Taylor linearized technique: "LinHB"= Hajek-Berger estimator, "LinHH" = Hansen-Hurwitz estimator, "LinHTSRS" = Horvitz-Thompson estimator under simple random sampling without replacement, and "LinHT" = Horvitz-Thompson estimator or a resampling technique: "bootstrapSRS" = bootstrap variance estimator under simple random sampling without replacement. The default is "LinHB".
#' @param pi2 A square matrix of the joint inclusion probabilities.  Needed for the "LinHT" variance estimator.
#' @param B The number of bootstrap samples if computing the bootstrap variance estimator.  Default is 1000.
#' @param fpc Default to TRUE, logical for whether or not the variance calculation should include a finite population correction when calculating the "LinHTSRS" or the "SRSbootstrap" variance estimator.
#' @param messages A logical indicating whether to output the messages internal to mase. Default is TRUE.
#' 
#' @examples 
#' library(dplyr)
#' data(IdahoSamp)
#' data(IdahoPop)
#' xsample <- filter(IdahoSamp, COUNTYFIPS == "16055")
#' xpop <- filter(IdahoPop, COUNTYFIPS == "16055") 
#' 
#' horvitzThompson(y = xsample$BA_TPA_ADJ,
#'                 N = xpop$npixels,
#'                 var_est = TRUE,
#'                 var_method = "LinHTSRS")
#'
#' @references
#' \insertRef{hor52}{mase}
#' 
#' @return 
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
#' @export horvitzThompson
#' @import boot
#' @importFrom Rdpack reprompt
#' @include varMase.R
#' @include bootHelpers.R
#'
horvitzThompson <- function(y,
                            pi = NULL,
                            N = NULL,
                            pi2 = NULL,
                            var_est = FALSE,
                            var_method = "LinHB",
                            B = 1000,
                            fpc = TRUE,
                            messages = TRUE) {

  if (!is.element(var_method, c("LinHB", "LinHH", "LinHTSRS", "LinHT", "bootstrapSRS"))) {
    stop("Variance method input incorrect. It has to be \"LinHB\", \"LinHH\", \"LinHT\", \"LinHTSRS\", or \"bootstrapSRS\".")
  }

  if (!(typeof(y) %in% c("numeric", "integer", "double"))) {
    stop("Must supply numeric y.  For binary variable, convert to 0/1's.")
  }
  
  if (is.null(pi) && is.null(N)) {
    stop("Must supply either ", sQuote("pi")," or ", sQuote("N"))
    
  }
  
  if (is.null(pi)) {
    if (messages) {
      message("Assuming simple random sampling") 
    }
  }  
  
  if (is.null(pi)) {
    pi <- rep(length(y)/N, length(y))
  }
  
  weight <- pi^(-1)
  
  #Sample size
  n <- length(y)
  
  # total population
  total <- as.vector(t(y) %*% weight)
  
  if (is.null(N)) {
    N <- sum(weight)
    mu <- as.vector(total * (1/N))
  } else {
    mu <- as.vector((total/N))
  }
  
  if (var_est == TRUE) {
    
    if (var_method != "bootstrapSRS") {
      varEst <- varMase(y = y, pi = pi,pi2 = pi2,method = var_method, N = N)
    } else {
      #Find bootstrap variance
      dat <- cbind(y,pi)
      #Bootstrap total estimates
      t_boot <- boot(data = dat, statistic = htt, R = B)
      
      if (fpc == TRUE) {
        varEst <- var(t_boot$t)*n/(n-1)*(N-n)/(N-1) 
      }
      if (fpc == FALSE) {
        varEst <- var(t_boot$t)*n/(n-1)
      }
      
    }
    
    varEstMu <- varEst*N^(-2)
    
    return(list(pop_total = total,
                pop_mean = mu, 
                pop_total_var = varEst,
                pop_mean_var = varEstMu))
    
  } else {
    return(list(pop_total = total,
                pop_mean = mu))
  }
}
