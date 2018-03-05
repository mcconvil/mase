#' Compute the Horvitz-Thompson Estimator
#' 
#' Calculate the Horvitz-Thompson Estimator for a finite population mean/proportion or total based on sample data collected from a complex sampling design.  
#'  
#' @param y  A numeric vector or a matrix with one column of the sampled response variable. 
#' @param pi Default assumes equal probability/simple random sampling. If unequal probability, requires vector of first-order inclusion probabilities of same length as y
#' @param N population size, if not provided estimate with the sum of the inverse inclusion probabilities
#' @param var_est Default is FALSE, logical indicating whether or not to compute a variance estimator
#' @param var_method Method to use when computing the variance estimator.  Options are "HB"= Hajek-Berger estimator, "HH" = Hansen-Hurwitz estimator, "HTSRS" = Horvitz-Thompson estimator under simple random sampling without replacement, "HT" = Horvitz-Thompson estimator, "bootstrapSRS" = bootstrap variance estimator under simple random sampling without replacement
#' @param pi2 = a square matrix of the joint inclusion probabilities.  Needed for the "HT" variance estimator
#' @param B number of bootstrap samples if computing the bootstrap variance estimator.  Default is 1000.
#'
#'@references 
#'\insertRef{hor52}{mase}
#'
#' @return List of output containing:
#' \itemize{
#' \item{pop_total:}{Estimate of population total}
#' \item{pop_mean:}{Estimate of population mean}
#' \item{pop_total_var:}{ Estimated variance of population total estimate}
#' \item{pop_mean_var:}{ Estimated variance of population mean estimate}
#' }
#'
#' @export horvitzThompson
#' @import boot
#' @include varMase.R
#' @include htt.R
#'
horvitzThompson <- function(y, pi = NULL, N = NULL, pi2 = NULL, var_est =FALSE, var_method="HB", B = 1000) {

  ### INPUT VALIDATION ###
  
  #Make sure the var_method is valid
  if(!is.element(var_method, c("HB", "HH", "HTSRS", "HT", "bootstrapSRS"))){
    message("Variance method input incorrect. It has to be \"HB\", \"HH\", \"HT\", \"HTSRS\", or \"bootstrapSRS\".")
    return(NULL)
  }

  #Check that y is numeric
  if(!(typeof(y) %in% c("numeric", "integer"))){
    stop("Must supply numeric y.  For binary variable, convert to 0/1's.")
  }
  
  
  
  if(is.null(pi) && is.null(N)){
    stop("Must supply either ", sQuote("pi")," or ", sQuote("N"))
    
  }
  
  if(is.null(pi)){
    message("Assuming simple random sampling")
  }  
  
    
  # convert pi into diagonal matrix format
  if (is.null(pi)) {
    pi <- rep(length(y)/N, length(y))
  }
  
  #weight: inverse first order inclusion probabilities
  weight <- pi^(-1)
  
  #Sample size
  n <- length(y)
  
  
##########################
  
  # total population
  total <- t(y) %*% weight
  
  # defining estimate for population size if N is unknown, otherwise use
  # known N.
  if (is.null(N)) {
    N <- sum(weight)
    mu <- total * (1/N)
  } else {
    mu <- (total/N)
  }
  
  if(var_est==TRUE){
    if(var_method!="bootstrapSRS"){
    varEst <- varMase(y = y,pi = pi,pi2 = pi2,method = var_method, N = N)

    }else{
      #Find bootstrap variance
      dat <- cbind(y,pi)
      #Bootstrap total estimates
      t_boot <- boot(data = dat, statistic = htt, R = B)
      
      #Adjust for bias and without replacement sampling
      varEst <- var(t_boot$t)*n/(n-1)*(N-n)/(N-1)
    }
    
    #return estimates and variance estimates
    varEstMu <- varEst*N^(-2)
    return(list(pop_total = total, pop_mean = mu, pop_total_var=varEst, pop_mean_var = varEstMu))
    
  }else{
  
  # return estimates
  return(list(pop_total = total, pop_mean = mu))
  }
}
