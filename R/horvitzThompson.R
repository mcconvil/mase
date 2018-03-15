#' Compute the Horvitz-Thompson Estimator
#' 
#' Calculate the Horvitz-Thompson Estimator for a finite population mean/proportion or total based on sample data collected from a complex sampling design.  
#'  
#' @param y  A numeric vector of the sampled response variable.
#' @param pi A numeric vector of inclusion probabilities for each sampled unit in y.  If NULL, then simple random sampling without replacement is assumed.
#' @param N A numeric value of the population size. If NULL, it is estimated with the sum of the inverse of the pis.
#' @param var_est A logical indicating whether or not to compute a variance estimator.  Default is FALSE.
#' @param var_method The method to use when computing the variance estimator.  Options are a Taylor linearized technique: "LinHB"= Hajek-Berger estimator, "LinHH" = Hansen-Hurwitz estimator, "LinHTSRS" = Horvitz-Thompson estimator under simple random sampling without replacement, and "LinHT" = Horvitz-Thompson estimator or a resampling technique: "bootstrapSRS" = bootstrap variance estimator under simple random sampling without replacement.
#' @param pi2 A square matrix of the joint inclusion probabilities.  Needed for the "LinHT" variance estimator.
#' @param B The number of bootstrap samples if computing the bootstrap variance estimator.  Default is 1000.
#' 
#' @examples 
#' library(survey)
#' data(api)
#' horvitzThompson(y = apisrs$api00, pi = apisrs$pw^(-1))
#' horvitzThompson(y = apisrs$api00, pi = apisrs$pw^(-1), var_est = TRUE, var_method = "LinHTSRS")
#'
#'@references{
#'\insertRef{hor52}{mase}
#'}
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
#' @importFrom Rdpack reprompt
#' @include varMase.R
#' @include htt.R
#'
horvitzThompson <- function(y, pi = NULL, N = NULL, pi2 = NULL, var_est =FALSE, var_method="LinHB", B = 1000) {

  ### INPUT VALIDATION ###
  
  #Make sure the var_method is valid
  if(!is.element(var_method, c("LinHB", "LinHH", "LinHTSRS", "LinHT", "bootstrapSRS"))){
    message("Variance method input incorrect. It has to be \"LinHB\", \"LinHH\", \"LinHT\", \"LinHTSRS\", or \"bootstrapSRS\".")
    return(NULL)
  }

  #Check that y is numeric
  if(!(typeof(y) %in% c("numeric", "integer", "double"))){
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
  total <- as.vector(t(y) %*% weight)
  
  # defining estimate for population size if N is unknown, otherwise use
  # known N.
  if (is.null(N)) {
    N <- sum(weight)
    mu <- as.vector(total * (1/N))
  } else {
    mu <- as.vector((total/N))
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
