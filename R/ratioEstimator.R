#' Compute a ratio estimator
#'
#' Calculates a ratio estimator for a finite population mean/proportion or total based on sample data collected from a complex sampling design and auxiliary population data.  
#' 
#' @param y vector of response variable of length n, or a matrix with 
#' dimension n * 1
#' @param xsample predictor vector of sampled data of length n
#' @param xpop Dataframe of population level auxiliary information.  Must come in the form of raw data, population totals or population means.
#' @param datatype Default to "raw", takes values "raw", "totals", and "means" for whether the user is providing the raw population values for x, the population totals for x, or the population means for x
#' @param pi Default to assume equal probability/simple random sampling, if unequal probability, requires vector of first-order inclusion probabilities of same length as y
#' @param N population size, if not provided estimate with the sum of the inverse inclusion probabilities
#' @param var_est Default to FALSE, logical for whether or not to compute estimate of variance
#' @param var_method Method to use when computing the variance estimate.  Options are "HB"= Hajek-Berger estimator, "HH" = Hansen-Hurwitz estimator, "HTSRS" = Horvitz-Thompson estimator under simple random sampling, "HT" = Horvitz-Thompson estimator, "bootstrapSRS" = bootstrap variance estimator under simple random sampling without replacement
#' @param pi2 = a n * n matrix of the joint inclusion probabilities.  Needed for the "HT" variance estimator
#' @param B number of bootstrap samples if computing the bootstrap variance estimator.  Default is 1000.
#'
#'@references 
#'\insertRef{coc77}{mase} 
#'\insertRef{sar92}{mase}
#'
#' @return List of output containing:
#' \itemize{
#' \item{pop_total:}{Estimate of population total}
#' \item{pop_mean:}{Estimate of population mean}
#' \item{pop_total_var:}{ Estimated variance of population total estimate}
#' \item{pop_mean_var:}{ Estimated variance of population mean estimate}
#' }
#' 
#' @export ratioEstimator
#' @import boot
#' @include ratioEstimatort.R
#' @include varMase.R


ratioEstimator <- function(
  y, xsample, xpop, datatype = "raw", pi = NULL, N = NULL, pi2 = NULL, var_est = FALSE, var_method = "HB", B = 1000) {


    ### INPUT VALIDATION ###
  
  #Make sure the var_method is valid
  if(!is.element(var_method, c("HB", "HH", "HTSRS", "HT", "bootstrapSRS"))){
    message("Variance method input incorrect. It has to be \"HB\", \"HH\", \"HT\", \"HTSRS\", or \"bootstrapSRS\".")
    return(NULL)
  }
  
  
  #Check on inclusion probabilities and create weight=inverse inclusion probabilities
  if(is.null(pi)){
    message("Assuming simple random sampling")
  }  
  
  #Determine N if not provided
  if (is.null(N)) {
    if( datatype == "raw") {
      N <- dim(data.frame(xpop))[1]
    }else {
      N <- sum(pi^(-1))
      warning("Assume N can be approximated by the sum of the inverse inclusion probabilities.")
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
  

  #Find population total of x
  if (datatype == "raw"){
    tau_x <- sum(xpop)
  }
  if (datatype == "totals"){
    tau_x <- xpop
  }
  if (datatype == "means"){
    tau_x <- xpop*N
  }
  
 #Mean of x
  mu_x <- tau_x/N
  
  
  #calculate estimates
  mu_r <- r * mu_x
  tau_r <- r * tau_x
  
  #Estimate the variance
  if(var_est==TRUE){
    
    if(var_method!="bootstrapSRS"){
    y.hat <- r*as.vector(xsample)
    varEst <- varMase(y = (y-y.hat),pi = pi,pi2 = pi2,method = var_method, N = N)
    varEstMu <- varEst*N^(-2)
    }
    
    if(var_method ==  "bootstrapSRS") {

      #Find bootstrap variance
      dat <- cbind(y,pi, xsample)
      #Bootstrap total estimates
      t_boot <- boot(data = dat, statistic = ratioEstimatort, R = B, tau_x = tau_x)
      
      #Adjust for bias and without replacement sampling
      n <- length(y)
      varEst <- var(t_boot$t)*n/(n-1)*(N-n)/(N-1)
      varEstMu <- varEst*N^(-2)
    }
    
    return(list(pop_total = tau_r, pop_mean = mu_r, pop_total_var=varEst, pop_mean_var = varEstMu))
    
  }else{
    
    return(list(pop_total = tau_r, pop_mean = mu_r))
  }
  
}
  
