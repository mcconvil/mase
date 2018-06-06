#' Compute a ratio estimator
#'
#' Calculates a ratio estimator for a finite population mean/proportion or total based on sample data collected from a complex sampling design and auxiliary population data.  
#' 
#' @inheritParams horvitzThompson
#' @param x_sample A numeric vector of the sampled auxiliary variable.
#' @param x_pop A numeric vector of population level auxiliary information.  Must come in the form of raw data, population total or population mean.
#' @param data_type A string that specifies the form of population auxiliary data. The possible values are "raw", "total" or "mean".  If data_type = "raw", then x_pop must contain a numeric vector of the auxiliary variable for each unit in the population. If data_type = "total" or "mean", then contains either the population total or population mean for the auxiliary variable.
#' 
#' @examples 
#' library(survey)
#' data(api)
#'ratioEstimator(y = apisrs$api00, x_sample = apisrs$meals, 
#'x_pop = sum(apipop$meals), data_type = "total", pi = apisrs$pw^(-1), 
#'N = dim(apipop)[1])
#' 
#'@references 
#'\insertRef{coc77}{mase} 
#'\insertRef{sar92}{mase}
#'
#' @return A list of output containing:
#' \itemize{
#' \item{pop_total: Estimate of population total}
#' \item{pop_mean: Estimate of the population mean}
#' \item{pop_total_var: Estimated variance of population total estimate}
#' \item{pop_mean_var: Estimated variance of population mean estimate}
#' \item{weights: Survey weights produced by ratio estimator}
#' }
#' 
#' @export ratioEstimator
#' @import boot
#' @include ratioEstimatort.R
#' @include varMase.R
#' 
#' @seealso \code{\link{greg}} for a linear or logistic regression model.


ratioEstimator <- function(
  y, x_sample, x_pop, data_type = "raw", pi = NULL, N = NULL, pi2 = NULL, var_est = FALSE, var_method = "lin_HB", B = 1000, strata = NULL) {


    ### INPUT VALIDATION ###
  
  #Check that y is numeric
  if(!(typeof(y) %in% c("numeric", "integer", "double"))){
    stop("Must supply numeric y.  For binary variable, convert to 0/1's.")
  }
  
  #Make sure the var_method is valid
  if(!is.element(var_method, c("lin_HB", "lin_HH", "lin_HTSRS", "lin_HT", "bootstrap_SRS"))){
    message("Variance method input incorrect. It has to be \"lin_HB\", \"lin_HH\", \"lin_HT\", \"lin_HTSRS\", or \"bootstrap_SRS\".")
    return(NULL)
  }
  
  #Need to provide either data_type="raw", N, or pi.  Give warning if not
  if(data_type %in% c("means", "totals") & is.null(N) & is.null(pi)){
    message("Must supply N, pi, or raw population data so that we can estimate N.")
    return(NULL)
  }
  
  
  #Check on inclusion probabilities
  if(is.null(pi)){
    message("Assuming simple random sampling")
  }  
  
  #Determine N if not provided
  if (is.null(N)) {
    if( data_type == "raw") {
      N <- dim(data.frame(x_pop))[1]
    }else {
      N <- sum(pi^(-1))
      warning("Assuming N can be approximated by the sum of the inverse inclusion probabilities.")
    }
  }
  
  #Check length of x_pop
  if(data_type == "raw" & length(x_pop) != N){
    message("Population size, N, and length of x_pop must match.")
    return(NULL)
  }
  if(data_type %in% c("mean", "total") & length(x_pop) != 1){
    message("If data_type is \"mean\" or \"total\", then x_pop must have length 1 since it is the corresponding population mean or total for x.")
    return(NULL)
  }
  
  
  # create equal SRS pis if not given
  if (is.null(pi)) {
    pi <- rep(length(y)/N, length(y))
  }
  
  
  #calculate Horvitz-Thompson estimator for y and x
  tyHT <- horvitzThompson(y=y,pi=pi)$pop_total
  txHT <- horvitzThompson(y=x_sample,pi=pi)$pop_total
  r <- as.vector(tyHT/txHT)
  

  #Find population total of x
  if (data_type == "raw"){
    tau_x <- sum(x_pop)
  }
  if (data_type %in% c("total", "totals")){
    tau_x <- x_pop
  }
  if (data_type %in% c("mean", "means")){
    tau_x <- x_pop*N
  }
  
 #Mean of x
  mu_x <- tau_x/N
  
  
  #calculate estimates
  mu_r <- r * mu_x
  tau_r <- r * tau_x
  
  #weights
  weights <- tau_x/txHT/pi
  
  #Estimate the variance
  if(var_est==TRUE){
    
    if(var_method!="bootstrap_SRS"){
    y_hat <- r*as.vector(x_sample)
    varEst <- varMase(y = (y-y_hat),pi = pi,pi2 = pi2,method = var_method, N = N, strata = strata)
    varEstMu <- varEst*N^(-2)
    }
    
    if(var_method ==  "bootstrap_SRS") {

      #Find bootstrap variance
      dat <- cbind(y,pi, x_sample)
      #Bootstrap total estimates
      t_boot <- boot(data = dat, statistic = ratioEstimatort, R = B, tau_x = tau_x)
      
      #Adjust for bias and without replacement sampling
      n <- length(y)
      varEst <- var(t_boot$t)*n/(n-1)*(N-n)/(N-1)
      varEstMu <- varEst*N^(-2)
    }
    
    return(list(pop_total = tau_r, 
                pop_mean = mu_r, 
                pop_total_var=varEst, 
                pop_mean_var = varEstMu, 
                weights = weights))
    
  }else{
    
    return(list(pop_total = tau_r, 
                pop_mean = mu_r, 
                weights = weights))
  }
  
}
  
