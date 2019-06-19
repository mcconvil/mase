#' 
#' Compute a generalized regression estimator
#' 
#' Summarizes population estimates given precomputed predicted sample and population values.
#'  
#' @inheritParams horvitzThompson
#' @param data_type A string that specifies the form of population auxiliary data. The possible values are "raw", "totals" or "means" for whether the user is providing population estimates at the unit level, aggregated to totals, or aggregated to means.  Default is "raw".
#' @param y_hat_sample A vector of predicted sample estimates from a survey model.
#' @param y_hat_pop If data_type is "raw", y_hat_pop must be a vector of predicted population estimates from a survey model. Otherwise it must be either the total or mean population estimate. 
#' 
#' 
#' @return A greg object containing:
#' \itemize{
#' \item{pop_total: Estimate of population total}
#' \item{pop_mean: Estimate of the population mean}
#' \item{pop_total_var: Estimated variance of population total estimate}
#' \item{pop_mean_var: Estimated variance of population mean estimate}
#' \item{weights: Survey weights produced by greg (linear model only)}
#' \item{coefficients: Survey-weighted model coefficients}
#' }
#' 
#' @export gregEstimate
#' @import survey
#' @import glmnet
#' @import boot
#' @importFrom stats model.matrix predict quasibinomial var
#' @include varMase.R
#' @include gregt.R
#' @include gregObject.R
#' 
gregEstimate <- function(y, y_hat_sample, pi = NULL, y_hat_pop, pi2 = NULL,
                         N = NULL, data_type = "raw", var_est =FALSE, var_method="lin_HB",
                         B = 1000, strata = NULL) {
  
  
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
  
  if(!is.element(model, c("linear","logistic"))){
    message("Method input incorrect, has to be either \"linear\" or \"logistic\"")
    return(NULL)
  }
  
  #Need to provide either data_type="raw", N, or pi.  Give warning if not
  if(data_type %in% c("means", "totals") & is.null(N) & is.null(pi)){
    message("Must supply N, pi, or raw population data so that we can estimate N.")
    return(NULL)
  }

  #Format y
  y <- as.vector(y)
  y_hat_sample <- as.vector(y_hat_sample)
  y_hat_pop <- as.vector(y_hat_pop)
  n <- length(y)
  
  
  #Make sure y's are complete
  if(NA %in% y){
    message("Must supply complete cases for y.")
    return(NULL)
  }
  if(NA %in% y_hat_sample){
    message("Must supply complete cases for y_hat_sample.")
    return(NULL)
  }
  if(NA %in% y_hat_pop){
    message("Must supply complete cases for y_hat_pop.")
    return(NULL)
  }
  #Check that the number of observations are consistent
  if(length(y_hat_sample) != length(y)){
    message("y and y_hat_sample must have same number of observations.")
    return(NULL)
  }
  #Need to get N if not provided
  if(is.null(N)){
    if(data_type=="raw"){
      N <- length(y_hat_pop)
    }else{
      N <- sum(pi^(-1))
      message("Assuming N can be approximated by the sum of the inverse inclusion probabilities.")
    }
  }
  #Check on inclusion probabilities and create weight=inverse inclusion probabilities
  if(is.null(pi)){
    message("Assuming simple random sampling")
    pi <- rep(length(y)/N, length(y))
  }  
  
  #weight: inverse first order inclusion probabilities
  weight <- as.vector(pi^(-1))
  #Format x_pop to be a vector of pop totals
  if (data_type=="totals"){
    #Make sure to only take the values which are also in x_sample
    if (length(y_hat_pop) != 1){
      message("y_hat_pop must be a single value if data_type == 'totals'")
      return(NULL)
    }
  }
  if (data_type=="means"){
    if (length(y_hat_pop) != 1){
      message("y_hat_pop must be a single value if data_type == 'means'")
      return(NULL)
    }
    y_hat_pop <- y_hat_pop * N
  }
    
  #Total estimate
  t <- sum((y-y_hat_sample)*weight) + sum(y_hat_pop)
  
  if ( var_est == TRUE ) {
    if ( var_method != "bootstrap_SRS") {
      varEst <- varMase(y = (y-y_hats_s), pi = pi, pi2 = pi2, method = var_method, N = N, strata = strata)
      }
    if ( var_method == "bootstrap_SRS"){
      #Find bootstrap variance
      #Sample data
      dat <- cbind(y,pi, x_sample_d)
      #Bootstrap total estimates
      t_boot <- boot(data = dat, statistic = gregElasticNett, R = B, x_pop_d = x_pop_d, alpha=alpha, lambda = lambda_opt, parallel = "multicore", ncpus = 2)
      #Adjust for bias and without replacement sampling
      varEst <- var(t_boot$t)*n/(n-1)*(N-n)/(N-1)
    }
  }
  if(var_est==TRUE){
    return(list( pop_total = as.numeric(t), 
                 pop_mean = as.numeric(t)/N, 
                 pop_total_var=varEst, 
                 pop_mean_var=varEst/N^2,
                 y_hat_sample = y_hats_sample,
                 y_hat_pop = y_hat_pop %>% as.vector()) %>%
             gregify())
  }
  else{
    return(list( pop_total = as.numeric(t), 
                 pop_mean = as.numeric(t)/N, 
                 y_hat_sample = y_hats_sample,
                 y_hat_pop = y_hat_pop %>% as.vector()) %>%
             gregify())
    
  }
}
  