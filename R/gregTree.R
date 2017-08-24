#' Compute a regression tree estimator
#' 
#' Calculates a regression tree estimator for a finite population mean/proportion or total based on sample data collected from a complex sampling design and auxiliary population data.  
#' 
#' @param y  vector or a matrix with one column of the sampled response variable
#' @param pi Default to assume equal probability/simple random sampling, if unequal probability, requires vector of first-order inclusion probabilities of same length as y
#' @param xsample matrix of auxiliary data in the sample.  Number of rows should match length of y.
#' @param xpop Dataframe of population level auxiliary information.
#' @param var_est Default to FALSE, logical for whether or not to compute estimate of variance
#' @param var_method Method to use when computing the variance estimate.  Options are "HB"= Hajek-Berger estimator, "HH" = Hansen-Hurwitz estimator, "HTSRS" = Horvitz-Thompson estimator under simple random sampling, "HT" = Horvitz-Thompson estimator, "bootstrapSRS" = bootstrap variance estimator under simple random sampling without replacement
#' @param pi2 a square matrix of the joint inclusion probabilities.  Needed for the "HT" variance estimator
#' @param N population size, if not provided estimated to be the sum of the inverse inclusion probabilities
#' @param B number of bootstrap samples if computing the bootstrap variance estimator.  Default is 1000.
#' 
#'@references 
#'\insertRef{mcc17b}{mase}

#' 
#' @return A list of output containing:
#' \itemize{
#' \item{pop_total:}{Estimate of population total}
#' \item{pop_mean:}{Estimate of the population mean (or proportion)}
#' \item{weights.greg:}{Survey weights produced by GREG (linear model only)}
#' \item{pop_total_var:}{Estimated variance of population total estimate}
#' \item{pop_mean_var:}{Estimated variance of population mean estimate}
#' }
#' 
#' @export gregTree
#' @import rpms
#' @import boot
#' @include varMase.R
#' @include gregt.R

# gregTree  <- function(y, xsample, xpop, pi = NULL,  pi2 = NULL, var_est = FALSE, var_method="HB", N = NULL, B = 1000){
# 
#   
#   
# ### INPUT VALIDATION ###
#   
#   
#   #Need to get N if not provided
#   if(is.null(N)){
#       N <- sum(pi^(-1))
#       message("Assuming N can be approximated by the sum of the inverse inclusion probabilities.")
#     }
# 
#   
#   #Convert y to a vector
#   y <- as.vector(y)
#   
#   #sample size
#   n <- length(y)
# 
#   
#   #Check on inclusion probabilities and create weight=inverse inclusion probabilities
#   if(is.null(pi)){
#     message("Assuming simple random sampling")
#   }  
#   
#   # convert pi into diagonal matrix format
#   if (is.null(pi)) {
#     pi <- rep(length(y)/N, length(y))
#   }
#   
#   #weight: inverse first order inclusion probabilities
#   weight <- as.vector(pi^(-1))
#   
#   #Fit model
#   
#   #Calculate weights  
#   # w <- 
#   
#   #calculating the total estimate for y
#   t <- w %*% y
#   
# 
#   #NOTE: check that weights times x's should equal total of x's to check for correct weight values
#   
#   #Coefficients
#   coefs <- solve(xsample.dt %*% diag(weight) %*% xsample.d) %*% (xsample.dt) %*% diag(weight) %*% y
#   
#   
#   if(var_est==TRUE){
#     if(var_method!="bootstrapSRS"){
#     y.hat <- xsample.d%*%solve(xsample.dt %*% diag(weight) %*% xsample.d) %*% (xsample.dt) %*% diag(weight)%*%y
#     e <- y-y.hat
#     varEst <- varMase(y = e,pi = pi,pi2 = pi2,method = var_method, N = N)
#     
#     }else if(var_method=="bootstrapSRS"){
#       #Find bootstrap variance
#       dat <- cbind(y,pi, xsample.d)
#       #Bootstrap total estimates
#       t_boot <- boot(data = dat, statistic = gregt, R = B, xpopd = xpop_d, parallel = "multicore", ncpus = 2)
#       
#       #Adjust for bias and without replacement sampling
#       varEst <- var(t_boot$t)*n/(n-1)*(N-n)/(N-1)
#     }
#     return(list( pop_total = as.numeric(t), 
#                  pop_mean = as.numeric(t)/N,
#                  pop_total_var=varEst, 
#                  pop_mean_var=varEst/N^2, 
#                  weights = as.vector(w),
#                  coefficients =  coefs))
#   }else{
#     return(list( pop_total = as.numeric(t), 
#                  pop_mean = as.numeric(t)/N,           
#                  weights = as.vector(w),
#                  coefficients =  coefs))      
#     
#   }
#   
#   }
#   
# }
# 
