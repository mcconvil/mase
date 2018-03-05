#' Compute an elastic net regression estimator
#' 
#' Calculates a lasso, ridge or elastic net generalized regression estimator for a finite population mean/proportion or total based on sample data collected from a complex sampling design and auxiliary population data.  
#' 
#' 
#' @param y vector of response variable of length n, or a matrix with 
#' dimension n * 1
#' @param xsample vectors of observed data of length n
#' @param xpop Dataframe of population level auxiliary information.  Must come in the form of raw data, population totals or population means.  
#' @param pi Default to assume equal probability/simple random sampling, if unequal probability, requires vector of first-order inclusion probabilities of same length as y
#' @param alpha mixing parameter for the lasso and ridge penalties in the elastic net.  When alpha = 1, uses only a lasso penalty.  When alpha=0, uses only a ridge penalty.
#' @param model regression model to utilize. User must choose 'linear' or 'logistic'
#' @param var_est Default to FALSE, logical for whether or not to compute estimate of variance
#' @param var_method Method to use when computing the variance estimate.  Options are "HB"= Hajek-Berger estimator, "HH" = Hansen-Hurwitz estimator, "HTSRS" = Horvitz-Thompson estimator under simple random sampling, "HT" = Horvitz-Thompson estimator, "bootstrapSRS" = bootstrap variance estimator under simple random sampling without replacement
#' @param pi2 = a n * n matrix of the joint inclusion probabilities.  Needed for the "HT" variance estimator
#' @param datatype Default to "raw", takes values "raw", "totals", and "means" for whether the user is providing the raw population values for x, the population totals for x, or the population means for x
#' @param N population size, if not provided estimated to be the sum of the inverse inclusion probabilities
#' @param lambda Default to "lambda.min", takes values "lambda.min", which is the lambda value associated with the minimum cross validation error or "lambda.1se", which is the lamabda value which is one standard error away from the minimizing lambda and produces a sparser fit
#' @param B number of bootstrap samples if computing the bootstrap variance estimator.  Default is 1000.
#' @param cvfolds number of folds for cross-validation to pick lambda.
#' 
#' @references 
#'\insertRef{mcc17}{mase}

#'
#' @return A list of output containing:
#' \itemize{
#' \item{pop_total:}{Estimate of population total}
#' \item{coefficients:} {Elastic net coefficient estimates}
#' \item{pop_mean:}{Estimate of the population mean (or proportion)}
#' \item{pop_total_var:}{Estimated variance of population total estimate}
#' \item{pop_mean_var:}{Estimated variance of population mean estimate}
#' } 
#' @import boot
#' @import glmnet
#' @importFrom stats model.matrix predict quasibinomial var
#' @export gregElasticNet
#' @include varMase.R
#' @include gregElasticNett.R

gregElasticNet  <- function(
  y, xsample, xpop, pi = NULL, alpha=1, model="linear", pi2 = NULL, var_est =FALSE, var_method="HB", datatype = "raw", N = NULL, lambda = "lambda.min", B = 1000, cvfolds = 10){
  
  
  ### INPUT VALIDATION ###
  #Check that y is numeric
  if(!(typeof(y) %in% c("numeric", "integer", "double"))){
    stop("Must supply numeric y.  For binary variable, convert to 0/1's.")
  }
  
  #Make sure the var_method is valid
  if(!is.element(var_method, c("HB", "HH", "HTSRS", "HT", "bootstrapSRS"))){
    message("Variance method input incorrect. It has to be \"HB\", \"HH\", \"HT\", \"HTSRS\", or \"bootstrapSRS\".")
    return(NULL)
  }
  
  if(!is.element(model, c("linear","logistic"))){
    message("Method input incorrect, has to be either \"linear\" or \"logistic\"")
    return(NULL)
  }
  

  #Need to get N if not provided
  if(is.null(N)){
    if(datatype=="raw"){
      N <- dim(as.matrix(xpop))[1]
    }else{
      N <- sum(pi^(-1))
      message("Assuming N can be approximated by the sum of the inverse inclusion probabilities.")
    }
  }
  
  
  #create design matrix, x matrix and transpose design matrix
  xsample.d <- model.matrix(~., data = data.frame(xsample))
  xsample <- data.frame(xsample.d[,-1])
  xsample.dt <- t(xsample.d) 
  
  #Format y
  y <- as.vector(y)
  n <- length(y)
  
  #Check on inclusion probabilities and create weight=inverse inclusion probabilities
  if(is.null(pi)){
    message("Assuming simple random sampling")
  }  
  
  
  # convert pi into diagonal matrix format
  if (is.null(pi)) {
    pi <- rep(length(y)/N, length(y))
  }
  
  #weight: inverse first order inclusion probabilities
  weight <- as.vector(pi^(-1))
  
  #Cross-validation to find lambdas
  if(model=="linear"){
    fam <- "gaussian"
  } else{
    fam <- "binomial"
  } 
  
  cv <- cv.glmnet(x = as.matrix(xsample), y = y, alpha = alpha, weights = weight, nfolds = cvfolds,family=fam, standardize=FALSE)
  
  if(lambda=="lambda.min"){
    lambda_opt <- cv$lambda.min
  }
  if(lambda=="lambda.1se"){
    lambda_opt <- cv$lambda.1se
  }
  
  
  ## MODEL SELECTION COEFFICIENTS ##
  pred.mod <- glmnet(x = as.matrix(xsample), y = y, alpha = alpha, family=fam, standardize = FALSE, weights=weight)
  elasticNet_coef <- predict(pred.mod,type = "coefficients",s = lambda_opt)[1:dim(xsample.d)[2],]
  
  #Estimated y values in sample
  y.hats.s <- predict(cv,newx = as.matrix(xsample), s = lambda_opt, type="response")

if (model == "logistic") {
  if (datatype != "raw"){
    message("For the Logistic Elastic Net Estimator, user must supply all x values for population.  Populations totals or means for x are not enough.")
    return(NULL)
  }
  
  #Population matrix
  xpop <- data.frame(model.matrix(~., data = xpop))[,-1]
  #Make sure to only take the columns which are also in xsample
  xpop <- dplyr::select_(xpop, .dots=names(xsample))
  xpop_d <- model.matrix(~., data = xpop)
  
  #Total estimate
  y.hats.U <- predict(cv,newx = xpop_d[,-1], s = lambda_opt, type = "response")
  t <- sum(y.hats.U) + t(y-y.hats.s)%*%pi^(-1)
  
  if ( var_est == TRUE){
    if (var_method != "bootstrapSRS") {
      varEst <- varMase(y = (y - y.hats.s),pi = pi,pi2 = pi2,method = var_method, N = N)
      
    }
    
    if(var_method == "bootstrapSRS"){
      #FILL IN: logistic, raw data!
      
      #Sample data
      dat <- cbind(y,pi, xsample.d)
      #Bootstrap total estimates
      t_boot <- boot(data = dat, statistic = logisticGregElasticNett, R = B, xpopd = xpop_d, alpha=alpha, lambda=lambda_opt, parallel = "multicore", ncpus = 2)
      
      #Adjust for bias and without replacement sampling
      varEst <- var(t_boot$t)*n/(n-1)*(N-n)/(N-1)
    } 
    
  }

  
}

if (model == "linear") {
  
  #Format xpop to be a vector of pop totals
  if (datatype=="raw"){
    xpop <- data.frame(model.matrix(~., data = xpop))[,-1]
    #Make sure to only take the columns which are also in xsample
    xpop <- dplyr::select_(xpop, .dots=names(xsample))
    xpop_d <- model.matrix(~., data = xpop)
    xpop_d <- apply(xpop_d,2,sum)
  }
  if (datatype=="totals"){
    #Make sure to only take the values which are also in xsample
    xpop_d <- unlist(c(N,xpop[names(xsample)]))
  }
  if (datatype=="means"){
    #Make sure to only take the values which are also in xsample
    xpop_d <- unlist(c(N,xpop[names(xsample)]*N))
  }
  
  #Total estimate
  t <- elasticNet_coef %*% (xpop_d) + t(y-y.hats.s)%*%pi^(-1)
  
  
  if ( var_est == TRUE ) {
    if ( var_method != "bootstrapSRS") {
      varEst <- varMase(y = (y-y.hats.s),pi = pi,pi2 = pi2,method = var_method, N = N)
      
    }
    
    if ( var_method == "bootstrapSRS"){
        #Find bootstrap variance
        
        #Sample data
        dat <- cbind(y,pi, xsample.d)
        #Bootstrap total estimates
        t_boot <- boot(data = dat, statistic = gregElasticNett, R = B, xpopd = xpop_d, alpha=alpha, lambda=lambda_opt, parallel = "multicore", ncpus = 2)
        
        #Adjust for bias and without replacement sampling
        varEst <- var(t_boot$t)*n/(n-1)*(N-n)/(N-1)
      }

  }
  
}   
  
  if(var_est==TRUE){

    return(list( pop_total = as.numeric(t), 
                 pop_mean = as.numeric(t)/N, 
                 pop_total_var=varEst, 
                 pop_mean_var=varEst/N^2,
                 coefficients = elasticNet_coef))
  }else{
    
    return(list( pop_total = as.numeric(t), 
                 pop_mean = as.numeric(t)/N, 
                 coefficients = elasticNet_coef))
    
  }
}
