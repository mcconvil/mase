#' Compute an elastic net regression estimator
#' 
#' Calculates a lasso, ridge or elastic net generalized regression estimator for a finite population mean/proportion or total based on sample data collected from a complex sampling design and auxiliary population data.  
#' 
#' 
#' @inheritParams horvitzThompson
#' @inheritParams greg
#' @param alpha A numeric value between 0 and 1 which signifies the mixing parameter for the lasso and ridge penalties in the elastic net.  When alpha = 1, only a lasso penalty is used.  When alpha = 0, only a ridge penalty is used. Default is alpha = 1. 
#' @param lambda A string specifying how to tune the lambda hyper-parameter.  Only used if modelselect = TRUE and defaults to "lambda.min". The possible values are "lambda.min", which is the lambda value associated with the minimum cross validation error or "lambda.1se", which is the lambda value associated with a cross validation error that is one standard error away from the minimum, resulting in a smaller model.
#' @param cvfolds The number of folds for the cross-validation to tune lambda.
#' @param weights_method A string specifying which method to use to calculate survey weights. Currently, "ridge" is the only option. The "ridge" method uses a ridge regression approximation to calculate weights (see McConville et al (2017), section 3.2 for details). Support for "calibration" to come soon, which employs the model calibration method of Wu and Sitter (2001).
#' @param eta A small positive number. Defaults to 0.0001. See McConville et al (2017), section 3.2 for details. 
#' 
#' @examples 
#' library(dplyr)
#' data(IdahoPop)
#' data(IdahoSamp)
#' xsample <- filter(IdahoSamp, COUNTYFIPS == "16055")
#' xpop <- filter(IdahoPop, COUNTYFIPS == "16055")
#' gregElasticNet(y = xsample$BA_TPA_ADJ,
#'                N = xpop$npixels,
#'                xsample = xsample[c("tcc", "elev", "ppt", "tmean")],
#'                xpop = xpop[c("tcc", "elev", "ppt", "tmean")],
#'                var_est = TRUE,
#'                var_method = "LinHB",
#'                datatype = "means",
#'                alpha = 0.5)

#' 
#' @references 
#' \insertRef{mcc17}{mase}

#'
#' @return 
#' A list of output containing:
#' 
#' * pop_total: Estimate of population total.
#' 
#' * coefficients: Elastic net coefficient estimates.
#' 
#' * pop_mean: Estimate of the population mean (or proportion).
#' 
#' * pop_total_var: Estimated variance of population total estimate.
#' 
#' * pop_mean_var:Estimated variance of population mean estimate.
#' 
#' @import boot
#' @import glmnet
#' @importFrom stats model.matrix predict quasibinomial var
#' @export gregElasticNet
#' @include varMase.R
#' @include bootHelpers.R

gregElasticNet  <- function(y,
                            xsample,
                            xpop,
                            pi = NULL,
                            alpha = 1,
                            model = "linear",
                            pi2 = NULL,
                            var_est = FALSE,
                            var_method="LinHB",
                            datatype = "raw",
                            N = NULL,
                            lambda = "lambda.min",
                            B = 1000,
                            cvfolds = 10,
                            weights_method = "ridge",
                            eta = 0.0001,
                            fpc = TRUE,
                            messages = TRUE){
  
  if (!(typeof(y) %in% c("numeric", "integer", "double"))) {
    stop("Must supply numeric y.  For binary variable, convert to 0/1's.")
  }
  
  if (!is.element(var_method, c("LinHB", "LinHH", "LinHTSRS", "LinHT", "bootstrapSRS"))) {
    stop("Variance method input incorrect. It has to be \"LinHB\", \"LinHH\", \"LinHT\", \"LinHTSRS\", or \"bootstrapSRS\".")
  }
  
  if (!is.element(model, c("linear","logistic"))) {
    stop("Method input incorrect, has to be either \"linear\" or \"logistic\"")
  }
  
  if (is.null(N)) {
    if (datatype == "raw") {
      N <- dim(as.matrix(xpop))[1]
    } else {
      N <- sum(pi^(-1))
      if (messages) {
        message("Assuming N can be approximated by the sum of the inverse inclusion probabilities.") 
      }
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
  if (is.null(pi)) {
    if (messages) {
      message("Assuming simple random sampling") 
    }
  }  
  
  if (is.null(pi)) {
    pi <- rep(length(y)/N, length(y))
  }
  
  weight <- as.vector(pi^(-1))
  
  #Cross-validation to find lambdas
  if (model == "linear") {
    fam <- "gaussian"
  } else{
    fam <- "binomial"
  } 
  
  cv <- cv.glmnet(x = as.matrix(xsample),
                  y = y,
                  alpha = alpha,
                  weights = weight,
                  nfolds = cvfolds,
                  family = fam,
                  standardize = FALSE)
  
  if (lambda == "lambda.min") {
    lambda_opt <- cv$lambda.min
  }
  if (lambda == "lambda.1se") {
    lambda_opt <- cv$lambda.1se
  }
  
  pred.mod <- glmnet(x = as.matrix(xsample),
                     y = y,
                     alpha = alpha,
                     family = fam,
                     standardize = FALSE,
                     weights = weight)
  
  elasticNet_coef <- predict(pred.mod, type = "coefficients", s = lambda_opt)[1:dim(xsample.d)[2],]
  
  #Estimated y values in sample
  y.hats.s <- predict(cv, newx = as.matrix(xsample), s = lambda_opt, type = "response")

  if (model == "logistic") {
    if (datatype != "raw") {
      stop("For the Logistic Elastic Net Estimator, user must supply all x values for population.  Populations totals or means for x are not enough.")
    }
    
    xpop <- data.frame(model.matrix(~., data = xpop))[,-1]
    xpop <- dplyr::select(xpop, names(xsample))
    xpop_d <- model.matrix(~., data = xpop)
    
    #Total estimate
    y.hats.U <- predict(cv,newx = xpop_d[,-1], s = lambda_opt, type = "response")
    t <- sum(y.hats.U) + t(y-y.hats.s)%*%pi^(-1)
    
    if ( var_est == TRUE){
      if (var_method != "bootstrapSRS") {
        varEst <- varMase(y = (y - y.hats.s),pi = pi,pi2 = pi2,method = var_method, N = N, fpc = fpc)
      }
      
      if(var_method == "bootstrapSRS"){
        #FILL IN: logistic, raw data!
        
        #Sample data
        dat <- cbind(y,pi, xsample.d)
        #Bootstrap total estimates
        t_boot <- boot(data = dat, statistic = logisticGregElasticNett, R = B, xpopd = xpop_d, alpha=alpha, lambda=lambda_opt, parallel = "multicore", ncpus = 2)
        
        #Adjust for bias and without replacement sampling
        if (fpc == T) {
          varEst <- var(t_boot$t)*n/(n-1)*(N-n)/(N-1) 
        }
        if (fpc == F) {
          varEst <- var(t_boot$t)*n/(n-1)
        }
      } 
      
    }
  
  }

  if (model == "linear") {
  
    if (datatype == "raw") {
      xpop <- data.frame(model.matrix(~., data = xpop))[, -1]
      xpop <- dplyr::select(xpop, names(xsample))
      xpop_d <- model.matrix(~., data = xpop)
      xpop_d <- apply(xpop_d,2,sum)
    }
    if (datatype == "totals") {
      xpop_d <- unlist(c(N, xpop[names(xsample)]))
    }
    if (datatype=="means") {
      xpop_d <- unlist(c(N, xpop[names(xsample)]*N))
    }
    
    if (weights_method == "ridge") {
      abs_beta_hat <- abs(c(0, elasticNet_coef[2:length(elasticNet_coef)]))
      Q_inverse <- diag(abs_beta_hat + eta)
      w <- as.matrix(
        1 + t(as.matrix(xpop_d) - xsample.dt %*% weight) %*% 
          solve(xsample.dt %*% diag(weight) %*% xsample.d + Q_inverse) %*% 
          (xsample.dt)
        ) %*% 
        diag(weight)
    }
  
    #Total estimate
    t <- elasticNet_coef %*% (xpop_d) + t(y-y.hats.s) %*% pi^(-1)
    
    if ( var_est == TRUE) {
      
      if (var_method != "bootstrapSRS") {
        varEst <- varMase(y = (y - y.hats.s),
                          pi = pi,
                          pi2 = pi2,
                          method = var_method,
                          N = N,
                          fpc = fpc)
      } else if (var_method == "bootstrapSRS") {

        dat <- cbind(y,pi, xsample.d)
        #Bootstrap total estimates
        t_boot <- boot(data = dat,
                       statistic = gregElasticNett,
                       R = B,
                       xpopd = xpop_d,
                       alpha = alpha,
                       lambda = lambda_opt,
                       parallel = "multicore",
                       ncpus = 2)
        
        #Adjust for bias and without replacement sampling
        if (fpc == T) {
          varEst <- var(t_boot$t)*n/(n-1)*(N-n)/(N-1)
        }
        if (fpc == F) {
          vvarEst <- var(t_boot$t)*n/(n-1)
        }
      }
    }
  }   
  
  if (var_est == TRUE) {

    return(list(pop_total = as.numeric(t), 
                pop_mean = as.numeric(t)/N, 
                pop_total_var=varEst, 
                pop_mean_var=varEst/N^2,
                weights = as.vector(w),
                coefficients = elasticNet_coef))
  } else {
    
    return(list(pop_total = as.numeric(t), 
                pop_mean = as.numeric(t)/N, 
                weights = as.vector(w),
                coefficients = elasticNet_coef))
    
  }
}
