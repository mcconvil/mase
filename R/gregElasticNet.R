#' Compute an elastic net regression estimator
#' 
#' Calculates a lasso, ridge or elastic net generalized regression estimator for a finite population mean/proportion or total based on sample data collected from a complex sampling design and auxiliary population data.  
#' 
#' 
#' @inheritParams horvitzThompson
#' @inheritParams greg
#' @param alpha A numeric value between 0 and 1 which signifies the mixing parameter for the lasso and ridge penalties in the elastic net.  When alpha = 1, only a lasso penalty is used.  When alpha = 0, only a ridge penalty is used. Default is alpha = 1. 
#' @param lambda A string specifying how to tune the lambda hyper-parameter.  Only used if model_select = TRUE and defaults to "lambda.min". The possible values are "lambda.min", which is the lambda_value associated with the minimum cross validation error or "lambda.1se", which is the lambda value associated with a cross validation error that is one standard error away from the minimum, resulting in a smaller model.
#' @param cvfolds The number of folds for the cross-validation to tune lambda.
#' @param hetero Set `TRUE` if model residuals show linear heteroskedasticity (non constant variance, "fanning", etc.). Alternatively, if an x variable is known to be heteroskedastic with y, provide a character vector of variable names. When used, the model weights will be divided by model estimated residual variance at each datum in x_sample.
#' 
#' @examples 
#' library(survey)
#' data(api)
#' gregElasticNet(y = apisrs$api00, 
#' x_sample = apisrs[c("col.grad", "awards", "snum", "dnum", "cnum", "pcttest", "meals", "sch.wide")], 
#' x_pop = apipop[c("col.grad", "awards", "snum", "dnum", "cnum", "pcttest", "meals", "sch.wide")], 
#' pi = apisrs$pw^(-1), var_est = TRUE, alpha = .5)
#' 
#' @references 
#'\insertRef{mcc17}{mase}

#'
#' @return A greg object containing:
#' \itemize{
#' \item{pop_total: Estimate of population total}
#' \item{pop_mean: Estimate of the population mean}
#' \item{pop_total_var: Estimated variance of population total estimate}
#' \item{pop_mean_var: Estimated variance of population mean estimate}
#' \item{coefficients: Survey-weighted model coefficients}
#' \item{formula: Model formula}
#' \item{model: glmnet model object}
#' \item{y_hat_sample: Response estimates for sample}
#' \item{y_hat_pop: Response estimates for population}
#' }
#' @import boot
#' @import glmnet
#' @import foreach
#' @import Matrix
#' @importFrom stats model.matrix predict quasibinomial var
#' @export gregElasticNet
#' @include varMase.R
#' @include gregElasticNett.R
#' @include gregObject.R
#' 
#' @seealso \code{\link{greg}} for a linear or logistic regression model.


gregElasticNet  <- function(
  y, x_sample, x_pop, pi = NULL, alpha = 1, model = "linear", pi2 = NULL, var_est = FALSE, var_method = "lin_HB", 
  data_type = "raw", N = NULL, lambda = "lambda.min", B = 1000, cvfolds = 10, strata = NULL, standardize = FALSE,
  hetero = FALSE){
  
  
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
  
  if(class(x_sample) != "data.frame"){
    message("x_sample must be a data.frame.")
    return(NULL)
  }
  
  #Need to provide either data_type="raw", N, or pi.  Give warning if not
  if(data_type %in% c("means", "totals") & is.null(N) & is.null(pi)){
    message("Must supply N, pi, or raw population data so that we can estimate N.")
    return(NULL)
  }

  #Need to get N if not provided
  if(is.null(N)){
    if(data_type=="raw"){
      N <- dim(as.matrix(x_pop))[1]
    }else{
      N <- sum(pi^(-1))
      message("Assuming N can be approximated by the sum of the inverse inclusion probabilities.")
    }
  }
  
  #Check standardization
  if(standardize == TRUE){
    x_sample <- base::scale(x_sample) %>% as.data.frame()
    x_pop <- base::scale(x_pop, center = colMeans(x_sample),
                         scale = apply(x_sample, 2, sd)) %>%
      as.data.frame()
  }
  
  #create design matrix, x matrix and transpose design matrix
  x_sample_d <- model.matrix(~., data = data.frame(x_sample))
  x_sample <- data.frame(x_sample_d[,-1])
  x_sample_dt <- t(x_sample_d) 
  
  #Format y
  y <- as.vector(y)
  n <- length(y)
  
  #Make sure y is complete
  if(NA %in% y){
    message("Must supply complete cases for y.")
    return(NULL)
  }
  #Check that the number of observations are consistent
  if(nrow(x_sample) != length(y)){
    message("y and x_sample must have same number of observations.")
    return(NULL)
  }
  #Check on inclusion probabilities and create weight=inverse inclusion probabilities
  if(is.null(pi)){
    message("Assuming simple random sampling")
    pi <- rep(length(y)/N, length(y))
  }  
  #Check for missing data:
  if(FALSE %in% complete.cases(x_sample) || FALSE %in% complete.cases(x_pop)){
    if(FALSE %in% complete.cases(x_sample)){
      message("Must supply complete cases for x_sample")
    }
    if(FALSE %in% complete.cases(x_pop)){
      message("Must supply complete cases for x_pop")
    }
    return(NULL)
  }
  
  #weight: inverse first order inclusion probabilities
  weight <- as.vector(pi^(-1))
  
  #Cross-validation to find lambda
  if(model=="linear"){
    fam <- "gaussian"
  } else{
    fam <- "binomial"
  } 
  
  cv <- cv.glmnet(x = as.matrix(x_sample), y = y, alpha = alpha, weights = weight,
                  nfolds = cvfolds,family=fam, standardize = FALSE)
  
  if(lambda == "lambda.min"){
    lambda_opt <- cv$lambda.min
  }
  if(lambda == "lambda.1se"){
    lambda_opt <- cv$lambda.1se
  }
  
  #formula
  f <- as.formula(paste("y ~ ", paste(names(x_sample), collapse= "+")))
  
  ## MODEL SELECTION COEFFICIENTS ##
  pred_mod <- glmnet(x = as.matrix(x_sample), y = y, alpha = alpha, family=fam,
                     standardize = FALSE, weights=weight)
  elasticNet_coef <- predict(pred_mod,type = "coefficients",
                             s = lambda_opt)[1:dim(x_sample_d)[2],]
  #Estimated y values in sample
  y_hats_s <- as.vector(predict(cv,newx = as.matrix(x_sample), s = lambda_opt, type="response"))
  #Heteroskedasticity weighting
  resid_mod <- NULL
  if(hetero != FALSE){
    selected <- names(x_sample)[which(elasticNet_coef != 0)-1]
    message("Adjusting sample weights for heteroskedasticity")
    vf <- as.formula(paste("e2 ~ ", paste(selected, collapse= "+")))
    if(typeof(hetero) == "character"){
      if(sum(hetero %in% selected) == length(hetero)){
        vf <- as.formula(paste("e2 ~ ", paste(hetero, collapse= "+")))
      }
      else{
        message(paste0("Heteroskedastic variable(s) supplied was not selected by ElasticNet: ",
                   paste(selected), ". Using selected variables instead."))
      }
    }
    e2 <- ((y-y_hats_s)^2) #fit model on squared residuals
    resid_mod <- lm(vf, data = cbind(e2, x_sample))
    #adjust weights, truncate if not positive
    weight <- weight/ifelse(resid_mod$fitted.values <= 0,
                            sd(resid_mod$fitted.values)/n,
                            resid_mod$fitted.values)
    #refit model
    pred_mod <- glmnet(x = as.matrix(x_sample), y = y, alpha = alpha, family=fam,
                       standardize = FALSE, weights=weight)
    elasticNet_coef <- predict(pred_mod,type = "coefficients",
                               s = lambda_opt)[1:dim(x_sample_d)[2],]
    y_hats_s <- as.vector(predict(cv,newx = as.matrix(x_sample), s = lambda_opt, type="response"))
  }
if (model == "logistic") {
  if (data_type != "raw"){
    message("For the Logistic Elastic Net Estimator, user must supply all x values for population.  Populations totals or means for x are not enough.")
    return(NULL)
  }
  #Population matrix
  x_pop <- data.frame(model.matrix(~., data = x_pop))[,-1]
  #Make sure to only take the columns which are also in x_sample
  x_pop <- dplyr::select_(x_pop, .dots=names(x_sample))
  x_pop_d <- model.matrix(~., data = x_pop)
  #Estimated y values in population
  y_hats_U <- predict(cv,newx = x_pop_d[,-1], s = lambda_opt, type = "response")
  #Total estimate
  t <- sum(y_hats_U) + t(y-y_hats_s)%*%pi^(-1)
  
  if ( var_est == TRUE){
    if (var_method != "bootstrap_SRS") {
      varEst <- varMase(y = (y - y_hats_s), pi = pi, pi2 = pi2, method = var_method, N = N, strata = strata)
      
    }
    
    if(var_method == "bootstrap_SRS"){
      #FILL IN: logistic, raw data!
      
      #Sample data
      dat <- cbind(y,pi, x_sample_d)
      #Bootstrap total estimates
      t_boot <- boot(data = dat, statistic = logisticGregElasticNett, R = B, x_pop_d = x_pop_d, alpha=alpha, lambda = lambda_opt, parallel = "multicore", ncpus = 2)
      
      #Adjust for bias and without replacement sampling
      varEst <- var(t_boot$t)*n/(n-1)*(N-n)/(N-1)
    } 
    
  }

  
}

if (model == "linear") {
  
  #Format x_pop to be a vector of pop totals
  if (data_type=="raw"){
    x_pop <- data.frame(model.matrix(~., data = x_pop))[,-1]
    #Make sure to only take the columns which are also in x_sample
    x_pop <- dplyr::select_(x_pop, .dots=names(x_sample))
    x_pop_d <- model.matrix(~., data = x_pop)
    #Estimated y values in population
    y_hats_U <- predict(cv,newx = x_pop_d[,-1], s = lambda_opt, type = "response")
    #Column sums
    x_pop_d <- apply(x_pop_d,2,sum)
  }
  if (data_type=="totals"){
    #Make sure to only take the values which are also in x_sample
    x_pop_d <- unlist(c(N,x_pop[names(x_sample)]))
    y_hats_U <- NULL
  }
  if (data_type=="means"){
    #Make sure to only take the values which are also in x_sample
    x_pop_d <- unlist(c(N,x_pop[names(x_sample)]*N))
    y_hats_U <- NULL
  }
  
  #Total estimate
  t <- elasticNet_coef %*% (x_pop_d) + t(y-y_hats_s)%*%pi^(-1)
  
  
  if ( var_est == TRUE ) {
    if ( var_method != "bootstrap_SRS") {
      varEst <- varMase(y = (y-y_hats_s), pi = pi, pi2 = pi2, method = var_method, N = N, strata = strata)
      
    }
    
    if ( var_method == "bootstrap_SRS"){
      #Find bootstrap variance via residual bootstrap
      dat <- cbind(y-y_hats_s, pi)
      #Bootstrap total estimates
      t_boot <- boot(data = dat, statistic = gregElasticNett, R = B, x_pop_d = x_pop_d,
                     y = y, x_sample_d = x_sample_d, alpha=alpha, lambda = lambda_opt,
                     parallel = "multicore", ncpus = 2)
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
                 formula = f,
                 coefficients = elasticNet_coef,
                 model = pred_mod,
                 y_hat_sample = y_hats_s,
                 y_hat_pop = y_hats_U %>% as.vector(),
                 resid_mod = resid_mod) %>%
             gregify())
  }else{
    
    return(list( pop_total = as.numeric(t), 
                 pop_mean = as.numeric(t)/N, 
                 formula = f,
                 coefficients = elasticNet_coef,
                 model = pred_mod,
                 y_hat_sample = y_hats_s,
                 y_hat_pop = y_hats_U %>% as.vector(),
                 resid_mod = resid_mod) %>%
             gregify())
    
  }
}
