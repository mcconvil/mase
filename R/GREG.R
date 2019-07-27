#' Compute a generalized regression estimator
#' 
#' Calculates a generalized regression estimator for a finite population mean/proportion or total based on sample data collected from a complex sampling design and auxiliary population data.  
#' 
#' @inheritParams horvitzThompson
#' @param x_sample A data frame of the auxiliary data in the sample.
#' @param x_pop A data frame of population level auxiliary information.  It must contain the same names as x_sample.  If data_type = "raw", must contain unit level data.  If data_type = "totals" or "means", then contains one row of aggregated, population totals or means for the auxiliary data. Default is "raw".
#' @param data_type A string that specifies the form of population auxiliary data. The possible values are "raw", "totals" or "means" for whether the user is providing population data at the unit level, aggregated to totals, or aggregated to means.  Default is "raw".
#' @param model A string that specifies the regression model to utilize. Options are "linear" or "logistic".
#' @param model_select A logical for whether or not to run lasso regression first and then fit the model using only the predictors with non-zero lasso coefficients. Default is FALSE.  
#' @param lambda A string specifying how to tune the lasso hyper-parameter.  Only used if model_select = TRUE and defaults to "lambda.min". The possible values are "lambda.min", which is the lambda value associated with the minimum cross validation error or "lambda.1se", which is the lambda value associated with a cross validation error that is one standard error away from the minimum, resulting in a smaller model.
#' @param hetero Set `TRUE` if model residuals show linear heteroskedasticity (non constant variance, "fanning", etc.). Alternatively, if an x variable is known to be heteroskedastic with y, provide a character vector of variable names. When used, the model weights will be divided by model estimated residual variance at each datum in x_sample. Not applicable for logistic regression.
#' @param ncpus Integer denoting the number of cpu cores to use in parallel for bootstrap variance.
#'
#' @examples 
#' library(survey)
#' data(api)
#' greg(y = apisrs$api00, x_sample = apisrs[c("col.grad", "awards")], 
#' x_pop = apipop[c("col.grad", "awards")], pi = apisrs$pw^(-1), 
#' var_est = TRUE)
#' 
#' #To estimate a proportion
#' y <- 0 + (apisrs$both == "Yes")
#' greg(y = y, x_sample = apisrs[c("col.grad")], 
#' x_pop = apipop[c("col.grad")], pi = apisrs$pw^(-1), 
#' var_est = TRUE, model = "logistic")
#' 
#'@references 
#'\insertRef{cas76}{mase}
#'
#'\insertRef{sar92}{mase}
#' 
#' @return A greg object containing:
#' \itemize{
#' \item{pop_total: Estimate of population total}
#' \item{pop_mean: Estimate of the population mean}
#' \item{pop_total_var: Estimated variance of population total estimate}
#' \item{pop_mean_var: Estimated variance of population mean estimate}
#' \item{weights: Survey weights produced by greg (linear model only)}
#' \item{coefficients: Survey-weighted model coefficients}
#' \item{logistic_model: glmnet logistic regression model object (logistic regression model only)}
#' \item{y_hat_sample: Response estimates for sample}
#' \item{y_hat_pop: Response estimates for population}
#' }
#' 
#' @export greg
#' @import survey
#' @import glmnet
#' @import boot
#' @importFrom stats model.matrix predict quasibinomial var
#' @include varMase.R
#' @include gregt.R
#' @include gregObject.R
#' 
#' @seealso \code{\link{gregElasticNet}} for a penalized regression model.

greg  <- function(y, x_sample, x_pop, pi = NULL, model = "linear",  pi2 = NULL,
                  var_est = FALSE, var_method = "lin_HB", data_type = "raw",
                  N = NULL, model_select = FALSE, lambda = "lambda.min", B = 1000,
                  strata = NULL, hetero = FALSE, ncpus = 1){

  
  
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
  
  if(model == "logistic" & data_type != "raw"){
    message("Must supply the raw population data to fit the logistic regression estimator.")
    return(NULL)
    
  }
  
  if(!is.element(data_type, c("raw","totals", "means"))){
    message("data_type input incorrect, has to be either \"raw\", \"totals\" or \"means\"")
    return(NULL)
  }
  
  if(class(x_sample) != "data.frame"){
    message("x_sample must be a data.frame.")
    return(NULL)
  }
  
  if(model_select == TRUE & dim(x_sample)[2] < 2){
    message("Only conduct model selection if you have at least two predictors. Set model_select = FALSE.")
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
  
  #Convert y to a vector
  y <- as.vector(y)
  
  #sample size
  n <- length(y)
  
  #create design matrix, x matrix and transpose design matrix
  x_sample_d <- model.matrix(~., data = data.frame(x_sample))
  x_sample <- data.frame(x_sample_d[,-1, drop = FALSE])
  x_sample_dt <- t(x_sample_d) 
 
  
  #Check on inclusion probabilities and create weight=inverse inclusion probabilities
  if(is.null(pi)){
    message("Assuming simple random sampling")
  }  
  
  
  # create pi if not given
  if (is.null(pi)) {
    pi <- rep(length(y)/N, length(y))
  }
  
  #weight: inverse first order inclusion probabilities
  weight <- as.vector(pi^(-1))
  
  if(model_select == TRUE){
    
    #Cross-validation to find lambda
    if(model == "linear"){
      fam <- "gaussian"
    } else{
      fam <- "binomial"
    } 
    
    #Run cv to find optimal lambda
    cv <- cv.glmnet(x = as.matrix(x_sample), y = y, alpha = 1, weights = weight,
                    nfolds = 10, family = fam, standardize = FALSE)
    
    #Pick lambda
    if(lambda == "lambda.min"){
      lambda_opt <- cv$lambda.min
    }
    if(lambda == "lambda.1se"){
      lambda_opt <- cv$lambda.1se
    }
    
    ## MODEL SELECTION COEFFICIENTS ##
    pred_mod <- glmnet(x = as.matrix(x_sample), y = y, alpha = 1, family = fam, standardize = FALSE, weights=weight)
    lasso_coef <- predict(pred_mod,type = "coefficients", s = lambda_opt)[1:dim(x_sample_d)[2],]
    
    #Collect the names of the variables with non-zero coefficients
    coef_select <- names(lasso_coef[lasso_coef != 0])[-1]
    
    #If select zero predictors, then fit a HT
    if(length(coef_select) == 0){
      message("No variables selected in the model selection stage.  Fitting a HT estimator.")
    
    if(var_est == TRUE){
      
      HT <- horvitzThompson(y = y, pi = pi, N = N, pi2 = pi2, var_est = TRUE, var_method = var_method)
      
      return(list( pop_total = HT$pop_total, 
                   pop_mean = HT$pop_total/N,
                   pop_total_var = HT$pop_total_var, 
                   pop_mean_var = HT$pop_total_var/N^2, 
                   weights = as.vector(pi^{-1})) %>%
               gregify())
      }
      else {
      
      HT <- horvitzThompson(y = y, pi = pi, N = N, pi2 = pi2, var_est = FALSE)
      
      return(list( pop_total = HT$pop_total, 
                   pop_mean = HT$pop_total/N,
                   weights = as.vector(pi^{-1})) %>%
               gregify())
      }
    }  
     
    #Create a new x_sample with only the columns in coef_select
    x_sample <- x_sample[,coef_select, drop = FALSE]#dplyr::select_(x_sample, .dots=coef_select)
    x_sample_d <- model.matrix(~., data = x_sample)
    x_sample_dt <- t(x_sample_d) 
    x_sample <- data.frame(x_sample_d[,-1, drop=FALSE])
    
  }
  
### LOGISTIC REGRESSION ###
  if (model == "logistic"){
    
    #Error if y has more than two categories 
    if(length(levels(as.factor(y)))!=2){
      message("Function can only handle categorical response with two categories.")
      return(NULL)
    }
    
    if (data_type=="raw"){
      x_pop <- data.frame(model.matrix(~., data = x_pop))[,-1, drop = FALSE]
      #Make sure to only take the columns which are also in x_sample
      x_pop <- dplyr::select(x_pop, dplyr::one_of(colnames(x_sample)))
      x_pop_d <- model.matrix(~., data = x_pop)
    }
    
    #Fit model using survey's svyglm since glm doesn't handle weights as we want
    dat <- data.frame(y, weight, x_sample)
    colnames(dat) <- c("y", "weight", colnames(x_sample))
    f <- paste(names(dat)[1], "~", paste(names(dat)[-c(1,2)], collapse=" + "))
    s_design <-survey::svydesign(ids=~1,weights=~weight,data=dat)
    mod <- survey::svyglm(f, design=s_design,family=quasibinomial())

    y_hats_U <- as.matrix(predict(mod, newdata = data.frame(x_pop_d[,-1, drop = FALSE]), type = "response", family = quasibinomial()))
    y_hats_s <- as.matrix(predict(mod, type = "response", family = quasibinomial()))
    #Estimator of total
    t <- t(y-y_hats_s)%*%pi^(-1) + sum(y_hats_U)
    
    #Coefficients
    coefs <- mod$coefficients
    
    #No weights
    w <- NULL
    
    if(var_est==TRUE){
      if(var_method!="bootstrap_SRS"){
        e <- y-y_hats_s
        varEst <- varMase(y = e,pi = pi,pi2 = pi2,method = var_method, N = N, strata = strata)
        }
      else if(var_method=="bootstrap_SRS"){
  
        #Find bootstrap variance
 
        #Bootstrap total estimates
        dat <- data.frame(y, weight, x_sample)
        colnames(dat) <- c("y", "weight", colnames(x_sample))
        
        t_boot <-  boot(data = dat, statistic = logisticGregt, R = B, x_pop_d = x_pop_d, parallel = "multicore", ncpus = ncpus)
    
       
        #Adjust for bias and without replacement sampling
        varEst <- var(t_boot$t)*n/(n-1)*(N-n)/(N-1)
      }
      return(list( pop_total = as.numeric(t),
                   pop_mean = as.numeric(t)/N,
                   pop_total_var=varEst, 
                   pop_mean_var=varEst/N^2,
                   coefficients = coefs,
                   logistic_model = mod,
                   y_hat_sample = y_hats_s,
                   y_hat_pop = y_hats_U) %>%
               gregify())
      }else{
        return(list( pop_total = as.numeric(t), 
                   pop_mean = as.numeric(t)/N,
                   coefficients =  coefs,
                   logistic_model = mod,
                   y_hat_sample = y_hats_s,
                   y_hat_pop = y_hats_U) %>%
               gregify())
      }
    }
### LINEAR REGRESSION ###
  else if (model == "linear"){
    
    #population design matrix, check whether its population totals, means or raw data
    
    if (data_type=="raw"){
      x_pop <- data.frame(model.matrix(~.-1, data = data.frame(x_pop)))
      #Make sure to only take the columns which are also in x_sample
      x_pop <- dplyr::select(x_pop, dplyr::one_of(colnames(x_sample)))
      x_pop_d <- model.matrix(~., data = x_pop)
      x_pop_d <- apply(x_pop_d,2,sum)
      }
    if (data_type=="totals"){
      #Make sure to only take the values which are also in x_sample
      x_pop_d <- unlist(c(N,x_pop[names(x_sample)]))
      }
    if (data_type=="means"){
      #Make sure to only take the values which are also in x_sample
      x_pop_d <- unlist(c(N,x_pop[names(x_sample)]*N))
      }
    
    w <- as.matrix(1 + t(as.matrix(x_pop_d) - x_sample_dt %*% weight) %*% solve(x_sample_dt %*% diag(weight) %*% x_sample_d) %*% (x_sample_dt)) %*% diag(weight)
  
  #calculating the total estimate for y
  t <- w %*% y
  

  #NOTE: check that weights times x's should equal total of x's to check for correct weight values
  
  #Coefficients
  coefs <- solve(x_sample_dt %*% diag(weight) %*% x_sample_d) %*% (x_sample_dt) %*% diag(weight) %*% y
  
  #model estimates
  y_hat <- (x_sample_d%*%solve(x_sample_dt %*% diag(weight) %*% x_sample_d) %*%
    (x_sample_dt) %*% diag(weight)%*%y) %>% as.vector()
  #heteroskedasticity
  resid_mod <- NULL
  if(hetero != FALSE){
    selected <- coef_select
    message("Adjusting sample weights for heteroskedasticity")
    vf <- as.formula(paste("e2 ~ ", paste(selected, collapse= "+")))
    if(typeof(hetero) == "character"){
      if(sum(hetero %in% selected) == length(hetero)){
        vf <- as.formula(paste("e2 ~ ", paste(hetero, collapse= "+")))
      }
      else{
        message(paste0("Heteroskedastic variable(s) supplied was not selected by LASSO: ",
                       paste(selected), ". Using selected variables instead."))
      }
    }
    e2 <- ((y-y_hat)^2) #fit model on squared residuals
    resid_mod <- lm(vf, data = cbind(e2, x_sample))
    #adjust weights, truncate if not positive
    weight <- weight/sqrt(ifelse(resid_mod$fitted.values <= 0,
                                 0, #truncate to zero and then add small constant
                                 resid_mod$fitted.values) + median(e2)/n)
    #refit original model
    cv <- cv.glmnet(x = as.matrix(cbind(rep(1, nrow(x_sample)), x_sample)),
                    y = y, alpha = 1, weights = weight, nfolds = 10, family = fam,
                    standardize = FALSE)
    if(lambda == "lambda.min"){
      lambda_opt <- cv$lambda.min
    }
    if(lambda == "lambda.1se"){
      lambda_opt <- cv$lambda.1se
    }
    pred_mod <- glmnet(x = as.matrix(cbind(rep(1, nrow(x_sample)), x_sample)),
                       y = y, alpha = 1, family = fam, standardize = FALSE, weights=weight)
    lasso_coef <- predict(pred_mod,type = "coefficients", s = lambda_opt)[1:dim(x_sample_d)[2],]
    coef_select <- names(lasso_coef[lasso_coef != 0])[-1]
    w <- as.matrix(1 + t(as.matrix(x_pop_d) - x_sample_dt %*% weight) %*% solve(x_sample_dt %*% diag(weight) %*% x_sample_d) %*% (x_sample_dt)) %*% diag(weight)
    t <- w %*% y
    coefs <- solve(x_sample_dt %*% diag(weight) %*% x_sample_d) %*% (x_sample_dt) %*% diag(weight) %*% y
    y_hat <- (x_sample_d%*%solve(x_sample_dt %*% diag(weight) %*% x_sample_d) %*%
                (x_sample_dt) %*% diag(weight)%*%y) %>% as.vector()
  }
  y_hat_pop <- NULL
  if(data_type == "raw") {
    y_hat_pop <- (model.matrix(~., data = x_pop) %*% coefs) %>% as.vector()
  }
  
  if(var_est==TRUE){
    if(var_method!="bootstrap_SRS"){
      e <- y-y_hat
      varEst <- varMase(y = e,pi = pi, pi2 = pi2, method = var_method, N = N, strata = strata)
      }
    else if(var_method=="bootstrap_SRS"){
      #Find bootstrap variance via residual bootstrap
      dat <- cbind(y-y_hat, pi)
      #Bootstrap total estimates
      t_boot <- boot(data = dat, statistic = gregt, R = B, y_hat = y_hat, x_pop_d = x_pop_d,
                     x_sample_d = x_sample_d, parallel = "multicore", ncpus = 2)
      #Adjust for bias and without replacement sampling
      varEst <- var(t_boot$t)*n/(n-1)*(N-n)/(N-1)
    }
    return(list( pop_total = as.numeric(t), 
                 pop_mean = as.numeric(t)/N,
                 pop_total_var=varEst, 
                 pop_mean_var=varEst/N^2, 
                 weights = as.vector(w),
                 coefficients =  coefs,
                 y_hat_sample = y_hat,
                 y_hat_pop = y_hat_pop,
                 skedastic_mod = resid_mod) %>%
             gregify())
    }
  else{
    return(list( pop_total = as.numeric(t), 
                 pop_mean = as.numeric(t)/N,           
                 weights = as.vector(w),
                 coefficients =  coefs,
                 y_hat_sample = y_hat,
                 y_hat_pop = y_hat_pop,
                 skedastic_mod = resid_mod) %>%
             gregify())      
    }
  }
}

