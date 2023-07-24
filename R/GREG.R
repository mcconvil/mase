#' Compute a generalized regression estimator
#' 
#' Calculates a generalized regression estimator for a finite population mean/proportion or total based on sample data collected from a complex sampling design and auxiliary population data.  
#' 
#' @inheritParams horvitzThompson
#' @param xsample A data frame of the auxiliary data in the sample.
#' @param xpop A data frame of population level auxiliary information.  It must contain the same names as xsample.  If datatype = "raw", must contain unit level data.  If datatype = "totals" or "means", then contains one row of aggregated, population totals or means for the auxiliary data. Default is "raw".
#' @param datatype A string that specifies the form of population auxiliary data. The possible values are "raw", "totals" or "means" for whether the user is providing population data at the unit level, aggregated to totals, or aggregated to means.  Default is "raw".
#' @param model A string that specifies the regression model to utilize. Options are "linear" or "logistic".
#' @param modelselect A logical for whether or not to run lasso regression first and then fit the model using only the predictors with non-zero lasso coefficients. Default is FALSE.  
#' @param lambda A string specifying how to tune the lasso hyper-parameter.  Only used if modelselect = TRUE and defaults to "lambda.min". The possible values are "lambda.min", which is the lambda value associated with the minimum cross validation error or "lambda.1se", which is the lambda value associated with a cross validation error that is one standard error away from the minimum, resulting in a smaller model.
#' 
#' @examples 
#' library(survey)
#' data(api)
#' greg(y = apisrs$api00, xsample = apisrs[c("col.grad", "awards")], 
#' xpop = apipop[c("col.grad", "awards")], pi = apisrs$pw^(-1), 
#' var_est = TRUE)
#' 
#'@references 
#'\insertRef{cas76}{mase}
#'
#'\insertRef{sar92}{mase}
#' 
#' @return A list of output containing:
#' \itemize{
#' \item{pop_total:}{Estimate of population total}
#' \item{pop_mean:}{Estimate of the population mean (or proportion)}
#' \item{weights:}{Survey weights produced by GREG (linear model only)}
#' \item{pop_total_var:}{Estimated variance of population total estimate}
#' \item{pop_mean_var:}{Estimated variance of population mean estimate}
#' }
#' 
#' @export greg
#' @import survey
#' @import glmnet
#' @import boot
#' @importFrom stats model.matrix predict quasibinomial var
#' @include varMase.R
#' @include gregt.R

greg  <- function(y,
                  xsample,
                  xpop,
                  pi = NULL,
                  model = "linear",
                  pi2 = NULL, 
                  var_est = FALSE,
                  var_method = "LinHB",
                  datatype = "raw",
                  N = NULL,
                  modelselect = FALSE,
                  lambda = "lambda.min",
                  B = 1000,
                  messages = T){

  
  
### INPUT VALIDATION ###
  
  #Check that y is numeric
  if(!(typeof(y) %in% c("numeric", "integer", "double"))){
    stop("Must supply numeric y.  For binary variable, convert to 0/1's.")
  }
  
  
  #Make sure the var_method is valid
  if(!is.element(var_method, c("LinHB", "LinHH", "LinHTSRS", "LinHT", "bootstrapSRS"))){
    stop("Variance method input incorrect. It has to be \"LinHB\", \"LinHH\", \"LinHT\", \"LinHTSRS\", or \"bootstrapSRS\".")
  }
  
  
  if(!is.element(model, c("linear","logistic"))){
    stop("Method input incorrect, has to be either \"linear\" or \"logistic\"")
  }
  
  if(model == "logistic" & datatype != "raw"){
    stop("Must supply the raw population data to fit the logistic regression estimator.")
  }
  
  if(!is.element(datatype, c("raw","totals", "means"))){
    stop("datatype input incorrect, has to be either \"raw\", \"totals\" or \"means\"")
  }
  
  
  #Need to get N if not provided
  if(is.null(N)){
    if(datatype=="raw"){
      N <- dim(as.matrix(xpop))[1]
    }else{
      N <- sum(pi^(-1))
      if (messages) {
        message("Assuming N can be approximated by the sum of the inverse inclusion probabilities.") 
      }
    }
  }
  
  #Convert y to a vector
  y <- as.vector(y)
  
  #sample size
  n <- length(y)
  
  #create design matrix, x matrix and transpose design matrix
  xsample.d <- model.matrix(~., data = data.frame(xsample))
  xsample <- data.frame(xsample.d[,-1, drop = FALSE])
  xsample.dt <- t(xsample.d) 
 
  
  #Check on inclusion probabilities and create weight=inverse inclusion probabilities
  if(is.null(pi)){
    if (messages) {
      message("Assuming simple random sampling") 
    }
  }  
  
  
  # convert pi into diagonal matrix format
  if (is.null(pi)) {
    pi <- rep(length(y)/N, length(y))
  }
  
  #weight: inverse first order inclusion probabilities
  weight <- as.vector(pi^(-1))
  
  if(modelselect == TRUE){
    
    #Cross-validation to find lambdas
    if(model == "linear"){
      fam <- "gaussian"
    } else{
      fam <- "binomial"
    } 
    
    #Run cv to find optimal lambda
    cv <- cv.glmnet(x = as.matrix(xsample), y = y, alpha = 1, weights = weight, nfolds = 10, family = fam, standardize = FALSE)
    
    #Pick lambda
    if(lambda=="lambda.min"){
      lambda.opt <- cv$lambda.min
    }
    if(lambda=="lambda.1se"){
      lambda.opt <- cv$lambda.1se
    }
    

    
    
    ## MODEL SELECTION COEFFICIENTS ##
    pred.mod <- glmnet(x = as.matrix(xsample), y = y, alpha = 1, family = fam, standardize = FALSE, weights=weight)
    lasso_coef <- predict(pred.mod,type = "coefficients",s = lambda.opt)[1:dim(xsample.d)[2],]
    #Collect the names of the variables with non-zero coefficients
    coef_select <- names(lasso_coef[lasso_coef != 0])[-1]
    
    #If select zero predictors, then fit a HT
    if(length(coef_select) == 0){
      
      if (messages) {
        message("No variables selected in the model selection stage.  Fitting a HT estimator.")  
      }
      
      if(var_est == TRUE){
        
        HT <- horvitzThompson(y = y, pi = pi, N = N, pi2 = pi2, var_est = TRUE, var_method = var_method)
        
        return(list( pop_total = HT$pop_total, 
                     pop_mean = HT$pop_total/N,
                     pop_total_var = HT$pop_total_var, 
                     pop_mean_var = HT$pop_total_var/N^2, 
                     weights = as.vector(pi^{-1})))
      } else {
        
        HT <- horvitzThompson(y = y, pi = pi, N = N, pi2 = pi2, var_est = FALSE)
        
        return(list( pop_total = HT$pop_total, 
                     pop_mean = HT$pop_total/N,
                     weights = as.vector(pi^{-1})))
       
        
      }
    }  
     
    
    
    #Create a new xsample with only the columns in coef_select
    xsample <- xsample[,coef_select, drop = FALSE]#dplyr::select_(xsample, .dots=coef_select)
    xsample.d <- model.matrix(~., data = xsample)
    xsample.dt <- t(xsample.d) 
    xsample <- data.frame(xsample.d[,-1, drop=FALSE])
    
  }
  
### LOGISTIC REGRESSION ###
  if (model == "logistic"){
    
    #Error if y has more than two categories (can only handle two right now)
    if(length(levels(as.factor(y))) != 2){
      stop("Function can only handle categorical response with two categories.")
    }
    
    if (datatype=="raw"){
      xpop <- data.frame(model.matrix(~., data = xpop))[,-1]
      #Make sure to only take the columns which are also in xsample
      xpop <- dplyr::select(xpop, one_of(colnames(xsample)))
      xpop_d <- model.matrix(~., data = xpop)
    }
    
    #Fit model using survey's svyglm since glm doesn't handle weights as we want
    dat <- data.frame(y, weight, xsample)
    colnames(dat) <- c("y", "weight", colnames(xsample))
    f <- paste(names(dat)[1], "~", paste(names(dat)[-c(1,2)], collapse=" + "))
    s.design<-survey::svydesign(ids=~1,weights=~weight,data=dat)
    mod <- survey::svyglm(f, design=s.design,family=quasibinomial())

    y.hats.U <- as.matrix(predict(mod,newdata=data.frame(xpop_d[,-1]),type="response",family=quasibinomial()))
    y.hats.s <- as.matrix(predict(mod,type="response",family=quasibinomial()))
    #Estimator of total
    t <- t(y-y.hats.s)%*%pi^(-1) + sum(y.hats.U)
    
    #Coefficients
    coefs <- mod$coefficients
    
    #No weights
    w <- NULL
    
    if(var_est==TRUE){
      if(var_method!="bootstrapSRS"){
        e <- y-y.hats.s
        varEst <- varMase(y = e,pi = pi,pi2 = pi2,method = var_method, N = N)
      
      }else if(var_method=="bootstrapSRS"){
  
        #Find bootstrap variance
 
        #Bootstrap total estimates
        dat <- data.frame(y, weight, xsample)
        colnames(dat) <- c("y", "weight", colnames(xsample))
        
       # system.time(t_boot <-  boot(data = dat, statistic = logisticGregt, R = B, xpopd = xpop_d, weights = pi^{-1}, parallel = "snow", cl=cluster))
        
        t_boot <-  boot(data = dat, statistic = logisticGregt, R = B, xpopd = xpop_d, parallel = "multicore", ncpus = 2)
    
       
        #Adjust for bias and without replacement sampling
        varEst <- var(t_boot$t)*n/(n-1)*(N-n)/(N-1)
      }
      return(list( pop_total = as.numeric(t),
                   pop_mean = as.numeric(t)/N,
                   pop_total_var=varEst, 
                   pop_mean_var=varEst/N^2,
                   coefficients =  coefs))
    }else{
      return(list( pop_total = as.numeric(t), 
                   pop_mean = as.numeric(t)/N,
                   coefficients =  coefs))     
      
    }
  }
### LINEAR REGRESSION ###
  else if (model == "linear"){
    
    #population design matrix, check whether its population totals, means or raw data
    
    if (datatype=="raw"){
      xpop <- data.frame(model.matrix(~.-1, data = data.frame(xpop)))
      #Make sure to only take the columns which are also in xsample
      xpop <- dplyr::select(xpop, one_of(colnames(xsample)))
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
    
    
    w <- as.matrix(1 + t(as.matrix(xpop_d) - xsample.dt %*% weight) %*% solve(xsample.dt %*% diag(weight) %*% xsample.d) %*% (xsample.dt)) %*% diag(weight)
  
  
  #calculating the total estimate for y
  t <- w %*% y
  

  #NOTE: check that weights times x's should equal total of x's to check for correct weight values
  
  #Coefficients
  coefs <- solve(xsample.dt %*% diag(weight) %*% xsample.d) %*% (xsample.dt) %*% diag(weight) %*% y
  
  
  if(var_est==TRUE){
    if(var_method!="bootstrapSRS"){
    y.hat <- xsample.d %*% (solve(xsample.dt %*% diag(weight) %*% xsample.d) %*% (xsample.dt) %*% diag(weight)%*%y)
    e <- y-y.hat
    varEst <- varMase(y = e,pi = pi, pi2 = pi2, method = var_method, N = N)
    
    }else if(var_method=="bootstrapSRS"){
      #Find bootstrap variance
      dat <- cbind(y,pi, xsample.d)
      #Bootstrap total estimates
      t_boot <- boot(data = dat, statistic = gregt, R = B, xpopd = xpop_d, parallel = "multicore", ncpus = 2)
      
      #Adjust for bias and without replacement sampling
      varEst <- var(t_boot$t)*n/(n-1)*(N-n)/(N-1)
    }
    return(list( pop_total = as.numeric(t), 
                 pop_mean = as.numeric(t)/N,
                 pop_total_var=varEst, 
                 pop_mean_var = varEst/N^2, 
                 weights = as.vector(w),
                 coefficients =  coefs))
  }else{
    return(list( pop_total = as.numeric(t), 
                 pop_mean = as.numeric(t)/N,           
                 weights = as.vector(w),
                 coefficients =  coefs))      
    
  }
  
  }
  
}

