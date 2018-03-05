#' Compute a generalized regression estimator
#' 
#' Calculates a generalized regression estimator for a finite population mean/proportion or total based on sample data collected from a complex sampling design and auxiliary population data.  
#' 
#' @param y  survey variable of interest. Should be a numeric vector or factor with two levels.
#' @param pi first order inclusion probabilities of same length as y. If not supplied, simple random sampling is assumed.
#' @param xsample dataframe of auxiliary data in the sample.  Number of rows should match length of y. 
#' @param xpop dataframe of population level auxiliary information.  Must contain same names as xsample.  Can supply either unit level data (denoted as datatype = "raw") or aggregate data (denoted as either datatype = "totals" or datatype = "means")
#' @param N population size, if not provided estimated to be the sum of the inverse inclusion probabilities
#' @param datatype form of population auxiliary data. Takes values "raw", "totals" or "means" for whether the user is providing population data at the unit level, aggregated to totals, or aggregated to means.
#' @param model regression model to utilize. Options are "linear" or "logistic".
#' @param var_est a logical indicating whether or not to compute a variance estimate
#' @param var_method method to use when computing the variance estimate.  Options are "HB"= Hajek-Berger estimator, "HH" = Hansen-Hurwitz estimator, "HTSRS" = Horvitz-Thompson estimator under simple random sampling, "HT" = Horvitz-Thompson estimator, "bootstrapSRS" = bootstrap variance estimator under simple random sampling without replacement.
#' @param B number of bootstrap samples if computing the bootstrap variance estimator.  Default is 1000.
#' @param pi2 a square matrix of the joint inclusion probabilities.  Needed for the "HT" variance estimator.
#' @param modelselect default to FALSE.  If TRUE, the Lasso is run and only the predictors with a non-zero Lasso coefficient are included in the model.
#' @param lambda Default to "lambda.min".  Only use if modelselect = TRUE. Takes values "lambda.min", which is the lambda value associated with the minimum cross validation error or "lambda.1se", which is the lamabda value which is one standard error away from the minimizing lambda and produces a sparser fit
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

greg  <- function(y, xsample, xpop, pi = NULL, model = "linear",  pi2 = NULL, var_est = FALSE, var_method="HB", datatype = "raw", N = NULL, modelselect = FALSE, lambda="lambda.min", B = 1000){

  
  
### INPUT VALIDATION ###
  
  #Check that y is numeric
  if(!(typeof(y) %in% c("numeric", "integer"))){
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
  
  if(model == "logistic" & datatype != "raw"){
    message("Must supply the raw population data to fit the logistic regression estimator.")
    return(NULL)
    
  }
  
  if(!is.element(datatype, c("raw","totals", "means"))){
    message("datatype input incorrect, has to be either \"raw\", \"totals\" or \"means\"")
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
    message("Assuming simple random sampling")
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
    if(length(coef_select)==0){
      message("No variables selected in the model selection stage.  Fitting a HT estimator.")
    
    if(var_est==TRUE){
      
      HT <- horvitzThompson(y = y, pi = pi, N = N, pi2 = pi2, var_est = TRUE, var_method = var_method)
      
      return(list( pop_total = HT$pop_total, 
                   pop_mean = HT$pop_total/N,
                   pop_total_var = HT$pop_total_var, 
                   pop_mean_var = HT$pop_total_var/N^2, 
                   weights = as.vector(pi^{-1})))
    }else {
      
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
    if(length(levels(as.factor(y)))!=2){
      message("Function can only handle categorical response with two categories.")
      return(NULL)
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
      xpop <- data.frame(model.matrix(~., data = xpop))[,-1]
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
    y.hat <- xsample.d%*%solve(xsample.dt %*% diag(weight) %*% xsample.d) %*% (xsample.dt) %*% diag(weight)%*%y
    e <- y-y.hat
    varEst <- varMase(y = e,pi = pi,pi2 = pi2,method = var_method, N = N)
    
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
                 pop_mean_var=varEst/N^2, 
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

