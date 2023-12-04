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
#' library(dplyr)
#' data(IdahoPop)
#' data(IdahoSamp)
#' 
#' xsample <- filter(IdahoSamp, COUNTYFIPS == "16055")
#' xpop <- filter(IdahoPop, COUNTYFIPS == "16055")
#' 
#' greg(y = xsample$BA_TPA_ADJ,
#'      N = xpop$npixels,
#'      xsample = xsample[c("tcc", "elev")],
#'      xpop = xpop[c("tcc", "elev")],
#'      var_est = TRUE,
#'      var_method = "LinHB",
#'      datatype = "means")
#' 
#' @references 
#' \insertRef{cas76}{mase}
#'
#' \insertRef{sar92}{mase}
#' 
#' @return 
#' A list of output containing:
#' 
#' * pop_total: Estimate of population total.
#' 
#' * pop_mean: Estimate of the population mean (or proportion).
#' 
#' * weights: Survey weights produced by GREG (linear model only).
#' 
#' * pop_total_var: Estimated variance of population total estimate.
#' 
#' * pop_mean_var: Estimated variance of population mean estimate.
#' 
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
  
  if (model == "logistic" & datatype != "raw") {
    stop("Must supply the raw population data to fit the logistic regression estimator.")
  }
  
  if (!is.element(datatype, c("raw","totals", "means"))) {
    stop("datatype input incorrect, has to be either \"raw\", \"totals\" or \"means\"")
  }
  
  if (is.null(N)) {
    if (datatype=="raw") {
      N <- dim(as.matrix(xpop))[1]
    } else {
      N <- sum(pi^(-1))
      if (messages) {
        message("Assuming N can be approximated by the sum of the inverse inclusion probabilities.") 
      }
    }
  }
  
  y <- as.numeric(y)
  
  #sample size
  n <- length(y)
  
  #create design matrix, x matrix and transpose design matrix
  xsample.d <- model.matrix(~., data = data.frame(xsample))
  xsample <- data.frame(xsample.d[,-1, drop = FALSE])
  xsample.dt <- t(xsample.d) 
 
  if (is.null(pi)) {
    if (messages) {
      message("Assuming simple random sampling") 
    }
  }  
  
  if (is.null(pi)) {
    pi <- rep(length(y)/N, length(y))
  }
  
  weight <- as.vector(pi^(-1))
  
  if (modelselect == TRUE) {
    
    #Cross-validation to find lambdas
    if(model == "linear"){
      fam <- "gaussian"
    } else{
      fam <- "binomial"
    } 
    
    #Run cv to find optimal lambda
    cv <- cv.glmnet(x = as.matrix(xsample), y = y, alpha = 1, weights = weight, nfolds = 10, family = fam, standardize = FALSE)
    
    #Pick lambda
    if (lambda=="lambda.min") {
      lambda.opt <- cv$lambda.min
    }
    if (lambda=="lambda.1se") {
      lambda.opt <- cv$lambda.1se
    }
  
    # model selection coefficients
    pred.mod <- glmnet(x = as.matrix(xsample), y = y, alpha = 1, family = fam, standardize = FALSE, weights=weight)
    lasso_coef <- predict(pred.mod,type = "coefficients",s = lambda.opt)[1:dim(xsample.d)[2],]
    #Collect the names of the variables with non-zero coefficients
    coef_select <- names(lasso_coef[lasso_coef != 0])[-1]
    
    #If select zero predictors, then fit a HT
    if(length(coef_select) == 0){
      
      if (messages) {
        message("No variables selected in the model selection stage.  Fitting a HT estimator.")  
      }
      
      if (var_est == TRUE) {
        
        HT <- horvitzThompson(y = y,
                              pi = pi,
                              N = N,
                              pi2 = pi2,
                              var_est = TRUE,
                              var_method = var_method,
                              fpc = fpc)
        
        return(list(pop_total = HT$pop_total, 
                    pop_mean = HT$pop_total/N,
                    pop_total_var = HT$pop_total_var, 
                    pop_mean_var = HT$pop_total_var/N^2, 
                    weights = as.vector(pi^{-1}),
                    estimator_used = "HT"))
      } else {
        
        HT <- horvitzThompson(y = y,
                              pi = pi,
                              N = N,
                              pi2 = pi2,
                              var_est = FALSE)
        
        return(list(pop_total = HT$pop_total, 
                    pop_mean = HT$pop_total/N,
                    weights = as.vector(pi^{-1}),
                    estimator_used = "HT"))
        
      }
    }  
     
    #Create a new xsample with only the columns in coef_select
    xsample <- xsample[ ,coef_select, drop = FALSE]
    xsample.d <- model.matrix(~., data = xsample)
    xsample.dt <- t(xsample.d) 
    xsample <- data.frame(xsample.d[ ,-1, drop = FALSE])
    
  }
  
  if (model == "logistic"){
    
    #Error if y has more than two categories (can only handle two right now)
    if (length(levels(as.factor(y))) != 2) {
      stop("Function can only handle categorical response with two categories.")
    }
    
    if (datatype == "raw") {
      xpop <- data.frame(model.matrix(~., data = xpop))[,-1]
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
    
    t <- t(y-y.hats.s)%*%pi^(-1) + sum(y.hats.U)
    
    #Coefficients
    coefs <- mod$coefficients
    
    #No weights
    w <- NULL
    
    if (var_est==TRUE) {
      
      if (var_method!="bootstrapSRS") {
        e <- y - y.hats.s
        varEst <- varMase(y = e,
                          pi = pi,
                          pi2 = pi2,
                          method = var_method,
                          N = N,
                          fpc = fpc)
      } else if (var_method == "bootstrapSRS") {
  
        #Bootstrap total estimates
        dat <- data.frame(y, weight, xsample)
        colnames(dat) <- c("y", "weight", colnames(xsample))
        t_boot <-  boot(data = dat,
                        statistic = logisticGregt,
                        R = B, 
                        xpopd = xpop_d,
                        parallel = "multicore",
                        ncpus = 2)
    
        if (fpc == T) {
          varEst <- var(t_boot$t)*n/(n-1)*(N-n)/(N-1) 
        }
        if (fpc == F) {
          varEst <- var(t_boot$t)*n/(n-1)
        }
        
      }
      return(list(pop_total = as.numeric(t),
                  pop_mean = as.numeric(t)/N,
                  pop_total_var = as.numeric(varEst), 
                  pop_mean_var = as.numeric(varEst)/N^2,
                  coefficients =  coefs))
    } else {
      return(list(pop_total = as.numeric(t), 
                  pop_mean = as.numeric(t)/N,
                  coefficients =  coefs))     
      
    }
  }

  else if (model == "linear"){
    
    if (datatype == "raw") {
      xpop <- data.frame(model.matrix(~.-1, data = data.frame(xpop)))
      xpop <- dplyr::select(xpop, one_of(colnames(xsample)))
      xpop_d <- model.matrix(~., data = xpop)
      xpop_d <- apply(xpop_d, 2, sum)
    }
    if (datatype == "totals") {
      xpop_d <- unlist(c(N, xpop[names(xsample)]))
    }
    if (datatype == "means") {
      xpop_d <- unlist(c(N, xpop[names(xsample)]*N))
    }
    
    one_mat <- matrix(rep(1, times = nrow(xsample.d)), nrow = 1)
    xpop_cpp <- as.matrix(xpop_d)
    weight_mat <- diag(weight)
    
    w <- get_weights_greg(xpop_cpp, xsample.d, weight_mat, one_mat)

  #calculating the total estimate for y
  t <- sum(as.numeric(w) * y)
  # t <- w %*% y

  #Coefficients
  coefs <- get_coefs(xsample.d, as.vector(y), weight_mat)
  names(coefs) <- colnames(xsample.d)
  
  if (var_est == TRUE) {
    
    if (var_method!="bootstrapSRS") {
      y.hat <- xsample.d %*% coefs
      e <- y - y.hat
      varEst <- varMase(y = e,
                        pi = pi,
                        pi2 = pi2,
                        method = var_method,
                        N = N,
                        fpc = fpc)
    } else if (var_method == "bootstrapSRS") {
      #Find bootstrap variance
      dat <- cbind(y,pi, xsample.d)
      #Bootstrap total estimates
      t_boot <- boot(data = dat,
                     statistic = gregt,
                     R = B,
                     xpopd = xpop_d,
                     parallel = "multicore",
                     ncpus = 2)
      
      if (fpc == T) {
        varEst <- var(t_boot$t)*n/(n-1)*(N-n)/(N-1)
      }
      if (fpc == F) {
        varEst <- var(t_boot$t)*n/(n-1)
      }

    }
    return(list(pop_total = as.numeric(t), 
                pop_mean = as.numeric(t)/N,
                pop_total_var=varEst, 
                pop_mean_var = varEst/N^2, 
                weights = as.vector(w),
                coefficients =  coefs))
  } else {
    return(list(pop_total = as.numeric(t), 
                pop_mean = as.numeric(t)/N,           
                weights = as.vector(w),
                coefficients =  coefs))      
    
  }
  
  }
  
}

