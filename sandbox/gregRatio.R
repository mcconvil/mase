#' Compute an estimator of a ratio with a ratio of generalized regression estimators
#' 
#' Calculates an estimator for a finite population ratio of totals where the numerator and denominator are estimated by generalized regression estimators. The estimator is 
#' based on sample data collected from a complex sampling design and a single auxiliary population variable.  
#' 
#' @param y_num  A numeric vector of the sampled response variable for the numerator estimator.
#' @param y_den  A numeric vector of the sampled response variable for the denominator estimator.
#' @param xsample A data frame of the auxiliary data in the sample for the numerator and denominator estimators.
#' @param xpop A data frame of population level auxiliary information for the numerator and denominator estimators.  It must contain the same names as xsample.  If datatype = "raw", must contain unit level data.  If datatype = "totals" or "means", then contains one row of aggregated, population totals or means for the auxiliary data. Default is "raw".
#' @param datatype A string that specifies the form of population auxiliary datar. The possible values are "raw", "totals" or "means" for whether the user is providing population data at the unit level, aggregated to totals, or aggregated to means.  Default is "raw".
#' @param model A string that specifies the regression model to utilize for the numerator and denominator estimators. Options are "linear" or "logistic".
#' @param pi A numeric vector of inclusion probabilities for each sampled unit.  If NULL, then simple random sampling without replacement is assumed.
#' @param N A numeric value of the population size. If NULL, it is estimated with the sum of the inverse of the pis.
#' @param var_est A logical indicating whether or not to compute a variance estimator.  Default is FALSE.
#' @param var_method A string that specifies the method to use when computing the variance estimator.  Currently the options are the following Taylor linearized techniques: "LinHTgWeighted" = Horvitz-Thompson estimator with g weights and "LinHT" = Horvitz-Thompson estimator.
#' @param pi2 A square matrix of the joint inclusion probabilities.  Needed for the "LinHT" variance estimator.
#' 
#' @examples 
#' 
#'@references 
#'
#'\insertRef{sar92}{mase}
#' 
#' @return A list of output containing:
#' \itemize{
#' \item{pop_ratio:}{Estimate of population ratio of totals}
#' \item{pop_ratio_var:}{Estimated variance of estimate of ratio}
#' }
#' 
#' @export gregRatio
#' @import survey
#' @import boot
#' @importFrom stats model.matrix predict quasibinomial var
#' @include varMase.R
#' @include gregt.R

gregRatio  <- function(y_num, y_den, xsample, xpop, pi = NULL, model = "linear", 
                       pi2 = NULL, var_est = FALSE, datatype = "raw", N = NULL,
                       var_method = "HTLin"){

  # Define variables
  # Sample size
  n <- length(y_num)
  
### INPUT VALIDATION ###
  
  # Check that y_num and y_den have same length
  # Proxy for checking that they are from the same sample
  if(length(y_num) != length(y_den)){
    message("The response variables must be from the same sample and so must have the same length.")
    return(NULL)
  }
  
  
  #Check that y is numeric
  if(!(typeof(y_num) %in% c("numeric", "integer", "double"))){
    stop("Must supply numeric y_num.  For binary variable, convert to 0/1's.")
  }
  if(!(typeof(y_den) %in% c("numeric", "integer", "double"))){
    stop("Must supply numeric y_den.  For binary variable, convert to 0/1's.")
  }
  
  # #Make sure the var_method is valid
  # if(!is.element(var_method, c("LinHB", "LinHH", "LinHTSRS", "LinHT", "bootstrapSRS"))){
  #   message("Variance method input incorrect. It has to be \"LinHB\", \"LinHH\", \"LinHT\", \"LinHTSRS\", or \"bootstrapSRS\".")
  #   return(NULL)
  # }
  # 
  
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
  y_num <- as.vector(y_num)
  y_den <- as.vector(y_den)
  

  
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
    pi <- rep(n/N, n)
  }
  
  #weight: inverse first order inclusion probabilities
  weight <- as.vector(pi^(-1))
  
  
  
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
        varEst <- varMase(y = e,pi = pi,pi2 = pi2, method = var_method, N = N)
      
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
  
  
  #calculating the total estimates for y's
  t_num <- w %*% y_num
  t_den <- w %*% y_den
  
  # Estimate of the ratio
  rat <- t_num/t_den

  #NOTE: check that weights times x's should equal total of x's to check for correct weight values
  
  #Coefficients
  coefs_num <- solve(xsample.dt %*% diag(weight) %*% xsample.d) %*% (xsample.dt) %*% diag(weight) %*% y_num
  coefs_den <- solve(xsample.dt %*% diag(weight) %*% xsample.d) %*% (xsample.dt) %*% diag(weight) %*% y_den
  
  if(var_est == TRUE){
#    if(var_method!="bootstrapSRS"){
    y_hat_num <- xsample.d%*%solve(xsample.dt %*% diag(weight) %*% xsample.d) %*% (xsample.dt) %*% diag(weight)%*%y_num
    y_hat_den <- xsample.d%*%solve(xsample.dt %*% diag(weight) %*% xsample.d) %*% (xsample.dt) %*% diag(weight)%*%y_den
    e_num <- y_num - y_hat_num
    e_den <- y_den - y_hat_den    

    if(var_method == "LinHT"){
      e_ratio <- e_num - as.vector(rat)*e_den
      a <- (pi2 - pi%*%t(pi))*pi2^(-1)*(pi%*%t(pi))^(-1)*e_ratio%*%t(e_ratio)
      varEst <- sum(a)*t_den^(-2)
    }
    
    if(var_method == "LinHTgWeighted"){
      g <- as.vector(w) * pi
      e_ratio <- e_num - as.vector(rat)*e_den
      a <- (pi2 - pi%*%t(pi))*pi2^(-1)*(pi%*%t(pi))^(-1)*(g*e_ratio)%*%t((g*e_ratio))
      varEst <- sum(a)*t_den^(-2)   
    }
    
 #   }
    
    # if(var_method=="bootstrapSRS"){
    #   #Find bootstrap variance
    #   dat <- cbind(y,pi, xsample.d)
    #   #Bootstrap total estimates
    #   t_boot <- boot(data = dat, statistic = gregt, R = B, xpopd = xpop_d, parallel = "multicore", ncpus = 2)
    #   
    #   #Adjust for bias and without replacement sampling
    #   varEst <- var(t_boot$t)*n/(n-1)*(N-n)/(N-1)
    # }
    return(list( pop_ratio = rat, 
                 pop_ratio_var = varEst))
  }else{
    return(list( pop_ratio = rat))      
    
  }
  
  }
  
}

