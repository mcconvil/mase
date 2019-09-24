#' Compute a neural network model assisted estimator
#' 
#' Calculates a dense feed forward neural network model assisted estimator for a finite population mean or total based on sample data collected from a complex sampling design and auxiliary population data. Requires installation of keras.  
#' 
#' @inheritParams horvitzThompson
#' @param x_sample A data frame of the auxiliary data in the sample.
#' @param x_pop A data frame of population level auxiliary information. It must contain the same names as x_sample.
#' @param standardize A logical indicating whether or not to standardize predictors before making model.
#' @param lambda scale parameter of gaussian kernel (unnormalized). If NULL, this will be numerically optimized until the condition number for interpolation matrix is near 10e12.
#' @param PCA If TRUE, principal component analysis will be done and the first two principal components will be used to support the RBF.
#' @param bag If TRUE, the weights will be bagged.
#' @param B Number of times to bag the weights
#' @examples
#' library(survey)
#' data(api)
#' gregNeuralNet(y = apisrs$api00, 
#' x_sample = apisrs[c("col.grad", "awards", "snum", "dnum", "cnum", "pcttest", "meals", "sch.wide")], 
#' x_pop = apipop[c("col.grad", "awards", "snum", "dnum", "cnum", "pcttest", "meals", "sch.wide")])
#' 
#'@references 
#'\insertRef{mcc17b}{mase}

#' 
#' @return A greg object containing:
#' \itemize{
#' \item{pop_total: Estimate of population total}
#' \item{pop_mean: Estimate of the population mean}
#' \item{pop_total_var: Estimated variance of population total estimate}
#' \item{pop_mean_var: Estimated variance of population mean estimate}
#' \item{formula: Model formula}
#' \item{model: feed forward neural network model}
#' \item{y_hat_sample: Response estimates for sample}
#' \item{y_hat_pop: Response estimates for population}
#' }
#' 
#' @export gregRBF
#' @import dplyr
#' @importFrom stats as.formula
#' @include varMase.R
#' @include gregt.R
#' @include gregObject.R
#' 
#' @seealso \code{\link{greg}} for a linear or logistic regression model.


gregRBF <- function(y, x_sample, x_pop, pi = NULL,  pi2 = NULL, var_est = FALSE,
                    var_method="lin_HB", strata = NULL, lambda = 1.666, PCA = FALSE,
                    standardize = FALSE, bag = FALSE, B = 100){
  #Make sure the var_method is valid
  if(!is.element(var_method, c("lin_HB", "lin_HH", "lin_HTSRS", "lin_HT", "bootstrap_SRS"))){
    message("Variance method input incorrect. It has to be \"lin_HB\", \"lin_HH\", \"lin_HT\", \"lin_HTSRS\", or \"bootstrap_SRS\".")
    return(NULL)
  }
  #Check that x_sample is a df
  if(class(x_sample) != "data.frame"){
    message("x_sample must be a data.frame.")
    return(NULL)
  }
  #Check that y is numeric
  if(!(typeof(y) %in% c("numeric", "integer", "double"))){
    stop("Must supply numeric y.  For binary variable, convert to 0/1's.")
  }
  #Convert y to a vector
  y <- as.vector(y)
  
  #Make sure y is complete
  if(NA %in% y){
    message("Must supply complete cases for y")
    return(NULL)
  }
  #Check that the number of observations are consistent
  if(nrow(x_sample) != length(y)){
    message("y and x_sample must have same number of observations.")
    return(NULL)
  }
  #sample size
  n <- length(y)
  
  #population size
  N <- dim(x_pop)[1]
  
  #Check on inclusion probabilities and create weight=inverse inclusion probabilities
  if(is.null(pi)){
    message("Assuming simple random sampling")
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
  # Create pi
  if (is.null(pi)) {
    pi <- rep(length(y)/N, length(y))
  }
  
  #weight: inverse first order inclusion probabilities
  weights <- as.vector(pi^(-1))
  
  #Make sure x_pop and x_sample have the same columns (in same order)
  x_pop <- x_pop[names(x_sample)]
  
  #Check standardization
  if(standardize == TRUE){
    cent <- colMeans(x_sample)
    scl <- apply(as.matrix(x_sample), 2, sd)
    if(0 %in% scl) {
      message("can not standardize a variable with zero variance.")
      return(NULL)
    }
    x_pop <- base::scale(x_pop, center = cent,
                         scale = scl) %>%
      as.data.frame()
    x_sample <- base::scale(x_sample) %>% as.data.frame()
    standardize <- list(center = cent, scale = scl)
  }
  
  #Create formula and model matrix
  f <- as.formula(paste("y ~ ", paste(names(x_sample), collapse= "+"), "+ 0"))
  
  #Do PCA, fight off the curse of dimensionality :)
  if(PCA == TRUE) {
    pca_obj <- princomp(~., x_sample)
    #use first two principal components
    model_mat <- pca_obj$scores[,1:2]
  }
  else{ #use everything
    pca_obj <- NULL
    model_mat <- model.matrix(f, x_sample)
  }
  
  #Make Radial Basis Function
  
  G0 <- model_mat %>% #make distance matrix
    dist(diag = T, upper = T) %>%
    as.matrix()
  #numerically tune lambda
  if(is.null(lambda)){
    find_lambda <- function(lambda) {
      G <- exp(-(lambda*G0)^2)
      (kappa(G) - 10e12)^2 #10e12 is the "suggested" condition number for interpolation matrix
    }
    lambda <- optim(5, find_lambda, lower = 0, upper = 200, method = "Brent")$par
  }
  #make final interpolation matrix
  G <- exp(-(lambda*G0)^2)
  rbf_weights <- as.vector(solve(G) %*% y) #find rbf_weights
  #bag
  if(bag == TRUE){
    smp <- round(2/3*length(y))
    bagged_weights <- lapply(1:B, function(i){
      #sample 2/3 in bag
      samp <- sample(rep(c(NA, 1), c(length(y) - smp, smp)))
      ind <- as.vector(na.omit(1:length(y)*samp))
      y_bag <- y[ind]
      #apply kernel
      G <- exp(-(lambda*(G0[ind, ind]))^2)
      return(as.vector(solve(G) %*% y_bag) * samp) #return as length y vector
    })
    rbf_weights <- bagged_weights %>% Reduce(f = "rbind") %>% colMeans(na.rm = TRUE)
  }
  
  #create rbf model
  phi <- function(x) { # x is a vector of predictors for obs i
    centers <- model_mat
    z <- (t(centers) - x) %>% apply(2, function(col_vec) sqrt(sum(col_vec^2)))
    sum(rbf_weights * exp(-(lambda*z)^2))
  }
  rbf_obj <- (list(phi = function(x) phi(x), rbf_weights = rbf_weights, lambda = lambda,
                   pca_obj = pca_obj, condition_number = kappa(G)))
  class(rbf_obj) <- "greg_rbf"
  ###find estimates
  y_hat_samp <- sapply(1:length(y), function(i) phi(model_mat[i,]))
  message("Fitting population values")
  y_hat_pop <- predict(rbf_obj, x_pop)
  
  #calculating the total estimate for y
  t <- sum(weights * (y - y_hat_samp)) + sum(y_hat_pop)
  
  if(var_est==TRUE){
    if(var_method != "bootstrap_SRS"){
      e <- y - y_hat_samp
      varEst <- varMase(y = e, pi = pi, pi2 = pi2, method = var_method, N = N,
                        strata = strata)
    }
    return(list( pop_total = as.numeric(t),
                 pop_mean = as.numeric(t)/N,
                 pop_total_var=varEst,
                 pop_mean_var=varEst/N^2,
                 formula = f,
                 model = rbf_obj,
                 y_hat_sample = y_hat_samp,
                 y_hat_pop = y_hat_pop,
                 standardize = standardize) %>%
             gregify())
  }
  else{
    return(list( pop_total = as.numeric(t),
                 pop_mean = as.numeric(t)/N,
                 formula = f,
                 model = rbf_obj,
                 y_hat_sample = y_hat_samp,
                 y_hat_pop = y_hat_pop,
                 standardize = standardize) %>%
             gregify())
  }
}

predict.greg_rbf <- function(obj, newdata) {
  newdata <- as.matrix(newdata)
  if(!is.null(obj$model$pca_obj)){ 
    #If PCA used, project onto first two loadings
    newdata <- (newdata %*% pca_obj$loadings)[,1:2]
  }
  sapply(1:nrow(newdata),function(i) obj$phi(x = newdata[i,]))
}