#' Compute a random forest estimator
#' 
#' Calculates a random forest estimator for a finite population mean or total based on sample data collected from a complex sampling design and auxiliary population data.  
#' 
#' @inheritParams horvitzThompson
#' @param x_sample A data frame of the auxiliary data in the sample.
#' @param x_pop A data frame of population level auxiliary information. It must contain the same names as x_sample.
#' @param p_value Designated p-value level to reject null hypothesis in permutation test used to fit the regression tree. Default value is 0.05.
#' @param perm_reps An integer specifying the number of permutations for each permutation test run to fit the regression tree. Default value is 500.
#' @param bin_size An integer specifying the minimum number of observations in each node.
#' @param mtry A number specifying how many predictors each tree uses. Default is the square root of total number of predictors.
#' @param ntrees An integer specifying the number of trees to train in forest.
#' @param forest_fun A string specifying which randomForest algorithm to use. Choices are "rpms" or "randomForest", default is mase modified rpms random forest.
#' @param oob_pred Logical. If TRUE, make predictions using only trees out of bag.
#' @param cores An integer specifying the number of cores to use in parallel if > 1 (not implemented)

#' 
#' @examples
#' library(survey)
#' data(api)
#' gregTree(y = apisrs$api00, 
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
#' \item{forest: rpms random forest object}
#' \item{oob_error: Out of bag error, average tree out of bag MSE or misclassification rate}
#' \item{y_hat_sample: Response estimates for sample}
#' \item{y_hat_pop: Response estimates for population}
#' }
#' 
#' @export gregForest
#' @import dplyr
#' @import rpms
#' @import boot
#' @importFrom stats as.formula
#' @include varMase.R
#' @include gregt.R
#' @include gregObject.R
#' @include rpmsForestt.R
#' 
#' @seealso \code{\link{greg}} for a linear or logistic regression model.


gregForest <- function(y, x_sample, x_pop, pi = NULL,  pi2 = NULL, var_est = FALSE,
                      var_method="lin_HB", B = 1000, p_value = 0.05, perm_reps = 500,
                      bin_size = NULL, strata = NULL, mtry = NULL, ntrees = 100,
                      forest_fun = "default", oob_pred = TRUE, cores = 1){
    
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
    if(is.null(mtry)){
      mtry <- round(sqrt(ncol(x_pop)))
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
    
    #Fit model
    dat <- data.frame(y, x_sample, weights = weights)
    #Create formula for rpms equation
    f <- as.formula(paste("y ~ ", paste(names(x_sample), collapse= "+")))
    if(tolower(forest_fun) == "rpms"){ #use original rpms
      forest <- rpms::rpms_forest(rp_equ = f, data = dat, weights = weights,
                            pval = p_value, perm_reps = perm_reps, 
                            bin_size = bin_size, f_size = ntrees, cores = cores)
      y_hat_sample <- predict(obj = forest, newdata = x_sample)
      y_hat_pop <- predict(obj = forest, newdata = x_pop)
    }
    else if(tolower(forest_fun) == "randomforest"){
      message("using randomForest")
      forest <- randomForest::randomForest(formula = f, data = dat, nodesize = bin_size, mtry = mtry, ntree = ntrees,
                            cores = cores)
      y_hat_sample <- predict(obj = forest, newdata = x_sample)
      y_hat_pop <- predict(obj = forest, newdata = x_pop)
    }
    else{ #use mase edited rpms forest function
      forest <- rpmsForestt(rp_equ = f, data = dat, weights = weights, pval = p_value,
                            perm_reps = perm_reps, bin_size = bin_size, mtry = mtry, f_size = ntrees,
                            cores = cores)
      y_hat_sample <- predict(obj = forest, newdata = x_sample, oob = oob_pred)
      y_hat_pop <- predict(obj = forest, newdata = x_pop)
    }
    #calculating the total estimate for y
    t <- sum(weights * (y - y_hat_sample)) + sum(y_hat_pop)
    
    if(var_est==TRUE){
      if(var_method != "bootstrap_SRS"){
        e <- y - y_hat_sample
        varEst <- varMase(y = e, pi = pi, pi2 = pi2, method = var_method, N = N,
                          strata = strata)
      }
      return(list( pop_total = as.numeric(t),
                   pop_mean = as.numeric(t)/N,
                   pop_total_var=varEst,
                   pop_mean_var=varEst/N^2,
                   formula = f,
                   forest = forest,
                   oob_error = forest$oob_error,
                   var_imp = forest$importance,
                   y_hat_sample = y_hat_sample,
                   y_hat_pop = y_hat_pop) %>%
               gregify())
    }
    else{
      return(list( pop_total = as.numeric(t),
                   pop_mean = as.numeric(t)/N,
                   formula = f,
                   forest = forest,
                   oob_error = forest$oob_error,
                   var_imp = forest$importance,
                   y_hat_sample = y_hat_sample,
                   y_hat_pop = y_hat_pop) %>%
               gregify())
    }
  }