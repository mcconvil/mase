#' Compute a regression tree estimator
#' 
#' Calculates a regression tree estimator for a finite population mean or total based on sample data collected from a complex sampling design and auxiliary population data.  
#' 
#' @inheritParams horvitzThompson
#' @param x_sample A data frame of the auxiliary data in the sample.
#' @param x_pop A data frame of population level auxiliary information. It must contain the same names as x_sample.
#' @param p_value Designated p-value level to reject null hypothesis in permutation test used to fit the regression tree. Default value is 0.05.
#' @param perm_reps An integer specifying the number of permutations for each permutation test run to fit the regression tree. Default value is 500.
#' @param bin_size An integer specifying the minimum number of observations in each node.
#' @param ncpus Integer denoting the number of cpu cores to use in parallel for bootstrap variance.
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
#' \item{weights: Survey weights produced by regression tree}
#' \item{formula: Model formula}
#' \item{tree: rpms object}
#' \item{y_hat_sample: Response estimates for sample}
#' \item{y_hat_pop: Response estimates for population}
#' }
#' 
#' @export gregTree
#' @import rpms
#' @import boot
#' @importFrom stats as.formula
#' @include varMase.R
#' @include gregt.R
#' @include gregObject.R
#' 
#' @seealso \code{\link{greg}} for a linear or logistic regression model.


gregTree  <- function(y, x_sample, x_pop, pi = NULL,  pi2 = NULL, var_est = FALSE,
                      var_method="lin_HB", B = 1000, p_value = 0.05, perm_reps = 500,
                      bin_size = NULL, strata = NULL, ncpus = 1){

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

  #Fit model
  dat <- data.frame(y, x_sample, weights = weights)
  #Create formula for rpms equation
  f <- as.formula(paste("y ~ ", paste(names(x_sample), collapse= "+")))
  tree <- rpms(rp_equ = f, data = dat, weights = weights, pval = p_value,
               perm_reps = perm_reps, bin_size = bin_size)
  
  
  #Calculate weights
  #Make sure x_pop and x_sample have the same columns (in same order)
  x_pop <- x_pop[names(x_sample)]
  #Design matrix for population
  x_pop_tree <- treeDesignMatrix(splits = tree$ln_split, data = x_pop)
  #Design matrix for sample
  x_sample_tree <- treeDesignMatrix(splits = tree$ln_split, data = x_sample)
  
  #weights
  w <- (1 + t(colSums(x_pop_tree)- colSums(x_sample_tree*pi^(-1)))%*%solve(t(x_sample_tree)%*%diag(pi^(-1))%*%as.matrix(x_sample_tree))%*%t(x_sample_tree))%*%diag(pi^(-1)) 
  
  #calculating the total estimate for y
  t <- w %*% y
  #find 
  y_hat_sample <- predict(object = tree, newdata = x_sample)
  y_hat_pop <- predict(object = tree, newdata = x_pop)

  #NOTE: check that weights times x's should equal total of x's to check for correct weight values
  # w %*% x_sample_tree[,4]
  # colSums(x_pop_tree)
  

  if(var_est==TRUE){
    if(var_method != "bootstrap_SRS"){
    e <- y - y_hat_sample
    varEst <- varMase(y = e, pi = pi, pi2 = pi2, method = var_method, N = N,
                      strata = strata)
    }
    else if(var_method == "bootstrap_SRS"){
      #Find bootstrap variance
      dat <- cbind(y, pi, x_sample)
      #Bootstrap total estimates
      t_boot <- boot(data = dat, statistic = gregTreet, R = B, x_pop = x_pop,
                     p_value= p_value, perm_reps = perm_reps, bin_size = bin_size,
                     parallel = "multicore", ncpus = ncpus)

      #Adjust for bias and without replacement sampling
      varEst <- var(t_boot$t)*n/(n-1)*(N-n)/(N-1)
    }
    return(list( pop_total = as.numeric(t),
                 pop_mean = as.numeric(t)/N,
                 pop_total_var=varEst,
                 pop_mean_var=varEst/N^2,
                 weights = as.vector(w),
                 formula = f,
                 tree = tree,
                 y_hat_sample = y_hat_sample,
                 y_hat_pop = y_hat_pop) %>%
           gregify())
  }
  else{
    return(list( pop_total = as.numeric(t),
                 pop_mean = as.numeric(t)/N,
                 weights = as.vector(w), 
                 formula = f,
                 tree = tree,
                 y_hat_sample = y_hat_sample,
                 y_hat_pop = y_hat_pop) %>%
             gregify())
    }
  }



