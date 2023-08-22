#' Compute a regression tree estimator
#' 
#' Calculates a regression tree estimator for a finite population mean/proportion or total based on sample data collected from a complex sampling design and auxiliary population data.  
#' 
#' @inheritParams horvitzThompson
#' @inheritParams greg
#' @param pval Designated p-value level to reject null hypothesis in permutation test used to fit the regression tree. Default value is 0.05.
#' @param perm_reps An integer specifying the number of permutations for each permutation test run to fit the regression tree. Default value is 500.
#' @param bin_size A integer specifying the minimum number of observations in each node.
#' 
#' @examples
#' data(IdahoPop)
#' data(IdahoSamp)
#' 
#' xsample <- select(IdahoSamp, tcc, elev) %>% as.data.frame()
#' xpop <- select(IdahoPop, names(xsample))
#' 
#' gregTree(y = xsample$BA_TPA_ADJ,
#'      N = xpop$npixels,
#'      xsample = xsample,
#'      xpop = xpop,
#'      var_est = TRUE,
#'      datatype = "means")
#' 
#'@references 
#'\insertRef{mcc17b}{mase}

#' 
#' @return A list of output containing:
#' \itemize{
#' \item{pop_total:}{Estimate of population total}
#' \item{pop_mean:}{Estimate of the population mean (or proportion)}
#' \item{weights:}{Survey weights produced by gregTree}
#' \item{pop_total_var:}{Estimated variance of population total estimate}
#' \item{pop_mean_var:}{Estimated variance of population mean estimate}
#' }
#' 
#' @export gregTree
#' @import rpms
#' @import boot
#' @importFrom stats as.formula
#' @include varMase.R
#' @include gregt.R

gregTree  <- function(y,
                      xsample,
                      xpop,
                      pi = NULL,
                      pi2 = NULL,
                      var_est = FALSE, 
                      var_method="LinHB",
                      B = 1000,
                      pval = 0.05,
                      perm_reps = 500,
                      bin_size = NULL,
                      fpc = T,
                      messages = T){

  if (!is.element(var_method, c("LinHB", "LinHH", "LinHTSRS", "LinHT", "bootstrapSRS"))) {
    stop("Variance method input incorrect. It has to be \"LinHB\", \"LinHH\", \"LinHT\", \"LinHTSRS\", or \"bootstrapSRS\".")
  }
  
  if (!(typeof(y) %in% c("numeric", "integer", "double"))) {
    stop("Must supply numeric y.  For binary variable, convert to 0/1's.")
  }
  
  #Convert y to a vector
  y <- as.vector(y)

  #sample size
  n <- length(y)

  #population size
  N <- dim(xpop)[1]

  #Check on inclusion probabilities and create weight=inverse inclusion probabilities
  if (is.null(pi)) {
    if (messages) {
      message("Assuming simple random sampling") 
    }
  }

  # convert pi into diagonal matrix format
  if (is.null(pi)) {
    pi <- rep(length(y)/N, length(y))
  }

  #weight: inverse first order inclusion probabilities
  weights <- as.vector(pi^(-1))

  #Fit model
  dat <- data.frame(y, xsample, weights = weights)
  #Create formula for rpms equation
  f <- as.formula(paste("y ~ ", paste(names(xsample), collapse= "+")))
  tree <- rpms(rp_equ = f, data = dat, weights = weights, pval = pval, perm_reps = perm_reps, bin_size = bin_size)
  
  #Make sure xpop and xsample have the same columns (in same order)
  xpop <- xpop[names(xsample)]
  #Design matrix for population
  xpop_tree <- box_ind(tree, xpop)
  #Design matrix for sample
  xsample_tree <- box_ind(tree, xsample)
  w <- (1 + t(colSums(xpop_tree)- colSums(xsample_tree*pi^(-1)))%*%solve(t(xsample_tree)%*%diag(pi^(-1))%*%as.matrix(xsample_tree))%*%t(xsample_tree))%*%diag(pi^(-1)) 

  #calculating the total estimate for y
  t <- w %*% y

  if (var_est == TRUE) {
    
    if (var_method != "bootstrapSRS") { 
      y.hat <- predict(object = tree, newdata = xsample)
      e <- y-y.hat
      varEst <- varMase(y = e,pi = pi,pi2 = pi2,method = var_method, N = N, fpc = fpc)
    } else if (var_method == "bootstrapSRS") {
      #Find bootstrap variance
      dat <- cbind(y, pi, xsample)
      #Bootstrap total estimates
      t_boot <- boot(data = dat,
                     statistic = gregTreet, 
                     R = B,
                     xpop = xpop,
                     pval= pval,
                     perm_reps = perm_reps, 
                     bin_size = bin_size,
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
                pop_mean_var=varEst/N^2,
                weights = as.vector(w),
                tree = tree))
  } else {
    
    return(list(pop_total = as.numeric(t),
                pop_mean = as.numeric(t)/N,
                weights = as.vector(w), 
                tree = tree))

  }

}



