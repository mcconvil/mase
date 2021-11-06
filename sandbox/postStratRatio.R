#' Compute an estimator of a ratio with a ratio of post-stratified estimators
#' 
#' Calculates an estimator for a finite population ratio of totals where the numerator and denominator are estimated by post-stratified estimators. The estimator is 
#' based on sample data collected from a complex sampling design and a single auxiliary population variable.  
#' 
#' @param y_num  A numeric vector of the sampled response variable for the numerator of the ratio.
#' @param y_den A numeric vector of the sampled response variable for the denominator of the ratio.
#' @param xsample A vector containing the post-stratum for each sampled unit.
#' @param xpop A vector or data frame, depending on datatype.  If datatype = "raw", then a vector containing the post-stratum for each population unit.  If datatype = "totals" or "means", then a data frame, where the first column lists the possible post-strata and the second column contains the population total or proportion in each post-stratum.  
#' @param datatype Default to "raw", takes values "raw", "totals" or "means" for whether the user is providing the raw population stratum memberships, the population totals of each stratum, or the population proportions of each stratum.
#' @param pi A numeric vector of inclusion probabilities for each sampled unit.  If NULL, then simple random sampling without replacement is assumed.
#' @param N A numeric value of the population size. If NULL, it is estimated with the sum of the inverse of the pis.
#' @param var_est Default to FALSE, logical for whether or not to compute estimate of variance
#' @param var_method A string that specifies the method to use when computing the variance estimator.  Currently the only option is "SRSunconditional" = simple random sampling variance estimator which accounts for random strata.
#' @param fpc Default to TRUE, logical for whether or not the variance calculation should include a finite population correction.
#' 
#' @examples 
#' 
#'@references 
#'\insertRef{coc77}{mase} 
#' 
#'\insertRef{sar92}{mase}
#'
#' @return A list of output containing:
#' \itemize{
#' \item{pop_ratio:}{Estimate of population ratio of totals}
#' \item{pop_total_var:}{Estimated variance of population ratio of totals estimate}
#' }

#' @import dplyr
#' @import magrittr
#' @export postStratRatio
#' @include postStratt.R
#' @include varMase.R


postStratRatio <- function(
  y_num, y_den, xsample, xpop, pi = NULL, N = NULL, var_est = FALSE, var_method = "SRSunconditional", datatype= "raw", fpc = TRUE){

  
  #Define variables
  x = N_h = N_h_hats = ps_h = poptotal_h = strat_pop_total = NULL  
  n <- length(y_num)
  
  ### INPUT VALIDATION ###
  
  #Check that y is numeric
  if(!(typeof(y_num) %in% c("numeric", "integer", "double"))){
    stop("Must supply numeric y for the numerator.  For binary variable, convert to 0/1's.")
  }

  if(!(typeof(y_den) %in% c("numeric", "integer", "double"))){
    stop("Must supply numeric y for the denominator.  For binary variable, convert to 0/1's.")
  }  
  
  #Make sure the var_method is valid
  if(!is.element(var_method, c("SRSunconditional"))){
    message("Variance method input incorrect. Currently only \"SRSunconditional\" is supported.")
    return(NULL)
  }

  
  #Need to provide either datatype="raw", N, or pi.  Give warning if not
  if(datatype=="means" & is.null(N) & is.null(pi)){
    message("Must supply N, pi, raw population data or population totals so that we can estimate N.")
    return(NULL)
  }
  
  
  #Check on inclusion probabilities and create weight=inverse inclusion probabilities
  if(is.null(pi)){
    message("Assuming simple random sampling")
  }
  
  #Check is unequal sampling (since can't handle this for var est yet)
  if(length(unique(pi) > 1) & var_est == TRUE){
    message("Currently variance estimators assume simple random sampling.")
  }
  
  # Check that y_num and y_den have same length
  # Proxy for checking that they are from the same sample
  if(length(y_num) != length(y_den)){
    message("The response variables must be from the same sample and so must have the same length.")
    return(NULL)
  }
  
  
  #Need to get N if not provided
  if(is.null(N)){
    if(datatype=="raw"){
      N <- dim(as.matrix(xpop))[1]
    }
    if(datatype=="totals"){
      N <- sum(xpop[,2])#apply(xpop,2,sum)[2]
    }
    if(datatype=="means"){
      N <- sum(pi^(-1))
      message("Assuming N can be approximated by the sum of the inverse inclusion probabilities.")
    }
  }
  
  # convert pi into vector
  if (is.null(pi)) {
    pi <- rep(n/N, n)
  }
  

  
  #Need to get N_h (stratum pop sizes)
  #population design matrix, check whether its population totals, means or raw data
  
  if (datatype=="raw"){
    xpop_tab <- data.frame(table(xpop))
    colnames(xpop_tab) <- c("x","N_h")
  }
  if (datatype=="totals"){
    xpop_tab <- xpop
    colnames(xpop_tab) <- c("x","N_h")    
  }
  if (datatype=="means"){
    xpop_tab <- xpop
    colnames(xpop_tab) <- c("x","N_h")
    xpop_tab$N_h <-  xpop_tab$N_h*N
    
  }
  
  #Make sure that x_pop_tab$x is a factor
  if(!is.factor(xpop_tab$x)){
    xpop_tab$x <- as.factor(xpop_tab$x)
  }
  
  #Make xsample have same column name as xpop_tab
  xsample <- data.frame(xsample)
  colnames(xsample) <- "x"
  
  #Make sure that xsample$x is a factor
  if(!is.factor(xsample$x)){
    xsample$x <- as.factor(xsample$x)
  }
  
  #Compute estimator
  y_num_ps <- postStrat(y = y_num, xsample = xsample,
                        xpop = xpop, var_est =  TRUE, 
                        var_method = var_method,
                        fpc = fpc, datatype = datatype, N = N)
  y_den_ps <- postStrat(y = y_den, xsample = xsample,
                        xpop = xpop, var_est =  TRUE, 
                        var_method = var_method, 
                        fpc = fpc, datatype = datatype, N = N)
  
  #Ratio: Point Estimate
  rat <- y_num_ps$pop_total/y_den_ps$pop_total
  
  #Compute variance estimator
  if(var_est==TRUE){
    
    if(var_method == "SRSunconditional"){
      # Compute components
      xsample_pi_num_den <- data.frame(xsample, pi, y_num, y_den)
      tab <- xsample_pi_num_den %>%
        group_by(x) %>%
        summarize(poptotal_h_num = y_num%*%pi^(-1),
                  poptotal_h_den = y_den%*%pi^(-1),
                  N_h_hats = sum(pi^(-1))) %>%
        inner_join(xpop_tab, by=c("x")) %>%
        mutate(ps_h = N_h/N_h_hats, 
               strat_pop_total_num = ps_h*poptotal_h_num,
               strat_pop_mean_num = strat_pop_total_num/N_h,
               strat_pop_total_den = ps_h*poptotal_h_den,
               strat_pop_mean_den = strat_pop_total_den/N_h)
      
      cov_piece <- xsample_pi_num_den %>%
        left_join(tab, by=c("x")) %>%
        group_by(x) %>%
        summarize(n_h = n(), strat_pop_mean_num = first(strat_pop_mean_num),
                  strat_pop_mean_den = first(strat_pop_mean_den),
                  strat_pop_total_den = first(strat_pop_total_den),
                  strat_pop_total_num = first(strat_pop_total_num),
                  Q = (sum(y_num*y_den) - n_h*strat_pop_mean_num*strat_pop_mean_den)/
                    (n_h*(n_h - 1)))
      
      if(fpc == TRUE){
        cov_est <- (1-n/N)/n*t(tab$N_h/N)%*%(cov_piece$n_h*cov_piece$Q) +
          (1-n/N)/n^2*t(1-tab$N_h/N)%*%(cov_piece$n_h*cov_piece$Q)
      }
      
      if(fpc == FALSE){
        cov_est <- 1/n*t(tab$N_h/N)%*%(cov_piece$n_h*cov_piece$Q) +
          1/n^2*t(1-tab$N_h/N)%*%(cov_piece$n_h*cov_piece$Q)
      }
      
      # Estimate
      var_est <- y_den_ps$pop_mean^(-2)*(y_num_ps$pop_mean_var + 
                                            rat^2*y_den_ps$pop_mean_var - 2*rat*cov_est)
      }
    
    # if(var_method=="bootstrapSRS"){
    #       #Find bootstrap variance
    # 
    #   #Bootstrap total estimates
    #   t_boot <- boot(data = xsample_pi_num_den,
    #                  statistic = postStratRatiot, R = B,
    #                  xpop = xpop)
    #   
    #   #Adjust for bias and without replacement sampling
    #   var_est <- var(t_boot$t)*n/(n-1)*(N-n)/(N-1)
    # 
    #     
    # }
    
    return(list( pop_ratio = rat, 
                 pop_ratio_var = var_est))
  }else{
    return(list( pop_ratio = rat))
    
  }

}
