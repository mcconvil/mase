#' Compute a post-stratified estimator
#' 
#' Calculates a generalized regression estimator for a finite population mean/proportion or total based on sample data collected from a complex sampling design and auxiliary population data.  
#' 
#' @inheritParams horvitzThompson
#' @inheritParams greg
#' @param x_sample A vector containing the post-stratum for each sampled unit.
#' @param x_pop A vector or data frame, depending on data_type.  If data_type = "raw", then a vector containing the post-stratum for each population unit.  If data_type = "totals" or "means", then a data frame, where the first column lists the possible post-strata and the second column contains the population total or proportion in each post-stratum.  
#' @param data_type Default to "raw", takes values "raw", "totals" or "means" for whether the user is providing the raw population stratum memberships, the population totals of each stratum, or the population proportions of each stratum.
#' @param var_est Default to FALSE, logical for whether or not to compute estimate of variance
#' @param var_method The method to use when computing the variance estimator.  Options are a Taylor linearized technique: "lin_HB"= Hajek-Berger estimator, "lin_HH" = Hansen-Hurwitz estimator, "lin_HTSRS" = Horvitz-Thompson estimator under simple random sampling without replacement, and "lin_HT" = Horvitz-Thompson estimator or a resampling technique: "bootstrap_SRS" = bootstrap variance estimator under simple random sampling without replacement, "unconditional_SRS" = simple random sampling variance estimator which accounts for random strata.
#' 
#' @examples 
#' library(survey)
#' data(api) 
#' postStrat(y = apisrs$api00, x_sample = apisrs$awards, 
#' x_pop = data.frame(table(apipop$awards)), data_type = "totals", 
#' pi = apisrs$pw^(-1))
#' 
#'@references 
#'\insertRef{coc77}{mase} 
#' 
#'\insertRef{sar92}{mase}
#'
#' @return A list of output containing:
#' \itemize{
#' \item{pop_total: Estimate of population total}
#' \item{pop_mean: Estimate of the population mean}
#' \item{pop_total_var: Estimated variance of population total estimate}
#' \item{pop_mean_var: Estimated variance of population mean estimate}
#' \item{strat_ests: Table of total and mean estimates for each strata}
#' \item{weights: Survey weights produced by post stratification}
#' }
#' @import dplyr
#' @import magrittr
#' @export postStrat
#' @include postStratt.R
#' @include varMase.R
#' @include gregObject.R
#' 
#' @seealso \code{\link{greg}} for a linear or logistic regression model.


postStrat <- function(
  y, x_sample, x_pop, pi = NULL, N=NULL, var_est = FALSE, var_method="lin_HB", pi2 = NULL, data_type= "raw", B=1000, strata = NULL){

  
  #Define variables
  x= N_h = N_h_hats = ps_h = poptotal_h = strat_pop_total = NULL  
  
  
  ### INPUT VALIDATION ###
  
  #Check that y is numeric
  if(!(typeof(y) %in% c("numeric", "integer", "double"))){
    stop("Must supply numeric y.  For binary variable, convert to 0/1's.")
  }
  
  
  #Make sure the var_method is valid
  if(!is.element(var_method, c("lin_HB", "lin_HH", "lin_HT", "lin_HTSRS", "bootstrap_SRS", "unconditional_SRS"))){
    message("Variance method input incorrect. It has to be \"lin_HB\", \"lin_HH\", \"lin_HT\", \"lin_HTSRS\", \"unconditional_SRS\" or \"bootstrap_SRS\".")
    return(NULL)
  }
  
  
  #Need to provide either data_type="raw", N, or pi.  Give warning if not
  if(data_type=="means" & is.null(N) & is.null(pi)){
    message("Must supply N, pi, raw population data or population totals so that we can estimate N.")
    return(NULL)
  }
  
  
  #Check on inclusion probabilities and create weight=inverse inclusion probabilities
  if(is.null(pi)){
    message("Assuming simple random sampling")
  }
  
  #Need to get N if not provided
  if(is.null(N)){
    if(data_type=="raw"){
      N <- dim(as.matrix(x_pop))[1]
    }
    if(data_type=="totals"){
      N <- sum(x_pop[,2])#apply(x_pop,2,sum)[2]
    }
    if(data_type=="means"){
      N <- sum(pi^(-1))
      message("Assuming N can be approximated by the sum of the inverse inclusion probabilities.")
    }
  }
  
  # convert pi into vector
  if (is.null(pi)) {
    pi <- rep(length(y)/N, length(y))
  }
  

  
  #Need to get N_h (stratum pop sizes)
  #population design matrix, check whether its population totals, means or raw data
  
  if (data_type=="raw"){
    x_pop_tab <- data.frame(table(x_pop))
    colnames(x_pop_tab) <- c("x","N_h")
  }
  if (data_type=="totals"){
    x_pop_tab <- x_pop
    colnames(x_pop_tab) <- c("x","N_h")    
  }
  if (data_type=="means"){
    x_pop_tab <- x_pop
    colnames(x_pop_tab) <- c("x","N_h")
    x_pop_tab$N_h <-  x_pop_tab$N_h*N
    
  }
  
  #Make sure that x_pop_tab$x is a factor
  if(!is.factor(x_pop_tab$x)){
    x_pop_tab$x <- as.factor(x_pop_tab$x)
  }
  
  #Make x_sample have same column name as x_pop_tab
  x_sample <- data.frame(x_sample)
  colnames(x_sample) <- "x"
  
  #Make sure that x_sample$x is a factor
  if(!is.factor(x_sample$x)){
    x_sample$x <- as.factor(x_sample$x)
  }
  
  #Compute estimator
  x_sample_pi_y <- data.frame(x_sample, pi,y)
  tab <- x_sample_pi_y %>%
    group_by(x) %>%
    summarize(poptotal_h = y%*%pi^(-1),N_h_hats = sum(pi^(-1)), var_h = var(y)) %>%
    inner_join(x_pop_tab, by=c("x")) %>%
    mutate(ps_h = N_h/N_h_hats, strat_pop_total = ps_h*poptotal_h, strat_pop_mean = strat_pop_total/N_h)
  
  #Estimates by strata  
  strat_ests <- tab[,c("x", "strat_pop_total", "strat_pop_mean")]
  #Estimate of population total
  pop_total <- sum(tab$strat_pop_total)
  
  #Create dataset with strata info
  dat_s <- left_join(x_sample_pi_y, tab, by="x") 
  
  #Survey weights  
  wts <- 1/dat_s$pi*dat_s$ps_h
  
  #Compute variance estimator
  if(var_est==TRUE){
    n <- length(y)
    if(var_method=="unconditional_SRS"){
      varEst <- N^2*(1-n/N)/n*t(tab$N_h/N)%*%tab$var_h + N^2*(1-n/N)/n^2*t(1-tab$N_h/N)%*%tab$var_h
    }
    
    if(var_method=="bootstrap_SRS"){
          #Find bootstrap variance

      #Bootstrap total estimates
      t_boot <- boot(data = x_sample_pi_y, statistic = postStratt, R = B, x_poptab = x_pop_tab)
      
      #Adjust for bias and without replacement sampling
      varEst <- var(t_boot$t)*n/(n-1)*(N-n)/(N-1)

        
    }
    
    if(is.element(var_method, c("lin_HB", "lin_HH", "lin_HT", "lin_HTSRS"))){
      
      y_hat <- dat_s$strat_pop_mean
      varEst <- varMase(y = (y-y_hat),pi = pi,pi2 = pi2, method = var_method, N = N, strata = strata)
      
    }
    return(list( pop_total = pop_total, 
                 pop_mean = pop_total/N,
                 pop_total_var = varEst, 
                 pop_mean_var = varEst/N^2, 
                 strat_ests = strat_ests,
                 weights = wts) %>% gregify())
  }else{
    return(list( pop_total = pop_total, 
                 pop_mean = pop_total/N,
                 strat_ests = strat_ests,
                 weights = wts) %>% gregify())
    
  }

}
