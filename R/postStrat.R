#' Compute a post-stratified estimator
#' 
#' Calculates a generalized regression estimator for a finite population mean/proportion or total based on sample data collected from a complex sampling design and auxiliary population data.  
#' 
#' @inheritParams horvitzThompson
#' @inheritParams greg
#' @param xsample A vector containing the post-stratum for each sampled unit.
#' @param xpop A vector or data frame, depending on datatype.  If datatype = "raw", then a vector containing the post-stratum for each population unit.  If datatype = "totals" or "means", then a data frame, where the first column lists the possible post-strata and the second column contains the population total or proportion in each post-stratum.  
#' @param datatype Default to "raw", takes values "raw", "totals" or "means" for whether the user is providing the raw population stratum memberships, the population totals of each stratum, or the population proportions of each stratum.
#' @param var_est Default to FALSE, logical for whether or not to compute estimate of variance
#' @param var_method The method to use when computing the variance estimator.  Options are a Taylor linearized technique: "LinHB"= Hajek-Berger estimator, "LinHH" = Hansen-Hurwitz estimator, "LinHTSRS" = Horvitz-Thompson estimator under simple random sampling without replacement, and "LinHT" = Horvitz-Thompson estimator or a resampling technique: "bootstrapSRS" = bootstrap variance estimator under simple random sampling without replacement, "SRSunconditional" = simple random sampling variance estimator which accounts for random strata.
#' 
#' @examples 
#' 
#' postStrat(y = apisrs$api00, xsample = apisrs$awards, xpop = data.frame(table(apipop$awards)), datatype = "totals", pi = apisrs$pw^(-1))
#' 
#'@references 
#'\insertRef{coc77}{mase} 
#' 
#'\insertRef{sar92}{mase}
#'
#' @return A list of output containing:
#' \itemize{
#' \item{pop_total:}{Estimate of population total}
#' \item{pop_mean:}{Estimate of the population mean (or proportion)}
#' \item{pop_total_var:}{Estimated variance of population total estimate}
#' \item{pop_mean_var:}{Estimated variance of population mean estimate}
#' \item{strat_ests:}{Table of total and mean estimates for each strata}
#' \item{weights.ps:}{Survey weights produced by PS}
#' }

#' @import dplyr
#' @import magrittr
#' @export postStrat
#' @include postStratt.R
#' @include varMase.R


postStrat <- function(
  y, xsample, xpop, pi = NULL, N=NULL, var_est =FALSE, var_method="HB", pi2 = NULL, datatype= "raw", B=1000){

  
  #Define variables
  x= N_h = N_h_hats = ps_h = poptotal_h = strat_pop_total = NULL  
  
  
  ### INPUT VALIDATION ###
  
  #Check that y is numeric
  if(!(typeof(y) %in% c("numeric", "integer", "double"))){
    stop("Must supply numeric y.  For binary variable, convert to 0/1's.")
  }
  
  
  #Make sure the var_method is valid
  if(!is.element(var_method, c("HB", "HH", "HTSRS", "HT", "bootstrapSRS", "SRSunconditional"))){
    message("Variance method input incorrect. It has to be \"HB\", \"HH\", \"HT\", \"HTSRS\", \"SRSunconditional\" or \"bootstrapSRS\".")
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
    pi <- rep(length(y)/N, length(y))
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
  xsample_pi_y <- data.frame(xsample, pi,y)
  tab <- xsample_pi_y %>%
    group_by(x) %>%
    summarize(poptotal_h = y%*%pi^(-1),N_h_hats = sum(pi^(-1)), var_h = var(y)) %>%
    inner_join(xpop_tab, by=c("x")) %>%
    mutate(ps_h = N_h/N_h_hats, strat_pop_total = ps_h*poptotal_h, strat_pop_mean = strat_pop_total/N_h)
  
  #Estimates by strata  
  strat_ests <- tab[,c("x", "strat_pop_total", "strat_pop_mean")]
  #Estimate of population total
  pop_total <- sum(tab$strat_pop_total)
  
  #Create dataset with strata info
  dat_s <- left_join(xsample_pi_y, tab, by="x") 
  
  #Survey weights  
  wts <- dat_s$pi*dat_s$ps_h
  
  #Compute variance estimator
  if(var_est==TRUE){
    n <- length(y)
    if(var_method=="SRSunconditional"){
      var_est <- N^2*(1-n/N)/n*t(tab$N_h/N)%*%tab$var_h + N^2*(1-n/N)/n^2*t(1-tab$N_h/N)%*%tab$var_h
    }
    
    if(var_method=="bootstrapSRS"){
          #Find bootstrap variance

      #Bootstrap total estimates
      t_boot <- boot(data = xsample_pi_y, statistic = postStratt, R = B, xpoptab = xpop_tab)
      
      #Adjust for bias and without replacement sampling
      var_est <- var(t_boot$t)*n/(n-1)*(N-n)/(N-1)

        
    }
    
    if(is.element(var_method, c("HB", "HH", "HT", "HTSRS"))){
      
      y_hat <- dat_s$strat_pop_mean
      var_est <- varMase(y = (y-y_hat),pi = pi,pi2 = pi2, method = var_method, N = N)
      
    }
    return(list( pop_total = pop_total, 
                 pop_mean = pop_total/N,
                 pop_total_var = var_est, 
                 pop_mean_var = var_est/N^2, 
                 strat_ests = strat_ests,
                 weights = wts))
  }else{
    return(list( pop_total = pop_total, 
                 pop_mean = pop_total/N,
                 strat_ests = strat_ests,
                 weights = wts))
    
  }

}
