#' Compute a post-stratified estimator
#' 
#' Calculates a post-stratified estimator for a finite population mean/proportion or total based on sample data collected from a complex sampling design and a single, categorical auxiliary population variable.  
#' 
#' @inheritParams horvitzThompson
#' @inheritParams greg
#' @param xsample A vector containing the post-stratum for each sampled unit.
#' @param xpop A vector or data frame, depending on datatype.  If datatype = "raw", then a vector containing the post-stratum for each population unit.  If datatype = "totals" or "means", then a data frame, where the first column lists the possible post-strata and the second column contains the population total or proportion in each post-stratum.  
#' @param datatype Default to "raw", takes values "raw", "totals" or "means" for whether the user is providing the raw population stratum memberships, the population totals of each stratum, or the population proportions of each stratum.
#' @param var_est Default to FALSE, logical for whether or not to compute estimate of variance
#' @param var_method The method to use when computing the variance estimator.  Options are a Taylor linearized technique: "LinHB"= Hajek-Berger estimator, "LinHH" = Hansen-Hurwitz estimator, "LinHTSRS" = Horvitz-Thompson estimator under simple random sampling without replacement, and "LinHT" = Horvitz-Thompson estimator or a resampling technique: "bootstrapSRS" = bootstrap variance estimator under simple random sampling without replacement, "SRSunconditional" = simple random sampling variance estimator which accounts for random strata.
#' @param fpc Default to TRUE, logical for whether or not the variance calculation should include a finite population correction when calculating the "LinHTSRS", "SRSunconditional", or the "SRSbootstrap" variance estimator.

#' @examples 
#' library(tidyr)
#' library(dplyr)
#' 
#' data(IdahoPop)
#' data(IdahoSamp)
#' 
#' xsample <- filter(IdahoSamp, COUNTYFIPS == "16055")
#' xpop <- filter(IdahoPop, COUNTYFIPS == "16055")
#' 
#' pop <- xpop[c("tnt.1", "tnt.2")] |> 
#'     pivot_longer(everything(), names_to = "tnt", values_to = "prop") |>
#'     mutate(tnt = as.numeric(gsub("\\D", "", tnt)))
#'     
#' postStrat(y = xsample$BA_TPA_ADJ,
#'           N = xpop$npixels,
#'           xsample = xsample$tnt,
#'           xpop = pop,
#'           datatype = "means",
#'           var_est = TRUE,
#'           var_method = "SRSunconditional")
#'           
#' @references 
#' \insertRef{coc77}{mase} 
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
#' * pop_total_var: Estimated variance of population total estimate.
#' 
#' * pop_mean_var: Estimated variance of population mean estimate.
#' 
#' * strat_ests: Table of total and mean estimates for each strata.
#' 
#' * weights: Survey weights produced by PS.
#' 
#' @import dplyr
#' @import tidyr
#' @export postStrat
#' @include bootHelpers.R
#' @include varMase.R


postStrat <- function(y,
                      xsample,
                      xpop,
                      pi = NULL,
                      N = NULL, 
                      var_est = FALSE, 
                      var_method = "LinHB", 
                      pi2 = NULL, 
                      datatype = "raw",
                      B = 1000,
                      fpc = TRUE,
                      messages = TRUE){

  
  #Define variables
  x <- N_h <- N_h_hats <- ps_h <- poptotal_h <- strat_pop_total <- NULL  

  if(!(typeof(y) %in% c("numeric", "integer", "double"))) {
    stop("Must supply numeric y.  For binary variable, convert to 0/1's.")
  }
  
  if (!is.element(var_method, c("LinHB", "LinHH", "LinHTSRS", "LinHT", "bootstrapSRS", "SRSunconditional"))) {
    stop("Variance method input incorrect. It has to be \"LinHB\", \"LinHH\", \"LinHT\", \"LinHTSRS\", \"SRSunconditional\", or \"bootstrapSRS\".")
  }

  if (datatype == "means" & is.null(N) & is.null(pi)) {
    stop("Must supply N, pi, raw population data or population totals so that we can estimate N.")
  }
  
  if (is.null(pi)) {
    if (messages) {
      message("Assuming simple random sampling") 
    }
  }
  
  #Need to get N if not provided
  if (is.null(N)) {
    if (datatype =="raw") {
      N <- dim(as.matrix(xpop))[1]
    }
    if (datatype == "totals") {
      N <- sum(xpop[,2])#apply(xpop,2,sum)[2]
    }
    if (datatype == "means") {
      N <- sum(pi^(-1))
      if (messages) {
        message("Assuming N can be approximated by the sum of the inverse inclusion probabilities.") 
      }
    }
  }
  
  if (is.null(pi)) {
    pi <- rep(length(y)/N, length(y))
  }
  
  #Need to get N_h (stratum pop sizes)
  #population design matrix, check whether its population totals, means or raw data
  
  if (datatype=="raw") {
    xpop_tab <- data.frame(table(xpop))
    colnames(xpop_tab) <- c("x","N_h")
  }
  if (datatype=="totals") {
    xpop_tab <- xpop
    colnames(xpop_tab) <- c("x","N_h")    
  }
  if (datatype=="means") {
    xpop_tab <- xpop
    colnames(xpop_tab) <- c("x","N_h")
    xpop_tab$N_h <-  xpop_tab$N_h*N
    
  }
  
  #Make sure that x_pop_tab$x is a factor
  if (!is.factor(xpop_tab$x)) {
    xpop_tab$x <- as.factor(xpop_tab$x)
  }
  
  #Make xsample have same column name as xpop_tab
  xsample <- data.frame(xsample)
  colnames(xsample) <- "x"
  
  #Make sure that xsample$x is a factor
  if (!is.factor(xsample$x)) { 
    xsample$x <- as.factor(xsample$x)
  }
  
  #Compute estimator
  xsample_pi_y <- data.frame(xsample, pi,y)
  tab <- xsample_pi_y |>
    group_by(x) |>
    summarize(poptotal_h = sum(y * pi^(-1)), N_h_hats = sum(pi^(-1)), var_h = var(y)) |>
    inner_join(xpop_tab, by=c("x")) |>
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
  if (var_est == TRUE) {
    n <- length(y)
    if (var_method == "SRSunconditional") {
      
      if (fpc == TRUE) {
        var_est <- N^2*(1-n/N)/n*t(tab$N_h/N)%*%tab$var_h + N^2*(1-n/N)/n^2*t(1-tab$N_h/N)%*%tab$var_h
      }
      
      if (fpc == FALSE) {
        var_est <- N^2/n*t(tab$N_h/N)%*%tab$var_h + N^2/n^2*t(1-tab$N_h/N)%*%tab$var_h
      }
            
    }
    
    if (var_method == "bootstrapSRS") {

      #Bootstrap total estimates
      t_boot <- boot(data = xsample_pi_y,
                     statistic = postStratt,
                     R = B, 
                     xpoptab = xpop_tab)
      
      if (fpc == TRUE) {
        var_est <- var(t_boot$t)*n/(n-1)*(N-n)/(N-1)
      }
      if (fpc == FALSE) {
        var_est <- var(t_boot$t)*n/(n-1)
      }
        
    }
    
    if (is.element(var_method, c("LinHB", "LinHH", "LinHT", "LinHTSRS"))) {
      
      y_hat <- dat_s$strat_pop_mean
      var_est <- varMase(y = (y - y_hat),
                         pi = pi,
                         pi2 = pi2,
                         method = var_method,
                         N = N,
                         fpc = fpc)
      
    }
    return(list(pop_total = pop_total, 
                pop_mean = pop_total/N,
                pop_total_var = var_est, 
                pop_mean_var = var_est/N^2, 
                strat_ests = strat_ests,
                weights = wts))
  } else { 
    return(list(pop_total = pop_total, 
                pop_mean = pop_total/N,
                strat_ests = strat_ests,
                weights = wts))
    
  }

}
