#varMase
#helper function to compute an approximate variance estimator
#Make flexible enough to work with and without the joint inclusion probabilities

varMase <- function(y, pi, pi2 = NULL, method = "lin_HB", N = NULL, strata = NULL){
  
  #Need sample size
  n <- length(y)
  
  
  #If strata = NULL, place all in same stratum
  if(is.null(strata)){
    strata <- rep("a", n)
  }

  #Make sure include pi2 if method is lin_HT
  if(method=="lin_HT" & is.null(pi2)){
    message("For lin_HT variance estimator, need to provide second order inclusion probabilities matrix.")
    return(NULL)
  }
  
  #Ensure pi2 is the correct size
  if(method == "lin_HT" & !is.null(pi2)){
    if(dim(pi2)[1] != n | dim(pi2)[2] != n){
      message("pi2 must be an n by n matrix of secord order inclusion probabilities where n is the length of y.")
      return(NULL)
    }
  }
    
  
  dat <- data.frame(y = y, pi = pi, strata = strata)

  if(method == "lin_HB"){
    
    varEst <- dplyr::group_by(dat, strata) %>%
      dplyr::mutate(n_h = n(), a = n_h/(n_h-1)*(1-pi)) %>%
      dplyr::summarize(varEst_h = sum(a*(as.vector(pi^(-1)*y) - c(sum(a)^(-1)*(pi^(-1)*a)%*%y))^2)) %>%
      dplyr::ungroup() %>%
      dplyr::summarize(varEst = sum(varEst_h)) %>%
      magrittr::extract2(1)
    
  }
  
  if(method == "lin_HH"){
    
    varEst <- dplyr::group_by(dat, strata) %>%
      dplyr::summarize(n_h = n(), t_h = sum(pi^(-1)*y), 
                       varEst_h = 1/(n_h*(n_h-1))*t(as.vector(n_h*y*pi^(-1)) - as.numeric(t_h))%*%(as.vector(n_h*y*pi^(-1)) - as.numeric(t_h))) %>%
      dplyr::ungroup() %>%
      dplyr::summarize(varEst = sum(varEst_h)) %>%
      magrittr::extract2(1)
      
  }
  
  if(method == "lin_HTSRS"){
    if(is.null(N)){
      N <- sum(pi^(-1))
    }
    varEst <- (N-n)*(N/n)*var(y)
  }
  
  if(method== "lin_HT"){
    a <- (pi2 - pi%*%t(pi))*pi2^(-1)*(pi%*%t(pi))^(-1)*y%*%t(y)
    varEst <- sum(a)
  }
  #lin_HT needs joint inclusion probabilities; in the event of lack, use lin_HB and lin_HH as substitutes; 
  #However, if have simple random, then could just use lin_HTSRS or the bootstrap simple random sampling
  #lin_HB Hajek-berger estimator, lin_HH Hansen Hurwitz
  
  
  return(varEst)
}
