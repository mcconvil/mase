#varMase
#helper function to compute an approximate variance estimator
#Make flexible enough to work with and without the joint inclusion probabilities

varMase <- function(y, pi, pi2 = NULL, method="HB", N = NULL){
  
  #Need sample size
  n <- length(y)
  
  

  #Make sure include pi2 if method is HT
  if(method=="HT" & is.null(pi2)){
    message("For HT variance estimator, need to provide second order inclusion probabilities matrix.")
    return(NULL)
  }
  
  #Need to add warning messages for if pi2 is not given for HT or is wrong dim!
  
  if(method == "HB"){
    a <- n/(n-1)*(1-pi)
    e <- as.vector(pi^(-1)*y) - c(sum(a)^(-1)*(pi^(-1)*a)%*%y)
    varEst <- sum(a*e^2)
    
  }
  if(method == "HH"){
    t <- pi^(-1)%*%y
    varEst <- 1/(n*(n-1))*t(as.vector(n*y*pi^(-1)) - as.numeric(t))%*%(as.vector(n*y*pi^(-1)) - as.numeric(t))
  }
  if(method == "HTSRS"){
    if(is.null(N)){
      N <- sum(pi^(-1))
    }
    varEst <- (N-n)*(N/n)*var(y)
  }
  if(method== "HT"){
    a <- (pi2 - pi%*%t(pi))*pi2^(-1)*(pi%*%t(pi))^(-1)*y%*%t(y)
    varEst <- sum(a)
  }
  
  
  return(varEst)
}
