#varMase
#helper function to compute an approximate variance estimator
#Make flexible enough to work with and without the joint inclusion probabilities

varMase <- function(y, pi, pi2 = NULL, method = "LinHB", N = NULL, fpc = T){
  
  #Need sample size
  n <- length(y)
  
  #Make sure include pi2 if method is LinHT
  if(method=="LinHT" & is.null(pi2)){
    stop("For LinHT variance estimator, need to provide second order inclusion probabilities matrix.")
  }
  
  #Need to add warning messages for if pi2 is not given for LinHT or is wrong dim!
  
  if(method == "LinHB"){
    a <- n/(n-1)*(1-pi)
    e <- as.vector(pi^(-1)*y) - c(sum(a)^(-1)*(pi^(-1)*a)%*%y)
    varEst <- sum(a*e^2)
    
  }
  if(method == "LinHH"){
    t <- pi^(-1)%*%y
    varEst <- 1/(n*(n-1))*t(as.vector(n*y*pi^(-1)) - as.numeric(t))%*%(as.vector(n*y*pi^(-1)) - as.numeric(t))
  }
  if(method == "LinHTSRS"){
    if(is.null(N)){
      N <- sum(pi^(-1))
    }
    
    #varEst <- (N-n)*(N/n)*var(y)
    
    if(fpc == TRUE){
      varEst <- (N-n)*(N/n)*var(y)
    }
    if(fpc == FALSE | is.null(fpc)){
      varEst <- (N^2/n)*var(y)
    }
        
  }
  if(method== "LinHT"){
    a <- (pi2 - pi%*%t(pi))*pi2^(-1)*(pi%*%t(pi))^(-1)*y%*%t(y)
    varEst <- sum(a)
  }
  #LinHT needs joint inclusion probabilities; in the event of lack, use LinHB and LinHH as substitutes; 
  #However, if have simple random, then could just use LinHTSRS or the bootstrap simple random sampling
  #LinHB Hajek-berger estimator, LinHH Hansen Hurwitz
  
  
  return(varEst)
}
