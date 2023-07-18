#Helper function to compute gregtreet total for bootstrapping
library(rpms)

gregTreet <- function(data, xpop, pval= pval, perm_reps = perm_reps, bin_size = bin_size, indices){
  #data: 1st column:y, 2nd column:pis, rest: xsample
  d <- data[indices,]
  
  #y
  y <- d[,1]
  
  #pis 
  pis <- d[,2]
  
  #Length of xsample_d
  p <- dim(d)[2] - 2
  #xsample
  xsample <- d[, 3:(p + 2)]
  
  #create tree
  f <- as.formula(paste("y ~ ", paste(names(xsample), collapse= "+")))
  treet <- rpms(rp_equ = f, data = d, weights = as.vector(pis^(-1)), pval = pval, perm_reps = perm_reps, bin_size = bin_size)
  
  
  #Calculate weights
  #Make sure xpop and xsample have the same columns (in same order)
  xpop <- xpop[names(xsample)]
  #Design matrix for population
  xpop_treet <- box_ind(treet, xpop)
  #Design matrix for sample
  xsample_treet <- box_ind(treet, xsample)
  w <- (1 + t(colSums(xpop_treet)- colSums(xsample_treet*pis^(-1))) %*%
          solve(t(xsample_treet) %*% diag(pis^(-1)) %*% as.matrix(xsample_treet)) %*%
          t(xsample_treet)) %*% diag(pis^(-1)) 
  
  #calculating the total estimate for y
  t <- w %*% y

  return(t)
}

