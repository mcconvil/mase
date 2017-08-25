#Helper function to compute gregtreet total for bootstrapping
library(rpms)

gregTreet <- function(data, xpop, indices){
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
  treet <- rpms(rp_equ = f, data=d, weights=pis, pval= pval, perm_reps = perm_reps, bin_size = bin_size)
  
  
  #Calculate weights
  #Make sure xpop and xsample have the same columns (in same order)
  xpop <- xpop[names(xsample)]
  #Design matrix for population
  xpop_treet <- treeDesignMatrix(splits = treet$ln_split, data = xpop)
  #Design matrix for sample
  xsample_treet <- treeDesignMatrix(splits = treet$ln_split, data = xsample)
  w <- (1 + t(colSums(xpop_treet)- colSums(xsample_treet*pis^(-1)))%*%solve(t(xsample_treet)%*%diag(pis^(-1))%*%as.matrix(xsample_treet))%*%t(xsample_treet))%*%diag(pis^(-1)) 
  
  #calculating the total estimate for y
  t <- w %*% y

  return(t)
}

