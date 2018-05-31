#Helper function to compute gregtreet total for bootstrapping
library(rpms)

gregTreet <- function(data, x_pop, p_value= p_value, perm_reps = perm_reps, bin_size = bin_size, indices){
  #data: 1st column:y, 2nd column:pis, rest: x_sample
  d <- data[indices,]
  
  #y
  y <- d[,1]
  
  #pis 
  pis <- d[,2]
  
  #Length of x_sample_d
  p <- dim(d)[2] - 2
  #x_sample
  x_sample <- d[, 3:(p + 2)]
  
  #create tree
  f <- as.formula(paste("y ~ ", paste(names(x_sample), collapse= "+")))
  treet <- rpms(rp_equ = f, data = d, weights = as.vector(pis^(-1)), pval= p_value, perm_reps = perm_reps, bin_size = bin_size)
  
  
  #Calculate weights
  #Make sure x_pop and x_sample have the same columns (in same order)
  x_pop <- x_pop[names(x_sample)]
  #Design matrix for population
  x_pop_treet <- treeDesignMatrix(splits = treet$ln_split, data = x_pop)
  #Design matrix for sample
  x_sample_treet <- treeDesignMatrix(splits = treet$ln_split, data = x_sample)
  w <- (1 + t(colSums(x_pop_treet)- colSums(x_sample_treet*pis^(-1)))%*%solve(t(x_sample_treet)%*%diag(pis^(-1))%*%as.matrix(x_sample_treet))%*%t(x_sample_treet))%*%diag(pis^(-1)) 
  
  #calculating the total estimate for y
  t <- w %*% y

  return(t)
}

