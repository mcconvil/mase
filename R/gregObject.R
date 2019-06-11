#An R script that contains all of the GREG S3 class methods

#Helper function to modify object class in place
gregify <- function(obj) {
  class(obj) <- "greg"
  return(obj)
}

#Print method for greg object
print.greg <- function(obj) {
  if(!identical(obj$formula, NULL)){
    f <- as.character(obj$formula)
    cat("formula:", "\n", f[2], f[1], f[3], "\n")
  }
  if(!identical(obj$pop_total_var, NULL)){
    cat("population estimates:", "\n", "total:", obj$pop_total, "\n",
        "total variance:", obj$pop_total_var, "\n", 
        "mean:", obj$pop_mean, "\n",
        "mean variance:", obj$pop_mean_var, "\n")
  }
  else{
    cat("population estimates:", "\n", "total:", obj$pop_total, "\n",
        "mean:", obj$pop_mean, "\n")
  }
} 

summary.greg <- function(obj) {
  cat("variable importance:")
}

predict.greg <- function(obj, new_data) {
  if(!identical(obj$coefficients, NULL)){
    X <- as.matrix(new_data)
    return(X %*% obj$coefficients)
  }
  if(!identical(obj$tree, NULL)){
    return(predict(obj$tree, new_data))
  }
  if(!identical(obj$forest, NULL)){
    return(predict(obj$forest, new_data))
  }
  else{
    return(NULL)
  }
}

plot.greg <- function(obj, spatial1, spatial2) {
  ggplot()
}