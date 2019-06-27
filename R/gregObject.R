#An R script that contains all of the GREG S3 class methods

#helper function to modify object class in place.
gregify <- function(obj) {
  class(obj) <- "greg"
  return(obj)
}

#print method
#' @export
print.greg <- function(obj) {
  if(!identical(obj$formula, NULL)) {
    f <- as.character(obj$formula)
    cat("formula:", "\n", f[2], f[1], f[3], "\n")
  }
  if(!identical(obj$pop_total_var, NULL)) {
    cat("population estimates:", "\n", "total:", obj$pop_total, "\n",
        "total variance:", obj$pop_total_var, "\n", 
        "mean:", obj$pop_mean, "\n",
        "mean variance:", obj$pop_mean_var, "\n")
  }
  else {
    cat("population estimates:", "\n", "total:", obj$pop_total, "\n",
        "mean:", obj$pop_mean, "\n")
  }
  if(!identical(obj$strat_ests, NULL)) {
    cat("stata estimates:", obj$strat_ests, "\n")
  }
  if(!identical(obj$oob_error, NULL)) {
    cat("out of bag error:", obj$oob_error, "\n")
  }
} 

#predict method
#' @export
predict.greg <- function(obj, new_data) {
  if(class(obj) == "greg"){
    if(!identical(obj$logistic_model, NULL)){
      return(predict(obj$logistic_model, new_data))
      }    
    if(!identical(obj$tree, NULL)){
      return(predict(obj$tree, new_data))
      }
    if(!identical(obj$forest, NULL)){
      return(predict(obj$forest, new_data))
      }
    if(!identical(obj$coefficients, NULL)){
      if(!identical(names(obj$coefficients), NULL)){
        var_names <- (obj$coefficients) %>% names()
      }
      else{
        var_names <- (obj$coefficients)[,1] %>% names()
      }
      if("(Intercept)" %in% var_names){
        var_names <- var_names[-1]
      }
      if(length(var_names) != base::ncol(new_data)){
        new_data <- new_data %>% dplyr::select(var_names)
      }
      design <- model.matrix(~., data = new_data)
      return((design %*% (as.vector(obj$coefficients))) %>%
               as.vector())
      }
    else{
      message("This greg object has no predict method.")
      return(NULL)
    }
  }
  else{
    message("object class is not 'greg'.")
    return(NULL)
  }
}

#plot method
#' Compute a geospatial heat map to summarize a generalized regression estimator
#' 
#' Computes a geospatial heat map of fitted values of a 'greg' object
#'
#' @param obj Object of class 'greg'.
#' @param x_sample A data frame of the sample level data with two spatial variables. Must have the same indicies as the original data frame used to train model.
#' @param x_pop A data frame of the population level data with two spatial variables. Must have the same indicies as the original population data frame used to make model.
#' @param spatial_1 A string for the name of variable in `x_sample` and `x_pop` to plot as x axis (column name of longitude). 
#' @param spatial_2 A string for the name of variable in `x_sample` and `x_pop` to plot as y axis (column name of latitude). 
#' @param bins An integer specifying the number of bins to aggregate estimated values
#' @param hexagon Boolean. Specify if you want hexagonal shaped bins or not.
#' @param fun Summary function to apply within each bin to be plotted. Default is `mean(x)`. See `stat_summary`.
#' 
#' @export plot.greg
#' @import dplyr
#' @import ggplot2
#' 
plot.greg <- function(obj, x_sample, x_pop, spatial_1, spatial_2, bins = 300,
                      fun = function(x) mean(x), hexagon = FALSE) {
  #Check if population estimates are available
  if(identical(obj$y_hat_pop, NULL)){
    message("individual population estimates not available without 'raw' data type")
    obj$y_hat_pop <- rep(NA, nrow(x_pop))
  }
  #create dataframe
  d <- rbind(x_sample, x_pop) %>%
    dplyr::select(x = spatial_1,y = spatial_2) %>%
    mutate(y_hat = c(obj$y_hat_sample, obj$y_hat_pop),
           plot_key = c(rep("sample", nrow(x_sample)),
                        rep("pop", nrow(x_pop)))) %>%
    arrange(plot_key)
  #summary kernel
  if(hexagon == TRUE){
    kern <- stat_summary_hex(data = d %>% filter(plot_key == "pop"),
                             bins = bins,
                             fun = fun)
  }
  else{
    kern <- stat_summary_2d(data = d %>% filter(plot_key == "pop"),
                            bins = bins,
                            fun = fun)
  }
  #create ggplot
  plt <- ggplot(d, aes(x = x, y = y, z = y_hat)) +
    kern +
    geom_point(data = d %>% filter(plot_key == "sample")) +
    labs(fill = "y_hat") + theme_minimal()
  #return plot using viridis scale, if client has it
  try(return(plt + scale_fill_viridis_c()), silent = TRUE)
  return(plt)
}
