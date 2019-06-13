#An R script that contains all of the GREG S3 class methods

#helper function to modify object class in place.
gregify <- function(obj) {
  class(obj) <- "greg"
  return(obj)
}

#print method
#' @export
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

#plot method (work in progress)
plot.greg <- function(obj, x_sample, x_pop, spatial_x, spatial_y, bins = 300,
                      fun = function(x) mean(x)) {
  #create dataframe
  d <- x_sample %>%
    mutate(plot_key = "sample") %>%
    rbind(x_pop %>%
            mutate(plot_key = "pop"))
  #predict using GREG object
  pred <- predict.greg(obj = obj,
                       new_data = d %>% dplyr::select(-LAT_PUBLIC))
  d <- d %>% mutate(y_hat = pred)
  
  #create ggplot
  plt <- ggplot(d %>% arrange(plot_key),
                    aes(x = LON_PUBLIC, y = LAT_PUBLIC, z = nlcd11)) +
    stat_summary_hex(data = d %>% filter(key == "pop"),
                     position = "identity",
                     bins = bins,
                     fun = fun) +
    scale_fill_viridis_c() +
    geom_point(data = d %>% filter(key == "plot")) +
    labs(title = "Daggett County Survey Plot Locations", x = "longitude", y = "latitude",
         fill = "")
  
  return(plt)
}
