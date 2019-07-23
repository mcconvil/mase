#' Compute a neural network model assisted estimator
#' 
#' Calculates a dense feed forward neural network model assisted estimator for a finite population mean or total based on sample data collected from a complex sampling design and auxiliary population data. Requires installation of keras.  
#' 
#' @inheritParams horvitzThompson
#' @param x_sample A data frame of the auxiliary data in the sample.
#' @param x_pop A data frame of population level auxiliary information. It must contain the same names as x_sample.
#' @param H An integer for the size of the hidden layers.
#' @param lr The learning rate of the gradient descent algorithm.
#' @param epochs The number of training epochs to run.
#' @param batch_size Integer or NULL. Number of samples per gradient update. If unspecified, batch_size will default to 32.
#' @param validation_split Proportion of sample to make validation set for each epoch. Default is 0.2
#' @param optimizer Choice of gradient descent algorithm, default is `keras::optimizer_rmsprop()`.
#' @param loss Choice of loss function, default is MSE. see `keras::compile.keras`
#' @param verbose Verbosity mode (0 = silent, 1 = progress bar, 2 = one line per epoch).
#' 
#' @examples
#' library(survey)
#' data(api)
#' gregNeuralNet(y = apisrs$api00, 
#' x_sample = apisrs[c("col.grad", "awards", "snum", "dnum", "cnum", "pcttest", "meals", "sch.wide")], 
#' x_pop = apipop[c("col.grad", "awards", "snum", "dnum", "cnum", "pcttest", "meals", "sch.wide")])
#' 
#'@references 
#'\insertRef{mcc17b}{mase}

#' 
#' @return A greg object containing:
#' \itemize{
#' \item{pop_total: Estimate of population total}
#' \item{pop_mean: Estimate of the population mean}
#' \item{pop_total_var: Estimated variance of population total estimate}
#' \item{pop_mean_var: Estimated variance of population mean estimate}
#' \item{formula: Model formula}
#' \item{model: feed forward neural network model}
#' \item{y_hat_sample: Response estimates for sample}
#' \item{y_hat_pop: Response estimates for population}
#' }
#' 
#' @export gregNeuralNet
#' @import dplyr
#' @importFrom stats as.formula
#' @include varMase.R
#' @include gregt.R
#' @include gregObject.R
#' 
#' @seealso \code{\link{greg}} for a linear or logistic regression model.


gregNeuralNet <- function(y, x_sample, x_pop, pi = NULL,  pi2 = NULL, var_est = FALSE,
                       var_method="lin_HB", strata = NULL, H = 256, lr = 0.001,
                       epochs = 10, batch_size = NULL, validation_split = 0.2,
                       optimizer = NULL, loss = "MSE", verbose = 1){
  #load keras
  require(keras)
  #Make sure the var_method is valid
  if(!is.element(var_method, c("lin_HB", "lin_HH", "lin_HTSRS", "lin_HT", "bootstrap_SRS"))){
    message("Variance method input incorrect. It has to be \"lin_HB\", \"lin_HH\", \"lin_HT\", \"lin_HTSRS\", or \"bootstrap_SRS\".")
    return(NULL)
  }
  #Check that x_sample is a df
  if(class(x_sample) != "data.frame"){
    message("x_sample must be a data.frame.")
    return(NULL)
  }
  #Check that y is numeric
  if(!(typeof(y) %in% c("numeric", "integer", "double"))){
    stop("Must supply numeric y.  For binary variable, convert to 0/1's.")
  }
  #Convert y to a vector
  y <- as.vector(y)
  
  #Make sure y is complete
  if(NA %in% y){
    message("Must supply complete cases for y")
    return(NULL)
  }
  #Check that the number of observations are consistent
  if(nrow(x_sample) != length(y)){
    message("y and x_sample must have same number of observations.")
    return(NULL)
  }
  #sample size
  n <- length(y)
  
  #population size
  N <- dim(x_pop)[1]
  
  #Check on inclusion probabilities and create weight=inverse inclusion probabilities
  if(is.null(pi)){
    message("Assuming simple random sampling")
  }
  #Check for missing data:
  if(FALSE %in% complete.cases(x_sample) || FALSE %in% complete.cases(x_pop)){
    if(FALSE %in% complete.cases(x_sample)){
      message("Must supply complete cases for x_sample")
    }
    if(FALSE %in% complete.cases(x_pop)){
      message("Must supply complete cases for x_pop")
    }
    return(NULL)
  }    
  # Create pi
  if (is.null(pi)) {
    pi <- rep(length(y)/N, length(y))
  }
  
  #weight: inverse first order inclusion probabilities
  weights <- as.vector(pi^(-1))
  
  #Make sure x_pop and x_sample have the same columns (in same order)
  x_pop <- x_pop[names(x_sample)]
  
  #Create formula and model matrix
  f <- as.formula(paste("y ~ ", paste(names(x_sample), collapse= "+")))
  model_mat <- model.matrix(f, x_sample)
  
  #build feed forward NN
  ffnn <- keras_model_sequential() %>%
    layer_dense(units = H, activation = 'relu', input_shape =  (ncol(x_sample)+1),
                kernel_initializer = initializer_glorot_normal()) %>%
    layer_dense(units = H, activation = 'relu',
                kernel_initializer = initializer_glorot_normal()) %>%
    layer_dense(units = H, activation = 'relu',
                kernel_initializer = initializer_glorot_normal()) %>%
    layer_dense(units = 1)
  #compile neural network
  ffnn <- ffnn %>% compile(loss = loss,
                           optimizer = optimizer_rmsprop(lr = lr, rho = 0.9),
                           metrics = list("mean_absolute_error"))
  #train neural network
  fit(ffnn, x = model_mat, y = y, epochs = epochs,
      batch_size = ifelse(is.null(batch_size), nrow(x_sample), batch_size), 
      validation_split = validation_split, verbose = verbose)
  
  #find estimates
  y_hat_sample <- predict(ffnn, model_mat) %>% as.vector()
  if(verbose != 0) {
    cat("fitting population level data", "\n")
  }
  y_hat_pop <- predict(ffnn, cbind(1, as.matrix(x_pop))) %>% as.vector()
  
  #calculating the total estimate for y
  t <- sum(weights * (y - y_hat_sample)) + sum(y_hat_pop)
  
  if(var_est==TRUE){
    if(var_method != "bootstrap_SRS"){
      e <- y - y_hat_sample
      varEst <- varMase(y = e, pi = pi, pi2 = pi2, method = var_method, N = N,
                        strata = strata)
    }
    return(list( pop_total = as.numeric(t),
                 pop_mean = as.numeric(t)/N,
                 pop_total_var=varEst,
                 pop_mean_var=varEst/N^2,
                 formula = f,
                 model = ffnn,
                 y_hat_sample = y_hat_sample,
                 y_hat_pop = y_hat_pop) %>%
             gregify())
  }
  else{
    return(list( pop_total = as.numeric(t),
                 pop_mean = as.numeric(t)/N,
                 formula = f,
                 model = ffnn,
                 y_hat_sample = y_hat_sample,
                 y_hat_pop = y_hat_pop) %>%
             gregify())
  }
}