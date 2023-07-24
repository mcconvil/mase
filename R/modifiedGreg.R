#' Compute a modified generalized regression estimator
#' 
#' Calculates a modified generalized regression estimator for a finite population mean/proportion or total based on sample data collected from a complex sampling design and auxiliary population data.  
#' 
#' @param y A vector of the response values from the sample
#' @param xsample A data frame of the auxiliary data in the sample.
#' @param xpop A data frame of population level auxiliary information.  It must contain all of the names from xsample. If datatype = "raw", must contain unit level data.  If datatype = "totals" or "means", then contains one row of aggregated, population totals or means for the auxiliary data and must include a column labeled N with the population sizes for each domain. Default is "raw".
#' @param domains A vector of the specific domain that each row of xsample belongs to.
#' @param pi First order inclusion probabilities.
#' @param pi2 Second order inclusion probabilities.
#' @param datatype A string that specifies the form of population auxiliary data. The possible values are "raw", "totals" or "means" for whether the user is providing population data at the unit level, aggregated to totals, or aggregated to means.  Default is "raw".
#' @param model A string that specifies the regression model to utilize. Options are "linear" or "logistic".
#' @param var_est A logical value that specifies whether variance estimation should be performed.
#' @param var_method A string that specifies the variance method to utilize. 
#' @param modelselect A logical for whether or not to run lasso regression first and then fit the model using only the predictors with non-zero lasso coefficients. Default is FALSE.  
#' @param lambda A string specifying how to tune the lasso hyper-parameter.  Only used if modelselect = TRUE and defaults to "lambda.min". The possible values are "lambda.min", which is the lambda value associated with the minimum cross validation error or "lambda.1se", which is the lambda value associated with a cross validation error that is one standard error away from the minimum, resulting in a smaller model.
#' @param domain_col_name A string that specifies the name of the column that contains the domain values in xpop.
#' @param estimation_domains A vector of domain values over which to produce estimates. If NULL, estimation will be performed over all of the domains included in xpop.
#' @param N The total population size.
#' 
#' @export modifiedGreg
#' @import survey
#' @import glmnet
#' @import boot
#' @importFrom stats model.matrix predict quasibinomial var aggregate 
#' @include varMase.R

modifiedGreg <- function(y,
                         xsample,
                         xpop,
                         domains,
                         pi = NULL, 
                         pi2 = NULL,
                         datatype = "raw",
                         model = "linear",
                         var_est = F,
                         var_method = "LinHB", 
                         modelselect = FALSE,
                         lambda = "lambda.min",
                         domain_col_name = NULL,
                         estimation_domains = NULL,
                         N = NULL,
                         B = 1000) {

  if (!(typeof(y) %in% c("numeric", "integer", "double"))) {
    stop("Must supply numeric y.  For binary variable, convert to 0/1's.")
  }

  if (!is.element(var_method, c("LinHB", "LinHH", "LinHTSRS", "LinHT", "bootstrapSRS"))) {
    message("Variance method input incorrect. It has to be \"LinHB\", \"LinHH\", \"LinHT\", \"LinHTSRS\", or \"bootstrapSRS\".")
    return(NULL)
  }

  if (!is.element(model, c("linear","logistic"))) {
    message("Method input incorrect, has to be either \"linear\" or \"logistic\"")
    return(NULL)
  }

  if (model == "logistic" && datatype != "raw") {
    message("Must supply the raw population data to fit the logistic regression estimator.")
    return(NULL)

  }

  if (!is.element(datatype, c("raw","totals", "means"))) {
    message("datatype input incorrect, has to be either \"raw\", \"totals\" or \"means\"")
    return(NULL)
  }
  
  if (is.null(N)) {
    if (datatype == "raw") {
      N <- nrow(xpop)
    } else {
      N <- sum(xpop$N)
    }
  }
  
  if (!is.null(N) && datatype != "raw"  && sum(xpop$N) != N) {
    stop("User inputted N does not equal the sum of the domain level population sizes in xpop.")
  }
  
  if (datatype != "raw" && !("N" %in% names(xpop))) {
    stop("xpop must contain a column for population size by domain called 'N' when datatype != raw.")
  }
  
  if (!all(names(xsample) %in% names(xpop))) {
    stop("All of the column names in `xsample` must exist in `xpop`.")
  }
  
  # xpop should only have either one or two more columns than xsample
  ncol_diff <- ncol(xpop) - ncol(xsample)
  if (datatype == "raw" && (ncol_diff != 1)) {
    stop("Incorrect number of columns in either xpop or xsample. When datatype = \"raw\" xpop should only contain columns with the same names as xsample and a column with the domains")
  } else if (is.element(datatype, c("means", "totals")) && (ncol_diff != 2)) {
    stop("Incorrect number of columns in either xpop or xsample. When datatype != \"raw\" xpop should only contain columns with the same names as xsample as well as column with the domains and a column with the population sizes.")
  }
  

  if (is.null(domain_col_name)) {
    
    if (datatype == "raw") {
      domain_col_name <- setdiff(names(xpop), names(xsample))
    } else {
      domain_col_name <- setdiff(names(xpop), c(names(xsample), "N"))
    }
    
    message(paste0("domain_col_name is not directly specified. ", domain_col_name, " is being used."))
    
  }
  
  pop_unique_domains <- unique(xpop[[domain_col_name]])
  samp_unique_domains <- unique(domains)
  
  if (!setequal(pop_unique_domains, samp_unique_domains)) {
    stop("`domains` must contain all the same unique domain values as xpop  ")
  }

  if (is.null(pi)) {
    message("Assuming simple random sampling")
  } 
  
  if (is.null(pi)) {
    pi <- rep(length(y)/N, length(y))
  }
  
  weight <- as.vector(pi^(-1))
  
  if (is.null(estimation_domains)) {
    estimation_domains <- pop_unique_domains
  }
  
  # creating a vector of common auxiliary variable names
  common_pred_vars <- intersect(names(xsample), names(xpop))
  
  # creating the design matrix for entire xsample
  xsample_d <- model.matrix(~., data = data.frame(xsample[common_pred_vars]))
  xsample <- cbind(data.frame(xsample_d[,-1, drop = FALSE]), domains)
  names(xsample) <- c(colnames(xsample_d[,-1, drop = FALSE]), domain_col_name)
  xsample_dt <- t(xsample_d) 
  
  # variable selection
  
  if (modelselect == TRUE) {
    
    if(model == "linear"){
      
      fam <- "gaussian"
      
    } else{
      
      fam <- "binomial"
      
    } 
    
    x <- xsample[ , -which(names(xsample) == domain_col_name)]
    cv <- cv.glmnet(x = as.matrix(x), y = y, alpha = 1, weights = weight, nfolds = 10, family = fam, standardize = FALSE)
    
    if(lambda=="lambda.min"){
      
      lambda.opt <- cv$lambda.min
      
    }
    if(lambda=="lambda.1se"){
      
      lambda.opt <- cv$lambda.1se
      
    }
    
    pred_mod <- glmnet(x = as.matrix(x), y = y, alpha = 1, family = fam, standardize = FALSE, weights = weight)
    lasso_coef <- predict(pred_mod, type = "coefficients", s = lambda.opt)[1:dim(xsample_d)[2],]
    coef_select <- names(lasso_coef[lasso_coef != 0])[-1]
    
    if(length(coef_select) == 0){
      
      message("No variables selected in the model selection stage.  Fitting a HT estimator.")
      
      if(var_est == TRUE){
        
        HT <- horvitzThompson(y = y, pi = pi, N = N, pi2 = pi2, var_est = TRUE, var_method = var_method)
        
        return(list(pop_total = HT$pop_total, 
                    pop_mean = HT$pop_total/N,
                    pop_total_var = HT$pop_total_var, 
                    pop_mean_var = HT$pop_total_var/N^2, 
                    weights = as.vector(pi^(-1))))
      } else {
        
        HT <- horvitzThompson(y = y, pi = pi, N = N, pi2 = pi2, var_est = FALSE)
        
        return(list(pop_total = HT$pop_total, 
                    pop_mean = HT$pop_total/N,
                    weights = as.vector(pi^(-1))))
        
      }
      
    } else {
      
      xsample <- xsample[ , c(coef_select, domain_col_name), drop = FALSE]
      xsample_d <- model.matrix(~., data = xsample[ , coef_select, drop = FALSE])
      xsample_dt <- t(xsample_d) 
      xsample <- cbind(data.frame(xsample_d[,-1, drop=FALSE]), domains)
      names(xsample) <- c(colnames(xsample_d[,-1, drop = FALSE]), domain_col_name)
      
      
    }  
    
  }
  
  
  if (model == "linear") {
    
    if (datatype == "raw"){
      
      # only using columns that occur in xsample too
      xpop_subset <- xpop[common_pred_vars]
      
      # expand factors and interaction terms and then re-add the domain column
      xpop <- cbind(data.frame(model.matrix(~. -1, data = data.frame(xpop_subset))), xpop[domain_col_name])
      xpop$`N` <- 1
      
      # sum auxiliary pop values by domain
      bydomain_formula <- as.formula(paste0(". ~", domain_col_name))
      xpop_d <- aggregate(bydomain_formula, xpop, FUN = sum)
      # move N column to the front (where it would normally be)
      xpop_d <- xpop_d[, c(ncol(xpop_d), 1:(ncol(xpop_d) - 1))]
      
    }
    if (datatype == "totals"){
      
      # need N for each specific domain
      # assume they are there?
      xpop_d <- xpop[ ,c("N", common_pred_vars, domain_col_name)]
      
    }
    if (datatype == "means"){
      
      # need N for each specific domain
      # assume they are there?
      xpop[common_pred_vars] <- lapply(xpop[common_pred_vars], function(x) x*(xpop$N))
      xpop_d <- cbind(N = xpop$N, xpop[ ,!(names(xpop) %in% "N")])
      
    }
    
    # computing the pieces that remain the same across all domains
    constant_component1 <- solve(xsample_dt %*% diag(weight) %*% xsample_d)
    constant_component2 <- t(weight * xsample_d)
    betas <- solve(xsample_dt %*% diag(weight) %*% xsample_d) %*% (xsample_dt) %*% diag(weight) %*% y
    
    # internal function to compute estimates by domain
    by_domain_linear <- function(domain_id) {

      domain_indic_vec <- as.integer(xsample[domain_col_name] == domain_id)

      xpop_domain <- xpop_d[xpop_d[domain_col_name] == domain_id, , drop = FALSE]
      xpop_d_domain <- unlist(xpop_domain[-which(names(xpop_domain) == domain_col_name)])
      xsample_domain <- xsample[xsample[domain_col_name] == domain_id, , drop = FALSE]
      xsample_d_domain <- model.matrix(~., data = data.frame(xsample_domain[common_pred_vars]))
      xsample_dt_domain <- t(xsample_d_domain)
      weights_domain <- weight[which(domain_indic_vec == 1)]

      w <- as.matrix(
        weight*domain_indic_vec + (
        t(as.matrix(xpop_d_domain) - xsample_dt_domain %*% weights_domain) %*%
          constant_component1
        ) %*%
        constant_component2
        )
      
      t <- w %*% y
      
      domain_N <- unlist(xpop_domain["N"])
    
      if(var_est == TRUE) {
        
        if(var_method != "bootstrapSRS") {
          
          y_hat <- xsample_d_domain %*% betas
          y_domain <- y[which(domain_indic_vec == 1)]
          e <- y_domain - y_hat
          varEst <- varMase(y = e, pi = pi[which(domain_indic_vec == 1)], pi2 = pi2, method = var_method, N = domain_N)
          
        } else if (var_method == "bootstrapSRS"){
          
          dat <- cbind(as.data.frame(cbind(y, pi, xsample_d)), xsample[[domain_col_name]])
          names(dat) <- c("y", "pi", colnames(xsample_d), domain_col_name)
          print(str(dat))
          t_boot <- boot(data = dat,
                         statistic = modifiedGregt,
                         R = B,
                         strata = as.factor(dat[ , ncol(dat)]),
                         xpopd = xpop_d_domain,
                         weight = weight,
                         domain = domain_id,
                         domain_col_name = domain_col_name,
                         parallel = "multicore",
                         ncpus = 2)

          varEst <- var(t_boot$t)
          
        }
        
        return(list(
          domain = domain_id,
          domain_total = as.numeric(t),
          domain_mean = as.numeric(t)/as.numeric(domain_N),
          domain_total_var = as.numeric(varEst),
          domain_mean_var = as.numeric(varEst)/as.numeric(domain_N^2) 
        ))
        
      } else {
        
        return(list(
          domain = domain_id,
          domain_total = as.numeric(t),
          domain_mean = as.numeric(t)/as.numeric(domain_N)
        ))
        
      }
      
    }

    # run by_domain function over domain_labels argument
    res <- lapply(estimation_domains, FUN = by_domain_linear)

  } else if (model == "logistic") {
    
    if(length(levels(as.factor(y))) != 2){
      stop("Function can only handle categorical response with two categories.")
    }
    
    xpop_subset <- xpop[common_pred_vars]
    xpop <- cbind(data.frame(model.matrix(~., data = data.frame(xpop_subset))), xpop[domain_col_name])
    
    # need N by domain
    bydomain_formula <- as.formula(paste0(names(xpop)[1], "~", domain_col_name))
    xpop_sums <- aggregate(bydomain_formula, xpop, FUN = sum)
    
    # preparing logistic model
    xsample_preds <- xsample[ , !(names(xsample) %in% domain_col_name)]
    dat <- data.frame(y, weight, xsample_preds)
    colnames(dat) <- c("y", "weight", names(xsample_preds))
    f <- paste(names(dat)[1], "~", paste(names(dat)[-c(1,2)], collapse = " + "))
    s_design <- survey::svydesign(ids = ~1, weights = ~weight, data = dat)
    mod <- survey::svyglm(f, design = s_design, family = quasibinomial())
    
    by_domain_logistic <- function(domain_id) {
      
      domain_indic_vec <- as.integer(xsample[domain_col_name] == domain_id)
      
      xpop_domain <- xpop[xpop[domain_col_name] == domain_id, , drop = FALSE]
      xsample_domain <- xsample[xsample[domain_col_name] == domain_id, , drop = FALSE]
      domain_N <- xpop_sums[xpop_sums[domain_col_name] == domain_id, ,drop = FALSE][names(xpop)[1]]
      
      y_hats_U <- as.matrix(predict(mod, newdata = xpop_domain[ , -1], type = "response", family = quasibinomial()))
      y_hats_s <- as.matrix(predict(mod, newdata = xsample_domain, type = "response", family = quasibinomial()))
      
      y_domain <- y[which(domain_indic_vec == 1)]
      weights_domain <- weight[which(domain_indic_vec == 1)]
      
      t <- t(y_domain - y_hats_s) %*% weights_domain + sum(y_hats_U)
      
      if (var_est == TRUE) {
        if (var_method != "bootstrapSRS") {
          
          e <- y_domain - y_hats_s
          varEst <- varMase(y = e, pi = pi[which(domain_indic_vec == 1)], pi2 = pi2, method = var_method, N = domain_N)
          
        } else if (var_method == "bootstrapSRS") {
          
          # xpop_sums, xpop_domain
          
          dat <- cbind(as.data.frame(cbind(y, pi, xsample_d)), xsample[[domain_col_name]])
          names(dat) <- c("y", "pi", colnames(xsample_d), domain_col_name)
          
          t_boot <- boot(dat,
               modifiedLogisticGregt,
               R = B,
               strata = as.factor(dat[ , ncol(dat)]),
               xpopd = xpop_domain,
               weight = weight,
               domain = domain_id,
               domain_col_name = domain_col_name,
               lab = names(xpop)[1],
               parallel = "multicore",
               ncpus = 2)
          
          varEst <- var(t_boot$t)
          
        }
        
        return(list(
          domain = domain_id,
          domain_total = as.numeric(t),
          domain_mean = as.numeric(t)/as.numeric(domain_N),
          domain_total_var = as.numeric(varEst),
          domain_mean_var = as.numeric(varEst)/as.numeric(domain_N^2) 
        ))
        
      } else {
        
        return(list(
          domain = domain_id,
          domain_total = as.numeric(t),
          domain_mean = as.numeric(t)/as.numeric(domain_N)
        ))
        
      }
      
    }
    
    res <- lapply(estimation_domains, FUN = by_domain_logistic)
    
  }
  
  pop_res <- do.call(
    rbind, lapply(
      res, FUN = function(x) data.frame(
        pop_total = x$domain_total,
        pop_totalvar = x$domain_total_var
      )
    )
  ) |>
    colSums()
  
  return(list(res, pop_res))
  
}



  
