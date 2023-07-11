#' Compute a modified generalized regression estimator
#' 
#' Calculates a modified generalized regression estimator for a finite population mean/proportion or total based on sample data collected from a complex sampling design and auxiliary population data.  
#' 
#' @param xsample A data frame of the auxiliary data in the sample.
#' @param xpop A data frame of population level auxiliary information.  It must contain all of the names from xsample. If datatype = "raw", must contain unit level data.  If datatype = "totals" or "means", then contains one row of aggregated, population totals or means for the auxiliary data and must include a column labeled N with the population sizes for each domain. Default is "raw".
#' @param domains A vector of the specific domain that each row of xsample belongs to.
#' @param domain_name A string that specifies the name of the column that contains the domain values in xpop.
#' @param estimation_domains A vector of domain values over which to produce estimates. If NULL, estimation will occur over all of the domains included in xpop.
#' @param pi First order inclusion probabilities.
#' @param pi2 Second order inclusion probabilities.
#' @param datatype A string that specifies the form of population auxiliary data. The possible values are "raw", "totals" or "means" for whether the user is providing population data at the unit level, aggregated to totals, or aggregated to means.  Default is "raw".
#' @param model A string that specifies the regression model to utilize. Options are "linear" or "logistic".
#' @param N The total population size.
#' 
#' @export modified_greg
#' @import survey
#' @import glmnet
#' @import boot
#' @importFrom stats model.matrix predict quasibinomial var
#' @include varMase.R
#' 
#' 
modified_greg <- function(y,
                          xsample,
                          xpop,
                          domains,
                          domain_name,
                          estimation_domains = NULL,
                          pi = NULL, 
                          pi2 = NULL,
                          datatype = "raw",
                          model = "linear",
                          var_est = F,
                          var_method = "LinHB", 
                          N = NULL,
                          B = 1000) {
  
  if (datatype != "raw" && !("N" %in% names(xpop))) {
    stop("xpop must contain a column for population size by domain called 'N' when datatype != raw.")
  }

  if (is.null(N)) {
    stop("Must supply total population size N")
  }

  if(!(typeof(y) %in% c("numeric", "integer", "double"))){
    stop("Must supply numeric y.  For binary variable, convert to 0/1's.")
  }

  if(!is.element(var_method, c("LinHB", "LinHH", "LinHTSRS", "LinHT", "bootstrapSRS"))){
    message("Variance method input incorrect. It has to be \"LinHB\", \"LinHH\", \"LinHT\", \"LinHTSRS\", or \"bootstrapSRS\".")
    return(NULL)
  }

  if(!is.element(model, c("linear","logistic"))){
    message("Method input incorrect, has to be either \"linear\" or \"logistic\"")
    return(NULL)
  }

  if(model == "logistic" && datatype != "raw"){
    message("Must supply the raw population data to fit the logistic regression estimator.")
    return(NULL)

  }

  if(!is.element(datatype, c("raw","totals", "means"))){
    message("datatype input incorrect, has to be either \"raw\", \"totals\" or \"means\"")
    return(NULL)
  }
  
  pop_unique_domains <- unique(xpop[[domain_name]])
  samp_unique_domains <- unique(domains)
  
  if(!setequal(pop_unique_domains, samp_unique_domains)) {
    stop("`domains` must contain all the same unique domain values as xpop  ")
  }
  
  # if domain_labels is null run over everything
  # probably need to drop NAs from xpop in a preprocessing step
  # need to check that all domain_labels exist in both xpop and xsample

  if(is.null(pi)){
    message("Assuming simple random sampling")
  } 
  
  if (is.null(pi)) {
    pi <- rep(length(y)/N, length(y))
  }
  
  weight <- as.vector(pi^(-1))
  
  y <- as.vector(y)
  
  
  if (is.null(estimation_domains)) {
    estimation_domains <- pop_unique_domains
  }
  
  # creating a vector of common auxiliary variable names
  common_pred_vars <- intersect(names(xsample), names(xpop))
  
  # creating the design matrix for entire xsample
  xsample_d <- model.matrix(~., data = data.frame(xsample[common_pred_vars]))
  xsample <- cbind(data.frame(xsample_d[,-1, drop = FALSE]), domains)
  names(xsample) <- c(colnames(xsample_d[,-1, drop = FALSE]), domain_name)
  xsample_dt <- t(xsample_d) 
  
  if (model == "linear") {
    
    if (datatype == "raw"){
      
      # only using columns that occur in xsample too
      xpop_subset <- xpop[common_pred_vars]
      
      # expand factors and interaction terms and then re-add the domain column
      xpop <- cbind(data.frame(model.matrix(~. -1, data = data.frame(xpop_subset))), xpop[domain_name])
      xpop$`N` <- 1
      
      # sum auxiliary pop values by domain
      bydomain_formula <- as.formula(paste0(". ~", domain_name))
      xpop_d <- aggregate(bydomain_formula, xpop, FUN = sum)
      # move N column to the front (where it would normally be)
      xpop_d <- xpop_d[, c(ncol(xpop_d), 1:(ncol(xpop_d) - 1))]
      
    }
    if (datatype == "totals"){
      
      # need N for each specific domain
      # assume they are there?
      xpop_d <- xpop[ ,c("N", common_pred_vars, domain_name)]
      
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

      domain_indic_vec <- as.integer(xsample[domain_name] == domain_id)

      xpop_domain <- xpop_d[xpop_d[domain_name] == domain_id, , drop = FALSE]
      xpop_d_domain <- unlist(xpop_domain[-which(names(xpop_domain) == domain_name)])
      xsample_domain <- xsample[xsample[domain_name] == domain_id, , drop = FALSE]
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
          
        } else if (var_method == "boostrapSRS"){
          
          # need to implement
          
        }
        
        return(list(
          domain = domain_id,
          pop_total = as.numeric(t),
          pop_mean = as.numeric(t)/as.numeric(domain_N),
          pop_total_var = as.numeric(varEst),
          pop_mean_var = as.numeric(varEst)/as.numeric(domain_N^2) 
        ))
        
      } else {
        
        return(list(
          domain = domain_id,
          pop_total = as.numeric(t),
          pop_mean = as.numeric(t)/as.numeric(domain_N)
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
    xpop <- cbind(data.frame(model.matrix(~., data = data.frame(xpop_subset))), xpop[domain_name])
    
    # need N by domain
    bydomain_formula <- as.formula(paste0(names(xpop)[1], "~", domain_name))
    xpop_sums <- aggregate(bydomain_formula, xpop, FUN = sum)
    
    # preparing logistic model
    xsample_preds <- xsample[ , !(names(xsample) %in% domain_name)]
    dat <- data.frame(y, weight, xsample_preds)
    colnames(dat) <- c("y", "weight", names(xsample_preds))
    f <- paste(names(dat)[1], "~", paste(names(dat)[-c(1,2)], collapse = " + "))
    s_design <- survey::svydesign(ids = ~1, weights = ~weight, data = dat)
    mod <- survey::svyglm(f, design = s_design, family = quasibinomial())
    
    by_domain_logistic <- function(domain_id) {
      
      domain_indic_vec <- as.integer(xsample[domain_name] == domain_id)
      
      xpop_domain <- xpop[xpop[domain_name] == domain_id, , drop = FALSE]
      xsample_domain <- xsample[xsample[domain_name] == domain_id, , drop = FALSE]
      domain_N <- xpop_sums[xpop_sums[domain_name] == domain_id, ,drop = FALSE][names(xpop)[1]]
      
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
          
          # need to implement
          
        }
        
        return(list(
          domain = domain_id,
          pop_total = as.numeric(t),
          pop_mean = as.numeric(t)/as.numeric(domain_N),
          pop_total_var = as.numeric(varEst),
          pop_mean_var = as.numeric(varEst)/as.numeric(domain_N^2) 
        ))
        
      } else {
        
        return(list(
          domain = domain_id,
          pop_total = as.numeric(t),
          pop_mean = as.numeric(t)/as.numeric(domain_N)
        ))
        
      }
      
    }
    
    res <- lapply(domain_labels, FUN = by_domain_logistic)
    
  }
  
  return(res)
  
}



  
t <- modified_greg(
  y = y,
  xsample = xsample,
  xpop = xpop,
  domains = domains,
  domain_name = "COUNTYFIPS",
  model = "linear",
  datatype = "means",
  N = sum(xpop$N),
  var_est = T
)
  

library(survey)
data(api)

modified_greg(y = as.numeric(apisrs$awards) - 1,
    xsample = apisrs[c("col.grad", "api00")],
    xpop = apipop[c("col.grad", "api00", "stype")],
    domains = apisrs$stype,
    domain_labels = c("H", "E"),
    model = "logistic",
    pi = apisrs$pw^(-1),
    var_est = T,
    N = nrow(apipop))


greg(y = as.numeric(apisrs$awards) - 1,
     xsample = apisrs[c("col.grad", "api00")],
     xpop = apipop[c("col.grad", "api00")],
     pi = apisrs$pw^(-1),
     model = "logistic",
     var_est = TRUE)

