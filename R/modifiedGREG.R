library(tidyverse)
library(mase)

# data preprocessing -----------------------------------------------------------

unitzonal <- read_rds("/Users/joshuayamamoto/Downloads/Oregon/unitzonal.rds")
tsumdatp <- read_rds("/Users/joshuayamamoto/Downloads/Oregon/tsumdatp.rds") 
pltassign <- read_rds("/Users/joshuayamamoto/Downloads/Oregon/pltassgn.rds") 

sampx <- tsumdatp %>% 
  left_join(pltassign, by = "CN") %>% 
  select(tcc16, elev, COUNTYCD) %>% 
  drop_na() %>% 
  mutate(COUNTYCD = as.character(COUNTYCD)) %>% 
  mutate(COUNTYFIPS = case_when(
    str_length(COUNTYCD) < 2 ~ paste0("4100", COUNTYCD),
    T ~ paste0("410", COUNTYCD)
  ))

domains <- sampx$COUNTYFIPS

xsample <- sampx %>% 
  select(tcc16, elev)

xpop <- unitzonal %>% 
  select(tcc16, elev, COUNTYFIPS, npixels) %>% 
  drop_na() %>% 
  rename(N = npixels)

y_aoi <- tsumdatp %>% 
  left_join(pltassign, by = "CN") %>% 
  select(tcc16, elev, COUNTYCD, BA_TPA_live_ADJ) %>%  
  drop_na() %>% 
  filter(COUNTYCD == "25") %>% 
  select(BA_TPA_live_ADJ) %>% 
  pull()

y <- tsumdatp %>% 
  left_join(pltassign, by = "CN") %>% 
  select(tcc16, elev, COUNTYCD, BA_TPA_live_ADJ) %>%  
  drop_na() %>% 
  select(BA_TPA_live_ADJ) %>% 
  pull()


## MODIFIED GREG FUNCTION CODE--------------------------------------------------



modified_greg <- function(y,
                          xsample,
                          xpop,
                          domains,
                          domain_labels = NULL,
                          pi = NULL, 
                          pi2 = NULL,
                          model = "linear",
                          var_est = F,
                          var_method = "LinHB", 
                          datatype = "raw",
                          N = NULL,
                          B = 1000) {
  
  if (datatype != "raw" && !("N" %in% names(xpop))) {
    stop("xpop must contain a column for population size by domain called 'N' when datatype != NULL.")
  }


  if (is.null(N)) {
    stop("Must supply total population size N")
  }

  if(!(typeof(y) %in% c("numeric", "integer", "double"))){
    stop("Must supply numeric y.  For binary variable, convert to 0/1's.")
  }


  #Make sure the var_method is valid
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
  
  # if domain_labels is null run over everything
  # probably need to drop NAs from xpop in a preprocessing step
  # need to check that all domain_labels exist in both xpop and xsample

  
  if(is.null(pi)){
    message("Assuming simple random sampling")
  } 
  
  # convert pi into diagonal matrix format
  if (is.null(pi)) {
    pi <- rep(length(y)/N, length(y))
  }
  
  weight <- as.vector(pi^(-1))
  
  y <- as.vector(y)
  
  # extracting the name of the domain column
  # potentially just ask for this as an argument to the function
  domain <- names(
    which(
      apply(xpop, 2, function(x) sum(unique(domains) %in% x)) == length(unique(domains))
    )
  )
  
  # creating a vector of common auxiliary variable names
  common_pred_vars <- intersect(names(xsample), names(xpop))
  
  # creating the design matrix for entire xsample
  xsample_d <- model.matrix(~., data = data.frame(xsample[common_pred_vars]))
  xsample <- cbind(data.frame(xsample_d[,-1, drop = FALSE]), domains)
  names(xsample) <- c(colnames(xsample_d[,-1, drop = FALSE]), domain)
  xsample_dt <- t(xsample_d) 
  
  if (model == "linear") {
    
    if (datatype == "raw"){
      
      # only using columns that occur in xsample too
      xpop_subset <- xpop[common_pred_vars]
      
      # expand factors and interaction terms and then re-add the domain column
      xpop <- cbind(data.frame(model.matrix(~. -1, data = data.frame(xpop_subset))), xpop[domain])
      xpop$`N` <- 1
      
      # sum auxiliary pop values by domain
      bydomain_formula <- as.formula(paste0(". ~", domain))
      xpop_d <- aggregate(bydomain_formula, xpop, FUN = sum)
      # move N column to the front (where it would normally be)
      xpop_d <- xpop_d[, c(ncol(xpop_d), 1:(ncol(xpop_d) - 1))]
      
    }
    if (datatype == "totals"){
      
      # need N for each specific domain
      # assume they are there?
      xpop_d <- xpop[ ,c("N", common_pred_vars, domain)]
      
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

      domain_indic_vec <- as.integer(xsample[domain] == domain_id)

      xpop_domain <- xpop_d[xpop_d[domain] == domain_id, ,drop = FALSE]
      xpop_d_domain <- unlist(xpop_domain[-which(names(xpop_domain) == domain)])
      xsample_domain <- xsample[xsample[domain] == domain_id, ,drop = FALSE]
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
      
      domain_N <- as.numeric(unlist(xpop_domain["N"]))
    
      if(var_est == TRUE) {
        if(var_method != "bootstrapSRS") {
          
          y_hat <- xsample_d_domain %*% betas
          y_domain <- y[which(domain_indic_vec == 1)]
          e <- y_domain - y_hat
          varEst <- varMase(y = e, pi = pi[which(domain_indic_vec == 1)], pi2 = pi2, method = var_method, N = aoi_N)
          
        } else if (var_method == "boostrapSRS"){
          
          # need to implement
          
        }
        
        return(list(
          domain = domain_id,
          pop_total = as.numeric(t),
          pop_mean = as.numeric(t)/domain_N,
          pop_total_var = as.numeric(varEst),
          pop_mean_var = as.numeric(varEst)/as.numeric(domain_N^2) 
        ))
        
        
      } else {
        
        return(list(
          domain = domain_id,
          pop_total = as.numeric(t),
          pop_mean = as.numeric(t)/domain_N
        ))
        
      }
      
    }

    # run by_domain function over domain_labels argument
    res <- lapply(domain_labels, FUN = by_domain_linear)

  } else if (model == "logistic") {
    
    if(length(levels(as.factor(y))) != 2){
      stop("Function can only handle categorical response with two categories.")
    }
    
    xpop_subset <- xpop[common_pred_vars]
    xpop <- cbind(data.frame(model.matrix(~., data = data.frame(xpop_subset))), xpop[domain])
    
    # need N by domain
    bydomain_formula <- as.formula(paste0(names(xpop)[1], "~", domain))
    xpop_sums <- aggregate(bydomain_formula, xpop, FUN = sum)
    
    # preparing logistic model
    xsample_preds <- xsample[ ,!(names(xsample) %in% domain)]
    dat <- data.frame(y, weight, xsample_preds)
    colnames(dat) <- c("y", "weight", names(xsample_preds))
    f <- paste(names(dat)[1], "~", paste(names(dat)[-c(1,2)], collapse = " + "))
    s_design <- survey::svydesign(ids = ~1, weights = ~weight, data = dat)
    mod <- survey::svyglm(f, design = s_design, family = quasibinomial())
    
    by_domain_logistic <- function(domain_id) {
      
      domain_indic_vec <- as.integer(xsample[domain] == domain_id)
      
      xpop_domain <- xpop[xpop[domain] == domain_id, , drop = FALSE]
      xsample_domain <- xsample[xsample[domain] == domain_id, , drop = FALSE]
      domain_N <- xpop_sums[xpop_sums[domain] == domain_id, ,drop = FALSE][names(xpop)[1]]
      
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
  model = "linear",
  datatype = "means",
  domain_labels = c("41025", "41037"),
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


## additive test (hold for totals but not variances) ---------------------------------------------------------------

full <- modified_greg(
  y = y,
  xsample = xsample,
  xpop = xpop,
  datatype = "means",
  domain = domains,
  domain_labels = xpop$COUNTYFIPS,
  N = sum(xpop$N),
  var_est = T
)

full %>% 
  map_df(.f = function(x) data.frame(total = x$pop_total, var = x$pop_total_var)) %>% 
  summarise(across(everything(), sum))

xpop_totals <- xpop %>% 
  mutate(tcc16 = tcc16*N, elev = elev*N) %>% 
  summarise(tcc16 = sum(tcc16), elev = sum(elev))

greg_est <- greg(
  y = y,
  xsample = xsample[c("tcc16", "elev")],
  xpop = xpop_totals,
  datatype = "totals",
  N = sum(xpop$N),
  var_est = T
)

data.frame(
  total = greg_est$pop_total,
  var = greg_est$pop_total_var
)


## speed test ------------------------------------------------------------------

test_speed <- function(domain_labels, n_domains) {
  
  a <- system.time(
    modified_greg(
      y = y,
      xsample = sampx,
      xpop = popx,
      datatype = "means",
      domain = "COUNTYFIPS",
      domain_labels = domain_labels,
      N = sum(popx$N),
      var_est = T
    )
  )
  
  return(
    data.frame(
      time = a[3],
      n_domains = n_domains
    )
  )
  
}

set.seed(11)
full_speed_res <- data.frame()

for(j in 1:15) {
  order <- sample(popx$COUNTYFIPS, length(popx$COUNTYFIPS))
  domain_label_sets <- list()
  
  for(i in 1:length(popx$COUNTYFIPS)) {
    domain_label_sets[[i]] <- order[1:i]
  }
  speed_res <-  map2_df(.x = domain_label_sets, .y = 1:36, .f =  ~ test_speed(.x, .y))
  rownames(speed_res) <- NULL
  full_speed_res <- rbind(full_speed_res, speed_res)
  print(paste0("done with ", j))
}



# over 15 "sim" runs

full_speed_res %>% 
  group_by(n_domains) %>% 
  summarise(avg_time = mean(time)) %>% 
  ggplot(aes(x = n_domains, y = avg_time)) +
  geom_line() +
  geom_point(color = "cyan4") +
  theme_bw() +
  ylim(c(0, 3))









