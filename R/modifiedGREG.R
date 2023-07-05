library(tidyverse)
library(mase)
unitzonal <- read_rds("/Users/joshuayamamoto/Downloads/Oregon/unitzonal.rds")
tsumdatp <- read_rds("/Users/joshuayamamoto/Downloads/Oregon/tsumdatp.rds") 
pltassign <- read_rds("/Users/joshuayamamoto/Downloads/Oregon/pltassgn.rds") 

sampx <- tsumdatp %>% 
  left_join(pltassign, by = "CN") %>% 
  select(tcc16, elev, COUNTYCD) %>% 
  rename(COUNTYFIPS = COUNTYCD) %>% 
  drop_na() %>% 
  mutate(COUNTYFIPS = paste0("410", COUNTYFIPS))

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
  

popx <- unitzonal %>% 
  select(tcc16, elev, COUNTYFIPS, npixels) %>% 
  drop_na() %>% 
  rename(N = npixels)




# notes
# xsample and xpop must contain a label/column that tells us which domain each point belongs to
# these labels must match up with the labels given in the domain_labels argument
# the harder case is when the datatype is raw, but can still be managed

modified_greg <- function(y, xsample, xpop, domain, domain_labels, pi = NULL, model = "linear", datatype = "raw", N = NULL) {
  
  # probably need to drop NAs from xpop in a preprocessing step
  # need to check that all domain_labels exist in both xpop and xsample
  # xpop and xsample must be dataframes 
  
  # convert pi into diagonal matrix format
  if (is.null(pi)) {
    pi <- rep(length(y)/N, length(y))
  }
  
  #weight: inverse first order inclusion probabilities
  weight <- as.vector(pi^(-1))
  
  # creating a vector of common auxiliary variable names
  common_vars <- intersect(names(xsample), names(xpop))
  common_pred_vars <- common_vars[!common_vars %in% domain]
  
  # only continue using domains we want to estimate
  #xpop <- xpop[xpop[domain] %in% domain_labels]
  
  # creating the design matrix for entire xsample
  xsample_d <- model.matrix(~., data = data.frame(xsample[common_pred_vars]))
  xsample <- cbind(data.frame(xsample_d[,-1, drop = FALSE]), xsample[domain])
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
      xpop_d <- xpop[ , c("N", common_pred_vars, domain)]
      
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

    by_domain <- function(y, weight, xpop_d, domain_id) {

      domain_indic_vec <- as.integer(xsample[domain] == domain_id)

      xpop_aoi <- xpop_d[xpop_d[domain] == domain_id, ,drop = FALSE]
      xpop_d_aoi <- unlist(xpop_aoi[-which(names(xpop_aoi) == domain)])
      xsample_aoi <- xsample[xsample[domain] == domain_id, ,drop = FALSE]
      xsample_d_aoi <- model.matrix(~., data = data.frame(xsample_aoi[common_pred_vars]))
      xsample_dt_aoi <- t(xsample_d_aoi)
      weights_aoi <- weight[domain_indic_vec]

      w <- as.matrix(
        weight*domain_indic_vec + (
        t(as.matrix(xpop_d_aoi) - xsample_dt_aoi %*% weights_aoi) %*%
          constant_component1
        ) %*%
        constant_component2
        )
      
      pop_total <- w %*% y
      pop_mean <- pop_total/unlist(xpop_aoi["N"])
      
      return(list(
        pop_total = as.numeric(pop_total),
        pop_mean = as.numeric(pop_mean)
      ))

    }

    # run by_domain function over domain_labels argument

    res <- by_domain(y, weight, xpop_d, domain_labels[1])

    #res <- lapply(domain_labels, FUN = by_domain(y, weights, xpop_d, x))
    
    
  }
  return(res)
  
}

t <- modified_greg(
  y = y,
  xsample = sampx,
  xpop = popx,
  datatype = "means",
  domain = "COUNTYFIPS",
  domain_labels = c("41025"),
  N = sum(popx$N)
  )

greg_est <- greg(
  y = y_aoi,
  xsample = sampx[sampx$COUNTYFIPS == "41025", ][c("tcc16", "elev")],
  xpop = popx[popx$COUNTYFIPS == "41025", ][c("tcc16", "elev")],
  datatype = "means",
  N = popx[popx$COUNTYFIPS == "41025", ]$N
  )

greg_est$pop_mean




