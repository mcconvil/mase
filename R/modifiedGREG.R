library(tidyverse)
unitzonal <- read_rds("/Users/joshuayamamoto/Downloads/Oregon/unitzonal.rds")
tsumdatp <- read_rds("/Users/joshuayamamoto/Downloads/Oregon/tsumdatp.rds") 
pltassign <- read_rds("/Users/joshuayamamoto/Downloads/Oregon/pltassgn.rds") 

sampx <- tsumdatp %>% 
  left_join(pltassign, by = "CN") %>% 
  select(tcc16, elev, COUNTYCD) %>% 
  rename(COUNTYFIPS = COUNTYCD) %>% 
  drop_na()

y <- tsumdatp$BA_TPA_live_ADJ

popx <- unitzonal %>% 
  select(tcc16, elev, COUNTYFIPS, npixels) %>% 
  drop_na() %>% 
  rename(N = npixels)




# notes
# xsample and xpop must contain a label/column that tells us which domain each point belongs to
# these labels must match up with the labels given in the domain_labels argument
# the harder case is when the datatype is raw, but can still be managed

modified_greg <- function(y, xsample, xpop, domain = "COUNTYFIPS", domain_labels, pi = NULL, model = "linear", datatype = "raw", N = NULL) {
  
  # probably need to drop NAs from xpop in a preprocessing step
  # xpop and xsample must be dataframes 
  
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
    constant_component1 <- solve(xsample.dt %*% diag(weights) %*% xsample.d)
    constant_component2 <- t(weights * xsample.d)

    by_domain <- function(weights, domain_indic_vec, xpop_d, domain_id) {
      
      xpop_aoi <- xpop_d[xpop_d[domain] == domain_id, drop = FALSE]
      xpop_d_aoi <- as.matrix(xpop_aoi[-which(names(xpop_aoi) == domain)])
      xsample_d_aoi <- as.matrix(xsample[xsample[domain] == domain_id, drop = FALSE])
      xsample_dt_aoi <- t(xsample_d_aoi)
      weights_aoi <- weights[domain_indic_vec]
      
      w <- as.matrix(weights*domain_indic_vec + (
        t(as.matrix(xpop_d_aoi) - xsample_dt_aoi %*% weights_aoi) %*%
          constant_component1
        ) %*%
        constant_component2)
      
      return(w)

    }

    
    
    
    
    # perform this computation by domain based on domain_labels argument
    # preferably use some functional programming style (e.g map, apply) or figure out a matrix algebra way to do it
    
    # w <- as.matrix(
    #   weights*ids + 
    #     (t(as.matrix(xpop_d_aoi) - xsample.dt_aoi %*% weights_aoi) %*%
    #        solve(xsample.dt %*% diag(weights) %*% xsample.d)) %*% 
    #     t(weights * xsample.d)
    #   )
    # 
    # t <- w %*% y
    
    
  }
  return(xpop_d)
  
}

modified_greg(y, sampx, popx, datatype = "means")

a <- data.frame(
  N = c(100, 200, 300),
  id = c(1,2,3),
  x1 = c(12, 14, 15),
  x2 = c(50, 60, 70)
)

a[a["id"] == 1, ]




