library(tidyverse)
library(mase)
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

popx <- unitzonal %>% 
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

# notes
# xsample and xpop must contain a label/column that tells us which domain each point belongs to
# these labels must match up with the labels given in the domain_labels argument
# the harder case is when the datatype is raw, but can still be managed

modified_greg <- function(y, xsample, xpop, domain, domain_labels = NULL, pi = NULL, model = "linear", pi2 = NULL, var_est = F, var_method = "LinHB", datatype = "raw", N = NULL) {
  
  # need N by domain if datatype != raw
  # force xpop to have N by domain
  # if domain_labels is null run over everything
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
  
  # only continue using domains we want to estimate on
  xpop <- xpop[xpop[[domain]] %in% domain_labels, ]
  
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
    betas <- solve(xsample_dt %*% diag(weight) %*% xsample_d) %*% (xsample_dt) %*% diag(weight) %*% y

    by_domain <- function(domain_id) {

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
      
      t <- w %*% y
      
      aoi_N <- as.numeric(unlist(xpop_aoi["N"]))
    
      if(var_est == TRUE) {
        if(var_method != "bootstrapSRS") {
          y_hat <- xsample_d_aoi %*% betas
          y_aoi <- y[domain_indic_vec]
          e <- y_aoi - y_hat
          varEst <- varMase(y = e, pi = pi[domain_indic_vec], pi2 = pi2, method = var_method, N = aoi_N)
        }
      }
      
      return(list(
        domain = domain_id,
        pop_total = as.numeric(t),
        pop_mean = as.numeric(t)/aoi_N,
        pop_total_var = as.numeric(varEst),
        pop_mean_var = as.numeric(varEst/aoi_N^2) 
      ))

    }
    

    # run by_domain function over domain_labels argument
    res <- lapply(domain_labels, FUN = by_domain)

  }
  return(res)
  
}


  
t <- modified_greg(
  y = y,
  xsample = sampx,
  xpop = popx,
  datatype = "means",
  domain = "COUNTYFIPS",
  domain_labels = c("41025", "41037"),
  N = sum(popx$N),
  var_est = T
)
  




## additive test ---------------------------------------------------------------

full <- modified_greg(
  y = y,
  xsample = sampx,
  xpop = popx,
  datatype = "means",
  domain = "COUNTYFIPS",
  domain_labels = popx$COUNTYFIPS,
  N = sum(popx$N),
  var_est = T
)

total <- 0
var <- 0
for(i in 1:36) {
  total <- total + full[[i]]$pop_total
  var <- var + full[[i]]$pop_total_var
}

test <- popx %>% 
  mutate(tcc16 = tcc16*N, elev = elev*N) %>% 
  summarise(tcc16 = sum(tcc16), elev = sum(elev))

greg_est <- greg(
  y = y,
  xsample = sampx[c("tcc16", "elev")],
  xpop = test,
  datatype = "totals",
  N = sum(popx$N),
  var_est = F
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









# making sure answers are sensible with regular greg ---------------------------

greg_est <- greg(
  y = y_aoi,
  xsample = sampx[sampx$COUNTYFIPS == "41037", ][c("tcc16", "elev")],
  xpop = popx[popx$COUNTYFIPS == "41037", ][c("tcc16", "elev")],
  datatype = "means",
  N = popx[popx$COUNTYFIPS == "41037", ]$N,
  var_est = T
  )






