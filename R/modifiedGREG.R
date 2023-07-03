library(tidyverse)
unitzonal <- read_rds("/Users/joshuayamamoto/Downloads/Oregon/unitzonal.rds")
tsumdatp <- read_rds("/Users/joshuayamamoto/Downloads/Oregon/tsumdatp.rds") 
pltassign <- read_rds("/Users/joshuayamamoto/Downloads/Oregon/pltassgn.rds") 

sampx <- tsumdatp %>% 
  left_join(pltassign, by = "CN") %>% 
  select(tcc16, elev, COUNTYCD, BA_TPA_live_ADJ) %>% 
  drop_na()

popx <- unitzonal %>% 
  select(tcc16, elev, COUNTYFIPS, npixels) %>% 
  drop_na()


xpop <- data.frame(model.matrix(~.-1, data = data.frame(popx_ex)))
xpop <- dplyr::select(xpop, all_of(colnames(xsample)))
xpop_d <- model.matrix(~., data = xpop)

aggregate( . ~ domain, xpop, FUN = sum)

xpop <- popx_ex
xsample <- sampx

# notes
# xsample and xpop must contain a label/column that tells us which domain each point belongs to
# these labels must match up with the labels given in the domain_labels argument
# the harder case is when the datatype is raw, but can still be managed

modified_greg <- function(y, xsample, xpop, domain, domain_labels, pi = NULL, model = "linear", datatype = "raw", N = NULL) {
  

  if (model == "linear") {
    
    if (datatype == "raw"){
      # sum auxiliary pop values by domain
    }
    if (datatype == "totals"){
      # check length
    }
    if (datatype == "means"){
      # check length
    }
    
    # perform this computation by domain based on domain_labels argument
    # preferably use some functional programming style (e.g map, apply) or figure out a matrix algebra way to do it
    w <- as.matrix(
      weights*ids + 
        (t(as.matrix(xpop_d) - xsample.dt_aoi %*% weights_aoi) %*%
           solve(xsample.dt %*% diag(weights) %*% xsample.d)) %*% 
        t(weights * xsample.d)
      )
    
    t <- w %*% y
    
    
    
    
  }
  
}
