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

modified_greg <- function(y, xsample, xpop, domain, domain_labels, pi = NULL, model = "linear", datatype = "raw", N = NULL) {
  
  # need xpop by domain var
  
  if (model == "linear") {
    
    if (datatype == "raw"){
      # sum by domain var
      xpop <- data.frame(model.matrix(~.-1, data = data.frame(xpop)))
      xpop <- xpop[intersect(names(xpop), names(xsample))]
      xpop_d <- cbind(as.data.frame(model.matrix(~., data = xpop)), domain = xpop_domain)
      aggregate( . ~ domain, xpop_d, FUN = sum)
      xpop_d <- apply(xpop_d, 2, sum)
      
    }
    if (datatype == "totals"){
      # check length
      xpop_d <- unlist(c(N, xpop[names(xsample)]))
    }
    if (datatype == "means"){
      # check length
      xpop_d <- unlist(c(N,xpop[names(xsample)]*N))
    }
    
  }
  
}
