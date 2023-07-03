library(tidyverse)
library(mase)

# county level population auxiliary data means
unitzonal <- read_rds("/Users/joshuayamamoto/Downloads/Oregon/unitzonal.rds")

# county level acreage
unitarea <- read_rds("/Users/joshuayamamoto/Downloads/Oregon/unitarea.rds")

# plot level forest attribute data
tsumdatp <- read_rds("/Users/joshuayamamoto/Downloads/Oregon/tsumdatp.rds") 

# plot level forest attribute data
tsumdatc <- read_rds("/Users/joshuayamamoto/Downloads/Oregon/tsumdatc.rds")

# strata weights
stratalut <- read_rds("/Users/joshuayamamoto/Downloads/Oregon/stratalut.rds")

# sample auxiliary data
pltassign <- read_rds("/Users/joshuayamamoto/Downloads/Oregon/pltassgn.rds") 



# using greg
sampx <- tsumdatp %>% 
  left_join(pltassign, by = "CN") %>% 
  select(tcc16, elev, COUNTYCD, BA_TPA_live_ADJ) %>% 
  drop_na()

sampx_aoi <- sampx %>% 
  filter(COUNTYCD == 1) %>% 
  select(BA_TPA_live_ADJ, tcc16, elev)

popx <- unitzonal %>% 
  select(tcc16, elev, COUNTYFIPS, npixels) %>% 
  drop_na()

popx_aoi <- popx %>% 
  filter(COUNTYFIPS == "41001") %>% 
  select(tcc16, elev, npixels)


greg(y = sampx_aoi$BA_TPA_live_ADJ, xsample = sampx_aoi[c("tcc16", "elev")],
     xpop = popx_aoi[c("tcc16", "elev")], datatype = "means", N = popx_aoi$npixels)$pop_mean




# rough bones of modified greg
# first test by doing manual fitting

mod_full <- lm(BA_TPA_live_ADJ ~ tcc16 + elev, data = sampx)

for_aoi <- function(model, popdata_aoi, sampdata) {
  
  est <- (1/nrow(sampdata))*sum(sampdata$BA_TPA_live_ADJ - predict(model, sampdata)) +
    (1/nrow(popdata_aoi)*sum(predict(model, popdata_aoi)))
  
  return(est)
  
}

for_aoi(mod_full, popx_aoi, sampx_aoi)


## modified greg matrix style --------------------------------------------------

n <- dim(sampx)[1]
N <- sum(popx$npixels)
pi <- rep(n/N, n)
weights <- as.vector(pi^(-1))
ids <- as.integer(sampx$COUNTYCD == 1)

# for whole sample
xsample.d <- model.matrix(~., data = data.frame(sampx[c("tcc16", "elev")]))
xsample <- data.frame(xsample.d[,-1, drop = FALSE])
xsample.dt <- t(xsample.d) 


# for specific aoi
xsample.d_aoi <- model.matrix(~., data = data.frame(sampx_aoi[c("tcc16", "elev")]))
xsample_aoi <- data.frame(xsample.d_aoi[,-1, drop = FALSE])
xsample.dt_aoi <- t(xsample.d_aoi) 
weights_aoi <- rep(weights[1], nrow(sampx_aoi))

# for specific aoi
#xpop <- data.frame(model.matrix(~.-1, data = data.frame(popx_aoi)))
xpop <- popx_aoi[1:2]*popx_aoi$npixels
#xpop_d <- model.matrix(~., data = xpop)
xpop_d <- unlist(c(popx_aoi$npixels,xpop))





w <- as.matrix(weights*ids + 
            t(as.matrix(xpop_d) - xsample.dt_aoi %*% weights_aoi) %*% 
            solve(xsample.dt %*% diag(weights) %*% xsample.d) %*% 
            (xsample.dt %*% diag(weights))
          )

t <- w %*% sampx$BA_TPA_live_ADJ

t/popx_aoi$npixels

































