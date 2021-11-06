# Playing and building the ratio estimator
# GREG estimators for numerator and denominator
library(tidyverse)
library(survey)
data(api)
devtools::load_all()

# Minimal Viable Product Work
# Start with PS

# Define terms
dat_s <- apisrs %>%
  drop_na(full, enroll)

xpop = data.frame(stype = apipop$stype)
xsample = data.frame(stype = dat_s$stype)
datatype = "raw"
y_num = dat_s$api.stu
y_den = dat_s$enroll
pi = NULL
N = NULL
var_est = TRUE
var_method = NULL
datatype= "raw"
model = "linear"
pi2 <- matrix(length(y_num)*(length(y_num) - 1)/nrow(xpop)/(nrow(xpop) - 1),
              nrow = length(y_num), ncol = length(y_num))
diag(pi2) <- length(y_num)/nrow(xpop)
head(pi2)


# -------------------------------------------
# Testing for bias in the estimators and variance estimators

store <- list()
store$B <- 1000
store$dat <- data.frame(ests = rep(NA, store$B),
                        var_ests_HTLin = rep(NA, store$B),
                        var_ests_gwtd = rep(NA, store$B))


#Remove missing values from the popdata
apipop  <- apipop %>%
  drop_na(stype, api.stu, enroll)
store$apipop <- apipop
store$n <- 100
store$pi2 <- matrix(store$n*(store$n - 1)/nrow(store$apipop )/(nrow(store$apipop ) - 1),
                    nrow = store$n, ncol = store$n)
diag(store$pi2) <- store$n/nrow(store$apipop )


for(i in 1:store$B){
  samp <- apipop %>%
    sample_n(store$n)
  
  store$dat$ests[i] <- gregRatio(xpop = data.frame(stype = apipop$stype),
                                      xsample = data.frame(stype = samp$stype),
                                      datatype = "raw",
                                      y_num = samp$api.stu,
                                      y_den = samp$enroll,
                                      pi = NULL,
                                      N = NULL,
                                      var_est = FALSE)$pop_ratio
  store$dat$var_ests_HTLin[i] <-
    gregRatio(xpop = data.frame(stype = apipop$stype),
              xsample = data.frame(stype = samp$stype),
              datatype = "raw",
              y_num = samp$api.stu,
              y_den = samp$enroll,
              pi = NULL,
              N = NULL,
              pi2 = store$pi2,
                   var_est = TRUE,
                   var_method = "LinHT")$pop_ratio_var
  
  store$dat$var_ests_gwtd[i] <-
    gregRatio(xpop = data.frame(stype = apipop$stype),
              xsample = data.frame(stype = samp$stype),
              datatype = "raw",
              y_num = samp$api.stu,
              y_den = samp$enroll,
              pi = NULL,
              N = NULL,
              pi2 = store$pi2,
              var_est = TRUE,
              var_method = "LinHTgWeighted")$pop_ratio_var
  
  # store$dat$var_ests_boot[i] <-
  #   postStratRatio(xpop = apipop$stype,
  #                  xsample = samp$stype,
  #                  datatype = "raw",
  #                  y_num = samp$api.stu,
  #                  y_den = samp$enroll,
  #                  pi = NULL,
  #                  N = NULL,
  #                  var_est = TRUE,
  #                  var_method = "bootstrapSRS")$pop_ratio_var
  
  print(i)
}

# save(store, file = "store_PS_ratio_2020_08_17.Rda")

# Results
store$dat %>%
  ggplot(aes(x = var_ests_HTLin)) +
  geom_histogram()

store$dat %>%
  ggplot(aes(x = ests)) +
  geom_histogram()


store$dat %>%
  summarize(var_emp = var(ests), 
            prb = (mean(ests) - sum(apipop$api.stu)/sum(apipop$enroll))/
              sum(apipop$api.stu)/sum(apipop$enroll)*100,
            var_exp = mean(var_ests_HTLin),
            var_exp_gwtd = mean(var_ests_gwtd),
            prb_var = (var_exp - var_emp)/var_emp*100,
            prb_var_gwtd = (var_exp_gwtd - var_emp)/var_emp*100,
            mse_emp = mean((ests - sum(apipop$api.stu)/sum(apipop$enroll))^2))

