# Playing and building the ratio estimator
# Post-stratified estimators for numerator and denominator
library(tidyverse)
library(survey)
data(api)

# Want total number of tested/total number of students enrolled
# = prop of students tested
# Why don't we just take number tested/number enrolled and 
# then find mean? -- different stat
dat_s <- apisrs %>%
  drop_na(full, enroll)


# Test it
thing1 <- postStratRatio(xpop = apipop$stype,
               xsample = dat_s$stype,
               datatype = "raw",
               y_num = dat_s$api.stu,
               y_den = dat_s$enroll,
               pi = NULL,
               N = NULL,
               var_est = TRUE,
               var_method = "SRSunconditional",
               B = 1000)
# thing2 <- postStratRatio(xpop = apipop$stype,
#                xsample = dat_s$stype,
#                datatype = "raw",
#                y_num = dat_s$api.stu,
#                y_den = dat_s$enroll,
#                pi = NULL,
#                N = NULL,
#                var_est = TRUE,
#                var_method = "bootstrapSRS",
#                B = 1000)
thing1
thing2



# Estimate ratio of totals
# Total students who took the exam/Total number enrolled
# Proportion of students who took exam
sum(dat_s$api.stu)/sum(dat_s$enroll)
# Not the same: Average proportion per school
mean(dat_s$api.stu/dat_s$enroll)

# -------------------------------------------

# Minimal Viable Product Work
# Start with PS

# Define terms
xpop = apipop$stype
xsample = dat_s$stype
datatype = "raw"
y_num = dat_s$api.stu
y_den = dat_s$enroll
pi = NULL
N = NULL
var_est = TRUE
var_method = "SRSunconditional" 
datatype= "raw"
B = 1000


y_num_ps <- postStrat(y = dat_s$api.stu, xsample = dat_s$stype,
          xpop = apipop$stype, var_est =  TRUE, 
          var_method = "SRSunconditional")
y_den_ps <- postStrat(y = dat_s$enroll, xsample = dat_s$stype,
          xpop = apipop$stype, var_est =  TRUE, 
          var_method = "SRSunconditional")

#Ratio: Point Estimate
Rat <- y_num_ps$pop_total/y_den_ps$pop_total

# Var est
# Use 4.17 in Green Book
# COV 
# Need:
# N = Pop size
# n = sample size
# N_h = proportion of pop in stratum h
# n_h = proportion of sample in stratum h
# Estimates for each stratum
# y for both y_num and y_den



# Note: Need to check that length of 
# y_num sample and y_den sample are equal
#Compute estimator

if (is.null(pi)) {
  pi <- rep(length(y)/N, length(y))
}

if (datatype=="raw"){
  xpop_tab <- data.frame(table(xpop))
  colnames(xpop_tab) <- c("x","N_h")
}
#Make sure that x_pop_tab$x is a factor
if(!is.factor(xpop_tab$x)){
  xpop_tab$x <- as.factor(xpop_tab$x)
}

xsample = dat_s$stype 
#Make xsample have same column name as xpop_tab
xsample <- data.frame(xsample)
colnames(xsample) <- "x"

# Compute components
xsample_pi_num_den <- data.frame(xsample, pi, y_num, y_den)
tab <- xsample_pi_num_den %>%
  group_by(x) %>%
  summarize(poptotal_h_num = y_num%*%pi^(-1),
            poptotal_h_den = y_den%*%pi^(-1),
            N_h_hats = sum(pi^(-1))) %>%
  inner_join(xpop_tab, by=c("x")) %>%
  mutate(ps_h = N_h/N_h_hats, 
         strat_pop_total_num = ps_h*poptotal_h_num,
         strat_pop_mean_num = strat_pop_total_num/N_h,
         strat_pop_total_den = ps_h*poptotal_h_den,
         strat_pop_mean_den = strat_pop_total_den/N_h)

cov_piece <- xsample_pi_num_den %>%
  left_join(tab, by=c("x")) %>%
  group_by(x) %>%
  summarize(n_h = n(), strat_pop_mean_num = first(strat_pop_mean_num),
            strat_pop_mean_den = first(strat_pop_mean_den),
            Q = (sum(y_num*y_den) - n_h*strat_pop_mean_num*strat_pop_mean_den)/
              n_h/(n_h - 1))
  
cov_est <- N^2*(1-n/N)/n*t(tab$N_h/N)%*%cov_piece$Q +
  N^2*(1-n/N)/n^2*t(1-tab$N_h/N)%*%cov_piece$Q


var_est <- y_den_ps$pop_total^(-2)*(y_num_ps$pop_total_var + 
               Rat^2*y_den_ps$pop_total_var - 2*Rat*cov_est)
se_est <- sqrt(var_est)
               

# -------------------------------------------
# Testing for bias in the estimators and variance estimators

store <- list()
store$B <- 1000
store$dat <- data.frame(ests = rep(NA, store$B),
                        var_ests_uncond = rep(NA, store$B))
store$n <- 100


#Remove missing values from the popdata
apipop  <- apipop %>%
  drop_na(stype, api.stu, enroll)
store$apipop <- apipop

for(i in 1:B){
  samp <- apipop %>%
    sample_n(store$n)
  
  store$dat$ests[i] <- postStratRatio(xpop = apipop$stype,
                                      xsample = samp$stype,
                                      datatype = "raw",
                                      y_num = samp$api.stu,
                                      y_den = samp$enroll,
                                      pi = NULL,
                                      N = NULL,
                                      var_est = FALSE)$pop_ratio
  store$dat$var_ests_uncond[i] <-
    postStratRatio(xpop = apipop$stype,
                   xsample = samp$stype,
                   datatype = "raw",
                   y_num = samp$api.stu,
                   y_den = samp$enroll,
                   pi = NULL,
                   N = NULL,
                   var_est = TRUE,
                   var_method = "SRSunconditional")$pop_ratio_var
  
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
  ggplot(aes(x = var_ests_uncond)) +
  geom_histogram()

store$dat %>%
  ggplot(aes(x = ests)) +
  geom_histogram()


store$dat %>%
  summarize(var_emp = var(ests), 
            prb = (mean(ests) - sum(apipop$api.stu)/sum(apipop$enroll))/
              sum(apipop$api.stu)/sum(apipop$enroll)*100,
            var_exp = mean(var_ests_uncond),
            prb_var = (var_exp - var_emp)/var_emp*100,
            mse_emp = mean((ests - sum(apipop$api.stu)/sum(apipop$enroll))^2))


# -------------------------------------------

# FIESTA Comparison
library(tidyverse)
tdomdat <- readRDS("~/Research/mase/sandbox/tdomdattot_volcf_ndead_dlive.rds")
#unit_totest <- readRDS("~/Research/mase/sandbox/ratio_unit_totest.rds")
unitlut <- readRDS("~/Research/mase/sandbox/unitlut.rds")
stratalut <- readRDS("~/Research/mase/sandbox/stratalut.rds")
strata_info <- stratalut %>%
  select(STRATUMCD, strwt)

devtools::load_all()

system.time({
thing1 <- postStratRatio(xpop = strata_info,
                         xsample = tdomdat$STRATUMCD,
                         datatype = "means",
                         y_num = tdomdat$VOLCFNET_TPA_ADJ,
                         y_den = tdomdat$VOLCFNET_TPA_ADJ.d,
                         pi = NULL,
                         N = unitlut$npixels,
                         var_est = TRUE,
                         var_method = "SRSunconditional",
                         fpc = FALSE)
thing1
})

thing2 <- postStratRatio(xpop = strata_info,
                         xsample = tdomdat$STRATUMCD,
                         datatype = "means",
                         y_num = tdomdat$BA_TPA_ADJ,
                         y_den = tdomdat$BA_TPA_ADJ.d,
                         pi = NULL,
                         N = unitlut$npixels,
                         var_est = TRUE,
                         var_method = "SRSunconditional",
                         fpc = FALSE)
