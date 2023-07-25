library(survey)
library(mase)
data(api) 

out <- ratio(y_num = apisrs$api.stu,
             y_den = apisrs$enroll,
             xsample = apisrs$stype,
             xpop = apipop$stype,
             pi = apisrs$pw^(-1), 
             estimator = "postStrat", 
             var_est = TRUE,
             var_method = "LinHB",
             datatype = "raw")

test_that("ratio.estimates", {
  
  expect_equal(out$ratio_est, 0.8258652)
  expect_equal(out$variance_est, 8.851525e-05)
  
})
