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
  
  expect_snapshot(out)
  
})
