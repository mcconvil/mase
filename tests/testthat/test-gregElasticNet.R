library(survey)
library(mase)
data(api)

set.seed(1)
out <- gregElasticNet(y = apisrs$api00, 
                      xsample = apisrs[c("col.grad", "awards", "snum", "dnum", "cnum", "pcttest", "meals", "sch.wide")], 
                      xpop = apipop[c("col.grad", "awards", "snum", "dnum", "cnum", "pcttest", "meals", "sch.wide")], 
                      pi = apisrs$pw^(-1),
                      model = "linear",
                      datatype = "raw",
                      var_est = TRUE,
                      var_method = "LinHB",
                      lambda = "lambda.min",
                      alpha = 1)

test_that("gregElasticNet.estimates", {
  
  expect_equal(out$pop_total, 4096571)
  expect_equal(out$pop_mean, 661.3774)
  expect_equal(out$pop_total_var, 803234237)
  expect_equal(out$pop_mean_var, 20.93629)
  
})
