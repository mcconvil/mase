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
  
  expect_snapshot(out)
  
})
