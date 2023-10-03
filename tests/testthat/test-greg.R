library(mase)
library(survey)
data(api)

out_lin <- greg(y = apisrs$api00,
                xsample = apisrs[c("col.grad", "awards")], 
                xpop = apipop[c("col.grad", "awards")],
                pi = apisrs$pw^(-1), 
                model = "linear",
                var_est = TRUE,
                datatype = "raw",
                var_method = "LinHB")

out_log <- greg(y = as.integer(apisrs$awards) - 1,
                xsample = apisrs[c("col.grad", "meals")],
                xpop = apipop[c("col.grad", "meals")],
                pi = apisrs$pw^(-1),
                model = "logistic",
                var_est = TRUE,
                datatype = "raw", 
                var_method = "LinHTSRS")


test_that("greg.estimates", {
  
  expect_snapshot(out_lin)
  expect_snapshot(out_log)
  
 
})

test_that("greg.weights", {
  
  expect_equal(sum(out_lin$weights), nrow(apipop))
  
})
