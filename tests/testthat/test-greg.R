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
  
  expect_equal(out_lin$pop_total, 4116418)
  expect_equal(out_lin$pop_mean, 664.5815)
  expect_equal(out_lin$pop_total_var, 1547721897)
  expect_equal(out_lin$pop_mean_var, 40.34136)
  
  expect_equal(out_log$pop_total, 3834.216)
  expect_equal(out_log$pop_mean, 0.619021)
  expect_equal(as.numeric(out_log$pop_total_var), 43380.63)
  expect_equal(as.numeric(out_log$pop_mean_var),  0.001130716)
  
 
})

test_that("greg.weights", {
  
  expect_equal(sum(out_lin$weights), nrow(apipop))
  
})
