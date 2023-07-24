library(mase)
library(survey)
data(api)


out <- horvitzThompson(y = apisrs$api00,
                       pi = apisrs$pw^(-1),
                       var_est = TRUE,
                       var_method = "LinHTSRS")


test_that("horvitzThompson.estimates", {
  
  expect_equal(out$pop_total, 4066887)
  expect_equal(out$pop_mean, 656.585)
  expect_equal(out$pop_total_var, 3282462447)
  expect_equal(out$pop_mean_var, 85.55736)
  
})
