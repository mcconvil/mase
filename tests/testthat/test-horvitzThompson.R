library(mase)
library(survey)
data(api)


out <- horvitzThompson(y = apisrs$api00,
                       pi = apisrs$pw^(-1),
                       var_est = TRUE,
                       var_method = "LinHTSRS")


test_that("horvitzThompson.estimates", {
  
  expect_snapshot(out)
  
})
