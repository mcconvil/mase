library(mase)
library(survey)
data(api)


out <- postStrat(y = apisrs$api00,
          xsample = apisrs$awards, 
          xpop = data.frame(table(apipop$awards)),
          pi = apisrs$pw^(-1),
          datatype = "totals",
          var_est = TRUE,
          var_method = "LinHB")

test_that("postStrat.dim", {
  
  expect_equal(dim(out$strat_ests), c(2, 3))
  
})

test_that("postStrat.estimates", {
  
  expect_snapshot(out)

})
