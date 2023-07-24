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

strat_expect <- data.frame(
  x = factor(c("No", "Yes")),
  strat_pop_total = c(1258153.57, 2827645.55),
  strat_pop_mean = c(620.69737, 678.58065)
)


test_that("postStrat.dim", {
  
  expect_equal(dim(out$strat_ests), c(2, 3))
  
})

test_that("postStrat.estimates", {
  
  expect_equal(out$pop_total, 4085799)
  expect_equal(out$pop_mean, 659.6382)
  expect_equal(out$pop_total_var, 3135191742)
  expect_equal(out$pop_mean_var, 81.71875)
  
  expect_equal(as.data.frame(out$strat_ests), strat_expect)
  
})
