library(survey)
library(mase)
data(api)
 
out <- gregTree(y = apisrs$api00, 
         xsample = apisrs[c("col.grad", "awards", "snum", "dnum", "cnum", "pcttest", "meals", "sch.wide")], 
         xpop = apipop[c("col.grad", "awards", "snum", "dnum", "cnum", "pcttest", "meals", "sch.wide")],
         var_est = T,
         var_method = "LinHB")


test_that("gregTree.estimates", {
  
  expect_equal(out$pop_total, 4118359)
  expect_equal(out$pop_mean, 664.8949)
  expect_equal(out$pop_total_var, 1125912522)
  expect_equal(out$pop_mean_var, 29.3469)
  
})

