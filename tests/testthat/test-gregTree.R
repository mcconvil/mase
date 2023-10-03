library(survey)
library(mase)
data(api)
 
out <- gregTree(y = apisrs$api00, 
         xsample = apisrs[c("col.grad", "awards", "snum", "dnum", "cnum", "pcttest", "meals", "sch.wide")], 
         xpop = apipop[c("col.grad", "awards", "snum", "dnum", "cnum", "pcttest", "meals", "sch.wide")],
         var_est = T,
         var_method = "LinHB")


test_that("gregTree.estimates", {
  
  expect_snapshot(out)
  
})

