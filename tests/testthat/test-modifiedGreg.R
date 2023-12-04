library(dplyr)
data(IdahoPop)
data(IdahoSamp)

res <- modifiedGreg(y = IdahoSamp$BA_TPA_ADJ,
                    xsample = IdahoSamp[c("tcc", "elev")],
                    xpop = IdahoPop[c("COUNTYFIPS","tcc", "elev", "npixels")] |> rename(N = npixels),
                    domains = IdahoSamp$COUNTYFIPS,
                    datatype = "means",
                    N = sum(IdahoPop$npixels),
                    var_est = TRUE)


test_that("modGreg.estimates", {
  
  expect_snapshot(res)
  
})


matrix(1:3, nrow = 3) - matrix(rep(1, 3*18), nrow = 3) %*% (1:18)
