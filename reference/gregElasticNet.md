# Compute an elastic net regression estimator

Calculates a lasso, ridge or elastic net generalized regression
estimator for a finite population mean/proportion or total based on
sample data collected from a complex sampling design and auxiliary
population data.

## Usage

``` r
gregElasticNet(
  y,
  xsample,
  xpop,
  pi = NULL,
  alpha = 1,
  model = "linear",
  pi2 = NULL,
  var_est = FALSE,
  var_method = "LinHB",
  datatype = "raw",
  N = NULL,
  lambda = "lambda.min",
  B = 1000,
  cvfolds = 10,
  weights_method = "ridge",
  eta = 1e-04,
  fpc = TRUE,
  messages = TRUE,
  parallel = "multicore",
  ncpus = 2
)
```

## Arguments

- y:

  A numeric vector of the sampled response variable.

- xsample:

  A data frame of the auxiliary data in the sample.

- xpop:

  A data frame of population level auxiliary information. It must
  contain the same names as xsample. If datatype = "raw", must contain
  unit level data. If datatype = "totals" or "means", then contains one
  row of aggregated, population totals or means for the auxiliary data.
  Default is "raw".

- pi:

  A numeric vector of inclusion probabilities for each sampled unit
  in y. If NULL, then simple random sampling without replacement is
  assumed.

- alpha:

  A numeric value between 0 and 1 which signifies the mixing parameter
  for the lasso and ridge penalties in the elastic net. When alpha = 1,
  only a lasso penalty is used. When alpha = 0, only a ridge penalty is
  used. Default is alpha = 1.

- model:

  A string that specifies the regression model to utilize. Options are
  "linear" or "logistic".

- pi2:

  A square matrix of the joint inclusion probabilities. Needed for the
  "LinHT" variance estimator.

- var_est:

  A logical indicating whether or not to compute a variance estimator.
  Default is FALSE.

- var_method:

  The method to use when computing the variance estimator. Options are a
  Taylor linearized technique: "LinHB"= Hajek-Berger estimator, "LinHH"
  = Hansen-Hurwitz estimator, "LinHTSRS" = Horvitz-Thompson estimator
  under simple random sampling without replacement, and "LinHT" =
  Horvitz-Thompson estimator or a resampling technique: "bootstrapSRS" =
  bootstrap variance estimator under simple random sampling without
  replacement. The default is "LinHB".

- datatype:

  A string that specifies the form of population auxiliary data. The
  possible values are "raw", "totals" or "means" for whether the user is
  providing population data at the unit level, aggregated to totals, or
  aggregated to means. Default is "raw".

- N:

  A numeric value of the population size. If NULL, it is estimated with
  the sum of the inverse of the pis.

- lambda:

  A string specifying how to tune the lambda hyper-parameter. Only used
  if modelselect = TRUE and defaults to "lambda.min". The possible
  values are "lambda.min", which is the lambda value associated with the
  minimum cross validation error or "lambda.1se", which is the lambda
  value associated with a cross validation error that is one standard
  error away from the minimum, resulting in a smaller model.

- B:

  The number of bootstrap samples if computing the bootstrap variance
  estimator. Default is 1000.

- cvfolds:

  The number of folds for the cross-validation to tune lambda.

- weights_method:

  A string specifying which method to use to calculate survey weights.
  Currently, "ridge" is the only option. The "ridge" method uses a ridge
  regression approximation to calculate weights (see McConville et al
  (2017), section 3.2 for details). Support for "calibration" to come
  soon, which employs the model calibration method of Wu and Sitter
  (2001).

- eta:

  A small positive number. Defaults to 0.0001. See McConville et al
  (2017), section 3.2 for details.

- fpc:

  Default to TRUE, logical for whether or not the variance calculation
  should include a finite population correction when calculating the
  "LinHTSRS" or the "SRSbootstrap" variance estimator.

- messages:

  A logical indicating whether to output the messages internal to mase.
  Default is TRUE.

- parallel:

  The type of parallel processing to use, if any. Can take values of
  "multicore", "snow", or "no". Only used if var_method =
  "bootstrapSRS". See boot::boot() for more details on parallel
  processing types.

- ncpus:

  The number of parallel workers to use. Only used if var_method =
  "bootstrapSRS" and parallel = "snow" or "multicore".

## Value

A list of output containing:

\* pop_total: Estimate of population total.

\* coefficients: Elastic net coefficient estimates.

\* pop_mean: Estimate of the population mean (or proportion).

\* pop_total_var: Estimated variance of population total estimate.

\* pop_mean_var:Estimated variance of population mean estimate.

## References

McConville K~S, Breidt F~J, Lee T~C~M, Moisen G~G (2017).
“Model-Assisted Survey Regression Estimation with the Lasso.” *Journal
of Survey Statistics and Methodology*, **5**, 131-158.

## Examples

``` r
library(dplyr)
data(IdahoPop)
data(IdahoSamp)
xsample <- filter(IdahoSamp, COUNTYFIPS == "16055")
xpop <- filter(IdahoPop, COUNTYFIPS == "16055")
gregElasticNet(y = xsample$BA_TPA_ADJ,
               N = xpop$npixels,
               xsample = xsample[c("tcc", "elev", "ppt", "tmean")],
               xpop = xpop[c("tcc", "elev", "ppt", "tmean")],
               var_est = TRUE,
               var_method = "LinHB",
               datatype = "means",
               alpha = 0.5)
#> Assuming simple random sampling
#> $pop_total
#> [1] 39824537
#> 
#> $pop_mean
#> [1] 95.13222
#> 
#> $pop_total_var
#> [1] 6.36173e+12
#> 
#> $pop_mean_var
#> [1] 36.30188
#> 
#> $weights
#>   [1] 8110.6960 8127.9596 7941.2650 3921.2833 7408.5365 4513.9805 5072.4348
#>   [8] 3113.4666 2668.6354 1624.0108 3050.8954 5767.8385 3309.6885 4758.3397
#>  [15] 3515.1740 1072.4100 1341.8433 2575.2740 4324.6775 7854.7047 1764.1328
#>  [22] 2033.4285 5607.9361 4334.7522 6112.1528 1717.4420 2122.6873 3394.7072
#>  [29] 1673.3116 7415.5079 4197.2586 6329.4901 2163.0174 3216.2893  738.0288
#>  [36] 1196.9901  665.2475 2882.4305 7690.6379 5571.6106 6321.0567  883.0486
#>  [43] 3980.2542 4728.0694 6818.1577 2608.9368 3721.9647 2126.2434 1576.9905
#>  [50] 4366.7803 4596.6650 4106.1462 3914.2026 5396.3184 1239.4076 7226.7119
#>  [57] 1828.1824 6284.2792 1678.8441 6388.1890 2120.5595 4024.6626 6659.0979
#>  [64] 6361.2053 4558.0869 7180.3791 1872.7467 3622.3401 3478.5790 4049.6880
#>  [71] 5161.6505 5505.3941 1062.8080 1378.3264 2591.6583  636.4388 3864.2962
#>  [78] 5134.7709 1522.9424 5719.7012 5138.4440 4183.3825 7971.2081 3122.3591
#>  [85] 7943.2118 4054.2819 2670.7986 2655.2078 3870.2713 2620.7725 6439.1774
#>  [92] 6255.7970 3504.0819 3620.3362 9988.3242 4310.8087 5048.3191 8485.6856
#>  [99] 6652.4722 2892.1072
#> 
#> $coefficients
#>  (Intercept)          tcc         elev          ppt        tmean 
#> -25.76086919   0.43723829   0.03828208   0.07555344   0.00000000 
#> 
```
