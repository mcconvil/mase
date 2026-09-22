# Compute a generalized regression estimator

Calculates a generalized regression estimator for a finite population
mean/proportion or total based on sample data collected from a complex
sampling design and auxiliary population data.

## Usage

``` r
greg(
  y,
  xsample,
  xpop,
  pi = NULL,
  model = "linear",
  pi2 = NULL,
  var_est = FALSE,
  var_method = "LinHB",
  datatype = "raw",
  N = NULL,
  modelselect = FALSE,
  lambda = "lambda.min",
  B = 1000,
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

- modelselect:

  A logical for whether or not to run lasso regression first and then
  fit the model using only the predictors with non-zero lasso
  coefficients. Default is FALSE.

- lambda:

  A string specifying how to tune the lasso hyper-parameter. Only used
  if modelselect = TRUE and defaults to "lambda.min". The possible
  values are "lambda.min", which is the lambda value associated with the
  minimum cross validation error or "lambda.1se", which is the lambda
  value associated with a cross validation error that is one standard
  error away from the minimum, resulting in a smaller model.

- B:

  The number of bootstrap samples if computing the bootstrap variance
  estimator. Default is 1000.

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

\* pop_mean: Estimate of the population mean (or proportion).

\* weights: Survey weights produced by GREG (linear model only).

\* pop_total_var: Estimated variance of population total estimate.

\* pop_mean_var: Estimated variance of population mean estimate.

## References

Cassel C~M, Sarndal C~E, Wretman J~H (1976). “Some results on
generalized difference estimation and generalized regression estimation
for finite populations.” *Biometrika*, **63**, 615–620.

Sarndal C~E, Swensson B, Wretman J (1992). *Model Assisted Survey
Sampling*. Springer-Verlag, New York.

## Examples

``` r
library(dplyr)
#> 
#> Attaching package: ‘dplyr’
#> The following objects are masked from ‘package:stats’:
#> 
#>     filter, lag
#> The following objects are masked from ‘package:base’:
#> 
#>     intersect, setdiff, setequal, union
data(IdahoPop)
data(IdahoSamp)

xsample <- filter(IdahoSamp, COUNTYFIPS == "16055")
xpop <- filter(IdahoPop, COUNTYFIPS == "16055")

greg(y = xsample$BA_TPA_ADJ,
     N = xpop$npixels,
     xsample = xsample[c("tcc", "elev")],
     xpop = xpop[c("tcc", "elev")],
     var_est = TRUE,
     var_method = "LinHB",
     datatype = "means")
#> Assuming simple random sampling
#> $pop_total
#> [1] 39521557
#> 
#> $pop_mean
#> [1] 94.40847
#> 
#> $pop_total_var
#> [1] 6.478447e+12
#> 
#> $pop_mean_var
#> [1] 36.9679
#> 
#> $weights
#>   [1] 8060.960 6489.037 6929.206 3145.983 8041.396 3737.945 5687.480 3577.979
#>   [9] 3020.678 1690.956 2427.990 6469.919 3123.403 4879.079 2639.280 2624.305
#>  [17] 1954.575 1186.989 3207.814 8886.670 3196.904 1876.082 4186.524 3406.140
#>  [25] 5158.757 2823.357 2449.285 4110.968 1848.589 8638.148 3758.716 4917.963
#>  [33] 2738.744 2695.673 1852.252 2086.525 2920.164 2799.124 7854.180 4707.118
#>  [41] 5730.753 1262.825 4115.680 4101.833 7274.657 2607.800 2055.045 2659.727
#>  [49] 1947.451 4639.413 3681.351 4315.134 3403.483 4823.892 1674.373 7991.121
#>  [57] 2524.517 6087.918 1799.399 5632.293 2457.738 3312.270 6021.943 5538.423
#>  [65] 3995.847 6102.368 3174.927 4497.681 4476.307 2101.377 6093.312 5503.605
#>  [73] 1478.582 2305.018 2687.543 2518.520 3852.988 5259.550 2145.174 5314.536
#>  [81] 5218.411 3132.459 7656.256 3177.261 8595.077 3461.808 2174.599 2428.715
#>  [89] 3508.262 3680.625 5900.982 5897.923 3362.667 2789.342 9305.020 5559.150
#>  [97] 5000.197 8538.360 7337.170 2925.479
#> 
#> $coefficients
#>  (Intercept)          tcc         elev 
#> -30.26159045   0.79210774   0.09679431 
#> 
```
