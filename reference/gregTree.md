# Compute a regression tree estimator

Calculates a regression tree estimator for a finite population
mean/proportion or total based on sample data collected from a complex
sampling design and auxiliary population data.

## Usage

``` r
gregTree(
  y,
  xsample,
  xpop,
  pi = NULL,
  pi2 = NULL,
  var_est = FALSE,
  var_method = "LinHB",
  B = 1000,
  pval = 0.05,
  perm_reps = 500,
  bin_size = NULL,
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

- B:

  The number of bootstrap samples if computing the bootstrap variance
  estimator. Default is 1000.

- pval:

  Designated p-value level to reject null hypothesis in permutation test
  used to fit the regression tree. Default value is 0.05.

- perm_reps:

  An integer specifying the number of permutations for each permutation
  test run to fit the regression tree. Default value is 500.

- bin_size:

  A integer specifying the minimum number of observations in each node.

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

\* weights: Survey weights produced by gregTree.

\* pop_total_var: Estimated variance of population total estimate.

\* pop_mean_var: Estimated variance of population mean estimate.

## References

McConville K~S, Toth D (2018). “Automated selection of post-strata using
a model-assisted regression tree estimator.” *Scandinavian Journal of
Statistics*.

## Examples

``` r
library(dplyr)
data(IdahoPop)
data(IdahoSamp)

xsample <- filter(IdahoSamp, COUNTYFIPS == "16055") 
xpop <- filter(IdahoSamp, COUNTYFIPS == "16055")

gregTree(y = xsample$BA_TPA_ADJ,
         xsample = xsample[c("tcc", "elev")],
         xpop = xpop[c("tcc", "elev")],
         var_est = TRUE)
#> Assuming simple random sampling
#> $pop_total
#> [1] 10722.31
#> 
#> $pop_mean
#> [1] 107.2231
#> 
#> $pop_total_var
#> [1] NaN
#> 
#> $pop_mean_var
#> [1] NaN
#> 
#> $weights
#>   [1] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#>  [38] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#>  [75] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#> 
#> $tree
#> 
#> RPMS Recursive Partitioning Equation 
#> y ~ tcc + elev
#> 
#> Estimating Equation 
#> y ~ 1
#> 
#> 
#> [1] "unequal probability of selection, sample design"
#> [1] "R-squared of model: 0.189980380370146"
#> 
#> ===================== Tree Model =================== 
#>  
#>    Splits         
#> sp   elev <= 934.5
#> 
#>     coefficients
#> node         1
#>    2  83.03977
#>    3 143.49798
#> 
#>  
#> 
```
