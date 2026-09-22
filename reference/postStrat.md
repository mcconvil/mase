# Compute a post-stratified estimator

Calculates a post-stratified estimator for a finite population
mean/proportion or total based on sample data collected from a complex
sampling design and a single, categorical auxiliary population variable.

## Usage

``` r
postStrat(
  y,
  xsample,
  xpop,
  pi = NULL,
  N = NULL,
  var_est = FALSE,
  var_method = "LinHB",
  pi2 = NULL,
  datatype = "raw",
  B = 1000,
  fpc = TRUE,
  messages = TRUE
)
```

## Arguments

- y:

  A numeric vector of the sampled response variable.

- xsample:

  A vector containing the post-stratum for each sampled unit.

- xpop:

  A vector or data frame, depending on datatype. If datatype = "raw",
  then a vector containing the post-stratum for each population unit. If
  datatype = "totals" or "means", then a data frame, where the first
  column lists the possible post-strata and the second column contains
  the population total or proportion in each post-stratum.

- pi:

  A numeric vector of inclusion probabilities for each sampled unit
  in y. If NULL, then simple random sampling without replacement is
  assumed.

- N:

  A numeric value of the population size. If NULL, it is estimated with
  the sum of the inverse of the pis.

- var_est:

  Default to FALSE, logical for whether or not to compute estimate of
  variance

- var_method:

  The method to use when computing the variance estimator. Options are a
  Taylor linearized technique: "LinHB"= Hajek-Berger estimator, "LinHH"
  = Hansen-Hurwitz estimator, "LinHTSRS" = Horvitz-Thompson estimator
  under simple random sampling without replacement, and "LinHT" =
  Horvitz-Thompson estimator or a resampling technique: "bootstrapSRS" =
  bootstrap variance estimator under simple random sampling without
  replacement, "SRSunconditional" = simple random sampling variance
  estimator which accounts for random strata.

- pi2:

  A square matrix of the joint inclusion probabilities. Needed for the
  "LinHT" variance estimator.

- datatype:

  Default to "raw", takes values "raw", "totals" or "means" for whether
  the user is providing the raw population stratum memberships, the
  population totals of each stratum, or the population proportions of
  each stratum.

- B:

  The number of bootstrap samples if computing the bootstrap variance
  estimator. Default is 1000.

- fpc:

  Default to TRUE, logical for whether or not the variance calculation
  should include a finite population correction when calculating the
  "LinHTSRS", "SRSunconditional", or the "SRSbootstrap" variance
  estimator.

- messages:

  A logical indicating whether to output the messages internal to mase.
  Default is TRUE.

## Value

A list of output containing:

\* pop_total: Estimate of population total.

\* pop_mean: Estimate of the population mean (or proportion).

\* pop_total_var: Estimated variance of population total estimate.

\* pop_mean_var: Estimated variance of population mean estimate.

\* strat_ests: Table of total and mean estimates for each strata.

\* weights: Survey weights produced by PS.

## References

Cochran W~G (1977). *Sampling Techniques*, 3rd edition. John Wiley &
Sons, New York.

Sarndal C~E, Swensson B, Wretman J (1992). *Model Assisted Survey
Sampling*. Springer-Verlag, New York.

## Examples

``` r
library(tidyr)
library(dplyr)

data(IdahoPop)
data(IdahoSamp)

xsample <- filter(IdahoSamp, COUNTYFIPS == "16055")
xpop <- filter(IdahoPop, COUNTYFIPS == "16055")

pop <- xpop[c("tnt.1", "tnt.2")] |> 
    pivot_longer(everything(), names_to = "tnt", values_to = "prop") |>
    mutate(tnt = as.numeric(gsub("\\D", "", tnt)))
    
postStrat(y = xsample$BA_TPA_ADJ,
          N = xpop$npixels,
          xsample = xsample$tnt,
          xpop = pop,
          datatype = "means",
          var_est = TRUE,
          var_method = "SRSunconditional")
#> Assuming simple random sampling
#> $pop_total
#> [1] 40148687
#> 
#> $pop_mean
#> [1] 95.90655
#> 
#> $pop_total_var
#>              [,1]
#> [1,] 6.724438e+12
#> 
#> $pop_mean_var
#>          [,1]
#> [1,] 38.37159
#> 
#> $strat_ests
#> # A tibble: 2 × 3
#>   x     strat_pop_total strat_pop_mean
#>   <fct>           <dbl>          <dbl>
#> 1 1           35004878.          114. 
#> 2 2            5143809.           46.1
#> 
#> $weights
#>   [1] 0.0006368499 0.0001946594 0.0006368499 0.0001946594 0.0001946594
#>   [6] 0.0001946594 0.0001946594 0.0001946594 0.0001946594 0.0001946594
#>  [11] 0.0001946594 0.0001946594 0.0001946594 0.0001946594 0.0001946594
#>  [16] 0.0001946594 0.0001946594 0.0001946594 0.0001946594 0.0006368499
#>  [21] 0.0001946594 0.0001946594 0.0001946594 0.0001946594 0.0001946594
#>  [26] 0.0001946594 0.0001946594 0.0001946594 0.0001946594 0.0006368499
#>  [31] 0.0001946594 0.0001946594 0.0001946594 0.0001946594 0.0001946594
#>  [36] 0.0001946594 0.0001946594 0.0001946594 0.0001946594 0.0001946594
#>  [41] 0.0001946594 0.0001946594 0.0001946594 0.0001946594 0.0001946594
#>  [46] 0.0001946594 0.0001946594 0.0001946594 0.0001946594 0.0001946594
#>  [51] 0.0001946594 0.0001946594 0.0001946594 0.0001946594 0.0001946594
#>  [56] 0.0006368499 0.0001946594 0.0001946594 0.0001946594 0.0001946594
#>  [61] 0.0001946594 0.0001946594 0.0006368499 0.0001946594 0.0006368499
#>  [66] 0.0001946594 0.0001946594 0.0001946594 0.0001946594 0.0001946594
#>  [71] 0.0006368499 0.0001946594 0.0001946594 0.0001946594 0.0001946594
#>  [76] 0.0001946594 0.0001946594 0.0001946594 0.0001946594 0.0006368499
#>  [81] 0.0001946594 0.0001946594 0.0001946594 0.0001946594 0.0001946594
#>  [86] 0.0001946594 0.0001946594 0.0001946594 0.0001946594 0.0001946594
#>  [91] 0.0001946594 0.0001946594 0.0001946594 0.0001946594 0.0006368499
#>  [96] 0.0001946594 0.0001946594 0.0001946594 0.0001946594 0.0001946594
#> 
          
```
