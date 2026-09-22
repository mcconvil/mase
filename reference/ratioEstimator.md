# Compute a ratio estimator

Calculates a ratio estimator for a finite population mean/proportion or
total based on sample data collected from a complex sampling design and
auxiliary population data.

## Usage

``` r
ratioEstimator(
  y,
  xsample,
  xpop,
  datatype = "raw",
  pi = NULL,
  N = NULL,
  pi2 = NULL,
  var_est = FALSE,
  var_method = "LinHB",
  B = 1000,
  fpc = TRUE,
  messages = TRUE
)
```

## Arguments

- y:

  A numeric vector of the sampled response variable.

- xsample:

  A numeric vector of the sampled auxiliary variable.

- xpop:

  A numeric vector of population level auxiliary information. Must come
  in the form of raw data, population total or population mean.

- datatype:

  A string that specifies the form of population auxiliary data. The
  possible values are "raw", "total" or "mean". If datatype = "raw",
  then xpop must contain a numeric vector of the auxiliary variable for
  each unit in the population. If datatype = "total" or "mean", then
  contains either the population total or population mean for the
  auxiliary variable.

- pi:

  A numeric vector of inclusion probabilities for each sampled unit
  in y. If NULL, then simple random sampling without replacement is
  assumed.

- N:

  A numeric value of the population size. If NULL, it is estimated with
  the sum of the inverse of the pis.

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

- fpc:

  Default to TRUE, logical for whether or not the variance calculation
  should include a finite population correction when calculating the
  "LinHTSRS" or the "SRSbootstrap" variance estimator.

- messages:

  A logical indicating whether to output the messages internal to mase.
  Default is TRUE.

## Value

List of output containing:

\* pop_total: Estimate of population total.

\* pop_mean: Estimate of population mean.

\* pop_total_var: Estimated variance of population total estimate.

\* pop_mean_var: Estimated variance of population mean estimate.

## References

Cochran W~G (1977). *Sampling Techniques*, 3rd edition. John Wiley &
Sons, New York. Sarndal C~E, Swensson B, Wretman J (1992). *Model
Assisted Survey Sampling*. Springer-Verlag, New York.

## Examples

``` r
library(dplyr)
data(IdahoPop)
data(IdahoSamp)

xsample <- filter(IdahoSamp, COUNTYFIPS == "16055")
xpop <- filter(IdahoPop, COUNTYFIPS == "16055")

ratioEstimator(y = xsample$BA_TPA_ADJ,
               xsample = xsample$tcc,
               xpop = xpop$tcc,
               datatype = "means",
               N = xpop$npixels)
#> Assuming simple random sampling
#> $pop_total
#> [1] 36891107
#> 
#> $pop_mean
#> [1] 88.12489
#> 
```
