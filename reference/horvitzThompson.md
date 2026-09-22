# Compute the Horvitz-Thompson Estimator

Calculate the Horvitz-Thompson Estimator for a finite population
mean/proportion or total based on sample data collected from a complex
sampling design.

## Usage

``` r
horvitzThompson(
  y,
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

Horvitz DG, Thompson DJ (1952). “A generalization of sampling without
replacement from a finite universe.” *Journal of the American
Statistical Association*, **47**, 663-685.

## Examples

``` r
library(dplyr)
data(IdahoSamp)
data(IdahoPop)
xsample <- filter(IdahoSamp, COUNTYFIPS == "16055")
xpop <- filter(IdahoPop, COUNTYFIPS == "16055") 

horvitzThompson(y = xsample$BA_TPA_ADJ,
                N = xpop$npixels,
                var_est = TRUE,
                var_method = "LinHTSRS")
#> Assuming simple random sampling
#> $pop_total
#> [1] 44886038
#> 
#> $pop_mean
#> [1] 107.2231
#> 
#> $pop_total_var
#> [1] 8.171847e+12
#> 
#> $pop_mean_var
#> [1] 46.63093
#> 
```
