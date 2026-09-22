# Compute a ratio of two estimators

Compute a ratio of two estimators

## Usage

``` r
ratio(
  y_num,
  y_den,
  xsample,
  xpop,
  pi = NULL,
  pi2 = NULL,
  N = NULL,
  estimator = NULL,
  var_est = F,
  var_method = "LinHTSRS",
  datatype = "raw",
  fpc = TRUE,
  messages = TRUE,
  ...
)
```

## Arguments

- y_num:

  A vector containing the response value for each sampled unit in the
  numerator

- y_den:

  A vector containing the response value for each sampled unit in the
  denominator

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

- N:

  A numeric value of the population size. If NULL, it is estimated with
  the sum of the inverse of the pis.

- estimator:

  A string containing the name of the estimators of which you are taking
  a ratio of. The names follow the same format as the functions
  independently do in mase. Options are "horvitzThompson", "postStrat",
  and "greg".

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

  Default to "raw", takes values "raw", "totals" or "means" for whether
  the user is providing the raw population stratum memberships, the
  population totals of each stratum, or the population proportions of
  each stratum.

- fpc:

  Default to TRUE, logical for whether or not the variance calculation
  should include a finite population correction when calculating the
  "LinHTSRS" or the "SRSbootstrap" variance estimator.

- messages:

  A logical indicating whether to output the messages internal to mase.
  Default is TRUE.

- ...:

  Any additional arguments that can be passed to mase::horvitzThompson,
  mase::greg, and mase::postStrat

## Value

A list of output containing:

\* ratio_est: Estimate of the ratio of the population totals/means of
the two estimators.

\* ratio_var_est: Estimate of the variance of the ratio of two
estimators.

## References

Cochran W~G (1977). *Sampling Techniques*, 3rd edition. John Wiley &
Sons, New York. Sarndal C~E, Swensson B, Wretman J (1992). *Model
Assisted Survey Sampling*. Springer-Verlag, New York.

## Examples

``` r
library(survey)
#> Loading required package: grid
#> Loading required package: Matrix
#> 
#> Attaching package: ‘Matrix’
#> The following objects are masked from ‘package:tidyr’:
#> 
#>     expand, pack, unpack
#> Loading required package: survival
#> 
#> Attaching package: ‘survey’
#> The following object is masked from ‘package:graphics’:
#> 
#>     dotchart
data(api) 
ratio(y_num = apisrs$api.stu, y_den = apisrs$enroll, xsample = apisrs$stype,
xpop = apipop$stype, pi = apisrs$pw^(-1), estimator = "postStrat",
var_est = TRUE, var_method = "LinHB", datatype = "raw")
#> $ratio_est
#> [1] 0.8258652
#> 
#> $ratio_var_est
#> [1] 8.851525e-05
#> 
```
