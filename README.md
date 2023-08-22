
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Development Mode

`mase` is still under development. Please use at your own risk!

# mase

`mase` contains a collection of model-assisted generalized regression
estimators (the post-stratification estimator, the ratio estimator, the
linear and logistic regression estimator, the modified linear and
logistic regression estimator, the elastic net regression estimator, and
the regression tree estimator) for finite population estimation of a
total or mean from a single stage, unequal probability without
replacement design. It also contains the Horvitz-Thompson estimator and
several variance estimators.

## Installation

You can install `mase` from github with:

``` r
# install.packages("devtools")
devtools::install_github("mcconvil/mase")
```

## Example

Here’s an example of fitting the Horvitz-Thompson estimator using
Forestry data in Idaho. The data comes from the Forestry Inventory &
Analysis department (FIA).

``` r
library(mase)
library(dplyr)

data(IdahoSamp)
data(IdahoPop)

samp <- filter(IdahoSamp, COUNTYFIPS == 16055) 
pop <- filter(IdahoPop, COUNTYFIPS == 16055) 

horvitzThompson(y = samp$BA_TPA_ADJ,
                N = pop$npixels,
                var_est = TRUE,
                var_method = "LinHTSRS")
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
```

We can also fit a linear regression estimator using that same data:

``` r
xsample <- select(samp, tcc, elev) %>%
  as.data.frame()

xpop <- select(pop, names(xsample))

greg_est <- greg(y = samp$BA_TPA_ADJ,
     N = pop$npixels,
     xsample = xsample,
     xpop = xpop,
     var_est = TRUE,
     var_method = "LinHB",
     datatype = "means")
```

We still get the population total and mean estimates along with their
variance estimates:

``` r
greg_est[1:4]
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
```

But with this estimator we also get the weights

``` r
greg_est[5] |> head(7)
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
```

and the coefficients for the model

``` r
greg_est[6]
#> $coefficients
#>                     [,1]
#> (Intercept) -30.26159045
#> tcc           0.79210774
#> elev          0.09679431
```
