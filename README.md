
# mase <a><img src='figs/mase_hex.png' align="right" height="280" /></a>

#### (Model Assisted Survey Estimation)

<!-- README.md is generated from README.Rmd. Please edit that file -->

<!-- badges: start -->
[![R-CMD-check](https://github.com/mcconvil/mase/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/mcconvil/mase/actions/workflows/R-CMD-check.yaml)
[![CRAN
status](https://www.r-pkg.org/badges/version/mase)](https://CRAN.R-project.org/package=mase)
<!-- badges: end -->

### Development Mode

`mase` is still under development. Please use at your own risk!

### Overview

`mase` contains a collection of model-assisted generalized regression
estimators for finite population estimation of a total or mean from a
single stage, unequal probability without replacement design. It also
contains several variance estimators.

The available estimators are currently:

- generalized regression: `greg()`
- hotvitz-thompson: `horvitzThompson()`
- post-stratification: `postStrat()`
- elastic net generalized regression: `gregElasticNet()`
- regression tree: `gregTree()`
- modified generalized regression: `modifiedGreg()`
- ratio estimator: `ratioEstimator()`
- estimate of a ratio estimator: `ratio()`

The available variance estimation techniques are:

- LinHB
- LinHH
- LinHTSRS
- LinHT
- bootstrapSRS

See `mase/inst/REFERENCES.bib` for sources related to each variance
estimator.

### Installation

Install the latest CRAN release with:

``` r
pak::pkg_install("mase")
```

You can also install the developmental version of `mase` from GitHub
with:

``` r
# install.packages("pak")
pak::pkg_install("mcconvil/mase")
```

### Usage

##### Horvitz-Thompson

Here’s an example of fitting the Horvitz-Thompson estimator using
Forestry data in Idaho. The data comes from the Forestry Inventory &
Analysis (FIA) program.

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

##### Linear Regression Estimator

We can also fit a linear regression estimator using that same data:

``` r
xsample <- select(samp, c(tcc, elev, ppt, tmean))

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
greg_est[c('pop_total','pop_mean', 'pop_total_var', 'pop_mean_var')]
#> $pop_total
#> [1] 39106643
#> 
#> $pop_mean
#> [1] 93.41733
#> 
#> $pop_total_var
#> [1] 6.328655e+12
#> 
#> $pop_mean_var
#> [1] 36.11314
```

But with this estimator we also get the weights

``` r
greg_est["weights"]
#> $weights
#>   [1] 8110.6960 8127.9599 7941.2651 3921.2834 7408.5365 4513.9805 5072.4347
#>   [8] 3113.4665 2668.6353 1624.0109 3050.8955 5767.8383 3309.6885 4758.3397
#>  [15] 3515.1741 1072.4099 1341.8432 2575.2742 4324.6776 7854.7045 1764.1326
#>  [22] 2033.4284 5607.9363 4334.7522 6112.1528 1717.4419 2122.6873 3394.7071
#>  [29] 1673.3117 7415.5078 4197.2586 6329.4902 2163.0174 3216.2894  738.0286
#>  [36] 1196.9899  665.2472 2882.4305 7690.6378 5571.6106 6321.0567  883.0485
#>  [43] 3980.2541 4728.0695 6818.1577 2608.9368 3721.9650 2126.2434 1576.9905
#>  [50] 4366.7802 4596.6651 4106.1462 3914.2027 5396.3184 1239.4076 7226.7119
#>  [57] 1828.1823 6284.2791 1678.8441 6388.1890 2120.5596 4024.6627 6659.0981
#>  [64] 6361.2053 4558.0869 7180.3791 1872.7464 3622.3400 3478.5788 4049.6881
#>  [71] 5161.6503 5505.3940 1062.8079 1378.3263 2591.6583  636.4387 3864.2963
#>  [78] 5134.7709 1522.9424 5719.7012 5138.4440 4183.3826 7971.2083 3122.3592
#>  [85] 7943.2118 4054.2819 2670.7987 2655.2078 3870.2713 2620.7724 6439.1774
#>  [92] 6255.7971 3504.0819 3620.3363 9988.3242 4310.8084 5048.3191 8485.6856
#>  [99] 6652.4721 2892.1071
```

and the coefficients for the model

``` r
greg_est["coefficients"]
#> $coefficients
#>   (Intercept)           tcc          elev           ppt         tmean 
#> -3.355552e+01  6.515276e-01  4.215046e-02  6.647252e-02  2.984714e-04
```

##### Variable Selection

All of the mase regression estimators can also perform variable
selection internally using the parameter `modelselect`

``` r
greg_select <- greg(y = samp$BA_TPA_ADJ,
                    N = pop$npixels,
                    xsample = xsample,
                    xpop = xpop,
                    modelselect = TRUE,
                    var_est = TRUE,
                    var_method = "LinHB",
                    datatype = "means")
```

And we can examine which predictors were chosen:

``` r
greg_select["coefficients"]
#> $coefficients
#>  (Intercept)          tcc         elev          ppt 
#> -33.24787647   0.65151379   0.04209371   0.06643125
```
