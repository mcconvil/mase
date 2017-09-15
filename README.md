
<!-- README.md is generated from README.Rmd. Please edit that file -->
Development Mode
================

mase is still under development. Please use at your own risk!

mase
====

mase contains a collection of model-assisted generalized regression estimators (the post-stratification estimator, the ratio estimator, the linear and logistic regression estimator, the elastic net regression estimator, and the regression tree estimator) for finite population estimation of a total or mean from a single stage complex sample design. It also contains the Horvitz-Thompson estimator and several variance estimators.

Installation
------------

You can install mase from github with:

``` r
# install.packages("devtools")
devtools::install_github("Swarthmore-Statistics/mase")
```

Example
-------

Here's an example fitting the Horvitz-Thompson estimator:

``` r

library(mase)

## Estimates the mean and total of the api00 variable using the apisrs dataset in the survey package
library(survey)
#> Loading required package: grid
#> Loading required package: Matrix
#> Loading required package: survival
#> 
#> Attaching package: 'survey'
#> The following object is masked from 'package:graphics':
#> 
#>     dotchart
data(api)
horvitzThompson(y = apisrs$api00, pi = apisrs$pw^(-1), var_est = TRUE, var_method = "HTSRS")
#> $pop_total
#>         [,1]
#> [1,] 4066887
#> 
#> $pop_mean
#>         [,1]
#> [1,] 656.585
#> 
#> $pop_total_var
#> [1] 3282462447
#> 
#> $pop_mean_var
#> [1] 85.55736
```

Please check out the vignette for more examples.
