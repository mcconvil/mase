# MASE
Model-Assisted Survey Estimators


## Description

Package: MASE

Type: Package

Title: Model-Assisted Survey Estimators

Version: 0.1

Date: 2016-02-14

Author: Kelly McConville, Becky Tang, George Zhu

Maintainer: Kelly McConville <kmcconv1@swarthmore.edu>


## Reference manual

[Vignette](https://github.swarthmore.edu/xzhu1/MASE/blob/master/vignettes/Model-Assisted%20Survey%20Estimators.Rmd)


## Version Notes (Need to revise)

### 0.1

* Initial commit 

* Added Horvitz-Thompson method

* Added GREG (Linear and Ridge)

* Added LASSO (Lasso and Ridge) 

* Added Ratio method

* Added Post-stratification method

* Added vignette


## Known Issues

### 0.1

* Issue on stratified weights for postStratification, g * weight overestimates population total


## Future Improvements

### 0.1

* Make inputs for postStrat more user friendly 

* Create a helper function to convert population totals dataframe to the nested list structure

* Design new search algorithm using Binary Search for Ridge Method in GREG in weights cacluation

* Add more testfiles in tests folder

* Modularize certain code chunks for abstraction

* Create a variance estimation function that can be used for all estimators instead of being estimator specific.  Allow users to decide when they run the function if they also want a variance estimator.

* Consider adding the g-weighted variance estimator of Sarndal, Swensson and Wretman (1989)

* Change weight to weights.  Better yet.  Make weight = pi and then make sure to invert pi in functions.

* Need to change horvitzThompson to reflect changes in var.horvitzThompson.

* Need to deal with issue of no N AND no pis in var.horvitzThompson.

* Need to make var.MASE give warnings for issues with pi2.


