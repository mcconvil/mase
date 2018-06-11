## Test Environments
* local OS X El Capitan, R 3.5.0
* win-builder, R-devel (R Under development (unstable) (2018-06-05 r74852))

## R CMD check results for local OS X
0 errors | 0 warnings | 0 notes

## R CMD check results for win-builder
There were no errors or warnings.  There was one note:

* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Kelly McConville <mcconville@reed.edu>'

New submission

****************************************************

#Feedback from initial submission:

Am 06.06.2018 um 22:39 schrieb CRAN submission:
> [This was generated from CRAN.R-project.org/submit.html]
> 
> The following package was uploaded to CRAN:
> ===========================================
> 
> Package Information:
> Package: mase
> Version: 0.1.1
> Title: Model-Assisted Survey Estimators
> Author(s): Kelly McConville, Becky Tang, George Zhu, Sida Li, Shirley
>    Cheung
> Maintainer: Kelly McConville <mcconville@reed.edu>
> Depends: R (>= 3.1)
> Suggests: NHANES, roxygen2, testthat, knitr, rmarkdown, readr, ggplot2
> Description: A set of model-assisted survey estimators and corresponding
>    variance estimators for single stage, unequal probability
>    sampling designs without replacement. All of the estimators
>    can be written as a generalized  regression estimator.

Thanks, please add a reference for the methods in the 'Description' 
field of your DESCRIPTION file in the form
authors (year) <doi:...>
authors (year) <arXiv:...>
authors (year, ISBN:...)
with no space after 'doi:', 'arXiv:' and angle brackets for auto-linking.

We are missing Daniell Toth in the authors list. Please add all authors 
and copyright holders in the Authors@R field with the appropriate roles.

 From the CRAN policies you agreed to:

"The ownership of copyright and intellectual property rights of all 
components of the package must be clear and unambiguous (including from 
the authors specification in the DESCRIPTION file). Where code is copied 
(or derived) from the work of others (including from R itself), care 
must be taken that any copyright/license statements are preserved and 
authorship is not misrepresented.

Preferably, an ‘Authors@R’ would be used with ‘ctb’ roles for the 
authors of such code. Alternatively, the ‘Author’ field should list 
these authors as contributors.

Where copyrights are held by an entity other than the package authors, 
this should preferably be indicated via ‘cph’ roles in the ‘Authors@R’ 
field, or using a ‘Copyright’ field (if necessary referring to an 
inst/COPYRIGHTS file)."

