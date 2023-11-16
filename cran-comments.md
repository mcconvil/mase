# Resubmission

This fixes the issues found in the previous submission which updated mase from 0.1.4 to 0.1.5

* I updated some of the formatting of the function documentation which had become outdated in Webian and Windows

The package has been re-tested on the following systems

* macos-latest (release)
* windows-latest (release)
* ubuntu-latest (devel)
* ububtu-latest (release)
* ubuntu-latest (oldrel-1)

And got 

0 errors | 0 warnings | 0 notes

We checked the strong reverse dependency: FIESTAutils, comparing R CMD check results across CRAN and dev versions of this package.

 * We saw 0 new problems
 * We failed to check 0 packages



# Resubmission

This is an update to mase 0.1.4. The new version is mase 0.1.5.

The package as tested on the following systems

* macos-latest (release)
* windows-latest (release)

And got 

0 errors | 0 warnings | 0 notes

As well as

* ubuntu-latest (devel)
* ububtu-latest (release)
* ubuntu-latest (oldrel-1)

And got 

And got

0 errors | 0 warnings | 1 note

Where the note had to do with the size of the installed package.



# Resubmission

This is an update to mase 0.1.3. The new version is mase 0.1.4.

Notably, the primary author and maintainer, Kelly McConville, has moved from Reed College to Harvard University and so her email has been changed from mcconville@reed.edu to kmcconville@fas.harvard.edu. Kelly's profile at Harvard can be found at the following link: https://statistics.fas.harvard.edu/people/kelly-mcconville and here biography that details this move can be found in the biography section of here personal website: https://mcconville.rbind.io/.

The package was tested on the following systems

* local MacOS Ventura (13.4.1), R 4.3.1
* win builder, R Under development (unstable) (2023-08-23 r85001 ucrt)

And got

0 errors | 0 warnings | 1 note

Where the note regarded the change in maintainer email.

# Resubmission

This is an update to mase 0.1.2.  The new version is mase 0.1.3.

Received the following message from CRAN:

 Version: 0.1.2
Check: LazyData
Result: NOTE
     'LazyData' is specified without a 'data' directory
Flavors: r-devel-linux-x86_64-debian-clang, r-devel-linux-x86_64-debian-gcc, r-devel-linux-x86_64-fedora-clang, r-devel-linux-x86_64-fedora-gcc, r-devel-windows-x86_64, r-devel-windows-x86_64-gcc10-UCRT, r-patched-linux-x86_64, r-patched-solaris-x86, r-release-linux-x86_64, r-release-macos-arm64, r-release-macos-x86_64, r-release-windows-ix86+x86_64

Version: 0.1.2
Check: examples
Result: ERROR
    Running examples in 'mase-Ex.R' failed
    The error most likely occurred in:
    
    > base::assign(".ptime", proc.time(), pos = "CheckExEnv")
    > ### Name: gregTree
    > ### Title: Compute a regression tree estimator
    > ### Aliases: gregTree
    >
    > ### ** Examples
    >
    > library(survey)
    Loading required package: grid
    Loading required package: Matrix
    Loading required package: survival
    
    Attaching package: 'survey'
    
    The following object is masked from 'package:graphics':
    
     dotchart
    
    > data(api)
    > gregTree(y = apisrs$api00,
    + x_sample = apisrs[c("col.grad", "awards", "snum", "dnum", "cnum", "pcttest", "meals", "sch.wide")],
    + x_pop = apipop[c("col.grad", "awards", "snum", "dnum", "cnum", "pcttest", "meals", "sch.wide")])
    Assuming simple random sampling
    Error in h(simpleError(msg, call)) :
     error in evaluating the argument 'x' in selecting a method for function 't': 'x' must be an array of at least two dimensions
    Calls: gregTree -> t -> colSums -> colSums -> <Anonymous>
    Execution halted
Flavors: r-devel-linux-x86_64-debian-clang, r-devel-linux-x86_64-debian-gcc, r-patched-linux-x86_64, r-release-linux-x86_64 


----------------------------------------------

* I removed the LazyData option from DESCRIPTION.
* I updated the functions that depend on rpms to now be compatiable with rpms.


----------------------------------------------



I tested the package on the following systems:

macOS 10.13.6 High Sierra, R-release, brew
Debian Linux, R-devel, clang, ISO-8859-15 locale
Windows Server 2008 R2 SP1, R-release, 32/64 bit

And got:

0 errors ✓ | 0 warnings ✓ | 0 notes ✓


-------------------------------------------------------------------

# Resubmission

This is an update to mase 0.1.1.  The new version is mase 0.1.2.

Received the following message from CRAN:

"This concerns packages

GGUM Rdpack erhcv ggtern idefix mase partitionComparison pdSpecEst

'check' in R-devel now checks Rd files after evaluating \Sexpr{}
expressions.  As the results pages on CRAN now show (for Fedora so far),
there are problems either with the Rd code inserted or with evaluating
the expression.  Note in particular that inserting non-ASCII text
requires an encoding to be declared (for the Rd file or the package).

Unfortunately, the line number currently reported refers to the \Sexpr
evaluation and not the text of the unexpanded .Rd file.

Please correct ASAP and before Oct 16 to safely retain the package on CRAN."

I have added "Encoding: UTF-8" to the DESCRIPTION file and that fixed the issue.  

Again I tested the package on two systems:

* local OS High Sierra (10.13.6), R 3.5.1
* win builder, R Under development (unstable) (2018-10-10 r75427)

With results:

R CMD check results
0 errors | 0 warnings | 0 notes




# Resubmission

This is a resubmission. Thank you for the helpful feedback to my initial submission. In this version I have:

* Again tested the package on two systems:
    + local OS Sierra (10.12.6), R 3.5.0
    + win-builder, R-devel (R Under development (unstable) (2018-06-10 r74877))

* Used Authors@R in the DESCRIPTION file.  I added Daniell Toth with the roles "ctb" and "cph" and the comment "Author and copyright holder of treeDesignMatrix helper function".
* I have added citations to the Description in the DESCRIPTION file and used the appropriate syntax.  
    + When I checked the package, I got the following message.  All of these are correctly spelled.
    
Possibly mis-spelled words in DESCRIPTION:
  Horvitz (19:35)
  Mashreghi (27:50)
  McConville (23:24, 24:44)
  Sarndal (20:41)
  Tille (26:28)
  Toth (24:59)
  al (20:52, 23:38, 27:63)
  et (20:49, 23:35, 27:60)
  
* Again, there is only one note:

* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Kelly McConville <mcconville@reed.edu>'

New submission

* I have included the initial submission note and feedback from CRAN reviewer below.


#Initial Submission

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



