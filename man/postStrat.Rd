% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/postStrat.R
\name{postStrat}
\alias{postStrat}
\title{Compute a post-stratified estimator}
\usage{
postStrat(
  y,
  xsample,
  xpop,
  pi = NULL,
  N = NULL,
  var_est = FALSE,
  var_method = "HB",
  pi2 = NULL,
  datatype = "raw",
  B = 1000,
  fpc = TRUE
)
}
\arguments{
\item{y}{A numeric vector of the sampled response variable.}

\item{xsample}{A vector containing the post-stratum for each sampled unit.}

\item{xpop}{A vector or data frame, depending on datatype.  If datatype = "raw", then a vector containing the post-stratum for each population unit.  If datatype = "totals" or "means", then a data frame, where the first column lists the possible post-strata and the second column contains the population total or proportion in each post-stratum.}

\item{pi}{A numeric vector of inclusion probabilities for each sampled unit in y.  If NULL, then simple random sampling without replacement is assumed.}

\item{N}{A numeric value of the population size. If NULL, it is estimated with the sum of the inverse of the pis.}

\item{var_est}{Default to FALSE, logical for whether or not to compute estimate of variance}

\item{var_method}{The method to use when computing the variance estimator.  Options are a Taylor linearized technique: "LinHB"= Hajek-Berger estimator, "LinHH" = Hansen-Hurwitz estimator, "LinHTSRS" = Horvitz-Thompson estimator under simple random sampling without replacement, and "LinHT" = Horvitz-Thompson estimator or a resampling technique: "bootstrapSRS" = bootstrap variance estimator under simple random sampling without replacement, "SRSunconditional" = simple random sampling variance estimator which accounts for random strata.}

\item{pi2}{A square matrix of the joint inclusion probabilities.  Needed for the "LinHT" variance estimator.}

\item{datatype}{Default to "raw", takes values "raw", "totals" or "means" for whether the user is providing the raw population stratum memberships, the population totals of each stratum, or the population proportions of each stratum.}

\item{B}{The number of bootstrap samples if computing the bootstrap variance estimator.  Default is 1000.}

\item{fpc}{Default to TRUE, logical for whether or not the variance calculation should include a finite population correction when calculating the "LinHTSRS", "SRSunconditional" or the "SRSbootstrap" variance estimator.}
}
\value{
A list of output containing:
\itemize{
\item{pop_total:}{Estimate of population total}
\item{pop_mean:}{Estimate of the population mean (or proportion)}
\item{pop_total_var:}{Estimated variance of population total estimate}
\item{pop_mean_var:}{Estimated variance of population mean estimate}
\item{strat_ests:}{Table of total and mean estimates for each strata}
\item{weights:}{Survey weights produced by PS}
}
}
\description{
Calculates a post-stratified estimator for a finite population mean/proportion or total based on sample data collected from a complex sampling design and a single, categorical auxiliary population variable.
}
\examples{
library(survey)
data(api) 
postStrat(y = apisrs$api00, xsample = apisrs$awards, 
xpop = data.frame(table(apipop$awards)), datatype = "totals", 
pi = apisrs$pw^(-1))

}
\references{
\insertRef{coc77}{mase} 

\insertRef{sar92}{mase}
}