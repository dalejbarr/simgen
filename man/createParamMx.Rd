% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/onefactor.R
\name{createParamMx}
\alias{createParamMx}
\title{Create Population Parameter Matrix for One-Factor Design}
\usage{
createParamMx(nexp = 1e+05, simparam.env = getParamRanges(),
  firstseed = 666, h0 = TRUE, outfile = NULL)
}
\arguments{
\item{nexp}{Number of experiments to run (default 100000).}

\item{simparam.env}{An \code{\link{environment}} object containing the
parameter ranges; see \code{\link{getParamRanges}} (default) for the
appropriate format.}

\item{firstseed}{First seed to start off generation of the matrix.  Included
for replicability of results.  A parameter matrix of a given size with a
given seed will always be identical.}

\item{h0}{Status of the null hypothesis: TRUE or FALSE.}

\item{outfile}{Name of a file where the matrix will be stored; if NULL
(default), no data will be stored.}
}
\value{
A matrix, with the following columns:

\item{int}{(grand mean) intercept}

\item{eff}{treatment effect}

\item{err}{error variance}

\item{miss}{proportion of response values missing}

\item{t00}{by-subject intercept variance}

\item{t11}{by-subject slope variance}

\item{rsub}{by-subject intercept/slope correlation}

\item{w00}{by-item intercept variance}

\item{w11}{by-item slope variance}

\item{ritm}{by-item intercept/slope correlation}

\item{seed}{random number seed for creating the dataframe (see}

\code{\link{mkDf}})
}
\description{
Create a matrix of population parameters for Monte Carlo simulation
with one-factor design.
}
\details{
Each row of the matrix returned from this function represents the
population parameters for a given experiment.  As of version 1.6 of
simgen, This function has been superceded by
\code{\link{genParamRanges}} and \code{\link{randParams}}, but is
preserved for backwards compatibility.
}
\examples{

# using defaults
createParamMx(10)

# now let's change one of the ranges
p.env <- getParamRanges()
get("t11.range", env=p.env)
assign("t11.range", c(5,10), env=p.env)
createParamMx(10, simparam.env=p.env)

}
\seealso{
\code{\link{getParamRanges}}, \code{\link{mkDf}}
}
