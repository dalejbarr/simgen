\name{createParamMx}
\alias{createParamMx}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create Population Parameter Matrix
}
\description{
Create a matrix of population parameters for Monte Carlo simulations.
}
\usage{
createParamMx(nexp = 100000, simparam.env = getParamRanges(), firstseed = 666, h0 = TRUE, outfile = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nexp}{
Number of experiments to run (default 100000).
}
  \item{simparam.env}{
An \code{\link{environment}} object containing the parameter ranges;
see \code{\link{getParamRanges}} (default) for the appropriate format.
}
  \item{firstseed}{
First seed to start off generation of the matrix.  Included for
replicability of results.  A parameter matrix of a given size with a
given seed will always be identical.
}
  \item{h0}{
Status of the null hypothesis: TRUE or FALSE.
}
  \item{outfile}{
Name of a file where the matrix will be stored; if NULL (default), no
data will be stored.
}
}
\details{
Each row of the matrix returned from this function represents the
population parameters for a given experiment.  Using
\code{firstseed}=666 and \code{nexp}=100000 will yield the same
parameter matrix used in the manuscript.
}
\value{
A matrix, with the following columns:
%%  ~Describe the value returned
%%  If it is a LIST, use
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
\item{seed}{random number seed for creating the dataframe (see \code{\link{mkDf}})}
%% ...
}
\seealso{
\code{\link{getParamRanges}}, \code{\link{mkDf}}
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
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
