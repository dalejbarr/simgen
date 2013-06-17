\name{fitstepwise.bestpath}
\alias{fitstepwise.bestpath}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fit lmer model using model selection
}
\description{
Analyzes data using \code{\link{lmer}}, using model selection to test
significance of random slope terms in the model (likelihood ratio
tests).  Does forward or backward selection.  Unlike \code{\link{fitstepwise}}, before
taking a step forward (or backward) both possible slopes are tested, and
the "best" (most conservative) step is taken.
}
\usage{
fitstepwise.bestpath(mcr.data, forward, crit = c(.01,.05,seq(.1,.8,.1)) )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mcr.data}{
A dataframe formatted as described in \code{\link{mkDf}}.  Named thus to interface with the function \code{\link{mcRun}}.
}
  \item{forward}{
Should a forward model (TRUE) or backward model (FALSE) be run?
}
  \item{crit}{
alpha level for each likelihood-ratio test of slope variance.
}
}
\details{
  It only makes sense to run this with data from a within-subject/within-items design (since
  for within-subject/between-items data, there is only one slope to be
  tested, and the procedure is therefore equivalent to a simple stepwise
  algorithm.

  Note that for purposes of efficiency, forward and backward models are
  not run simultaneously as with \code{\link{fitstepwise}}.
}
\value{
  A single row of a dataframe with number of fields depending on \code{crit}.
  Stepwise lmer models output six values for each alpha level (i.e.,
  level of \code{crit}.  These values are:
  \item{fm}{the model that was selected.  4 means that no model converged.  For forward stepping models, 3
  = maximal model was selected; 2 = model includes only first slope; 1 =
model is random intercept only.  For backward stepping models, 1 =
maximal model, 2 = model minus one slope, 3 = random intercept model.}
  \item{t}{t-statistic for the treatment effect}
  \item{chi}{chi-square statistic for the likelihood ratio test (1 df)}
  \item{pt}{p-value for the t-statistic (normal distribution)}
  \item{pchi}{p-value for the chi-square statistic}
  The values are given once for the specified direction (forward or backward).  Thus, if 
  there are two values of crit, there will be 6 (value)
  x 2 (levels of crit) = 12 values in each row of the dataframe.  To
  assemble a dataframe of results from a file into a three-dimensional array, see \code{\link{reassembleStepwiseFile}}.
}
% \references{
%% ~put references to the literature/web site here ~
% }
% \author{
%%  ~~who you are~~
% }
% \note{
%%  ~~further notes~~
% }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{fitlmer}}, \code{\link{modSpace}},
\code{\link{fitstepwise}}, \code{\link{reassembleStepwiseFile}}
}
\examples{
pmx <- createParamMx(10)

x8 <- mkDf(nsubj=24, nitem=24, pmx[8,], wsbi=FALSE)

# forward
fitstepwise.bestpath(mcr.data=x8, forward=TRUE, crit=.05)

# backward
fitstepwise.bestpath(mcr.data=x8, forward=FALSE, crit=.05)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
