\name{reassembleStepwiseFile}
\alias{reassembleStepwiseFile}
\title{
Parse and assemble data in a stepwise file into an array.
}
\description{
Convert file containing values from stepwise lmer models to a
three-dimensional array in the workspace.
}
\usage{
reassembleStepwiseFile(fname)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fname}{
Name of the file containing output from \code{\link{mcRun}} calling
either \code{\link{fitstepwise}} or \code{\link{fitstepwise.bestpath}}.
}
}
%\details{
%}
\value{
Returns either a three dimension array (if file is the result of
\code{\link{fitstepwise.bestpath}}) or a list of two three dimensional
arrays if \code{\link{fitstepwise}} was used, one for forward and one
for backward selection.  Each array has the following dimensions:
\item{run}{the Monte Carlo run}
\item{crit}{the alpha level for selection}
\item{params}{the resulting values from individual calls to \code{fitstepwise}
or \code{fitstepwise.bestpath}}
}

\seealso{
\code{\link{fitstepwise}}, \code{\link{fitstepwise.bestpath}}.
}
\examples{
pmx <- createParamMx(10, h0=TRUE)  # 20 experimental runs

mf <- modSpace(wsbi=FALSE) # figure out set of possible models
mf.sfirst <- c(mf[3], mf[[2]][1], mf[1]) # test subject slope first

mcRun("fitstepwise", mcr.outfile="test.txt",
      mcr.xdatFnc="mkDf", mcr.varying=pmx, mf=mf.sfirst, wsbi=FALSE)

ff <- reassembleStepwiseFile("test.txt")
ff$forward[,2,]  # alpha-level=.05
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
