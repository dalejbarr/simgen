\name{reassembleStepwiseFile}
\alias{reassembleStepwiseFile}
\title{Parse and assemble data in a stepwise file into an array.}
\usage{
  reassembleStepwiseFile(fname)
}
\arguments{
  \item{fname}{Name of the file containing output from
  \code{\link{mcRun}} calling either
  \code{\link{fitstepwise}} or
  \code{\link{fitstepwise.bestpath}}.}
}
\value{
  Returns either a three dimension array (if file is the
  result of \code{\link{fitstepwise.bestpath}}) or a list
  of two three dimensional arrays if
  \code{\link{fitstepwise}} was used, one for forward and
  one for backward selection.  Each array has the following
  dimensions: \item{run}{the Monte Carlo run}
  \item{crit}{the alpha level for selection}
  \item{params}{the resulting values from individual calls
  to} \code{fitstepwise} or \code{fitstepwise.bestpath}
}
\description{
  Convert file containing values from stepwise lmer models
  to a three-dimensional array in the workspace (one-factor
  design only).
}
\examples{
nmc <- 10 # number of monte carlo simulations
# create the parameter matrix
pmx <- cbind(randParams(genParamRanges(), nmc, 1001), seed=mkSeeds(nmc, 1001))

# mf contains model formulae for the various models to fit
mf <- list(max=Resp ~ Cond + (1 + Cond | SubjID) + (1 + Cond | ItemID),
            mid.srs=Resp ~ Cond + (1 + Cond | SubjID) + (1 | ItemID),
            min=Resp ~ Cond + (1 | SubjID) + (1 | ItemID))

mcRun("fitstepwise", mcr.outfile="test.txt",
      mcr.fnArgs=list(wsbi=FALSE, mf=mf),
      mcr.datFn="mkDf", mcr.datArgs=list(wsbi=FALSE),
      mcr.varying=pmx)

ff <- reassembleStepwiseFile("test.txt")
ff$forward[,2,]  # alpha-level=.05
}
\seealso{
  \code{\link{fitstepwise}},
  \code{\link{fitstepwise.bestpath}}.
}

