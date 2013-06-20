\name{psig}
\alias{psig}
\title{Proportion of runs achieving statistical significance}
\usage{
  psig(x, alpha = 0.05)
}
\arguments{
  \item{x}{Vector or matrix where each row is a single run}

  \item{alpha}{alpha level (default=.05)}
}
\value{
  A vector with the proportion of significant runs.
}
\description{
  Proportion of runs achieving statistical significance
}
\examples{
nmc <- 10
pmx <- cbind(randParams(genParamRanges(), nmc, 1001), seed=mkSeeds(nmc, 1001))

res <- mcRun("fitanova", mcr.fnArgs=list(wsbi=FALSE),
             mcr.datFn="mkDf", mcr.datArgs=list(wsbi=FALSE),
             mcr.varying=pmx)

psig(res[,c("p1","p2","pmf","pmax")])
}
\seealso{
  \code{\link{pconv}}
}

