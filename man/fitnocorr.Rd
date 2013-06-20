\name{fitnocorr}
\alias{fitnocorr}
\title{Run mixed-effects model using (\code{\link{lmer}})}
\usage{
  fitnocorr(mcr.data, wsbi = FALSE)
}
\arguments{
  \item{mcr.data}{A dataframe formatted as described in
  \code{\link{mkDf}}.}

  \item{wsbi}{Whether the design is between-items (TRUE) or
  within-items (FALSE).}
}
\value{
  A vector with elements: \item{fm} Code for the model that
  converged: (1) dropped one slope; (2) dropped two slopes;
  (3) main model did not converge; (4) comparision model
  for likelihood ratio test did not converge. \item{t}
  t-statistic for the treatment effect \item{chi}
  chi-square statistic for the likelihood ratio test (1 df)
  \item{pt} p-value for the t-statistic (normal
  distribution) \item{pchi} p-value for the chi-square
  statistic
}
\description{
  Fits a LMEM with random slopes and random intercepts but
  no random correlations (one-factor design only).
}
\examples{
nmc <- 10
pmx <- cbind(randParams(genParamRanges(), nmc, 1001), seed=mkSeeds(nmc, 1001))

# between-items dataset
x.bi <- mkDf(nsubj=24, nitem=24, mcr.params=pmx[1,], wsbi=TRUE)

# within-items dataset
x.wi <- mkDf(nsubj=24, nitem=24, mcr.params=pmx[1,], wsbi=FALSE)

# maximal model
fitnocorr(x.bi, wsbi=TRUE)
fitnocorr(x.wi, wsbi=FALSE)
}
\seealso{
  \code{\link{mkDf}}, \code{\link{fitlmer.mcmc}}
}

