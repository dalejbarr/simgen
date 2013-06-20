\name{fitlmer.mcmc}
\alias{fitlmer.mcmc}
\title{Run mixed-effects model (\code{\link{lmer}}) with MCMC p-value}
\usage{
  fitlmer.mcmc(mcr.data, wsbi, nmcmc = 10000)
}
\arguments{
  \item{mcr.data}{A dataframe formatted as described in
  \code{\link{mkDf}}.}

  \item{wsbi}{Whether the design is between-items (TRUE) or
  within-items (FALSE); has no effect because the model is
  random-intercepts only, but was included for consistency
  with \code{\link{fitlmer}} and \code{\link{fitanova}}.}

  \item{nmcmc}{Number of Markov-Chain Monte Carlo
  simulations (default = 10000).}
}
\value{
  A single element vector with the mcmc p-value.
}
\description{
  Uses \code{\link{lmer}} to fit a random-intercepts only
  model, and then uses \code{\link{mcmcsamp}} to obtain
  p-values (one-factor design only).
}
\details{
  Tries to fit the most complex model requested, and if
  that fails to converge, then tries progressively simpler
  models, before calculating mcmc p-value.  If no model
  converges, returns NAs.  The MCMC procedure is based on
  Baayen's \code{pvals.fnc} in package \code{languageR}.
}
\note{
  This function ONLY fits random intercept models.
}
\examples{
nmc <- 10
pmx <- cbind(randParams(genParamRanges(), nmc, 1001), seed=mkSeeds(nmc, 1001))

# between-items dataset
x.bi <- mkDf(nsubj=24, nitem=24, mcr.params=pmx[1,], wsbi=TRUE)

fitlmer.mcmc(x.bi, wsbi=FALSE, nmcmc=1000) # maximal
}
\seealso{
  \code{\link{fitlmer}}
}

