\name{fitlmer.mcmc}
\alias{fitlmer.mcmc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Run mixed-effects model (\code{\link{lmer}}) with MCMC p-value
}
\description{
Uses \code{\link{lmer}} to fit a random-intercepts only model, and then uses \code{\link{mcmcsamp}}
  to obtain p-values.
}
\usage{
fitlmer.mcmc(mcr.data, wsbi, nmcmc = 10000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mcr.data}{
A dataframe formatted as described in \code{\link{mkDf}}.
}
  \item{wsbi}{
Whether the design is between-items (TRUE) or within-items (FALSE); has no effect because the model is
random-intercepts only, but was included for consistency with \code{\link{fitlmer}} and \code{\link{fitanova}}.
}
  \item{nmcmc}{
Number of Markov-Chain Monte Carlo simulations (default = 10000).
}
}
\details{
Tries to fit the most complex model requested, and if that fails to
converge, then tries progressively simpler models, before calculating
mcmc p-value.  If no model converges, returns NAs.  The MCMC procedure
is based on Baayen's \code{pvals.fnc} in package \code{languageR}.
}
\value{
  A single element vector with the mcmc p-value.
}
% \references{
%% ~put references to the literature/web site here ~
% }
% \author{
%%  ~~who you are~~
%}
\note{
This function ONLY fits random intercept models.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{fitlmer}}
}
\examples{
pmx <- createParamMx(10)

# between-items dataset
x.bi <- mkDf(nsubj=24, nitem=24, mcr.params=pmx[1,], wsbi=TRUE)

fitlmer.mcmc(x.bi, wsbi=FALSE, nmcmc=1000) # maximal
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
