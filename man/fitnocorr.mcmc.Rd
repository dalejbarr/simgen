\name{fitnocorr.mcmc}
\alias{fitnocorr.mcmc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Run mixed-effects model (\code{\link{lmer}}) with MCMC p-value
}
\description{
Uses \code{\link{lmer}} to fit a random-intercepts only model, and then uses \code{\link{mcmcsamp}}
  to obtain p-values.
}
\usage{
fitnocorr.mcmc(mcr.data, wsbi, nmcmc = 10000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mcr.data}{
A dataframe formatted as described in \code{\link{mkDf}}.
}
  \item{wsbi}{
Whether the design is between-items (TRUE) or within-items (FALSE).
}
  \item{nmcmc}{
Number of Markov-Chain Monte Carlo simulations (default = 10000).
}
}
\details{
If the model does not converge, returns \code{NA}.  The MCMC procedure
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
This function ONLY fits models with independent random intercepts and slopes.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{fitnocorr}}
}
\examples{
pmx <- createParamMx(10)

# between-items dataset
x.bi <- mkDf(nsubj=24, nitem=24, mcr.params=pmx[1,], wsbi=TRUE)

# NB: small number of MCMC runs so that the example runs quickly
# increase the number of runs for stable results
fitnocorr.mcmc(x.bi, wsbi=TRUE, nmcmc=1000) 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
