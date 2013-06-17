\name{fitlmer}
\alias{fitlmer}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Run mixed-effects model using (\code{\link{lmer}})
}
\description{
Runs either a random-intercepts only model or a maximal random-effects
model (by-subject and by-item random intercepts, by-subject random
slope, and by-item random slope for within-item design).
}
\usage{
fitlmer(mcr.data, ri.only = FALSE, wsbi = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mcr.data}{
A dataframe formatted as described in \code{\link{mkDf}}.
}
  \item{ri.only}{
Whether the random effects specification is to be random-intercepts only (TRUE) or
maximal random-effects (FALSE).
}
  \item{wsbi}{
Whether the design is between-items (TRUE) or within-items (FALSE).
}
}
\details{
  \code{fitlmer} will attempt to fit the model specified by the user,
  and will progressively simplify the model as needed to get it to
  converge.  If no model converges, it returns NAs.  \code{fitlmer}
  performs a likelihood ratio test for the treatment effect, as well as
  returns a p-value for the t-statistic using an approximation from the
  normal distribution.
}
\value{
  \item{fm}{Code for the model that converged: (1) dropped one slope; (2)
    dropped two slopes; (3) main model did not converge; (4) comparision
    model for likelihood ratio test did not converge.}
  \item{t}{t-statistic for the treatment effect}
  \item{chi}{chi-square statistic for the likelihood ratio test (1 df)}
  \item{pt}{p-value for the t-statistic (normal distribution)}
  \item{pchi}{p-value for the chi-square statistic}
}
% \references{
%% ~put references to the literature/web site here ~
% }
%\author{
%%  ~~who you are~~
%}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{mkDf}}, \code{\link{fitlmer.mcmc}}
}
\examples{
pmx <- createParamMx(10)

# between-items dataset
x.bi <- mkDf(nsubj=24, nitem=24, mcr.params=pmx[1,], wsbi=TRUE)

# within-items dataset
x.wi <- mkDf(nsubj=24, nitem=24, mcr.params=pmx[1,], wsbi=FALSE)

# maximal model
fitlmer(x.bi, wsbi=TRUE)

fitlmer(x.wi, wsbi=FALSE, ri.only=FALSE) # maximal
fitlmer(x.wi, wsbi=FALSE, ri.only=TRUE) # random intercepts only
}
