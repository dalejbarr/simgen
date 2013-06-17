\name{fitnocorr}
\alias{fitnocorr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Run mixed-effects model using (\code{\link{lmer}})
}
\description{
Fits a LMEM with random slopes and random intercepts but no random correlations.
}
\usage{
fitnocorr(mcr.data, wsbi)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mcr.data}{
A dataframe formatted as described in \code{\link{mkDf}}.
}
  \item{wsbi}{
Whether the design is between-items (TRUE) or within-items (FALSE).
}
}
%\details{
%}
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
fitnocorr(x.bi, wsbi=TRUE)
fitnocorr(x.wi, wsbi=FALSE)
}
