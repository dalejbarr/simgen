\name{fitanova}
\alias{fitanova}
\title{
Run ANOVA analysis
}
\description{
Given a dataset \code{xd}, calculates F1, F2, F'min, and F1+F2, and
associated p-values.
}
\usage{
fitanova(mcr.data, wsbi)
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
%%%  ~~ If necessary, more details than the description above ~~
%}
\value{
  A vector, with:
  \item{F1}{the by-subjects F value}
  \item{F2}{the by-items F value}
  \item{minF}{the min F' value}
  \item{p1}{the p-value associated with F1}
  \item{p2}{the p-value associated with F2}
  \item{pmf}{the p-value associated with min F'}
  \item{pmax}{the p-value associated with F1+F2 ( =max(c(p1,p2)) )}
}
%\references{
%% ~put references to the literature/web site here ~
%}
%\author{
%%  ~~who you are~~
%}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{mkDf}}
}
\examples{
pmx <- createParamMx(10)

# between-items dataset
x.bi <- mkDf(nsubj=24, nitem=24, mcr.params=pmx[1,], wsbi=TRUE)

# within-items dataset
x.wi <- mkDf(nsubj=24, nitem=24, mcr.params=pmx[1,], wsbi=FALSE)

fitanova(x.bi, wsbi=TRUE)
fitanova(x.wi, wsbi=FALSE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
