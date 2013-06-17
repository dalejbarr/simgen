\name{tryFit}
\alias{tryFit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fit model with error-trapping
}
\description{
Try to fit an lmer model, catching any errors/warnings
}
\usage{
tryFit(mf, xdat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mf}{The model formula}
  \item{xdat}{The data frame}
}
%\details{
%}
\value{
  %%  ~Describe the value returned
%%  If it is a LIST, use
\item{value}{fitted model object from lmer}
\item{converged}{whether or not the model converged}
%% ...
}
\seealso{
\code{\link{fitlmer}}, \code{\link{modSpace}}
}
\examples{
mf.wswi <- modSpace(FALSE) # WSWI design

# random intercepts only
ff <- tryFit(mf.wswi$min, mkDf(mcr.params=createParamMx(1)[1,]))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
