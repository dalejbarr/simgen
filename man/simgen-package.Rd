\docType{package}
\name{simgen-package}
\alias{simgen}
\alias{simgen-package}
\title{Monte Carlo Simulation of Simultaneous Generalization to Subjects and Items}
\usage{
  tryFit(tf.formula, tf.data, ...)
}
\arguments{
  \item{mf}{The model formula}

  \item{xdat}{The data frame}
}
\value{
  %% ~Describe the value returned %% If it is a LIST, use

  %% ...

  \item{value} fitted model object from lmer
  \item{converged} whether or not the model converged
}
\description{
  Routines for running Monte Carlo simulations testing the
  generalizability of various ways of analyzing datasets
  with crossed random effects.  This accompanies a paper by
  Barr, Levy, Scheepers, & Tily (under revision).  See
  http://idiom.ucsd.edu/~rlevy/papers/barr-levy-scheepers-tily-8-aug-2011-submitted.pdf.
}
\details{
  \tabular{ll}{ Package: \tab simgen\cr Type: \tab
  Package\cr Version: \tab 1.5\cr Date: \tab 2012-02-20\cr
  License: \tab GPL (>= 2)\cr LazyLoad: \tab yes\cr } To
  get started, see \code{\link{createParamMx}},
  \code{\link{mkDf}}, and \code{\link{mcRun}}.
}
\examples{
# random intercepts only
ff <- tryFit(Resp~Cond + (1+Cond|SubjID)+(1+Cond|ItemID), mkDf())
}
\author{
  Dale J. Barr <dale.barr@glasgow.ac.uk>
}
\seealso{
  \code{\link{fitlmer}}, \code{\link{modSpace}}
}
\keyword{Fit}
\keyword{Try}
\keyword{an}
\keyword{any}
\keyword{catching}
\keyword{error-trapping}
\keyword{errors/warnings}
\keyword{fit}
\keyword{lmer}
\keyword{model}
\keyword{model,}
\keyword{package}
\keyword{to}
\keyword{with}

