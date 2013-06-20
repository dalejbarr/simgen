\name{tryFit}
\alias{tryFit}
\title{Fit model with error-trapping}
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

  \item{value}{fitted model object from lmer}
  \item{converged}{whether or not the model converged}
}
\description{
  Try to fit an lmer model, catching any errors/warnings
}
\examples{
# random intercepts only
ff <- tryFit(Resp~Cond + (1+Cond|SubjID)+(1+Cond|ItemID), mkDf())
}
\seealso{
  \code{\link{fitlmer}}, \code{\link{modSpace}}
}

