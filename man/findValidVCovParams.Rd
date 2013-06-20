\name{findValidVCovParams}
\alias{findValidVCovParams}
\title{Find valid parameters for a variance-covariance matrix}
\usage{
  findValidVCovParams(vars,
    corrs = as.list(rep(0, sum(1:(length(vars) - 1)))),
    maxtries = 10000)
}
\arguments{
  \item{vars}{(named) list of ranges of variances; if not
  named, given names V1..VN}

  \item{corrs}{(named) list of correlations (range -1 to
  1); if not named, given names C1..CN.  List them in the
  following order: take the upper triangle of the matrix,
  list all elements in the first row (left to right), then
  the second row, etc.}

  \item{maxtries}{the number of attempts before giving up}
}
\value{
  A vector of variances and correlations
}
\description{
  To generate data from a multivariate normal distribution
  using \code{mvrnorm} in the \code{MASS} package, it is
  necessary to specify a square variance-covariance matrix
  that is symmetric and 'positive definite.'
}
\details{
  This function takes data-generating parameter ranges
  (variances and covariances) and randomly generates values
  for a single population. Use a brute-force method, so to
  avoid hanging, set \code{maxtries} to a reasonable value.

  A 3x3 matrix corresponds to 3 variances on the diagonal,
  and 6 variances (forming the upper or lower triangle).
  Stated generally, an NxN matrix requires
  \code{sum(1:(N-1))} covariances.
}
\examples{
pars <- findValidVCovParams(vars=list(A=c(1,3),B=c(1,3), C=c(0,2)),
                            corrs=list(AB=c(-.5,.5),AC=c(-.5,.5),BC=c(-.5,.5)))

formVCovMx(pars[1:3], pars[4:6])
}
\seealso{
  \code{\link{formVCovMx}}
}

