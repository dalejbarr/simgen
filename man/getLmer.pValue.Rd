\name{getLmer.pValue}
\alias{getLmer.pValue}
\title{Get p-value from model comparison of lmer objects}
\usage{
  getLmer.pValue(m1, m2)
}
\arguments{
  \item{m1}{A model fitted by \code{\link{tryFit}}}

  \item{m2}{A model fitted by \code{\link{tryFit}}}
}
\value{
  A vector with: \item{chisq}{observed chi-square
  (difference in -2LL for models)} \item{df}{degrees of
  freedom for the comparison} \item{p}{p-value}
}
\description{
  Derive p-value by comparing two models fitted by
  \code{\link{tryFit}}.
}
\examples{
d1 <- mkDf()

m1 <- tryFit(Resp ~ Cond + (1 + Cond | SubjID) + (1 + Cond | ItemID), d1, REML=FALSE,
              na.action=na.omit)
m2 <- tryFit(Resp ~ (1 + Cond | SubjID) + (1 + Cond | ItemID), d1, REML=FALSE,
              na.action=na.omit)

getLmer.pValue(m1, m2)
}

