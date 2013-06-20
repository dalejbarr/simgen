\name{lrCompare}
\alias{lrCompare}
\title{Fit mixed-models and obtain p-values from model comparison}
\usage{
  lrCompare(mcr.data, lrc.mods, lrc.modCompare, ...)
}
\arguments{
  \item{mcr.data}{Data to be analyzed}

  \item{lrc.mods}{List of (named) model formulae}

  \item{lrc.modCompare}{List of comparisons to perform
  (named)}

  \item{...}{Arguments to be passed along to
  \code{\link{tryFit}} (and thus to \code{lmer}).}
}
\value{
  A vector with the p-values from tests specified by
  \code{lrc.modCompare}. If a test could not be performed
  (e.g., because one or both of the corresponding models
  did not converge), \code{NA} is returned.
}
\description{
  A set of model formula are fed to \code{lrCompare} as
  named elements in a list.  These models are fit using
  \code{\link{tryFit}}.  Models are compared as specified
  in the argument \code{lrc.modCompare} using LLR
  (likelihood-ratio) tests. This function is useful for
  performing multiple tests on relatively complex designs
  (e.g., all main effects and interactions in a factorial
  design).
}
\examples{
dat <- mkDf()

lrCompare(dat,
          lrc.mods=list(
              ri=    Resp~Cond+(1|SubjID)+     (1|ItemID),
              ri.noC= Resp~     (1|SubjID)+     (1|ItemID),
              max=    Resp~Cond+(1+Cond|SubjID)+(1+Cond|ItemID),
              max.noC=Resp~     (1+Cond|SubjID)+(1+Cond|ItemID)),
          lrc.modCompare=list(
              ri=c("ri",  "ri.noC"),
              max=c("max","max.noC")),
          REML=FALSE, na.action=na.omit) # pass to tryFit/lmer

dat <- mkDf.facMixedAB(nitem=24)

lrCompare(dat,
          lrc.mods=list(
              m=   Y~Ad*Bd+   (1+Bd|SubjID)+(1+Ad*Bd|ItemID),
              m.A= Y~Bd+Ad:Bd+(1+Bd|SubjID)+(1+Ad*Bd|ItemID),
              m.B= Y~Ad+Ad:Bd+(1+Bd|SubjID)+(1+Ad*Bd|ItemID),
              m.AB=Y~Ad+Bd+   (1+Bd|SubjID)+(1+Ad*Bd|ItemID)),
          lrc.modCompare=list(
              A=c("m","m.A"), B=c("m","m.B"), AB=c("m","m.AB")),
          REML=FALSE, na.action=na.omit)
}

