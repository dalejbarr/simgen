\name{heatmap2}
\alias{heatmap2}
\title{Create a heatmap from a given loess fit.}
\usage{
  heatmap2(lofit, labelX, labelY, title,
    mxseq = seq(0.01, 2.99, length.out = 300),
    myseq = seq(0.01, 2.99, length.out = 300),
    powermap = FALSE)
}
\arguments{
  \item{lofit}{fitted loess object}

  \item{labelX}{label for the x-axis variable.  This can be
  a character string or an expression.  If one of "t00",
  "t11", "w00", "w11", will be properly formatted using
  \code{\link{plotmath}}.}

  \item{labelY}{same as \code{labelX}, but for the y axis.}

  \item{title}{title for the plot}

  \item{mxseq}{sequence for the x variable, determining
  points at which predictions are calculated}

  \item{myseq}{same as \code{mxseq}, but for the y
  variable}

  \item{powermap}{if plotting power (TRUE) or type I error
  (FALSE); allows for different grayscale sensitivities.}
}
\description{
  Creates a graph of estimated p-values from a loess fit to
  stored simulation data.  The loess fit should be created
  using \code{\link{loessPred}}.
}
\examples{
pmx <- createParamMx(100)

result <- apply(pmx, 1, function(x) {fitanova(mkDf(24,24,mcr.params=x,FALSE),FALSE)})
lpred <- loessPred(pval=result["pmax",],
                paramX=pmx[,"t11"], paramY=pmx[,"w11"])

heatmap2(lpred, "t11", "w11", "F1+F2")
}
\seealso{
  \code{\link{loessPred}}
}

