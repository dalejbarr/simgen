\name{loessPred}
\alias{loessPred}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fit a loess model to the p-value data for later use in generating heatmaps.
}
\description{
Given a vector of p-values and two population parameters, this will fit
a \code{\link{loess}} model to the data, with the model of the form
z ~ x * y, where z is the predicted p-value.  This can later be used to
make a heatmap, and is separated from the \code{\link{heatmap2}} function because of the long amount of time it
can take for this function to run.
}
\usage{
loessPred(pval, paramX, paramY, minx = 0.01, maxx = 2.99, miny = 0.01,
maxy = 2.99, svec = rep(TRUE, length(paramX)), span = .9)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pval}{
vector of p-values from simulated runs
}
  \item{paramX}{
population parameter to condition \code{z} on (e.g., a column from a
matrix created using \code{\link{createParamMx}}).  This will appear on
the x-axis of the heatmap.
}
  \item{paramY}{
same as \code{paramX}, but for the Y-axis
}
  \item{minx}{
minimum x value
}
  \item{maxx}{
maximum x value
}
  \item{miny}{
minimum y value
}
  \item{maxy}{
maximum y value
}
  \item{svec}{
a logical vector to select out only certain rows from \code{pval},
\code{paramX} and \code{paramY}
}
\item{span}{
  a value to be passed along to the \code{\link{loess}} function that
  determines the smoothing level
  }
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
A fitted \code{\link{loess}} model.
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
\code{\link{loess}}, \code{\link{heatmap2}}
}
\examples{
pmx <- createParamMx(100)

result <- apply(pmx, 1, function(x) {fitanova(mkDf(24,24,mcr.params=x,FALSE),FALSE)})
lpred <- loessPred(pval=result["pmax",],
                paramX=pmx[,"t11"], paramY=pmx[,"w11"])

heatmap2(lpred, "t11", "w11", "F1+F2")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
