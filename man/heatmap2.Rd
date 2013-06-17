\name{heatmap2}
\alias{heatmap2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create a heatmap from a given loess fit.
}
\description{
Creates a graph of estimated p-values from a loess fit to stored
simulation data.  The loess fit should be created using
\code{\link{loessPred}}.
}
\usage{
heatmap2(lofit, labelX, labelY, title, mxseq = seq(0.01, 2.99, length.out = 300), myseq = seq(0.01, 2.99, length.out = 300), powermap = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lofit}{
fitted loess object
}
  \item{labelX}{
label for the x-axis variable.  This can be a character string or an
expression.  If one of "t00", "t11", "w00", "w11", will be
properly formatted using \code{\link{plotmath}}.
}
  \item{labelY}{
same as \code{labelX}, but for the y axis.
}
  \item{title}{
title for the plot
}
  \item{mxseq}{
sequence for the x variable, determining points at which predictions are calculated
}
  \item{myseq}{
same as \code{mxseq}, but for the y variable
}
  \item{powermap}{
if plotting power (TRUE) or type I error (FALSE); allows for different
grayscale sensitivities.
}
}
%\details{
%Fits a loess model using function \code{\link{loess}} with \code{span}=.9.
%}
%\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
%}
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
\code{\link{loessPred}}
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
