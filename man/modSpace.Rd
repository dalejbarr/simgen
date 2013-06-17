\name{modSpace}
\alias{modSpace}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Derive model formulas
}
\description{
Derive all random effects structures for a given design, from maximal
to random intercepts only.
}
\usage{
modSpace(wsbi)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{wsbi}{
Whether the design is within items (FALSE) or between items (TRUE).
}
}
%\details{
%}
\value{
A list of model formulas.  The model formulas can be used as arguments to lmer.
  %%  ~Describe the value returned
%%  If it is a LIST, use
\item{max}{'maximal' model}
\item{mid}{'middle' model with a single slope, either by-subject
  (mid$srs) or by-item (mid$irs) (WSWI only)}
\item{min}{random intercepts only model}
%% ...
}
\seealso{
\code{\link{fitlmer}}
}
\examples{
mf.wswi <- modSpace(FALSE) # WSWI design

# subject random slope
ff <- tryFit(mf.wswi$mid$srs, mkDf(mcr.params=createParamMx(1)[1,]))
ff$converged # did it converge?

# random intercepts only
ff <- tryFit(mf.wswi$min, mkDf(mcr.params=createParamMx(1)[1,]))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
