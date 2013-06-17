\name{getParamRanges}
\alias{getParamRanges}
\title{
Get the parameter ranges used in the Monte Carlo simulations.
}
\description{
Retrieve the parameter ranges used in the simulations.
}
\usage{
getParamRanges()
}
%- maybe also 'usage' for other objects documented here.
\details{
This function is also useful for generating parameter ranges that can
then be modified by the user.
}
\value{
An environment object, containing the following variables:
\item{evar.range}{range of error variance (default [0,3]) }
\item{icept.range}{range of intercept (default [-3, 3]) }
\item{pmissing.range}{range of proportion missing obs (default [0, .05])
}
\item{r01.item.range}{range of by-subject intercept/slope correlation
  (default [-.8, .8])}
\item{r01.subj.range}{range of by-item intercept/slope correlation
  (default [-.8, .8])}
\item{slope}{treatment effect, when h0 TRUE, 0; when h0 FALSE, .8}
\item{t00.range}{range of by-subject random intercept variance (default [0,
3]}
\item{t11.range}{range of by-subject random slope variance (default [0, 3])}
\item{w00.range}{range of by-item random intercept variance (default [0,
3])}
\item{w11.range}{range of by-item random slope variance (default [0, 3])}
}
% \references{
%% ~put references to the literature/web site here ~
% }
% \author{
%%  ~~who you are~~
% }
% \note{
%%  ~~further notes~~
% }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{createParamMx}}
}
\examples{
# display all the parameter ranges
print(as.list(getParamRanges()))

# example of how to alter one of the values
p.env <- getParamRanges()
p.env$pmissing.range <- c(0, .25)
assign("pmissing.range", c(0, .25), env=p.env) # alternate way
print(as.list(p.env))
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
