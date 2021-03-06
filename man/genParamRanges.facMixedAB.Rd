% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/factorial.R
\name{genParamRanges.facMixedAB}
\alias{genParamRanges.facMixedAB}
\title{Data-generating parameter ranges}
\usage{
genParamRanges.facMixedAB()
}
\value{
A list containing ranges for the following parameters, default values in parentheses:

\item{t_00}{by-subject random intercept variance (1 to 3)}
\item{t_11}{by-subject random slope variance (1 to 3)}
\item{r_01}{by-subject random intercept/slope correlation (-.9,.9)}
\item{evar}{error variance (fixed at 3)}
\item{mu}{intercept}
\item{A}{effect of A}
\item{B}{effect of B}
\item{AB}{effect of AB}
\item{w_00}{by-item random intercept (1 to 3)}
\item{w_11}{by-item random slope variance for A (1 to 3)}
\item{w_22}{by-item random slope variance for B (1 to 3)}
\item{w_33}{by-item random slope variance for AB (1 to 3)}
\item{w_01}{by-item random intercept/A-slope variance (-.9,.9)}
\item{w_02}{by-item random intercept/B-slope variance (-.9,.9)}
\item{w_03}{by-item random intercept/AB-slope variance (-.9,.9)}
\item{w_12}{by-item random A/B-slope variance (-.9,.9)}
\item{w_13}{by-item random A/AB-slope variance (-.9,.9)}
\item{w_23}{by-item random B/AB-slope variance (-.9,.9)}
}
\description{
Generates parameter ranges for two-way factorial design, with A
between-subjects and B within-subjects, and A,B within-items.
}
\seealso{
\code{\link{genParamRanges.facMixedAB}}, \code{\link{randParams}}, \code{\link{mkDf}}
}
