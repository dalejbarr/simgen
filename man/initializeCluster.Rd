\name{initializeCluster}
\alias{initializeCluster}
\title{Load in package simgen on all available clusters}
\usage{
  initializeCluster(cl)
}
\arguments{
  \item{cl}{cluster object (created e.g. by
  \code{\link{makeCluster}}; see package parallel for
  details.}
}
\description{
  This function is called by \code{\link{mcRun}} to load
  the package functions across all clusters.  Typically
  users should not need to call this function.
}

