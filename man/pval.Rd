\name{pval}
\alias{pval}

\title{Example p-Values}

\description{
  This data set contains 4,289 p-values.  These data are used to
  illustrate the functionality of the functions \code{\link{fdr.control}}
  and \code{\link{fdr.estimate.eta0}}.
}

\usage{
data(pval)
}

\format{
  \code{pval} is a vector with 4,289 p-values. 
}


\examples{
# load fdrtool library
library("fdrtool")

# load data set
data(pval)

# plot distribution of p-values
density.pr(pval, df=30, plot=TRUE)
}
\keyword{datasets}

