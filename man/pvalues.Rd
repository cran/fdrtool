\name{pvalues}
\alias{pvalues}

\title{Example p-Values}

\description{
  This data set contains 4,289 p-values.  These data are used to
  illustrate the functionality of the functions \code{\link{fdrtool}}
  and \code{\link{pval.estimate.eta0}}.
}

\usage{
data(pvalues)
}

\format{
  \code{pvalues} is a vector with 4,289 p-values. 
}


\examples{
# load fdrtool library
library("fdrtool")

# load data set
data(pvalues)

# estimate density and distribution function, 
# and compute corresponding (local) false discovery rates
fdrtool(pvalues, statistic="pvalue")
}
\keyword{datasets}

