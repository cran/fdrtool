\name{censored.fit}
\alias{censored.fit}
\alias{fndr.cutoff}

\title{Fit Null Distribution To Censored Data by Maximum Likelihood}

\usage{
censored.fit(x, cutoff, statistic=c("normal", "correlation", "pvalue", "studentt"))
fndr.cutoff(x, statistic=c("normal", "correlation", "pvalue", "studentt"))
}

\description{
  \code{censored.fit} fits  a null distribution 
  to censored data.

  \code{fndr.cutoff} finds a suitable cutoff point based on the 
   (approximate) false non-discovery rate (FNDR).
}

\arguments{
  \item{x}{vector of test statistics.}
  \item{cutoff}{truncation point (this may a single value or a vector).}
  \item{statistic}{type of statistic - normal, correlation, or student t.}
}
\details{

   As null model truncated normal, truncated student t or a truncated
   correlation density is assumed.  The truncation point is specified
   by the cutoff parameter.  All data points whose absolute value
   are large than the cutoff point are ignored when fitting the truncated
   null model via maximum likelihood.  The total number of data points is
   only used to estimate the fraction of null values eta0.
}
\value{
  \code{censored.fit} returns a matrix whose rows contain the estimated parameters and corresponding errors
  for each cutoff point. 

  \code{fndr.cutoff} returns a tentative cutoff point.
}


\seealso{\code{\link{fdrtool}}.}

\examples{
# load "fdrtool" library
library("fdrtool")

# simulate normal data
sd.true = 2.232
n = 5000
z = rnorm(n, sd=sd.true)
censored.fit(z, c(2,3,5), statistic="normal")


# simulate contaminated mixture of correlation distribution
r = rcor0(700, kappa=10)
u1 = runif(200, min=-1, max=-0.7)
u2 = runif(200, min=0.7, max=1)
rc = c(r, u1, u2)

censored.fit(r, 0.7, statistic="correlation")
censored.fit(rc, 0.7, statistic="correlation")

# pvalue example
data(pvalues)
co = fndr.cutoff(pvalues, statistic="pvalue")
co
censored.fit(pvalues, cutoff=co, statistic="pvalue")
}
\keyword{htest}
