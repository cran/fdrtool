\name{censored.fit}
\alias{censored.fit}

\title{Fit Null Distribution To Censored Data}

\usage{
censored.fit(x, statistic=c("normal", "correlation", "studentt"), 
  pct0=3/4,  method=c("MM", "ML"))
}

\description{
  \code{censored.fit} fits the scale parameter of a null distribution with mean zero
  to censored data.
}

\arguments{
  \item{x}{vector of test statistics.}
  \item{statistic}{type of statistic - normal, correlation, or student t.}
  \item{pct0}{fraction of data used to fit null distribution.}
  \item{method}{procedure used to fit null distibution, either of MM (=median matching, the default)
        or ML (maximum likelihood).}
}
\details{
   For the default method (median matching), the median of the truncated 
   half-normal distribution (or the equivalent for correlation and student t distribution) 
   is matched to the median of the absolute values of the censored sample.

   For maximum likelihood the truncated density is fit to the censored    
   sample - this can be quite time consuming for large data set (and there
   may be convergence problems).
}
\value{

  The estimate scale parameter of the null distribution. This depends
  on the type of statistic:
  \item{normal}{sd value.} 
  \item{correlation}{kappa value (see \code{\link{dcor0}}).}
  \item{studentt}{degrees of freedom (df).}
}


\seealso{\code{\link{dcor0}}.}

\examples{
# load "fdrtool" library
library("fdrtool")

# simulate normal data
sd.true = 2.232
n = 5000
z = rnorm(n, sd=sd.true)
censored.fit(z, statistic="normal")


# simulate contaminated mixture of correlation distribution
r <- rcor0(700, kappa=10)
u1 <- runif(200, min=-1, max=-0.7)
u2 <- runif(200, min=0.7, max=1)
rc <- c(r, u1, u2)

censored.fit(r, statistic="correlation")
censored.fit(rc, statistic="correlation")
}
\keyword{htest}
