\name{grenader}
\alias{grenander}
\alias{plot.grenander}

\title{Grenander Estimator of a Decreasing or Increasing Density}
\description{
  The function \code{grenander} computes the Grenander estimator
  of a one-dimensional decreasing or increasing density. 
}
\usage{
grenander(F, type=c("decreasing", "increasing"))
}
\arguments{
  \item{F}{an \code{\link{ecdf}} containing the empirical cumulative density.}
  \item{type}{specifies whether the distribution is decreasing (the default) or increasing.}
}
\details{
  The Grenander (1956) density estimator is given by the slopes of the 
  least concave majorant (LCM) of the empirical distribution function (ECDF).
  It is a decreasing piecewise-constant function and can be shown to be the 
  non-parametric maximum likelihood estimate (NPMLE) under the assumption
  of a decreasing density (note that the ECDF is the NPMLE without
  this assumption).  Similarly, an increasing density function is obtained
  by using the greatest convex minorant (GCM) of the ECDF.
}
\value{
  A list of class \code{grenander} with the following components:
  
  \item{F}{the empirical distribution function specified as input.   } 
  \item{x.knots}{x locations of the knots of the least concave majorant of the ECDF.}
  \item{F.knots}{the corresponding y locations of the least concave majorant of the ECDF.}
  \item{f.knots}{the corresponding slopes (=density).}

}
\author{
  Korbinian Strimmer (\url{http://strimmerlab.org}).
}
\references{

  Grenander, U. (1956). On the theory of mortality measurement, Part II.  
  \emph{Skan. Aktuarietidskr},
  \bold{39}, 125--153.
}
\seealso{\code{\link{ecdf}}, \code{\link{gcmlcm}}, \code{\link{density}}.}

\examples{
# load "fdrtool" library
library("fdrtool")

# samples from random exponential variable 
z = rexp(30,1)
e = ecdf(z)
g = grenander(e)
g

# plot ecdf, concave cdf, and Grenander density estimator (on log scale)
plot(g, log="y") 

# for comparison the kernel density estimate
plot(density(z)) 

# area under the Grenander density estimator 
sum( g$f.knots[-length(g$f.knots)]*diff(g$x.knots) )
}
\keyword{univar}
