\name{density.pr}
\alias{density.pr}
\title{Density Estimation Via Poisson Regression}
\description{
  The function \code{density.pr} estimates a one-dimensional probability density distribution 
  by a Poisson regression fit to histogram counts. 
}
\usage{
density.pr(x, ncells=100, trim.lo=1/1000, trim.up=1/1000, df=7, plot=FALSE, \dots)
}
\arguments{
  \item{x}{the raw input data (missing and infinite values will be automatically removed)}
  \item{ncells}{number of cells of the underlying histogram}
  \item{trim.lo}{fraction of data points to be removed from the lower tail of the
                 distribution (this prevents convergence problems)}
  \item{trim.up}{fraction of data points to be removed from the upper tail of the
                 distribution (this prevents convergence problems)}
  \item{df}{degree of freedom of the natural spline used as regression function}		 
  \item{plot}{plot estimated density and the underlying histogram}
  \item{\dots}{arguments passed to the plot function if \code{plot=TRUE}}
}
\details{
  Density estimation using Poisson regression is advocated in Efron and Tibshirani (1996),
  and these authors trace the origins of this method back to Lindsey (1974a, b).
  
  The algorithm in \code{density.pr} proceeds in two steps. First, a histogram is constructed 
  as prelimary density estimate. Subsequently, a Poisson regression is employed to fit
  the histogram counts, where a natural spline is used as regression function. 
}
\value{
  A list of class \code{density} with the following component:
  
  \item{x}{ the coordinates of the points where the density is estimated.   } 
  \item{y}{ the estimated density values.   }
  \item{bw}{ the bin width used in the underlying histogram.
  }
  \item{n}{  the original sample size (before elimination of missing and infinite values).   }
 \item{dispersion}{  estimated dispersion (a warning is issued if dispersion > 1.5 - this can usually
                      be fixed by increasing the \code{df} parameter).   }  
  \item{call}{ the call which produced the result.}
  \item{data.name}{ the deparsed name of the input data.   }
  \item{has.na}{ \code{TRUE} if input data contains missing data or infinite values.}
  \item{histogram}{ the underlying histogram.   }

}
\author{
  Korbinian Strimmer (\url{http://strimmerlab.org}).
  
  Part of the code of the \code{density.pr} function was adopted from the \code{locfdr} package.
}
\references{

  Efron, B., and Tibshirani, R. (1996). Using specially designed exponential 
  families for density estimation.  \emph{Annals of Statistics},
  \bold{24}, 2431--2461.
  
  Lindsey, J. (1974a). Comparison of probability distributions.
  \emph{JRSS B}, \bold{36}, 38-47.
  
   Lindsey, J. (1974b). Construction and comparison of statistical models.
  \emph{JRSS B}, \bold{36}, 418-425.
 
}
\seealso{\code{\link{density}}.}

\examples{
# load fdrtools library
library("fdrtool")

# load data 
data(faithful)
z <- faithful[,1]
rz <- range(z)

# estimate density 
d1a <- density.pr(z, plot=TRUE)
d1a

# discretization is not critical in this algorithm
d1b <- density.pr(z, plot=TRUE, ncells=80)
d1c <- density.pr(z, plot=TRUE, ncells=40)

plot(d1a, col=2)
lines(d1b, col=2)
lines(d1c, col=2)

# comparison with kernel density estimate
d2a <- density(z, from=rz[1], to=rz[2])
d2b <- density(z, from=rz[1], to=rz[2], bw="bcv")

plot(d1a, col=2, xlim=rz)
lines(d2a)
lines(d2b)
}
\keyword{univar}
