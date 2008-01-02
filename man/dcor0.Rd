\name{dcor0}
\alias{dcor0}
\alias{pcor0}
\alias{rcor0}
\alias{qcor0}

\title{Distribution of the Vanishing Correlation Coefficient (rho=0) and Related Functions}

\usage{
dcor0(x, kappa, log=FALSE)
pcor0(q, kappa, lower.tail=TRUE, log.p=FALSE)
qcor0(p, kappa, lower.tail=TRUE, log.p=FALSE)
rcor0(n, kappa)
}

\arguments{
	\item{x,q}{vector of sample correlations}
        \item{p}{vector of probabilities}
	\item{kappa}{the degree of freedom of the distribution (= inverse variance)}
	\item{n}{number of values to generate. If n is a vector, length(n)
                values will be generated}
	\item{log, log.p}{logical vector; if TRUE, probabilities p are given as log(p)}
	\item{lower.tail}{logical vector; if TRUE (default), probabilities are \eqn{P[R <= r]},
	                  otherwise, \eqn{P[R > r]}}
}

\description{
  The above functions describe the distribution of the Pearson correlation 
  coefficient \code{r} assuming that there is no correlation present (\code{rho = 0}).
 
  Note that the distribution has only a single parameter: the degree 
  of freedom \code{kappa}, which is equal to the inverse of the variance of the distribution.
  
  The theoretical value of  \code{kappa} depends both on the sample size \code{n} and the number 
  \code{p} of considered variables.  If a simple correlation coefficient between two
  variables  (\code{p=2})  is considered the degree of freedom equals \code{kappa = n-1}.
  However, if a partial correlation coefficient is considered (conditioned on \code{p-2} remaining 
  variables) the degree of freedom is \code{kappa = n-1-(p-2) = n-p+1}.
}

\details{
  For density and distribution functions as well as a corresponding random number generator
  of the correlation coefficient for arbitrary non-vanishing correlation \code{rho} please refer to the
  \code{SuppDists} package by  Bob Wheeler \email{bwheeler@echip.com} (available on CRAN).
  Note that the parameter \code{N} in his \code{dPearson} function corresponds to  \code{N=kappa+1}.
}

\value{
   \code{dcor0} gives the density,  \code{pcor0} 
  gives the distribution function, \code{qcor0} gives
  the quantile function, and \code{rcor0} generates random deviates.
}

\seealso{\code{\link{cor}}.}

\author{
   Korbinian Strimmer (\url{http://strimmerlab.org}).
}

\examples{
# load fdrtool library
library("fdrtool")

# distribution of r for various degrees of freedom
x = seq(-1,1,0.01)
y1 = dcor0(x, kappa=7)
y2 = dcor0(x, kappa=15)
plot(x,y2,type="l", xlab="r", ylab="pdf",
  xlim=c(-1,1), ylim=c(0,2))
lines(x,y1)

# simulated data
r = rcor0(1000, kappa=7)
hist(r, freq=FALSE, 
  xlim=c(-1,1), ylim=c(0,5))
lines(x,y1,type="l")

# distribution function
pcor0(-0.2, kappa=15)
}
\keyword{distribution}
