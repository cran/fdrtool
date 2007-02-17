\name{halfnormal}
\alias{halfnormal}

\alias{dhalfnorm}
\alias{phalfnorm}
\alias{qhalfnorm}
\alias{rhalfnorm}
\alias{sd2theta}
\alias{theta2sd}


\title{The Half-Normal Distribution}
\description{
  Density, distribution function, quantile function and random
  generation for the half-normal distribution with parameter \code{theta}.
}
\usage{
dhalfnorm(x, theta=sqrt(pi/2), log = FALSE)
phalfnorm(q, theta=sqrt(pi/2), lower.tail = TRUE, log.p = FALSE)
qhalfnorm(p, theta=sqrt(pi/2), lower.tail = TRUE, log.p = FALSE)
rhalfnorm(n, theta=sqrt(pi/2))
sd2theta(sd)
theta2sd(theta)
}
\arguments{
  \item{x,q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations. If \code{length(n) > 1}, the length
    is taken to be the number required.}
  \item{theta}{parameter of half-normal distribution.}
  \item{log, log.p}{logical; if TRUE, probabilities p are given as log(p).}
  \item{lower.tail}{logical; if TRUE (default), probabilities are
    \eqn{P[X \le x]}{P[X <= x]}, otherwise, \eqn{P[X > x]}{P[X > x]}.}
  \item{sd}{standard deviation of the zero-mean normal distribution 
            that corresponds to the half-normal with parameter \code{theta}. }
}
\value{
  \code{dhalfnorm} gives the density,
  \code{phalfnorm} gives the distribution function,
  \code{qhalfnorm} gives the quantile function, and
  \code{rhalfnorm} generates random deviates.
  \code{sd2theta} computes a \code{theta} parameter. 
  \code{theta2sd} computes a \code{sd} parameter. 
}
\details{
 \code{x = abs(z)} follows a half-normal distribution with
  if \code{z} is a normal variate with zero mean. 
 The half-normal distribution has density
  \deqn{
    f(x) =
    \frac{2 \theta}{\pi} e^{-x^2 \theta^2/\pi}}{
    f(x) = 2*theta/pi e^-(x^2*theta^2/pi)
  }
  It has mean  \eqn{E(x) = \frac{1}{\theta}}{E(x) = 1/theta} and variance 
  \eqn{Var(x) = \frac{\pi-2}{2 \theta^2}}{Var(x) = (pi-2)/(2*theta^2)}.
  
 
  The parameter \eqn{\theta}{theta} is related to the
  standard deviation \eqn{\sigma}{sigma} of the corresponding 
  zero-mean normal distribution by the equation
  \eqn{\theta = \sqrt{\pi/2}/\sigma}{theta = sqrt(pi/2)/sigma}.

  If \eqn{\theta}{theta} is not specified in the above functions 
  it assumes the default values of \eqn{\sqrt{\pi/2}}{sqrt(pi/2)}, 
  corresponding to \eqn{\sigma=1}{sigma=1}.  
}
\seealso{
  \code{\link{Normal}}.
}

\examples{
# load "fdrtool" library
library("fdrtool")


## density of half-normal compared with a corresponding normal
par(mfrow=c(1,2))

sd.norm = 0.64
x  = seq(0, 5, 0.01)
x2 = seq(-5, 5, 0.01)
plot(x, dhalfnorm(x, sd2theta(sd.norm)), type="l", xlim=c(-5, 5), lwd=2,
   main="Probability Density", ylab="pdf(x)")
lines(x2, dnorm(x2, sd=sd.norm), col=8 )


plot(x, phalfnorm(x, sd2theta(sd.norm)), type="l", xlim=c(-5, 5),  lwd=2,
   main="Distribution Function", ylab="cdf(x)")
lines(x2, pnorm(x2, sd=sd.norm), col=8 )

legend("topleft", 
c("half-normal", "normal"), lwd=c(2,1),
col=c(1, 8), bty="n", cex=1.0)

par(mfrow=c(1,1))


## distribution function

integrate(dhalfnorm, 0, 1.4, theta = 1.234)
phalfnorm(1.4, theta = 1.234)

## quantile function
qhalfnorm(-1) # NaN
qhalfnorm(0)
qhalfnorm(.5)
qhalfnorm(1)
qhalfnorm(2) # NaN

## random numbers
theta = 0.72
hz = rhalfnorm(10000, theta)
hist(hz, freq=FALSE)
lines(x, dhalfnorm(x, theta))

mean(hz) 
1/theta  # theoretical mean

var(hz)
(pi-2)/(2*theta*theta) # theoretical variance


## relationship with two-sided normal p-values
z = rnorm(1000)

# two-sided p-values
pvl = 1- phalfnorm(abs(z))
pvl2 = 2 - 2*pnorm(abs(z)) 
sum(pvl-pvl2)^2 # equivalent
hist(pvl2, freq=FALSE)  # uniform distribution

# back to half-normal scores
hz = qhalfnorm(1-pvl)
hist(hz, freq=FALSE)
lines(x, dhalfnorm(x))
}
\keyword{distribution}
