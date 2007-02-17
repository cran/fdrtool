\name{gcmlcm}
\alias{gcmlcm}

\title{Greatest Convex Minorant and Least Concave Majorant}

\usage{
gcmlcm(x, y, type=c("gcm", "lcm"))
}

\description{
  \code{gcmlcm} computes the greatest convex minorant (GCM) or the
  least concave majorant (LCM) of a piece-wise linear function.
}

\arguments{
  \item{x, y}{coordinate vectors of the piece-wise linear function. Note
   that the x values need to be unique and be arranged in sorted order.}
  \item{type}{specifies whether to compute the greatest convex 
        minorant (\code{type="gcm"}, the default) or the
        least concave majorant (\code{type="lcm"}).}
}

\details{

 The GCM is obtained by isotonic regression of the raw slopes,
 whereas the LCM is obtained by antitonic regression.
 See Robertson et al. (1988).
}

\value{
  A list with the following entries:
  
  \item{x.knots}{the x values belonging to the knots of the LCM/GCM curve}
  \item{y.knots}{the corresponding y values}
  \item{slope.knots}{the slopes of the corresponding line segments}
}

\references{
 
 Robertson, T., F. T. Wright, and R. L. Dykstra. 1988.  Order restricted
 statistical inference. John Wiley and Sons.

}
\author{
  Korbinian Strimmer (\url{http://strimmerlab.org}).
}
\seealso{\code{\link{monoreg}}.}

\examples{
# load "fdrtool" library
library("fdrtool")

# generate some data
x = 1:20
y = rexp(20)
plot(x, y, type="l", lty=3, main="GCM (red) and LCM (blue)")
points(x, y)

# greatest convex minorant (red)
gg = gcmlcm(x,y)
lines(gg$x.knots, gg$y.knots, col=2, lwd=2)

# least concave majorant (blue)
ll = gcmlcm(x,y, type="lcm")
lines(ll$x.knots, ll$y.knots, col=4, lwd=2)

}
\keyword{smooth}
