\name{monoreg}
\alias{monoreg}
\alias{plot.monoreg}
\alias{fitted.monoreg}
\alias{residuals.monoreg}

\title{Monotone Regression: Isotonic Regression and Antitonic Regression}

\usage{
monoreg(x, y=NULL, w=rep(1, length(x)), type=c("isotonic", "antitonic"))
}

\description{
  \code{monoreg} performs monotone regression (either isotonic
  or antitonic) with weights. 
}

\arguments{
  \item{x, y}{coordinate vectors of the regression
    points.  Alternatively a single \dQuote{plotting} structure can be
    specified: see \code{\link{xy.coords}}.}
  \item{w}{data weights (default values: 1).}
    
 
  \item{type}{fit a monotonely increasing ("isotonic") or 
  monotonely decreasing ("antitonic") function.}
}

\details{

 \code{monoreg} is similar to \code{\link{isoreg}}, with the addition
 that \code{monoreg} accepts weights.

 If several identical \code{x} values are given as input, the 
 corresponding \code{y} values and the
 weights \code{w} are automatically merged, and a warning is issued.
 
 The \code{plot.monoreg} function optionally plots the cumulative
 sum diagram with the greatest convex minorant (isotonic regression)
 or the least concave majorant (antitonic regression), see the
 examples below. 
}

\value{
  A list with the following entries:
  
  \item{x}{the sorted and unique x values}
  \item{y}{the corresponding y values}
  \item{w}{the corresponding weights}
  \item{yf}{the fitted y values}
  \item{type}{the type of monotone regression ("isotonic" or "antitonic"}
  \item{call}{the function call}
}

\references{
 
 Robertson, T., F. T. Wright, and R. L. Dykstra. 1988.  Order restricted
 statistical inference. John Wiley and Sons.

}
\author{
  Korbinian Strimmer (\url{http://strimmerlab.org}).
  
  Part of this function is C code that has been ported from
  R code originally written by Kaspar Rufibach. 
}
\seealso{\code{\link{isoreg}}.}

\examples{
# load "fdrtool" library
library("fdrtool")


# an example with weights

# Example 1.1.1. (dental study) from Robertson, Wright and Dykstra (1988)
age = c(14, 14, 8, 8, 8, 10, 10, 10, 12, 12, 12)
size = c(23.5, 25, 21, 23.5, 23, 24, 21, 25, 21.5, 22, 19)

mr = monoreg(age, size)

# sorted x values
mr$x # 8 10 12 14
# weights and merged y values
mr$w  # 3 3 3 2
mr$y #  22.50000 23.33333 20.83333 24.25000
# fitted y values
mr$yf # 22.22222 22.22222 22.22222 24.25000
fitted(mr)
residuals(mr)

plot(mr, ylim=c(18, 26))  # this shows the averaged data points
points(age, size, pch=2)  # add original data points


###

y = c(1,0,1,0,0,1,0,1,1,0,1,0)
x = 1:length(y)
mr = monoreg(y)

# plot with greatest convex minorant
plot(mr, plot.type="row.wise")  

# this is the same
mr = monoreg(x,y)
plot(mr)

# antitonic regression and least concave majorant
mr = monoreg(-y, type="a")
plot(mr, plot.type="row.wise")  

# the fit yf is independent of the location of x and y
plot(monoreg(x + runif(1, -1000, 1000), 
             y +runif(1, -1000, 1000)) )

###

y = c(0,0,2/4,1/5,2/4,1/2,4/5,5/8,7/11,10/11)
x = c(5,9,13,18,22,24,29,109,120,131)

mr = monoreg(x,y)
plot(mr, plot.type="row.wise")

# the fit (yf) only depends on the ordering of x
monoreg(1:length(y), y)$yf 
monoreg(x, y)$yf 


}
\keyword{regression}
\keyword{smooth}
