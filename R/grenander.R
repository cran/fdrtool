### grenander.R  (2007-02-14)
###
###     Grenander Density Estimator
###
### Copyright 2006 Korbinian Strimmer 
###
###
### This file is part of the `fdrtool' library for R and related languages.
### It is made available under the terms of the GNU General Public
### License, version 2, or at your option, any later version,
### incorporated herein by reference.
### 
### This program is distributed in the hope that it will be
### useful, but WITHOUT ANY WARRANTY; without even the implied
### warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
### PURPOSE.  See the GNU General Public License for more
### details.
### 
### You should have received a copy of the GNU General Public
### License along with this program; if not, write to the Free
### Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
### MA 02111-1307, USA


grenander = function(F)
{
  if( !any(class(F) == "ecdf") ) stop("ecdf object required as input!")

  # find least concave majorant of ECDF
  ll = gcmlcm(environment(F)$x, environment(F)$y, type="lcm")

  f.knots = ll$slope.knots
  f.knots = c(f.knots, f.knots[length(f.knots)])

  g = list(F=F,
       x.knots=ll$x.knots,
       F.knots=ll$y.knots,
       f.knots=f.knots)

  class(g) <- "grenander"
  
  return(g)
}

plot.grenander <- function(x, ...)
{
  par(mfrow=c(1,2))

  plot(x$x.knots, x$f.knots, type="s", xlab="x", ylab="fn(x)",
     main="Grenander Decreasing Density", col=4, lwd=2, log="y", ...)
 
  plot(x$F, do.points=FALSE)
  lines(x$x.knots, x$F.knots, type='l', col=4, lwd=2)


  par(mfrow=c(1,1))
}



