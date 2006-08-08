### density.pr.R  (2005-06-05)
###
###    Density estimation via Poisson regression
###
### Copyright 2005 Korbinian Strimmer
###
### This functions is based in part on code from the "locfdr" package 
###
###
### This file is part of the `GeneTS' library for R and related languages.
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


# Reference:
#
# B. Efron and R. Tibshirani. 1996. Using specially designed exponential 
# families for density estimation.  Annals of Statistics 24:2431-2461


density.pr <- function(x, ncells=100,
    trim.lo=1/1000, trim.up=1/1000, df=7, plot=FALSE, ...)
{     
    # original data name and length
    name <- deparse(substitute(x))
    l.x <- length(x)
    
    # remove missing and infinite values
    has.na <- any(is.na(x))
    x <- x[is.finite(x)]
  
    # determine lower and upper boundaries of histogram
    # (i.e remove outliers that may lead to convergence problems)
    v <- quantile(x, c(trim.lo, 1 - trim.up))
    lo <- v[1]
    up <- v[2]
    x <- x[x > lo & x < up]
       
    # Note: in the original locfdr procedure the outliers are all put into
    # the outer cells:
    #       x <- pmax(pmin(x, up), lo)
    #
    # Later, just before performing the Poisson regression the inflated
    # counts at the borders are then corrected via
    #      if (trim.lo > 0) yy[1] <- min(yy[1], 1)
    #      if (trim.up > 0) yy[ncells] <- min(yy[ncells], 1)    
    #
    # This works well *if* the density is vanishing near the borders.
    # However, if this is not the case (and there are many counter examples!)
    # then setting yy[] to 0 or 1 is not really a good idea.  Instead, it 
    # generally seems safer to simply eliminate the outliers.
 
    # determine cell breaks and cell (bin) width
    breaks <- seq(lo, up, length=(ncells+1))
    bw <- (up-lo)/ncells
          
    # determine cell counts
    hist.x <- hist(x, breaks = breaks, plot = FALSE)
    xx <- hist.x$mids
    yy <- hist.x$counts
    
    # scale factor (= area under the histogram bars)
    # this is used turn f below into a density
    scale <- sum(diff(breaks)*yy)    
            
    # fit Poisson model using a natural spline
    require("splines")
    f <- glm(yy ~ ns(xx, df = df), poisson)$fit
    
    # some other possibilities
    #f <- glm(yy ~ poly(xx, df = df), poisson)$fit
    #f <- glm(yy ~ bs(xx, df = df), poisson)$fit

  
    # check fit of Poisson model 
    SSR <- (yy - f)^2/(f + 1) # squared studentized residuals
    Phi <- sum(SSR)/(ncells - df) # estimated dispersion
    if(Phi > 1.5)
    {
      warning(paste("Estimated dispersion =", round(Phi, 2),
       "- increase df to improve fit of Poisson model?"))
    }

            
    # density object
    dx <- list(
              x = xx,
              y = f/scale,
	      bw = bw,
	      n = l.x,
	      dispersion = Phi,
	      call = match.call(),
	      data.name = name,
	      has.na = has.na,
	      histogram = hist.x # return histogram as well
              )
    class(dx) <- "density"
    
    # if desired plot histogram and estimated density    
    if (plot)
    {
      plot(hist.x, freq=FALSE, ...)
      lines(dx, lwd = 3, col = 4)
    }
    
    return(dx)
}

