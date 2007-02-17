### censored.fit.R  (2007-02-07)
###
###     Fit Null Distribution To Censored Data
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


# estimate parameters of null distribution using censored data

# available null distributions
# - normal (with mean zero)
# - correlation (with rho zero)
# - student t


# Methods:
# - maximum likelihood (ML) fit 
# - MM median matching  
# truncated density to censored sample


# pct0 = seq(0.7, 1.0, 0.02)
# pct0 = 3/4


censored.fit <- function(x, 
   statistic=c("normal", "correlation", "studentt"),
   pct0=3/4,  method=c("MM", "ML"))
{
    statistic <- match.arg(statistic)
    method <- match.arg(method)

    if ( !is.vector(x) ) stop("x needs to be a vector!")
    if ( length(x) < 100 ) warning("estimates may be unreliable as length(x) = ", length(x))

    diagnostic.plot=TRUE  # only if smooth = TRUE

    if (length(pct0) > 1) smooth=TRUE 
    else smooth=FALSE

    if (smooth)
    {
      x0 = quantile(abs(x), probs=pct0)
      sc.vec = rep(NA, length(x0))
      for(i in 1:length(x0))
      {
        if (method=="ML")
          sc.vec[i] <- pvt.censored.fit1(x, x0[i], statistic=statistic)$param
        else # method=="MM"
          sc.vec[i] <- pvt.censored.fit2(x, x0[i], statistic=statistic)$param
      }
      sc.spline <- smooth.spline(pct0, sc.vec,df=3)
      sc.param <- predict(sc.spline, x=min(pct0))$y
      
      if (diagnostic.plot)
      {
        get(getOption("device"))() # open new plot window
        plot(pct0, sc.vec, main="Smoothing Curve Employed For Estimating scale parameter", 
          xlab="pct0", ylab="estimated sd")
        lines( sc.spline )
        points( min(pct0), sc.param, pch=20, col=2 )
      }
    }
    else # length(pct0) = 1
    {
      x0 = quantile(abs(x), probs=pct0)
      if (method=="ML")
        sc.param <- pvt.censored.fit1(x, x0, statistic=statistic)$param
      else # method=="MM"
        sc.param <- pvt.censored.fit2(x, x0, statistic=statistic)$param
    }


    ## the optimization is done on the level of "sd" parameter
    ## convert back to natural parameter

    if (statistic=="normal")
    {  
      attr(sc.param, "names") <- "sd"
    }

    if (statistic=="correlation")
    {  
      sc.param = 1/sc.param^2    # kappa
      attr(sc.param, "names") <- "kappa"
    }

    if (statistic=="studentt")
    {  
      sc.param = sc.param*sc.param  # var
      sc.param = 2*sc.param/(sc.param-1)  # df
    
      attr(sc.param, "names") <- "df"
    }

    return(sc.param)
}



