### censored.fit2.R  (2007-02-07)
###
###     Fit Null Distribution Using To Data
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


# private function for MM method (median matching)


# estimate scale parameter of null component: 
# fit truncated density to censored sample by matching empirical
# and theoretical medians

pvt.censored.fit2 <- function(x, x0,
   statistic=c("normal", "correlation", "studentt") )
{

######################
# perhaps better to use this directly as argument

  ax = abs(x)
  vx = var(x)
######################


  x0 = abs(x0)
  ax.cens <- ax[ ax <= x0 ]

  if (length(ax.cens) < 10) 
    warning(paste("Adjust threshold - censored sample has only size", 
      length(ax.cens), "!"), call.=FALSE)

  if (length(ax.cens) < 2) 
    stop(paste("Adjust threshold - censored sample has only size", 
      length(ax.cens), "!"), call.=FALSE)

  ################
 
  if (statistic=="normal")
  {  
    # theoretical median for given truncation point and scale parameter
    m.theory = function(sd)
    {
      #th = sd2theta(sd)
      #pp = phalfnorm(x0, theta=th)
      #qq = qhalfnorm(1/2*pp, theta=th)

      pp = 2*pnorm(x0, sd=sd)-1
      qq = qnorm((1/2*pp+1)/2,  sd=sd)

      return(qq)
    }
  }

  if (statistic=="correlation")
  {  
    # theoretical median for given truncation point and scale parameter
    m.theory = function(sd)
    {
      kappa = 1/(sd*sd)
      
      pp = 2*pcor0(x0, kappa=kappa)-1
      qq = qcor0((1/2*pp+1)/2,  kappa=kappa)

      return(qq)
    }
  }

  if (statistic=="studentt")
  {  
    # theoretical median for given truncation point and scale parameter
    m.theory = function(sd)
    {
      if (sd < 1) stop("sd must be >= 1")
      v = sd*sd 
      df = 2*v/(v-1)
      
      pp = 2*pt(x0, df=df)-1
      qq = qt((1/2*pp+1)/2,  df=df)

      return(qq)
    }
  }


  ###############


  # median of truncated sample
  m.observed = median(ax.cens)  
  minp = 0
  if (statistic=="studentt") minp=1.001002 # df = 1000
  maxp = 2*sqrt(vx) 

  #cat("DEBUG: minp=", minp, ",  maxp=", maxp, "\n")

  #f <- function(y) m.theory(y)-m.observed
  #param <- uniroot(f, lower=minp, upper=maxp)$root

  f <- function(y) (m.theory(y)-m.observed)^2
  param <- optimize(f, lower=minp, upper=maxp)$minimum

  # param -> sd !

  return(
    list(param=param) 
  )
}

