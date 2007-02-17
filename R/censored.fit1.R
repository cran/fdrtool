### censored.fit1.R  (2007-02-07)
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


# private function for ML method (maximum likelihood fit)


# estimate scale parameter of null component: 
# fit truncated density to censored sample

pvt.censored.fit1 <- function(x, x0, 
   statistic=c("normal", "correlation", "studentt") )
{
  x0 = abs(x0)
  x.cens <- x[ abs(x) <= x0 ]

  vx = var(x)

  if (length(x.cens) < 10) 
    warning(paste("Adjust threshold - censored sample has only size", 
      length(x.cens), "!"), call.=FALSE)

  if (length(x.cens) < 2) 
    stop(paste("Adjust threshold - censored sample has only size", 
      length(x.cens), "!"), call.=FALSE)


  if (statistic=="studentt")
  {
    # log-density of truncated t  
    logDens <- function(xc, x0, sd)
    {
       if (sd < 1) stop("sd must be >= 1")
       v = sd*sd 
       df = 2*v/(v-1)

       # note: argument xc must be from censored sample,
       #       otherwise d=-Inf !
       m <- log( 2*pt(x0, df=df)-1 )
       d <- dt(xc, df=df, log=TRUE)-m

       return(d)
    }
   }


  if (statistic=="normal")
  {
    # log-density of truncated normal  
    logDens <- function(xc, x0, sd)
    {
       # note: argument xc must be from censored sample,
       #       otherwise d=-Inf !
       m <- log( 2*pnorm(x0, sd=sd)-1 )
       d <- dnorm(xc, sd=sd, log=TRUE)-m

       return(d)
    }
  }

  if (statistic=="correlation")
  {
    # log-density of truncated normal  
    logDens <- function(xc, x0, sd)
    {
       kappa = 1/(sd*sd)

       # note: argument xc must be from censored sample,
       #       otherwise d=-Inf !
       m <- log( 2*pcor0(x0, kappa=kappa)-1 )
       d <- dcor0(xc, kappa=kappa, log=TRUE)-m

       return(d)
    }
  }

  # find ML estimate of sd parameter
  minp = 0
  if (statistic=="studentt") minp=1.001002 # df = 1000
  maxp = 2*sqrt(vx) 

  # log-likelihood function
  logL <- function(sd)
  {
    sum(logDens(x.cens, x0, sd))
  }

  opt.out <- optimize(logL, lower=minp, upper=maxp, maximum=TRUE)


  return(
    list(logL=opt.out$objective, 
         param=opt.out$maximum) 
  )
}
