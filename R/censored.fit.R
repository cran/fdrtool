### censored.fit.R  (2008-02-25)
###
###     Fit Null Distribution To Censored Data by Maximum Likelihood
###
### Copyright 2006-08 Korbinian Strimmer 
###
###
### This file is part of the `fdrtool' library for R and related languages.
### It is made available under the terms of the GNU General Public
### License, version 3, or at your option, any later version,
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


# estimate parameters of null distribution 
# using using truncated distributions
# available null distributions
# - normal (with mean zero)
# - correlation (with rho zero)
# - student t
# - uniform


censored.fit = function(x, cutoff,
   statistic=c("normal", "correlation", "pvalue", "studentt"))
{
    statistic = match.arg(statistic)
    cutoff = abs(cutoff)

    if ( !is.vector(x) ) stop("x needs to be a vector!")
    
    
    if (statistic=="pvalue")
    {
      result = matrix(nrow=length(cutoff), ncol=4)
      colnames(result)= c("cutoff", "N0", "eta0", "eta0.SE")
    }
    else
    {
      result = matrix(nrow=length(cutoff), ncol=6)
      colnames(result)= c("cutoff", "N0", "eta0", "eta0.SE", "sd", "sd.SE")
    }
    if (statistic=="correlation") colnames(result)[5] = "kappa"
    if (statistic=="studentt") colnames(result)[5] = "df"
    if (statistic=="correlation") colnames(result)[6] = "kappa.SE"
    if (statistic=="studentt") colnames(result)[6] = "df.SE"

    
    for (i in 1:length(cutoff))
    {
      x0 = cutoff[i]
      result[i,1] = x0

      out = pvt.fit.nullmodel(x, x0, statistic=statistic)
      result[i,2] = out$N0
      result[i,3] = out$eta0
      result[i,4] = out$eta0.SE

      if (statistic!="pvalue")
      {
         result[i,5] = out$param
         result[i,6] = out$param.SE
      }
    }

    return(result)
}

### helper functions


# Richardson extrapolation approximation 
# for numerical computation of curvature
num.curv = function(x, fun) 
{
  macheps = .Machine$double.eps
  h = max( 1e-4, macheps^(1/4)*abs(x) )

  w = c(-1/12,4/3,-5/2,4/3,-1/12) 
  xd = x + h*c(-2,-1,0,1,2)
 
  return( sum(w*fun(xd))/h^2 )
}



### internal functions 

pvt.fit.nullmodel = function(x, x0, statistic)
{
  N = length(x)

  if (statistic=="pvalue")
    x.cens = x[ x >= x0 ]
  else
    x.cens = x[ abs(x) <= x0 ] 

  N0 = length(x.cens) 
  if (N0 > N) stop("N must be larger or equal to the size of the censored sample!")
  if (N0 < 10) 
    warning(paste("Adjust threshold - censored sample has only size", 
      length(x.cens), "!"), call.=FALSE)

  if (N0 < 2) 
    stop(paste("Adjust threshold - censored sample has only size", 
      length(x.cens), "!"), call.=FALSE)


  ##############

  nm = get.nullmodel(statistic)

  
  # negative log-likelihood function (truncated density)
  nlogL = function(pp) 
  {
    out = rep(0, length(pp))
    for (i in 1:length(pp))
    {
      out[i] = length(x.cens)*log(1-nm$get.pval(x0, pp[i]))-
               sum(nm$f0(x.cens, pp[i], log=TRUE))
    }
    return(out)
  } 

  ##############
  
  # estimate parameters of null model
  if (statistic!="pvalue")
  {
    start = iqr.fit(x.cens) # start value for scale parameter
    
    #sup = nm$get.support()
    #opt.out = nlminb( start, nlogL, lower=sup[1], upper=sup[2] )
    #sc.param = opt.out$par[1]
   
    lo = start/1000
    up = start*1000
    sc.param = optimize(nlogL, lower=lo, upper=up)$minimum

    sc.var = 1/num.curv(sc.param,nlogL) # inverse curvature of negative logL
    if(sc.var < 0) 
    {
       sc.var = 0
       warning("Variance of scale parameter set to zero due to numerical problems")
    }
    sc.SE = sqrt(sc.var)
  }
  else
  {
    sc.param = NULL # no scale parameter
    sc.SE = NULL    
  }

  # ML estimate of eta0
  m = 1-nm$get.pval(x0, sc.param)
  th = N0/N
  eta0 = min(1, th / m )
  #eta0 = th / m 
  eta0.SE = sqrt( th*(1-th)/(N*m*m) )

  rm(x.cens)

  return(
    list(N0=N0,
         eta0=eta0,              # proportion
         eta0.SE=eta0.SE,        # corresponding standard error
         param=sc.param,         # scale parameter
         param.SE=sc.SE          # corresponding standard error
        )
  )
}

