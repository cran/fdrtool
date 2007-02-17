### fdrtool.R  (2007-02-17)
###
###    Estimate (Local) False Discovery Rates For Diverse Test Statistics
###    
###
### Copyright 2007 Korbinian Strimmer
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



 
fdrtool <- function(x, 
  statistic=c("normal", "correlation", "pvalue", "studentt"),
  plot=TRUE, censored.fit.arg=NULL, pval.estimate.eta0.arg=NULL)
{
  statistic <- match.arg(statistic)
  ax = abs(x) 

  if (statistic=="pvalue") ax = 1-ax  # reverse p-val plot 

  
#### step 1 ####

  cat("Step 1... fit null distribution\n")

  # determine parameters of null distribution
  if (statistic != "pvalue")
  {
    cf.param <- do.call("censored.fit", c(list(x=x, statistic=statistic),
                  censored.fit.arg ) ) 
  }

#### step 2 ####

  cat("Step 2... compute p-values and estimate eta0\n")

  nf <- pvt.nullfunction(statistic=statistic, cf.param=cf.param)
  pval = 1- nf$F0(ax)

  # determine eta0
  eta0 <- do.call("pval.estimate.eta0", c(list(p=pval, diagnostic.plot=FALSE),
             pval.estimate.eta0.arg ) ) 
 
#### step 3 ####

  cat("Step 3... estimate empirical PDF/CDF of p-values (this may take a while!)\n")

  # determine cumulative empirical distribution function (pvalues)
  ee <- ecdf.pval(pval, eta0=eta0)

  g.pval <- grenander(ee)

  #cat("DEBUG: Grenander eta0=", g.pval$f.knots[length(g.pval$f.knots)], "\n")
  #cat("DEBUG: estimated eta0=", eta0 , "\n\n")

  # mixture density and CDF  
  f.pval = approxfun( g.pval$x.knots,  g.pval$f.knots, method="constant" )
  f0.pval = function(x) return( ifelse(x > 1 | x < 0, 0, rep(1, length(x))) )  

  F.pval = approxfun( g.pval$x.knots,  g.pval$F.knots, method="linear", 
           yleft=0, yright=g.pval$F.knots[length(g.pval$F.knots)])
  F0.pval = function(x) return( ifelse(x > 1, 1, ifelse(x < 0, 0, x )) )

  fdr.pval = function(p) pmin( eta0 / f.pval(p), 1)   # eta0*f0/ f
  Fdr.pval = function(p) pmin( eta0*p / F.pval(p), 1) # eta0*F0/ F
  

#### step 4 ####

  cat("Step 4... compute q-values and local fdr for each case\n")

  qval <- Fdr.pval(pval) 
  lfdr <- fdr.pval(pval)


 
#### return results ####

  if (statistic=="pvalue") 
    param = c(eta0)
  else
    param = c(eta0, cf.param)
  if (statistic=="studentt") attr(param, "names") <- c("eta0", "df")
  if (statistic=="normal") attr(param, "names") <- c("eta0", "sd")
  if (statistic=="correlation") attr(param, "names") <- c("eta0", "kappa")
  if (statistic=="pvalue") attr(param, "names") <- c("eta0")
  
  nm = list(pval=pval, qval=qval, lfdr=lfdr, 
             statistic=statistic, param=param)

  if (plot)
  {
    cat("Step 5... prepare for plotting\n")

    fdr <- function(ax) return( fdr.pval( 1- nf$F0(ax) ) )
    Fdr <- function(ax) return( Fdr.pval( 1- nf$F0(ax) ) )
     
    F  = function(ax) 1-eta0*(1-nf$F0(ax))/Fdr(ax) 
    FA = function(ax) (F(ax)-eta0*nf$F0(ax))/(1-eta0)		

    f = function(ax) eta0*(nf$f0(ax))/fdr(ax) 
    fA = function(ax) (f(ax)-eta0*nf$f0(ax))/(1-eta0)		

    xxx = seq(0, max(ax), length.out=200)
    ll = pvt.plotlabels(statistic, cf.param, eta0)

    par(mfrow=c(3,1))
    hist(ax, freq=FALSE, #bre=30,
      main=ll$main, xlab=ll$xlab, cex.main=1.8)
    lines(xxx, eta0*nf$f0(xxx), col=2, lwd=2, lty=3 )
    #lines(xxx, f(xxx), col=1, lwd=1 ) # show histogram instead
    lines(xxx, (1-eta0)*fA(xxx), col=4, lwd=2 )
    if (statistic=="pvalue") 
      pos1 = "topleft" else pos1="topright"
    legend(pos1, 
      c("Mixture", "Null Component", "Alternative Component"), 
      lwd=c(1, 2, 2), col=c(1, 2, 4), lty=c(1,3,1), bty="n", cex=1.5)
 
    plot(xxx, F(xxx), col=1, lwd=1, type="l", ylim=c(0,1),
      main="Density (first row) and Distribution Function (second row)",
      xlab=ll$xlab, ylab="CDF", cex.main=1.5)
    lines(xxx, eta0*nf$F0(xxx), col=2, lwd=2, lty=3)
    lines(xxx, (1-eta0)*FA(xxx), col=4, lwd=2)

    plot(xxx, Fdr(xxx), type="l", lwd=2, ylim=c(0,1),
      main="(Local) False Discovery Rate", ylab="Fdr and fdr",
      xlab=ll$xlab, col=1, lty=3, cex.main=1.5)
    lines(xxx, fdr(xxx), lwd=2, col=1)
    if (eta0 > 0.98) 
      pos2 = "bottomleft" else pos2="topright"
    legend(pos2, 
      c("fdr (density-based)", "Fdr (tail area-based)"), 
      lwd=c(2,2), col=c(1,1), lty=c(1,3), bty="n", cex=1.5)


    par(mfrow=c(1,1))

  }

  return(nm)
}


#####

## create labels for plots
pvt.plotlabels <- function(statistic, cf.param, eta0)
{
   if (statistic=="pvalue")
   {
     main = paste("Type of Statistic: p-Value (eta0 = ", round(eta0, 4), ")", sep="")
     xlab ="1-pval"
   }

   if (statistic=="studentt")
   {
     df = cf.param 
     main = paste("Type of Statistic: t-Score (df = ", round(df,3),
                       ", eta0 = ", round(eta0, 4), ")", sep="")
     xlab = "abs(t)"
   }

   if (statistic=="normal")
   {
     sd = cf.param 
     main = paste("Type of Statistic: z-Score  (sd = ", round(sd,3),
                       ", eta0 = ", round(eta0, 4), ")", sep="")
     xlab = "abs(z)"
   }

   if (statistic=="correlation")
   {
     kappa =cf.param      
     main = paste("Type of Statistic: Correlation (kappa = ", round(kappa,1),
                       ", eta0 = ", round(eta0, 4), ")", sep="")
     xlab = "abs(r)"
   }

   return(list(main=main, xlab=xlab))
}

pvt.nullfunction <- function(statistic, cf.param)
{
   if (statistic!="pvalue") attr(cf.param, "names") <- NULL

   if (statistic=="pvalue")
   {
     f0 = function(x) rep(1, length(x))
     F0 = function(x) x
   }

   if (statistic=="studentt")
   {
     df = cf.param       
     f0 = function(x) 2*dt(x, df=df)
     F0 = function(x) 2*pt(x, df=df)-1  
   }

   if (statistic=="normal")
   {
     sd = cf.param
     f0 = function(x) dhalfnorm(x, theta=sd2theta(sd))
     F0 = function(x) phalfnorm(x, theta=sd2theta(sd))
   }

   if (statistic=="correlation")
   {
     kappa =cf.param      
     f0 = function(x) 2*dcor0(x, kappa=kappa)
     F0 = function(x) 2*pcor0(x, kappa=kappa)-1 
   }

   checkf0 = function(ax) pmax(f0(ax),0) # make sure f0 is in [0,Inf]
   checkF0 = function(ax) pmax(pmin(F0(ax),1),0) # make sure F0 is in [0,1]

   return(list(f0=checkf0, F0=checkF0))
 }

# empirical cumulative distribution of p-values,
# constrained such that the known fraction of null p-values is taken into account
ecdf.pval <- function (x, eta0=1) 
{
    x <- sort(x)
    n <- length(x)
    if (n < 1) 
        stop("'x' must have 1 or more non-missing values")
    vals <- sort(unique(x))
    F.raw <- cumsum(tabulate(match(x, vals)))/n
    
    # this is the bit that makes sure that the gradient of 1-F(1-pval) is >= eta0 
    F.raw <- pmin(F.raw, 1-eta0*(1-vals) ) 
    
    rval <- approxfun(vals, F.raw, 
        method = "constant", yleft = 0, yright = 1, f = 0, ties = "ordered")
    class(rval) <- c("ecdf", "stepfun", class(rval))
    attr(rval, "call") <- sys.call()
    rval
}


