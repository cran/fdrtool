### fdr.estimate.eta0.R  (2006-08-11)
###
###     Estimating the Proportion of Null p-Values
###
### Copyright 2003-06 Korbinian Strimmer 
###
### Parts of this code is adapted from 
### S-PLUS code (c) by Y. Benjamini (available from
### http://www.math.tau.ac.il/~roee/FDR_Splus.txt ) 
### and from R code (c) J.D. Storey (available from 
### http://faculty.washington.edu/~jstorey/qvalue/ )
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


#Input
#=============================================================================
#p: a vector of p-values 
#method:  method for computing eta0
#lambda:  optional tuning parameter (vector, needed for "bootstrap" and "smoothing") 
#
#  conservative: Benjamini and Hochberg (1995) JRSSB
#  adaptive:     Benjamini and Hochberg (2000) J. Behav. Educ. Statist.
#  bootstrap:    Storey (2002) JRSSB
#  smoother:     Storey and Tibshirani (2003) PNAS

#Output
#=============================================================================
#eta0: an estimate of the proportion of null p-values

fdr.estimate.eta0 <- function(p,
    method=c("smoother", "bootstrap", "conservative", "adaptive"),
    lambda=seq(0,0.9,0.05), diagnostic.plot=TRUE)
{
    method <- match.arg(method)
        
    
    ########
    
    if (method == "conservative") # Benjamini and Hochberg (1995)
    {
        eta0 = 1.0
    }

    ########
   
    if (method == "adaptive") # Benjamini and Hochberg (2000)
    {
        m <- length(p)
	sortp <- sort(p)
	s <- sort(1 - sortp)/(1:m)
	
	m0raw <- m
        i <- m
        while(i > 1 && s[i] <= s[i - 1]) i <- i - 1
        if(i > 1)
            m0raw <- 1/s[i - 1]
        else m0raw <- 1/s[1]
        
	m0 <- min(floor(1 + m0raw), m)
        
	eta0 <- m0/m      
    }

    ########
    if(method == "bootstrap" || method == "smoother")
    {
    
      # for the remaining methods we require a lambda vector
      if (length(lambda)<4)
        stop("At least 4 values in lambda tuning vector required")
    
      eta0 <- rep(0,length(lambda))
      for(i in 1:length(lambda))
      {
        eta0[i] <- mean(p >= lambda[i])/(1-lambda[i])
      }
    }
        
    ########
   
    if(method == "bootstrap") # Storey (2002) JRSSB
    {
            m <- length(p)
	    mineta0 <- min(eta0)
            mse <- rep(0,length(lambda))
            eta0.boot <- rep(0,length(lambda))
            for(i in 1:100) {
                p.boot <- sample(p,size=m,replace=TRUE)
                for(i in 1:length(lambda)) {
                    eta0.boot[i] <- mean(p.boot>lambda[i])/(1-lambda[i])
                }
                mse <- mse + (eta0.boot-mineta0)^2
            }
            eta0 <- min(eta0[mse==min(mse)])
            eta0 <- min(eta0,1)
    }    
 
    ########
 
    if(method == "smoother") # Storey and Tibshirani (2003) PNAS
    {
            seta0 <- smooth.spline(lambda,eta0,df=3)
            eta0 <- predict(seta0,x=max(lambda))$y
            eta0 <- min(eta0,1)
    } 
    
    #######
    
    if (diagnostic.plot)
    {
      get(getOption("device"))() # this more platform-independent than X11()
  
      info0 = paste("method =", method)
      info1 = paste("eta0 =", round(eta0, 4))
  
      dpr = density.pr(p, plot=TRUE, df=30, 
                    main="Diagnostic Plot:  Distribution of p-Values and Estimated eta0", xlim=c(0,1), 
                    xlab="p-values")
   
  
      y0 = dunif(dpr$x)*eta0
      lines(dpr$x, y0, col=2)
      maxy = max(dpr$histogram$density)
  
      text(0.5, 5/6*maxy, info0) 
      text(0.5, 4/6*maxy, info1)  
    }
    
    
    #######
    
    
    
    return(eta0)
   
}
