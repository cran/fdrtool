### fdr.control.R  (2006-08-08)
###
###     Controlling False Discovery Rate in Multiple Testing
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


# FDR controlling procedures.
#
# The procedure below controls the False Discovery Rate (FDR) at a
# given level Q using the algorithms described in Benjamini and Hochberg (1995)
# and Storey (2002).  The FDR is the expected proportion
# of false positives (erroneous rejections) among the significant tests (rejections).
# For a given vector of p-values and the desired FDR level Q the corresponding p-value
# cut-off and the q-values for each hypothesis (see Storey (2002) ) are computed.
#
# Notes:  
# -the default settings correspond to the step-up procedure to control the FDR 
#    by Benjamini and Hochberg (1995)
# -q-values for each hypothesis are computed as defined in Storey (2002) JRSSB
# -small sample approximation for q-value (robust=TRUE) is from Storey (2002) JRSSB.
# -default eta0=0 is safe but most conservative choice (for other possibilities
#    see estimate.eta0)
#
# References:
#
# Benjamini, Y., and Y. Hochberg (1995)  Controlling the false
# discovery rate: a practical and powerful approach to multiple testing.
# J. Roy. Statist. Soc. B. 57:289-300
#
# Storey, J.D. (2002) A direct approach to false discovery rates. 
# J. Roy. Statist. Soc. B, 64: 479-498


#Input
#=============================================================================
#p:      a vector of p-values 
#Q:      a level at which to control the FDR (default: 0.05)
#eta0     an estimate of the proportion of null p-values (default 1)
#robust: an indicator of whether it is desired to make the estimate of q-values
#         more robust for small p-values (default: FALSE)
#Output
#=============================================================================
#qvalues          q-values for each hypothesis - see Storey (2002)
#significant      a vector with TRUE/FALSE value for each hypothesis
#num.significant  number of significant hypothesis
#pvalue.cutoff    corresponding p-value cut-off (all p <= pvalue.cutoff significant)

fdr.control <- function(p, Q=0.05, eta0, robust=FALSE, ...)
{ 
    if(min(p)<0 || max(p)>1)
    {
       stop("p-values not in valid range")
    }
    
    if( missing(eta0) ) eta0 <- fdr.estimate.eta0(p, ...)
    
      
    m <- length(p)
       
    # compute q-values
    u <- order(p)
    v <- rank(p)
    qvalue <- eta0*m*p/v
    if(robust)
    {
        qvalue <- eta0*m*p/(v*(1-(1-p)^m))
    }
    qvalue[u[m]] <- min(qvalue[u[m]],1)
    for(i in (m-1):1)
    {
        qvalue[u[i]] <- min(qvalue[u[i]],qvalue[u[i+1]],1)
    }

    # test hypothesis and compute p-value cutoff
    rej <- (qvalue <= Q)
    if ( sum(rej) == 0 )
    {
      cutoff <- 0
    }
    else
    {
      cutoff <- max(p[rej])
    }
              
    return( list(qvalues=qvalue, significant=rej,
                   num.significant = sum(rej), 
                   pvalue.cutoff=cutoff,
		   qvalue.cutoff=Q,
		   eta0=eta0)
           )
}
