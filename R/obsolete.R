
### obsolete.R  (2007-02-17)
###
###     Obsolete Functions
###
### Copyright 2007 Korbinian Strimmer 
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



# this file contains contains warning functions that are activated
# when an obsolete function is used.

fdr.estimate.eta0 <- function(p, 
  method=c("smoother", "bootstrap", "conservative", "adaptive"),
   lambda=seq(0,0.9,0.05), diagnostic.plot=TRUE )
{
  stop(paste(
      "You are calling an obsolete function!!",  
      "Instead of 'fdr.estimate.eta0()' please use the function",
      "'pval.estimate.eta0()'."))
}

fdr.control <- function(p, Q=0.05, eta0, robust=FALSE, ...)
{
  stop(paste(
      "You are calling an obsolete function!!",  
      "Instead of 'fdr.control()' please use the function",
      "'fdrtool()'."))
}

