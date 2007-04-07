\name{fdrtool}
\alias{fdrtool}
\title{Estimate (Local) False Discovery Rates For Diverse Test Statistics}
\usage{
fdrtool(x, statistic=c("normal", "correlation", "pvalue", "studentt"),
  plot=TRUE, verbose=TRUE, censored.fit.arg=NULL, pval.estimate.eta0.arg=NULL, use.locfdr=TRUE)
}
\description{
  \code{fdrtool} takes a vector of z-scores (or of correlations, p-values,
  or t-statistics), and estimates for each case both the tail area-based Fdr
  as well as the density-based fdr (=q-value resp. local false discovery rate).
  The parameters of the null distribution are 
  estimated adaptively from the data (except for the case of p-values where
  this is not necessary).
}

\arguments{
  \item{x}{vector of the observed test statistics.}
  \item{statistic}{one of ""normal" (default), "correlation", "pvalue", 
                   "studentt".  This species the null model.}
  \item{plot}{plot a figure with estimated densities, distribution functions, 
             and (local) false discovery rates.}
  \item{verbose}{print out status messages.}
  \item{censored.fit.arg}{list with arguments passed on to 
                              \code{\link{censored.fit}}.}
  \item{pval.estimate.eta0.arg}{list with arguments passed on to 
                                 \code{\link{pval.estimate.eta0}}.}
  \item{use.locfdr}{employ the locfdr package for estimation eta0 and 
        scale parameter (only for "normal" and "correlation" statistic).}

}
\details{
  The algorithm implemented in this function proceeds as follows:

 \enumerate{
   \item  The parameters of the null model are estimated from the 
          data using \code{\link{censored.fit}}.
   \item  Subsequently the corresponding p-values are computed.
   \item  The proportion of null-p-values \code{eta0} is estimated using 
          \code{\link{pval.estimate.eta0}}.
   \item  A modified \code{\link{grenander}} algorithm is employed
          to obtain the overall density and distribution function 
          (note that this respects the estimated \code{eta0}).
   \item  Finally, q-values and local fdr values are computed for each case.
  }

  The assumed null models all have (except for p-values) one free
  scale parameter.  Note that the z-scores are assumed to have zero mean. 

}
\value{
  A list with the following components:

  \item{pval}{a vector with p-values for each case.} 
  \item{qval}{a vector with q-values (Fdr) for each case.}
  \item{lfdr}{a vector with local fdr values for each case.}
  \item{statistic}{the specified type of null model.} 
  \item{param}{a vector containing the estimated parameters (the null 
               proportion \code{eta0}  and the free parameter of the null model).}  
}

\author{
  Korbinian Strimmer (\url{http://strimmerlab.org}).
}
\seealso{\code{\link{pval.estimate.eta0}}, \code{\link{censored.fit}}.}

\examples{
# load "fdrtool" library and p-values
library("fdrtool")
data(pvalues)


# estimate fdr and Fdr from p-values

data(pvalues)
fdr <- fdrtool(pvalues, statistic="pvalue")
fdr$qval # estimated Fdr values 
fdr$lfdr # estimated local fdr 


# estimate fdr and Fdr from z-scores

sd.true = 2.232
n = 500
z = rnorm(n, sd=sd.true)
z = c(z, runif(30, 5, 10)) # add some contamination
fdr <- fdrtool(z)

# you may change some parameters of the underlying functions
fdr <- fdrtool(z, censored.fit.arg=list( pct0=0.9)) 
}
\keyword{htest}
