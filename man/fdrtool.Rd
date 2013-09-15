\name{fdrtool}
\alias{fdrtool}
\title{Estimate (Local) False Discovery Rates For Diverse Test Statistics}
\usage{
fdrtool(x, statistic=c("normal", "correlation", "pvalue"),
  plot=TRUE, color.figure=TRUE, verbose=TRUE, 
  cutoff.method=c("fndr", "pct0", "locfdr"),
  pct0=0.75)
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
  \item{statistic}{one of "normal" (default), "correlation", "pvalue".  
         This species the null model.}
  \item{plot}{plot a figure with estimated densities, distribution functions, 
             and (local) false discovery rates.}
  \item{verbose}{print out status messages.}
  \item{cutoff.method}{one of "fndr" (default), "pct0", "locfdr".}
  \item{pct0}{fraction of data used for fitting null model - only if \code{cutoff.method}="pct0"}
  \item{color.figure}{determines whether a color figure or a black and white
        figure is produced (defaults to "TRUE", i.e. to color figure).}

}
\details{
  The algorithm implemented in this function proceeds as follows:

 \enumerate{
   \item  A suitable cutoff point is determined.  If \code{cutoff.method}
          is "fndr" then first an approximate null model is fitted and
          subsequently a cutoff point is sought with false nondiscovery
          rate as small as possible (see \code{\link{fndr.cutoff}}). 
          If \code{cutoff.method} is "pct0"
          then a specified quantile (default value: 0.75) of the data
          is used as the cutoff point.  If \code{cutoff.method} equals
          "locfdr" then the heuristic of the "locfdr" package (version 1.1-6)
          is employed to find the cutoff (z-scores and correlations only).
   \item  The parameters of the null model are estimated from the 
          data using \code{\link{censored.fit}}. This results
          in estimates for scale parameters und and proportion
          of null values (\code{eta0}). 
   \item  Subsequently the corresponding p-values are computed, and
          a modified \code{\link{grenander}} algorithm is employed
          to obtain the overall density and distribution function 
          (note that this respects the estimated \code{eta0}).
   \item  Finally, q-values and local fdr values are computed for each case.
  }

  The assumed null models all have (except for p-values) one free
  scale parameter.  Note that the z-scores and the correlations
  are assumed to have zero mean. 

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
\references{

  Strimmer, K. (2008a).   A unified approach to false discovery 
  rate estimation. BMC Bioinformatics 9: 303.
  Available from \url{http://www.biomedcentral.com/1471-2105/9/303/}.

  Strimmer, K. (2008b). fdrtool: a versatile R package for estimating 
  local and tail area- based false discovery rates.
  Bioinformatics 24: 1461-1462.
  Available from \url{http://bioinformatics.oxfordjournals.org/cgi/content/abstract/24/12/1461}.

}

\seealso{\code{\link{pval.estimate.eta0}}, \code{\link{censored.fit}}.}

\examples{
# load "fdrtool" library and p-values
library("fdrtool")
data(pvalues)


# estimate fdr and Fdr from p-values

data(pvalues)
fdr = fdrtool(pvalues, statistic="pvalue")
fdr$qval # estimated Fdr values 
fdr$lfdr # estimated local fdr 

# the same but with black and white figure  
fdr = fdrtool(pvalues, statistic="pvalue", color.figure=FALSE)


# estimate fdr and Fdr from z-scores

sd.true = 2.232
n = 500
z = rnorm(n, sd=sd.true)
z = c(z, runif(30, 5, 10)) # add some contamination
fdr = fdrtool(z)

# you may change some parameters of the underlying functions
fdr = fdrtool(z, cutoff.method="pct0", pct0=0.9) 
}
\keyword{htest}
