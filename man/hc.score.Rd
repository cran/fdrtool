\name{hc.score}
\alias{hc.score}
\alias{hc.thresh}

\title{Compute Empirical Higher Criticism Scores and Corresponding Decision Threshold From p-Values}

\usage{
hc.score(pval)
hc.thresh(pval, alpha0=1, plot=TRUE)
}

\description{
  \code{hc.score} computes the empirical higher criticism (HC) scores from p-values.

  \code{hc.thresh} determines the HC decision threshold by searching for the p-value with the maximum HC score.
}

\arguments{
  \item{pval}{vector of p-values.}
  \item{alpha0}{look only at a fraction \code{alpha0} of the p-values (default: 1, i.e. all p-values).}
  \item{plot}{show plot with HC decision threshold.}
}
\details{

   Higher Criticism (HC) provides an alternative means to determine decision thresholds for
   signal identification, especially if the signal is rare and weak.

   See Donoho and Jin (2008) for details of this approach 
   and Klaus and Strimmer (2012) for a review and connections with FDR methdology. 
}
\value{
  \code{hc.score} returns a vector with the HC score corresponding to each p-value. 

  \code{hc.thresh} returns the p-value corresponding to the maximum HC score.
}

\author{
 Bernd Klaus and  Korbinian Strimmer (\url{http://strimmerlab.org}).
}
\references{

Donoho, D. and J. Jin. (2008). Higher criticism thresholding: optimal feature selection
when useful features are rare and weak. Proc. Natl. Acad. Sci. USA 105:14790-15795.

Klaus, B., and K. Strimmer (2013). Signal identification for rare and
weak features: higher criticism or false discovery rates?
Biostatistics 14: 129-143.  Preprint available from \url{http://arxiv.org/abs/1112.2615}.
}

\seealso{\code{\link{fdrtool}}.}

\examples{
# load "fdrtool" library
library("fdrtool")

# some p-values
data(pvalues)

# compute HC scores
hc.score(pvalues)

# determine HC threshold
hc.thresh(pvalues)
}
\keyword{htest}
