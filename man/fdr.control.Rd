\name{fdr.control}
\alias{fdr.control}
\title{Controlling the False Discovery Rate in Multiple Testing}
\usage{
fdr.control(p, Q=0.05, eta0, robust=FALSE, \dots)
}
\description{
  \code{fdr.control} controls the False Discovery Rate (FDR) at a
  given level Q.
}

\arguments{
  \item{p}{vector of p-values}
  \item{Q}{desired FDR level}
  \item{eta0}{proportion of null p-values. If value is not specified (the default) it is estimated
              via the function \code{\link{fdr.estimate.eta0}}.}
  \item{robust}{use small sample approximation for estimating q-values (default value: \code{FALSE})}
  \item{\dots}{arguments passed on to \code{\link{fdr.estimate.eta0}}.}
}
\details{
  The FDR is the expected proportion
  of false positives (erroneous rejections) among the significant tests (rejections).
 \code{fdr.control} uses the algorithms described in Benjamini and Hochberg (1995)
  and Storey (2002).  
  For a given vector of p-values and the desired FDR level Q the corresponding p-value
  cut-off and the q-values for each hypothesis (see Storey, 2002) are computed. 

Notes:
\enumerate{
 \item the default settings correspond to the step-up procedure to control the FDR 
    by Benjamini and Hochberg (1995).
    
 \item q-values for each hypothesis are computed as defined in Storey (2002).
 
 \item small sample approximation for q-value (robust=TRUE) is from Storey (2002).
 \item default eta0=0 is safe but also most conservative choice (for other possibilities
    see \code{\link{fdr.estimate.eta0}}).
}
}
\value{

  A list object with the following components:
  \item{qvalues}{a vector with the q-values for each hypothesis.} 
  \item{significant}{a vector with a TRUE/FALSE value for each hypothesis}
  \item{num.significant}{number of significant hypotheses.}
  \item{pvalue.cutoff}{cutoff level for the individual p-values to obtain the
        desired control of FDR.
        Hypotheses whose corresponding p-values are below or equal to this
        cut-off level are rejected (i.e. significant).} 
  \item{qvalue.cutoff}{the specified q-value cut-off (Q)}
  \item{eta0}{the fraction of null p-values used to estimate the cut-offs}  
}

\references{
 
 Benjamini, Y., and Y. Hochberg (1995)  Controlling the false
 discovery rate: a practical and powerful approach to multiple testing.
 \emph{J. Roy. Statist. Soc. B}, \bold{57}, 289--300.

 Storey, J. D. (2002) A direct approach to false
 discovery rates.
 \emph{J. Roy. Statist. Soc. B.}, \bold{64}, 479--498.
}
\author{
  Korbinian Strimmer (\url{http://strimmerlab.org}).
  
  Adapted in part from S-PLUS code by  Y. Benjamini (\url{http://www.math.tau.ac.il/~roee/FDR_Splus.txt})
  and R code from J.D. Storey (\url{http://faculty.washington.edu/~jstorey/}). 
}
\seealso{\code{\link{fdr.estimate.eta0}}.}

\examples{
# load fdrtool library and p-values
library("fdrtool")
data(pval)


# Proportion of null p-values for different methods
fdr.estimate.eta0(pval, method="conservative")
fdr.estimate.eta0(pval, method="adaptive")
fdr.estimate.eta0(pval, method="bootstrap")
fdr.estimate.eta0(pval, method="smoother")


# FDR test on the level 0.05 using the default "smoother" method
fdr.control(pval, Q = 0.05)$num.significant 

# FDR test on the level 0.05 using the conservative choice eta0=1
fdr.control(pval, Q = 0.05, method="conservative")$num.significant 
}
\keyword{htest}
