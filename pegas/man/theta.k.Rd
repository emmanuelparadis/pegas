\name{theta.k}
\alias{theta.k}
\title{Population Parameter THETA using Expected Number of Alleles}
\usage{
theta.k(x, n = NULL, k = NULL)
}
\arguments{
  \item{x}{a vector or a factor.}
  \item{n}{a numeric giving the sample size.}
  \item{k}{a numeric giving the number of alleles.}
}
\description{
  This function computes the population parameter THETA using the
  expected number of alleles.
}
\value{
  A numeric vector of length one with the estimated theta.
}
\details{
  This function can be used in two ways: either with a vector giving the
  individual genotypes from which the sample size and number of alleles
  are derived (e.g., \code{theta.k(x)}), or giving directly these two
  quantities (e.g., \code{theta.k(n = 50, k = 5)}).

  The argument \code{x} can be either a factor or a vector. If it is a
  factor, then it is taken to give the individual alleles in the
  population. If it is a numeric vector, then its values are taken to be
  the numbers of each allele in the population. If it is a non-numeric
  vector, it is a coerced as a factor.

  Both arguments \code{n} and \code{k} must be single numeric values.
}
\note{
  For the moment, no standard-error or confidence interval is computed.
}
\references{
  Ewens, W. J. (1972) The sampling theory of selectively neutral
  alleles. \emph{Theoretical Population Biology}, \bold{3}, 87--112.
}
\author{Emmanuel Paradis}
\seealso{
  \code{\link{theta.h}}, \code{\link{theta.s}}, \code{\link{theta.tree}}
}
\examples{
require(adegenet)
data(nancycats)
## convert the data and compute frequencies:
S <- summary(as.loci(nancycats))
## compute THETA for all loci:
sapply(S, function(x) theta.k(x$allele))
}
\keyword{manip}
\keyword{univar}
