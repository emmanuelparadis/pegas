\name{heterozygosity}
\alias{heterozygosity}
\alias{H}
\alias{H.loci}
\alias{H.default}
\title{Heterozygosity at a Locus Using Gene Frequencies}
\description{
  Thes functions compute the mean heterozygosity(ies) from gene
  frequencies, and return optionally the associated variance(s).
}
\usage{
H(x, ...)
\method{H}{loci}(x, variance = FALSE, observed = FALSE, ...)
\method{H}{default}(x, variance = FALSE, ...)

heterozygosity(x, variance = FALSE)
}
\arguments{
  \item{x}{an object of class \code{"loci"}, or vector or a factor.}
  \item{variance}{a logical indicating whether the variance of the
    estimated heterozygosity should be returned (\code{TRUE}), the
    default being \code{FALSE}.}
  \item{observed}{a logical specifying whether to calculate the observed
    heterozygosity.}
  \item{...}{currently unused}
}
\details{
  The argument \code{x} can be either a factor or a vector. If it is a
  factor, then it is taken to give the individual alleles in the
  population. If it is a numeric vector, then its values are taken to be
  the numbers of each allele in the population. If it is a non-numeric
  vector, it is a coerced as a factor.

  The mean heterozygosity is estimated with:

  \deqn{\hat{H} = \frac{n}{n-1} \left(1 - \sum_{i=1}^k p_i^2 \right)}{%
    H = n(1 - SUM (FROM i=1 TO k) p_i^2)/(n - 1)}

  where \eqn{n} is the number of genes in the sample, \eqn{k} is the
  number of alleles, and \eqn{p_i} is the observed (relative) frequency
  of the \eqn{i}th allele.
}
\value{
  For the default method: a numeric vector of length one with the
  estimated mean heterozygosity (the default), or of length two if the
  variance is returned.

  For the \code{"loci"} method: a numeric matrix with one, two, or three
  columns with a row for each locus and the values of heterozygosity as
  columns.
}
\references{
  Nei, M. (1987) \emph{Molecular evolutionary genetics}. New York:
  Columbia University Press.
}
\author{Emmanuel Paradis}
\seealso{
  \code{\link{theta.s}}
}
\examples{
data(jaguar)
H(jaguar, TRUE, TRUE)
## use the (old) default method:
## convert the data and compute frequencies:
S <- summary(jaguar)
## compute H for all loci:
sapply(S, function(x) H(x$allele))
## ... and its variance
sapply(S, function(x) H(x$allele, variance = TRUE))
}
\keyword{manip}
\keyword{univar}
