\name{dist.snp}
\alias{dist.snp}
\title{Allelic Sharing Distance with SNPs}
\description{
  This function computes the allelic sharing distance (ASD) with
  (strict) SNPs.
}
\usage{
dist.snp(x)
}
\arguments{
  \item{x}{an object of class \code{"loci"}.}
}
\details{
  The ASD between two diploid genotypes is (Gao and Martin, 2009):

  \deqn{\frac{1}{L}\sum_{j=1}^L d_j}{(1/L)sum d_j, j = 1, ... L}

  where \eqn{L} is the number loci, \eqn{d_j} is the value for the
  \eqn{j}th locus: 0 if both genotypes are identical, 1 if they have one
  allele in common, or 2 if they have no allele in common.

  This function works only with biallelic, unphased, diploid
  genotypes. (These do not actually need to be SNPs: for instance, they
  could be A/A, A/a, or a/a, or else). Note that the required conditions
  are not checked by the present function: see the functions below.
}
\value{
  an object of class \code{"dist"}.
}
\references{
  Gao, X. and Martin, E. R. (2009) Using allele sharing distance for
  detecting human population stratification. \emph{Human Hederity},
  \bold{68}, 182--191.
}
\author{Emmanuel Paradis}
\seealso{
  \code{\link{is.snp}}, \code{\link{is.phased}}, \code{\link{getPloidy}}
}
%\examples{
%}
\keyword{manip}
\keyword{multivariate}
\keyword{cluster}
