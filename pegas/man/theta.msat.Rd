\name{theta.msat}
\alias{theta.msat}
\title{
  Population Parameter THETA From Micro-Satellites
}
\description{
  This function estimates the population parameter \eqn{\theta}{THETA}
  using micro-satellite data with three different estimators.
}
\usage{
theta.msat(x)
}
\arguments{
  \item{x}{an object of class \code{"loci"}.}
}
\details{
  The data must be micro-satellites, so the allele names must be the
  repeat counts (see the example).

  The three estimators are based on (i) the variance of the number of
  repeats, (ii) the expected homozygosity (both described in Kimmel et
  al., 1998), and (iii) the mean allele frequencies (Haasl and Payseur,
  2010).
}
\value{
  a numeric matrix with loci as rows and the three estimates of
  \eqn{\theta}{THETA} as columns.
}
\references{
  Kimmel, M., Chakraborty, R., King, J. P., Bamshad, M., Watkins,
  W. S. and Jorde, L. B. (1998) Signatures of population expansion in
  microsatellite repeat data. \emph{Genetics}, \bold{148}, 1921--1930.

  Haasl, R. J. and Payseur, B. A. (2010) The number of alleles at a
  microsatellite defines the allele frequency spectrum and facilitates
  fast accurate estimation of \eqn{\theta}{THETA}. \emph{Molecular
  Biology and Evolution}, \bold{27}, 2702--2715.
}
\author{Emmanuel Paradis}
\seealso{
  \code{\link{theta.h}}, \code{\link{theta.tree}}
}
\examples{
require(adegenet)
data(nancycats)
x <- as.loci(nancycats)
theta.msat(x)
}
\keyword{models}
