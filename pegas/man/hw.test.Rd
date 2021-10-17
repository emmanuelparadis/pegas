\name{hw.test}
\alias{hw.test}
\alias{hw.test.loci}
\alias{hw.test.genind}
\title{Test of Hardy--Weinberg Equilibrium}
\description{
  This function tests, for a series of loci, the hypothesis that
  genotype frequencies follow the Hardy--Weinberg equilibrium.
  \code{hw.test} is a generic with methods for the classes
  \code{"\link{loci}"} and \code{genind}. Note that the latter
  replaces \code{HWE.test.genind} in the \pkg{adegenet} package.
}
\usage{
hw.test(x, B = 1000, \dots)
\method{hw.test}{loci}(x, B = 1000, \dots)
\method{hw.test}{genind}(x, B = 1000, \dots)
}
\arguments{
  \item{x}{an object of class \code{"\link{loci}"} or \code{genind}.}
  \item{B}{the number of replicates for the Monte Carlo procedure; for
    the regular HW test, set B = 0 (see details).}
  \item{\dots}{further arguments to be passed.}
}
\details{
  This test can be performed with any level of ploidy. Two versions
  of the test are available: the classical \eqn{\chi^2}{chi^2}-test based
  on the expected genotype frequencies calculated from the allelic
  frequencies, and an exact test based on Monte Carlo permutations of
  alleles (Guo and Thompson 1992). For the moment, the latter version is
  available only for diploids. Set \code{B = 0} if you want to skip the
  second test.
}
\value{
  A matrix with three or four columns with the \eqn{\chi^2}{chi^2}-value,
  the number of degrees of freedom, the associated \emph{P}-value, and
  possibly the \emph{P}-value from the Monte Carlo test. The rows of
  this matrix are the different loci in \code{x}.
}
\author{
  Main code by Emmanuel Paradis; wrapper for \code{genind}
  objects by Thibaut Jombart.
}
\references{
  Guo, S. W. and Thompson, E. A. (1992) Performing the exact test of
  Hardy--Weinberg proportion for multiple alleles. \emph{Biometrics},
  \bold{48}, 361--372.
}
\examples{
\dontrun{
require(adegenet)

## load data
data(nancycats)

## test on genind object, no permutation
hw.test(nancycats, B=0)

## test on loci object
x <- as.loci(nancycats)
hw.test(x)
}
data(jaguar)
hw.test(jaguar)
}
\keyword{htest}
