\name{hw.test}
\alias{hw.test}
\title{Test of Hardy--Weinberg Equilibrium}
\description{
  This function tests, for a series of loci, the hypothesis that
  genotype frequencies follow the Hardy--Weinberg equilibrium.
}
\usage{
hw.test(x, B = 1000)
}
\arguments{
  \item{x}{an object of class \code{"loci"}.}
  \item{B}{the number of replicates for the Monte Carlo procedure.}
}
\details{
  This test can be performed with any level of ploidy. Two versions
  of the test are available: the classical \eqn{\chi^2}{chi^2}-test based
  on the expected genotype frequencies calculated from the allelic
  frequencies, and an exact test based on Monte Carlo permutations of
  alleles. For the moment, the latter version is available only for
  diploids. Set \code{B = 0} if you want to skip the second test.
}
\value{
  A matrix with three or four columns with the \eqn{\chi^2}{chi^2}-value,
  the number of degrees of freedom, the associated \emph{P}-value, and
  possibly the \emph{P}-value from the Monte Carlo test. The rows of
  this matrix are the different loci in \code{x}.
}
\author{Emmanuel Paradis}
\examples{
require(adegenet)
data(nancycats)
x <- as.loci(nancycats)
hw.test(x)
}
\keyword{htest}
