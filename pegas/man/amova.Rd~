\name{amova}
\alias{amova}
\alias{print.amova}
\title{Analysis of Molecular Variance}
\description{
  This function performs a hierarchical analysis of molecular variance
  as described in Excoffier et al. (1992). This implementation accepts
  any number of hierarchical levels.
}
\usage{
amova(formula, data = NULL, nperm = 1000, is.squared = FALSE)
\method{print}{amova}(x, ...)
}
\arguments{
  \item{formula}{a formula giving the AMOVA model to be fitted with the
    distance matrix on the left-hand side of the \code{~}, and the
    population, region, etc, levels on its right-hand side (see details).}
  \item{data}{an optional data frame where to find the hierarchical
    levels; by default they are searched for in the user's workspace.}
  \item{nperm}{the number of permutations for the tests of hypotheses
    (1000 by default). Set this argument to 0 to skip the tests and
    simply estimate the variance components.}
  \item{is.squared}{a logical specifying whether the distance matrix has
    already been squared.}
  \item{x}{an object of class \code{"amova"}.}
  \item{\dots}{unused (here for compatibility.}
}
\details{
  The formula must be of the form \code{d ~ A/B/...} where \code{d} is a
  distance object, and \code{A}, \code{B}, etc, are the hierarchical
  levels from the highest to the lowest one. Any number of levels is
  accepted, so specifying \code{d ~ A} will simply test for population
  differentiation.

  It is assumed that the rows of the distance matrix are in the same
  order than the hierarchical levels (which may be checked by the user).
}
\note{
  If there are more than three levels, approximate formulae are used to
  estimate the variance components.
}
\value{
  An object of class \code{"amova"} which is a list with a table of sums
  of square deviations (SSD), mean square deviations (MSD), and the
  number of degrees of freedom, and a vector of variance components.
}
\author{Emmanuel Paradis}
\references{
  Excoffier, L., Smouse, P. E. and Quattro, J. M. (1992) Analysis of
  molecular variance inferred from metric distances among DNA
  haplotypes: application to human mitochondrial DNA restriction
  data. \emph{Genetics}, \bold{131}, 479--491.
}
\seealso{
  \code{\link[ade4]{amova}} in \pkg{ade4} for an implementation of
  the original Excoffier et al.'s model;  \code{\link[vegan]{adonis}} in
  \pkg{vegan} for a general (multivariate) implementation of an ANOVA
  framework with distances.
}
\examples{
### All examples below have 'nperm = 100' for faster execution times.
### The default 'nperm = 1000' is recommended.
require(ape)
data(woodmouse)
d <- dist.dna(woodmouse)
g <- factor(c(rep("A", 7), rep("B", 8)))
p <- factor(c(rep(1, 3), rep(2, 4), rep(3, 4), rep(4, 4)))
amova(d ~ g/p, nperm = 100) # 2 levels
amova(d ~ p, nperm = 100) # 1 level
amova(d ~ g, nperm = 100)

## 3 levels (quite slow):
\dontrun{
pop <- gl(64, 5, labels = paste("pop", 1:64))
region <- gl(16, 20, labels = paste("region", 1:16))
conti <- gl(4, 80, labels = paste("conti", 1:4))
dd <- as.dist(matrix(runif(320^2), 320))
amova(dd ~ conti/region/pop, nperm = 100)
}}
\keyword{models}
