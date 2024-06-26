\name{amova}
\alias{amova}
\alias{print.amova}
\alias{getPhi}
\alias{write.pegas.amova}
\title{Analysis of Molecular Variance}
\description{
  This function performs a hierarchical analysis of molecular variance
  as described in Excoffier et al. (1992). This implementation accepts
  any number of hierarchical levels.
}
\usage{
amova(formula, data = NULL, nperm = 1000, is.squared = FALSE)
\method{print}{amova}(x, ...)
getPhi(sigma2)
write.pegas.amova(x, file = "")
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
  \item{sigma2}{a named vector of variance components.}
  \item{file}{a file name.}
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

  The function \code{getPhi()} is a convenience function for extracting a
  table of hierarchical Phi-statistics for reporting. This will be an N+1
  by N matrix where N is the number of hierarchcial levels and GLOBAL is
  always the first row of the matrix. The matrix can read as COLUMN in ROW.

  If the variance components passed to getPhi() are not named, they will
  be reported as "level 1", "level 2", etc.
}
\note{
  If there are more than three levels, approximate formulae are used to
  estimate the variance components.

  If there is an error message like this:

  \code{Error in FUN(X[[1L]], ...) : 'bin' must be numeric or a factor}

  it may be that the factors you use in the formula were not read
  correctly. You may convert them with the function \code{factor}, or,
  before reading your data files, do this command (in case this option
  was modified):

  \code{options(stringsAsFactors = TRUE)}
}
\value{
  An object of class \code{"amova"} which is a list with a table of sums
  of square deviations (SSD), mean square deviations (MSD), and the
  number of degrees of freedom, and a vector of variance components.
}
\author{Emmanuel Paradis, Zhian N. Kamvar, and Brian Knaus}
\references{
  Excoffier, L., Smouse, P. E. and Quattro, J. M. (1992) Analysis of
  molecular variance inferred from metric distances among DNA
  haplotypes: application to human mitochondrial DNA restriction
  data. \emph{Genetics}, \bold{131}, 479--491.
}
\seealso{
  \code{amova} in \pkg{ade4} for an implementation of the original
  Excoffier et al.'s model; \code{adonis} in \pkg{vegan} for a general
  (multivariate) implementation of an ANOVA framework with distances.
}
\examples{
### All examples below have 'nperm = 100' for faster execution times.
### The default 'nperm = 1000' is recommended.
require(ape)
data(woodmouse)
d <- dist.dna(woodmouse)
g <- factor(c(rep("A", 7), rep("B", 8)))
p <- factor(c(rep(1, 3), rep(2, 4), rep(3, 4), rep(4, 4)))
(d_gp <- amova(d ~ g/p, nperm = 100)) # 2 levels
sig2 <- setNames(d_gp$varcomp$sigma2, rownames(d_gp$varcomp))
getPhi(sig2) # Phi table
amova(d ~ p, nperm = 100) # 1 level
amova(d ~ g, nperm = 100)

## 3 levels (quite slow):
\dontrun{
pop <- gl(64, 5, labels = paste0("pop", 1:64))
region <- gl(16, 20, labels = paste0("region", 1:16))
conti <- gl(4, 80, labels = paste0("conti", 1:4))
dd <- as.dist(matrix(runif(320^2), 320))
(dd_crp <- amova(dd ~ conti/region/pop, nperm = 100))
sig2 <- setNames(dd_crp$varcomp$sigma2, rownames(dd_crp$varcomp))
getPhi(sig2)
}
}
\keyword{models}
