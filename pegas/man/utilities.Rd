\name{utilities}
\alias{getPloidy}
\alias{getAlleles}
\alias{getGenotypes}
\alias{expand.genotype}
\alias{proba.genotype}
\alias{is.snp}
\alias{is.phased}
\title{Utily Functions for pegas}
\description{
  The first three functions extract information on loci,
  \code{expand.genotype} creates a table of all possible genotypes given
  a set of alleles, \code{proba.genotype} calculates expected
  probabilities of genotypes under Hardy--Weinberg equilibrium, and the
  last two functions test whether a locus is a SNP or whether a
  gentotype is phased.
}
\usage{
getPloidy(x)
getAlleles(x)
getGenotypes(x)
expand.genotype(n, alleles = NULL, ploidy = 2, matrix = FALSE)
proba.genotype(alleles = c("1", "2"), p, ploidy = 2)
is.snp(x)
is.phased(x)
}
\arguments{
  \item{x}{an object of class \code{"loci"}.}
  \item{n}{an integer giving how many alleles to consider (ignored if
    \code{alleles} is used).}
  \item{alleles}{the allele names as a vector of mode character.}
  \item{ploidy}{an integer giving the ploidy level (either 2 or 4 for
    the moment).}
  \item{matrix}{a logical specifying whether to return the genotypes in
    a matrix or as a character vector.}
  \item{p}{a vector of allele probabilities; if missing, equal
    probabilities are assumed.}
}
\details{
  \code{expand.genotype} and \code{proba.genotype} accept any level of
  ploidy and any number of alleles.

  For \code{is.snp}, a locus is defined as a SNP if it has two alleles
  and their labels are made of a single character (e.g., A and T, or 1
  and 2, but not A and AT).
}
\value{
  \code{getPloidy} returns the ploidy level of all loci in an object of
  class \code{"loci"} as a numeric vector.

  \code{getAlleles} and \code{getGenotypes} return the alleles and
  genotypes, respectively, observed in all loci in an object of class
  \code{"loci"} as a list.

  \code{expand.genotype} returns a character vector (the default) or a
  matrix where the rows are the genotypes and the columns are the
  alleles. The matrix is numeric by default, or character if the
  argument \code{alleles} is given.

  \code{proba.genotype} returns a numeric vector with names set as the
  genotypes.

  \code{is.snp} returns a logical vector specifying whether each locus
  is a SNP.

  \code{is.phased} returns a matrix of the same size than the original
  data specifying whether each genotype is phased or not.
}
\author{Emmanuel Paradis}
\examples{
data(jaguar)
X <- jaguar[, 1:2]
getAlleles(X)
getGenotypes(X)
expand.genotype(2)
expand.genotype(2, LETTERS[1:3])
expand.genotype(3, ploidy = 4)
proba.genotype() # classical HWE with 2 alleles
## an octoploid with a six-allele locus (1287 possible genotypes):
length(p <- proba.genotype(alleles = LETTERS[1:6], ploidy = 8))
max(p) # ~ 0.006
## back to the jaguar data:
s <- summary(X)
## allele counts from the first locus:
p <- s[[1]]$allele
## expected probabilities for the 136 possible genotypes...
proba.genotype(names(p), p/sum(p))
## ... to be compared with s[[1]]$genotype
}
\keyword{manip}
