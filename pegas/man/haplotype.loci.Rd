\name{haplotype.loci}
\alias{haplotype.loci}
\alias{plot.haplotype.loci}
\alias{dist.haplotype.loci}
\title{Haplotype Extraction and Frequencies From Allelic Data}
\description{
  This function extracts haplotypes from phased genotypes.
}
\usage{
\method{haplotype}{loci}(x, locus = 1:2, quiet = FALSE, compress = TRUE,
           check.phase = TRUE, ...)
\method{plot}{haplotype.loci}(x, ...)
dist.haplotype.loci(x)
}
\arguments{
  \item{x}{an object of class \code{"loci"} or of class
    \code{"haplotype.loci"}.}
  \item{locus}{a vector of integers giving the loci to analyse.}
  \item{quiet}{a logical value specifying whether to not print the
    progress of the analysis (\code{FALSE} by default).}
  \item{compress}{by default only the unique haplotypes are returned
    with their frequencies. If \code{compress = FALSE}, a matrix with
    all observed haplotypes is returned (with the number of columns
    equals to the number of individuals times the ploidy level).}
  \item{check.phase}{a logical value specifying whether to check if the
    individual genotypes are phased.}
  \item{\dots}{arguments passed to and from methods.}
}
\details{
  The individuals with at least one unphased genotype are ignored with a
  warning.

  \code{dist.haplotype.loci} computes pairwise distances among
  haplotypes by counting the number of different alleles.

  Checking whether the genotypes are phased can be time consuming with
  very big data sets. It may be useful to set \code{check.phase = FALSE}
  if several analyses are done on the same data and no warning was
  issued after the first scan, or you are sure that the genotypes are phased.
}
\note{
  \code{haplotype} is a generic function with methods for objects of
  class \code{"DNAbin"} and of class \code{"loci"}. Note that the class
  returned by these methods is different: \code{c("haplotype", "DNAbin")}
  and \code{"haplotype.loci"}, respectively. This and other details are
  likely to change in the future.
}
\value{
  \code{haplotype} returns a matrix of mode character with the loci as
  rows and the haplotypes as columns. The attribute \code{"freq"} gives
  the counts of each haplotype and the class is \code{"haplotype.loci"}.

  \code{dist.haplotype.loci} returns an object of class \code{"dist"}.
}
\author{Emmanuel Paradis}
\seealso{ \code{\link{haplotype}}, \code{\link{LD}}}
\keyword{manip}
