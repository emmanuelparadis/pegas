\name{alleles2loci}
\alias{alleles2loci}
\alias{loci2alleles}
\title{Build Loci Object From Matrix of Alleles}
\description{
  These functions transform a matrix of alleles into an object of class
  \code{"loci"}, or the reverse operation.
}
\usage{
alleles2loci(x, ploidy = 2, rownames = NULL, population = NULL,
             phased = FALSE)
loci2alleles(x)
}
\arguments{
  \item{x}{a matrix or a data frame where each column is an allele, or
    an object of class \code{"loci"}.}
  \item{ploidy}{an integer specifying the level of ploidy.}
  \item{rownames}{an integer giving the column number to be used as
    rownames of the output.}
  \item{population}{an integer giving the column number to be used as
    population (if any).}
  \item{phased}{a logical specifying whether the genotypes should be
    output as phased. By default, they are unphased.}
}
\details{
  Genetic data matrices are often arranged with one allele
  in each column of the matrix (particularly for micro-satellites), so
  that the number of columns is equal to the number of loci times the
  level of ploidy. \code{alleles2loci} transforms such matrices into a
  \code{"loci"} object.

  If the rownames of the input matrix are already set, they are used in
  the output. Alternatively, it is possible to specify which column to
  use as rownames (this column will be deleted before creating the
  genotypes).

  If the input matrix has colnames, then the names of the first column of
  each genotype is used as names of the output loci (see examples).

  \code{loci2alleles} checks that all individuals have the ploidy for a
  given locus (if not an error occurs), but ploidy can vary among loci.
}
\value{
  an object of class \code{"loci"} or a matrix.
}
\author{Emmanuel Paradis}
\seealso{
  \code{\link{read.loci}}, \code{\link{as.loci}}

  The vignette \dQuote{ReadingFiles} explains how to read such a data
  set from Dryad (\url{https://datadryad.org/stash}).
}
\examples{
x <- matrix(c("A", "A", "A", "a"), 2)
colnames(x) <- c("Loc1", NA)
y <- alleles2loci(x)
print(y, details = TRUE)
loci2alleles(y)
}
\keyword{IO}
