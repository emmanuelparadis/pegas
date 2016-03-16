\name{na.omit.loci}
\alias{na.omit.loci}
\title{Missing Allelic Data}
\description{
  This is a method of the generic function \code{\link[stats]{na.omit}}.
}
\usage{
\method{na.omit}{loci}(object, na.alleles = c("0", "."), ...)
}
\arguments{
  \item{object}{an object of class \code{"loci"}.}
  \item{na.alleles}{a vector of charater strings giving the alleles
    to be treated as missing data.}
  \item{\dots}{unused.}
}
\details{
  The side effect of this function is to drop the rows (individuals)
  with unclearly identified genotypes, i.e., with at least one allele
  among \code{na.alleles}.

  Other variables in the data table are eventually checked and levels
  with no observation (e.g., population) are dropped.
}
\value{
  an object of class \code{"loci"}.
}
\author{Emmanuel Paradis}
\examples{
data(jaguar)
nrow(jaguar)
nrow(na.omit(jaguar))
}
\keyword{manip}
