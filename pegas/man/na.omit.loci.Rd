\name{na.omit.loci}
\alias{na.omit.loci}
\alias{nullAlleles2NA}
\title{Missing Allelic Data}
\description{
  The first function is a method of the generic function
  \code{\link[stats]{na.omit}}.

  \code{nullAlleles2NA} changes all genotypes with at least one `null'
  allele (that is among the values in \code{na.alleles}) into \code{NA}.
}
\usage{
\method{na.omit}{loci}(object, na.alleles = c("0", "."), ...)

nullAlleles2NA(object, na.alleles = c("0", "."))
}
\arguments{
  \item{object}{an object of class \code{"loci"}.}
  \item{na.alleles}{a vector of character strings giving the alleles
    to be treated as missing data.}
  \item{\dots}{(unused)}
}
\details{
  The side effect of \code{na.omit} is to drop the rows (individuals)
  with unclearly identified genotypes, i.e., with at least one allele
  among \code{na.alleles}.

  Other variables in the data table are eventually checked and levels
  with no observation (e.g., population) are dropped.

  \code{nullAlleles2NA} does not remove any observation but changes
  these genotypes into \code{\link[base]{NA}}.
}
\value{an object of class \code{"loci"}.}
\author{Emmanuel Paradis}
\examples{
data(jaguar)
nrow(jaguar)
nrow(na.omit(jaguar))
nrow(nullAlleles2NA(jaguar))
}
\keyword{manip}
