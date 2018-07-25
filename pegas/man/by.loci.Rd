\name{by.loci}
\alias{by.loci}
\title{Summary by Population of Other Factor}
\description{
  This is an implementation of the generic \code{\link[base]{by}}
  function which applies a function to some data for a each level of a
  categorical factor.
}
\usage{
\method{by}{loci}(data, INDICES = data$population, FUN = NULL, ..., simplify = TRUE)
}
\arguments{
  \item{data}{an object of class \code{"loci"}.}
  \item{INDICES}{a vector of the same length as the number of rows in \code{data}.}
  \item{FUN}{a function}
  \item{\dots}{(currently unused).}
  \item{simplify}{(currently unused).}
}
\details{
  The default \code{FUN = NULL} calculates allele frequencies for each
  population in \code{data}.
}
\value{
  a list by default indexed by locus.
}
\author{Emmanuel Paradis}
\seealso{
  \code{\link[base]{by}}
}
\examples{
data(jaguar)
by(jaguar)
by(na.omit(jaguar))
}
\keyword{manip}
