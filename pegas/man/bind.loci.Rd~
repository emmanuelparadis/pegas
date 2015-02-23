\name{bind.loci}
\alias{bind.loci}
\alias{rbind.loci}
\alias{cbind.loci}
\title{Bind Loci Objects}
\description{
  These functions combine objects of class \code{"loci"} by binding
  their rows or their columns.
}
\usage{
\method{rbind}{loci}(...)
\method{cbind}{loci}(...)
}
\arguments{
  \item{\dots}{some object(s) of class \code{"loci"}, separated with
    commas.}
}
\details{
  These two methods call \code{[rc]bind.data.frame} and take care to
  respect the attribute ``locicol'' of the returned object.

  You can pass a data frame in the \code{\dots}, but then you should
  bypass the generic by calling \code{cbind.loci} directly. Do not try
  to pass a vector: this will mess the ``locicol'' attribute. Instead,
  make a data frame with this vector (see examples).
}
\value{
  An object of class \code{"loci"}.
}
\author{Emmanuel Paradis}
\seealso{[.loci}
\examples{
a <- as.loci(data.frame(x = "A/a", y = 1), col.loci = 1)
b <- as.loci(data.frame(y = 2, x = "A/A"), col.loci = 2)
## rbind.loci reorders the columns if necessary:
str(rbind(a, b))
## cbind sets "locicol" correctly:
str(cbind(a, b))
str(cbind(b, a))
## Unexpected result...
str(cbind(a, data.frame(z = 10)))
## ... bypass the generic:
str(cbind.loci(a, data.frame(z = 10)))
## Here "locicol" is not correct...
str(cbind.loci(z = 10, a))
## ... instead
str(cbind.loci(data.frame(z = 10), a))
}
\keyword{manip}
