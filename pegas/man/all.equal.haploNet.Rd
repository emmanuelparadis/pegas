\name{all.equal.haploNet}
\alias{all.equal.haploNet}
\title{Compare Two Haplotype Networks}
\description{
  This function compares two haplotype networks and returns either
  \code{TRUE} or a description of the differences.
}
\usage{
\method{all.equal}{haploNet}(target, current, use.steps = TRUE, ...)
}
\arguments{
  \item{target, current}{two objects of class \code{"haplotype"}.}
  \item{use.steps}{a logical value: whether to consider the number of
    steps (or length) in each link.}
  \item{\dots}{(unused).}
}
\details{
  This function should return \code{TRUE} if the two networks are
  identical even if the links are ordered differently. In all other
  situations, a vector of character strings describing the differences
  is returned.

  As usual with the \code{\link[base]{all.equal}} function, this cannot
  be used directly to return a \code{TRUE}/\code{FALSE} value (see
  examples).
}
\value{
  either a logical value (\code{TRUE}), or a vector of mode character.
}
\author{Emmanuel Paradis}
\seealso{\code{\link{haploNet}}, \code{\link{mst}}}
\examples{
data(woodmouse)
d <- dist.dna(woodmouse, "n")
nt1 <- mst(d)
nt2 <- msn(d)
(comp <- all.equal(nt1, nt2)) # clearly different

## how to use all.equal to return TRUE/FALSE:
isTRUE(comp) # FALSE
}
\keyword{models}
