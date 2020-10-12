\name{diffHaplo}
\alias{diffHaplo}
\title{Comparison Between Two Haplotypes}
\description{
  This function compares two haplotypes and returns a summary of the
  differences.
}
\usage{
diffHaplo(h, a = 1, b = 2, strict = FALSE, trailingGapsAsN = TRUE)
}
\arguments{
  \item{h}{an object of class \code{"haplotype"}.}
  \item{a, b}{two integers (or character strings) giving the indices (or
    labels) of the two haplotypes to be compared.}
  \item{strict}{a logical value; if \code{TRUE}, ambiguities and gaps in
    the sequences are ignored and treated as separate characters.}
  \item{trailingGapsAsN}{a logical value; if \code{TRUE} (the default),
    the leading and trailing alignment gaps are considered as unknown
    bases (i.e., N). This option has no effect if \code{strict = TRUE}.}
}
\details{
  The options \code{strict} and \code{trailingGapsAsN} are passed to
  \code{\link[ape]{seg.sites}}.
}
\value{
  a data frame with three columns named \code{pos} (position of the
  differences) and the labels of the two haplotypes compared.
}
\author{Emmanuel Paradis}
\seealso{
  \code{\link{haploNet}}, \code{\link{haplotype}}
}
\examples{
data(woodmouse)
h <- haplotype(woodmouse)
diffHaplo(h) # compares the 1st and 2nd haplotypes
diffHaplo(h, 1, 3)
diffHaplo(h, "I", "III") # same than above but using labels
}
\keyword{manip}
