\name{diffHaplotype}
\alias{diffHaplotype}
\title{ a Haplotype}
\description{
  This function compares two haplotypes and returns a summary of the
  differences.
}
\usage{
diffHaplotype(h, a = 1, b = 2)
}
\arguments{
  \item{h}{an object of class \code{"haplotype"}.}
  \item{a, b}{two integers (or character strings) giving the indices (or
    labels) of the two haplotypes to be compared.}
}
\details{
  This function prints the number of transitions and transversions
  between both sequences, and returns a data frame with three columns
  giving the positions of the differences and the nucleotides in each
  sequence at these positions.
}
\value{
  a data frame with three columns named \code{pos} and with the labels
  of the two haplotypes compared.
}
\author{Emmanuel Paradis}
\seealso{
  \code{\link{haploNet}}, \code{\link{haplotype}}
}
\examples{
data(woodmouse)
h <- haplotype(woodmouse)
diffHaplotype(h) # compares the 1st and 2nd haplotypes
diffHaplotype(h, 14, 15)
diffHaplotype(h, "XIV", "XV") # same than above but using labels
}
\keyword{manip}
