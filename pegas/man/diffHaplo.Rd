\name{diffHaplo}
\alias{diffHaplo}
\title{Comparison Between Two Haplotypes}
\description{
  This function compares two haplotypes and returns a summary of the
  differences.
}
\usage{
diffHaplo(h, a = 1, b = 2)
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
