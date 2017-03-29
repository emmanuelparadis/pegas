\name{dist.hamming}
\alias{dist.hamming}
\title{Hamming Distance}
\description{
  This function implements a general purpose Hamming distance.
}
\usage{
dist.hamming(x)
}
\arguments{
  \item{x}{a matrix or a data frame.}
}
\details{
  This function should work for a wide range of data types. A typical
  usage would be with an object of class \code{c("haplotype",
    "character")}.

  For objects of class \code{c("haplotype", "DNAbin")}, it is better to
  use \code{dist.dna(x, "n")} to compute the Hamming distances.
}
\value{
  an object of class \code{"dist"}.
}
\author{Emmanuel Paradis}
\seealso{
  \code{\link{haplotype}}, \code{\link{dist.haplotype.loci}}
}
\keyword{model}
