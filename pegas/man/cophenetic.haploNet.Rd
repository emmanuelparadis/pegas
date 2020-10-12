\name{cophenetic.haploNet}
\alias{cophenetic.haploNet}
\title{Cophenetic Matrix on Haplotype Networks}
\description{
  This function calculates the cophenetic distance on a network. The
  output can be used to find nodes with short distances to most nodes.
}
\usage{
\method{cophenetic}{haploNet}(x)
}
\arguments{
  \item{x}{an object of class \code{"haploNet"}.}
}
\details{
  The results of the function are likely to be approximate in most
  cases with reticulations in the network. In the case of MSTs, the
  results are exact.
}
\value{
  a numeric matrix with colnames and rownames set to the labels of the
  network nodes.
}
\author{Emmanuel Paradis}
\seealso{
  \code{\link[ape]{cophenetic.phylo}} in \pkg{ape},
  \code{\link[stats]{cophenetic}} for the generic function
}
\examples{
example(mst)
coph <- cophenetic(r)
rowSums(coph)
}
\keyword{manip}
