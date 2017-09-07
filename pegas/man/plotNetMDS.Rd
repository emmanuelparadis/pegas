\name{plotNetMDS}
\alias{plotNetMDS}
\title{Plot Networks With MDS Layout}
\description{
  This function plots a haplotype network using a layout calculated from
  an MDS performed on the pairwise distance matrix. The haplotypes have
  always the same positions for different networks.
}
\usage{
plotNetMDS(net, d, k = 2, show.mutation = FALSE, col = NULL, font = 2, cex = 1)
}
\arguments{
  \item{net}{an object of class \code{"haploNet"}.}
  \item{d}{an object of class \code{"dist"} (or a matrix).}
  \item{k}{the number of dimensions of the plot (2 or 3).}
  \item{show.mutation}{a logical value: if \code{TRUE}, the number of
    steps is printed on the links.}
  \item{col}{the colours of the links; by default, semi-transparent
    green.}
  \item{font}{the font used to print the labels; bold by default.}
  \item{cex}{the character expansion of the labels.}
}
\value{NULL}
\author{Emmanuel Paradis}
\references{
  Paradis, E. (2017) Analysis of haplotype networks: the randomized
  minimum spanning tree method. Manuscript.
}
\seealso{
  \code{\link{haploNet}}
}
\examples{
data(woodmouse)
d <- dist.dna(woodmouse, "n")
net <- rmst(d)
plotNetMDS(net, d)
}
\keyword{hplot}
