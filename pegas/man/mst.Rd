\name{mst}
\alias{mst}
\description{
  This function computes a minimum spanning tree using Kruskal's
  algorithm.
}
\title{Minimum Spanning Tree}
\usage{
mst(d)
}
\arguments{
  \item{d}{a distance matrix, either as an object of class \code{"dist"},
    or a (square symmetric) matrix.}
}
\value{
  an object of class \code{"\link{haploNet}"}.
}
\references{
  Kruskal, J. B., Jr. (1956) On the shortest spanning subtree of a graph
  and the traveling salesman problem. \emph{Proceedings of the American
    Mathematical Society}, \bold{7}, 48--50.
}
\author{Emmanuel Paradis}
\seealso{
  \code{\link{haploNet}}
}
\examples{
data(woodmouse)
d <- dist.dna(woodmouse, "n")
(r <- mst(d))
plot(r)
}
\keyword{models}
