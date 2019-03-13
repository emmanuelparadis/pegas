\name{mst}
\alias{mst}
\alias{msn}
\alias{rmst}
\title{Minimum Spanning Tree and Network}
\description{
  Computes a minimum spanning tree using Kruskal's algorithm, the
  minimum spanning network using Bandelt et al.'s algorithm, or the
  randomized minimum spanning tree (Paradis 2017).
}
\usage{
mst(d)
msn(d)
rmst(d, B = 100)
}
\arguments{
  \item{d}{a distance matrix, either as an object of class \code{"dist"},
    or a (square symmetric) matrix.}
  \item{B}{number of randomizations.}
}
\value{
  an object of class \code{"\link{haploNet}"}.
}
\references{
  Bandelt, H. J., Forster, P. and Rohl, A. (1999) Median-joining networks
  for inferring intraspecific phylogenies. \emph{Molecular Biology and
  Evolution}, \bold{16}, 37--48.

  Kruskal, J. B., Jr. (1956) On the shortest spanning subtree of a graph
  and the traveling salesman problem. \emph{Proceedings of the American
    Mathematical Society}, \bold{7}, 48--50.

  Paradis, E. (2018) Analysis of haplotype networks: the randomized
  minimum spanning tree method. \emph{Methods in Ecology and Evolution},
  \bold{9}, 000--000. DOI: 10.1111/2041-210X.12969.
}
\author{Emmanuel Paradis}
\seealso{
  \code{\link{haploNet}}, \code{\link{mjn}}
}
\note{
  \pkg{ape} has a function named \code{mst} which is older (and used by
  other packages) and returns its results in a different form. The
  present version is more efficient. If you want to use the older
  version after loading \pkg{pegas}, use \code{ape::mst} since \pkg{ape}
  will certainly always be loaded before \pkg{pegas}.
}
\examples{
data(woodmouse)
d <- dist.dna(woodmouse, "n")
(r <- mst(d))
plot(r)

## a case where the RMST and the MJN are identical:
x <- c(">A", "TAAGTGCAT", ">B", "TAAATGCAT", ">C", "TAGGTGCAT", ">D", "TAAGTACAT",
       ">E", "TAAGTGTAT", ">F", "TAAGTACAC", ">G", "TAAGTACGT", ">H", "CAAGTACAC",
       ">I", "CAAGCACAC", ">J", "CAAGTACAT", ">K", "CGAGTACAT", ">L", "TAAGTACGC",
       ">M", "CAAGCACAT")
fl <- tempfile()
cat(x, file = fl, sep = "\n")
x <- read.dna(fl, "f")
tr <- rmst(dist.dna(x, "n"))
ts <- mjn(x)
stopifnot(all.equal(tr, ts))
unlink(fl)
}
\keyword{models}
