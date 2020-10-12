\name{mutations}
\alias{mutations}
\title{Plot Mutations on Networks}
\description{
  \code{mutations} draws annotations about mutations related to the link
  of a haplotype network.
}
\usage{
mutations(haploNet, link, x, y, data = NULL, style = "table", POS, SEQLEN, ...)
}
\arguments{
  \item{haploNet}{an object of class \code{"haploNet"} which should be
    plotted beforehand.}
  \item{link}{the link number; can be left missing in which case the
    list of links in the network is printed and the function exits.}
  \item{x, y}{the coordinates where to draw the annotations; can be left
    missing: the user is then asked to click where to draw them and the
    chosen coordinates are printed.}
  \item{data}{the sequence data; can be left missing if the data are
    attached to the network (for a MJN network output by
    \code{\link{mjn}}.)}
  \item{style}{the type annotations. There two possible choices:
    \code{"table"} (default) and \code{"sequence"} (can be abbreviated).}
  \item{POS, SEQLEN}{a vector of genomic positions and the sequence
    length in case \code{data} is of class \code{"haplotype.loci"}.}
  \item{\dots}{options}
}
\details{
  The easiest way to use this function is with an output from
  \code{\link{mjn}} since the data are attached to the network. In other
  cases, the sequence data must given to the argument \code{data} or
  attached to the network as an attribute named \code{"data"}.
}
\value{none}
\author{Emmanuel Paradis}
\seealso{\code{\link{plot.haploNet}}, \code{\link{haplotype.loci}}}
\examples{
## simple example
x <- as.DNAbin(matrix(c("a", "g"), 2, 1))
rownames(x) <- paste("Ind", 1:2, sep = "_")
nt <- mst(dist.dna(x, "N"))
plot(nt)
mutations(nt, link = 1, x = 2, y = 2, data = x)

example(mjn)
plot(ntz, xlim = c(-5, 20))
mutations(ntz, 6, 10, 0, style = "s")
mutations(ntz, 8, 10, -2, style = "s")
}
\keyword{lplot}