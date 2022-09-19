\name{haploNet}
\alias{haploNet}
\alias{print.haploNet}
\alias{plot.haploNet}
\alias{as.network.haploNet}
\alias{as.igraph.haploNet}
\alias{as.phylo.haploNet}
\alias{as.evonet.haploNet}
\title{Haplotype Networks}
\description{
  \code{haploNet} computes a haplotype network. There is a plot method
  and two conversion functions towards other packages.
}
\usage{
haploNet(h, d = NULL, getProb = TRUE)
\method{print}{haploNet}(x, ...)
\method{plot}{haploNet}(x, size = 1, col, bg, col.link, lwd, lty,
     shape = "circles", pie = NULL, labels, font, cex, scale.ratio,
     asp = 1, legend = FALSE, fast = FALSE, show.mutation,
     threshold = c(1, 2), xy = NULL, ...)
\method{as.network}{haploNet}(x, directed = FALSE, altlinks = TRUE, ...)
\method{as.igraph}{haploNet}(x, directed = FALSE, use.labels = TRUE,
        altlinks = TRUE, ...)
\method{as.phylo}{haploNet}(x, quiet, ...)
\method{as.evonet}{haploNet}(x, ...)
}
\arguments{
  \item{h}{an object of class \code{"haplotype"}.}
  \item{d}{an object giving the distances among haplotypes (see
    details).}
  \item{getProb}{a logical specifying whether to calculate Templeton's
    probabilities (see details).}
  \item{x}{an object of class \code{"haploNet"}.}
  \item{size}{a numeric vector giving the diameter of the circles
    representing the haplotypes: this is in the same unit than the links
    and eventually recycled.}
  \item{col}{a character vector specifying the colours of the circles;
    eventually recycled.}
  \item{bg}{a character vector (or a function) specifying either the
    colours of the background of the symbols (if \code{pie = NULL}), or
    the colours of the slices of the pies (could be a function);
    eventually recycled.}
  \item{col.link}{a character vector specifying the colours of the links;
    eventually recycled.}
  \item{lwd}{a numeric vector giving the width of the links; eventually
    recycled.}
  \item{lty}{idem for the line types.}
  \item{shape}{the symbol shape used for the haplotypes (eventually
    recycled): \code{"circles"}, \code{"squares"}, \code{"diamonds"}
    (can be abbreviated).}
  \item{pie}{a matrix used to draw pie charts for each haplotype; its
    number of rows must be equal to the number of haplotypes.}
  \item{labels}{a logical specifying whether to identify the haplotypes
    with their labels (the default).}
  \item{font}{the font used for these labels (bold by default); must be
    an integer between 1 and 4.}
  \item{cex}{a numerical specifying the character expansion of the labels.}
  \item{scale.ratio}{the ratio of the scale of the links representing
    the number of steps on the scale of the circles representing the
    haplotypes. It may be needed to give a value greater than one to
    avoid overlapping circles.}
  \item{asp}{the aspect ratio of the plot. Do not change the default
    unless you want to distort your network.}
  \item{legend}{a logical specifying whether to draw the legend, or a
    vector of length two giving the coordinates where to draw the
    legend; \code{FALSE} by default. If \code{TRUE}, the user is asked
    to click where to draw the legend.}
  \item{fast}{a logical specifying whether to optimize the spacing of
    the circles; \code{FALSE} by default.}
  \item{show.mutation}{an integer value: if 0, nothing is drawn on the
    links; if 1, the mutations are shown with small segments on the
    links; if 2, they are shown with small dots; if 3, the number of
    mutations are printed on the links.}
  \item{threshold}{a numeric vector with two values (or 0) giving the
    lower and upper numbers of mutations for alternative links to be
    displayed. If \code{threshold = 0}, alternative links are not drawn
    at all.}
  \item{directed}{a logical specifying whether the network is directed
    (\code{FALSE} by default).}
  \item{use.labels}{a logical specifying whether to use the original
    labels in the returned network.}
  \item{altlinks}{whether to output the alternative links when
    converting to another class; \code{TRUE} by default.}
  \item{quiet}{whether to give a warning when reticulations are dropped
    when converting a network into a tree.}
  \item{xy}{the coordinates of the nodes (see \code{\link{replot}}).}
  \item{\dots}{further arguments passed to \code{plot}.}
}
\details{
  By default, the haplotype network is built using an infinite site
  model (i.e., uncorrected or Hamming distance) of DNA sequences and
  pairwise deletion of missing data (see \code{\link[ape]{dist.dna}}).
  Users may specify their own distance with the argument \code{d}. There
  is no check of labels, so the user must make sure that the distances
  are ordered in the same way than the haplotypes.

  The probabilities calculated with Templeton et al.'s (1992) method may
  give non-finite values with very divergent sequences, resulting in an
  error from \code{haploNet}. If this happens, it may be better to use
  \code{getProb = FALSE}.

  If two haplotypes are very different, \code{haploNet} will likely fail
  (error during integration due to non-finite values).
}
\note{
  Plotting haplotype networks is a difficult task. There is a vignette
  in \pkg{pegas} (see \code{vignette("PlotHaploNet")}) giving some
  information on this isseu. You may also see two posts on r-sig-genetics
  (July 2022) that give some tricks in the situation when one haplotype
  is abundant and the others are in low frequencies (the symbols are
  likely to overlap a lot by default):

  \url{https://stat.ethz.ch/pipermail/r-sig-genetics/2022-July/000237.html}

  \url{https://stat.ethz.ch/pipermail/r-sig-genetics/2022-July/000238.html}

  The first post explains how to use the package \pkg{network} in
  combination with \pkg{pegas}, and the second one gives a trick that
  works with \pkg{pegas} only for a similar result.
}
\value{
  \code{haploNet} returns an object of class \code{"haploNet"} which is
  a matrix where each row represents a link in the network, the first
  and second columns give the numbers of the linked haplotypes, the
  third column, named \code{"step"}, gives the number of steps in this
  link, and the fourth column, named \code{"Prob"}, gives the
  probability of a parsimonious link as given by Templeton et
  al. (1992). There are three additional attributes: \code{"freq"}, the
  absolute frequencies of each haplotype, \code{"labels"}, their labels,
  and \code{"alter.links"}, the alternative links of the network.

  \code{as.network} and \code{as.igraph} return objects of the
  appropriate class.
}
\author{Emmanuel Paradis, Klaus Schliep}
\references{
  Templeton, A. R., Crandall, K. A. and Sing, C. F. (1992) A cladistic
  analysis of phenotypic association with haplotypes inferred from
  restriction endonuclease mapping and DNA sequence data. III. Cladogram
  estimation. \emph{Genetics}, \bold{132}, 619--635.
}
\seealso{
  \code{\link{haplotype}}, \code{\link{haploFreq}}, \code{\link{replot}},
  \code{\link{diffHaplo}}, \code{\link{mst}}, \code{\link{mjn}}
}
\examples{
## generate some artificial data from 'woodmouse':
data(woodmouse)
x <- woodmouse[sample(15, size = 110, replace = TRUE), ]
h <- haplotype(x)
(net <- haploNet(h))
plot(net)
## symbol sizes equal to haplotype sizes:
plot(net, size = attr(net, "freq"), fast = TRUE)
plot(net, size = attr(net, "freq"))
plot(net, size = attr(net, "freq"), scale.ratio = 2, cex = 0.8)
}
\keyword{models}
\keyword{hplot}
