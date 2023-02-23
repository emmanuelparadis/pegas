\name{replot}
\alias{replot}
\title{Edit the Layout of a Haplotype Network}
\description{
  This function makes possible to change the layout of a haplotype
  network interactively or with specified coordinates.
}
\usage{
replot(xy = NULL, col.identifier = "purple", ...)
}
\arguments{
  \item{xy}{an optional list with vectors names \code{x} and \code{y}
    (or \code{xx} and \code{yy}) giving the coordinates of the nodes.}
  \item{col.identifier}{the colour used to identify the node to be moved.}
  \item{\dots}{further arguments passed to \code{plot}.}
}
\details{
  This function can be used in two ways. By default (i.e.,
  \code{replot()}), the user can edit a plotted haplotype network by
  clicking with the mouse on the graphical window: a message is printed
  asking to click once close to the node to move and then clicking again
  where this node should be placed (careful: two separate single
  clicks). Editing is stopped with a right click.

  The second possible use is to specify the new coordinates of the nodes
  with the argument \code{xy}, typically, from a previous call to
  \code{replot} (see examples).

  Since \pkg{pegas} 1.0, these coordinates can be used directly in
  \code{\link{plot.haploNet}} making possible to combine networks with
  other graphics (which not possible with \code{replot} because the
  network is replotted).
}
\note{
  For users of RStudio: the function does not work within this
  application. It seems the best is to run R from a shell (or maybe
  opening a new graphical device with \code{\link[grDevices]{X11}}).
}
\value{
  a named list with two numeric vertors (\code{x} and \code{y}).
}
\author{Emmanuel Paradis}
\seealso{
  \code{\link{haploNet}}, \code{\link{haploFreq}}
}
\examples{
## a non-interactive example:
example(mjn)
layout(matrix(1:2, 1))
plot(ntz, labels = TRUE)
## it is possible plot this network with no line-crossing
## with these coordinates:
xy <- list(x = c(3.2, -2.6, -6.6, -7.2, 0, 3.5, 2.6, -2.9, -0.3, 3.4, -3.4),
           y = c(3.4, 4.4, 1.3, -3.9, -5.5, -10.9, 0.1, -0.8, -2.3,  -7.9, -8.1))
replot(ntz, xy = xy) # or plot(ntz, xy = xy, labels = TRUE)
layout(1)

## an interactive example:
\dontrun{
data(woodmouse)
net <- haploNet(haplotype(woodmouse))
plot(net)
o <- replot() # interactive
## click to rearrange the network at will...
## then do a different plot using the same coordinates:
plot(net, bg = "red", labels = FALSE, show.mutation = 2)
replot(o) # not interactive
}}
\keyword{hplot}
