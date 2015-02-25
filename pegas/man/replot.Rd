\name{replot}
\alias{replot}
\title{Edit the Layout of a Haplotype Network}
\description{
  This function makes possible to change the layout of a haplotype
  network interactively or with specified coordinates.
}
\usage{
replot(xy = NULL, ...)
}
\arguments{
  \item{xy}{an optional list with vectors names \code{x} and \code{y}
    (or \code{xx} and \code{yy}) giving the coordinates of the nodes.}
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
}
\value{
  a named list with two numeric vertors (\code{x} and \code{y}).
}
\author{Emmanuel Paradis}
\seealso{
  \code{\link{haploNet}}, \code{\link{haploFreq}}
}
\examples{
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
