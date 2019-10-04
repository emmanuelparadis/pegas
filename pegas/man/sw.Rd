\name{sw}
\alias{sw}
\alias{sw.default}
\alias{sw.DNAbin}
\alias{plot.sw}
\title{Sliding Windows}
\description{
  Applies a function over a matrix or a vector using sliding
  windows. \code{sw} is a generic function with a method for
  \code{"DNAbin"} matrices.
}
\usage{
sw(x, width, step, ...)
\method{sw}{default}(x, width = 100, step = 50, POS = NULL,
   FUN = mean, out.of.pos = NA_real_, na.rm = TRUE, L = NULL, ...)
\method{sw}{DNAbin}(x, width = 100, step = 50, FUN = GC.content,
   rowAverage = FALSE, quiet = TRUE, ...)
\method{plot}{sw}(x, type = "l", xlab = "Position", x.scaling = 1,
     show.ranges = FALSE, col.ranges = "blue",
     lty.ranges = 1, lwd.ranges = 1, ...)
}
\arguments{
  \item{x}{a vector or a matrix.}
  \item{width}{an integer giving the window width.}
  \item{step}{an integer giving the step separating successive windows.}
  \item{POS}{a numeric vector giving the positions of the sites.}
  \item{FUN}{the function to be applied to the windows.}
  \item{rowAverage}{a logical value: if \code{TRUE}, \code{FUN} is
    applied over all rows of \code{x}; if \code{FALSE} (the default)
    \code{FUN} is applied to each row of \code{x}.}
  \item{out.of.pos}{the values used for the sites which are not in
    \code{POS}.}
  \item{na.rm}{option passed to \code{FUN}.}
  \item{L}{the length of the chromosome (or sequence). If not given,
    this is largest value in \code{POS} or the length of \code{x} if
    \code{POS} is not given.}
  \item{quiet}{a logical value: if \code{FALSE}, the progress of the
    calculations is printed.}
  \item{type}{the type of plotting (see \code{\link[graphics]{plot.default}}).}
  \item{xlab}{the label under the x-axis.}
  \item{x.scaling}{the scaling of the x-axis.}
  \item{show.ranges}{a logical value specifying whether to show the
    ranges of the windows with horizontal segments (ignored with a
    warning if \code{x} is a matrix).}
  \item{col.ranges, lty.ranges, lwd.ranges}{arguments to modify the
    appearance of the above segments (see \code{\link[graphics]{segments}}).}
  \item{\dots}{further arguments passed to and from methods.}
}
\details{
  \code{FUN} should return a single value.

  \code{x} should be a matrix for the \code{"DNAbin"} method, or a
  vector for the default one.

  For the default method, the vector \code{x} is expanded into a vector
  of length \code{L} (see above on how this value is found) and the
  positions which are not in \code{POS} are filled with the value given
  in \code{out.of.pos}. The resulting vector is then analysed with the
  function \code{FUN} which must have an option \code{na.rm}. If the
  function you want to use does not have this option, you can use
  something like \code{FUN = function(x, na.rm = TRUE)
  foo(x[!is.na(x)])}, replacing `foo' by the name of your function. You
  may also include more control on the handling of missing data.
}
\value{
  a matrix or a vector (if \code{rowAverage = TRUE}).
}
\author{Emmanuel Paradis}
\examples{
data(woodmouse)
sw(woodmouse)
sw(woodmouse, 200, 200)
sw(woodmouse, 200, 200, rowAverage = TRUE)

## to get the proportions of G:
foo <- function(x) base.freq(x)["g"]
sw(woodmouse, 200, 200, FUN = foo, rowAverage = TRUE)

## a simulated example with the default method:
x <- runif(100)
pos <- sort(sample(1e6, 100))
resx <- sw(x, w = 2e4, s = 5e3, POS = pos, L = 1e6)
plot(resx, show.ranges = TRUE, x.scaling = 1e6, xlab = "Position (Mb)")
}
\keyword{manip}
