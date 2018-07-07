\name{sw}
\alias{sw}
\alias{sw.DNAbin}
\title{Sliding Windows}
\description{
  Applies a function over a matrix using sliding windows. \code{sw} is a
  generic function with a method for \code{"DNAbin"} matrices.
}
\usage{
sw(x, width, step, ...)
\method{sw}{DNAbin}(x, width = 100, step = 50, FUN = GC.content,
   rowAverage = FALSE, quiet = TRUE, ...)
}
\arguments{
  \item{x}{a matrix.}
  \item{width}{an integer giving the window width.}
  \item{step}{an integer giving the step separating successive windows.}
  \item{FUN}{the function to be applied to the windows.}
  \item{rowAverage}{a logical value: if \code{TRUE}, \code{FUN} is
    applied over all rows of \code{x}; if \code{FALSE} (the default)
    \code{FUN} is applied to each row of \code{x}.}
  \item{quiet}{a logical value: if \code{FALSE}, the progress of the
    calculations is printed.}
  \item{\dots}{further arguments passed to and from methods.}
}
\details{
  \code{FUN} should return a single value, and \code{x} should be a matrix.
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
}
\keyword{manip}
