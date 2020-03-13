\encoding{utf8}
\name{geoTrans}
\alias{geoTrans}
\alias{geoTrans2}
\title{Manipulate Geographical Coordinates}
\description{
  \code{geoTrans} transforms geographical coordinates in degrees,
  minutes and seconds input as characters (or a factor) into numerical
  values in degrees. \code{geoTrans2} does the reverse operation.
}
\usage{
geoTrans(x, degsym = NULL, minsym = "'", secsym = "\"")
geoTrans2(lon, lat = NULL, degsym = NULL, minsym = "'",
          secsym = "\"", dropzero = FALSE, digits = 3,
          latex = FALSE)
}
\arguments{
  \item{x}{a vector of character strings storing geographical
    coordinates; this can be a factor with the levels correctly set.}
  \item{degsym, minsym, secsym}{a single character giving the symbol
    used for degrees, minutes and seconds, respectively.}
  \item{lon}{either a vector of numeric values with the longitudes in
    degrees, or, if \code{lat = NULL}, a matrix (or a data frame) giving
    the longitudes in the first column and the latitudes in the second
    column.}
  \item{lat}{a vector with the latitudes.}
  \item{dropzero}{a logical value: if \code{TRUE}, the number of
    arc-seconds is dropped if it is zero; similarly for the number of
    arc-minutes if the number of arc-seconds is also zero.}
  \item{digits}{an integer used for rounding the number of arc-seconds.}
  \item{latex}{a logical value: if \code{TRUE}, the returned character
    is formatted with LaTeX code.}
}
\details{
  \code{geoTrans} should be robust to any pattern of spacing around the
  values and the symbols (see examples). If the letter S, W, or O is
  found is the coordinate, the returned value is negative. Note that
  longitude and latitude should not be mixed in the same character
  strings.

  \code{geoTrans2} can be used with \code{\link[base]{cat}} (see
  examples).

  The default for \code{degsym} (\code{NULL}) is because the degree
  symbol (°) is coded differently in different character encodings.
  By default, the function will use the appropriate character depending
  on the system and encoding used.
}
\value{
  \code{geoTrans} returns a numeric vector with the coordinates in
  degrees (eventually as decimal values). \code{geoTrans2} returns a
  character vector.
}
\author{Emmanuel Paradis}
\seealso{\code{\link{geod}}}
\examples{
coord <- c("N 43°27'30\"", "N43°27'30\"", "43°27'30\"N",
           "43° 27' 30\" N", "43 ° 27 ' 30 \" N",
           "43°27'30\"", "43°27.5'")
cat(coord, sep = "\n")
geoTrans(coord)
geoTrans("43 D 27.5'", degsym = "D")
geoTrans("43° 27' 30\" S")

XL <- c(100.6417, 102.9500)
YL <- c(11.55833, 14.51667)
cat(geoTrans2(XL, YL, dropzero = TRUE), sep = "\n")
cat(geoTrans2(XL, YL, latex = TRUE), sep = "\\\\\n")
}
\keyword{manip}
