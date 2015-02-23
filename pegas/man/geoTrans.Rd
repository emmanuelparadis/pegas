\encoding{utf8}
\name{geoTrans}
\alias{geoTrans}
\title{Manipulate Geographic Coordinates}
\description{
  This function transforms standard geographical coordinates in degrees,
  minutes and seconds input as characters (or a factor) into numerical
  values in degrees.
}
\usage{
geoTrans(x, degsym = NULL, minsym = "'", secsym = "\"")
}
\arguments{
  \item{x}{a vector of character strings storing geographical
    coordinates; this can be a factor with the levels correctly set.}
  \item{degsym, minsym, secsym}{a single character giving the symbol
    used for degrees, minutes and seconds, respectively.}
}
\details{
  This function should be robust to any pattern of spacing around the
  values and the symbols (see examples).

  If the letter S, W, or O is found is the coordinate, the returned
  value is negative.

  Note that longitude and latitude should not be mixed in the same
  character strings.

  The default for \code{degsym} (\code{NULL}) is because the degree
  symbol (°) is coded differently in different character encodings.
  By default, the function will use the appropriate character depending
  on the system and encoding used.
}
\value{
  a numeric vector with the coordinates in degrees (eventually as
  decimal values).
}
\author{Emmanuel Paradis}
\seealso{geod}
\examples{
coord <- c("N 43°27'30\"", "N43°27'30\"", "43°27'30\"N",
           "43° 27' 30\" N", "43 ° 27 ' 30 \" N",
           "43°27'30\"", "43°27.5'")
cat(coord, sep = "\n")
geoTrans(coord)
geoTrans("43 D 27.5'", degsym = "D")
geoTrans("43° 27' 30\" S")
}
\keyword{manip}
