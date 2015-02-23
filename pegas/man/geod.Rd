\name{geod}
\alias{geod}
\title{Geodesic Distances}
\description{
  This function calculates geodesic (or great-circle) distances
  between pairs of points with their longitudes and latitudes given in
  (decimal) degrees.
}
\usage{
geod(lon, lat = NULL, R = 6371)
}
\arguments{
  \item{lon}{either a vector of numeric values with the longitudes in
    degrees, or, if \code{lat = NULL}, a matrix giving the longitudes
    (first column) and the latitudes (second column).}
  \item{lat}{a vector with the latitudes.}
  \item{R}{the mean radius of the Earth (see details).}
}
\details{
  The default value of \code{R} is the mean radius of the Earth which is
  slightly smaller than the radius at the equator (6378.1 km).
}
\value{
  a numeric symmetric matrix with the distances between pairs of points
  in kilometres.
}
\author{Emmanuel Paradis}
\references{
  \url{http://en.wikipedia.org/wiki/Great-circle_distance}

  \url{http://en.wikipedia.org/wiki/Earth}
}
\seealso{geoTrans, as.dist}
\examples{
## the distance between 0N 0E and 0N 180E...
geod(c(0, 180), c(0, 0)) # ~ 20015.09 km
## ... the same using the radius of the Earth at the equator:
geod(c(0, 180), c(0, 0), 6378.1) # ~ 20037.39 km
## The same comparison for two points 5 degrees apart:
geod(c(0, 5), c(0, 0)) # ~ 555.9746 km
geod(c(0, 5), c(0, 0), 6378.1) # ~ 556.5942 km
}
\keyword{manip}
