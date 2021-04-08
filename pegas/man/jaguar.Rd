\name{jaguar}
\alias{jaguar}
\title{Jaguar Micro-Satellites}
\description{
  Fifty nine jaguars (\emph{Panthera onca}) from four populations
  genotyped at thirteen micro-satellites by Haag et al. (2010).
}
\usage{
data(jaguar)
}
\format{
  An object of class \code{"loci"} with 59 rows and 14 columns.
}
\source{
  Haag, T., Santos, A. S., Sana, D. A., Morato, R. G., Cullen, Jr., L.,
  Crawshaw, Jr., P. G., De Angelo, C., Di Bitetti, M. S., Salzano,
  F. M. and Eizirik, E. (2010) The effect of habitat fragmentation on
  the genetic structure of a top predator: loss of diversity and high
  differentiation among remnant populations of Atlantic Forest jaguars
  (\emph{Panthera onca}). \emph{Molecular Ecology}, \bold{22},
  4906--4921.

  Haag, T., Santos, A. S., Sana, D. A., Morato, R. G., Cullen, Jr., L.,
  Crawshaw, Jr., P. G., De Angelo, C., Di Bitetti, M. S., Salzano,
  F. M. and Eizirik, E. (2010) Data from: The effect of habitat
  fragmentation on the genetic structure of a top predator: loss of
  diversity and high differentiation among remnant populations of
  Atlantic Forest jaguars (Panthera onca). Dryad Digital
  Repository. \doi{10.5061/dryad.1884}
}
\seealso{
  \code{\link{loci}}, \code{\link{alleles2loci}}

  The vignette \dQuote{ReadingFiles} explains how to read data like
  these from Dryad (\url{https://datadryad.org/stash}).
}
\examples{
data(jaguar)
str(jaguar)
s <- summary(jaguar)
\dontrun{
## works if the device is large enough:
plot(s, layout = 30, las = 2)
}}
\keyword{datasets}
