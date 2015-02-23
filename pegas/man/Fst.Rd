\name{Fst}
\alias{Fst}
\title{F-Statistics}
\description{
  This function computes the \eqn{F_{IT}}{Fit}, \eqn{F_{ST}}{Fst} and
  \eqn{F_{IS}}{Fis} for each locus in the data.
}
\usage{
Fst(x, pop = NULL)
}
\arguments{
  \item{x}{an object of class \code{"loci"}.}
  \item{pop}{a vector or factor giving the population assignment of each
    row of \code{x}, or a single numeric value specifying which column
    of \code{x} to use as population indicator. By default, the column
    labelled \code{"population"} is used.}
}
\details{
  The formulae in Weir and Cockerham (1984) are used for each allele,
  and then averaged within each locus over the different alleles as
  suggested by these authors.
}
\value{
  A matrix with genes (loci) as rows and the three \emph{F}-statistics
  as columns.
}
\note{
  Programming bugs have been fixed in version 0.3-2 of
  \pkg{pegas}. Further tests and feed-back are still welcome.
}
\references{
  Weir, B. S. and Cockerham, C. C. (1984) Estimating \emph{F}-statistics
  for the analysis of population structure. \emph{Evolution}, \bold{38},
  1358--1370.

  Weir, B. S. and Hill, W. G. (2002) Estimating \emph{F}-statistics.
  \emph{Annu Review of Genetics}, \bold{36}, 721--750.
}
\author{Emmanuel Paradis}
\seealso{
  \code{\link[adegenet]{fstat}} in \pkg{adegenet}; package \pkg{dirmult}
  on CRAN that implements various estimators of the
  Dirichlet-multinomial distribution, including maximum likekihood and
  the moments estimator of Weir and Hill (2002); \code{Fst} in
  \pkg{Biodem} that caculates \eqn{F_{ST}}{Fst} from a ``kinship
  matrix''.
}
\examples{
require(adegenet)
data(nancycats)
x <- as.loci(nancycats)
Fst(x)
}
\keyword{htest}
