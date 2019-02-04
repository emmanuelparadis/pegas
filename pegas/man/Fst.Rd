\name{Fst}
\alias{Fst}
\alias{Rst}
\title{F-Statistics}
\description{
  \code{Fst} computes the \eqn{F_{IT}}{Fit}, \eqn{F_{ST}}{Fst} and
  \eqn{F_{IS}}{Fis} for each locus in the data. \code{Rst} computes the
  \eqn{R_{ST}}{Rst} for microsatellites.
}
\usage{
Fst(x, pop = NULL, quiet = TRUE, na.alleles = "")
Rst(x, pop = NULL, quiet = TRUE, na.alleles = "")
}
\arguments{
  \item{x}{an object of class \code{"loci"}.}
  \item{pop}{a vector or factor giving the population assignment of each
    row of \code{x}, or a single numeric value specifying which column
    of \code{x} to use as population indicator. By default, the column
    labelled \code{"population"} is used.}
  \item{quiet}{a logical value: should calculations be quiet?}
  \item{na.alleles}{by default, only genotypes coded as
    \link[base]{NA} are considered as missing data. This option
    is to specify if some alleles code for missing data.}
}
\details{
  \code{Fst} uses the formulae in Weir and Cockerham (1984) for each
  allele, and then averaged within each locus over the different alleles
  as suggested by these authors.

  \code{Rst} uses the formulae in Slatkin (1995).
}
\value{
  A matrix with genes (loci) as rows and the three \emph{F}-statistics
  as columns.
}
\references{
  Slatkin, M. (1995) A measure of population subdivision based on
  microsatellite allele frequencies.  \emph{Genetics}, \bold{139},
  457--462.

  Weir, B. S. and Cockerham, C. C. (1984) Estimating \emph{F}-statistics
  for the analysis of population structure. \emph{Evolution}, \bold{38},
  1358--1370.

  Weir, B. S. and Hill, W. G. (2002) Estimating \emph{F}-statistics.
  \emph{Annual Review of Genetics}, \bold{36}, 721--750.
}
\author{Emmanuel Paradis}
\seealso{
  \code{fstat} in package \pkg{hierfstat}; package \pkg{dirmult}
  on CRAN that implements various estimators of the
  Dirichlet-multinomial distribution, including maximum likekihood and
  the moments estimator of Weir and Hill (2002); \code{Fst} in
  \pkg{Biodem} that caculates \eqn{F_{ST}}{Fst} from a ``kinship
  matrix''.
}
\examples{
data(jaguar)
Fst(jaguar)
Rst(jaguar)
}
\keyword{htest}
