\name{dist.asd}
\alias{dist.asd}
\alias{dist.snp}
\title{Allelic Sharing Distance}
\description{
  This function computes the allelic sharing distance (ASD) for diploid
  genotypes.
}
\usage{
dist.asd(x, scaled = TRUE)
dist.snp(x, scaled = TRUE)
}
\arguments{
  \item{x}{an object of class \code{"loci"}.}
  \item{scaled}{a logical value specifying whether the distances should
    be scaled by the number of loci.}
}
\details{
  The ASD between two diploid genotypes is (Gao and Martin, 2009):

  \deqn{\frac{1}{L}\sum_{j=1}^L d_j}{(1/L)sum d_j, j = 1, ... L}

  where \eqn{L} is the number loci, \eqn{d_j} is the value for the
  \eqn{j}th locus: 0 if both genotypes are identical, 1 if they have one
  allele in common, or 2 if they have no allele in common.

  \code{dist.asd} works for all diploid genotypes (phased or unphased,
  with two alleles or more). Note that the required conditions are not
  checked by the present function: see the functions below.
}
\value{
  an object of class \code{"dist"}.
}
\note{
  \code{dist.snp} will be removed soon; it currently calls
  \code{dist.asd} with a warning.
}
\references{
  Gao, X. and Martin, E. R. (2009) Using allele sharing distance for
  detecting human population stratification. \emph{Human Hederity},
  \bold{68}, 182--191.
}
\author{Emmanuel Paradis}
\seealso{
  \code{\link{is.snp}}, \code{\link{is.phased}}, \code{\link{getPloidy}}
}
\examples{
data(jaguar)
## ASD for micro-satellites:
d <- dist.asd(jaguar)
co <- rainbow(nlevels(jaguar$pop))
plot(nj(d), "u", tip.color = co[jaguar$pop], font = 2, lab4 = "a")
legend("topleft", legend = levels(jaguar$pop), text.col = co, text.font = 2)
}
\keyword{manip}
\keyword{multivariate}
\keyword{cluster}
