\name{theta.h}
\alias{theta.h}
\title{Population Parameter THETA using Homozygosity}
\usage{
theta.h(x, standard.error = FALSE)
}
\arguments{
  \item{x}{a vector or a factor.}
  \item{standard.error}{a logical indicating whether the standard error
    of the estimated theta should be returned (\code{TRUE}), the default
    being \code{FALSE}.}
}
\description{
  This function computes the population parameter THETA using the
  homozygosity (or mean heterozygosity) from gene frequencies.
}
\value{
  A numeric vector of length one with the estimated theta (the default),
  or of length two if the standard error is returned
  (\code{standard.error = TRUE}).
}
\details{
  The argument \code{x} can be either a factor or a vector. If it is a
  factor, then it is taken to give the individual alleles in the
  population. If it is a numeric vector, then its values are taken to be
  the numbers of each allele in the population. If it is a non-numeric
  vector, it is a coerced as a factor.

  The standard error is computed with an approximation due to
  Chakraborty and Weiss (1991).
}
\references{
  Zouros, E. (1979) Mutation rates, population sizes and amounts of
  electrophoretic variation at enzyme loci in natural
  populations. \emph{Genetics}, \bold{92}, 623--646.

  Chakraborty, R. and Weiss, K. M. (1991) Genetic variation of the
  mitochondrial DNA genome in American Indians is at mutation-drift
  equilibrium. \emph{American Journal of Physical Anthropology},
  \bold{86}, 497--506.
}
\author{Emmanuel Paradis}
\seealso{
  \code{\link{heterozygosity}}, \code{\link{theta.s}},
  \code{\link{theta.k}}, \code{\link{theta.tree}}
}
\examples{
## similar to what is in ?H:
data(jaguar)
## convert the data and compute frequencies:
S <- summary(jaguar)
## compute THETA for all loci:
sapply(S, function(x) theta.h(x$allele))
}
\keyword{manip}
\keyword{univar}
