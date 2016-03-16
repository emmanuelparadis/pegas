\name{nuc.div}
\alias{nuc.div}
\alias{nuc.div.DNAbin}
\alias{nuc.div.haplotype}
\title{Nucleotide Diversity}
\description{
  This function computes the nucleotide diversity from a sample of DNA
  sequences or a set of haplotypes.
}
\usage{
nuc.div(x, ...)
\method{nuc.div}{DNAbin}(x, variance = FALSE, pairwise.deletion = FALSE, ...)
\method{nuc.div}{haplotype}(x, variance = FALSE, pairwise.deletion = FALSE, ...)
}
\arguments{
  \item{x}{a matrix or a list which contains the DNA sequences.}
  \item{variance}{a logical indicating whether to compute the variance
    of the estimated nucleotide diversity.}
  \item{pairwise.deletion}{a logical indicating whether to delete the
    sites with missing data in a pairwise way. The default is to delete
    the sites with at least one missing data for all sequences.}
  \item{\dots}{further arguments to be passed.}
}
\details{
  This is a generic function with methods for classes \code{"DNAbin"}
  and \code{"haplotype"}. The first method uses the sum of the number of
  differences between pairs of sequences divided by the number of
  comparisons (i.e. \eqn{n(n - 1)/2}, where \eqn{n} is the number of
  sequences). The second method uses haplotype frequencies. It could be
  that both methods give (slightly) different results because of missing
  or ambiguous nucleotides: this is generally solved by setting
  \code{pairwise.deletion = TRUE}.

  The variance of the estimated diversity uses formula (10.9) from Nei
  (1987). This applies only if all sequences are of the same lengths,
  and cannot be used if \code{pairwise.deletion = TRUE}. A bootstrap
  estimate may be in order if you insist on using the latter option.
}
\value{
  A numeric vector with one or two values if \code{variance = TRUE}.
}
\references{
  Nei, M. (1987) \emph{Molecular evolutionary genetics}. New York:
  Columbia University Press.
}
\author{Emmanuel Paradis}
\seealso{
  \code{\link[ape]{base.freq}}, \code{\link[ape]{GC.content}},
  \code{\link{theta.s}}, \code{\link[ape]{seg.sites}}
}
\examples{
data(woodmouse)
nuc.div(woodmouse)
nuc.div(woodmouse, TRUE)
nuc.div(woodmouse, FALSE, TRUE)
}
\keyword{manip}
\keyword{univar}
