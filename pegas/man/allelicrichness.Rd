\name{allelicrichness}
\alias{allelicrichness}
\alias{rarefactionplot}
\alias{rhost}
\title{Allelic Richness and Rarefaction Plots}
\description{
  These functions analyse allelic richness.
}
\usage{
allelicrichness(x, pop = NULL, method = "extrapolation", min.n = NULL)
rarefactionplot(x, maxn = nrow(x), type = "l", xlab = "Sample size",
                ylab = "Expected number of alleles", plot = TRUE, ...)
rhost(x, pop = NULL, method = "extrapolation")
}
\arguments{
  \item{x}{an object of class \code{"loci"}.}
  \item{pop}{a vector or factor giving the population assignment of each
    row of \code{x}, or a single numeric value specifying which column
    of \code{x} to use as population indicator. By default, the column
    labelled \code{"population"} is used.}
  \item{method}{a character string which should be one of
    ``extrapolation'', ``rarefaction'', ``raw'' or an unambiguous
    abbreviation of these.}
  \item{min.n}{the value of \eqn{n} used in the rarefaction method; by
    default, the smallest observed number of genotypes within a population.}
  \item{maxn}{the largest sample size used to calculate the rarefaction
    curve.}
  \item{type, xlab, ylab}{arguments passed to \code{plot}.}
  \item{plot}{a logical value specifying whether to do the rarefaction
    plot (\code{TRUE} by default).}
  \item{\dots}{arguments passed to and from methods.}
}
\details{
  \code{allelicrichness} computes for each locus in \code{x} the
  estimated allelic richness. Three methods are available: the
  extrapolation method (Foulley and Ollivier 2006), the rarefaction
  method (Hurlbert 1971), and the raw numbers of alleles.

  \code{rarefactionplot} computes the rarefaction curves of the number
  of alleles with respect to sample size using Hurlbert's (1971)
  method. A plot is made by default.
}
\value{
  \code{allelicrichness} returns a numeric matrix.

  \code{rarefactionplot} returns invisibly a list of matrices with the
  coordinates of the rarefaction plots for each locus.

  \code{rhost} returns a numeric vector.
}
\references{
  El Mousadik, A. and Petit, R. J. (1996) High level of genetic
  differentiation for allelic richness among populations of the argan
  tree [Argania spinosa (L. Skeels)] endemic to
  Morocco. \emph{Theoretical and Applied Genetics}, \bold{92},
  832--836.

  Foulley, J. L. and Ollivier, L. (2006) Estimating allelic richness and
  its diversity. \emph{Livestock Science}, \bold{101}, 150--158.

  Hurlbert, S. H. (1971) The nonconcept of species diversity: a critique
  and alternative parameters. \emph{Ecology}, \bold{52}, 577--586.
}
\author{Emmanuel Paradis}
\examples{
data(jaguar)
rarefactionplot(jaguar)
allelicrichness(jaguar)
rhost(jaguar)
}
\keyword{hplot}
