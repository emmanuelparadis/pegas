\name{site.spectrum}
\alias{site.spectrum}
\alias{plot.spectrum}
\title{Site Frequency Spectrum}
\description{
  \code{site.spectrum} computes the (un)folded site frequency spectrum
  of a set of aligned DNA sequences.
}
\usage{
site.spectrum(x, folded = TRUE, outgroup = 1)
\method{plot}{spectrum}(x, col = "red", main = NULL, ...)
}
\arguments{
  \item{x}{a set of DNA sequences (as an object of class
    \code{"DNAbin"}), or an object of class \code{"spectrum"}.}
  \item{folded}{a logical specifying whether to compute the folded site
    frequency spectrum (the default), or the unfolded spectrum if
    \code{folded = FALSE}.}
  \item{outgroup}{a single integer value giving which sequence is
    ancestral; ignored if \code{folded = TRUE}.}
  \item{col}{the colour of the barplot (red by default).}
  \item{main}{a character string for the title of the plot; a generic
    title is given by default (use \code{main = ""} to have no title).}
  \item{\dots}{further arguments passed to \code{\link[graphics]{barplot}}.}
}
\details{
  Under the infinite sites model of mutation, mutations occur on
  distinct sites, so every segregating (polymorphic) site defines a
  partition of the \eqn{n} sequences (see Wakeley, 2009). The \emph{site
  frequency spectrum} is a series of values where the \eqn{i}th element
  is the number of segregating sites defining a partition of \eqn{i} and
  \eqn{n - i} sequences. The \emph{unfolded} version requires to define
  an ancestral state with an external (outgroup) sequence, so \eqn{i}
  varies between 1 and \eqn{n - 1}. If no ancestral state can be
  defined, the \emph{folded} version is computed, so \eqn{i} varies
  between 1 and \eqn{n/2} or \eqn{(n - 1)/2}, for \eqn{n} even or odd,
  respectively.

  If \code{folded = TRUE}, sites with more than two states are ignored
  and a warning is returned giving how many were found.

  If \code{folded = FALSE}, sites with an ambiguous state at the
  external sequence are ignored and a warning is returned giving how
  many were found. Note that it is not checked if some sites have more
  than two states.
}
\value{
  \code{site.spectrum} returns an object of class \code{"spectrum"}
  which is a vector of integers (some values may be equal to zero) with
  the attribute \code{"folded"} (a logical value) indicating which
  version of the spectrum has been computed.
}
\author{Emmanuel Paradis}
\references{
  Wakeley, J. (2009) \emph{Coalescent Theory: An Introduction.}
  Greenwood Village, CO: Roberts and Company Publishers.
}
\seealso{
  \code{\link[ape]{DNAbin}} for manipulation of DNA sequences in R,
  \code{\link{haplotype}}
}
\examples{
require(ape)
data(woodmouse)
(sp <- site.spectrum(woodmouse))
plot(sp)
}
\keyword{manip}
\keyword{hplot}
