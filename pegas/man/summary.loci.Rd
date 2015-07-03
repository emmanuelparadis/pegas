\name{summary.loci}
\alias{summary.loci}
\alias{[.loci}
\alias{print.summary.loci}
\alias{plot.summary.loci}
\alias{print.loci}
\title{Print and Summaries of Loci Objects}
\description{
  These functions print and summarize table of alleles and loci (objects
  of class \code{"loci"}).
}
\usage{
\method{print}{loci}(x, details = FALSE, ...)
\method{summary}{loci}(object, ...)
\method{print}{summary.loci}(x, ...)
\method{[}{loci}(x, i, j, drop = TRUE)
\method{plot}{summary.loci}(x, loci, what = "both", layout = 1, col = c("blue", "red"), ...)
}
\arguments{
  \item{x, object}{an object of class \code{"loci"} or \code{"summary.loci"}.}
  \item{details}{a logical value: if \code{TRUE} the data are printed as
    a data frame; the default is \code{FALSE}.}
  \item{i, j}{indices of the rows and/or columns to select or to drop.
    They may be numeric, logical, or character (in the same way than for
    standard R objects).}
  \item{drop}{a logical specifying whether to returned an object of
    the smallest dimension possible, i.e., may return a vector or a
    factor if \code{drop = TRUE} (this is not the default).}
  \item{loci}{the loci (genes) to be plotted. By default, all loci are
    plotted.}
  \item{what}{the frequencies to be plotted. Three choices are possible:
    \code{"alleles"}, \code{"genotypes"}, and \code{"both"} (the default),
    or any unambiguous abbreviations.}
  \item{layout}{the number of graphs to be plotted simultaneously.}
  \item{col}{the colours used for the barplots.}
  \item{\dots}{further arguments to be passed to or from other methods.}
}
\details{
  Genotypes not observed in the data frame are not counted.

  When using the \code{[} method, if only one column is extracted or if
  the returned data frame has no `loci' column, then the class
  \code{"loci"} is dropped.

  An object of class \code{"loci"} can be edited in the R data editor
  with, e.g., \code{fix(x)} or \code{x <- edit(x)}.

  \code{summary.loci} computes the absolute frequencies (counts); see
  the examples on how to compute the relative frequencies (proportions).
}
\value{
  \code{summary.loci} returns a list with the genes as names and each
  element made a list with two vectors \code{"genotype"} and
  \code{"allele"} with the frequencies (numbers) of genotypes and
  alleles, respectively. The names of these two vectors are the observed
  genotypes and alleles.

  \code{print} and \code{plot} methods return NULL.
}
\author{Emmanuel Paradis}
\seealso{
  \code{\link{read.loci}}, \code{\link{getAlleles}}, \code{\link{edit.loci}}
}
\examples{
data(jaguar)
s <- summary(jaguar)
\dontrun{
## works if the device is large enough:
plot(s, layout = 30, las = 2)
layout(1)
}
## compute the relative frequencies:
rapply(s, function(x) x/sum(x), how = "replace")
}
\keyword{manip}
\keyword{hplot}
