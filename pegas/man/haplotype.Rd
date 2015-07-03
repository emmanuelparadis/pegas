\name{haplotype}
\alias{haplotype}
\alias{haplotype.DNAbin}
\alias{plot.haplotype}
\alias{print.haplotype}
\alias{plot.haplotype}
\alias{sort.haplotype}
\alias{[.haplotype}
\title{Haplotype Extraction and Frequencies}
\description{
  \code{haplotype} extracts the haplotypes from a set of DNA
  sequences. The result can be plotted with the appropriate function.
}
\usage{
haplotype(x, ...)
\method{haplotype}{DNAbin}(x, labels = NULL, ...)
\method{plot}{haplotype}(x, ...)
\method{print}{haplotype}(x, ...)
\method{sort}{haplotype}(x,
     decreasing = ifelse(what == "frequencies", TRUE, FALSE),
     what = "frequencies", ...)
\method{[}{haplotype}(x, ...)
}
\arguments{
  \item{x}{a set of DNA sequences (as an object of class
    \code{"DNAbin"}), or an object of class \code{"haplotype"}.}
  \item{labels}{a vector of character strings used as names for the rows
    of the returned object. By default, Roman numerals are given.}
  \item{\dots}{further arguments passed to
    \code{\link[graphics]{barplot}} (unused in \code{print} and
    \code{sort}).}
  \item{decreasing}{a logical value specifying in which order to sort
    the haplotypes; by default this depends on the value of
    \code{what}.}
  \item{what}{a character specifying on what feature the haplotypes
    should be sorted: this must be \code{"frequencies"} or
    \code{"labels"}, or an unambiguous abbreviation of these.}
}
\details{
The \code{sort} method sorts the haplotypes in decreasing frequencies
(the default) or in alphabetical order of their labels (if \code{what =
"labels"}). Note that if these labels are Roman numerals (as assigned by
\code{haplotype}), their alphabetical order may not be their numerical
one (e.g., IX is alphabetically before VIII).

From \pkg{pegas} 0.7, \code{haplotype} extracts haplotypes taking into
account base ambiguities.
}
\value{
  \code{haplotype} returns an object of class \code{c("haplotype",
    "DNAbin")} which is an object of class \code{"DNAbin"} with two
  additional attributes: \code{"index"} identifying the index of each
  observation that share the same haplotype, and \code{"from"} giving
  the name of the original data.

  \code{sort} returns an object of the same class respecting its
  attributes.
}
\author{Emmanuel Paradis}
\seealso{
  \code{\link{haploNet}},  \code{\link{haploFreq}},
  \code{\link{subset.haplotype}},
  \code{\link[ape]{DNAbin}} for manipulation of DNA sequences in R.

  The \code{haplotype} method for objects of class \code{"loci"} is
  documented separately: \code{\link{haplotype.loci}}.
}
\examples{
## generate some artificial data from 'woodmouse':
data(woodmouse)
x <- woodmouse[sample(15, size = 110, replace = TRUE), ]
(h <- haplotype(x))
## the indices of the individuals belonging to the 1st haplotype:
attr(h, "index")[[1]]
plot(sort(h))
}
\keyword{manip}
\keyword{hplot}
