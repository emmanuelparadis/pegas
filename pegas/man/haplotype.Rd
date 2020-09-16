\name{haplotype}
\alias{haplotype}
\alias{haplotype.DNAbin}
\alias{haplotype.character}
\alias{haplotype.numeric}
\alias{plot.haplotype}
\alias{print.haplotype}
\alias{summary.haplotype}
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
\method{haplotype}{DNAbin}(x, labels = NULL, strict = FALSE,
                  trailingGapsAsN = TRUE, ...)
\method{haplotype}{character}(x, labels = NULL, ...)
\method{haplotype}{numeric}(x, labels = NULL, ...)
\method{plot}{haplotype}(x, xlab = "Haplotype", ylab = "Number", ...)
\method{print}{haplotype}(x, ...)
\method{summary}{haplotype}(object, ...)
\method{sort}{haplotype}(x,
     decreasing = ifelse(what == "frequencies", TRUE, FALSE),
     what = "frequencies", ...)
\method{[}{haplotype}(x, ...)
}
\arguments{
  \item{x}{a set of DNA sequences (as an object of class
    \code{"DNAbin"}), or an object of class \code{"haplotype"}.}
  \item{object}{an object of class \code{"haplotype"}.}
  \item{labels}{a vector of character strings used as names for the rows
    of the returned object. By default, Roman numerals are given.}
  \item{strict}{a logical value; if \code{TRUE}, ambiguities and gaps in
    the sequences are ignored.}
  \item{trailingGapsAsN}{a logical value; if \code{TRUE} (the default),
    the leading and trailing alignment gaps are considered as unknown
    bases (i.e., N). This option has no effect if \code{strict = TRUE}.}
  \item{xlab, ylab}{labels for the x- and x-axes.}
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
The way ambiguities in the sequences are taken into account is explained
in a post to r-sig-phylo (see the examples below):

\url{https://www.mail-archive.com/r-sig-phylo@r-project.org/msg05541.html}

The \code{sort} method sorts the haplotypes in decreasing frequencies
(the default) or in alphabetical order of their labels (if \code{what =
"labels"}). Note that if these labels are Roman numerals (as assigned by
\code{haplotype}), their alphabetical order may not be their numerical
one (e.g., IX is alphabetically before VIII).

From \pkg{pegas} 0.7, \code{haplotype} extracts haplotypes taking into
account base ambiguities (see Note below).
}
\note{
The presence of ambiguous bases and/or alignment gaps in DNA sequences
can make the interpretation of haplotypes difficult. It is recommended
to check their distributions with \code{\link[ape]{image.DNAbin}} and
\code{\link[ape]{base.freq}} (using the options in both functions).

Comparing the results obtained playing with the options \code{strict}
and \code{trailingGapsAsN} of \code{haplotype.DNAbin} may be useful.
Note that the \pkg{ape} function \code{\link[ape]{seg.sites}} has the
same two options (as from \pkg{ape} 5.4) which may be useful to find the
relevant sites in the sequence alignment.
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
\note{
  There are cases where the algorithm that pools the different sequences
  into haplotypes has difficulties, although it seems to require a
  specific configuration of missing/ambiguous data. The last example
  below is one of them.
}
\examples{
## generate some artificial data from 'woodmouse':
data(woodmouse)
x <- woodmouse[sample(15, size = 110, replace = TRUE), ]
(h <- haplotype(x))
## the indices of the individuals belonging to the 1st haplotype:
attr(h, "index")[[1]]
plot(sort(h))
## get the frequencies in a named vector:
setNames(lengths(attr(h, "index")), labels(h))

## data posted by Hirra Farooq on r-sig-phylo (see link above):
cat(">[A]\nCCCGATTTTATATCAACATTTATTT------",
    ">[D]\nCCCGATTTT----------------------",
    ">[B]\nCCCGATTTTATATCAACATTTATTT------",
    ">[C]\nCCCGATTTTATATCACCATTTATTTTGATTT",
    file = "x.fas", sep = "\n")
x <- read.dna("x.fas", "f")
unlink("x.fas")

## show the sequences and the distances:
alview(x)
dist.dna(x, "N", p = TRUE)

## by default there are 3 haplotypes with a warning about ambiguity:
haplotype(x)

## the same 3 haplotypes without warning:
haplotype(x, strict = TRUE)

## if we remove the last sequence there is, by default, a single haplotype:
haplotype(x[-4, ])

## to get two haplotypes separately as with the complete data:
haplotype(x[-4, ], strict = TRUE)

## a simpler example:
y <- as.DNAbin(matrix(c("A", "A", "A", "A", "R", "-"), 3))
haplotype(y) # 1 haplotype
haplotype(y, strict = TRUE) # 3 haplotypes
haplotype(y, trailingGapsAsN = FALSE) # 2 haplotypes

## a tricky example with 4 sequences and 1 site:
z <- as.DNAbin(matrix(c("Y", "A", "R", "N"), 4))
alview(z, showpos = FALSE)

## a single haplotype is identified:
haplotype(z)
## 'Y' has zero-distance with (and only with) 'N', so they are pooled
## together; at a later iteration of this pooling step, 'N' has
## zero-distance with 'R' (and ultimately with 'A') so they are pooled

## if the sequences are ordered differently, 'Y' and 'A' are separated:
haplotype(z[c(4, 1:3), ])
}
\keyword{manip}
\keyword{hplot}
