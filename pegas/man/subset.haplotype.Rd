\name{subset.haplotype}
\alias{subset.haplotype}
\title{Subsetting and Filtering Haplotypes}
\description{
  This function selects haplotypes based on their (absolute) frequencies
  and/or proportions of missing nucleotides.
}
\usage{
\method{subset}{haplotype}(x, minfreq = 1, maxfreq = Inf, maxna = Inf, na = c("N", "?"), ...)
}

\arguments{
  \item{x}{an object of class \code{c("haplotype", "DNAbin")}.}
  \item{minfreq, maxfreq}{the lower and upper limits of (absolute)
    haplotype frequencies. By default, all haplotypes are selected
    whatever their frequency.}
  \item{maxna}{the maximum frequency (absolute or relative; see details)
    of missing nucleotides within a given haplotype.}
  \item{na}{a vector of mode character specifying which nucleotide
    symbols should be treated as missing data; by default, unknown
    nucleotide (N) and completely unknown site (?) (can be lower- or
    uppercase). There are two shortcuts: see details.}
  \item{\dots}{unused.}
}
\details{
  The value of \code{maxna} can be either less than one, or greater or
  equal to one. In the former case, it is taken as specifying the
  maximum proportion (relative frequency) of missing data within a given
  haplotype. In the latter case, it  is taken as the maximum number
  (absolute frequency).

  \code{na = "all"} is a shortcut for all ambiguous nucleotides
  (including N) plus alignment gaps and completely unknown site (?).

  \code{na = "ambiguous"} is a shortcut for only ambiguous nucleotides
  (including N).
}
\value{
  an object of class \code{c("haplotype", "DNAbin")}.
}
\author{Emmanuel Paradis}
\seealso{\code{\link{haplotype}}}
\examples{
data(woodmouse)
h <- haplotype(woodmouse)
subset(h, maxna = 20)
subset(h, maxna = 20/ncol(h)) # same thing than above
}
\keyword{manip}
