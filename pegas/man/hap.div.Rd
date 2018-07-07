\name{hap.div}
\alias{hap.div}
\alias{hap.div.haplotype}
\alias{hap.div.DNAbin}
\title{Haplotype Diversity}
\description{
  This function computes haplotype diversity from DNA
  sequences. This is a generic function.
}
\usage{
hap.div(x, ...)
\method{hap.div}{haplotype}(x, variance = FALSE, method = "Nei", ...)
\method{hap.div}{DNAbin}(x, variance = FALSE, method = "Nei", ...)
}
\arguments{
  \item{x}{an object with DNA data.}
  \item{variance}{a logical value specifying whether to calculate the
    variance of the estimated haplotype diversity.}
  \item{method}{(unused, see details).}
  \item{\dots}{further arguments passed to and from methods.}
}
\details{
  Currently, only Nei and Tajima's (1981) method is available.
}
\value{
  a numeric vector with one or two values (if \code{variance = TRUE}).
}
\author{Emmanuel Paradis}
\references{
  Nei, M. and Tajima, F. (1981) DNA polymorphism detectable by
  restriction endonuclease. \emph{Genetics}, \bold{97}, 145--163.
}
\seealso{
  \code{\link{nuc.div}}
}
\examples{
data(woodmouse)
hap.div(woodmouse) # all haplotypes are unique

## neuraminidase sequences from the 2009 H1N1 data (delivered with adegenet):
fl <- system.file("files/pdH1N1-NA.fasta", package = "adegenet")
H1N1.NA <- read.dna(fl, "fasta")
hap.div(H1N1.NA, TRUE)
}
\keyword{manip}
