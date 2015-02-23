\name{read.vcf}
\alias{read.vcf}
\title{Read Variant Calling Format Files}
\description{
  This function reads allelic data from VCF (variant calling format) files.
}
\usage{
read.vcf(file, nloci = 1000, skip = 0)
}
\arguments{
  \item{file}{a file name specified by either a variable of mode character,
    or a quoted string.}
  \item{nloci}{the number of loci to read.}
  \item{skip}{the number of loci to be skipped before reading the
    genetic data.}
}
\details{
  The VCF file can compressed (*.gz) or not.

  A TABIX file is not required.

  Because VCF files can be very big, only 1000 loci are read by default,
  but 5000 loci can be read rather quickly.

  In the VCF standard, missing data are represented by a dot and these
  are read ``as is'' by the present function without trying to
  substitute by \code{NA}.
}
\value{
  an object of class \code{c("loci", "data.frame")} with four additional
  attributes:
  \itemize{
    \item{CHR}{the chromosome number (integer values);}
    \item{POS}{the position of the locus (numeric values);}
    \item{QUAL}{the quality of the inferred locus (integer values);}
    \item{FILTER}{the filter (characters).}
  }
}
\author{Emmanuel Paradis}
\references{\url{http://www.1000genomes.org/node/101}}
\seealso{
  \code{\link{read.loci}}, \code{\link{read.gtx}},
  \code{\link{write.loci}}
}
\keyword{IO}
