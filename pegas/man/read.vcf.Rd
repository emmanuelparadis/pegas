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
  The VCF file can compressed (*.gz) or not, but compressed files cannot
  be read remotely (see examples).

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
    \item{CHR}{the chromosome number (characters);}
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
\examples{
\dontrun{
## Chr Y from the 1000 genomes:
a <- "ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502"
b <- "ALL.chrY.phase3_integrated_v1a.20130502.genotypes.vcf.gz"
## WARNING: the name of the file above may have changed
url <- paste(a, b, sep = "/")
## file is compressed, so we download first:
download.file(url, "tmp.vcf.gz")
## no need to uncompress to read now that the file is local:
(x <- read.vcf("tmp.vcf.gz", 4000)) # read only 4000 loci

SNP <- is.snp(x)
table(SNP) # how many loci are really SNPs?

op <- par(mfcol = c(4, 1), xpd = TRUE)
lim <- c(2.65e6, 2.95e6)
## distribution of SNP and non-SNP mutations along the Y chr:
plot(attr(x, "POS"), !SNP, "h", col = "red", main = "non-SNP mutations")
rect(lim[1], -0.1, lim[2], 1.1, lwd = 2, lty = 2)
plot(attr(x, "POS"), SNP, "h", col = "blue", main = "SNP mutations")
rect(lim[1], -0.1, lim[2], 1.1, lwd = 2, lty = 2)
par(xpd = FALSE)
## same focusing on a smaller portion of the chromosome:
plot(attr(x, "POS"), !SNP, "h", col = "red", xlim = lim)
plot(attr(x, "POS"), SNP, "h", col = "blue", xlim = lim)
par(op)

## get haplotypes for the first 10 loci:
h <- haplotype(x, 1:10)
## plot their frequencies:
op <- par(mar = c(3, 10, 1, 1))
plot(h, horiz=TRUE, las = 1)
par(op)

## read the next 4000 loci:
(x2 <- read.vcf("tmp.vcf.gz", 4000, 4000))
}
}
\keyword{IO}
