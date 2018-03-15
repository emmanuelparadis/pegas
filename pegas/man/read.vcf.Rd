\name{read.vcf}
\alias{read.vcf}
\title{Read Variant Calling Format Files}
\description{
  This function reads allelic data from VCF (variant calling format) files.
}
\usage{
read.vcf(file, from = 1, to = 10000, which.loci = NULL, quiet = FALSE)
}
\arguments{
  \item{file}{a file name specified by either a variable of mode character,
    or a quoted string.}
  \item{from, to}{the loci to read; by default, the first 10,000.}
  \item{which.loci}{an alternative way to specify which loci to read is
    to give their indices (see \code{link{VCFloci}} how to obtain them).}
  \item{quiet}{a logical: should the progress of the operation be printed?}
}
\details{
  The VCF file can be compressed (*.gz) or not. Since pegas 0.11,
  compressed remote files can be read (see examples).

  A TABIX file is not required (and will be ignored if present).

  In the VCF standard, missing data are represented by a dot and these
  are read ``as is'' by the present function without trying to
  substitute by \code{NA}.
}
\note{
  Like for \code{\link{VCFloci}}, the present function can read either
  compressed (*.gz) or uncompressed files. There should be no difference
  in performance between both types of files if they are relatively
  small (less than 1 Gb as uncompressed, equivalent to ~50 Mb when
  compressed). For bigger files, it is more efficient to uncompress them
  (if disk space is sufficient), especially if they have to be accessed
  several times during the same session.
}
\value{
  an object of class \code{c("loci", "data.frame")}.
}
\author{Emmanuel Paradis}
\references{
  \url{http://www.1000genomes.org/node/101}

  \url{https://github.com/samtools/hts-specs}
}
\seealso{
  \code{\link{VCFloci}}, \code{\link{read.loci}},
  \code{\link{read.gtx}}, \code{\link{write.loci}}
}
\examples{
\dontrun{
## Chr Y from the 1000 Genomes:
a <- "ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502"
b <- "ALL.chrY.phase3_integrated_v1b.20130502.genotypes.vcf.gz"
## WARNING: the name of the file above may change
url <- paste(a, b, sep = "/")
## Solution 1: download first
download.file(url, "chrY.vcf.gz")
## no need to uncompress:
(info <- VCFloci("chrY.vcf.gz"))
str(info) # show the modes of the columns
## Solution 2: read remotely (since pegas 0.11)
info2 <- VCFloci(url)
identical(info, info2)
rm(info2)

SNP <- is.snp(info)
table(SNP) # how many loci are SNPs?
## compare with:
table(getINFO(info, "VT"))

op <- par(mfcol = c(4, 1), xpd = TRUE)
lim <- c(2.65e6, 2.95e6)
## distribution of SNP and non-SNP mutations along the Y chr:
plot(info$POS, !SNP, "h", col = "red", main = "non-SNP mutations",
     xlab = "Position", ylab = "", yaxt = "n")
rect(lim[1], -0.1, lim[2], 1.1, lwd = 2, lty = 2)
plot(info$POS, SNP, "h", col = "blue", main = "SNP mutations",
     xlab = "Position", ylab = "", yaxt = "n")
rect(lim[1], -0.1, lim[2], 1.1, lwd = 2, lty = 2)
par(xpd = FALSE)
## same focusing on a smaller portion of the chromosome:
plot(info$POS, !SNP, "h", col = "red", xlim = lim, xlab = "Position",
     ylab = "", yaxt = "n")
plot(info$POS, SNP, "h", col = "blue", xlim = lim, xlab = "Position",
     ylab = "", yaxt = "n")
par(op)

## read both types of mutations separately:
X.SNP <- read.vcf("chrY.vcf.gz", which.loci = which(SNP))
X.other <- read.vcf("chrY.vcf.gz", which.loci = which(!SNP))

identical(rownames(X.SNP), VCFlabels("chrY.vcf.gz")) # TRUE
cat(VCFheader("chrY.vcf.gz"))

## get haplotypes for the first 10 loci:
h <- haplotype(X.SNP, 1:10)
## plot their frequencies:
op <- par(mar = c(3, 10, 1, 1))
plot(h, horiz=TRUE, las = 1)
par(op)
}
}
\keyword{IO}
