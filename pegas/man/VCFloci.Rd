\name{VCFloci}
\alias{VCFloci}
\alias{VCFheader}
\alias{VCFlabels}
\alias{print.VCFinfo}
\alias{is.snp.VCFinfo}
\alias{rangePOS}
\alias{selectQUAL}
\alias{getINFO}
\title{Information From VCF Files}
\description{
  These functions help to extract information from VCF files and to
  select which loci to read with \code{\link{read.vcf}}.
}
\usage{
VCFloci(file, what = "all", chunck.size = 1e9, quiet = FALSE)
\method{print}{VCFinfo}(x, \dots)
VCFheader(file)
VCFlabels(file)
\method{is.snp}{VCFinfo}(x)
rangePOS(x, from, to)
selectQUAL(x, threshold = 20)
getINFO(x, what = "DP", as.is = FALSE)
}
\arguments{
  \item{file}{file name of the VCF file.}
  \item{what}{a character specifying the information to be extracted
    (see details).}
  \item{chunck.size}{the size of data in bytes read at once.}
  \item{quiet}{a logical: should the progress of the operation be printed?}
  \item{x}{an object of class \code{"VCFinfo"}.}
  \item{from, to}{integer values giving the range of position values.}
  \item{threshold}{a numerical value indicating the minimum value of
    quality for selecting loci.}
  \item{as.is}{a logical. By default, \code{getINFO} tries to convert
    its output as numeric: if too many NA's are produced, the output is
    returned as character. Use \code{as.is = TRUE} to force the output
    to be in character mode.}
  \item{\dots}{further arguments passed to and from other methods.}
}
\details{
  The variant call format (VCF) is described in details in the
  References. Roughly, a VCF file is made of two parts: the header and
  the genotypes. The last line of the header gives the labels of the
  genotypes: the first nine columns give information for each locus and
  are (always) "CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER",
  "INFO", and "FORMAT". The subsequent columns give the labels
  (identifiers) of the individuals; these may be missing if the file
  records only the variants. Note that the data are arranged as the
  transpose of the usual way: the individuals are as columns and the
  loci are as rows.

  \code{VCFloci} is the main function documented here: it reads the
  information relative to each locus. The option \code{what} specifies
  which column(s) to read. By default, all of them are read. If the user
  is interested in only the locus positions, the option \code{what =
  "POS"} would be used.

  Since VCF files can be very big, the data are read in portions of
  \code{chunk.size} bytes. The default (1 Gb) should be appropriate in
  most situations. This value should not exceed 2e9.

  \code{VCFheader} returns the header of the VCF file (excluding the
  line of labels). \code{VCFlabels} returns the individual labels.

  The output of \code{VCFloci} is a data frame with as many rows as
  there are loci in the VCF file and storing the requested
  information. The other functions help to extract specific information
  from this data frame: their outputs may then be used to select which
  loci to read with \code{\link{read.vcf}}.

  \code{is.snp} tests whether each locus is a SNP (i.e., the reference
  allele, REF, is a single charater and the alternative allele, ALT,
  also). It returns a logical vector with as many values as there are
  loci. Note that some VCF files have the information VT (variant type)
  in the INFO column.

  \code{rangePOS} and \code{selectQUAL} select some loci with respect to
  values of position or quality. They return the indices (i.e., row
  numbers) of the loci satisfying the conditions.

  \code{getINFO} extracts a specific information from the INFO
  column. By default, these are the total depths (DP) which can be
  changed with the option \code{what}. The meaning of these information
  should be described in the header of the VCF file.
}
\note{
  \code{VCFloci} is able to read either compressed (*.gz) or
  uncompressed files.
}
\value{
  \code{VCFloci} returns an object of class \code{"VCFinfo"} which is a
  data frame with a specific print method.

  \code{VCFheader} returns a single character string which can be
  printed nicely with \code{cat}.

  \code{VCFlabels} returns a vector of mode character.

  \code{is.snp} returns a vector of mode logical.

  \code{rangePOS} and \code{selectQUAL} return a vector of mode
  numeric.

  \code{getINFO} returns a vector of mode character or numeric (see above).
}
\references{
  \url{http://www.1000genomes.org/node/101}

  \url{https://github.com/samtools/hts-specs}
}
\author{Emmanuel Paradis}
\seealso{\code{\link{read.vcf}}}
\examples{
## see ?read.vcf
}
\keyword{IO}
\keyword{manip}
