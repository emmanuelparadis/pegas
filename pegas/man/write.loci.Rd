\name{write.loci}
\alias{write.loci}
\title{Write Allelic Data Files}
\description{
  This function writes allelic data into a text file.
}
\usage{
write.loci(x, file = "", loci.sep = " ", allele.sep = "/|", ...)
}
\arguments{
  \item{x}{an object of class \code{"loci"}.}
  \item{file}{a file name specified by either a variable of mode character,
    or a quoted string. By default, the data are printed on the console.}
  \item{loci.sep}{the character(s) use to separate the loci (columns) in
    the file (a space by default).}
  \item{allele.sep}{the character(s) used to separate the alleles for
    each locus in the file (a slash by default).}
  \item{\dots}{further arguments passed to \code{write.table}.}
}
\value{
  NULL
}
\author{Emmanuel Paradis}
\seealso{
  \code{\link{read.loci}}, \code{\link[utils]{write.table}} for all its options
}
\examples{
data(jaguar)
x <- jaguar[1:10, 1:3] # take a small subset
write.loci(x)
## use of '...':
write.loci(x, loci.sep = "\t", quote = FALSE, col.names = FALSE)
}
\keyword{IO}
