\name{haploFreq}
\alias{haploFreq}
\title{Haplotype Frequencies With a Covariate}
\description{
  This utility function extracts the absolute frequencies of haplotypes
  with respect to a categorical variable (a factor). The output is
  useful when ploting haplotype networks.
}
\usage{
haploFreq(x, fac, split = "_", what = 2, haplo = NULL)
}
\arguments{
  \item{x}{a set of DNA sequences (as an object of class
    \code{"DNAbin"}).}
  \item{fac}{a factor giving the categorical variable (can be missing).}
  \item{split}{a single character (see details).}
  \item{what}{a single integer (see details).}
  \item{haplo}{an object of class \code{"haplotype"}.}
}
\details{
  The frequencies of each haplotype in \code{x} are counted with respect
  to a factor which is either specified with \code{fac}, or extracted
  from the labels of \code{x}. In the second case, these labels are
  split with respect to the character specified in \code{split} and the
  \code{what}'th substrings are extracted and taken as the categorical
  variable (see example).

  If \code{haplo} is specified, the haplotype frequencies are taken from
  it, otherwise they are calculated from \code{x}.
}
\value{
  a matrix of counts.
}
\author{Klaus Schliep and Emmanuel Paradis}
\seealso{
  \code{\link{haplotype}}, \code{\link{haploNet}}
}
\examples{
## generate some artificial data from 'woodmouse':
data(woodmouse)
x <- woodmouse[sample(15, size = 50, replace = TRUE), ]
## labels IdXXX_PopXXX_LocXXX
rownames(x) <- paste("Id", 1:50, "_Pop", 1:2, "_Loc", 1:5, sep = "")
head(labels(x))
h <- haplotype(x)
## frequencies of haplotypes wrt 'Pop':
f.pop <- haploFreq(x, haplo = h)
## frequencies of haplotypes wrt 'Loc':
f.loc <- haploFreq(x, what = 3, haplo = h)
nt <- haploNet(h)
fq <- attr(nt, "freq")
op <- par(mfcol = c(1, 2))
plot(nt, size = fq, pie = f.pop, labels = FALSE)
plot(nt, size = fq, pie = f.loc, labels = FALSE)
par(op)
}
\keyword{manip}
