\name{rr.test}
\alias{rr.test}
\title{Tajima Relative Rate Test of Molecular Clock}
\description{
  This function tests the hypothesis of a molecular evolutionary clock
  (i.e., a constant rate of molecular evolution) between two samples
  using an outgroup sample. It can be applied to both nucleotide and
  amino acid sequences.
}
\usage{
rr.test(x, y, out)
}
\arguments{
  \item{x, y}{a single DNA sequence (object class \code{"DNAbin"}).}
  \item{out}{a single DNA sequence to be used as outgroup.}
}
\value{
  a list with two numeric values: \code{Chi} (Chi-squared statistic) and
  \code{Pval} (the P-value).
}
\author{Alastair Potts \email{potts.a@gmail.com}}
\references{
  Tajima, F. (1993) Simple methods for testing molecular clock
  hypothesis. \emph{Genetics}, \bold{135}, 599--607. (Equation 4)
}
\examples{
require(ape)
data(woodmouse)
rr.test(x = woodmouse[2, ], y = woodmouse[3, ], out = woodmouse[1, ])

# Test all pairs in a sample:
outgroup <- woodmouse[1, ]
n <- nrow(woodmouse)
cc <- combn(2:n, 2)
FUN <- function(x)
    rr.test(woodmouse[x[1], ], woodmouse[x[2], ], outgroup)$Pval
OUT <- apply(cc, 2, FUN)
### two ways to arrange the output:
RES <- matrix(NA, n - 1, n - 1)
RES[row(RES) > col(RES)] <- OUT
RES <- t(RES)
RES[row(RES) > col(RES)] <- OUT
RES <- t(RES)
dimnames(RES) <- list(2:n, 2:n)
RES <- as.dist(RES)
### 2nd method:
class(OUT) <- "dist"
attr(OUT, "Labels") <- as.character(2:15)
attr(OUT, "Size") <- n - 1L
attr(OUT, "Diag") <- attr(OUT, "Upper") <- FALSE
### they are the same:
all(OUT == RES)
}
\keyword{htest}
