\name{mjn}
\alias{mjn}
\title{Median-Joining Network}
\description{
  This function computes the median-joining network (MJN) as described
  by Bandelt et al. (1999).
}
\usage{
mjn(x, epsilon = 0, max.n.cost = 10000, prefix = "median.vector_")
}
\arguments{
  \item{x}{a matrix (or data frame) of DNA sequences or binary 0/1 data.}
  \item{epsilon}{tolerance parameter.}
  \item{max.n.cost}{the maximum number of costs to be computed.}
  \item{prefix}{the prefix used to label the median vectors.}
}
\details{
  MJN is a network method where intermediate (unobserved) sequences (the
  median vectors) are reconstructed and included in the final
  network. Unlike \code{\link{mst}}, \code{rmst}, and \code{msn},
  \code{mjn} works with the original sequences, the distances being
  calculated internally using a Hamming distance method (with
  \code{dist(x, "manhattan")} for binary data or \code{dist.dna(x,
    "N")} for DNA sequences).

  The parameter \code{epsilon} controls how the search for new median
  vectors is performed: the larger this parameter, the wider the search
  (see the example with binary data).

  If the sequences are very divergent, the search for new median vectors
  can take a very long time. The argument \code{max.n.cost} controls how
  many such vectors are added to the network (the default value should
  avoid the function to run endlessly).
}
\note{
  This version still needs to be tested with large data sets. Bandelt et
  al. (1999) reported long computing times because of the need to
  compute a lot of median vectors.
}
\value{
  an object of class \code{"\link{haploNet}"} with an extra attribute
  (data) containing the original data together with the median vectors.
}
\references{
  Bandelt, H. J., Forster, P. and Rohl, A. (1999) Median-joining networks
  for inferring intraspecific phylogenies. \emph{Molecular Biology and
    Evolution}, \bold{16}, 37--48.
}
\author{Emmanuel Paradis}
\seealso{
  \code{\link{haploNet}}, \code{\link{mst}}
}
\examples{
## data in Table 1 of Bandelt et al. (1999):
x <- c(0, 0, 0, 0, 0, 0, 0, 0, 0,
       1, 1, 1, 1, 0, 0, 0, 0, 0,
       1, 0, 0, 0, 1, 1, 1, 0, 0,
       0, 1, 1, 1, 1, 1, 0, 1, 1)
x <- matrix(x, 4, 9, byrow = TRUE)
rownames(x) <- LETTERS[1:4]
(nt0 <- mjn(x))
(nt1 <- mjn(x, 1))
(nt2 <- mjn(x, 2))
plot(nt0)

\dontrun{
## same like in Fig. 4 of Bandelt et al. (1999):
plotNetMDS(nt2, dist(attr(nt2, "data"), "manhattan"), 3)
}

## data in Table 2 of Bandelt et al. (1999):
z <- list(c("g", "a", "a", "a", "a", "a", "a", "a", "a", "a", "a", "a"),
          c("a", "g", "g", "a", "a", "a", "a", "a", "a", "a", "a", "a"),
          c("a", "a", "a", "g", "a", "a", "a", "a", "a", "a", "g", "g"),
          c("a", "a", "a", "a", "g", "g", "a", "a", "a", "a", "g", "g"),
          c("a", "a", "a", "a", "a", "a", "a", "a", "g", "g", "c", "c"),
          c("a", "a", "a", "a", "a", "a", "g", "g", "g", "g", "a", "a"))
names(z) <- c("A1", "A2", "B1", "B2", "C", "D")
z <- as.matrix(as.DNAbin(z))
(ntz <- mjn(z, 2))

\dontrun{
## same like in Fig. 5 of Bandelt et al. (1999):
plotNetMDS(ntz, dist.dna(attr(ntz, "data"), "N"), 3)
}
}
\keyword{models}
