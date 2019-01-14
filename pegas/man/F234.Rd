\name{F4}
\alias{F2}
\alias{F3}
\alias{F4}
\title{F-Statistics From Patterson et al}
\description{
  These functions compute the \emph{F}-statistics developed by Patterson
  et al. (2012).
}
\usage{
F2(x, allele.freq = NULL, population = NULL, check.data = TRUE,
   pops = NULL, jackknife.block.size = 10, B = 1e4)
F3(x, allele.freq = NULL, population = NULL, check.data = TRUE,
   pops = NULL, jackknife.block.size = 10, B = 1e4)
F4(x, allele.freq = NULL, population = NULL, check.data = TRUE,
   pops = NULL, jackknife.block.size = 10, B = 1e4)
}
\arguments{
  \item{x}{an object of class \code{"loci"}.}
  \item{allele.freq}{alternatively, a list of allele (absolute)
    frequencies as output by \code{\link{by.loci}} (if this is used,
    \code{x} is ignored).}
  \item{population}{a column name or number giving which column of
    \code{x} should be treated as the population variable. By default,
    the column named ``population'' is used.}
  \item{check.data}{if \code{FALSE}, it is assumed that the user checked
    that all loci are strict SNPs. By default, the data are checked for
    the number of alleles and the non-SNP loci are dropped with a
    warning.}
  \item{pops}{a vector giving two, three, or four population names
    depending on the function. The order of these names is important
    (see Patterson et al. 2012). By default, the populations in the
    \code{x} are taken in the order they appear, and an error is
    returned if the number of populations does not match number required
    by the function.}
  \item{jackknife.block.size}{the size of the block used in the
    jackknife to assess the significance of the \emph{F}-statistic (this
    should be around one thousandth of the number of loci, or not less
    than 10.}
  \item{B}{the number of replications of the bootstrap used to assess
    the significance of the \emph{F}-statistic.}
}
\details{
  These functions are provisional versions.

  It is much better to compute the allele frequencies, and then use
  \code{allele.freq} with different combinations of \code{pops}.
}
\value{
  A vector with names.
}
\references{
  Patterson, N., Moorjani, P., Luo, Y., Mallick, S., Rohland, N., Zhan,
  Y., Genschoreck, T., Webster, T. and Reich, D. (2012) Ancient
  admixture in human history. \emph{Genetics}, \bold{192}, 1065--1093.
}
\author{Emmanuel Paradis}
\seealso{
  \code{\link{by.loci}}, \code{\link{Fst}}, the package
  \pkg{admixturegraph} that can draw graphs from the output of this
  function.
}
\keyword{htest}
