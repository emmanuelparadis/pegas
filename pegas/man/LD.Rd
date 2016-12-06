\name{LD}
\alias{LD}
\alias{LD2}
\title{Linkage Disequilibrium}
\description{
  These two functions analyse linkage disequilibrium in the case of
  phased (\code{LD}) or unphased (\code{LD2}) genotypes.
}
\usage{
LD(x, locus = c(1, 2), details = TRUE)
LD2(x, locus = c(1, 2), details = TRUE)
}
\arguments{
  \item{x}{an object of class \code{"loci"}.}
  \item{locus}{a vector of two integers giving the loci to analyse.}
  \item{details}{a logical value indicating whether to print the
    correlation matrix among alleles.}
}
\details{
  These functions consider a pair of loci and compute the correlations
  among pairs of alleles.

  \code{LD} first scans the data for unphased genotypes: all individuals
  with at least one unphased genotype are dropped with a warning. It is
  based on the observed frequencies of haplotypes (Zaykin et
  al. 2008). \code{LD2} is based on the observed frequencies of
  different genotypes (Schaid 2004).

  Both functions accept any number of alleles. \code{LD} can work with
  any level of ploidy; \code{LD2} works with diploid data.

  The present version does not test the significance of the \eqn{T_2}
  test (Zaykin et al. 2008) with permutations. These authors present
  simulation results suggesting that the chi-squared approximation has
  similar type I error rates and power than the test based on
  permutations even for small sample sizes. Furthermore, this test has
  better statistical properties than alternatives such as those reported
  here (LRT and Pearson's test).
}
\value{
  For both functions, if \code{details = FALSE}, only the T2 test is
  returned.

  For \code{LD}: if \code{details = TRUE}, a named list with the
  following elements:

  \item{Observed frequencies}{the counts of haplotypes in the data.}
  \item{Expected frequencies}{the expected frequencies of haplotypes
    computed from the observed proportions of alleles under the
    assumption of no linkage disequilibrium.}
  \item{Correlations among alleles}{the observed correlations among
    alleles from both loci.}
  \item{LRT (G-squared)}{the likelihood-ratio test of the null
    hypothesis of no linkage disequilibrium.}
  \item{Pearson's test (chi-squared)}{the chi-squared test based on
    haplotypes counts.}
  \item{T2}{the \eqn{T_2} test with its number of degrees of freedom
    (df).}

  For \code{LD2}: if \code{details = TRUE}, a named list with two
  elements:

  \item{Delta}{the correlations among alleles (denoted \eqn{Delta} in
    Schaid 2004).}
  \item{T2}{the \eqn{T_2} test with its number of degrees of freedom
    (df).}
}
\author{Emmanuel Paradis}
\references{
  Schaid, D. J. (2004) Linkage disequilibrium testing when linkage phase
  is unknown. \emph{Genetics}, \bold{166}, 505--512.

  Zaykin, D. V., Pudovkin, A. and  Weir, B. S. (2008) Correlation-based
  inference for linkage disequilibrium with multiple alleles.
  \emph{Genetics}, \bold{180}, 533--545.
}
\seealso{
  \code{\link{haplotype.loci}}, \code{\link{is.phased}},
  \code{\link{LDscan}}}
\examples{
data(jaguar)
LD2(jaguar, details = FALSE)
LD2(jaguar, locus = 8:9, details = FALSE)
}
\keyword{models}
