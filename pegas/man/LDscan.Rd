\name{LDscan}
\alias{LDscan}
\alias{LDscan.DNAbin}
\alias{LDscan.loci}
\alias{LDmap}
\title{Multi-Locus Linkage Disequilibrium}
\description{
  \code{LDscan} computes a matrix of pairwise linkage disequilibrium (LD)
  coefficients (\eqn{r^2}) from a set of loci (which must be bi-allelic;
  if not, the results are not guaranteed to be meaningful). The
  genotypes must be phased.

  \code{LDmap} plots a matrix of LD coefficients, optionally with the
  positions of the loci.
}
\usage{
LDscan(x, ...)

\method{LDscan}{DNAbin}(x, quiet = FALSE, ...)
\method{LDscan}{loci}(x, depth = NULL, quiet = FALSE, ...)

LDmap(d, POS = NULL, breaks = NULL, col = NULL, border = NA,
      angle = 0, asp = 1, cex = 1, scale.legend = 0.8, ...)
}
\arguments{
  \item{x}{an object of class \code{"loci"} with phased genotypes.}
  \item{depth}{a vector of integers giving the the depth(s) (or lags) at
    which the \eqn{r^2}'s are calculated. By default, all possible
    depths are considered.}
  \item{quiet}{a logical: should the progress of the operation be
    printed?}
  \item{d}{a correlation matrix (can be an object of class \code{"dist"}).}
  \item{POS}{an optional vector of locus positions (e.g., from a VCF
    file; see examples).}
  \item{breaks}{a vector of break intervals to count the values in
    \code{d}; by default, ten equally-sized intervals are used.}
  \item{col}{an optional vector of colours; a scale from lightyellow to
    red is used by default.}
  \item{border}{the border of the rectangles: the default is to have no
    border (this is not the same than default in
    \code{\link[graphics]{rect}}; see examples).}
  \item{angle}{value (in degrees) to rotate the graphic.}
  \item{asp}{the aspect ratio of the graphic; one by default so the
    elements are squares (not rectangles).}
  \item{cex}{the scaling of the labels and text.}
  \item{scale.legend}{the scaling of the legend rectangles.}
  \item{\dots}{further arguments passed to methods (\code{LDscan}) or to
    \code{plot.default} (\code{LDmap}).}
}
\details{
  The LD coefficient \eqn{r^2} is well defined when the two loci have
  only two alleles. In other cases, LD is well defined (see
  \code{\link{LD}}) but the definition of \eqn{r^2} is not clear.

  All levels of ploidy are accepted, but all loci should have the same
  ploidy level.

  If \code{depth} is used, the \eqn{r^2}'s are calculated only for the
  pairs of loci that are distant by these values in \code{x}, but
  necessarily on the chromosome. The returned list has names set with
  the values of \code{depth}.
}
\value{
  \code{LDscan} returns an object of class \code{"dist"} by default, or
  a list if \code{depth} is used.
}
\author{Emmanuel Paradis}
\seealso{\code{\link{LD}}, \code{\link{read.vcf}}}
\examples{
data(woodmouse)
d <- LDscan(woodmouse)
LDmap(d, seg.sites(woodmouse), seq(0, 1, .1))

\dontrun{
## Download the VCF file from Dryad:
## http://dx.doi.org/10.5061/dryad.446sv.2

## the VCF file should have this name:
fl <- "global.pop.GATK.SNP.hard.filters.V3.phased_all.pop.maf.05.recode.vcf.gz"

info.fly <- VCFloci(fl)

bks <- seq(0, 1, 0.2)

## LD map from the first 100 loci:
x <- read.vcf(fl, to = 100) # read only 100 loci
res <- LDscan(x)
LDmap(res, info.fly$POS[1:100], bks, scale.legend = 3)

## check the chromosomes:
table(info.fly$CHROM)

## LD map from 100 loci randomly distributed on the chromosome:
s <- ceiling(seq(1, 224253, length.out = 100))
xs <- read.vcf(fl, which.loci = s)
res2 <- LDscan(xs)
LDmap(res2, info.fly$POS[s], bks, scale.legend = 3)

## something simpler with 10 loci:
x10 <- x[, 1:10]
## the VCF file has no locus IDs, so we give some here:
names(x10) <- paste0("Loc", 1:10)
res10 <- LDscan(x10, quiet = TRUE)
LDmap(res10, angle = 45, border = NULL)
}
}
\keyword{models}
