\name{read.gtx}
\alias{read.gtx}
\title{Read Genetix Data Files}
\description{
  This function reads allelic data from a Genetix file (.gtx).
}
\usage{
read.gtx(file)
}
\arguments{
  \item{file}{a file name specified by either a variable of mode character
    or a quoted string.}
}
\value{
  A data frame with class \code{c("loci", "data.frame")}.
}
\author{Emmanuel Paradis}
\note{
  The package \pkg{adegenet} has a similar function,
  \code{\link[adegenet]{read.genetix}}, but it returns an object of
  class \code{"genind"}.
}
\references{
  Belkhir, K., Borsa, P., Chikhi, L., Raufaste, N. and Bonhomme,
  F. (1996--2004) GENETIX 4.05, logiciel sous Windows(TM) pour la
  genetique des populations. Laboratoire Genome, Populations,
  Interactions, CNRS UMR 5000, Universite de Montpellier II, Montpellier
  (France). \url{https://kimura.univ-montp2.fr/genetix/}
}
\seealso{
  \code{\link{read.loci}}, \code{\link{write.loci}},
  \code{\link{read.vcf}}, \code{\link[adegenet]{read.genetix}}
}
\examples{
require(adegenet)
(X <- read.gtx(system.file("files/nancycats.gtx", package = "adegenet")))
## compare with the example in ?read.genetix
}
\keyword{IO}
