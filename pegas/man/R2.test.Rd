\name{R2.test}
\alias{R2.test}
\title{Ramos-Onsins--Rozas Test of Neutrality}
\description{
  This function computes Ramos-Onsins and Rozas's test of neutrality
  for a set of DNA sequences.
}
\usage{
R2.test(x, B = 1000, theta = 1, plot = TRUE, quiet = FALSE, ...)
}
\arguments{
  \item{x}{a DNA matrix (object of class \code{"DNAbin"}).}
  \item{B}{the number of replicates used for the simulation procedure.}
  \item{theta}{the value of the \eqn{\theta}{theta} population parameter
    used in the simulation.}
  \item{plot}{a logical value specifying whether to plot the results
    (\code{TRUE} by default).}
  \item{quiet}{a logical value specifying whether to not display the
    progress of the simulations. The default is \code{FALSE} meaning
    that a progress bar is displayed by default.}
  \item{\dots}{further arguments passed to \code{hist}.}
}
\note{
  The simulation procedure probably needs to be tested and
  improved. However the results make sense so far.
}
\value{
  a list with two elements: \code{R2} the value of the test statistic
  \eqn{R_2}{R_2}, and \code{P.val} the associated \emph{P}-value. If
  \code{B = 0} a single value, the test statistic, is returned
}
\author{Emmanuel Paradis}
\references{
  Ramos-Onsins, R. and Rozas, R. (2002) Statistical properties of new
  neutrality tests against population growth. \emph{Molecular Biology
    and Evolution}, \bold{19}, 2092--2100.

  Sano, J. and Tachida, G. (2005) Gene genealogy and properties of test
  statistics of neutrality under population growth. \emph{Genetics},
  \bold{169}, 1687--1697.
}
\seealso{
  \code{\link[ape]{read.dna}}, \code{\link[ape]{dist.dna}}
}
\examples{
data(woodmouse)
R2.test(woodmouse, quiet = TRUE)
}
\keyword{htest}
