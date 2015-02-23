\name{tajima.test}
\alias{tajima.test}
\title{Test of the Neutral Mutation Hypothesis}
\description{
  This function tests the neutral mutation hypothesis with Tajima's
  \emph{D}.
}
\usage{
tajima.test(x)
}
\arguments{
  \item{x}{a set of DNA sequences (object of class \code{"DNAbin"}).}
}
\details{
}
\value{
  A list with three numeric values:
  \item{D}{Tajima's \emph{D} statistic.}
  \item{Pval.normal}{the p-value assuming that \emph{D} follows a normal
    distribution with mean zero and variance one.}
  \item{Pval.beta}{the p-value assuming that \emph{D} follows a beta
    distribution after rescaling on [0, 1] (Tajima, 1989).}
}
\author{Emmanuel Paradis}
\references{
  Tajima, F. (1989) Statistical method for testing the neutral mutation
  hypothesis by DNA polymorphism. \emph{Genetics}, \bold{123}, 595--595.

}
\examples{
require(ape)
data(woodmouse)
tajima.test(woodmouse)
}
\keyword{htest}
