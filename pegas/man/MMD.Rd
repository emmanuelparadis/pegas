\name{MMD}
\alias{MMD}
\title{Mismatch Distribution}
\description{
  This function draws a histogram of the frequencies of pairwise
  distances from a set of DNA sequences.
}
\usage{
MMD(x, xlab = "Distance", main = "", rug = TRUE, legend = TRUE,
    lcol = c("blue", "red"), lty = c(1, 1), bw = 2, ...)
}
\arguments{
  \item{x}{a set of DNA sequences (object of class \code{"DNAbin"}).}
  \item{xlab}{the label for the x-axis.}
  \item{main}{the title (none by default).}
  \item{rug}{a logical specifying whether to add a rug of the pairwise
    distances on the horizontal axis (see \code{\link[graphics]{rug}}).}
  \item{legend}{a logical specifying whether to draw a legend.}
  \item{lcol}{the colours used for the curves.}
  \item{lty}{the line types for the curves}
  \item{bw}{the bandwidth used for the empirical density curve (passed
    to \code{\link[stats]{density}}).}
  \item{\dots}{further arguments passed to \code{hist}.}
}
\details{
  The histogram shows the observed distribution of pairwise
  distances. The lines show an empirical density estimate (in blue) and
  the expected distribution under stable population (Rogers and
  Harpending 1992).
}
\author{Emmanuel Paradis and David Winter}
\references{
  Rogers, A. R. and Harpending, H. (1992) Population growth makes waves
  in the distribution of pairwise genetic-differences. \emph{Molecular
    Biology and Evolution}, \bold{9}, 552--569.
}
\examples{
data(woodmouse)
MMD(woodmouse, col = "grey")
MMD(woodmouse, breaks = 20, legend = FALSE)
MMD(woodmouse, lty = 1:2, lcol = rep("black", 2), col = "lightgrey")
}
\keyword{hplot}
