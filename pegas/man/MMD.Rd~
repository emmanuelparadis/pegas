\name{MMD}
\alias{MMD}
\title{Mismatch Distribution}
\description{
  This function draws a histogram of the frequencies of pairwise
  distances from a set of DNA sequences.
}
\usage{
MMD(x, xlab = "Distance", main = "", rug = TRUE,
    legend = TRUE, lcol = "blue", ...)

}
\arguments{
  \item{x}{a set of DNA sequences (object of class \code{"DNAbin"}).}
  \item{xlab}{the label for the x-axis.}
  \item{main}{the title (none by default).}
  \item{rug}{a logical specifying whether to add a rug of the pairwise
    distances on the horizontal axis (see \code{\link[graphics]{rug}}).}
  \item{legend}{a logical specifying whether to draw a legend.}
  \item{lcol}{the colour used for the empirical density curve.}
  \item{\dots}{further arguments passed to \code{hist}.}
}
\details{
  Currently only the observed distribution of pairwise distances are
  plotted, as well as an empirical density estimate. Rogers and
  Harpending (1992) give formulae to compute the expected curves under
  some demographic models.
}
\author{Emmanuel Paradis}
\references{
  Rogers, A. R. and Harpending, H. (1992) Population growth makes waves
  in the distribution of pairwise genetic-differences. \emph{Molecular
    Biology and Evolution}, \bold{9}, 552--569.
}
\examples{
data(woodmouse)
MMD(woodmouse, col = "grey")
MMD(woodmouse, breaks = 20, legend = FALSE)
}
\keyword{hplot}
