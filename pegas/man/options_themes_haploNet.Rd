\name{getHaploNetOptions}
\alias{getHaploNetOptions}
\alias{setHaploNetOptions}
\title{Options to Plot haploNet Objects}
\description{
  These functions change the graphical options to plot haplotype networks.
}
\usage{
getHaploNetOptions()
setHaploNetOptions(...)
}
\arguments{
  \item{...}{option(s) and value(s) to be changed (separated by commas
    if several).}
}
\details{
  The options are listed below with their default values. Most of these
  values use the standard \R graphical paramters (see
  \code{\link[graphics]{par}}).

\itemize{
  \item{bg = "transparent"}{the background colour of the plot.}
  \item{labels = TRUE}{whether to show the haplotype labels.}
  \item{labels.cex = 1}{size of the haplotype labels.}
  \item{labels.font = 2}{font of the haplotype labels.}
  \item{link.color = "black"}{colour of the links.}
  \item{link.type = 1}{type of line for the links.}
  \item{link.type.alt = 2}{type of lines for the alternative links.}
  \item{link.width = 1}{line width for the links.}
  \item{link.width.alt = 1}{line width for the alternative links.}
  \item{haplotype.inner.color = "white"}{colour used inside the
    haplotype symbols.}
  \item{haplotype.outer.color = "black"}{colour used for the border of
    the haplotype symbols.}
  \item{mutations.cex = 1}{size of the mutation annotations.}
  \item{mutations.font = 1}{font of the mutation annotations.}
  \item{mutations.frame.background = "#0000FF4D"}{transparent blue.}
  \item{mutations.frame.border = "black"}{}
  \item{mutations.text.color = 1}{colour of the mutation annotations.}
  \item{mutations.arrow.color = "black"}{colour of the arrow pointing to
    the link.}
  \item{mutations.arrow.type = "triangle"}{type of the above arrow.}
  \item{mutations.sequence.color = "#BFBFBF4D"}{transparent grey.}
  \item{mutations.sequence.end = "round"}{possible choices: \code{"round"},
    \code{"butt"}, or \code{"square"} (or alternatively 0, 1, or 2).}
  \item{mutations.sequence.length = 0.3}{the length of the segment
    showing the sequence as fraction of the graphical window.}
  \item{mutations.sequence.width = 5}{thickness of this segment.}
  \item{pie.outer.color = "black"}{}
  \item{pie.inner.segments.color = "black"}{}
  \item{pie.colors.function = rainbow}{functio used to define colours
    for the frequencies.}
  \item{scale.ratio = 1}{}
  \item{show.mutation = 1}{}
}}
\value{
  \code{getHaploNetOptions} returns a list of options. The other
  function returns nothing.
}
\author{Emmanuel Paradis}
\seealso{\code{\link{plot.haploNet}}, \code{\link{mutations}}}
\examples{
getHaploNetOptions()
}
\keyword{hplot}
\keyword{lplot}