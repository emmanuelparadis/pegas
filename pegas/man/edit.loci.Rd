\name{edit.loci}
\alias{edit.loci}
\title{Edit Allelic Data with R's Data Editor}
\description{
  This allows to edit a data frame of class \code{"loci"} with R's
  spreadsheet-like data editor.
}
\usage{
\method{edit}{loci}(name, edit.row.names = TRUE, ...)
}
\arguments{
  \item{name}{an object of class \code{"loci"}.}
  \item{edit.row.names}{a logical specifying to allow editing the
    rownames, \code{TRUE} by default (by contrast to data frames).}
  \item{\dots}{further arguments to be passed to or from other methods.}
}
\details{
  This `method' of the generic \code{edit} respects the class and the
  attribute \code{"locicol"} of the allelic data frame.
}
\value{
  A data frame with class \code{c("loci", "data.frame")}.
}
\author{Emmanuel Paradis}
\seealso{
  \code{\link{read.loci}}, \code{\link{summary.loci}}
}
\keyword{IO}
