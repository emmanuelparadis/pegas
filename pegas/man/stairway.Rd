\name{stairway}
\alias{stairway}
\alias{plot.stairway}
\alias{lines.stairway}
\title{The Stairway Plot}
\description{
  This function fits a model of population change using the site
  frequency spectrum (SFS). The default assumes \eqn{\Theta=1}{THETA =
    1}. A model of population change estimates the temporal changes in
  \eqn{\Theta}{THETA} with respect to the value of this parameter at
  present time. The model is specified by the user with the option
  \code{epoch}.
}
\usage{
stairway(x, epoch = NULL, step.min = 1e-6, step.max = 1e-3)
\method{plot}{stairway}(x, type = "S", xlab = "Coalescent intervals",
             ylab = expression(Theta), \dots)
\method{lines}{stairway}(x, type = "S", \dots)
}
\arguments{
  \item{x}{an object of class \code{\link{site.spectrum}} or of class
    \code{stairway}.}
  \item{epoch}{an optional vector of integers giving the periods of time
    (or epochs) with distinct \eqn{\Theta}{THETA}.}
  \item{step.min}{a single numeric value giving the smallest step size
    used during optimization.}
  \item{step.max}{id. for the largest step size (see
    \code{\link[stats]{nlminb}}).}
  \item{type}{the type of lines.}
  \item{xlab, ylab}{the default labels on the axes.}
  \item{\dots}{further arguments passed to other methods.}
}
\details{
  The basic method implemented in this function is similar to Polanski
  and Kimmel (2003). The temporal model with ``epochs'' is from Liu and
  Fu (2015).
}
\value{
  By default, a single numeric value with the null deviance. If
  \code{epoch} is used, a list with the following components:

\item{estimates}{the maximum likelihood estimates.}
\item{deviance}{the deviance of the fitted model.}
\item{null.deviance}{the deviance of the null model.}
\item{LRT}{the likelihood-ratio test comparing the null and the fitted models.}
\item{AIC}{the Akaike information criterion of the fitted model.}
}
\references{
Liu, X. M. and Fu, Y. X. (2015) Exploring population size changes using SNP frequency spectra. \emph{Nature Genetics}, \bold{47}, 555--559.

Polanski, A. and Kimmel, M. (2003) New explicit expressions for relative frequencies of single-nucleotide polymorphisms with application to statistical inference on population growth. \emph{Genetics}, \bold{165}, 427--436.
}
\author{Emmanuel Paradis}
\seealso{\code{\link{site.spectrum}}, \code{\link[stats]{nlminb}}}
\examples{
data(woodmouse)
sp <- site.spectrum(woodmouse)
stairway(sp, c(1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2))
}
\keyword{models}
