\name{theta.tree}
\alias{theta.tree}
\title{
  Population Parameter THETA Using Genealogy
}
\description{
  This function estimates the population parameter \eqn{\theta}{THETA}
  from a genealogy (coded a as phylogenetic tree) under the coalescent.
}
\usage{
theta.tree(phy, theta, fixed = FALSE, analytical = TRUE, log = TRUE)
}
\arguments{
  \item{phy}{an object of class \code{"phylo"}.}
  \item{theta}{a numeric vector.}
  \item{fixed}{a logical specifying whether to estimate \code{theta}
    (the default), or to return the likelihoods for all values in
    \code{theta}.}
  \item{analytical}{a logical specifying whether to use analytical
    formulae to estimate \code{theta} and its standard-error. If
    \code{FALSE}, a numerical optimisation of the likelihood is
    performed (this option is ignored if \code{fixed = TRUE})}.
  \item{log}{a logical specifying whether to return the likelihoods on a
    log scale (the default); ignored if \code{fixed = FALSE}.}
}
\details{
  The tree \code{phy} is considered as a genealogy, and therefore should
  be ultrametric. By default, \eqn{\theta}{THETA} is estimated by
  maximum likelihood and the value given in \code{theta} is used as
  starting value for the minimisation function (if several values are
  given as a vector the first one is used). If \code{fixed = TRUE},
  then the [log-]likelihood values are returned corresponding to each
  value in \code{theta}.

  The present implementation does a numerical optimisation of the
  log-likelihood function (with \code{\link[stats]{nlminb}}) with the
  first partial derivative as gradient. It is possible to solve the
  latter and have a direct analytical MLE of \eqn{\theta}{THETA} (and
  its standard-error), but this does not seem to be faster.
}
\value{
  If \code{fixed = FALSE}, a list with two elements:

  \item{theta}{the maximum likelihood estimate of \eqn{\theta}{THETA};}
  \item{logLik}{the log-likelihood at its maximum.}

  If \code{fixed = TRUE}, a numeric vector with the [log-]likelihood
  values.
}
\references{
  Kingman, J. F. C. (1982) The coalescent. \emph{Stochastic Processes
  and their Applications}, \bold{13}, 235--248.

  Kingman, J. F. C. (1982) On the genealogy of large
  populations. \emph{Journal of Applied Probability}, \bold{19A},
  27--43.

  Wakeley, J. (2009) \emph{Coalescent Theory: An Introduction.}
  Greenwood Village, CO: Roberts and Company Publishers.
}
\author{Emmanuel Paradis}
\seealso{
  \code{\link{theta.h}}, \code{\link{theta.s}}, \code{\link{theta.k}}
}
\examples{
tr <- rcoal(50) # assumes theta = 2
(o <- theta.tree(tr))
theta.tree(tr, 10, analytical = FALSE) # uses nlminb()
## profile log-likelihood:
THETA <- seq(0.5, 4, 0.01)
logLikelihood <- theta.tree(tr, THETA, fixed = TRUE)
plot(THETA, logLikelihood, type = "l")
xx <- seq(o$theta - 1.96 * o$se, o$theta + 1.96 * o$se, 0.01)
yy <- theta.tree(tr, xx, fixed = TRUE)
polygon(c(xx, rev(xx)), c(yy, rep(0, length(xx))),
        border = NA, col = "lightblue")
segments(o$theta, 0, o$theta, o$logLik, col = "blue")
abline(v = 2, lty = 3)
legend("bottomright", legend = expression("log-likelihood",
       "True " * theta, hat(theta) * " (MLE)", "95\%\ conf. interv."),
       lty = c(1, 3, 1, 1), lwd = c(1, 1, 1, 15),
       col = c("black", "black", "blue", "lightblue"))
}
\keyword{models}
