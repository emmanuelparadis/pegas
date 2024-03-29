## theta.R (2022-01-04)

##   Population Parameter THETA

## theta.h: using homozygosity
## theta.k: using expected number of alleles
## theta.s: using segregating sites in DNA sequences
## theta.tree: using a genealogy
## theta.tree.hetero: using a genealogy with heterochronous dates
## theta.msat: using micro-satellites

## Copyright 2002-2022 Emmanuel Paradis

## This file is part of the R-package `pegas'.
## See the file ../DESCRIPTION for licensing issues.

theta.h <- function(x, standard.error = FALSE)
{
    HE <- H(x, variance = TRUE)
    sdH <- HE[2]
    HE <- HE[1]
    f <- function(th) HE - th * (1 + (2 * (1 + th)) / ((2 + th) * (3 + th)))
    th <- uniroot(f, interval = c(0, 1))$root
    if (standard.error) {
        SE <- (2 + th)^2 * (2 + th)^3 * sdH /
            HE^2 * (1 + th) * ((2 + th) * (3 + th) * (4 + th) + 10 * (2 + th) + 4)
        th <- c(th, SE)
    }
    th
}

theta.k <- function(x, n = NULL, k = NULL)
{
    if (is.null(n)) {
        if (!is.factor(x)) {
            if (is.numeric(x)) {
                n <- sum(x)
                k <- length(x)
            }
            else x <- factor(x)
        }
        if (is.factor(x)) { # ne pas remplacer par `else'...
            n <- length(x)
            k <- nlevels(x)
        }
    }
    f <- function(th) th * sum(1 / (th + (0:(n - 1)))) - k
    uniroot(f, interval = c(1e-8, 100))$root
}

theta.s <- function(x, ...) UseMethod("theta.s")

theta.s.default <- function(x, n, variance = FALSE, ...)
{
    b <- 1:(n - 1)
    a1 <- sum(1/b)
    th <- x/a1
    if (variance) {
        a2 <- sum(1/b^2)
        var.th <- (th * a1 + th^2 * a2) / a1^2 # fix by Carsten Wiuf (2022-01-04)
        th <- c(th, var.th)
    }
    th
}

theta.s.DNAbin <- function(x, variance = FALSE, ...)
{
    s <- length(seg.sites(x))
    n <- nrow(x)
    theta.s.default(s, n, variance = variance)
}

theta.tree <-
    function(phy, theta, fixed = FALSE, analytical = TRUE, log = TRUE)
{
    ## coalescent intervals from the most recent to the oldest one:
    x <- diff(c(0, sort(branching.times(phy))))
    k <- length(phy$tip.label):2 # c(n, ..., 2)
    K <- length(k) # n - 1
    tmp <- (k * (k - 1))/2 # choose(k, 2)
    tmp2 <- sum(x * tmp)
    sltmp <- sum(log(tmp))
    if (fixed) {
        res <- sltmp - K * log(theta) - tmp2/theta
        if (!log) res <- exp(res)
    } else {
        if (analytical) {
            theta <- tmp2/K
            se <- sqrt(-1/(K/theta^2 - 2 * tmp2/theta^3))
            logLik <- sltmp - K * log(theta) - tmp2/theta
            res <- list(theta = theta, se = se, logLik = logLik)
        } else {
            minusLogLik <- function(theta) # vectorized on 'theta'
                -(sltmp - K*log(theta) - tmp2/theta)
            gr <- function(theta) K/theta - tmp2/theta^2
            out <- nlminb(theta[1], minusLogLik, gr,
                          lower = .Machine$double.eps, upper = Inf)
            res <- list(theta = out$par, logLik = -out$objective)
        }
    }
### alternative version based on L-BFGS-B
###out <- optim(theta[1], minusLogLik, gr, method = "L-BFGS-B",
###             lower = .Machine$double.eps, upper = Inf,
###             hessian = TRUE)
###res <- list(theta = out$par, se = sqrt(1/out$hessian[, ]),
###            logLik = -out$value)
### I prefered nlminb() because it is slightly faster and in most cases
### the hessian-based estimate of SE(theta) are not needed
    res
}

theta.tree.hetero <- function(phy, theta, fixed = FALSE, log = TRUE)
{
    n <- length(phy$tip.label)
    ROOT <- n + 1L
    times <- dist.nodes(phy)[, ROOT]
    ## rescale so the most recent sample has t=0 and the root has t=max(t):
    times <- max(times) - times
    o <- order(times)
    isSamp <- o <= n
    isCoal <- !isSamp
    trans <- as.integer(isSamp)
    trans[isCoal] <- -1L
    x <- cumsum(trans)
    coal.times <- c(0, times[o][isCoal]) # sorted
    x.coal <- x[which(isCoal) - 1] # number of lineages before the coalescence
    coal.ints <- diff(coal.times)
    tmp <- x.coal * (x.coal - 1)/2
    if (fixed) {
        res <- sum(lchoose(x.coal, 2) - log(theta) - tmp * coal.ints/theta)
        if (log) res else exp(res)
    } else {
        theta <- sum(tmp * coal.ints)/(n - 1)
        se <- sqrt(-1/((n - 1)/theta^2 - 2 * sum(tmp * coal.ints)/theta^3))
        logLik <- sum(lchoose(x.coal, 2) - log(theta) - tmp * coal.ints/theta)
        list(theta = theta, se = se, logLik = logLik)
    }
}

theta.msat <- function(x)
{
    s <- summary(x)
    getThetas <- function(x) {
        wi <- x$allele
        n <- sum(wi) # number of alleles
        ai <- as.numeric(names(x$allele))
        abar <- weighted.mean(ai, wi)
        fi <- wi/n
        H0 <- (n * sum(fi^2) - 1) / (n - 1)
        xbar <- mean(fi)
        c(2 * sum(wi * (ai - abar)^2)/(n - 1), # theta_va
          0.5 * (1/H0^2 - 1), # theta_h
          1/(8 * xbar^2) - .5) # theta_xbar
    }
    res <- t(sapply(s, getThetas))
    colnames(res) <- c("theta.v", "theta.h", "theta.x")
    res
}
