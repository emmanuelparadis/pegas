## stairway.R (2019-01-17)

##   Stairway Plot

## Copyright 2019 Emmanuel Paradis

## This file is part of the R-package `pegas'.
## See the file ../DESCRIPTION for licensing issues.

stairway <- function(x, epoch = NULL, step.min = 1e-6, step.max = 1e-3)
### "Your stairway lies on the whispering wind" (R. Plant)
{
    if (!inherits(x, "spectrum"))
        stop("x should be of class \"spectrum\"")
    if (!attr(x, "folded"))
        stop("site spectrum should be folded")

    n <- attr(x, "sample.size")

    ## eq. 12 in Polanski & Kimmel (2003, Genetics)
    ## (vectorized on j = 1, ..., n)
    V <- function(n, j) {
        A <- lfactorial(n)
        B <- lfactorial(n - 1)
        C <- lfactorial(n + j - 1)
        D <- lfactorial(n - j)
        res <- (4*j - 2) * exp(A + B - C - D)
        odd <- as.logical(j %% 2) # if j is odd, then (1 + (-1)^j) = 0
        if (any(odd)) res[odd] <- 0
        res
    }

    ## eqs. 13-15
    ## (vectorized on j = 1, ..., n with w[1] = 0)
    W <- function(n, b, j) {
        w <- numeric(n)
        w[2L] <- 6/(n + 1)
        w[3L] <- 30*(n - 2*b)/((n + 1)*(n + 2))
        for (i in 4:n) {
            k <- i - 2L
            X <- 3 + 2*k
            Y <- k * (n + k + 1)
            B <- (X * (n - 2*b))/Y
            A <- (1 + k) * X * (n - k)/((2*k -1) * Y)
        w[i] <- B * w[i - 1] - A * w[i - 2]
        }
        w[j]
    }

    makeVcache <- function(n) V(n, 1:n)
    makeWcache <- function(n) {
        res <- matrix(0, n, n)
        for (i in 1:n) res[i, ] <- W(n, i, 1:n)
        res
    }

    v.cache <- makeVcache(n)
    w.cache <- makeWcache(n)

    ## CDF of coalescent times -- eq. 3
    e <- function(n, theta) { # THETA can be a vector or a scalar
        k <- n:2
        tmp <- k*(k - 1)
        cumsum(theta/tmp)
    }

    ## probability to observe a locus of size 'b' in the spectrum
    ## among 'n' alleles -- eq. 8 in Polanski & Kimmel
    ## (vectorized on b, with b integers)
    q <- function(n, b, theta) {
        E <- e(n, theta)
        res <- numeric(lb <- length(b))
        for (i in 1:lb)
            res[i] <- sum(E * w.cache[b[i], n:2])/sum(E * v.cache[n:2])
            ##res[i] <- sum(E * w.cache[b[i], 2:n])/sum(E * v.cache[2:n])
        ## res[i] <- sum(E * W(n, b[i], n:2))/sum(E * V(n, n:2))
        res
    }

    ll <- function(theta) {
        b <- 1:length(x)
        ## eq. 22 in Polanski & Kimmel (next 2 lines):
        p <- q(n, b, theta)
        if (any(p < 0)) p <- abs(p)
        if (any(m <- b == n/2)) p[m] <- 2 * p[m]
        ##sum(x * log(p[x > 0])) # eq. 24
        sum(x * log(p)) # eq. 24
    }

    if (!is.null(epoch)) {
        if (length(epoch) != n - 1)
            stop("vector 'epoch' should have ", n - 1, " elements")
        epoch <- as.integer(epoch)
        tmp <- sort(unique(epoch))
        if (!all(diff(tmp) == 1))
            stop("values in 'epoch' should be numbered 1, 2, ... without gap")
        np <- length(tmp) - 1
        dev <- function(p) -2 * ll(c(1, p)[epoch])
        lo <- rep(1e-8, np)
        ip <- rep(1, np)
        ctr <- list(step.min = step.min, step.max = step.max)
        res <- nlminb(ip, dev, lower = lo, control = ctr)
    }

    dev0 <- -2 * ll(1)
    if (is.null(epoch)) return(dev0)
    res <- res[1:2]
    names(res) <- c("estimates", "deviance")
    res$null.deviance  <-  dev0
    chi2  <-  dev0 - res$deviance
    res$LRT <- c(chi2 = chi2, df = np, P.val = 1 - pchisq(chi2, np))
    res$AIC <- res$deviance + 2 * np
    res$epoch <- epoch
    class(res) <- "stairway"
    res
}

plot.stairway <- function(x, type = "S", xlab = "Coalescent intervals",
                          ylab = expression(Theta), ...)
{
    y <- c(1, x$estimates)[x$epoch]
    plot(0:length(y), c(NA, rev(y)), type = type, xaxt = "n",
         xlab = xlab, ylab = ylab, ...)
    xx <- pretty(1:length(y))
    axis(1, at = xx, labels = rev(xx))
}

lines.stairway <- function(x, type = "S", ...)
{
    y <- c(1, x$estimates)[x$epoch]
    lines.default(0:length(y), c(NA, rev(y)), type = type, ...)
}
