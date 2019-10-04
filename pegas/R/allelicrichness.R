## allelicrichness.R (2019-10-03)

##   F-Statistics

## Copyright 2019 Emmanuel Paradis

## This file is part of the R-package `pegas'.
## See the file ../DESCRIPTION for licensing issues.

.eq13.Hurlbert1971 <- function(x, ni)
{
    n <- sum(ni, na.rm = TRUE)
    length(ni) - sum(choose(n - ni, x))/choose(n, x)
}

rarefactionplot <-
    function(x, maxn = nrow(x), type = "l", xlab = "Sample size",
             ylab = "Expected number of alleles", plot = TRUE, ...)
{
    s <- summary.loci(x)
    nms <- names(s)
    s <- lapply(s, "[[", "allele")
    nloci <- length(s)
    xx <- 1:maxn
    YY <- vector("list", nloci)
    for (i in seq_len(nloci))
        YY[[i]] <- sapply(xx, .eq13.Hurlbert1971, ni = s[[i]])
    if (plot) {
        if (!par("ask") && names(dev.cur()) %in% deviceIsInteractive()) {
            par(ask = TRUE)
            on.exit(par(ask = FALSE))
        }
        for (i in seq_len(nloci)) {
            plot(xx, YY[[i]], type = type, xlab = xlab, ylab = ylab,
                 main = nms[i], ...)
        }
    }
    names(YY) <- nms
    invisible(YY)
}

allelicrichness <- function(x, pop = NULL, method = "extrapolation")
{
    method <- match.arg(method, c("extrapolation", "rarefaction", "raw"))
    s <- summary.loci(x)
    s <- lapply(s, "[[", "allele")
    nloci <- length(s)
    bypop <- by(x)
    nbypop <- tabulate(x$population)
    npop <- length(nbypop)
    res <- matrix(NA_real_, npop, nloci)
    switch(method, "extrapolation" = {
        p <- lapply(s, function(x) x/sum(x, na.rm = TRUE))
        for (j in seq_len(nloci)) {
            pj <- p[[j]]
            for (i in 1:npop) {
                tab <- bypop[[j]]
                nalleles <- ncol(bypop[[j]])
                absent <- tab[i, ] == 0
                tmp <- nalleles
                if (any(absent))
                    tmp <- tmp - sum((1 - pj[which(absent)])^nbypop[i])
                res[i, j] <- tmp
            }
        }
    }, "rarefaction" = {
        for (j in seq_len(nloci)) {
            tab <- bypop[[j]]
            for (i in 1:npop)
                res[i, j] <- .eq13.Hurlbert1971(sum(tab[i, ]), tab[i, ])
        }
    }, "raw" = {
        for (j in seq_len(nloci)) {
            tab <- bypop[[j]]
            for (i in 1:npop) res[i, j] <- sum(tab[i, ] > 0)
        }
    })
    dimnames(res) <- list(levels(x$population), names(s))
    res
}

rhost <- function(x, pop = NULL, method = "extrapolation")
{
    R <- allelicrichness(x, method = method)
    Rbar <- apply(R, 2, mean, na.rm = TRUE)
    1 - (Rbar - 1)/(ncol(R) - 1)
}
