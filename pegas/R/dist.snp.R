## dist.snp.R (2017-12-18)

##   Allelic Sharing Distance

## Copyright 2017 Emmanuel Paradis

## This file is part of the R-package `pegas'.
## See the file ../DESCRIPTION for licensing issues.

## if all loci are unphased, diploid, and biallelic (e.g., strict SNP)
dist.snp <- function(x, scaled = TRUE)
{
    n <- nrow(x)
    labs <- row.names(x)
    class(x) <- NULL # makes things MUCH faster
    locicol <- attr(x, "locicol")
    p <- length(locicol)
    y <- matrix(0L, n, p)
    for (j in 1:p) {
        tmp <- x[[locicol[j]]]
        attributes(tmp) <- NULL
        y[, j] <- tmp
    }
    D <- numeric(n*(n - 1)/2)
    k <- 1L
    for (i in 1:(n - 1)) {
        a <- y[i, ]
        for (j in (i + 1):n) {
            D[k] <- sum(abs(a - y[j, ]))
            k <- k + 1L
        }
    }
    if (scaled) D <- D/p
    attr(D, "Size") <- n
    attr(D, "Labels") <- labs
    attr(D, "Diag") <- attr(D, "Upper") <- FALSE
    attr(D, "call") <- match.call()
    class(D) <- "dist"
    D
}

dist.asd <- function(x, scaled = TRUE)
{
    labs <- row.names(x)
    locicol <- attr(x, "locicol")
    nloc <- length(locicol)
    n <- nrow(x)
    class(x) <- NULL # makes things MUCH faster
    foo <- function(x) {
        geno <- levels(x)
        ng <- length(geno)
        alle <- strsplit(geno, "[/|]")
        ualle <- lapply(alle, unique.default)
        if (length(ualle) == 1) return(0)
        m <- matrix(0, ng, ng)
        for (i in 1:(ng - 1)) {
            a <- ualle[[i]]
            for (j in (i + 1):ng)
                m[i, j] <- m[j, i] <- 2 - sum(outer(a, ualle[[j]], "=="))
        }
        x <- unclass(x)
        d <- numeric(n*(n - 1)/2)
        k <- 1L
        for (i in 1:(n - 1)) {
            a <- x[i]
            for (j in (i + 1):n) {
                d[k] <- m[a, x[j]]
                k <- k + 1L
            }
        }
        d
    }
    D <- 0
    for (j in locicol) D <- D + foo(x[[j]])
    if (scaled) D <- D/nloc
    attr(D, "Size") <- n
    attr(D, "Labels") <- labs
    attr(D, "Diag") <- attr(D, "Upper") <- FALSE
    attr(D, "call") <- match.call()
    class(D) <- "dist"
    D
}
