## nuc.div.R (2018-07-05)

##   Nucleotide and Haplotype Diversity

## Copyright 2009-2018 Emmanuel Paradis

## This file is part of the R-package `pegas'.
## See the file ../DESCRIPTION for licensing issues.

nuc.div <- function(x, ...) UseMethod("nuc.div")

nuc.div.DNAbin <- function(x, variance = FALSE, pairwise.deletion = FALSE, ...)
{
    if (pairwise.deletion && variance) {
        warning("cannot compute the variance of nucleotidic diversity\nwith pairwise deletion: try 'pairwise.deletion = FALSE' instead.")
        variance <- FALSE
    }
    if (is.list(x)) x <- as.matrix(x)
    d <- dim(x)
    n <- d[1]
    ans <- sum(dist.dna(x, "RAW", pairwise.deletion = pairwise.deletion))/
        (n*(n - 1)/2)
    if (variance) {
        var <- (n + 1)*ans/(3*(n - 1)*d[2]) + 2*(n^2 + n + 3)*ans^2/(9*n*(n - 1))
        ans <- c(ans, var)
    }
    ans
}

nuc.div.haplotype <- function(x, variance = FALSE, pairwise.deletion = FALSE, ...)
{
    if (pairwise.deletion && variance) {
        warning("cannot compute the variance of nucleotidic diversity\nwith pairwise deletion: try 'pairwise.deletion = FALSE' instead.")
        variance <- FALSE
    }
    D <- dist.dna(x, "RAW", pairwise.deletion = pairwise.deletion)
    f <- sapply(attr(x, "index"), length)
    n <- sum(f)
    ff <- outer(f, f)
    ff <- ff[lower.tri(ff)]
    ans <- 2 * sum(ff * D) / (n * (n - 1)) # scale with n(n - 1) to have the same results than with nuc.div.DNAbin()
    if (variance) {
        var <- (n + 1)*ans/(3*(n - 1)*dim(x)[2]) + 2*(n^2 + n + 3)*ans^2/(9*n*(n - 1))
        ans <- c(ans, var)
    }
    ans
}

hap.div <- function(x, ...) UseMethod("hap.div")

hap.div.haplotype <- function(x, variance = FALSE, method = "Nei", ...)
{
    f <- sapply(attr(x, "index"), length)
    n <- sum(f)
    p <- f/n
    sump2 <- sum(p^2)
    n1 <- n - 1L
    res <- (1 - sump2) * n / n1
    if (variance) {
        tmp <- sump2^2
        sump3 <- sum(p^3)
        var <- (sump2 - tmp + 4 * n1 * (sump3 - tmp))/(n * n1)
        res <- c(res, var)
    }
    res
}

hap.div.DNAbin <- function(x, variance = FALSE, method = "Nei", ...)
{
    h <- haplotype(x)
    hap.div.haplotype(h, variance = variance, method = method, ...)
}
