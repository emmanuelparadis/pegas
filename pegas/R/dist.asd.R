## dist.asd.R (2022-05-16)

##   Allelic Sharing Distance

## Copyright 2017-2022 Emmanuel Paradis

## This file is part of the R-package `pegas'.
## See the file ../DESCRIPTION for licensing issues.

dist.asd <- function(x, scaled = TRUE, pairwise.deletion = FALSE)
{
    labs <- row.names(x)
    locicol <- attr(x, "locicol")
    nloc <- length(locicol)
    n <- nrow(x)
    np <- n*(n - 1)/2

    ## check if all loci are diploid and biallelic
    ploidy <- getPloidy(x)
    if (any(ploidy != 2)) stop("all genotypes must be diploid")
    alleles <- getAlleles(x)
    FAST <- all(lengths(alleles) == 2)
    if (FAST && any(is.phased(x))) x <- unphase(x)

    ## if all genotypes are biallelic (FAST == TRUE), then there can
    ## only three genotypes which are coded 1, 2, 3 in the factor with
    ## 2 being the heterozygote, so the distance between two genotypes
    ## is equal to the difference of their respective codes

    class(x) <- NULL # makes things MUCH faster
    spd <- scaled & pairwise.deletion
    if (spd) NLOC <- integer(np)

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
        ## There should be no NA in 'm'
        x <- unclass(x)
        d <- numeric(np)
        k <- 1L
        for (i in 1:(n - 1)) {
            a <- x[i]
            if (is.na(a)) next
            for (j in (i + 1):n) {
                b <- x[j]
                if (is.na(b)) next
                ## 'a' and 'b' cannot be NA
                d[k] <- m[a, b]
                k <- k + 1L
                if (spd) NLOC[k] <<- NLOC[k] + 1L
            }
        }
        d
    }

    if (FAST) {
        y <- matrix(0L, n, nloc)
        for (j in 1:nloc) {
            tmp <- x[[locicol[j]]]
            attributes(tmp) <- NULL
            y[, j] <- tmp
        }
        D <- numeric(np)
        k <- 1L
        for (i in 1:(n - 1)) {
            a <- y[i, ]
            for (j in (i + 1):n) {
                tmp <- abs(a - y[j, ])
                D[k] <- sum(tmp, na.rm = pairwise.deletion)
                if (spd) NLOC[k] <- sum(!is.na(tmp))
                k <- k + 1L
            }
        }
    } else {
        D <- 0
        for (j in locicol) D <- D + foo(x[[j]])
    }
    if (scaled) {
        if (pairwise.deletion) nloc <- NLOC
        D <- D/nloc
    }
    attr(D, "Size") <- n
    attr(D, "Labels") <- labs
    attr(D, "Diag") <- attr(D, "Upper") <- FALSE
    attr(D, "call") <- match.call()
    class(D) <- "dist"
    D
}
