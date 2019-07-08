## heterozygosity.R (2019-07-08)

##   Heterozygosity at a Locus Using Gene Frequencies

## Copyright 2002-2018 Emmanuel Paradis

## This file is part of the R-package `pegas'.
## See the file ../DESCRIPTION for licensing issues.

H <- function(x, ...) UseMethod("H")

H.loci <- function(x, variance = FALSE, observed = FALSE, ...)
{
    if (any(.checkPloidy(x) != 2))
        stop("H() requires diploid data")
    n <- nrow(x)
    LOCI <- attr(x, "locicol")
    nloci <- length(LOCI)
    s <- summary(x)
    class(s) <- NULL
    a <- 2*n/(2*n - 1)
    ## get allele sample size for all loci:
    na <- unlist(lapply(s, function(x) sum(x$allele, na.rm = TRUE)))
    ## get sum of pi_^2 for all loci:
    p <- vector("list", nloci) # to store relative allele frequencies
    sp2 <- numeric(nloci)
    for (j in 1:nloci) {
        p[[j]] <- pj <- s[[j]][[2]]/na[j]
        sp2[j] <- sum(pj * pj, na.rm = TRUE)
    }
    res <- a * (1 - sp2)
    dim(res) <- c(nloci, 1)
    colnames(res) <- "Hs"

    if (variance) {
        b <- 2 * (2*n - 2)
        v <- numeric(nloci)
        for (j in 1:nloci) {
            sp3 <- sum(p[[j]]^3, na.rm = TRUE)
            sp2j2 <- sp2[j]^2
            v[j] <- 2 * (b * (sp3 - sp2j2) + sp2[j] - sp2j2)/a
        }
        res <- cbind(res, Var_Hs = v)
    }

    if (observed) {
        o <- numeric(nloci)
        for (j in 1:nloci) {
            geno <- s[[j]][[1]]
            ag <- strsplit(names(geno), "[/|]")
            homoz <- unlist(lapply(ag, anyDuplicated))
            o[j] <- 1 - sum(geno[homoz > 0], na.rm = TRUE)/sum(geno, na.rm = TRUE)
        }
        res <- cbind(res, Hi = o)
    }

    rownames(res) <- names(x)[LOCI]
    res
}

H.default <- function(x, variance = FALSE, ...)
{
    if (!is.factor(x)) {
        if (is.numeric(x)) {
            n <- sum(x)
            k <- length(x)
            freq <- x/n
        } else x <- factor(x)
    }
    if (is.factor(x)) { # ne pas remplacer par `else'...
        n <- length(x)
        k <- nlevels(x)
        freq <- table(x)/n
    }
    sp2 <- sum(freq^2)
    H <- n * (1 - sp2) / (n - 1)
    if (variance) {
        sp3 <- sum(freq^3)
        var.H <- 2 * (2 * (n - 2) * (sp3 - sp2^2) + sp2 - sp2^2) / (n * (n - 1))
        return(c(H, var.H))
    }
    else return(H)
}

heterozygosity <- function(x, variance = FALSE) H(x, variance)
