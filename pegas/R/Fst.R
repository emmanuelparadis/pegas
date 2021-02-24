## Fst.R (2021-02-24)

##   F-Statistics

## Copyright 2009-2021 Emmanuel Paradis

## This file is part of the R-package `pegas'.
## See the file ../DESCRIPTION for licensing issues.

Fst <- function(x, pop = NULL, quiet = TRUE, na.alleles = "")
{
    if (any(.checkPloidy(x) != 2))
        stop("Fst() requires diploid data")

    NAMESX <- names(x)
    if (is.null(pop)) {
        ipop <- which(NAMESX == "population")
        if (!length(ipop)) stop("no 'population' column in x")
    } else {
        if (is.numeric(pop) && length(pop) == 1) {
            ipop <- pop
        } else {
            x$populationforthisanalysis <- factor(pop)
            ipop <- length(x)
        }
    }

    LOCI <- attr(x, "locicol")
    nloci <- length(LOCI)

    ## 'p' is a matrix with alleles as columns and populations
    ##    as rows, and its entries are the counts
    ## 'h' is the same with the number of heterozygotes

    res <- matrix(0, nloci, 3)
    dimnames(res) <- list(NAMESX[LOCI], c("Fit", "Fst", "Fis"))

    for (j in 1:nloci) {
        if (!quiet) cat("\rAnalyzing locus", j, "/", nloci)
        Z <- x[, c(LOCI[j], ipop)]
        Z <- na.omit(Z, na.alleles = na.alleles) # all n's are calculated locus-wise (2018-04-20)
        N <- nrow(Z)
        nBYpop <- tabulate(Z$pop)
        r <- length(nBYpop) # number of pops
        nbar <- N/r
        nC <- if (r == 1) 0 else (N - sum(nBYpop^2)/N)/(r - 1)
        ALLELES <- getAlleles(Z)[[1]]
        h <- p <- matrix(0, r, length(ALLELES))
        for (i in 1:r) {
            s <- summary(Z[as.integer(Z$pop) == i, ])[[1]] # levels are preserved
            allel <- names(s$allele)
            genot <- names(s$genotype)
            p[i, ] <- s$allele
            for (k in seq_along(allel)) {
                for (l in seq_along(genot)) {
                    ag <- unlist(strsplit(genot[l], "/"))
                    if (sum(ag %in% allel[k]) == 1)
                        h[i, k] <- h[i, k] + s$genotype[l]
                }
            }
        }
        ptild <- p/(2 * nBYpop)
        pbar <- colSums(p)/(2 * N) # for each allele in the locus
        s2 <- if (r == 1) rep(0, length(ALLELES)) else colSums(nBYpop * (ptild - rep(pbar, each = r))^2)/((r - 1) * nbar)
        hbar <- colSums(h)/N # id.
        A <- pbar * (1 - pbar) - (r - 1) * s2/r
        a <- nbar * (s2 - (A - hbar/4)/(nbar - 1))/nC
        b <- nbar * (A - (2*nbar - 1) * hbar/(4*nbar))/(nbar - 1)
        c <- hbar/2
        res[j, 1] <- if (r == 1) NA else 1 - sum(c)/sum(a + b + c)
        res[j, 2] <- if (r == 1) NA else sum(a)/sum(a + b + c)
        res[j, 3] <- 1 - sum(c)/sum(b + c)
    }
    if (!quiet) cat("... Done.\n")
    res
}

Rst <- function(x, pop = NULL, quiet = TRUE, na.alleles = "")
{
    foo <- function(gr, g1, g2) {
        o <- outer(gr, gr, function(x, y)
            x == g1 & y == g2 | x == g2 & y == g1)
        o[lower.tri(o)] | o[upper.tri(o)]
    }

    NAMESX <- names(x)
    if (is.null(pop)) {
        pop <- x$population
        if (is.null(pop)) stop("no 'population' column in x")
    } else {
        if (length(pop) == 1) {
            pop <- x[[pop]]
        } else {
            pop <- factor(pop)
            if (length(pop) != nrow(x))
                stop("length of 'pop' does not match number of rows in 'x'")
        }
    }
    pop <- as.integer(pop)
    LOCI <- attr(x, "locicol")
    nloci <- length(LOCI)
    res <- numeric(nloci)
    names(res) <- NAMESX[LOCI]
    Z <- loci2alleles(x)
    storage.mode(Z) <- "double"
    if (all(is.na(Z))) stop("no allele interpretable as number of repeats")

    ploidy <- .checkPloidy(x)

    K <- length(unique(pop)) # Nb of pops
    one2K <- 1:K
    expr1 <- 2 / (K * (K - 1))

    gr.diploid <- c(pop, pop)

    with.na.alleles <- !identical(na.alleles, "")

    b <- 0L
    for (j in 1:nloci) {
        if (!quiet) cat("\rAnalyzing locus", j, "/", nloci)
        a <- b + 1L
        b <- a + ploidy[j] - 1L
        y <- Z[, a:b]
        dim(y) <- NULL
        gr <- if (ploidy[j] == 2) gr.diploid else rep(pop, ploidy[j])

        ## drop missing data (null alleles, ...)
        del <- is.na(y)
        if (with.na.alleles) del <- y %in% na.alleles | del
        if (any(del)) {
            y <- y[!del]
            gr <- gr[!del]
        }

        sqdf <- as.dist(outer(y, y, "-")^2) # squared differences
        Nbypop <- tabulate(gr, K)

        SW <- 0
        for (i in one2K) {
            if (!Nbypop[i]) next
            o <- foo(gr, i, i)
            SW <- SW + sum(sqdf[o]) / sum(o) # sum(o) is equal to (n*(n - 1))/2
        }
        SW <- SW / sum(Nbypop > 0)

        SB <- 0
        for (i1 in 1:(K - 1)) {
            if (!Nbypop[i1]) next
            for (i2 in (i1 + 1):K) {
                if (!Nbypop[i2]) next
                o <- foo(gr, i1, i2)
                SB <- SB + sum(sqdf[o]) / sum(o)
            }
        }
        SB <- SB * expr1

        nbar <- mean(Nbypop)
        denom <- sum(Nbypop) - 1
        Sbar <- SW * (nbar - 1) / denom + SB * (nbar * (K - 1)) / denom
        res[j] <- (Sbar - SW) / Sbar
    }
    if (!quiet) cat("... Done.\n")
    res
}
