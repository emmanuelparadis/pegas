## Fst.R (2018-07-07)

##   F-Statistics

## Copyright 2009-2018 Emmanuel Paradis

## This file is part of the R-package `pegas'.
## See the file ../DESCRIPTION for licensing issues.

Fst <- function(x, pop = NULL, quiet = TRUE)
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
        Z <- na.omit(Z) # all n's are calculated locus-wise (2018-04-20)
        N <- nrow(Z)
        nBYpop <- tabulate(Z$pop)
        r <- length(nBYpop) # number of pops
        nbar <- N/r
        nC <- (N - sum(nBYpop^2)/N)/(r - 1)
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
        s2 <- colSums(nBYpop * (ptild - rep(pbar, each = r))^2)/((r - 1) * nbar)
        hbar <- colSums(h)/N # id.
        A <- pbar * (1 - pbar) - (r - 1) * s2/r
        a <- nbar * (s2 - (A - hbar/4)/(nbar - 1))/nC
        b <- nbar * (A - (2*nbar - 1) * hbar/(4*nbar))/(nbar - 1)
        c <- hbar/2
        res[j, 1] <- 1 - sum(c)/sum(a + b + c)
        res[j, 2] <- sum(a)/sum(a + b + c)
        res[j, 3] <- 1 - sum(c)/sum(b + c)
    }
    if (!quiet) cat("... Done.\n")
    res
}
