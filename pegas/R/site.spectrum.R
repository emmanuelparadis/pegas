## site.spectrum.R (2019-01-15)

##   Site Frequency Spectrum

## Copyright 2009-2019 Emmanuel Paradis

## This file is part of the R-package `pegas'.
## See the file ../DESCRIPTION for licensing issues.

site.spectrum <- function(x, ...) UseMethod("site.spectrum")

site.spectrum.DNAbin <- function(x, folded = TRUE, outgroup = 1, ...)
{
    if (is.list(x)) x <- as.matrix(x)
    n <- dim(x)[1]
    if (n == 1 || is.null(n))
        stop("only one sequence in data set")

    ss <- seg.sites(x)

    if (folded) {
        more.than.two.states <- 0L
        spectrum <- integer(floor(n/2))
        for (i in ss) {
            bfi <- base.freq(x[, i], freq = TRUE)
            bfi <- bfi[bfi > 0]
            if (length(bfi) > 2)
                more.than.two.states <- more.than.two.states + 1L
            else {
                j <- min(bfi)
                spectrum[j] <- spectrum[j] + 1L
            }
        }
        if (more.than.two.states)
            warning(paste(more.than.two.states,
                          "sites with more than two states were ignored"))
    } else { # unfolded spectrum
        anc <- x[outgroup, ss, drop = TRUE]
        outgroup.state <- anc %in% as.raw(c(24, 40, 72, 136))
        if (s <- sum(!outgroup.state)) {
            warning(paste(s, "sites with ambiguous state were ignored"))
            ss <- ss[outgroup.state]
        }
        spectrum <- apply(x[, ss], 2, function(y) sum(y[outgroup] != y[-outgroup]))
        spectrum <- tabulate(spectrum, n - 1)
    }
    class(spectrum) <- "spectrum"
    attr(spectrum, "sample.size") <- n
    attr(spectrum, "folded") <- folded
    spectrum
}

plot.spectrum <- function(x, col = "red", main = NULL, ...)
{
    if (is.null(main)) {
        main <- "Site Frequency Spectrum"
        main <- if (attr(x, "folded")) paste("Folded", main) else paste("Unfolded", main)
    }

    barplot(as.numeric(x), names.arg = 1:length(x), col = col, main = main, ...)
}

site.spectrum.loci <- function(x, folded = TRUE, ancestral = NULL, ...)
{
    nonsnp <- !is.snp(x)
    p <- length(nonsnp)
    if (!folded) {
        if (is.null(ancestral))
            stop("need a vector of ancestral alleles to compute the unfolded spectrum")
        if (length(ancestral) != p)
            stop("length of 'ancestral' not equal to number of loci")
    }
    if (any(nonsnp)) {
        if (all(nonsnp)) stop("no SNP loci in data set")
        d <- attr(x, "locicol")[nonsnp]
        ld <- length(d)
        x <- x[, -d]
        warning(paste(ld, "non-SNP loci were dropped"))
        p <- p - ld
        if (!folded) ancestral <- ancestral[-d]
    }
    n <- nrow(x)
    s <- summary(x) # works in all situations
    if (folded) {
        f <- sapply(s, function(x) min(x[[2]]), USE.NAMES = FALSE)
        res <- tabulate(f, floor(n/2))
    } else {
        f <- numeric(p)
        for (i in 1:p) f[i] <- s[[i]][[2]][ancestral[i]]
        res <- tabulate(f, n - 1)
    }
    attr(res, "sample.size") <- n
    attr(res, "folded") <- folded
    class(res) <- "spectrum"
    res
}
