## site.spectrum.R (2011-07-19)

##   Site Frequency Spectrum

## Copyright 2009-2011 Emmanuel Paradis

## This file is part of the R-package `pegas'.
## See the file ../COPYING for licensing issues.

site.spectrum <- function(x, folded = TRUE, outgroup = 1)
{
    if (is.list(x)) x <- as.matrix(x)
    n <- dim(x)[1]
    if (n == 1 || is.vector(x))
        stop("only one sequence in the data set")
    if (folded) {
        more.than.two.states <- 0L
        spectrum <- integer(floor(n/2))
        ss <- seg.sites(x)
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
        ambiguous.outgroup.state <- 0L
        spectrum <- integer(n - 1)
        ss <- seg.sites(x[-outgroup, ])
        for (i in ss) {
            anc <- x[outgroup, i, drop = TRUE]
            if (!anc %in% as.raw(c(24, 40, 72, 136)))
                ambiguous.outgroup.state <- ambiguous.outgroup.state + 1L
            else {
                j <- sum(x[-outgroup, i, drop = TRUE] == anc)
                spectrum[j] <- spectrum[j] + 1L
            }
        }
        if (ambiguous.outgroup.state)
            warning(paste(ambiguous.outgroup.state,
                          "sites with ambiguous state were ignored"))
    }
    class(spectrum) <- "spectrum"
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
