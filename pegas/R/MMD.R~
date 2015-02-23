## MMD.R (2009-05-11)

##   Mismatch Distribution

## Copyright 2009 Emmanuel Paradis

## This file is part of the R-package `pegas'.
## See the file ../COPYING for licensing issues.

MMD <- function(x, xlab = "Distance", main = "", rug = TRUE,
                legend = TRUE, lcol = "blue", ...)
{
    d <- dist.dna(x, "N")
    hist(d, xlab = xlab, main = main, freq = FALSE, ...)
    lines(density(d, bw = 2), col = lcol)
    if (rug) rug(d)
    if (legend)
        legend("topleft", "Empirical", lty = 1, col = lcol, bg = "white", bty = "n")
}
