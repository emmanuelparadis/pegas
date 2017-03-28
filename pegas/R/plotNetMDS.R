## plotNetMDS.R (2017-03-28)

##   Plot Networks With MDS Layout

## Copyright 2017 Emmanuel Paradis

## This file is part of the R-package `pegas'.
## See the file ../DESCRIPTION for licensing issues.

plotNetMDS <- function(net, d, k = 2, col = NULL, font = 2, cex = 1)
{
    if (! k %in% 2:3) stop("'k' should be 2 or 3")
    n <- attr(d, "Size")
    if (k == 3) {
        if (n == 3) {
            k <- 2
            warning("you set k = 3 but there are only 3 observations: k was set to 2")
        } else {
            if (!requireNamespace("rgl", quietly = TRUE)) {
                k <- 2
                warning("package 'rgl' not available: k was set to 2")
            } else {
                require(rgl)
            }
        }
    }

    mds <- cmdscale(d, k)

    if (is.null(col)) col <- rgb(0, 1, 0, 0.5)
    mat <- rbind(net[, ], attr(net, "alt"))
    if (k == 2) {
        plot.default(mds, type = "n", xlab = "MDS axis 1", ylab = "MDS axis 2", xaxt = "n")
        for (i in 1:nrow(mat)) {
            a <- mat[i, 1]
            b <- mat[i, 2]
            segments(mds[a, 1], mds[a, 2], mds[b, 1], mds[b, 2], col = col)
        }
        text(mds, labels = rownames(mds), font = font, cex = cex)
    } else { # k == 3
        clear3d()
        for (i in 1:nrow(mat)) {
            a <- mat[i, 1:2]
            segments3d(mds[a, 1:3], col = col)
        }
        text3d(mds, texts = rownames(mds), font = font, cex = cex)
    }
}
