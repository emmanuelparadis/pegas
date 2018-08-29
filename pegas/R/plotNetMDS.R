## plotNetMDS.R (2018-08-29)

##   Plot Networks With MDS Layout

## Copyright 2017-2018 Emmanuel Paradis

## This file is part of the R-package `pegas'.
## See the file ../DESCRIPTION for licensing issues.

plotNetMDS <- function(net, d, k = 2, show.mutation = FALSE, col = NULL, font = 2, cex = 1)
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
            }
        }
    }

    mds <- cmdscale(d, k)

    if (is.null(col)) col <- rgb(0, 1, 0, 0.5)
    mat <- rbind(net[, ], attr(net, "alt"))
    a <- mat[, 1]
    b <- mat[, 2]
    if (k == 2) {
        plot.default(mds, type = "n", xlab = "MDS axis 1", ylab = "MDS axis 2")
        segments(mds[a, 1], mds[a, 2], mds[b, 1], mds[b, 2], col = col)
        if (show.mutation)
            .labelSegmentsHaploNet(mds[, 1], mds[, 2], mat[, 1:2], mat[, 3], 1, lwd = 1, col.link = "", 3)
        text(mds, labels = rownames(mds), font = font, cex = cex)
    } else { # k == 3
        rgl::clear3d()
        for (i in 1:nrow(mat)) {
            a <- mat[i, 1:2]
            rgl::segments3d(mds[a, 1:3], col = col)
            if (show.mutation) rgl::text3d(colMeans(mds[a, 1:3]), texts = mat[i, 3], col = "grey")
        }
        rgl::text3d(mds, texts = rownames(mds), font = font, cex = cex)
        rgl::axes3d()
        rgl::mtext3d("MDS axis 1", "X-+", line = 1)
        rgl::mtext3d("MDS axis 2", "Y-+", line = 1)
        rgl::mtext3d("MDS axis 3", "Z+-", line = 1)
    }
}
