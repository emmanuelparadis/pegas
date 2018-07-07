## MMD.R (2018-07-07)

##   Sliding Windows

## Copyright 2018 Emmanuel Paradis

## This file is part of the R-package `pegas'.
## See the file ../DESCRIPTION for licensing issues.

sw <- function(x, width, step, ...) UseMethod("sw")

sw.DNAbin <- function(x, width = 100, step = 50, FUN = GC.content,
                      rowAverage = FALSE, quiet = TRUE, ...)
{
    n <- nrow(x)
    s <- ncol(x)
    nout <- ceiling(if (step > width) s/step else 1 + (s - width)/step)
    step <- as.integer(step)
    width <- as.integer(width)
    out <- FUN(x[1, 1:width])
    if (length(out) > 1) stop("FUN should return a vector of length one")
    if (is.list(out)) stop("FUN should return a vector (not a list)")
    if (rowAverage) {
        res <- rep.int(NA, nout)
    } else {
        res <- rep.int(NA, n * nout)
        dim(res) <- c(n, nout)
        skip <- 0L
    }
    storage.mode(res) <- storage.mode(out)
    a <- 1L
    b <- width
    for (j in 1:nout) {
        if (!quiet) cat("\rSliding window:", j, "/", nout)
        z <- x[, a:b]
        if (rowAverage) {
            res[j] <- FUN(z)
        } else {
            for (i in 1:n)
                res[i + skip] <- FUN(z[i, ])
            skip <- skip + n
        }
        a <- a + step
        b <- a + width - 1L
        if (b > s) b <- s
    }
    if (!quiet) cat("... done.\n")
    a <- seq(1L, by = step, length.out = nout)
    b <- a + width - 1L
    b[b > s] <- s
    collabs <- paste0("[", a, ",", b, "]")
    if (rowAverage) names(res) <- collabs
    else dimnames <- list(rownames(x), collabs)
    res
}
