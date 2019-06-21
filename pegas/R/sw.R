## sw.R (2019-06-21)

##   Sliding Windows

## Copyright 2018-2019 Emmanuel Paradis

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
    else dimnames(res) <- list(rownames(x), collabs)
    class(res) <- "sw"
    res
}

sw.default <- function(x, width = 100, step = 50, POS = NULL, FUN = mean,
                       out.of.pos = NA_real_, na.rm = TRUE, L = NULL, ...)
{
    if (!is.null(dim(x))) {
        x <- as.vector(x)
        warning("'x' has non-null 'dim': converted as a vector")
    }
    lx <- length(x)
    if (is.null(L)) {
        L <- if (!is.null(POS)) max(POS) else lx
    }
    nout <- ceiling(if (step > width) L/step else 1 + (L - width)/step)
    step <- as.integer(step)
    width <- as.integer(width)
    z <- rep(out.of.pos, L)
    if (!is.null(POS)) {
        if (lx != length(POS))
            stop("lengths of 'x' and 'POS' not the same")
        z[POS] <- x
    } else z[1:lx] <- x
    res <- rep.int(NA_real_, nout)
    a <- 1L
    b <- width
    for (j in 1:nout) {
        res[j] <- FUN(z[a:b], na.rm = na.rm)
        a <- a + step
        b <- a + width - 1L
        if (b > L) b <- L
    }
    a <- seq(1L, by = step, length.out = nout)
    b <- a + width - 1L
    b[b > L] <- L
    collabs <- paste0("[", a, ",", b, "]")
    names(res) <- collabs
    class(res) <- "sw"
    res
}

plot.sw <- function(x, type = "l", xlab = "Position", x.scaling = 1,
                    show.ranges = FALSE, col.ranges = "blue",
                    lty.ranges = 1, lwd.ranges = 1, ...)
{
    BOUNDS <- if (is.matrix(x)) colnames(x) else names(x)
    BOUNDS <- gsub("^\\[|\\]$", "", BOUNDS)
    BOUNDS <- strsplit(BOUNDS, ",")
    LOWER <- as.numeric(sapply(BOUNDS, "[", 1)) / x.scaling
    UPPER <- as.numeric(sapply(BOUNDS, "[", 2)) / x.scaling
    xx <- (UPPER + LOWER)/2
    if (is.matrix(x)) {
        matplot(matrix(xx, nrow(x), ncol(x)), x, type = type, xlab = xlab, ...)
        if (show.ranges) warning("x is a matrix: option 'show.ranges' was ignored")
    } else {
        plot.default(xx, x, type = type, xlab = xlab, ...)
        if (show.ranges) segments(LOWER, x, UPPER, x, col = col.ranges,
                                  lty = lty.ranges, lwd = lwd.ranges)
    }
}
