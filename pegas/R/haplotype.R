## haplotype.R (2019-03-13)

##   Haplotype Extraction, Frequencies, and Networks

## Copyright 2009-2019 Emmanuel Paradis, 2013 Klaus Schliep

## This file is part of the R-package `pegas'.
## See the file ../DESCRIPTION for licensing issues.

rmst <- function(d, B = 100)
{
    if (!is.matrix(d) && !inherits(d, "dist"))
        stop("'d' must be a matrix or a 'dist' object")
    D <- as.matrix(d)
    n <- nrow(D)
    MAT <- matrix(NA_character_, 0, 2)
    for (rep in 1:B) {
        if (rep == B) d <- D
        else {
            ri <- sample(n)
            d <- D[ri, ri]
        }
        d <- as.dist(d)
        rnt <- mst(d)
        m <- rnt[, 1:2]
        dm <- dim(m)
        mat <- attr(rnt, "labels")[m]
        dim(mat) <- dm
        MAT <- rbind(MAT, t(apply(mat, 1, sort)))
    }
    X <- paste(MAT[, 1], MAT[, 2], sep = "\r")
    tab <- table(X)
    if (length(tab) == n - 1) {
        attr(rnt, "weight") <- rep(B, n - 1)
        return(rnt)
    }
    links <- names(tab)
    ## define the alternative links wrt the last MST:
    labs <- attr(rnt, "labels")
    MST.str <- labs[rnt[, 1:2]]
    dim(MST.str) <- c(n - 1, 2)
    MST.str <- apply(MST.str, 1, sort)
    MST.str <- paste(MST.str[1, ], MST.str[2, ], sep = "\r")
    k <- match(links, MST.str)
    nak <- is.na(k)
    tab <- c(tab[match(MST.str, links)], tab[nak])
    alt <- unlist(strsplit(links[nak], "\r"))
    alt <- match(alt, labs)
    alt <- matrix(alt, length(alt)/2, 2, byrow = TRUE)
    alt <- t(apply(alt, 1, sort))
    i <- alt[, 1]
    j <- alt[, 2]
    k <- n*(i - 1) - i*(i - 1)/2 + j - i
    alt <- cbind(alt, d[k])
    colnames(alt) <- c("", "", "step")
    attr(rnt, "alter.links") <- alt
    names(tab) <- gsub("\r", "--", names(tab))
    attr(rnt, "weight") <- tab
    rnt
}

msn <- function(d)
{
    getIandJ <- function(ij, n) {
        ## assumes a lower triangle, so i > j
        ## n must be > 1 (not checked)
        ## ij must be <= (n - 1)*n/2 (not checked too)
        j <- 1L
        N <- n - 1L
        while (ij > N) {
            j <- j + 1L
            N <- N + n - j
        }
        i <- n - (N - ij)
        c(j, i) # return the smaller index first
    }
    if (is.matrix(d)) d <- as.dist(d)
    MST <- mst(d)
    n <- attr(d, "Size")
    if (n < 3) {
        warning("less than 3 observations: function mst() was used")
        return(MST)
    }
    m <- matrix(NA_real_, 0, 3)
    forest <- 1:n
    ud <- sort(unique(d)) # unique distances in increasing order

    i <- 1L
    while (length(unique(forest)) > 1) {
        delta <- ud[i]
        for (iud in which(d == delta)) {
            p <- getIandJ(iud, n)
            f1 <- forest[p[1L]]
            f2 <- forest[p[2L]]
            m <- rbind(m, c(p, delta))
            if (f1 != f2) forest[forest == f2] <- f1
        }
        i <- i + 1L
    }
    ## define the alternative links wrt the MST:
    MST.str <- paste(MST[, 1], MST[, 2], sep = "\r")
    m.str <- paste(m[, 1], m[, 2], sep = "\r")
    k <- match(m.str, MST.str)
    nak <- is.na(k)
    if (any(nak)) {
        alt <- m[nak, , drop = FALSE]
        colnames(alt) <- c("", "", "step")
        m <- m[!nak, , drop = FALSE]
        attr(m, "alter.links") <- alt
    }
    colnames(m) <- c("", "", "step")
    attr(m, "labels") <- attr(d, "Labels")
    class(m) <- "haploNet"
    m
}

mst <- function(d)
{
    getIandJ <- function(ij, n) {
        ## assumes a lower triangle, so i > j
        ## n must be > 1 (not checked)
        ## ij must be <= (n - 1)*n/2 (not checked too)
        j <- 1L
        N <- n - 1L
        while (ij > N) {
            j <- j + 1L
            N <- N + n - j
        }
        i <- n - (N - ij)
        c(j, i) # return the smaller index first
    }
    if (is.matrix(d)) d <- as.dist(d)
    n <- attr(d, "Size")
    if (n < 2) stop("less than 2 observations in distance matrix")
    Nedge <- n - 1L
    m <- matrix(NA_real_, Nedge, 3)
    forest <- 1:n
    o <- order(d)
    p <- getIandJ(o[1L], n)
    m[1, ] <- c(p, d[o[1L]])
    forest[p[2L]] <- forest[p[1L]]
    i <- j <- 2L
    while (j <= Nedge) {
        p <- getIandJ(o[i], n)
        f1 <- forest[p[1L]]
        f2 <- forest[p[2L]]
        if (f2 != f1) {
            m[j, ] <- c(p, d[o[i]])
            forest[forest == f2] <- f1
            j <- j + 1L
        }
        i <- i + 1L
    }
    colnames(m) <- c("", "", "step")
    attr(m, "labels") <- attr(d, "Labels")
    class(m) <- "haploNet"
    m
}

.TempletonProb <- function(j, S, b = 2, r = 1)
{
    br <- b * r
    P <- numeric(max(j))
    L_jm <- function(q, j, m) {
        jm1 <- j - 1
        qonbr <- q/br
        (2*q)^jm1 * (1 - q)^(2*m + 1) * (1 - qonbr) *
            (2 - q*(br + 1)/br)^jm1 *
                (1 - 2*q*(1 - qonbr))
    }
    for (i in seq_along(P)) {
        M <- S - i
        denom <- integrate(L_jm, 0, 1, j = i, m = M)$value
        ## eq.7 from Templeton et al. 1992:
        out <- integrate(function(q) q*L_jm(q, j = i, m = M), 0, 1)$value/denom
        P[i] <- 1 - out
    }
    cumprod(P)[j]
}

haplotype <- function(x, ...) UseMethod("haplotype")

haplotype.DNAbin <- function(x, labels = NULL, ...)
{
    nms.x <- deparse(substitute(x))
    if (is.list(x)) x <- as.matrix(x)
    n <- nrow(x)
    s <- ncol(x)
    h <- .C(haplotype_DNAbin, x, n, s, integer(n), NAOK = TRUE)[[4]]
    u <- h == 0
    h[u] <- i <- which(u)
    obj <- x[i, ]
    if (is.null(labels))
        labels <- as.character(as.roman(seq_along(i)))
    rownames(obj) <- labels
    class(obj) <- c("haplotype", "DNAbin")
    attr(obj, "index") <- lapply(i, function(x) which(h == x))
    attr(obj, "from") <- nms.x
    obj
}

haplotype.character <- function(x, labels = NULL, ...)
{
    nms.x <- deparse(substitute(x))
    if (!is.matrix(x)) stop("x must be a matrix")
    h <- factor(apply(x, 1, paste, collapse = "\r"))
    h <- as.integer(h)
    obj <- x[which(!duplicated(h)), , drop = FALSE]
    N <- nrow(obj)
    if (is.null(labels))
        labels <- as.character(as.roman(1:N))
    rownames(obj) <- labels
    class(obj) <- c("haplotype", "character")
    attr(obj, "index") <- lapply(1:N, function(y) which(h == y))
    attr(obj, "from") <- nms.x
    obj
}

haplotype.numeric <- function(x, labels = NULL, ...)
    haplotype.character(x = x, labels = labels, ...)

haploFreq <- function(x, fac, split = "_", what = 2, haplo = NULL)
{
    if (missing(fac)) {
        fac <- strsplit(rownames(x), split)
        fac <- factor(sapply(fac, function(xx) xx[what]))
    } else if (length(fac) != nrow(x))
        stop("number of elements in 'fac' not the same than number of sequences")

    if (is.null(haplo)) haplo <- haplotype(x)
    h.index <- attr(haplo, "index")
    l <- nlevels(fac)
    if (l == 1) {
        res <- sapply(h.index, length)
        dim(res) <- c(length(res), 1L)
    } else {
        res <- sapply(h.index, function(xx) tabulate(fac[xx], l))
        res <- t(res)
    }
    colnames(res) <- levels(fac)
    res
}

diffHaplo <- function(h, a = 1, b = 2)
{
    s <- seg.sites(h[c(a, b), ])
    x <- h[c(a, b), s]
    cat("\n", dist.dna(x, "TS"), "transitions,",
        dist.dna(x, "TV"), "transversions\n\n")
    data.frame(pos = s, toupper(t(as.character(x))))
}

haploNet <- function(h, d = NULL, getProb = TRUE)
{
    if (!inherits(h, "haplotype"))
        stop("'h' must be of class 'haplotype'")
    freq <- sapply(attr(h, "index"), length)
    n <- length(freq) # number of haplotypes
    link <- matrix(0, 0, 3)
    if (is.null(d)) {
        ## d <- dist.dna(h, "N", pairwise.deletion = TRUE)
        d <- .C(distDNA_pegas, h, as.integer(n), ncol(h), numeric(n * (n - 1)/2), NAOK = TRUE)[[4]]
        attr(d, "Size") <- n
        attr(d, "Labels") <- rownames(h)
        attr(d, "Diag") <- attr(d, "Upper") <- FALSE
        class(d) <- "dist"
    }
    d <- as.matrix(d)
    d[col(d) >= row(d)] <- NA # put NA's in the diag and above-diag elts
    one2n <- seq_len(n)
    dimnames(d) <- list(one2n, one2n)
    step <- 1
    gr <- one2n
    altlink <- matrix(0, 0, 3) # alternative links
    while (length(unique(gr)) > 1) {
        newLinks <- which(d == step, TRUE)
        if (length(newLinks)) {
            del <- NULL
            for (i in seq_len(nrow(newLinks))) {
                h1 <- newLinks[i, 1]
                h2 <- newLinks[i, 2]
                a <- gr[h1]
                b <- gr[h2]
                ## if both nodes are already in the same
                ## subnet, then count the link as alternative
                if (a == b) {
                    del <- c(del, i)
                    if (d[h1, h2] <= step)
                        altlink <- rbind(altlink, c(newLinks[i, ], step))
                } else gr[which(gr == b)] <- a
            }
            if (!is.null(del)) newLinks <- newLinks[-del, , drop = FALSE]
            newLinks <- cbind(newLinks, rep(step, nrow(newLinks)))
            link <- rbind(link, newLinks)
        }
        step <- step + 1
    }
    Probs <- if (getProb) .TempletonProb(link[, 3], ncol(h)) else NA
    link <- cbind(link, Probs)
    dimnames(link) <- list(NULL, c("", "", "step", "Prob"))
    attr(link, "freq") <- freq
    attr(link, "labels") <- rownames(h)
    if (nrow(altlink)) {
        Probs <- if (getProb) .TempletonProb(altlink[, 3], ncol(h)) else NA
        altlink <- cbind(altlink, Probs)
        dimnames(altlink) <- dimnames(link)
        attr(link, "alter.links") <- altlink
    }
    class(link) <- "haploNet"
    link
}

print.haploNet <- function(x, ...)
{
    cat("Haplotype network with:\n")
    cat(" ", length(attr(x, "labels")), "haplotypes\n")
    N <- n <- nrow(x)
    altlinks <- attr(x, "alter.links")
    if (!is.null(altlinks)) N <- N + nrow(altlinks)
    cat(" ", N, if (N > 1) "links\n" else "link\n")
    cat("  link lengths between", x[1, 3], "and", x[n, 3], "steps\n\n")
    cat("Use print.default() to display all elements.\n")
}

.drawSymbolsHaploNet <- function(xx, yy, size, col, bg, pie)
{
    if (is.null(pie))
        symbols(xx, yy, circles = size/2, inches = FALSE,
                add = TRUE, fg = col, bg = bg)
    else {
        nc <- ncol(pie)
        co <- if (length(bg) == 1 && bg == "white") rainbow(nc) else rep(bg, length.out = nc)
        for (i in seq_along(xx))
            floating.pie.asp(xx[i], yy[i], pie[i, ], radius = size[i]/2, col = co)
    }
}

.drawAlternativeLinks <- function(xx, yy, altlink, threshold, show.mutation)
{
    s <- altlink[, 3] >= threshold[1] & altlink[, 3] <= threshold[2]
    if (!any(s)) return(NULL)
    xa0 <- xx[altlink[s, 1]]
    ya0 <- yy[altlink[s, 1]]
    xa1 <- xx[altlink[s, 2]]
    ya1 <- yy[altlink[s, 2]]
    segments(xa0, ya0, xa1, ya1, col = "grey", lty = 2)
    if (show.mutation) {
        n <- length(xx)
        .labelSegmentsHaploNet(xx, yy, altlink[s, 1:2, drop = FALSE],
                               altlink[s, 3, drop = FALSE], rep(1, n),
                               1, "black", show.mutation)
    }
}

.mutationRug <- function(x0, y0, x1, y1, n, space = 0.05, length = 0.2)
{
    ## convert inches into user-coordinates:
    l <- xinch(length)
    sp <- xinch(space)
    for (i in seq_along(x0)) {
        xstart <- seq(-(sp*(n[i] - 1))/2, by = sp, length.out = n[i])
        ystart <- rep(-l/2, n[i])
        xstop <- xstart
        ystop <- -ystart
        ## rotation if the line is not horizontal
        if (y0[i] != y1[i]) {
            theta <- atan2(y1[i] - y0[i], x1[i] - x0[i])
            tmpstart <- rect2polar(xstart, ystart)
            tmpstop <- rect2polar(xstop, ystop)
            Rstart <- tmpstart$r
            Rstop <- tmpstop$r
            THETAstart <- tmpstart$angle + theta
            THETAstop <- tmpstop$angle + theta
            xy.start <- polar2rect(Rstart, THETAstart)
            xy.stop <- polar2rect(Rstop, THETAstop)
            xstart <- xy.start$x
            ystart <- xy.start$y
            xstop <- xy.stop$x
            ystop <- xy.stop$y
        }
        ## translation:
        xm <- (x0[i] + x1[i])/2
        ym <- (y0[i] + y1[i])/2
        xstart <- xstart + xm
        ystart <- ystart + ym
        xstop <- xstop + xm
        ystop <- ystop + ym
        segments(xstart, ystart, xstop, ystop)
    }
}

.labelSegmentsHaploNet <- function(xx, yy, link, step, size, lwd, col.link, method)
{
### method: the way the segments are labelled
### 1: small segments
### 2: Klaus's method
### 3: show the number of mutations on each link

    l1 <- link[, 1]
    l2 <- link[, 2]
    switch(method, {
        .mutationRug(xx[l1], yy[l1], xx[l2], yy[l2], step)
    }, {
        ld1 <- step
        ld2 <- step# * scale.ratio
        for (i in seq_along(ld1)) {
            pc <- ((1:ld1[i]) / (ld1[i] + 1) * ld2[i] + size[l1[i]]/2) / (ld2[i] + (size[l1[i]] + size[l2[i]])/2)
            xr <- pc * (xx[l2[i]] - xx[l1[i]]) +  xx[l1[i]]
            yr <- pc * (yy[l2[i]] - yy[l1[i]]) +  yy[l1[i]]
            symbols(xr, yr, circles = rep(lwd/15, length(xr)), inches = FALSE,
                    add = TRUE, fg = col.link, bg = col.link)
        }
    }, {
        x <- (xx[l1] + xx[l2])/2
        y <- (yy[l1] + yy[l2])/2
        BOTHlabels(step, NULL, x, y, c(0.5, 0.5), "rect", NULL, NULL,
                   NULL, NULL, "black", "lightgrey", FALSE, NULL, NULL)
    })
}

replot <- function(xy = NULL, ...)
{
    on.exit(return(list(x = xx, y = yy)))

    Last.phn <- get("last_plot.haploNet", envir = .PlotHaploNetEnv)
    xx <- Last.phn$xx
    yy <- Last.phn$yy
    l1 <- Last.phn$link[, 1]
    l2 <- Last.phn$link[, 2]
    step <- Last.phn$step
    lwd <- Last.phn$lwd
    size <- Last.phn$size
    col <- Last.phn$col
    bg <- Last.phn$bg
    lty <- Last.phn$lty
    col.link <- Last.phn$col.link
    labels <- Last.phn$labels
    asp <- Last.phn$asp
    pie <- Last.phn$pie
    labels <- Last.phn$labels
    show.mutation <- Last.phn$show.mutation
    altlink <- Last.phn$alter.links
    threshold <- Last.phn$threshold

    if (is.character(labels)) {
        font <- Last.phn$font
        cex <- Last.phn$cex
    }

    getMove <- function() {
        p1 <- locator(1)
        if (is.null(p1)) return("done")
        p2 <- locator(1)
        ## find the closest node to p1...
        i <- which.min(sqrt((p1$x - xx)^2 + (p1$y - yy)^2))
        ## ... and change its coordinates:
        xx[i] <<- p2$x
        yy[i] <<- p2$y
    }

    if (is.null(xy)) {
        cat("Click close to the node you want to move, then click where to place it\n(only single left-clicks).\nRight-click to exit.\n")
        res <- getMove()
    } else {
        xx <- xy$x
        yy <- xy$y
        res <- "done"
    }

    repeat {
        plot(xx, yy, type = "n", xlab = "", ylab = "",
             axes = FALSE, bty = "n", asp = asp, ...)
        segments(xx[l1], yy[l1], xx[l2], yy[l2], lwd = lwd,
                 lty = lty, col = col.link)
        if (show.mutation)
            .labelSegmentsHaploNet(xx, yy, cbind(l1, l2), step, size, lwd,
                                   col.link, as.numeric(show.mutation))
        if (!is.null(altlink) && !identical(as.numeric(threshold), 0))
            .drawAlternativeLinks(xx, yy, altlink, threshold, show.mutation)
        .drawSymbolsHaploNet(xx, yy, size, col, bg, pie)
        if (is.character(labels))
            text(xx, yy, labels, font = font, cex = cex)
        assign("tmp", list(xx, yy), envir = .PlotHaploNetEnv)
        eval(quote(last_plot.haploNet[1:2] <- tmp), envir = .PlotHaploNetEnv)
        if (identical(res, "done")) break
        res <- getMove()
    }
}

plot.haploNet <-
    function(x, size = 1, col = "black", bg = "white",
             col.link = "black", lwd = 1, lty = 1, pie = NULL,
             labels = TRUE, font = 2, cex = 1, scale.ratio = 1,
             asp = 1, legend = FALSE, fast = FALSE, show.mutation = 1,
             threshold = c(1, 2), ...)
{
    par(xpd = TRUE)
    link <- x[, 1:2, drop = FALSE]
    l1 <- x[, 1]
    l2 <- x[, 2]
    ld <- x[, 3] * scale.ratio

    tab <- tabulate(link)
    n <- length(tab)
    show.mutation <- as.integer(show.mutation)

    ## adjust 'ld' wrt the size of the symbols:
    size <- rep(size, length.out = n)
    ld <- ld + (size[l1] + size[l2])/2

    if (n < 3) {
        xx <- c(-ld, ld)/2
        yy <- rep(0, 2)
        fast <- TRUE
    } else {
        xx <- yy <- angle <- theta <- r <- numeric(n)
        avlb <- !logical(length(ld))

        H <- vector("list", n) # the list of hierarchy of nodes...

        ## the recursive function to allocate quadrants
        foo <- function(i) {
            j <- integer() # indices of the haplotypes linked to 'i'
            for (ii in 1:2) { # look at both columns
                ll <- which(link[, ii] == i & avlb)
                if (length(ll)) {
                    newj <- link[ll, -ii]
                    r[newj] <<- ld[ll]
                    j <- c(j, newj)
                    avlb[ll] <<- FALSE
                }
            }
            if (nlink <- length(j)) {
                H[[i]] <<- j
                start <- theta[i] - angle[i]/2
                by <- angle[i]/nlink
                theta[j] <<- seq(start, by = by, length.out = nlink)
                angle[j] <<- by
                xx[j] <<- r[j] * cos(theta[j]) + xx[i]
                yy[j] <<- r[j] * sin(theta[j]) + yy[i]
                for (ii in j) foo(ii)
            }
        }

        ## start with the haplotype with the most links:
        central <- which.max(tab)
        angle[central] <- 2*pi
        foo(central)
    }

    if (!fast) {
        fCollect <- function(i) {
            ## find all nodes to move simultaneously
            ii <- H[[i]]
            if (!is.null(ii)) {
                j <<- c(j, ii)
                for (jj in ii) fCollect(jj)
            }
        }

        ## NOTE: moving this function outside of the body of plot.haploNet() is not more efficient
        energy <- function(xx, yy) {
            ## First, check line crossings
            nlink <- length(l1)
            ## rounding makes it work better (why???)
            x0 <- round(xx[l1])
            y0 <- round(yy[l1])
            x1 <- round(xx[l2])
            y1 <- round(yy[l2])
            ## compute all the slopes and intercepts:
            beta <- (y1 - y0)/(x1 - x0)
            if (any(is.na(beta))) return(Inf)
            intp <- y0 - beta*x0
            for (i in 1:(nlink - 1)) {
                for (j in (i + 1):nlink) {
                    ## in case they are parallel:
                    if (beta[i] == beta[j]) next
                    ## if both lines are vertical:
                    if (abs(beta[i]) == Inf && abs(beta[j]) == Inf) next
                    ## now find the point where both lines cross
                    if (abs(beta[i]) == Inf) { # in case the 1st line is vertical...
                        xi <- x0[i]
                        yi <- beta[j]*xi + intp[j]
                    } else if (abs(beta[j]) == Inf) { # ... or the 2nd one
                        xi <- x0[j]
                        yi <- beta[i]*xi + intp[i]
                    } else {
                        xi <- (intp[j] - intp[i])/(beta[i] - beta[j])
                        yi <- beta[i]*xi + intp[i]
                    }
                    xi <- round(xi) # rounding as above
                    yi <- round(yi)

                    if (x0[i] < x1[i]) {
                        if (xi <= x0[i] || xi >= x1[i]) next
                    } else {
                        if (xi >= x0[i] || xi <= x1[i]) next
                    }

                    if (y0[i] < y1[i]) {
                        if (yi <= y0[i] || yi >= y1[i]) next
                    } else {
                        if (yi >= y0[i] || yi <= y1[i]) next
                    }

                    ## if we reach here, the intersection point is on
                    ## the 1st segment, check if it is on the 2nd one
                    if (x0[j] < x1[j]) {
                        if (xi <= x0[j] || xi >= x1[j]) next
                    } else {
                        if (xi >= x0[j] || xi <= x1[j]) next
                    }

                    if (y0[i] < y1[j]) {
                        if (yi <= y0[j] || yi >= y1[j]) next
                    } else {
                        if (yi >= y0[j] || yi <= y1[j]) next
                    }
                    return(Inf)
                }
            }
            D <- dist(cbind(xx, yy))
            sum(1/c(D)^2, na.rm = TRUE)
        }

        Rotation <- function(rot, i, beta) {
            ## rot: indices of the nodes to rotate
            ## i: index of the node connected to 'rot' (= fixed rotation point)
            xx.rot <- xx[rot] - xx[i]
            yy.rot <- yy[rot] - yy[i]
            theta <- atan2(yy.rot, xx.rot) + beta
            h <- sqrt(xx.rot^2 + yy.rot^2)
            new.xx[rot] <<- h*cos(theta) + xx[i]
            new.yy[rot] <<- h*sin(theta) + yy[i]
        }

        OptimizeRotation <- function(node, rot) {
            ## test the direction 1st
            inc <- pi/90
            Rotation(rot, node, inc)
            new.E <- energy(new.xx, new.yy)
            if (new.E >= E) {
                inc <- -inc
                Rotation(rot, node, inc)
                new.E <- energy(new.xx, new.yy)
            }

            while (new.E < E) {
                xx <<- new.xx
                yy <<- new.yy
                E <<- new.E
                Rotation(rot, node, inc)
                new.E <- energy(new.xx, new.yy)
            }
        }

        E <- energy(xx, yy)
        new.xx <- xx
        new.yy <- yy

        nextOnes <- NULL
        for (i in H[[central]][-1]) {
            ## collect the nodes descending from 'i':
            j <- NULL # j must be initialized before calling fCollect
            fCollect(i)
            rot <- c(i, j) # index of the nodes that will be moved
            OptimizeRotation(central, rot)
            nextOnes <- c(nextOnes, i)
        }

        while (!is.null(nextOnes)) {
            newNodes <- nextOnes
            nextOnes <- NULL
            for (i in newNodes) {
                if (is.null(H[[i]])) next
                for (j in H[[i]]) {
                    fCollect(j)
                    rot <- j
                    OptimizeRotation(i, rot)
                    nextOnes <- c(nextOnes, rot)
                }
            }
        }
    }

    plot(xx, yy, type = "n", xlab = "", ylab = "",
         axes = FALSE, bty = "n", asp = asp, ...)
    segments(xx[l1], yy[l1], xx[l2], yy[l2], lwd = lwd,
             lty = lty, col = col.link)

    ## draw alternative links
    altlink <- attr(x, "alter.links")
    if (!is.null(altlink) && !identical(as.numeric(threshold), 0))
        .drawAlternativeLinks(xx, yy, altlink, threshold, show.mutation)

    if (show.mutation) {
        if (show.mutation == 1 && all(x[, 3] < 1)) {
            warning("link lengths < 1: cannot use the default for 'show.mutation' (changed to 3)")
            show.mutation <- 3
        }
        .labelSegmentsHaploNet(xx, yy, link, x[, 3], size, lwd, col.link,
                               show.mutation)
    }

    .drawSymbolsHaploNet(xx, yy, size, col, bg, pie)

    if (labels) {
        labels <- attr(x, "labels") # for export below
        text(xx, yy, labels, font = font, cex = cex)
    }
    if (legend[1]) {
        if (is.logical(legend)) {
            cat("Click where you want to draw the legend\n")
            xy <- unlist(locator(1))
        } else xy <- legend
        segments(xy[1], xy[2], xy[1] + scale.ratio, xy[2])
        text(xy[1] + scale.ratio, xy[2], " 1", adj = 0)
        if (length(unique(size)) > 1) {
            vspace <- strheight(" ")
            symbols(xy[1] + 0.5, xy[2] - 2*vspace, circles = 0.5,
                    inches = FALSE, add = TRUE)
            text(xy[1] + 0.5, xy[2] - 2*vspace, "  1", adj = 0)
        }
        if (!is.null(pie)) {
            TEXT <- paste(" ", colnames(pie))
            nc <- ncol(pie)
            co <- if (length(bg) == 1 && bg == "white") rainbow(nc) else rep(bg, length.out = nc)
            for (i in 1:nc) {
                Y <- xy[2] - 2 * (i + 1) * vspace
                symbols(xy[1] + 0.5, Y, circles = 0.5,
                        inches = FALSE, add = TRUE, bg = co[i])
                text(xy[1] + 0.5, Y, TEXT[i], adj = 0)
            }
        }
    }
    assign("last_plot.haploNet",
           list(xx = xx, yy = yy, link = link, step = x[, 3], size = size,
                col = col, bg = bg, lwd = lwd, lty = lty, col.link = col.link,
                labels = labels, font = font, cex = cex, asp = asp, pie = pie,
                show.mutation = show.mutation, alter.links = altlink,
                threshold = threshold),
           envir = .PlotHaploNetEnv)
}

plot.haplotype <- function(x, xlab = "Haplotype", ylab = "Number", ...)
{
    barplot(sapply(attr(x, "index"), length), xlab = xlab, ylab = ylab,
            names.arg = rownames(x), ...)
}

sort.haplotype <-
    function(x, decreasing = ifelse(what == "frequencies", TRUE, FALSE), what = "frequencies", ...)
{
    oc <- oldClass(x)
    from <- attr(x, "from")
    what <- match.arg(what, c("frequencies", "labels"))
    idx <- attr(x, "index")
    o <- switch(what,
                frequencies = order(sapply(idx, length), decreasing = decreasing),
                labels = order(rownames(x), decreasing = decreasing))
    x <- x[o, ]
    attr(x, "index") <- idx[o]
    class(x) <- oc
    attr(x, "from") <- from
    x
}

subset.haplotype <- function(x, minfreq = 1, maxfreq = Inf, maxna = Inf,
                             na = c("N", "?"), ...)
{
    oc <- oldClass(x)
    from <- attr(x, "from")
    idx <- attr(x, "index")
    f <- sapply(idx, length)
    s <- f <= maxfreq & f >= minfreq
    if (is.finite(maxna)) {
        na <- tolower(na)
        all <- c("r", "m", "w", "s", "k", "y", "v", "h", "d", "b", "n", "-", "?")
        if (identical(na, "all")) na <- all
        if (identical(na, "ambiguous")) na <- all[1:11]
        freq <- if (maxna < 1) FALSE else TRUE
        foo <- function(x) sum(base.freq(x, freq, TRUE)[na])
        count.na <- numeric(n <- nrow(x))
        for (i in seq_len(n)) count.na[i] <- foo(x[i, ]) # cannot use apply
        s <- s & count.na <= maxna
    }
    x <- x[s, ]
    attr(x, "index") <- idx[s]
    class(x) <- oc
    attr(x, "from") <- from
    x
}

summary.haplotype <- function(object, ...)
{
    res <- sapply(attr(object, "index"), length)
    names(res) <- rownames(object)
    res
}

print.haplotype <- function(x, ...)
{
    n <- (d <- dim(x))[1]
    DF <- summary.haplotype(x)
    cat("\nHaplotypes extracted from:", attr(x, "from"), "\n\n")
    cat("    Number of haplotypes:", n, "\n")
    cat("         Sequence length:", d[2], "\n\n")
    cat("Haplotype labels and frequencies:\n\n")
    if (n <= 40) print(DF)
    else {
        print(DF[1:40])
        cat("...\n(use summary() to print all)\n")
    }
}

"[.haplotype" <- function(x, ...)
{
    y <- NextMethod("[")
    class(y) <- "DNAbin"
    y
}

as.phylo.haploNet <- function(x, quiet = FALSE, ...)
{
    if (!is.null(attr(x, "alter.links")) && !quiet)
        warning("some links (edges) were dropped because of reticulations")

    foo <- function(xx) {
        NODES <<- c(NODES, xx)
        W <- which(mat == xx, TRUE)
        for (i in 1:nrow(W)) {
            k <- W[i, 1]
            if (DONE[k]) next
            DONE[k] <<- TRUE
            link <- mat[k, ]
            if (link[1] != xx) link <- rev(link)
            edge[ie, ] <<- link
            el[ie] <<- x[k, 3]
            ie  <<- ie + 1L
            neigh <- link[2]
            if (deg[neigh] == 1L) TIPS <<- c(TIPS, neigh) else foo(neigh)
        }
    }

    LABS <- attr(x, "labels")
    n <- length(LABS) # number of haplotype nodes
    mat <- x[, 1:2]
    deg <- tabulate(mat, n) # get the degree of each node
    deg1 <- deg == 1
    ntip <- sum(deg1)

    TIPS <- integer()
    NODES <- integer()
    edge <- mat
    edge[] <- 0L
    el <- numeric(nrow(mat))
    ie <- 1L
    DONE <- logical(nrow(mat))
    ROOT <- which(!deg1)[1]
    foo(ROOT)
    TIPNODE <- c(TIPS, NODES)
    edge <- match(edge, TIPNODE)
    dim(edge) <- dim(mat)
    phy <- list(edge = edge, edge.length = el, tip.label = LABS[TIPS],
                Nnode = length(NODES), node.label = LABS[NODES])
    class(phy) <- "phylo"
    phy
}

as.evonet.haploNet <- function(x, ...)
{
    res <- as.phylo.haploNet(x, quiet = TRUE)
    alt <- attr(x, "alter.links")
    if (!is.null(alt)) {
        LABS <- attr(x, "labels")
        o <- match(LABS, c(res$tip.label, res$node.label))
        alt <- o[alt[, 1:2]]
        dim(alt) <- c(length(alt)/2, 2)
        res$reticulation <- alt
        class(res) <- c("evonet", "phylo")
    }
    res
}

if (getRversion() >= "2.15.1") utils::globalVariables(c("network", "network.vertex.names<-"))
as.network.haploNet <- function(x, directed = FALSE, altlinks = TRUE, ...)
{
    res <- x[, 1:2]
    if (altlinks) res <- rbind(res, attr(x, "alter.links")[, 1:2])
    res <- network(res, directed = directed, ...)
    network.vertex.names(res) <- attr(x, "labels")
    res
}

if (getRversion() >= "2.15.1") utils::globalVariables("graph.edgelist")
as.igraph.haploNet <- function(x, directed = FALSE, use.labels = TRUE,
                               altlinks = TRUE, ...)
{
    y <- x[, 1:2]
    if (altlinks) y <- rbind(y, attr(x, "alter.links")[, 1:2])
    y <-
        if (use.labels) matrix(attr(x, "labels")[y], ncol = 2)
        else y - 1L
    graph.edgelist(y, directed = directed, ...)
}

haplotype.loci <- function(x, locus = 1:2, quiet = FALSE, compress = TRUE,
                           check.phase = TRUE, ...)
{
    x <- x[, attr(x, "locicol")[locus]]
    nloc <- ncol(x)
    n <- nrow(x)

    if (check.phase) {
        ## drop the rows with at least one unphased genotype:
        s <- apply(is.phased(x), 1, all)
        if (any(!s)) {
            x <- x[s, ]
            warning(paste("dropping", sum(!s), "individual(s) out of", n, "due to unphased genotype(s)"))
            n <- nrow(x)
        }
    }

### NOTE: trying to find identical rows first does not speed calculations

    ## initialise (works for all levels of ploidy)
    nh <- .checkPloidy(x[, 1, drop = FALSE]) # the number of haplotypes
    names(nh) <- NULL
    res <- matrix("", nloc, nh * n)
    class(x) <- "data.frame" # drop "loci"

### unlist(lapply()) is a bit faster than filling the matrix with 'for'
### thanks to use.names = FALSE in unlist(). The old code is:
###    y <- matrix(NA_integer_, n, nloc)
###    for (j in seq_len(nloc)) y[, j] <- as.integer(x[[j]])
    y <- unlist(lapply(x, as.character), FALSE, FALSE)
    dim(y) <- c(n, nloc)

    mysplit <- function(x)
        unlist(strsplit(x, "|", fixed = TRUE, useBytes = TRUE), FALSE, FALSE)

    k <- 1:nh
    for (i in seq_len(n)) { # loop along each individual
        if (!quiet && !(i %% 100)) cat("\rAnalysing individual no.", i, "/", n)
        tmp <- mysplit(y[i, ])
        dim(tmp) <- c(nh, nloc) # arrange the alleles in a matrix...
        res[, k] <- t(tmp)      # faster than using matrix(, byrow = TRUE)
        k <- k + nh
    }

### experimental code to run the above loop with parallel
### foo <- function(i) {
###     tmp <- mysplit(y[i, ])
###     dim(tmp) <- c(nh, nloc)
###     t(tmp)
### }
### res <- mclapply(seq_len(n), foo)
### res <- unlist(res, FALSE, FALSE)
### dim(res) <- c(nloc, n * nh)

    if (!compress) {
        rownames(res) <- colnames(x)
        return(res)
    }

    nc <- ncol(res)
    h <- .Call(unique_haplotype_loci, res, nloc, nc)
    u <- h == 0
    if (all(u)) freq <- rep(1L, nc) else {
        i <- which(u)
        res <- res[, i, drop = FALSE]
        h[u] <- i
        freq <- tabulate(h)
        freq <- freq[freq > 0]
    }

    if (!quiet) cat("\rAnalysing individual no.", n, "/", n, "\n")

    rownames(res) <- colnames(x)
    class(res) <- "haplotype.loci"
    attr(res, "freq") <- freq
    res
}

plot.haplotype.loci <- function(x, ...)
{
    y <- attr(x, "freq")
    names(y) <- apply(x, 2, paste, collapse = "-")
    barplot(y, ...)
}

dist.hamming <- function(x)
{
    n <- nrow(x)
    if (n < 2) stop("less than two haplotypes")
    d <- numeric(n * (n - 1)/2)
    k <- 1L
    for (i in 1:(n - 1)) {
        for (j in (i + 1):n) {
            d[k] <- sum(x[i, ] != x[j, ])
            k <- k + 1L
        }
    }
    attr(d, "Size") <- n
    attr(d, "Labels") <- rownames(x)
    attr(d, "Diag") <- attr(d, "Upper") <- FALSE
    attr(d, "call") <- match.call()
    attr(d, "method") <- "N"
    class(d) <- "dist"
    d
}

dist.haplotype.loci <- function(x)
{
    n <- ncol(x)
    if (n < 2) stop("less than two haplotypes")
    d <- numeric(n * (n - 1) / 2)
    k <- 1L
    for (i in 1:(n - 1)) {
        for (j in (i + 1):n) {
            d[k] <- sum(x[, i] != x[, j])
            k <- k + 1L
        }
    }
    attr(d, "Size") <- n
    attr(d, "Labels") <- colnames(x)
    attr(d, "Diag") <- attr(d, "Upper") <- FALSE
    attr(d, "call") <- match.call()
    attr(d, "method") <- "N"
    class(d) <- "dist"
    d
}

LD <- function(x, locus = c(1, 2), details = TRUE)
{
    if (length(locus) != 2)
        stop("you must specify two loci to compute linkage disequilibrium")
    hap <- haplotype.loci(x, locus = locus, quiet = TRUE)
    alleles1 <- unique(hap[1, ])
    alleles2 <- unique(hap[2, ])
    k <- length(alleles1)
    m <- length(alleles2)
    nij <- matrix(0L, k, m, dimnames = list(alleles1, alleles2))
    nij[t(hap)] <- attr(hap, "freq")
    N <- sum(nij)
    pij <- nij / N
    pi <- rowSums(pij)
    qj <- colSums(pij)
    eij <- pi %o% qj * N
    pi <- rep(pi, ncol(pij))
    qj <- rep(qj, each = nrow(pij))
    D <- pij - pi * qj
    rij <- D / sqrt(pi * (1 - pi) * qj * (1 - qj))
    df <- (k - 1) * (m - 1)
    T2 <- df * N * sum(rij^2) / (k * m)
    res <- c("T2" = T2, "df" = df, "P-val" = pchisq(T2, df, lower.tail = FALSE))
    if (details) {
        res <- list(nij, eij, rij, 2 * sum(nij * log(nij / eij)),
                    2 * sum((nij - eij)^2 / eij), res)
        names(res) <- c("Observed frequencies", "Expected frequencies", "Correlations among alleles",
                       "LRT (G-squared)", "Pearson's test (chi-squared)", "T2")
    }
    res
}

LD2 <- function(x, locus = c(1, 2), details = TRUE)
{
    if (length(locus) != 2)
        stop("you must specify two loci to compute linkage disequilibrium")
    x <- x[, attr(x, "locicol")[locus]]
    if (any(.checkPloidy(x) != 2))
        stop("linkage disequilibrium with unphased genotypes works only for diploid data")
    n <- nrow(x)
    s <- summary(x)
    alleles1 <- names(s[[1]]$allele)
    alleles2 <- names(s[[2]]$allele)
    genotypes1 <- names(s[[1]]$genotype)
    genotypes2 <- names(s[[2]]$genotype)

    ## get allele proportions:
    P <- lapply(s, "[[", "allele")
    P <- rapply(P, function(x) x/sum(x), how = "replace")

    ## compute HW disequilibrium coefficients:
    G <- lapply(s, "[[", "genotype")
    G <- rapply(G, function(x) x/sum(x), how = "replace")
    DHW <- vector("list", 2)
    for (i in 1:2) {
        alls <- if (i == 1) alleles1 else alleles2
        tmp <- G[[i]][paste(alls, alls, sep = "/")] # frequencies of all possible homozygotes
        tmp[is.na(tmp)] <- 0 # for the unobserved
        names(tmp) <- alls
        DHW[[i]] <- tmp - P[[i]]^2
    }

    Delta <- r <- numeric()

    x1 <- as.integer(x[, 1L, drop = TRUE])
    x2 <- as.integer(x[, 2L, drop = TRUE])

    for (a in alleles1) {
        Pa <- P[[1]][a]
        Da <- DHW[[1]][a]
        for (b in alleles2) {
            Pb <- P[[2]][b]
            Db <- DHW[[2]][b]

            ## find the homozygote genotypes:
            i <- match(paste(a, a, sep = "/"), genotypes1)
            j <- match(paste(b, b, sep = "/"), genotypes2)

            n1 <- n2 <- n3 <- n4 <- 1e100 # initialise
            ## if the homozygotes were not observed set the n's accordingly,
            ## else find them in the respective columns:
            if (is.na(i)) n1 <- n2 <- 0 else Ho1 <- x1 == i
            if (is.na(j)) n1 <- n3 <- 0 else Ho2 <- x2 == j

            ## count the homozygotes at both loci:
            if (as.logical(n1)) n1 <- sum(Ho1 & Ho2)

            ## the next command will find the genotypes with the 1st allele (homozygote OR heterozygote):
            tmp1 <- grep(paste0("^", a, "/|/", a, "$"), genotypes1)
            if (!is.na(i)) tmp1 <- tmp1[tmp1 != i] # remove the homozygote genotype if needed
            ## and get the heterozygotes:
            if (length(tmp1)) He1 <- x1 %in% tmp1 else n3 <- n4 <- 0

            ## ... and the same for the 2nd allele:
            tmp2 <- grep(paste0("^", b, "/|/", b, "$"), genotypes2)
            if (!is.na(j)) tmp2 <- tmp2[tmp2 != j]
            if (length(tmp2)) He2 <- x2 %in% tmp2 else n2 <- n4 <- 0

            ## if homozygotes were observed for each locus:
            if (as.logical(n2)) n2 <- sum(Ho1 & He2)
            if (as.logical(n3)) n3 <- sum(Ho2 & He1)
            if (as.logical(n4)) n4 <- sum(He1 & He2)

            delta <- (2 * n1 + n2 + n3 + 0.5 * n4)/n - 2 * Pa * Pb
            Delta <- c(Delta, delta)
            r <- c(r, delta^2 / ((Pa * (1 - Pa) + Da) * (Pb * (1 - Pb) + Db)))
        }
    }
    k <- length(alleles1)
    m <- length(alleles2)
    df <- (k - 1) * (m - 1)
    T2 <- df * n * sum(r) / (k * m)
    res <- c("T2" = T2, "df" = df, "P-val" = pchisq(T2, df, lower.tail = FALSE))
    if (details) {
        dim(Delta) <- c(k, m)
        dimnames(Delta) <- list(alleles1, alleles2)
        res <- list(Delta = Delta, T2 = res)
    }
    res
}

LDscan <- function(x, ...) UseMethod("LDscan")

LDscan.DNAbin <- function(x, quiet = FALSE, ...)
{
    if (!quiet) cat("Scanning haplotypes... ")
    ss <- seg.sites(x)
    nloci <- length(ss)
    hap <- t(as.character(x[, ss]))
    hap[hap == "n"] <- NA_character_
    if (!quiet) cat("done.\n")
    .LD <- function (hap, loc1, loc2) {
        nij <- table(hap[loc1, ], hap[loc2, ])
        if (any(dim(nij) != 2)) return(NA_real_)
        N <- sum(nij)
        pij <- nij/N
        pi <- rep(rowSums(pij), 2)
        qj <- rep(colSums(pij), each = 2)
        D <- pij - pi * qj
        rij <- D/sqrt(pi * (1 - pi) * qj * (1 - qj))
        abs(rij[1])
    }
    M <- nloci * (nloci - 1) / 2
    ldx <- numeric(M)
    k <- 0L
    for (i in 1:(nloci - 1)) {
        for (j in (i + 1):nloci) {
            k <- k + 1L
            ldx[k] <- .LD(hap, i, j)
            if (!quiet) cat("\r", round(100 * k / M), "%")
        }
    }
    if (!quiet) cat("\n")
    class(ldx) <- "dist"
    attr(ldx, "Size") <- nloci
    attr(ldx, "Labels") <- names(x)
    attr(ldx, "Diag") <- attr(ldx, "Upper") <- FALSE
    attr(ldx, "call") <- match.call()
    ldx
}

LDscan.loci <- function(x, depth = NULL, quiet = FALSE, ...)
{
    nloci <- length(attr(x, "locicol"))
    if (!quiet) cat("Scanning haplotypes... ")
    hap <- haplotype.loci(x, seq_len(nloci), TRUE, FALSE, FALSE)
    if (!quiet) cat("done.\n")
    .LD <- function (hap, loc1, loc2) {
        nij <- table(hap[loc1, ], hap[loc2, ])
        N <- sum(nij)
        pij <- nij/N
        pi <- rep(rowSums(pij), 2)
        qj <- rep(colSums(pij), each = 2)
        D <- pij - pi * qj
        rij <- D/sqrt(pi * (1 - pi) * qj * (1 - qj))
        abs(rij[1])
    }
    if (is.null(depth)) {
        M <- nloci * (nloci - 1) / 2
        ldx <- numeric(M)
        k <- 0L
        for (i in 1:(nloci - 1)) {
            for (j in (i + 1):nloci) {
                k <- k + 1L
                ldx[k] <- .LD(hap, i, j)
                if (!quiet) cat("\r", round(100 * k / M), "%")
            }
        }
        if (!quiet) cat("\n")
        class(ldx) <- "dist"
        attr(ldx, "Size") <- nloci
        attr(ldx, "Labels") <- names(x)
        attr(ldx, "Diag") <- attr(ldx, "Upper") <- FALSE
        attr(ldx, "call") <- match.call()
    } else {
        depth <- as.integer(depth)
        ldx <- vector("list", length(depth))
        k <- 0L
        for (o in depth) {
            vec <- numeric(nloci - o)
            j <- 0L
            for (i in 1:(nloci - o)) {
                j <- j + 1L
                vec[j] <- .LD(hap, i, i + o)
                if (!quiet) cat("\rScanning at depth", o, ":\t", round(100 * j / (nloci - o)), "%")
            }
            if (!quiet) cat("\n")
            k <- k + 1L
            ldx[[k]] <- vec
        }
        names(ldx) <- depth
    }
    ldx
}

LDmap <- function(d, POS = NULL, breaks = NULL, col = NULL, border = NA,
                  angle = 0, asp = 1, cex = 1, scale.legend = 0.8, ...)
{
    if (!is.null(POS) & angle != 0) {
        warning("'angle' set to 0 since 'POS' is provided")
        angle <- 0
    }
    if (is.matrix(d)) d <- as.dist(d)
    nloci <- attr(d, "Size")
    n <- nloci - 1
    m <- nloci * n /2
    nl <- if (is.null(col)) 10 else length(col) # Nb of colour levels
    if (is.null(breaks)) {
        rgd <- range(d, na.rm = TRUE)
        breaks <- seq(rgd[1], rgd[2], length.out = nl + 1)
    } else {
        nl <- length(breaks) - 1
        if (!is.null(col))
            col <- rep(col, length.out = nl)
    }
    if (is.null(col))
        col <- colorRampPalette(c("lightyellow", "red"))(nl)

    co <- col[cut(d, breaks, include.lowest = TRUE)]

    x.lab.loci <- 1:nloci - 0.75
    y.lab.loci <- 1:nloci - 0.25

    use.rect <- angle == 45

    ## assume angle = 45
    yb <- unlist(mapply(":", 1:n, n, USE.NAMES = FALSE))
    yt <- yb + 1
    xl <- unlist(mapply(rep, 0:(n - 1), n:1, USE.NAMES = FALSE))
    xr <- xl + 1

    if (!use.rect) {
        xx <- rbind(xl, xl, xr, xr, rep(NA, m))
        yy <- rbind(yb, yt, yt, yb, rep(NA, m))
        dim(xx) <- dim(yy) <- NULL
        tmp <- rect2polar(xx, yy)
        new.angle <- tmp$angle + 2 * pi * (angle - 45)/360
        xy <- polar2rect(tmp$r, new.angle)
        X <- range(xy$x, na.rm = TRUE)
        Y <- range(xy$y, na.rm = TRUE)
        ## adjust the coordinates of the labels of the loci:
        tmp <- rect2polar(x.lab.loci, y.lab.loci)
        new.angle <- tmp$angle + 2 * pi * (angle - 45)/360
        tmp <- polar2rect(tmp$r, new.angle)
        x.lab.loci <- tmp$x
        y.lab.loci <- tmp$y
    } else {
        X <- Y <- c(0, nloci)
    }

    if (!is.null(POS)) Y[1] <- -0.1 * Y[2]

    plot.default(X, Y, "n", axes = FALSE, xlab = "", ylab = "",
                 asp = asp, ...)
    if (use.rect) rect(xl, yb, xr, yt, col = co, border = border)
    else polygon(xy, col = co, border = border)

    psr <- par("usr")

    if (!is.null(POS)) {
        f <- function(x, mn, mx) {
            x <- x - mn # translate to the origin
            x <- x / mx # rescale to 1
            (psr[2] - psr[1]) * x + psr[1]
        }
        POS <- as.double(POS)
        mn <- POS[1]
        mx <- POS[nloci] - mn
        AT <- pretty(POS)
        at <- f(AT, mn, mx)
        segments(psr[1], Y[1], psr[2], Y[1], lwd = 2)
        segments(at, Y[1] - 1, at, Y[1])
        text(at, Y[1], AT, adj = c(0.5, 2), cex = cex)
        segments(f(POS, mn, mx), Y[1], x.lab.loci, y.lab.loci, lty = 3)
        mtext("Position", 1, 1.5, cex = cex)
    }

    x1 <- psr[1] + scale.legend
    y1 <- seq(psr[4] - scale.legend * nl, psr[4] - scale.legend,
              by = scale.legend)
    y2 <- y1 + scale.legend
    rect(psr[1], y1, x1, y2, col = col, border = border)
    text(x1, c(y1[1], y2), round(breaks, 3), adj = -0.1, xpd = TRUE,
         cex = cex)
    text(x.lab.loci, y.lab.loci, attr(d, "Labels"),
         srt = angle - 90, adj = 0, cex = cex)
}

all.equal.haploNet <- function(target, current, use.steps = TRUE, ...)
{
    nt1 <- sQuote(deparse(substitute(target)))
    nt2 <- sQuote(deparse(substitute(current)))
    if (identical(target, current)) return(TRUE)
    ## function to build a list to make comparisons easier
    foo <- function(x) {
        mat <- x[, 1:2]
        step <- x[, 3]
        if (!is.null(attr(x, "alter.links"))) {
            mat <- rbind(mat, attr(x, "alter.links")[, 1:2])
            step <- c(step, attr(x, "alter.links")[, 3])
        }
        dm <- dim(mat)
        mat <- attr(x, "labels")[mat]
        dim(mat) <- dm
        mat <- t(apply(mat, 1, sort))
        list(mat = mat, step = step)
    }
    ## function to arrange print of links:
    bar <- function(x) gsub("\r", "--", x)
    ## another one to print comparison of link lengths:
    bar2 <- function(x, y, z) paste0(bar(x), " (", y, ", ", z, ")")

    msg <- NULL

    labs1 <- attr(target, "labels")
    labs2 <- attr(current, "labels")
    comp12 <- is.na(match(labs1, labs2))
    comp21 <- is.na(match(labs2, labs1))
    if (all(comp12) && all(comp21))
        return(paste0("No common label between ", nt1, " and ", nt2))
    else {
        if (any(comp12))
            msg <- c(msg, paste0("Labels in ", nt2, " not in ", nt1, ":"),
                     paste(labs1[comp12], sep = ", "))
        if (any(comp21))
            msg <- c(msg, paste0("Labels in ", nt1, " not in ", nt2, ":"),
                     paste(labs2[comp21], sep = ", "))
    }

    X1 <- foo(target)
    X2 <- foo(current)
    links1 <- paste(X1$mat[, 1], X1$mat[, 2], sep = "\r")
    links2 <- paste(X2$mat[, 1], X2$mat[, 2], sep = "\r")
    comp12 <- match(links1, links2)
    if (length(links1) == length(links2)) {
        if (anyNA(comp12)) {
            comp21 <- match(links2, links1)
            msg <- c(msg, "Number of links equal",
                     paste0("Links in ", nt1, " not in ", nt2, ":"),
                     bar(links1[is.na(comp12)]),
                     paste0("Links in ", nt2, " not in ", nt1, ":"),
                     bar(links2[is.na(comp21)]))
        }
        if (use.steps) {
            tmp <- X2$step[comp12]
            test <- X1$step != tmp
            if (anyNA(test)) test[is.na(test)] <- FALSE
            if (any(test))
                msg <- c(msg, "Link lengths different (in target, in current):",
                         bar2(links1[test], X1$step[test], tmp[test]))
        }
    } else {
        msg <- c(msg, "Number of links different.")
        comp21 <- match(links2, links1)
        if (anyNA(comp12))
            msg <- c(msg, paste0("Links in ", nt1, " not in ", nt2, ":"),
                     bar(links1[is.na(comp12)]))
        if (anyNA(comp21))
            msg <- c(msg, paste0("Links in ", nt2, " not in ", nt1, ":"),
                     bar(links2[is.na(comp21)]))
        if (use.steps) {
            if (length(links1) > length(links2)) {
                tmp1 <- X1$step[!is.na(comp12)]
                tmp2 <- X2$step[comp21]
            } else {
                tmp1 <- X1$step[comp12]
                tmp2 <- X2$step[!is.na(comp21)]
            }
            test <- tmp1 != tmp2
            if (any(test))
                msg <- c(msg,
                         paste0("Links identical but of lengths different (in ", nt1, ", in ", nt2, "):"),
                         bar2(links1[tmp1][test], tmp1[test], tmp2[test]))
        }
    }
    if (is.null(msg)) TRUE else msg
}
