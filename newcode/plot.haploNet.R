diffHaplotype <- function(h, a = 1, b = 2)
{
    s <- seg.sites(h[c(a, b), ])
    x <- h[c(a, b), s]
    cat("\n", dist.dna(x, "TS"), "transitions,",
        dist.dna(x, "TV"), "transversions\n\n")
    data.frame(pos = s, toupper(t(as.character(x))))
}

.PlotHaploNetEnv <- new.env()
source("haploNet.R")

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
    if (length(threshold) == 2) {
        min <- threshold[1]
        max <- threshold[2]
    } else {
        min <- 1
        max <- threshold
    }
    s <- altlink[, 3] >= min & altlink[, 3] <= max
    xa0 <- xx[altlink[s, 1]]
    ya0 <- yy[altlink[s, 1]]
    xa1 <- xx[altlink[s, 2]]
    ya1 <- yy[altlink[s, 2]]
    segments(xa0, ya0, xa1, ya1, col = "grey", lty = 2)
    if (show.mutation)
        .labelSegmentsHaploNet(xx, yy, altlink[s, 1:2, drop = FALSE],
                               altlink[s, 3, drop = FALSE], NULL, 1, NULL, 2)
}

.labelSegmentsHaploNet <- function(xx, yy, link, step, size, lwd, col.link, method)
{
### method: the way the segments are labelled
###   1: Klaus's method
###   2: show the number of mutations on each link

    l1 <- link[, 1]
    l2 <- link[, 2]
    switch(method, {
        ld1 <- step
        ld2 <- step# * scale.ratio
        for (i in seq_along(ld1)) {
            pc <- ((1:ld1[i]) / (ld1[i] + 1) * ld2[i] + size[l1[i]]/2) / (ld2[i] + (size[l1[i]] + size[l2[i]])/2)
            xr <- pc * (xx[l2[i]] - xx[l1[i]]) +  xx[l1[i]]
            yr <- pc * (yy[l2[i]] - yy[l1[i]]) +  yy[l1[i]]
            symbols(xr, yr, circles = rep(lwd/15, length(xr)), inches = FALSE, add = TRUE,
                    fg = col.link, bg = col.link)
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
        if (!is.null(altlink) && !identical(threshold, 0))
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
             asp = 1, legend = FALSE, fast = FALSE, show.mutation = TRUE,
             threshold = 2, ...)
{
    par(xpd = TRUE)
    link <- x[, 1:2]
    l1 <- x[, 1]
    l2 <- x[, 2]
    ld <- x[, 3] * scale.ratio

    tab <- tabulate(link)
    n <- length(tab)
    xx <- yy <- angle <- theta <- r <- numeric(n)
    avlb <- !logical(length(ld))

    ## adjust 'ld' wrt the size of the symbols:
    size <- rep(size, length.out = n)
    ld <- ld + (size[l1] + size[l2])/2

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
    if (!is.null(altlink) && !identical(threshold, 0))
        .drawAlternativeLinks(xx, yy, altlink, threshold, show.mutation)

    if (show.mutation)
        .labelSegmentsHaploNet(xx, yy, link, x[, 3], size, lwd, col.link, as.numeric(show.mutation))

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
