## mjn.R (2020-10-09)

##   Median-Joining Network

## Copyright 2017-2020 Emmanuel Paradis

## This file is part of the R-package `pegas'.
## See the file ../DESCRIPTION for licensing issues.

mjn <- function(x, epsilon = 0, max.n.cost = 10000, prefix = "median.vector_", quiet = FALSE)
{
    if (is.data.frame(x)) x <- as.matrix(x)
    if (mode(x) == "numeric") {
        if (!all(unique.default(x) %in% 0:1))
            stop("mjn() requires binary 0/1 data")
    } else {
        if (!inherits(x, "DNAbin"))
            stop("mjn() requires DNA or binary 0/1 data")
        if (is.list(x)) x <- as.matrix(x)
    }

    getIndex <- function(i, j, n) {
        if (i < j) n*(i - 1) - i*(i - 1)/2 + j - i
        else n*(j - 1) - j*(j - 1)/2 + i - j
    }

    ## the reverse operation:
    getIandJ <- function(ij, n) {
        j <- 1L
        N <- n - 1L
        while (ij > N) {
            j <- j + 1L
            N <- N + n - j
        }
        i <- n - (N - ij)
        c(j, i)
    }

    purgeObsolete <- function()
    {
        res <- integer()
        if (nrow(x) > n) {
            tab <- tabulate(c(link1, link2))[(n + 1):nrow(x)]
            sel <- tab == 1 | tab == 2
            if (any(sel)) {
                res <- which(sel) + n
                x <<- x[-res, ]
            }
        }
        res
    }

    if (inherits(x, "DNAbin")) {
        getMedianVectors <- function(x) {
            ## returns 1 or 3 sequences
            p <- as.integer(ncol(x))
            #browser()
            ans <- .C(getMedianVectors_DNAbin_mjn, x, p, raw(3L * p), 1L)
            res <- ans[[3L]]
            dim(res) <- c(3L, p)
            if (ans[[4L]]) res <- res[1L, , drop = FALSE]
            res
        }
    } else {
        getMedianVectors <- function(x) { # for binary (0/1) data
            p <- ncol(x)
            res <- matrix(NA_real_, 1L, p)
            for (j in 1:p) res[, j] <- as.integer(sum(x[, j]) > 1)
            res
        }
    }

    funDist <-
        switch(mode(x),
               "raw" = function(x) dist.dna(x, "n", pairwise.deletion = TRUE),
               "numeric" = function(x) dist(x, "manhattan"))

    funVariableCol <-
        switch(mode(x),
               "raw" = seg.sites,
               "numeric" = function(x) {
            foo <- function(x) length(unique.default(x)) > 1
            which(apply(x, 2, foo))
        })

    alreadyIn <-
        switch(mode(x),
               "raw" = function(x, table) .Call(alreadyIn_mjn_DNAbin, x, table),
               "numeric" = function(x, table) any(apply(table, 1, function(y) all(y == x))))

    if (mode(x) == "raw") class(x) <- NULL

    ## initialization --
    n <- nrow(x)
    nextX <- n + 1L

    ## step 1 --
    Delta <- funDist(x)
    ## find the d-step components with an MST:
    MST <- mst(Delta)

    if (n < 3) {
        warning("less than 3 observations: function mst() was used")
        return(MST)
    }

    ## step 2 --
    link1 <- MST[, 1]
    link2 <- MST[, 2]
    weight <- MST[, 3]

    ## add extra links depending on the value of epsilon
    threshold <- max(weight) + epsilon
    minDelta <- min(Delta)
    extra <- which(Delta > minDelta & Delta <= minDelta + threshold)
    i <- 1L
    while (i <= length(extra)) {
        p <- getIandJ(extra[i], n = n)
        if (!any(p[1] == link1 & p[2] == link2)) {
            link1 <- c(link1, p[1])
            link2 <- c(link2, p[2])
            weight <- c(weight, Delta[extra[i]])
        }
        i <- i + 1L
    }

    n.cost <- 0L

    lambda <- Inf
    if (!quiet) cat("Adding median vectors:")
    repeat {
        ## step 3 --
        purged <- purgeObsolete()
        if (length(purged)) {
            sel <- purged == link1 | purged == link2
            link1 <- link1[!sel]
            link2 <- link2[!sel]
            weight <- weight[!sel]
        }
        triplet <- combn(1:nrow(x), 3)
        median.vector <- NULL
        for (i in 1:ncol(triplet)) {
            s <- triplet[, i]
            ## step 4 --
            U <- s[1]; V <- s[2]; W <- s[3]
            l1 <- any(U == link1 & V == link2)
            l2 <- any(U == link1 & W == link2)
            l3 <- any(V == link1 & W == link2)
            ## check if at least 2 links feasible:
            if (l1 + l2 + l3 < 2) next
            ## get the median vectors:
            subsamp <- x[s, ] # get the sequences for triplet 's'
            s4 <- funVariableCol(subsamp)
            if (!length(s4))
                stop("maybe there are duplicated sequences in your data")
            subsamp4 <- subsamp[, s4, drop = FALSE] # get only the variable cols for 's'
            X <- getMedianVectors(subsamp4)
            tmp <- subsamp[rep(1, nrow(X)), , drop = FALSE]
            tmp[, s4] <- X

            ## check that 'tmp' is not yet among the current sequence types:
            check.presence <- alreadyIn(tmp, x)
            if (all(check.presence)) next
            if (any(check.presence)) tmp <- tmp[!check.presence, , drop = FALSE]
            ## ... same thing with the sequences in median.vector:
            if (!is.null(median.vector)) {
                check.presence <- alreadyIn(tmp, median.vector)
                if (all(check.presence)) next
                if (any(check.presence)) tmp <- tmp[!check.presence, , drop = FALSE]
            }

            ## code for 0/1 data:
            ## if (alreadyIn(tmp, x)) next
            ## if (!is.null(median.vector)) if (alreadyIn(tmp, median.vector)) next

            m <- nrow(tmp)
            rownames(tmp) <- paste0(prefix, nextX:(nextX + m - 1))
            ## rownames(tmp) <- paste("median.vector", nextX, sep = "_")
            ## dist2X <- funDist(rbind(tmp, subsamp))[1:3]
            dist2X <- matrix(0, m, 3)
            for (k in 1:m)
                dist2X[k, ] <- funDist(rbind(tmp[k, ], subsamp))[1:3]

            connection.cost <- rowSums(dist2X)
            ## connection.cost <- sum(dist2X)

            n.cost <- n.cost + length(connection.cost)
            if (n.cost > max.n.cost) stop("too many connection costs computed")

            ## if (connection.cost < lambda) lambda <- connection.cost
            ## if (connection.cost <= lambda + epsilon) {
            ##     link1 <- c(link1, s)
            ##     link2 <- c(link2, rep(nextX, 3))
            ##     weight <- c(weight, dist2X)
            ##     median.vector <- if (is.null(median.vector)) tmp else rbind(median.vector, tmp)
            ##    nextX <- nextX + 1L
            ## }
            if (any(connection.cost < lambda)) lambda <- min(connection.cost)
            for (k in 1:m) {
                if (connection.cost[k] <= lambda + epsilon) {
                    link1 <- c(link1, s)
                    link2 <- c(link2, rep(nextX, 3))
                    weight <- c(weight, dist2X[k, ])
                    median.vector <-
                        if (is.null(median.vector)) tmp[k, , drop = FALSE]
                        else rbind(median.vector, tmp[k, , drop = FALSE])
                    nextX <- nextX + 1L
                }
            }
        }
        if (!quiet) cat("", nrow(median.vector))
        if (is.null(median.vector)) break else x <- rbind(x, median.vector)
    }

    ## step 5 --
    N <- nrow(x)
    d <- funDist(x)
    dnet <- numeric(N*(N - 1)/2) # distances inferred from the network
    m <- matrix(NA_real_, 0, 3)
    forest <- 1:N
    ud <- sort(unique(d)) # unique distances in increasing order

    ## build the network:
    i <- 1L
    if (!quiet) {
        cat("\nBuilding the network... number of clusters:")
        previousNC <- N
    }
    while (length(unique(forest)) > 1) {
        if (!quiet) {
            if (length(unique(forest)) < previousNC) {
                previousNC <- length(unique(forest))
                cat("", previousNC)
            }
        }
        delta <- ud[i]
        for (iud in which(d == delta)) {
            p <- getIandJ(iud, N)
            p1 <- p[1L]
            p2 <- p[2L]
            f1 <- forest[p1]
            f2 <- forest[p2]
            if (f1 != f2) {
                ## update dnet:
                for (k1 in which(forest == f1)) {
                    for (k2 in which(forest == f2)) {
                        A <- if (k1 == p1) 0 else dnet[getIndex(k1, p1, N)]
                        B <- if (k2 == p2) 0 else dnet[getIndex(k2, p2, N)]
                        dnet[getIndex(k1, k2, N)] <- A + B + delta
                    }
                }
                forest[forest == f2] <- f1
                m <- rbind(m, c(p, delta))
            } else {
                if (delta < dnet[iud]) {
                    m <- rbind(m, c(p, delta))
                    ## need to update dnet
                    n1 <- c(m[m[, 2] == p1, 1], m[m[, 1] == p1, 2])
                    n2 <- c(m[m[, 2] == p2, 1], m[m[, 1] == p2, 2])
                    for (k1 in c(p1, n1)) {
                        for (k2 in c(p2, n2)) {
                            if (k1 == k2) next
                            j <- getIndex(k1, k2, N)
                            A <- if (k1 == p1) 0 else dnet[getIndex(k1, p1, N)]
                            B <- if (k2 == p2) 0 else dnet[getIndex(k2, p2, N)]
                            tmp <- A + B + delta
                            if (dnet[j] > tmp) dnet[j] <- tmp
                        }
                    }
                }
            }
        }
        i <- i + 1L
        if (i == 2e3) stop("too many links")
    }
    ## define the alternative links wrt the MST:
    MST <- mst(funDist(x))
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
    dimnames(m) <- list(NULL, c("", "", "step"))
    attr(m, "labels") <- rownames(x)
    if (mode(x) == "raw") class(x) <- "DNAbin"
    attr(m, "data") <- x
    attr(m, "prefix") <- prefix
    class(m) <- c("mjn", "haploNet")
    if (!quiet) cat(" 1\n")
    m
}

plot.mjn <- function(x, shape = c("circles", "diamonds"),
                     bg = c("green", "slategrey"), labels = FALSE, ...)
{
    prefix <- attr(x, "prefix")
    labs <- labels(x)
    ## identify the median vectors among the labels (1 or 2):
    mv <- grepl(paste0("^", prefix), labs) + 1L
    plot.haploNet(x, bg = bg[mv], shape = shape[mv], labels = labels, ...)
}
