## F234.R (2018-11-17)

##   F-Statistics From Patterson et al. (2012)

## Copyright 2018 Emmanuel Paradis

## This file is part of the R-package `pegas'.
## See the file ../DESCRIPTION for licensing issues.

## quick function for unphased SNPs
.bySNP <- function(data) {
    LOCI <- attr(data, "locicol")
    p <- length(LOCI)
    INDICES <- factor(data$population)
    lv <- levels(INDICES)
    nlv <- length(lv)
    gr <- as.integer(INDICES)
    tmp <- vector("list", nlv)
    FUN <- function(x) {
        foo <- function(x) {
            tab <- tabulate(x, 3L)
            cbind(2 * tab[1] + tab[2], tab[2] + 2 * tab[3])
        }
        lapply(x, foo)
    }
    for (i in 1:nlv) tmp[[i]] <- FUN(data[gr == i, ])
    names(tmp) <- lv
    res <- vector("list", p)
    alls <- getAlleles(data)
    for (i in 1:p) {
        z <- do.call(rbind, lapply(tmp, "[[", i))
        dimnames(z) <- list(lv, alls[[i]])
        res[[i]] <- z
    }
    names(res) <- names(data)[LOCI]
    res
}

F2 <- function(x, allele.freq = NULL, population = NULL, check.data = TRUE,
               pops = NULL, jackknife.block.size = 10, B = 1e4)
{
    if (!is.null(pops)) {
        if (length(pops) != 2)
            stop("'pops' should give two population names")
    }

    s <- allele.freq

    if (is.null(s)) {
        if (!is.null(population)) {
            if (nrow(x) != length(population))
                stop("length of 'population' not equal to number of rows in 'x'")
        } else {
            if (is.null(x$population)) stop("no 'population' column in 'x'")
            population <- x$population
        }
        population <- factor(population)
        npop <- nlevels(population)
        if (npop < 2) stop("less than 2 populations")
        if (npop > 2) {
            if (is.null(pops))
                stop("more than 2 populations but none specified in 'pops'")
            sel <- population %in% pops
            x <- x[sel, ]
            population <- population[sel]
        }
        ## drop the other variables if any:
        p <- ncol(x)
        LOCI <- attr(x, "locicol")
        if (length(LOCI) < p) {
            x <- x[, LOCI]
            p <- ncol(x)
        }

        if (check.data) {
            snpx <- is.snp(x)
            if (!all(snpx)) {
                if (all(!snpx)) stop("no SNP loci in x")
                x <- x[, snpx]
                warning(paste(p - ncol(x), "non-SNP loci were dropped"))
                p <- ncol(x)
            }
        }
        ## append the population column:
        x$population <- population
        s <- if (any(is.phased(x))) by(x) else .bySNP(x)
    } else {
        p <- length(s)
        if (check.data) {
            nall <- unlist(lapply(s, ncol)) # number of alleles in each locus
            snpx <- nall == 2
            if (!all(snpx)) {
                if (all(!snpx)) stop("no SNP loci")
                s <- s[snpx]
                warning(paste(p - length(s), "non-SNP loci were dropped"))
                p <- length(s)
            }
        }
        npop <- nrow(s[[1]])
        if (npop < 2) stop("less than 2 populations")
        if (npop > 2) {
            if (is.null(pops))
                stop("more than 2 populations but none specified in 'pops'")
            s <- lapply(s, function(x) x[pops, ])
        }
    }

    s <- unlist(s, use.names = FALSE)
    ii <- seq(1, by = 4, length.out = p)

    na <- s[ii]
    nb <- s[ii + 1L]
    nap <- s[ii + 2L]
    nbp <- s[ii + 3L]
    sa <- na + nap
    sb <- nb + nbp
    ha <- na * nap / (sa * (sa - 1))
    hb <- nb * nbp / (sb * (sb - 1))
    a <- na/sa
    b <- nb/sb

    ## adjust if there are pops without some alleles
    if (any(tmp <- is.na(a))) a[tmp] <- 0
    if (any(tmp <- is.na(b))) b[tmp] <- 0
    if (any(tmp <- is.na(sa))) sa[tmp] <- 0
    if (any(tmp <- is.na(sb))) sb[tmp] <- 0
    if (any(tmp <- is.na(ha))) ha[tmp] <- 0
    if (any(tmp <- is.na(hb))) hb[tmp] <- 0

    F2 <- sum((a - b)^2 - ha/sa - hb/sb, na.rm = TRUE)

    if (jackknife.block.size) { # blocked jackknife
        one2p <- 1:p
        one2bs <- 1:jackknife.block.size
        Nblock <- p %/% jackknife.block.size
        rest <- p %% jackknife.block.size
        BLOCK <- lapply(0:(Nblock - 1),
                        function(i) one2p[one2bs + i * jackknife.block.size])
        if (rest) BLOCK <- c(BLOCK, list((p - rest + 1):p))
        f2_j <- sapply(BLOCK, function(i) sum((a[-i] - b[-i])^2 - ha[-i]/sa[-i] - hb[-i]/sb[-i]))
        f2jack <- sum(p * F2 - (p - 1) * f2_j)/p # jackknife estimate of F2
        varf2jack <- (Nblock - 1) * sum((f2_j - f2jack)^2)/Nblock
        Zvalue <- f2jack * sqrt(Nblock/varf2jack)
    }
    if (B) { # simple bootstrap
        fs <- numeric(B)
        for (i in 1:B) {
            is <- sample(p, size = p, replace = TRUE)
            fs[i] <- sum((a[is] - b[is])^2 - ha[is]/sa[is] - hb[is]/sb[is], na.rm = TRUE)
        }
        Pboot <- 2 * pnorm(abs(F2), sd = sd(fs), lower.tail = FALSE) # P-value
    }
    res <- c(F2 = F2)
    if (jackknife.block.size) res <- c(res, Z.value = Zvalue)
    if (B) res <- c(res, P.boot = Pboot)
    res
}

F3 <- function(x, allele.freq = NULL, population = NULL, check.data = TRUE,
               pops = NULL, jackknife.block.size = 10, B = 1e4)
{
    if (!is.null(pops)) {
        if (length(pops) != 3)
            stop("'pops' should give three population names")
    }

    s <- allele.freq

    if (is.null(s)) {
        if (!is.null(population)) {
            if (nrow(x) != length(population))
                stop("length of 'population' not equal to number of rows in 'x'")
        } else {
            if (is.null(x$population)) stop("no 'population' column in 'x'")
            population <- x$population
        }
        population <- factor(population)
        npop <- nlevels(population)
        if (npop < 3) stop("less than 3 populations")
        if (npop > 3) {
            if (is.null(pops))
                stop("more than 3 populations but none specified in 'pops'")
            sel <- population %in% pops
            x <- x[sel, ]
            population <- population[sel]
        }
        ## drop the other variables if any:
        p <- ncol(x)
        LOCI <- attr(x, "locicol")
        if (length(LOCI) < p) {
            x <- x[, LOCI]
            p <- ncol(x)
        }

        if (check.data) {
            snpx <- is.snp(x)
            if (!all(snpx)) {
                if (all(!snpx)) stop("no SNP loci in x")
                x <- x[, snpx]
                warning(paste(p - ncol(x), "non-SNP loci were dropped"))
                p <- ncol(x)
            }
        }
        ## append the population column:
        x$population <- population
        s <- if (any(is.phased(x))) by(x) else .bySNP(x)
    } else {
        p <- length(s)
        if (check.data) {
            nall <- unlist(lapply(s, ncol)) # number of alleles in each locus
            snpx <- nall == 2
            if (!all(snpx)) {
                if (all(!snpx)) stop("no SNP loci")
                s <- s[snpx]
                warning(paste(p - length(s), "non-SNP loci were dropped"))
                p <- length(s)
            }
        }
        npop <- nrow(s[[1]])
        if (npop < 3) stop("less than 3 populations")
        if (npop > 3) {
            if (is.null(pops))
                stop("more than 3 populations but none specified in 'pops'")
            s <- lapply(s, function(x) x[pops, ])
        }
    }

    s <- unlist(s, use.names = FALSE)
    ii <- seq(1, by = 6, length.out = p)

    na <- s[ii]
    nb <- s[ii + 1L]
    nc <- s[ii + 2L]
    nap <- s[ii + 3L]
    nbp <- s[ii + 4L]
    ncp <- s[ii + 5L]
    sa <- na + nap
    sb <- nb + nbp
    sc <- nc + ncp
    hc <- nc * ncp / (sc * (sc - 1))
    a <- na/sa
    b <- nb/sb
    c <- nc/sc

    ## adjust if there are pops without some alleles
    if (any(tmp <- is.na(a))) a[tmp] <- 0
    if (any(tmp <- is.na(b))) b[tmp] <- 0
    if (any(tmp <- is.na(c))) c[tmp] <- 0
    if (any(tmp <- is.na(sc))) sc[tmp] <- 0
    if (any(tmp <- is.na(hc))) hc[tmp] <- 0

    F3 <- (c - a)*(c - b) - hc/sc
    if (any(tmp <- is.na(F3))) F3[tmp] <- 0
    f3star <- sum(F3)/(2*sum(hc))

    if (jackknife.block.size) { # blocked jackknife
        one2p <- 1:p
        one2bs <- 1:jackknife.block.size
        Nblock <- p %/% jackknife.block.size
        rest <- p %% jackknife.block.size
        BLOCK <- lapply(0:(Nblock - 1),
                        function(i) one2p[one2bs + i * jackknife.block.size])
        if (rest) BLOCK <- c(BLOCK, list((p - rest + 1):p))
        f3_j <- sapply(BLOCK, function(i) sum(F3[-i])/(2*sum(hc[-i])))
        f3jack <- sum(p * f3star - (p - 1) * f3_j)/p # jackknife estimate of f3star
        varf3jack <- (Nblock - 1) * sum((f3_j - f3jack)^2)/Nblock
        Zvalue <- f3jack * sqrt(Nblock/varf3jack)
    }
    if (B) { # simple bootstrap
        fs <- numeric(B)
        for (i in 1:B) {
            is <- sample(p, size = p, replace = TRUE)
            fs[i] <- sum(F3[is])/(2*sum(hc[is]))
        }
        Pboot <- 2 * pnorm(abs(f3star), sd = sd(fs), lower.tail = FALSE) # P-value
    }
    res <- c(F3 = f3star)
    if (jackknife.block.size) res <- c(res, Z.value = Zvalue)
    if (B) res <- c(res, P.boot = Pboot)
    res
}

F4 <- function(x, allele.freq = NULL, population = NULL, check.data = TRUE,
               pops = NULL, jackknife.block.size = 10, B = 1e4)
{
    if (!is.null(pops)) {
        if (length(pops) != 4)
            stop("'pops' should give four population names")
    }

    s <- allele.freq

    if (is.null(s)) {
        if (!is.null(population)) {
            if (nrow(x) != length(population))
                stop("length of 'population' not equal to number of rows in 'x'")
        } else {
            if (is.null(x$population)) stop("no 'population' column in 'x'")
            population <- x$population
        }
        population <- factor(population)
        npop <- nlevels(population)
        if (npop < 4) stop("less than 4 populations")
        if (npop > 4) {
            if (is.null(pops))
                stop("more than 4 populations but none specified in 'pops'")
            sel <- population %in% pops
            x <- x[sel, ]
            population <- population[sel]
        }
        ## drop the other variables if any:
        p <- ncol(x)
        LOCI <- attr(x, "locicol")
        if (length(LOCI) < p) {
            x <- x[, LOCI]
            p <- ncol(x)
        }

        if (check.data) {
            snpx <- is.snp(x)
            if (!all(snpx)) {
                if (all(!snpx)) stop("no SNP loci in x")
                x <- x[, snpx]
                warning(paste(p - ncol(x), "non-SNP loci were dropped"))
                p <- ncol(x)
            }
        }
        ## append the population column:
        x$population <- population
        s <- if (any(is.phased(x))) by(x) else .bySNP(x)
    } else {
        p <- length(s)
        if (check.data) {
            nall <- unlist(lapply(s, ncol)) # number of alleles in each locus
            snpx <- nall == 2
            if (!all(snpx)) {
                if (all(!snpx)) stop("no SNP loci")
                s <- s[snpx]
                warning(paste(p - length(s), "non-SNP loci were dropped"))
                p <- length(s)
            }
        }
        npop <- nrow(s[[1]])
        if (npop < 4) stop("less than 4 populations")
        if (npop > 4) {
            if (is.null(pops))
                stop("more than 4 populations but none specified in 'pops'")
            s <- lapply(s, function(x) x[pops, ])
        }
    }

    s <- unlist(s, use.names = FALSE)
    ii <- seq(1, by = 8, length.out = p)

    na <- s[ii]
    nb <- s[ii + 1L]
    nc <- s[ii + 2L]
    nd <- s[ii + 3L]
    nap <- s[ii + 4L]
    nbp <- s[ii + 5L]
    ncp <- s[ii + 6L]
    ndp <- s[ii + 7L]
    sa <- na + nap
    sb <- nb + nbp
    sc <- nc + ncp
    sd <- nd + ndp
    ##ha <- na * nap / (sa * (sa - 1))
    ##hb <- nb * nbp / (sb * (sb - 1))
    ##hc <- nc * ncp / (sc * (sc - 1))
    ##hd <- nd * ndp / (sd * (sd - 1))
    a <- na/sa
    b <- nb/sb
    c <- nc/sc
    d <- nd/sd

    ## adjust if there are pops without some alleles
    if (any(tmp <- is.na(a))) a[tmp] <- 0
    if (any(tmp <- is.na(b))) b[tmp] <- 0
    if (any(tmp <- is.na(c))) c[tmp] <- 0
    if (any(tmp <- is.na(d))) d[tmp] <- 0

    F4 <- (a - b) * (c - d)
    sumF4 <- sum(F4)
    Den <- (a + b - 2*a*b) * (c + d - 2*c*d)
    D <- sum(F4)/sum(Den)
    if (jackknife.block.size) { # blocked jackknife
        one2p <- 1:p
        one2bs <- 1:jackknife.block.size
        Nblock <- p %/% jackknife.block.size
        rest <- p %% jackknife.block.size
        BLOCK <- lapply(0:(Nblock - 1),
                        function(i) one2p[one2bs + i * jackknife.block.size])
        if (rest) BLOCK <- c(BLOCK, list((p - rest + 1):p))
        D_j <- sapply(BLOCK, function(i) sum(F4[-i])/sum(Den[-i]))
        Dstar <- sum(p * D - (p - 1) * D_j)/p # jackknife estimate of D
        varDstar <- (Nblock - 1) * sum((D_j - Dstar)^2)/Nblock
        Zvalue <- Dstar * sqrt(Nblock/varDstar)
    }
    if (B) { # simple bootstrap
        Ds <- numeric(B)
        for (i in 1:B) {
            is <- sample(p, size = p, replace = TRUE)
            Ds[i] <- sum(F4[is])/sum(Den[is])
        }
        Pboot <- 2 * pnorm(abs(D), sd = sd(Ds), lower.tail = FALSE) # P-value
    }
    res <- c(D = D, F4 = sumF4)
    if (jackknife.block.size) res <- c(res, Z.value = Zvalue)
    if (B) res <- c(res, P.boot = Pboot)
    res
}
