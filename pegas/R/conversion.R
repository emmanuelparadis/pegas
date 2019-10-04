## conversion.R (2019-10-03)

##   Conversion Among Allelic Data Classes

## Copyright 2009-2019 Emmanuel Paradis

## This file is part of the R-package `pegas'.
## See the file ../DESCRIPTION for licensing issues.

loci2SnpMatrix <- function(x, checkSNP = TRUE)
{
    LOCI <- attr(x, "locicol")
    p <- length(LOCI)
    n <- nrow(x)
    if (checkSNP) {
        SNP <- is.snp(x)
        nonSNP <- !SNP
        NnonSNP <- sum(nonSNP)
        if (NnonSNP == p) {
            warning("no SNP found: returning NULL")
            return(NULL)
        }
        if (NnonSNP) {
            msg <- ifelse(NnonSNP == 1, "locus is not SNP: it was dropped",
                          "loci are not SNPs: they were dropped")
            warning(paste(NnonSNP, msg))
            LOCI <- LOCI[SNP]
            p <- length(LOCI)
        }
    }

    res <- matrix(as.raw(0), n, p)
    rownames(res) <- row.names(x)
    colnames(res) <- names(x)[LOCI]
    class(x) <- NULL
    for (j in 1:p) {
        y <- x[[LOCI[j]]]
        geno <- levels(y)
        ngeno <- length(geno)
        map <- integer(ngeno)
        ## take arbitrarily the first allele as the REF allele:
        REF <- charToRaw(geno[1])[1]
        for (i in 1:ngeno) {
            tmp <- charToRaw(geno[i])
            a1 <- tmp[1]
            a2 <- tmp[3]
            if (a1 != a2) {
                map[i] <- 2L
            } else {
                map[i] <- if (a1 == REF) 1L else 3L
            }
        }
        res[, j] <- as.raw(map[y])
    }
    requireNamespace("snpStats")
    new("SnpMatrix", res)
}

loci2genind <- function(x, ploidy = 2, na.alleles = c("0", "."), unphase = TRUE)
{
    ipop <- which(names(x) == "population")
    pop <- if (length(ipop)) x[[ipop]] else NULL

    if (unphase) x <- unphase(x)

    if (any(isDot <- na.alleles == ".")) na.alleles[isDot] <- "\\."
    pat <- c(paste0("^", na.alleles, "/"), paste0("/", na.alleles, "$"), paste0("/", na.alleles, "/"))
    pat <- paste(pat, collapse = "|")

    for (i in attr(x, "locicol")) {
        z <- x[[i]]
        if (length(na <- grep(pat, z))) {
            z[na] <- NA_integer_
            z <- factor(z)
        }
        x[[i]] <- z
    }

    adegenet::df2genind(as.matrix(x[, attr(x, "locicol"), drop = FALSE]),
                        sep = "/", pop = pop, ploidy = ploidy)
}

as.loci <- function(x, ...) UseMethod("as.loci")

as.loci.genind <- function(x, ...)
{
    obj <- adegenet::genind2df(x, sep = "/")
    icol <- 1:ncol(obj)
    pop <- which(names(obj) == "pop")
    if (length(pop)) {
        names(obj)[pop] <- "population"
        icol <- icol[-pop]
    }
    for (i in icol) obj[, i] <- factor(obj[, i] )
    class(obj) <- c("loci", "data.frame")
    attr(obj, "locicol") <- icol
    .check.order.alleles(obj)
}

genind2loci <- function(x) as.loci.genind(x)

## to be sure that alleles are sorted in their ASCII code
## (not in lexicographical order) whatever the locale
.sort.alleles <- function(x, index.only = FALSE)
{
    locale <- Sys.getlocale("LC_COLLATE")
    if (!identical(locale, "C")) {
        Sys.setlocale("LC_COLLATE", "C")
        on.exit(Sys.setlocale("LC_COLLATE", locale))
    }
    o <- sort.list(x) # faster than order(x)
    if (index.only) o else x[o] # faster than sort()
}

.check.order.alleles <- function(x)
{
    locale <- Sys.getlocale("LC_COLLATE")
    if (!identical(locale, "C")) {
        Sys.setlocale("LC_COLLATE", "C")
        on.exit(Sys.setlocale("LC_COLLATE", locale))
    }

    reorder.alleles <- function(x) {
        ALLOK <- TRUE
        for (i in seq_along(x)) {
            y <- x[i]
            if (!length(grep("/", y))) next # phased genotype (mixed with unphased ones)
            y <- strsplit(y, "/")[[1L]]
            z <- paste0(y[sort.list(y)], collapse = "/")
            if (!identical(y, z)) {
                ALLOK <- FALSE
                x[i] <- z
            }
        }
        if (ALLOK) return(ALLOK)
        x
    }

    LOCI <- attr(x, "locicol")
    if (is.null(LOCI)) return(x)
    oc <- oldClass(x)
    class(x) <- NULL

    for (k in LOCI) {
        y <- x[[k]]
        if (is.numeric(y)) { # haploid with alleles coded with numerics
            x[[k]] <- factor(y)
            next
        }
        lv <- levels(y) # get the genotypes of the k-th genotype
        if (!length(grep("/", lv))) next # if haploid or phased genotype
        a <- reorder.alleles(lv) # works with all levels of ploidy > 1
        if (!is.logical(a)) x[[k]] <- factor(a[as.numeric(y)])
    }
    class(x) <- oc
    x
}

as.loci.data.frame <-
    function(x, allele.sep = "/|", col.pop = NULL, col.loci = NULL, ...)
{
    if (is.null(col.pop)) {
        ipop <- which(tolower(names(x)) == "population")
        if (length(ipop)) col.pop <- ipop
    }
    if (is.character(col.pop))
        col.pop <- which(names(x) == col.pop)
    if (is.numeric(col.pop)) {
        names(x)[col.pop] <- "population"
        x[, col.pop] <- factor(x[, col.pop])
    }
    if (is.null(col.loci)) {
        col.loci <- 1:ncol(x)
        if (is.numeric(col.pop))
            col.loci <- col.loci[-col.pop]
    }
    if (is.character(col.loci))
        col.loci <- match(col.loci, names(x))
    if (allele.sep != "/|") {
        if (allele.sep == "")
            stop("alleles within a genotype must be separated")
        for (i in col.loci)
            levels(x[, i]) <- gsub(allele.sep, "/", levels(x[, i]))
    }
    class(x) <- c("loci", "data.frame")
    attr(x, "locicol") <- col.loci
    .check.order.alleles(x)
}

as.loci.factor <- function(x, allele.sep = "/|", ...)
    as.loci.data.frame(data.frame(x), allele.sep = allele.sep, ...)

as.loci.character <- function(x, allele.sep = "/|", ...)
    as.loci.data.frame(data.frame(factor(x)), allele.sep = allele.sep, ...)

alleles2loci <- function(x, ploidy = 2, rownames = NULL, population = NULL,
                         phased = FALSE)
{
    withPop <- !is.null(population)
    x <- as.data.frame(x)
    if (is.null(rownames)) {
        idx <- rownames(x)
        if (is.null(idx)) idx <- as.character(seq_len(nrow(x)))
    } else {
        idx <- as.character(x[[rownames]])
        x[[rownames]] <- NULL
        if (withPop && rownames < population)
            population <- population - 1
    }
    if (withPop) {
        pop <- x[, population]
        x <- x[, -population, drop = FALSE]
    }
    p <- ncol(x)
    if (ploidy == 1) {
        loci.nms <- colnames(x)
        obj <- vector("list", p)
        for (i in 1:p) obj[[i]] <- factor(x[, i])
    } else {
        if (p %% ploidy) stop("number of columns not a multiple of ploidy")
        nloci <- p / ploidy
        start <- seq(1, by = ploidy, length.out = nloci)
        end <- start + ploidy - 1
        loci.nms <- colnames(x)[start]
        obj <- vector("list", nloci)
        sep <- if (phased) "|" else "/"
        foo <- function(...) paste(..., sep = sep)
        for (i in seq_len(nloci))
            obj[[i]] <- factor(do.call(foo, x[, start[i]:end[i], drop = FALSE]))
    }
    names(obj) <- loci.nms
    obj <- as.data.frame(obj, row.names = idx)
    obj <- as.loci(obj)
    if (withPop) obj$population <- factor(pop)
    obj
}

loci2alleles <- function(x)
{
    ploidy <- .checkPloidy(x)
    if (any(ploidy == 0)) stop("ploidy not homogeneous within some loci")
    n <- nrow(x)
    LOCI <- attr(x, "locicol")
    x <- x[, LOCI]
    fl <- tempfile()
    on.exit(unlink(fl))
    write.loci(x, fl, allele.sep = " ", quote = FALSE,
               row.names = FALSE, col.names = FALSE)
    res <- scan(fl, what = "", quiet = TRUE)
    res <- matrix(res, n, length(res)/n, byrow = TRUE)
    rownames(res) <- rownames(x)
    colnames(res) <- paste(unlist(mapply(rep, colnames(x), each = ploidy)),
                           unlist(mapply(":", 1, ploidy)), sep = ".")
    res
}

na.omit.loci <- function(object, na.alleles = c("0", "."), ...)
{
    if (any(isDot <- na.alleles == ".")) na.alleles[isDot] <- "\\."
    pat <- c(paste0("^", na.alleles, "/"), paste0("/", na.alleles, "$"), paste0("/", na.alleles, "/"))
    pat <- paste(pat, collapse = "|")
    drop <- logical(nrow(object))
    M <- 1:ncol(object)
    for (i in attr(object, "locicol")) {
        x <- object[[i]]
        if (length(na <- grep(pat, x))) drop[na] <- TRUE
        if (any(na <- is.na(x))) drop[na] <- TRUE
    }
    object <- object[!drop, , drop = FALSE]
    for (i in M) {
        if (is.factor(x <- object[[i]])) {
            drop <- tabulate(x, nlevels(x)) == 0
            if (any(drop)) object[[i]] <- factor(x)
        }
    }
    object
}
