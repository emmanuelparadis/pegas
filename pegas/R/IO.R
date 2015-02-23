## IO.R (2014-07-08)

##   Input/Ouput

## Copyright 2009-2014 Emmanuel Paradis

## This file is part of the R-package `pegas'.
## See the file ../DESCRIPTION for licensing issues.

read.loci <-
    function(file, header = TRUE, loci.sep = "", allele.sep = "/|",
             col.pop = NULL, col.loci = NULL, ...)
{
    res <- read.table(file = file, header = header, sep = loci.sep, ...)
### the lines below are in case 'row.names' is used so the col#'s
### must be, possibly, decremented by one
    ddd <- list(...)
    row.nms <- ddd$row.names
    if (!is.null(row.nms)) {
        row.nms.idx <- NULL
        if (is.character(row.nms) && length(row.nms) == 1) {
            hdr <- strsplit(scan(file, what = "", n = 1), loci.sep)[[1]]
            row.nms.idx <- which(hdr == row.nms)
        }
        if (is.numeric(row.nms)) row.nms.idx <- row.nms
        if (is.numeric(col.loci) && is.numeric(row.nms.idx))
            col.loci[col.loci > row.nms.idx] <-
                col.loci[col.loci > row.nms.idx] - 1L
        if (is.numeric(col.pop) && !is.null(row.nms.idx))
            col.pop[col.pop > row.nms.idx] <-
                col.pop[col.pop > row.nms.idx] - 1L
    }
    as.loci.data.frame(res, allele.sep = allele.sep, col.pop = col.pop,
                       col.loci = col.loci)
}

read.vcf <- function(file, nloci = 1000, skip = 0)
{
    i <- 0
    repeat {
        X <- scan(file, what = "", sep = "\n", nlines = 1, skip = i, quiet = TRUE)
        i <- i + 1
        if (identical(substr(X, 1, 6), "#CHROM")) break
    }
    labs <- strsplit(X, "\t")[[1]]
    X <- scan(file, what = "", sep = "\n", nlines = nloci, skip = i + skip, quiet = TRUE)
    nloci <- length(X) # in case X is shorter than nloci
    obj <- vector("list", nloci)
    BUFFER <- matrix("", nloci, 7)
    for (j in seq_len(nloci)) {
        x <- strsplit(X[j], "\t")[[1]]
        BUFFER[j, ] <- x[1:7]
        obj[[j]] <- factor(gsub(":[[:graph:]]*$", "", x[-(1:9)]))
    }
    names(obj) <- BUFFER[, 3]
    obj <- as.data.frame(obj)
    rownames(obj) <- labs[-(1:9)]
    obj <- as.loci(obj)#, allele.sep = "[|]")
    ## substitute the allele names:
    REF <- BUFFER[, 4]
    ALT <- BUFFER[, 5]
    for (j in seq_len(nloci)) {
        lv <- levels(obj[[j]])
        tmp <- strsplit(ALT[j], ",")[[1]]
        ## we start from last allele in case there are more than 10 alleles...
        for (i in length(tmp):1)
            lv <- gsub(as.character(i), tmp[i], lv)
        ## ... and we finish with the main allele:
        lv <- gsub("0", REF[j], lv)
        levels(obj[[j]]) <- lv
    }
    attr(obj, "CHR") <- as.integer(BUFFER[, 1])
    attr(obj, "POS") <- as.numeric(BUFFER[, 2])
    attr(obj, "QUAL") <- as.integer(BUFFER[, 6])
    attr(obj, "FILTER") <- BUFFER[, 7]
    obj
}

read.gtx <- function(file)
{
    x <- scan(file, what = "", sep = "\n", quiet = TRUE)
    last <- x[length(x)]
    nloci <- length(strsplit(substr(last, 12, nchar(last)), " +")[[1]])
    npop <- as.integer(strsplit(x[2], " +")[[1]][1])
    ## the number of individuals is found easily:
    n <- length(x) - 2L * (1L + nloci + npop)
    loci.nms <- x[seq(from = 3, by = 2, length.out = nloci)]
    j <- 2L * nloci + 3L
    k <- 1L
    pop.nms <- character(npop)
    keep <- logical(length(x))
    pop <- integer(n)
    for (i in 1:npop) {
        pop.nms[i] <- x[j]
        m <- as.integer(x[j + 1L])
        pop[k:(k + m - 1L)] <- i
        k <- k + m
        j <- j + 2L
        keep[j:(j + m - 1L)] <- TRUE
        j <- j + m
    }
    x <- gsub("^ +", "", x[keep])
    x <- matrix(unlist(strsplit(x, " +")), n, nloci + 1L, byrow = TRUE)
    obj <- as.data.frame(x[, -1])
    for (i in 1:ncol(obj)) {
        levels(obj[, i]) <-
            paste(substr(levels(obj[, i]), 1, 3),
                  substr(levels(obj[, i]), 4, 6), sep = "/")
    }
    dimnames(obj) <- list(x[, 1], loci.nms)
    pop.nms <- gsub("^ +", "", pop.nms)
    pop.nms <- gsub(" +$", "", pop.nms)
    class(pop) <- "factor"
    levels(pop) <- pop.nms
    obj$population <- pop
    attr(obj, "locicol") <- 1:nloci
    obj <- .check.order.alleles(obj)
    class(obj) <- c("loci", "data.frame")
    obj
}

write.loci <- function(x, file = "", loci.sep = " ", allele.sep = "/|", ...)
{
    if (allele.sep != "/|") {
        for (i in attr(x, "locicol"))
            levels(x[[i]]) <- gsub("[/|]", allele.sep, levels(x[[i]]))
    }
    write.table(x, file = file, sep = loci.sep, ...)
}

edit.loci <- function(name, edit.row.names = TRUE, ...)
{
    oc <- oldClass(name)
    locicol <- attr(name, "locicol")
    class(name) <- "data.frame"
    name <- NextMethod("[", edit.row.names = edit.row.names)
    class(name) <- oc
    attr(name, "locicol") <- locicol
    name
}
