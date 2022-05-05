## IO.R (2022-05-05)

##   Input/Ouput

## Copyright 2009-2022 Emmanuel Paradis

## This file is part of the R-package `pegas'.
## See the file ../DESCRIPTION for licensing issues.

read.loci <-
    function(file, header = TRUE, loci.sep = "", allele.sep = "/|",
             col.pop = NULL, col.loci = NULL, ...)
{
    res <- read.table(file = file, header = header, sep = loci.sep,
                      stringsAsFactors = TRUE, ...)
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

read.vcf <- function(file, from = 1, to = 1e4, which.loci = NULL, quiet = FALSE)
{
    meta <- .getMETAvcf(file)
    f <- .VCFconnection(file)
    GZ <- if (inherits(f, "connection")) TRUE else FALSE

    if (is.null(which.loci)) which.loci <- from:to
    nLoci <- length(which.loci)

    labs <- strsplit(meta$LABELS, "\t")[[1]]
    nCol <- length(labs)
    n <- nCol - 9L
    hop <- 2L * nCol - 1L

    cache <- ls(envir = .cacheVCF, all.names = TRUE)
    if (! file %in% cache) {
        if (!quiet) cat("File apparently not yet accessed:\n")
        info <- VCFloci(file, what = "POS", quiet = quiet)
    }
    cache <- get(file, envir = .cacheVCF)

    FROM <- cache$FROM
    TO <- cache$TO

    nChunks <- nrow(cache)
    obj <- vector("list", nLoci)
    locnms <- character(nLoci)

    if (GZ) open(f)
    BACK <- 1L + as.integer(meta$DOS)

    ii <- 0L # number of loci read
    for (k in seq_len(nChunks)) {
        sel <- match(which.loci, FROM[k]:TO[k])
        sel <- sel[!is.na(sel)]
        ck <- cache$CHUNK.SIZES[k]
        if (GZ) Y <- readBin(f, "raw", ck)
        if (!length(sel)) next

        if (!GZ) {
            skip <- if (k == 1) 0L else sum(cache$CHUNK.SIZES[1L:(k - 1L)])
            Y <- .Call(read_bin_pegas, file, ck, skip)
        }

        skip <- if (k == 1) meta$position else 0L
        EOL <- .Call(findEOL_C, Y, skip, hop) # can multiply 'hop' by 2 if diploid

        for (i in sel) {
            start <- if (i == 1) skip + 1L else EOL[i - 1] + 1L
            end <- EOL[i] - BACK
            out <- .Call(build_factor_loci, Y[start:end], n)
            ii <- ii + 1L
            locnms[ii] <- out[[1L]]
            REF <- out[[2L]]
            ALT <- out[[3L]]
            geno <- out[[4L]]
            lv <- out[[5L]]

            ## substitute the allele names:
            tmp <- strsplit(ALT, ",")[[1L]]
            ## we start from last allele in case there are more than 10 alleles...
            for (j in length(tmp):1)
                lv <- gsub(as.character(j), tmp[j], lv)
            ## ... and we finish with the reference allele:
            lv <- gsub("0", REF, lv)

            attr(geno, "levels") <- lv
            class(geno) <- "factor"
            obj[[ii]] <- geno
            if (!quiet && !(ii %% 100)) cat("\rReading", ii, "/", nLoci, "loci")
        }
    }
    if (GZ) close(f)
    if (!quiet) cat("\rReading", ii, "/", ii, "loci.\nDone.\n")

    i2nLoci <- seq_len(ii)

    if (ii < nLoci) {
        obj[(ii + 1):nLoci] <- NULL
        locnms <- locnms[i2nLoci]
    }
    names(obj) <- locnms
    class(obj) <- c("loci", "data.frame")
    attr(obj, "locicol") <- i2nLoci
    rownames(obj) <- labs[-(1:9)]
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

write.vcf <- function(x, file, CHROM = NULL, POS = NULL, quiet = FALSE)
{
    if (!inherits(x, "loci")) stop("'x' should be of class \"loci\"")
    ALLELES <- summary(x)
    ALLELES <- lapply(ALLELES, "[[", "allele")
    LOCI <- attr(x, "locicol")
    nLOCI <- length(LOCI)
    if (is.null(CHROM)) {
        CHROM <- rep(".", nLOCI)
    } else {
        if (length(CHROM) != nLOCI)
            stop("length of 'CHROM' must be equal to the number of loci")
    }
    if (is.null(POS)) {
        POS <- rep(".", nLOCI)
    } else {
        if (length(POS) != nLOCI)
            stop("length of 'POS' must be equal to the number of loci")
    }
    NAMES <- names(x)
    if (file.exists(file)) unlink(file)
    con <- file(file, "ab")
    tmp <- paste0(c("##fileformat=VCFv4.1\n##File produced by pegas (",
                    date(), ")\n##QUAL=<missingValue=9999>\n"), collapse = "")
    writeBin(charToRaw(tmp), con)
    tmp <- paste0(c("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT",
                    row.names(x)), collapse = "\t")
    LF <- charToRaw("\n")
    writeBin(c(charToRaw(tmp), LF), con)
    class(x) <- NULL
    i <- 0L
    for (j in LOCI) {
        i <- i + 1L
        if (!quiet && j %% 100 == 0) cat("\r", j, "/", length(LOCI))
        alls <- sort.int(ALLELES[[j]], decreasing = TRUE, method = "quick")
        nalls <- length(alls)
        allele.names <- names(alls)
        withDotAllele <- any(s <- allele.names == ".")
        y <- x[[j]]

        ## if only missing data for this locus, no need to translate the genotypes
        if (!(withDotAllele && nalls == 1)) {
            if (withDotAllele) {
                alls <- alls[!s]
                allele.names <- allele.names[!s]
                nalls <- nalls - 1L
            }
            REF <- allele.names[1]
            ALT <- if (nalls == 1) "." else allele.names[-1]
            geno <- levels(y)
            ngeno <- length(geno)
            newgeno <- character(ngeno)
            o <- .C(translateGenotypesForVCF, geno, newgeno, ngeno, allele.names, nalls)
            newgeno <- o[[2]]
            y <- newgeno[as.integer(y)]
            if (nalls > 2) ALT <- paste0(ALT, collapse = ",")
        }

        tmp <- paste0(c(CHROM[i], POS[i], NAMES[j], REF, ALT,
                      "9999\t.\t.\tGT", y), collapse = "\t")
        writeBin(c(charToRaw(tmp), LF), con)
    }
    close(con)
    if (!quiet) cat("\n")
}

