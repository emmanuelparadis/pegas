## summary.loci.R (2019-10-03)

##   Print and Summaries of Loci Objects

## Copyright 2009-2019 Emmanuel Paradis

## This file is part of the R-package `pegas'.
## See the file ../DESCRIPTION for licensing issues.

## utility to check that all individuals have the same ploidy level for each locus:
## outputs a vector with ploidy level if all individuals have the same, 0 otherwise
.checkPloidy <- function(x)
{
    PLOIDY <- getPloidy(x)
    ploidy <- integer(ncol(PLOIDY))
    for (j in seq_along(ploidy)) {
        tmp <- unique(PLOIDY[, j])
        if (length(tmp) == 1) ploidy[j] <- tmp
    }
    ploidy
}

getPloidy <- function(x)
{
    foo <- function(x, n) {
        class(x) <- NULL
        res <- integer(n)
        ploidyGENO <- nchar(gsub("[^|/]", "", levels(x))) + 1L
        uniquePloidy <- unique(ploidyGENO)
        if (length(uniquePloidy) == 1L) {
            res[] <- uniquePloidy
            return(res)
        }
        for (i in uniquePloidy) {
            k <- which(ploidyGENO == i)
            res[x %in% k] <- i
        }
        res
    }
    n <- nrow(x)
    LOCI <- attr(x, "locicol")
    p <- length(LOCI)
    ans <- matrix(NA_integer_, n, p)
    colnames(ans) <- names(x)[LOCI]
    rownames(ans) <- row.names(x)
    class(x) <- NULL
    for (j in 1:p) ans[, j] <- foo(x[[LOCI[j]]], n)
    ans
}

getAlleles <- function(x)
    lapply(x[, attr(x, "locicol"), drop = FALSE],
           function(x) unique(unlist(strsplit(levels(x), "[/|]"))))

getGenotypes <- function(x)
    lapply(x[, attr(x, "locicol"), drop = FALSE], levels)

is.phased <- function(x)
{
    foo <- function(x) {
        phased <- grep("/", levels(x), invert = TRUE)
        if (length(phased) == 0) return(logical(length(x)))
        as.integer(x) %in% phased
    }
    sapply(x[, attr(x, "locicol")], foo)
}

unphase <- function(x)
{
    locale <- Sys.getlocale("LC_COLLATE")
    if (!identical(locale, "C")) {
        Sys.setlocale("LC_COLLATE", "C")
        on.exit(Sys.setlocale("LC_COLLATE", locale))
    }
    foo <- function(x) {
        f <- function(y) sort(strsplit(y, "|", fixed = TRUE)[[1]])
        geno <- levels(x)
        bar <- grep("|", geno, fixed = TRUE)
        if (!length(bar)) return(x)
        for (i in bar)
            geno[i] <- paste(f(geno[i]), collapse = "/")
        x <- geno[x]
        factor(x)
    }
    for (i in attr(x, "locicol")) x[, i] <- foo(x[, i, drop = TRUE])
    x
}

is.snp <- function(x) UseMethod("is.snp")

is.snp.loci <- function(x) sapply(lapply(getAlleles(x), nchar), sum) == 2

print.loci <- function(x, details = FALSE, ...)
{
    if (details) print.data.frame(x) else {
        n <- dim(x)
        nloci <- length(attr(x, "locicol"))
        cat("Allelic data frame:", n[1])
        cat(" individual")
        if (n[1] > 1) cat("s")
        cat("\n                   ", nloci)
        if (nloci == 1) cat(" locus\n") else cat(" loci\n")
        nav <- n[2] - nloci
        if (nav) {
            cat("                   ", nav, "additional variable")
            if (nav > 1) cat("s")
            cat("\n")
        }
    }
}

summary.loci <- function(object, ...)
{
    LOCI <- attr(object, "locicol")
    res <- .Call("summary_loci_pegas", object, LOCI)
    names(res) <- names(object)[LOCI]
    class(res) <- "summary.loci"
    res
}

print.summary.loci <- function(x, ...)
{
    nms <- names(x)
    for (i in 1:length(x)) {
        cat("Locus ", nms[i], ":\n", sep = "")
        cat("-- Genotype frequencies:\n")
        print(x[[i]][[1]])
        cat("-- Allele frequencies:\n")
        print(x[[i]][[2]])
        cat("\n")
    }
}

"[.loci" <- function (x, i, j, drop = FALSE)
{
    ## From base-R [.data.frame
    Narg  <- nargs()  # number of arg from x,i,j that were specified
    oc    <- oldClass(x)
    colnms.old <- names(x)
    names(x) <- colnms.new <- as.character(seq_len(ncol(x)))
    loci.nms <- names(x)[attr(x, "locicol")]
    if (Narg > 2) { # More than two arguments indicates that they used matrix-like subset
       x <- NextMethod("[", drop = drop)
    } else { # Two or fewer arguments: list-like subsetting
       x <- NextMethod("[")
    }
    ## restore the class and the "locicol" attribute only if there
    ## is at least 1 col *and* at least one loci returned:
    if (inherits(x, "data.frame")) {
        locicol <- match(loci.nms, names(x))
        locicol <- locicol[!is.na(locicol)]
        if (length(locicol)) {
            attr(x, "locicol") <- locicol
            class(x) <- oc
        }
        names(x) <- colnms.old[as.numeric(names(x))]
    }
    x
}

rbind.loci <- function(...) NextMethod("rbind")
## No need to drop the class cause it calls rbind.data.frame.
## The successive bindings eventually reorders the columns to
## agree with the 1st data frame, AND its "locicol" attribute
## is kept, so no need to change anything.

cbind.loci <- function(...)
{
    x <- list(...)
    n <- length(x)
    if (n == 1) return(x[[1]])
    NC <- unlist(lapply(x, ncol))
    LOCICOL <- lapply(x, attr, "locicol")
    offset.col <- cumsum(NC)
    for (i in 2:n) LOCICOL[[i]] <- LOCICOL[[i]] + offset.col[i - 1]
    for (i in 2:n) x[[1]] <- cbind.data.frame(x[[1]], x[[i]])
    x <- x[[1]]
    ## need to restore the class, so can't use NextMethod()
    class(x) <- c("loci", "data.frame")
    attr(x, "locicol") <- unlist(LOCICOL)
    x
}

by.loci <- function(data, INDICES = data$population, FUN = NULL, ...,
                    simplify = TRUE)
{
    if (is.null(INDICES))
        stop("no 'population' column in ", sQuote(deparse(substitute(data))))
    LOCI <- attr(data, "locicol")
    p <- length(LOCI)
    if (is.null(FUN))
        FUN <- function(x) lapply(summary(x), "[[", 2)
    INDICES <- factor(INDICES)
    lv <- levels(INDICES)
    nlv <- length(lv)
    gr <- as.integer(INDICES)
    tmp <- vector("list", nlv)
    for (i in 1:nlv) tmp[[i]] <- FUN(data[gr == i, ])
    names(tmp) <- lv
    res <- vector("list", p)
    for (i in 1:p)
        res[[i]] <- do.call(rbind, lapply(tmp, "[[", i))
    names(res) <- names(data)[LOCI]
    res
}
