## plot.loci.R (2018-08-29)

##   Plot Loci Frequencies

## Copyright 2009-2018 Emmanuel Paradis

## This file is part of the R-package `pegas'.
## See the file ../DESCRIPTION for licensing issues.

plot.summary.loci <-
    function(x, loci, what = "both", layout = 1, col = c("blue", "red"), ...)
{
    what <- match.arg(what, c("both", "alleles", "genotypes"))
    layout(matrix(1:layout, ceiling(sqrt(layout))))
    if (!devAskNewPage() && !names(dev.cur()) %in% c("pdf", "postscript")) {
        devAskNewPage(TRUE)
        on.exit(devAskNewPage(FALSE))
    }
    nms <- names(x)
    N <- if (!missing(loci)) match(loci, nms) else seq_along(x)
    if (what == "both") {
        for (i in N){
            barplot(x[[i]]$allele, main = paste(nms[i], "- alleles"),
                    col = col[1], ...)
            barplot(x[[i]]$genotype, main = paste(nms[i], "- genotypes"),
                    col = col[2], ...)
        }
    } else if (what == "alleles")
        for (i in N)
            barplot(x[[i]]$allele, main = nms[i], col = col[1], ...)
    else if (what == "genotypes")
        for (i in N)
            barplot(x[[i]]$genotype, main = nms[i], col = col[2], ...)
}
