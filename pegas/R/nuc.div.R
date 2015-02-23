## nuc.div.R (2012-01-19)

##   Nucleotide Diversity

## Copyright 2009-2012 Emmanuel Paradis

## This file is part of the R-package `pegas'.
## See the file ../DESCRIPTION for licensing issues.

nuc.div <- function(x, variance = FALSE, pairwise.deletion = FALSE)
{
    if (pairwise.deletion && variance)
      warning("cannot compute the variance of nucleotidic diversity\nwith pairwise deletion: try 'pairwise.deletion = FALSE' instead.")
    if (is.list(x)) x <- as.matrix(x)
    n <- dim(x)[1]
    ans <- sum(dist.dna(x, "raw", pairwise.deletion = pairwise.deletion))/
        (n*(n - 1)/2)
    if (variance) {
        var <- (n + 1)*ans/(3*(n - 1)*dim(x)[2]) + 2*(n^2 + n + 3)*ans^2/(9*n*(n - 1))
        ans <- c(ans, var)
    }
    ans
}
