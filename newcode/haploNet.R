haploNet <- function(h, d = NULL)
{
    if (!inherits(h, "haplotype"))
        stop("'h' must be of class 'haplotype'")
    freq <- sapply(attr(h, "index"), length)
    n <- length(freq) # number of haplotypes
    link <- matrix(0, 0, 3)
    if (is.null(d)) d <- dist.dna(h, "N", pairwise.deletion = TRUE)
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
    link <- cbind(link, pegas:::.TempletonProb(link[, 3], ncol(h)))
    dimnames(link) <- list(NULL, c("", "", "step", "Prob"))
    attr(link, "freq") <- freq
    attr(link, "labels") <- rownames(h)
    attr(link, "alter.links") <-
        altlink <- cbind(altlink, pegas:::.TempletonProb(altlink[, 3], ncol(h)))
    class(link) <- "haploNet"
    link
}
