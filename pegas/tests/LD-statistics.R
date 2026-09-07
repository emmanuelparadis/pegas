library(pegas)

# Construct phased diploid individuals with two identical haplotypes each.
make_loci <- function(counts) {
    stopifnot(all(counts %% 2 == 0))
    cells <- which(counts > 0, arr.ind = TRUE)
    rows <- rep(seq_len(nrow(cells)), counts[cells] / 2)
    a <- rownames(counts)[cells[rows, 1]]
    b <- colnames(counts)[cells[rows, 2]]
    x <- data.frame(L1 = factor(paste(a, a, sep = "|")),
                    L2 = factor(paste(b, b, sep = "|")))
    class(x) <- c("loci", "data.frame")
    attr(x, "locicol") <- 1:2
    x
}

check_statistics <- function(counts) {
    x <- make_loci(counts)
    ans <- LD(x)
    observed <- ans[["Observed frequencies"]]
    expected <- ans[["Expected frequencies"]]
    stopifnot(identical(unname(observed),
                        unname(counts[rownames(observed), colnames(observed)])))
    pearson <- unname(stats::chisq.test(observed, correct = FALSE)$statistic)
    positive <- observed > 0
    g2 <- 2 * sum(observed[positive] * log(observed[positive] / expected[positive]))
    stopifnot(isTRUE(all.equal(ans[["Pearson's test (chi-squared)"]], pearson)),
              isTRUE(all.equal(ans[["LRT (G-squared)"]], g2)),
              identical(ans$T2, LD(x, details = FALSE)))
    invisible(ans)
}

# Empty cells are valid and must not cause NaN in G-squared.
sparse <- matrix(c(20L, 0L, 0L, 20L), 2,
                 dimnames = list(c("A", "B"), c("C", "D")))
ans <- check_statistics(sparse)
stopifnot(isTRUE(all.equal(ans[["Pearson's test (chi-squared)"]], 40)),
          isTRUE(all.equal(ans[["LRT (G-squared)"]], 80 * log(2))),
          isTRUE(all.equal(unname(ans$T2["T2"]), 40)))

# Nonzero cells, including a multiallelic example.
check_statistics(matrix(c(20L, 10L, 10L, 20L), 2,
                        dimnames = list(c("A", "B"), c("C", "D"))))
check_statistics(matrix(c(20L, 10L, 10L, 20L, 20L, 10L), 2,
                        dimnames = list(c("A", "B"), c("C", "D", "E"))))
