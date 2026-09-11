library(pegas)

check_population <- function() {
    x <- data.frame(L1 = factor(c("A/A", "A/A", "B/B", "B/B")),
                    L2 = factor(c("C/C", "C/C", "D/D", "D/D")),
                    population = factor(c("one", "one", "two", "two")))
    class(x) <- c("loci", "data.frame")
    attr(x, "locicol") <- 1:2
    alternate <- factor(c("one", "two", "one", "two"))
    changed <- x
    changed$population <- alternate
    for (method in c("raw", "rarefaction", "extrapolation")) {
        stopifnot(isTRUE(all.equal(rhost(x, pop = alternate, method = method),
                                  rhost(changed, method = method))))
    }
    stopifnot(all(rhost(x, method = "raw") == 1),
              all(rhost(x, pop = alternate, method = "raw") == 0))
    # Numeric column selection, with no pre-existing population column.
    column <- changed
    names(column)[3] <- "group"
    stopifnot(isTRUE(all.equal(rhost(column, pop = 3, method = "raw"),
                              rhost(changed, method = "raw"))))
}

check_rarefaction <- function() {
    hurlbert <- getFromNamespace(".eq13.Hurlbert1971", "pegas")
    # Match the direct formula where binomial coefficients are finite.
    for (counts in list(c(3, 2, 1), c(4, 0, 2), c(6))) {
        n <- sum(counts)
        for (size in 0:n) {
            expected <- length(counts) - sum(choose(n - counts, size)) / choose(n, size)
            stopifnot(isTRUE(all.equal(hurlbert(size, counts), expected)))
        }
    }
    x <- data.frame(L1 = factor(rep(c("A/A", "B/B", "C/C", "D/D"), each = 250)),
                    population = factor(rep("one", 1000)))
    class(x) <- c("loci", "data.frame")
    attr(x, "locicol") <- 1L
    # Independent reference: probability of drawing no copies of each allele.
    expected <- sum(stats::phyper(0, m = rep(500, 4), n = rep(1500, 4),
                                 k = 1000, lower.tail = FALSE))
    result <- allelicrichness(x, method = "rarefaction", min.n = 1000)
    stopifnot(is.finite(result[1, 1]),
              isTRUE(all.equal(unname(result[1, 1]), expected)))
    curve <- rarefactionplot(x, maxn = 1000, plot = FALSE)[[1]]
    stopifnot(all(is.finite(curve)),
              isTRUE(all.equal(unname(curve[1000]), expected)))
}

check_population()
check_rarefaction()
