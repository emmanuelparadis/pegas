library(pegas)

m <- cbind(group = c("north", "south"), L1 = c("A:A", "A:B"),
           L2 = c("C:D", "D:D"), note = c("sample1", "sample2"))
for (selection in list(list(pop = 1, loci = 2:3),
                       list(pop = "group", loci = c("L1", "L2")))) {
    actual <- as.loci(m, allele.sep = ":", col.pop = selection$pop,
                      col.loci = selection$loci)
    expected <- as.loci(as.data.frame(m), allele.sep = ":",
                        col.pop = selection$pop, col.loci = selection$loci)
    stopifnot(identical(actual, expected),
              identical(as.character(actual$population), c("north", "south")),
              identical(as.character(actual$L1), c("A/A", "A/B")),
              identical(attr(actual, "locicol"), 2:3))
}

# Defaults and singleton dimensions must also match data-frame conversion.
for (x in list(m[, 2:3, drop = FALSE], m[1, 2:3, drop = FALSE],
               m[, 2, drop = FALSE])) {
    stopifnot(identical(as.loci(x), as.loci(as.data.frame(x))))
}

# Preserve arguments intended for the initial matrix-to-data-frame step.
actual <- as.loci(m, allele.sep = ":", col.pop = 1, col.loci = 2:3,
                  row.names = c("individual1", "individual2"))
expected <- as.loci(as.data.frame(m, row.names = c("individual1", "individual2")),
                    allele.sep = ":", col.pop = 1, col.loci = 2:3)
stopifnot(identical(actual, expected))
