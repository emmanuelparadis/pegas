## zzz.R (2024-03-06)

##   Library Loading

## Copyright 2015-2024 Emmanuel Paradis

## This file is part of the R-package `pegas'.
## See the file ../COPYING for licensing issues.

.PlotHaploNetEnv <- new.env()
.cacheVCF <- new.env()

assign("plotHaploNetOptions",
list(
    labels = TRUE,
    labels.cex = 1,
    labels.font = 2,
    labels.color = "black",
    link.color = "black",
    link.color.alt = "grey",
    link.type = 1,
    link.type.alt = 2,
    link.width = 1,
    link.width.alt = 1,
    haplotype.inner.color = "white",
    haplotype.outer.color = "black",
    mutations.cex = 1,
    mutations.font = 1,
    mutations.frame.background = "#0000FF4D",
    mutations.frame.border = "black",
    mutations.text.color = 1,
    mutations.arrow.color = "black",
    mutations.arrow.type = "triangle",
    mutations.sequence.color = "#BFBFBF4D",
    mutations.sequence.end = "round",
    mutations.sequence.length = 0.3,
    mutations.sequence.width = 5,
    pie.inner.segments.color = "black",
    pie.colors.function = rainbow,
    scale.ratio = 1,
    show.mutation = 1),
envir = .PlotHaploNetEnv)
