## options_themes_haploNet.R (2020-12-17)

##   Options and Themes to Plot "haploNet" Objects

## Copyright 2020 Emmanuel Paradis

## This file is part of the R-package `pegas'.
## See the file ../DESCRIPTION for licensing issues.

getHaploNetOptions <- function()
    get("plotHaploNetOptions", envir = .PlotHaploNetEnv)

setHaploNetOptions <- function(...)
{
    dots <- list(...)
    opts <- get("plotHaploNetOptions", envir = .PlotHaploNetEnv)
    if (length(dots)) {
        IN <- names(dots) %in% names(opts)
        if (any(IN)) {
            w <- names(dots)[IN]
            opts[w] <- dots[w]
            assign("plotHaploNetOptions", opts, envir = .PlotHaploNetEnv)
        }
        if (any(!IN))
            warning(paste("Option(s) ignored:",
                          paste(sQuote(names(dots)[!IN]), collapse = ", ")))
    }
}

#setHaploNetTheme <- function(theme)
#{
#    if (missing(theme)) {
#        cat('List of themes: "puma", "tiger", "ocean"\n')
#    } else {
#        if (! theme %in% c("puma", "tiger", "whale")) {
#            warning(paste("theme", dQuote(theme), "not found"))
#        } else {
#            L <- get(paste0(".theme.", theme))
#            par(bg = L[[1]])
#            do.call(setHaploNetOptions, L[-1])
#        }
#    }
#}
#
#.theme.puma <- list(
#    bg = "white",
#    labels = FALSE,
#    labels.cex = 1,
#    labels.font = 2,
#    link.color = "black",
#    link.type = 1,
#    link.type.alt = 2,
#    link.width = 1,
#    link.width.alt = 1,
#    haplotype.inner.color = "#CCCC4D",
#    haplotype.outer.color = "#CCCC4D",
#    mutations.cex = 1,
#    mutations.font = 1,
#    mutations.frame.background = "#0000FF4D",
#    mutations.frame.border = "black",
#    mutations.text.color = 1,
#    mutations.arrow.color = "grey",
#    mutations.arrow.type = "triangle",
#    mutations.sequence.color = "slategrey",
#    mutations.sequence.end = "round",
#    mutations.sequence.length = 0.3,
#    mutations.sequence.width = 5,
#    pie.outer.color = "black",
#    pie.inner.segments.color = NULL,
#    pie.colors.function = colorRampPalette(c("white", "#666626")),
#    scale.ratio = 1,
#    show.mutation = 3)
#
#.theme.tiger <- list(
#    bg = "yellow3",
#    labels = FALSE,
#    labels.cex = 1,
#    labels.font = 2,
#    link.color = "black",
#    link.type = 1,
#    link.type.alt = 2,
#    link.width = 1,
#    link.width.alt = 1,
#    haplotype.inner.color = "blue",
#    haplotype.outer.color = "blue",
#    mutations.cex = 1,
#    mutations.font = 1,
#    mutations.frame.background = "#0000FF4D",
#    mutations.frame.border = "black",
#    mutations.text.color = 1,
#    mutations.arrow.color = "grey",
#    mutations.arrow.type = "triangle",
#    mutations.sequence.color = "slategrey",
#    mutations.sequence.end = "round",
#    mutations.sequence.length = 0.3,
#    mutations.sequence.width = 5,
#    pie.outer.color = "black",
#    pie.inner.segments.color = NULL,
#    pie.colors.function = colorRampPalette(c("white", "orange")),
#    scale.ratio = 1,
#    show.mutation = 1)
#
#.theme.ocean <- list(
#    bg = "lightblue",
#    labels = FALSE,
#    labels.cex = 1,
#    labels.font = 2,
#    link.color = "white",
#    link.type = 1,
#    link.type.alt = 2,
#    link.width = 1,
#    link.width.alt = 1,
#    haplotype.inner.color = "navy",
#    haplotype.outer.color = "navy",
#    mutations.cex = 1,
#    mutations.font = 1,
#    mutations.frame.background = "#0000FF4D",
#    mutations.frame.border = "black",
#    mutations.text.color = 1,
#    mutations.arrow.color = "grey",
#    mutations.arrow.type = "triangle",
#    mutations.sequence.color = "slategrey",
#    mutations.sequence.end = "round",
#    mutations.sequence.length = 0.3,
#    mutations.sequence.width = 5,
#    pie.outer.color = "black",
#    pie.inner.segments.color = NULL,
#    pie.colors.function = colorRampPalette(c("white", "navy")),
#    scale.ratio = 1,
#    show.mutation = 1)
#
#
