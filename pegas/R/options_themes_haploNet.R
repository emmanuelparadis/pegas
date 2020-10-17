## options_themes_haploNet.R (2020-10-17)

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

setHaploNetTheme <- function(theme)
{
    if (missing(theme)) {
        cat("List of themes:\n")
        cat(.list.themes.haploNet, sep = "\n")
    } else {
        if (! theme %in% .list.themes.haploNet) {
            warning(paste("theme", theme, "not found"))
        } else {
            L <- get(paste0(".theme.", theme))
            par(bg = L[[1]])
            do.call(setHaploNetOptions, L[-1])
        }
    }
}

## themes ... maybe better to move them elsewhere?...

.list.themes.haploNet <- c("puma", "tiger")

.theme.puma <- list(
    bg = "slategrey",
    labels = FALSE,
    labels.cex = 1,
    labels.font = 2,
    link.color = "black",
    link.type = 1,
    link.type.alt = 2,
    link.width = 1,
    link.width.alt = 1,
    haplotype.inner.color = "orange",
    haplotype.outer.color = "orange",
    haplotype.shape = "circles",
    mutations.cex = 1,
    mutations.font = 1,
    mutations.frame.background = "#0000FF4D",
    mutations.frame.border = "black",
    mutations.text.color = 1,
    mutations.arrow.color = "grey",
    mutations.arrow.type = "triangle",
    mutations.sequence.color = "slategrey",
    mutations.sequence.end = "round",
    mutations.sequence.length = 0.3,
    mutations.sequence.width = 5,
    pie.outer.color = "black",
    pie.inner.segments.color = NULL,
    pie.colors.function = colorRampPalette(c("white", "orange")),
    scale.ratio = 1,
    show.mutation = 3)

.theme.tiger <- list(
    bg = "slategrey",
    labels = FALSE,
    labels.cex = 1,
    labels.font = 2,
    link.color = "black",
    link.type = 1,
    link.type.alt = 2,
    link.width = 1,
    link.width.alt = 1,
    haplotype.inner.color = "yellow3",
    haplotype.outer.color = "yellow3",
    haplotype.shape = "squares",
    mutations.cex = 1,
    mutations.font = 1,
    mutations.frame.background = "#0000FF4D",
    mutations.frame.border = "black",
    mutations.text.color = 1,
    mutations.arrow.color = "grey",
    mutations.arrow.type = "triangle",
    mutations.sequence.color = "slategrey",
    mutations.sequence.end = "round",
    mutations.sequence.length = 0.3,
    mutations.sequence.width = 5,
    pie.outer.color = "black",
    pie.inner.segments.color = NULL,
    pie.colors.function = colorRampPalette(c("white", "orange")),
    scale.ratio = 1,
    show.mutation = 1)


