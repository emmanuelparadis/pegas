The function has been improved (a bit faster) and gains a new option:

\arguments{
  \item{threshold}{one or two numeric values specifying the alternative
    links to draw: if a single value is given, then all links with a
    number of mutations equal or less than this value are drawn; if two
    values are given (as a vector), then the links with a number of
    mutations equal or greater than the first one and equal or less
    than the second one are drawn (so, use twice the same value if you
    want to represent only links with this number of mutations).}
}

The alternative links are shown with dashed grey lines and the number of
mutations printed on them (see next).

The argument \code{show.mutation} has been improved:

\arguments{
  \item{show.mutation}{a logical or a numeric value: if zero (or FALSE)
    nothing is drawn on the links; if 1 (or TRUE) the mutations are
    shown with small dots on the links; if 2 the number of mutations are
    printed on the links.}
}

See the new function \code{diffHaplotype} to print the details of the
differences between two haplotypes.