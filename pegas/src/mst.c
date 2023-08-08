/* mst.c    2023-08-08 */

/* Copyright 2023 Emmanuel Paradis */

/* This file is part of the R-package `pegas'. */
/* See the file ../DESCRIPTION for licensing issues. */

#include <R.h>
#include <Rinternals.h>

void _getIandJ_(int k, int n, int *i, int *j)
{
    /* algo taken from pegas 1.2 */
    double b, x, y;
    k++;

    b = n - 0.5;
    x = ceil(b - sqrt(b * b - 2 * k));
    y = n * (1 - x) + (x + 1) * x/2 + k;
    *i = (int)x - 1;
    *j = (int)y - 1;
}

typedef struct DATA {
    double x;
    int o;
} DATA;

static int comp__(const void *a, const void *b)
{
    double x, y;
    x = ((struct DATA *)a)->x;
    y = ((struct DATA *)b)->x;;
    return (x > y) - (x < y);
}

void order_(double *x, int n, int *o)
{
    int i;
    struct DATA *X;

    X = (DATA*)R_alloc(n, sizeof(DATA));

    for (i = 0; i < n; i++) {
	X[i].x = x[i];
	X[i].o = i;
    }

    qsort(X, n, sizeof(struct DATA), comp__);

    for (i = 0; i < n; i++) o[i] = X[i].o;
}

SEXP mst_C(SEXP D, SEXP size)
{
    int e, i, j, k, n, l, z, Nedge, *o, *forest, f1, f2;
    double *d, *m;
    SEXP res;

    PROTECT(D = coerceVector(D, REALSXP));
    PROTECT(size = coerceVector(size, INTSXP));
    d = REAL(D);
    l = LENGTH(D);
    n = INTEGER(size)[0];

    Nedge = n - 1;
    PROTECT(res = allocMatrix(REALSXP, Nedge, 3));
    m = REAL(res);

    forest = (int*)R_alloc(n, sizeof(int));
    for (i = 0; i < n; i++) forest[i] = i + 1;

    o = (int*)R_alloc(l, sizeof(int));
    order_(d, l, o);

    /* create the first edge */
    _getIandJ_(o[0], n, &i, &j);
    m[0] = (double) i + 1;
    m[Nedge] = (double) j + 1;
    m[2*Nedge] = d[o[0]];
    forest[j] = forest[i];

    k = 1; /* the next distance to be looked at */
    e = 1; /* the number of edges done */

    while (e < Nedge) {
	_getIandJ_(o[k], n, &i, &j);
	f1 = forest[i];
	f2 = forest[j];
	if (f2 != f1) {
            m[e] = (double) i + 1;
	    m[e + Nedge] = (double) j + 1;
	    m[e + 2*Nedge] = d[o[k]];
	    for (z = 0; z < n; z++)
		if (forest[z] == f2) forest[z] = f1;
            e++;
	}
	k++;
    }
    UNPROTECT(3);
    return res;
}
