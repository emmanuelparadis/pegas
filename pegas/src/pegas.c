/* pegas.c    2015-06-29 */

/* Copyright 2015 Emmanuel Paradis */

/* This file is part of the R-package `pegas'. */
/* See the file ../DESCRIPTION for licensing issues. */

#include <R.h>
#include <R_ext/Rdynload.h>
#include <Rinternals.h>

/* returns 1 if both bases are different surely, 0 otherwise */
#define DifferentBase(a, b) (a & b) < 16

void haplotype_DNAbin(unsigned char *x, int *n, int *s, int *haplo)
{
    int i1, i2, s1, s2, flag;

    i1 = 0;
    while (i1 < *n - 1) {
	if (!haplo[i1]) {
	    i2 = i1 + 1;
	    while (i2 < *n) {
		if (!haplo[i2])  {
		    s1 = i1;
		    s2 = i2;
		    flag = 1; /* initially the two sequences are considered identical */
		    while (s1 < i1 + *n * (*s - 1)) {
			if (x[s1] != x[s2]) { /* fix the fact that a and b could be gaps (-) or completely missing (?) */
			    if (DifferentBase(x[s1], x[s2])) {
				flag = 0;
				break;
			    }
			}
			s1 += *n;
			s2 += *n;
		    }
		    if (flag) haplo[i2] = i1 + 1;
		}
		i2++;
	    }
	}
	i1++;
    }
}

static R_CMethodDef C_entries[] = {
    {"haplotype_DNAbin", (DL_FUNC) &haplotype_DNAbin, 4},
    {NULL, NULL, 0}
};

SEXP read_bin_pegas(SEXP FILENAME, SEXP SIZE, SEXP SKIP);
SEXP findEOL_C(SEXP x, SEXP SKIP, SEXP HOP);
SEXP extract_POS(SEXP x, SEXP EOL, SEXP nTABtoSKIP);
SEXP extract_REF(SEXP x, SEXP EOL, SEXP nTABtoSKIP);
SEXP build_factor_loci(SEXP x, SEXP N);

static R_CallMethodDef Call_entries[] = {
    {"read_bin_pegas", (DL_FUNC) &read_bin_pegas, 3},
    {"findEOL_C", (DL_FUNC) &findEOL_C, 3},
    {"extract_POS", (DL_FUNC) &extract_POS, 3},
    {"extract_REF", (DL_FUNC) &extract_REF, 3},
    {"build_factor_loci", (DL_FUNC) &build_factor_loci, 2},
    {NULL, NULL, 0}
};

void R_init_pegas(DllInfo *info)
{
    R_registerRoutines(info, C_entries, Call_entries, NULL, NULL);
    R_useDynamicSymbols(info, FALSE);
}
