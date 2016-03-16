/* readVCFbin.c    2016-01-19 */

/* Copyright 2015-2016 Emmanuel Paradis */

/* This file is part of the R-package `pegas'. */
/* See the file ../DESCRIPTION for licensing issues. */

#include <R.h>
#include <Rinternals.h>

/* reads a binary stream of 'SIZE' bytes from file 'FILENAME'
   after skipping 'SKIP' bytes */
SEXP read_bin_pegas(SEXP FILENAME, SEXP SIZE, SEXP SKIP)
{
    SEXP res;
    const char *filename;
    FILE *fl;
    int sz;
    double skip;
    unsigned char *p;

    PROTECT(FILENAME = coerceVector(FILENAME, STRSXP));
    PROTECT(SIZE = coerceVector(SIZE, INTSXP)); /* OK to have INT cause biggest chunk is 1e9 */
    PROTECT(SKIP = coerceVector(SKIP, REALSXP)); /* must be REAL cause file size can be > 2 Gb */
    filename = CHAR(STRING_ELT(FILENAME, 0));
    sz = INTEGER(SIZE)[0];
    skip = REAL(SKIP)[0];
    PROTECT(res = allocVector(RAWSXP, sz));
    p = RAW(res);
    fl = fopen(filename, "r");
    fseek(fl, (long)skip, SEEK_SET);
    fread(p, 1, sz, fl);
    fclose(fl);
    UNPROTECT(4);
    return res;
}

/* change bytes into an integer */
static int raw2int(unsigned char *x, int a, int b)
{
	int i, k = 1, ans = 0;

	for (i = b; i >= a; i--, k *= 10)
		ans += ((int)x[i] - 48) * k;

	return ans;
}

void extract_substring(unsigned char *x, int a, int b, char *y)
{
	int i, j;

	for (i = a, j = 0; i <= b; i++, j++) y[j] = x[i];

	y[j] = '\0';
}

/* find the locations of end-of-lines (= linefeeds; LF) in a strem of bytes 'x'
   after skipping 'SKIP' bytes; once an LF is found, 'HOP' bytes are skipped */
SEXP findEOL_C(SEXP x, SEXP SKIP, SEXP HOP)
{
    int n, i, j, *p, *buf, nEOL, hop;
    unsigned char *xr;
    SEXP res;

    PROTECT(x = coerceVector(x, RAWSXP));
    PROTECT(SKIP = coerceVector(SKIP, INTSXP));
    PROTECT(HOP = coerceVector(HOP, INTSXP));
    n = LENGTH(x);
    xr = RAW(x);
    hop = INTEGER(HOP)[0];

    buf = (int*)R_alloc(n/hop, sizeof(int));

    i = INTEGER(SKIP)[0];
    nEOL = j = 0;
    while (i < n) {
	if (xr[i] == 0x0a) {
	    nEOL++;
	    buf[j] = i + 1;
	    j++;
	    i += hop;
	}
	i++;
    }

    PROTECT(res = allocVector(INTSXP, nEOL));
    p = INTEGER(res);

    for (i = 0; i < nEOL; i++) p[i] = buf[i];

    UNPROTECT(4);
    return res;
}

/* extract the (int) field of a VCF file given a stream
   of bytes 'x' with end-of-lines stored in 'EOL' after
   skipping 'nTABtoSKIP' TABs from the start of each line
   nTABtoSKIP = 1 -> POS
   nTABtoSKIP = 5 -> QUAL */
SEXP extract_POS(SEXP x, SEXP EOL, SEXP nTABtoSKIP)
{
    int n, i, j, k, a, *p, *eol;
    unsigned char *xr;
    SEXP res;

    PROTECT(x = coerceVector(x, RAWSXP));
    PROTECT(EOL = coerceVector(EOL, INTSXP));
    PROTECT(nTABtoSKIP = coerceVector(nTABtoSKIP, INTSXP));
    xr = RAW(x);
    n = LENGTH(EOL) - 1;
    eol = INTEGER(EOL);

    PROTECT(res = allocVector(INTSXP, n));
    p = INTEGER(res);

    for (i = 0; i < n; i++) {
	j = eol[i];
	for (k = 1; k <= INTEGER(nTABtoSKIP)[0]; k++) {
	    while (xr[j] != 0x09) j++;
	    j++;
	}
	a = j;
	while (xr[j] != 0x09) j++;
	p[i] = raw2int(xr, a, j - 1);
    }

    UNPROTECT(4);
    return res;
}

/* extract the (char) field of a VCF file given a stream
   of bytes 'x' with end-of-lines stored in 'EOL' after
   skipping 'nTABtoSKIP' TABs from the start of each line
   nTABtoSKIP = 0 -> CHROM
   nTABtoSKIP = 2 -> ID
   nTABtoSKIP = 3 -> REF
   nTABtoSKIP = 4 -> ALT
   nTABtoSKIP = 6 -> FILTER
   nTABtoSKIP = 7 -> INFO
   nTABtoSKIP = 8 -> FORMAT */
SEXP extract_REF(SEXP x, SEXP EOL, SEXP nTABtoSKIP)
{
    int n, i, j, k, a, *eol;
    unsigned char *xr;
    char str[10000];
    SEXP res;

    PROTECT(x = coerceVector(x, RAWSXP));
    PROTECT(EOL = coerceVector(EOL, INTSXP));
    PROTECT(nTABtoSKIP = coerceVector(nTABtoSKIP, INTSXP));
    xr = RAW(x);
    n = LENGTH(EOL) - 1;
    eol = INTEGER(EOL);

    PROTECT(res = allocVector(STRSXP, n));

    for (i = 0; i < n; i++) {
	j = eol[i];
	for (k = 1; k <= INTEGER(nTABtoSKIP)[0]; k++) {
	    while (xr[j] != 0x09) j++;
	    j++;
	}
	a = j;
	while (xr[j] != 0x09) j++;
	/* skip if the string is longer than 10,000 bytes */
	if (j - a > 10000) continue;
	extract_substring(xr, a, j - 1, str);
	SET_STRING_ELT(res, i, mkChar(str));
    }

    UNPROTECT(4);
    return res;
}

SEXP build_factor_loci(SEXP x, SEXP N)
{
    int Nind, n, i, i1, i2, j, k, nunique, done, a, *p, *buf;
    SEXP res, geno, locnms, REF, ALT, levels;
    unsigned char *xr, RIGHT;
    char str[1000];

    PROTECT(x = coerceVector(x, RAWSXP));
    PROTECT(N = coerceVector(N, INTSXP));
    Nind = INTEGER(N)[0];
    PROTECT(geno = allocVector(INTSXP, Nind));
    p = INTEGER(geno);
    PROTECT(locnms = allocVector(STRSXP, 1));
    PROTECT(REF = allocVector(STRSXP, 1));
    PROTECT(ALT = allocVector(STRSXP, 1));
    xr = RAW(x);
    n = LENGTH(x);

    i = 0;
    while (xr[i] != 0x09) i++; /* 1st TAB */
    i++;
    while (xr[i] != 0x09) i++; /* 2nd TAB */
    a = ++i;
    while (xr[i] != 0x09) i++; /* 3rd TAB */
    extract_substring(xr, a, i - 1, str);
    SET_STRING_ELT(locnms, 0, mkChar(str));
    a = ++i;
    while (xr[i] != 0x09) i++; /* 4th TAB */
    extract_substring(xr, a, i - 1, str);
    SET_STRING_ELT(REF, 0, mkChar(str));
    a = ++i;
    while (xr[i] != 0x09) i++; /* 5th TAB */
    extract_substring(xr, a, i - 1, str);
    SET_STRING_ELT(ALT, 0, mkChar(str));
    i++;

    for (k = 1; k < 5; k++) { /* 6-9th TABs */
	while (xr[i] != 0x09) i++;
	i++;
    }

    /* test if FORMAT == GT, then genotypes are delimited by two TABs,
       otherwise by a TAB on the left and a colon on the right */
    if (xr[i - 1] == 0x09 && xr[i - 2] == 0x54 && xr[i - 3] == 0x47 && xr[i - 4] == 0x09) RIGHT = 0x09;
    else RIGHT = 0x3a;

    nunique = 1;
    buf = (int*)R_alloc(Nind, sizeof(int));
    buf[0] = i;
    p[0] = 1;
    if (Nind > 1) {
	for (j = 1; j < Nind - 1; j++) { /* start at the 2nd individual */
	    while (xr[i] != 0x09) i++;
	    i++;
	    done = 0;
	    for (k = 0; k < nunique; k++) {
		for (i1 = i, i2 = buf[k]; ; i1++, i2++) {
		    if (xr[i1] != xr[i2]) break;
		    if (xr[i1] != RIGHT && xr[i1] != 0x09) continue;
		    p[j] = k + 1;
		    done = 1;
		    break;
		}
		if (done) break;
	    }
	    if (!done) {
		buf[nunique] = i;
		p[j] = ++nunique;
	    }
	    /* DELETE THIS LINE TO FIX A BUG (2016-01-19): i = i1; */
	}

	/* CHANGE THE LINE BELOW BY THE NEXT ONE, SAME FIX THAN ABOVE (2016-01-19): */
	/* if (RIGHT == 0x3a) while (xr[i] != 0x09) i++; */
	while (xr[i] != 0x09) i++;

	/* treat the last individual separately */
	done = 0;
	i++;
	for (k = 0; k < nunique; k++) {
	    for (i1 = i, i2 = buf[k]; ; i1++, i2++) {
		if (xr[i1] != xr[i2]) break;
		if (i1 == n - 1 || xr[i1] == RIGHT) {
		    p[j] = k + 1;
		    done = 1;
		    break;
		}
	    }
	    if (done) break;
	}
	if (!done) {
	    buf[nunique] = i;
	    p[j] = ++nunique;
	}
    }

    PROTECT(levels = allocVector(STRSXP, nunique));

    for (j = 0; j < nunique; j++) {
	k = a = buf[j];
	while (xr[k + 1] != RIGHT && xr[k + 1] != 0x09 && k < n - 1) k++;
	extract_substring(xr, a, k, str);
	SET_STRING_ELT(levels, j, mkChar(str));
    }

    PROTECT(res = allocVector(VECSXP, 5));
    SET_VECTOR_ELT(res, 0, locnms);
    SET_VECTOR_ELT(res, 1, REF);
    SET_VECTOR_ELT(res, 2, ALT);
    SET_VECTOR_ELT(res, 3, geno);
    SET_VECTOR_ELT(res, 4, levels);

    UNPROTECT(8);
    return res;
}
