/* summary_loci_pegas.c    2019-10-04 */

/* Copyright 2019 Emmanuel Paradis */

/* This file is part of the R-package `pegas'. */
/* See the file ../DESCRIPTION for licensing issues. */

#include <R.h>
#include <Rinternals.h>

#define EMPTY_STRING_LEN 100
#define MAX_NGENO 256
#define MAX_NALL 512

void tabulateAllelesGenotypes(const char **geno, int *ngeno, char **alleles, int *nall, int *taballgeno)
{
    int i = 0, j, k, l, a = 0, b = 1, newallele;
    *nall = 0;
    char buf[100];

    for (;;) {
	if (i >= *ngeno) break;
	// is it safer to start with b=0 (in case a genotype is "/...")?
	for (;;) {
	    if (geno[i][b] == '/' || geno[i][b] == '|' || geno[i][b] == '\0')
		break;
	    b++;
	}
	for (k = 0, l = a; l < b; k++, l++)
	    buf[k] = geno[i][l];
	buf[k] = '\0';
        if (*nall) {
	    newallele = 1;
	    for (j = 0; j < *nall; j++) {
		if (! strcmp(alleles[j], buf)) {
		    newallele = 0;
		    taballgeno[i + j * (*ngeno)]++;
		    break;
		}
	    }
	    if (newallele) {
		strcpy(alleles[*nall], buf);
		taballgeno[i + j * (*ngeno)]++;
		(*nall)++;
	    }
	} else {
	    strcpy(alleles[0], buf);
	    (taballgeno[0])++;
	    *nall = 1;
	}
	if (geno[i][b] == '\0') {
	    i++;
	    a = 0;
	    b = 1;
	} else {
	    a = b + 1;
	    b = a + 1;
	}
    }
}

SEXP summary_loci_pegas(SEXP x, SEXP LOCI)
{
    SEXP res, GENOTYPES, ALLELES, freqs, genotype_freqs, allele_freqs;
    const char **geno;
    char **alleles;
    int i, j, k, m, n, *loci, Nloci, *xp, ngeno, nall, *gp, *ga, a, *taballgeno, NA_count;

    SEXP TWO_NAMES;
    PROTECT(TWO_NAMES = allocVector(STRSXP, 2));
    SET_STRING_ELT(TWO_NAMES, 0, mkChar("genotype"));
    SET_STRING_ELT(TWO_NAMES, 1, mkChar("allele"));

    PROTECT(x = coerceVector(x, VECSXP));
    PROTECT(LOCI = coerceVector(LOCI, INTSXP));
    Nloci = LENGTH(LOCI);
    //    Rprintf("Nloci = %d\t", Nloci);
    loci = INTEGER(LOCI);
    n = LENGTH(VECTOR_ELT(x, 0));

    //    Rprintf("n = %d\n", n);
    PROTECT(res = allocVector(VECSXP, Nloci));

    /* fix the size of arrays passed to tabulateAllelesGenotypes() */
    geno = (const char**)R_alloc(MAX_NGENO, sizeof(const char*));
    alleles = (char**)R_alloc(MAX_NALL, sizeof(char*));
    for (k = 0; k < MAX_NALL; k++)
	alleles[k] = (char*)R_alloc(EMPTY_STRING_LEN, sizeof(char));
    /* the table genotype by allele */
    taballgeno = (int*)R_alloc(MAX_NGENO*MAX_NALL, sizeof(int));

    for (i = 0; i < Nloci; i++) {
	//	Rprintf("i = %d\t", i);
	j = loci[i] - 1;
	//	Rprintf("j = %d\t", j);
	xp = INTEGER(VECTOR_ELT(x, j));

	PROTECT(GENOTYPES = getAttrib(VECTOR_ELT(x, j), install("levels")));
	ngeno = LENGTH(GENOTYPES);

	// Rprintf("ngeno = %d\n", ngeno);

	PROTECT(genotype_freqs = allocVector(INTSXP, ngeno));
	setAttrib(genotype_freqs, R_NamesSymbol, GENOTYPES);
	gp = INTEGER(genotype_freqs);
	memset(gp, 0, ngeno * sizeof(int));
	/* count the genotypes: */
	NA_count = 0;
	for (k = 0; k < n; k++) {
	    a = xp[k];
	    (a == 0 || a == NA_INTEGER) ? NA_count++ : ++gp[a - 1]; // O's and NA's are counted together
	}
	// NOTE: for the moment, NA_count is not output

	/* copy the pointers to the genotype names */
	if (ngeno > MAX_NGENO) // resize if needed
	    geno = (const char**)R_alloc(ngeno, sizeof(const char*));
	for (k = 0; k < ngeno; k++)
	    geno[k] = CHAR(STRING_ELT(GENOTYPES, k));

	/* initialise for counting alleles */
	nall = 0;
	if (ngeno > MAX_NGENO) {
	    taballgeno = (int*)R_alloc(ngeno*MAX_NALL, sizeof(int)); // resize if needed
	    memset(taballgeno, 0, ngeno*MAX_NALL * sizeof(int)); // initialize
	} else {
	    memset(taballgeno, 0, MAX_NGENO*MAX_NALL * sizeof(int));
	}
	tabulateAllelesGenotypes(geno, &ngeno, alleles, &nall, taballgeno);

	/* copy the strings in alleles to the object ALLELES */
	PROTECT(ALLELES = allocVector(STRSXP, nall));
	for (k = 0; k < nall; k++)
	    SET_STRING_ELT(ALLELES, k, mkChar(alleles[k]));

	PROTECT(allele_freqs = allocVector(INTSXP, nall));
	ga = INTEGER(allele_freqs);
	memset(ga, 0, nall * sizeof(int));
	for (k = 0; k < nall; k++) {
	    for (m = 0; m < ngeno; m++)
		ga[k] += gp[m] * taballgeno[m + k*ngeno];
	}
	setAttrib(allele_freqs, R_NamesSymbol, ALLELES);
	PROTECT(freqs = allocVector(VECSXP, 2));
	SET_VECTOR_ELT(freqs, 0, genotype_freqs);
	SET_VECTOR_ELT(freqs, 1, allele_freqs);
	setAttrib(freqs, R_NamesSymbol, TWO_NAMES);
	SET_VECTOR_ELT(res, i, freqs);
	UNPROTECT(5);
    }

    UNPROTECT(4);
    return res;
}
