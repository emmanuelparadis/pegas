/* pegas.c    2021-04-08 */

/* Copyright 2015-2021 Emmanuel Paradis */

/* This file is part of the R-package `pegas'. */
/* See the file ../DESCRIPTION for licensing issues. */

#include <R.h>
#include <R_ext/Rdynload.h>
#include <Rinternals.h>

/* returns 8 if the base is known surely, 0 otherwise */
#define KnownBase(a) (a & 8)

/* returns 1 if both bases are different surely, 0 otherwise */
#define DifferentBase(a, b) (a & b) < 16

int identical_seqs(unsigned char *x, int i, int j, int n, int s)
{
    int k = i + n * (s - 1);
    for (;;) {
	if (x[i] != x[j]) return 0;
	i += n;
	j += n;
	if (i > k) break;
    }
    return 1;
}

/* *index must be initialised with 0's; on output contains the indices (1, ... n) of decreasing values of *x */
void order_int(int *x, int *index, int n)
{
    int i, k, m = 1;

    while (m <= n) {
	k = 0;
	while (index[k]) k++;
	if (m == n) {
	    index[k] = m;
	} else {
	    for (i = 0; i < n; i++)
		if (index[i] == 0 && x[i] < x[k]) k = i;
	    index[k] = m;
	}
	m++;
    }
}

int anyElementZero(int *x, int n)
{
    int i;
    for (i = 0; i < n; i++) {
	if (!x[i]) return 1;
    }
    return 0;
}

void update_dist_mat(int *D, int n, int j)
{
    int i;
    for (i = 0; i < n; i++) D[i + n * j] = 1;
    for (i = 0; i < n; i++) D[j + n * i] = 1;
}

/* Description of the algorihm:
   https://www.mail-archive.com/r-sig-phylo@r-project.org/msg05541.html */
void haplotype_DNAbin(unsigned char *x, int *n, int *s, int *haplo, int *warn, int *strict)
{
    int i, j, k, *ihap, Nhaplo, new_haplo;
    /* i: index of the sequence, j: index of the haplotype */

    /* step 1 */
    ihap = (int *)R_alloc(*n, sizeof(int)); /* safe allocation */
    /* ihap stores the C-indices of the haplotype sequences:
       ihap[i] = j means that the sequence of the i-th haplotype
       is given by the j-th sequence of x */

    /* the 1st sequence is always the 1st haplotype */
    ihap[0] = 0;
    Nhaplo = 1;

    for (i = 1; i < *n; i++) { /* start with the 2nd sequence */
	new_haplo = 1;
	for (k = 0; k < Nhaplo; k++) { /* loop among the haplotypes */
	    j = ihap[k];
	    if (identical_seqs(x, i, j, *n, *s)) {
		haplo[i] = j + 1;
		new_haplo = 0;
		break;
	    }
	}
	if (new_haplo) {
	    ihap[Nhaplo] = i;
	    Nhaplo++;
	}
    }

    if (*strict) return;

    /* step 2 is already processed at the R level */

    /* step 3 */
    int *D, ii, jj, S, Ndist = Nhaplo*Nhaplo;
    D = (int *)R_alloc(Ndist, sizeof(int));
    /* initialise the diagonal with 1's: */
    for (i = 0; i < Nhaplo; i++) D[i + Nhaplo * i] = 1;

    for (i = 0; i < Nhaplo - 1; i++) {
	for (j = i + 1; j < Nhaplo; j++) {
	    ii = ihap[i];
	    jj = ihap[j];
	    S = 0;
	    k = ii + *n * (*s - 1); /* the last site of seq i */
	    /* no need to compute the real distance,
	       just need to know if it is > 0 */
	    while (ii <= k) {
		if (DifferentBase(x[ii], x[jj])) {
		    S++;
		    break;
		}
		ii += *n;
		jj += *n;
	    }
	    D[i + Nhaplo * j] = D[j + Nhaplo * i] = S;
	}
    }

    if (! anyElementZero(D, Ndist)) return;

    int *NknownBases;
    NknownBases = (int *)R_alloc(Nhaplo, sizeof(int));
    memset(NknownBases, 0, Nhaplo * sizeof(int));
    for (i = 0; i < Nhaplo; i++) {
	j = ihap[i];
	k = j + *n * (*s - 1);
	while (j <= k) {
	    if (KnownBase(x[j])) (NknownBases[i])++;
	    j += *n;
	}
    }

    int *index;
    index = (int *)R_alloc(Nhaplo, sizeof(int));
    memset(index, 0, Nhaplo * sizeof(int));
    /* order of haplotypes with decreasing numbers of missing data: */
    order_int(NknownBases, index, Nhaplo);

    /* a buffer to store haplotype indices below: */
    int NdistZero, *buf;
    buf = (int *)R_alloc(Nhaplo, sizeof(int));

    /* step 4 */
    for (;;) {
	for (i = 0; i < Nhaplo; i++) {
	    ii = 0;
	    /* find the haplotypes with the largest number of missing data first */
	    while (index[ii] - 1 != i) ii++;
	    NdistZero = 0;
	    /* compare this haplotype (ii) with all the others */
	    for (j = 0; j < Nhaplo; j++) {
		if (j == ii) continue; /* skip the diagonal */
		if (!D[ii + Nhaplo * j]) {
		    buf[NdistZero] = j;
		    /* buf stores the haplotype indices that have 0-dist
		       with haplotype ii */
		    NdistZero++;
		}
	    }
	    /* step 5 */
	    if (NdistZero > 1) {
/* check that the 2+ haplotypes with dist=0 to 'ii' don't have also dist=0 among them;
   if yes then all are pooled in the same haplotype, if no keep all separate */
		int allOthersZero = 1;
		/* use 'j' and 'k' for the indices of the haplotype pair */
		for (jj = 0; jj < NdistZero - 1; jj++) {
		    j = buf[jj];
		    for (S = jj + 1; S < NdistZero; S++) {
			k = buf[S];
			if (D[j + Nhaplo * k]) allOthersZero = 0;
		    }
		}
		if (allOthersZero) NdistZero = 1; else warn[1] = 1;
/* just set NdistZero = 1 and the haplotypes will be pooled together at a later iteration */
	    }
	    if (NdistZero == 1) {
		int z = ihap[buf[0]] + 1;
		warn[0] = 1;
		k = ihap[ii];
		haplo[k] = z;
		k++;
		for (jj = 0; jj < *n; jj++)
		    if (haplo[jj] == k) haplo[jj] = z;
	    }
	    update_dist_mat(D, Nhaplo, ii);
	    if (!anyElementZero(D, Ndist)) return;
	}
    }
}

void distDNA_pegas(unsigned char *x, int *n, int *s, double *d)
{
    int i1, i2, s1, s2, target, Nd;
    target = 0;
    for (i1 = 1; i1 < *n; i1++) {
	for (i2 = i1 + 1; i2 <= *n; i2++) {
	    Nd = 0;
	    for (s1 = i1 - 1, s2 = i2 - 1; s1 < i1 + *n*(*s - 1); s1 += *n, s2 += *n) {
		if (x[s1] <= 0x07 || x[s2] <= 0x07) continue;
		if (DifferentBase(x[s1], x[s2])) Nd++;
	    }
	    d[target] = ((double) Nd);
	    target++;
	}
    }
}

SEXP unique_haplotype_loci(SEXP x, SEXP NROW, SEXP NCOL)
{
    SEXP res;
    int nr, nc, i, j, k, *H, flag;

    PROTECT(x = coerceVector(x, STRSXP));
    PROTECT(NROW = coerceVector(NROW, INTSXP));
    PROTECT(NCOL = coerceVector(NCOL, INTSXP));

    nr = INTEGER(NROW)[0];
    nc = INTEGER(NCOL)[0];

    PROTECT(res = allocVector(INTSXP, nc));
    H = INTEGER(res);
    memset(H, 0, nc * sizeof(int));

    j = 0;
    while (j < nc - 1) {
	if (!H[j]) {
	    k = j + 1;
	    while (k < nc) {
		if (!H[k]) {
		    i = 0;
		    flag = 1; /* initially the two haplotypes are considered identical */
		    while (i < nr) {
			if (strcmp(CHAR(STRING_ELT(x, i + j*nr)), CHAR(STRING_ELT(x, i + k*nr)))) { /* the strings are different */
			    flag = 0;
			    break;
			}
			i++;
		    }
		    if (flag) H[k] = j + 1;
		}
		k++;
	    }
	}
	j++;
    }

    UNPROTECT(4);
    return res;
}

SEXP alreadyIn_mjn_DNAbin(SEXP x, SEXP table)
{
    int n1, n2, s, i1, i2, j;
    unsigned char *xp, *tablep;
    SEXP res;
    int *y, diff_seq = 0;

    PROTECT(x = coerceVector(x, RAWSXP));
    PROTECT(table = coerceVector(table, RAWSXP));
    xp = RAW(x);
    tablep = RAW(table);

    n1 = nrows(x);
    n2 = nrows(table);
    s = ncols(x);

    PROTECT(res = allocVector(INTSXP, n1));
    y = INTEGER(res);
    memset(y, 0, n1 * sizeof(int));

    for (i1 = 0; i1 < n1; i1++) {
	for (i2 = 0; i2 < n2; i2++) {
	    for (j = 0; j < s; j++) {
		if (DifferentBase(xp[i1 + j * n1], tablep[i2 + j * n2])) {
		    diff_seq = 1;
		    break;
		}
	    }
	    if (diff_seq) {
		diff_seq = 0;
		continue;
	    } else {
		y[i1] = 1;
		break;
	    }
	}
	if (y[i1]) continue;
    }

    UNPROTECT(3);
    return res;
}

void getMedianVectors_DNAbin_mjn(unsigned char *x, int *p, unsigned char *res, int *flag)
{
    int j, k, s = p[0];
    unsigned char z;

    for (j = 0; j < s; j++) {
	k = j * 3; /* there are 3 seqs */
	z = x[k];
	if (z == x[k + 1] || z == x[k + 2]) {
	    res[k] = res[k + 1] = res[k + 2] = z;
	} else {
	    z = x[k + 1];
	    if (z == x[k + 2]) {
		res[k] = res[k + 1] = res[k + 2] = z;
	    } else {
		res[k] = x[k]; k++;
		res[k] = x[k]; k++;
		res[k] = x[k];
		*flag = 0;
	    }
	}
    }
}


#include <R.h>
#include <Rinternals.h>

/* fast tabulation (4 times faster than tabulate(); see write.vcf.R
   ===> code unused <===
SEXP C_tabulate_pegas(SEXP x, SEXP ncat) {
    SEXP res;
    PROTECT(x = coerceVector(x, INTSXP));
    PROTECT(ncat = coerceVector(ncat, INTSXP));
    int n = LENGTH(x);
    int nc = INTEGER(ncat)[0];
    PROTECT(res = allocVector(INTSXP, nc));
    int *y = INTEGER(res);
    memset(y, 0, nc * sizeof(int));
    int *xp = INTEGER(x);
    for (int i = 0; i < n; i++) ++y[xp[i]]; // O's are also counted
    UNPROTECT(3);
    return res;
    }*/

/* two functions to compare sprintf() and int2str():

void funfoo(int *x, char **y, int *n)
{
    for (int i = 0; i < *n; i++) sprintf(y[i], "%d", x[i]);
}

SEXP funbar(SEXP x) {
    SEXP res;
    PROTECT(x = coerceVector(x, INTSXP));
    int n = LENGTH(x);
    PROTECT(res = allocVector(STRSXP, n));
    int * xp = INTEGER(x);
    char y[100];
    for (int i = 0; i < n; i++) {
	int2str(xp[i], y);
	SET_STRING_ELT(res, i, mkChar(y));
    }
    UNPROTECT(2);
    return res;
}
*/

/* translates integers 0..999 into a char string
   => 8 times faster than using sprintf() */
int int2str(int x, char *str)
{
    static char digits[] = {'0', '1', '2', '3', '4', '5', '6', '7', '8', '9'};
    if (x < 10) {
	str[0] = digits[x];
	str[1] = '\0';
	return 1;
    }
    if (x < 100) {
	str[0] = digits[x/10];
	str[1] = digits[x % 10];
	str[2] = '\0';
	return 2;
    }
    if (x < 1000) {
	str[0] = digits[x/100];
	str[1] = digits[x/10 % 10];
	str[2] = digits[x % 10];
	str[3] = '\0';
	return 3;
    }
    return 0;
}

void translateGenotypesForVCF(char **old, char **new, int *ngeno, char **alleles, int *nall)
{
    int a = 0, b = 1, i = 0, j, k, l, in = 0, nd;
    char buf[100], c;

    for (;;) {
    if (i >= *ngeno) break;
	for (;;) {
	    c = old[i][b];
	    if (c == '/' || c == '|' || c == '\0')
		break;
	    b++;
	}
	for (k = 0, l = a; l < b; k++, l++)
	    buf[k] = old[i][l];
	buf[k] = '\0';
	if (strcmp(".\0", buf)) {
	    for (j = 0; j < *nall; j++) {
		if (! strcmp(alleles[j], buf)) break;
	    }
	    nd = int2str(j, buf);
	    for (k = 0; k < nd; k++, in++) new[i][in] = buf[k];
	} else {
	    new[i][in++] = '.';
	    new[i][in] = '\0';
	}
	c = old[i][b];
	new[i][in++] = c;
	if (c == '\0') {
	    i++;
	    in = 0;
	    a = 0;
	    b = 1;
	} else {
	    a = b + 1;
	    b = a + 1;
	}
    }
}

static R_CMethodDef C_entries[] = {
    {"haplotype_DNAbin", (DL_FUNC) &haplotype_DNAbin, 6},
    {"distDNA_pegas", (DL_FUNC) &distDNA_pegas, 4},
    {"getMedianVectors_DNAbin_mjn", (DL_FUNC) &getMedianVectors_DNAbin_mjn, 4},
    {"translateGenotypesForVCF", (DL_FUNC) &translateGenotypesForVCF, 5},
    {NULL, NULL, 0}
};

SEXP read_bin_pegas(SEXP FILENAME, SEXP SIZE, SEXP SKIP);
SEXP findEOL_C(SEXP x, SEXP SKIP, SEXP HOP);
SEXP extract_POS(SEXP x, SEXP EOL, SEXP nTABtoSKIP);
SEXP extract_REF(SEXP x, SEXP EOL, SEXP nTABtoSKIP);
SEXP build_factor_loci(SEXP x, SEXP N);
SEXP summary_loci_pegas(SEXP x, SEXP LOCI);

static R_CallMethodDef Call_entries[] = {
    {"read_bin_pegas", (DL_FUNC) &read_bin_pegas, 3},
    {"findEOL_C", (DL_FUNC) &findEOL_C, 3},
    {"extract_POS", (DL_FUNC) &extract_POS, 3},
    {"extract_REF", (DL_FUNC) &extract_REF, 3},
    {"build_factor_loci", (DL_FUNC) &build_factor_loci, 2},
    {"unique_haplotype_loci", (DL_FUNC) &unique_haplotype_loci, 3},
    {"summary_loci_pegas", (DL_FUNC) &summary_loci_pegas, 2},
    {"alreadyIn_mjn_DNAbin", (DL_FUNC) &alreadyIn_mjn_DNAbin, 2},
    {NULL, NULL, 0}
};

void R_init_pegas(DllInfo *info)
{
    R_registerRoutines(info, C_entries, Call_entries, NULL, NULL);
    R_useDynamicSymbols(info, FALSE);
}
