/* fftlib.f -- translated by f2c (version of 25 November 1991  18:53:50).
 * 6 August 1992 by Michael Pique (mp@scripps.edu)
 *  Used by permission of The Scripps Research Institute */
#include <cstdio>
#include <cmath>

#ifndef F2C_INCLUDE
#define F2C_INCLUDE

typedef long int integer;
typedef char* address;
typedef short int shortint;
typedef float real;
typedef double doublereal;
typedef struct
{
  real r, i;
} micomplex;
typedef struct
{
  doublereal r, i;
} doublecomplex;
typedef long int logical;
typedef short int shortlogical;

/* Extern is for use with -E */
#ifndef Extern
#define Extern extern
#endif

/* I/O stuff */

#ifdef f2c_i2
/* for -i2 */
typedef short flag;
typedef short ftnlen;
typedef short ftnint;
#else
typedef long flag;
typedef long ftnlen;
typedef long ftnint;
#endif


/*external read, write*/
typedef struct
{
  flag cierr;
  ftnint ciunit;
  flag ciend;
  char* cifmt;
  ftnint cirec;
} cilist;

/*internal read, write*/
typedef struct
{
  flag icierr;
  char* iciunit;
  flag iciend;
  char* icifmt;
  ftnint icirlen;
  ftnint icirnum;
} icilist;

/*open*/
typedef struct
{
  flag oerr;
  ftnint ounit;
  char* ofnm;
  ftnlen ofnmlen;
  char* osta;
  char* oacc;
  char* ofm;
  ftnint orl;
  char* oblnk;
} olist;

/*close*/
typedef struct
{
  flag cerr;
  ftnint cunit;
  char* csta;
} cllist;

/*rewind, backspace, endfile*/
typedef struct
{
  flag aerr;
  ftnint aunit;
} alist;

/* inquire */
typedef struct
{
  flag inerr;
  ftnint inunit;
  char* infile;
  ftnlen infilen;
  ftnint* inex;  /*parameters in standard's order*/
  ftnint* inopen;
  ftnint* innum;
  ftnint* innamed;
  char* inname;
  ftnlen innamlen;
  char* inacc;
  ftnlen inacclen;
  char* inseq;
  ftnlen inseqlen;
  char* indir;
  ftnlen indirlen;
  char* infmt;
  ftnlen infmtlen;
  char* inform;
  ftnint informlen;
  char* inunf;
  ftnlen inunflen;
  ftnint* inrecl;
  ftnint* innrec;
  char* inblank;
  ftnlen inblanklen;
} inlist;

#define VOID void

union Multitype     /* for multiple entry points */
{
  shortint h;
  integer i;
  real r;
  doublereal d;
  micomplex c;
  doublecomplex z;
};

typedef union Multitype Multitype;

typedef long Long;  /* No longer used; formerly in Namelist */

struct Vardesc      /* for Namelist */
{
  char* name;
  char* addr;
  ftnlen* dims;
  int type;
};
typedef struct Vardesc Vardesc;

struct Namelist
{
  char* name;
  Vardesc** vars;
  int nvars;
};
typedef struct Namelist Namelist;

#define abs(x) ((x) >= 0 ? (x) : -(x))
#define dabs(x) (doublereal)abs(x)

/* procedure parameter types for -A and -C++ */

#define F2C_proc_par_types 1

typedef int /* Unknown procedure type */ (*U_fp)(...);
typedef shortint (*J_fp)(...);
typedef integer (*I_fp)(...);
typedef real (*R_fp)(...);
typedef doublereal (*D_fp)(...), (*E_fp)(...);
typedef /* Complex */ VOID (*C_fp)(...);
typedef /* Double Complex */ VOID (*Z_fp)(...);
typedef logical (*L_fp)(...);
typedef shortlogical (*K_fp)(...);
typedef /* Character */ VOID (*H_fp)(...);
typedef /* Subroutine */ int (*S_fp)(...);

/* E_fp is for real functions when -R is not specified */
typedef VOID C_f;   /* complex function */
typedef VOID H_f;   /* character function */
typedef VOID Z_f;   /* double complex function */
typedef doublereal E_f; /* real function with -R not specified */

/* undef any lower-case symbols that your C compiler predefines, e.g.: */

#ifndef Skip_f2c_Undefs
#undef cray
#undef gcos
#undef mc68010
#undef mc68020
#undef mips
#undef pdp11
#undef sgi
#undef sparc
#undef sun
#undef sun2
#undef sun3
#undef sun4
#undef u370
#undef u3b
#undef u3b2
#undef u3b5
#undef unix
#undef vax
#endif
#endif
/* end of f2c.h */

/* Subroutine */
int cmplft_(real* x, real* y, integer* n, integer* d);
/* Subroutine */
int srfp_(integer* pts, integer* pmax, integer* twogrp, integer* factor, integer* sym, integer* psym, integer* unsym, logical* error);
/* Subroutine */
int mdftkd_(integer* n, integer* factor, integer* dim, real* x, real* y);
/* Subroutine */
int diprp_(integer* pts, integer* sym, integer* psym, integer* unsym, integer* dim, real* x, real* y);
/* Subroutine */
int r2cftk_(integer* n, integer* m, real* x0, real* y0, real* x1, real* y1, integer* dim);
/* Subroutine */
int r3cftk_(integer* n, integer* m, real* x0, real* y0, real* x1, real* y1, real* x2, real* y2, integer* dim);
/* Subroutine */
int r4cftk_(integer* n, integer* m, real* x0, real* y0, real* x1, real* y1, real* x2, real* y2, real* x3, real* y3, integer* dim);
/* Subroutine */
int r5cftk_(integer* n, integer* m, real* x0, real* y0, real* x1, real* y1, real* x2, real* y2, real* x3, real* y3, real* x4, real* y4,
            integer* dim);
/* Subroutine */
int r8cftk_(integer* n, integer* m, real* x0, real* y0, real* x1, real* y1, real* x2, real* y2, real* x3, real* y3, real* x4, real* y4, real* x5,
            real* y5, real* x6, real* y6, real* x7, real* y7, integer* dim);
/* Subroutine */
int rpcftk_(integer* n, integer* m, integer* p, integer* r, real* x, real* y, integer* dim);

/* Table of constant values */

/* Subroutine */
int realft_(real* even, real* odd, integer* n, integer* dim) {
  /* System generated locals */
  integer i__1, i__2, i__3, i__4, i__5;

  /* Local variables */
  static real twon, a, b, c, d, e, f;
  static integer i, j, k, l;
  static real angle;
  static integer d2, d3, d4, d5, i0, i1, i2;
  static real twopi;
  static integer nover2;
  static real co, si;
  static integer nt;


  /*     REAL FOURIER TRANSFORM */

  /*     GIVEN A REAL SEQUENCE OF LENGTH 2N THIS SUBROUTINE CALCULATES THE
   */
  /*     UNIQUE PART OF THE FOURIER TRANSFORM.  THE FOURIER TRANSFORM HAS */

  /*     N + 1 UNIQUE REAL PARTS AND N - 1 UNIQUE IMAGINARY PARTS.  SINCE */

  /*     THE REAL PART AT X(N) IS FREQUENTLY OF INTEREST, THIS SUBROUTINE */

  /*     STORES IT AT X(N) RATHER THAN IN Y(0).  THEREFORE X AND Y MUST BE
   */
  /*     OF LENGTH N + 1 INSTEAD OF N.  NOTE THAT THIS STORAGE ARRANGEMENT
   */
  /*     IS DIFFERENT FROM THAT EMPLOYED BY THE HERMITIAN FOURIER TRANSFORM
   */
  /*     SUBROUTINE. */

  /*     FOR CONVENIENCE THE DATA IS PRESENTED IN TWO PARTS, THE FIRST */
  /*     CONTAINING THE EVEN NUMBERED REAL TERMS AND THE SECOND CONTAINING
   */
  /*     THE ODD NUMBERED TERMS (NUMBERING STARTING AT 0).  ON RETURN THE */

  /*     REAL PART OF THE TRANSFORM REPLACES THE EVEN TERMS AND THE */
  /*     IMAGINARY PART OF THE TRANSFORM REPLACES THE ODD TERMS. */


  /* Parameter adjustments */
  --dim;
  --odd;
  --even;

  /* Function Body */
  twopi = (float)6.283185;
  twon = (real) (*n << 1);

  cmplft_(&even[1], &odd[1], n, &dim[1]);

  nt = dim[1];
  d2 = dim[2];
  d3 = dim[3];
  d4 = dim[4] - 1;
  d5 = dim[5];
  nover2 = *n / 2 + 1;

  if (nover2 < 2) {
    goto L400;
  }
  i__1 = nover2;
  for (i = 2; i <= i__1; ++i) {
    angle = twopi * (real) (i - 1) / twon;
    co = (float)cos(angle);
    si = (float)sin(angle);
    i0 = (i - 1) * d2 + 1;
    j = (*n + 2 - (i << 1)) * d2;
    i__2 = nt;
    i__3 = d3;
    for (i1 = i0; i__3 < 0 ? i1 >= i__2 : i1 <= i__2; i1 += i__3) {
      i2 = i1 + d4;
      i__4 = i2;
      i__5 = d5;
      for (k = i1; i__5 < 0 ? k >= i__4 : k <= i__4; k += i__5) {
        l = k + j;
        a = (even[l] + even[k]) / (float)2.;
        c = (even[l] - even[k]) / (float)2.;
        b = (odd[l] + odd[k]) / (float)2.;
        d = (odd[l] - odd[k]) / (float)2.;
        e = c * si + b * co;
        f = c * co - b * si;
        even[k] = a + e;
        even[l] = a - e;
        odd[k] = f - d;
        odd[l] = f + d;
        /* L100: */
      }
      /* L200: */
    }
    /* L300: */
  }

L400:
  if (*n < 1) {
    goto L600;
  }
  j = *n * d2;
  i__1 = nt;
  i__3 = d3;
  for (i1 = 1; i__3 < 0 ? i1 >= i__1 : i1 <= i__1; i1 += i__3) {
    i2 = i1 + d4;
    i__2 = i2;
    i__5 = d5;
    for (k = i1; i__5 < 0 ? k >= i__2 : k <= i__2; k += i__5) {
      l = k + j;
      even[l] = even[k] - odd[k];
      odd[l] = (float)0.;
      even[k] += odd[k];
      odd[k] = (float)0.;
      /* L500: */
    }
  }

L600:
  return 0;

} /* realft_ */

/* Subroutine */
int hermft(real* x, real* y, integer* n, integer* dim) {
  /* System generated locals */
  integer i__1, i__2, i__3, i__4, i__5;

  /* Local variables */
  static real twon, a, b, c, d, e, f;
  static integer i, j, k;
  static real angle;
  static integer d2, d3, d4, d5, i0, i1, i2, k1;
  static real twopi;
  static integer nover2;
  static real co, si;
  static integer nt;


  /*     HERMITIAN SYMMETRIC FOURIER TRANSFORM */

  /*     GIVEN THE UNIQUE TERMS OF A HERMITIAN SYMMETRIC SEQUENCE OF LENGTH
   */
  /*     2N THIS SUBROUTINE CALCULATES THE 2N REAL NUMBERS WHICH ARE ITS */
  /*     FOURIER TRANSFORM.  THE EVEN NUMBERED ELEMENTS OF THE TRANSFORM */
  /*     (0, 2, 4, . . ., 2N-2) ARE RETURNED IN X AND THE ODD NUMBERED */
  /*     ELEMENTS (1, 3, 5, . . ., 2N-1) IN Y. */

  /*     A FINITE HERMITIAN SEQUENCE OF LENGTH 2N CONTAINS N + 1 UNIQUE */
  /*     REAL NUMBERS AND N - 1 UNIQUE IMAGINARY NUMBERS.  FOR CONVENIENCE
   */
  /*     THE REAL VALUE FOR X(N) IS STORED AT Y(0). */


  /* Parameter adjustments */
  --dim;
  --y;
  --x;

  /* Function Body */
  twopi = (float)6.283185;
  twon = (real) (*n << 1);

  nt = dim[1];
  d2 = dim[2];
  d3 = dim[3];
  d4 = dim[4] - 1;
  d5 = dim[5];

  i__1 = nt;
  i__2 = d3;
  for (i0 = 1; i__2 < 0 ? i0 >= i__1 : i0 <= i__1; i0 += i__2) {
    i1 = i0 + d4;
    i__3 = i1;
    i__4 = d5;
    for (i = i0; i__4 < 0 ? i >= i__3 : i <= i__3; i += i__4) {
      a = x[i];
      b = y[i];
      x[i] = a + b;
      y[i] = a - b;
      /* L100: */
    }
  }

  nover2 = *n / 2 + 1;
  if (nover2 < 2) {
    goto L500;
  }
  i__4 = nover2;
  for (i0 = 2; i0 <= i__4; ++i0) {
    angle = twopi * (real) (i0 - 1) / twon;
    co = (float)cos(angle);
    si = (float)sin(angle);
    k = (*n + 2 - (i0 << 1)) * d2;
    k1 = (i0 - 1) * d2 + 1;
    i__3 = nt;
    i__2 = d3;

    for (i1 = k1; i__2 < 0 ? i1 >= i__3 : i1 <= i__3; i1 += i__2) {
      i2 = i1 + d4;
      i__1 = i2;
      i__5 = d5;
      for (i = i1; i__5 < 0 ? i >= i__1 : i <= i__1; i += i__5) {
        j = i + k;
        a = x[i] + x[j];
        b = x[i] - x[j];
        c = y[i] + y[j];
        d = y[i] - y[j];
        e = b * co + c * si;
        f = b * si - c * co;
        x[i] = a + f;
        x[j] = a - f;
        y[i] = e + d;
        y[j] = e - d;
        /* L200: */
      }
      /* L300: */
    }
    /* L400: */
  }

  cmplft_(&x[1], &y[1], n, &dim[1]);

L500:
  return 0;

} /* hermft */

/* Subroutine */
int sdiad_(real* x, real* y, integer* n, integer* dim) {

  /* System generated locals */
  integer i__1, i__2, i__3, i__4, i__5;


  /* Local variables */
  static logical fold;
  static real twon, a, c;
  static integer i, j, k, l, m;
  static real s, angle;
  static integer d1, d2, d3, d4, d5, k0, k1, k2;
  static real twopi;
  static integer nover2, mn, nn;




  /*     THIS SUBROUTINE COMPUTES HALF THE FOURIER SYNTHESIS ALONG A SCREW
   */
  /*     DIAD LYING ALONG A CRYSTALLOGRAPHIC AXIS GIVEN HALF THE FOURIER */
  /*     COEFFICIENTS.  THAT IS, IT ASSUMES THAT F(T) = CONJG(F(-T)) FOR T
   */
  /*     EVEN AND F(T) = -CONJG(F(-T)) FOR T ODD.  N IS THE LENGTH OF THE */

  /*     DESIRED HALF OF THE TRANSFORM.  THE LOCATION X(N+1) IS REQUIRED AS
   */
  /*     A SCRATCH LOCATION AND THEREFORE A VALUE IS ALSO RETURNED IN */
  /*     X(N+1) AND Y(N+1).  THE VALUE OF THE SECOND HALF OF THE TRANSFORM
   */
  /*     MAY BE GENERATED FROM THE FIRST HALF BY THE FORMULA X(N+T) = X(T),
   */
  /*     Y(N+T) = -Y(T).  IN OTHER WORDS, THE LAST HALF OF THE TRANSFORM IS
   */
  /*     THE COMPLEX CONJUGATE OF THE FIRST HALF. */

  /*     THE TRANSFORM IS CALCULATED BY FORMING THE SUM OF THE EVEN TERMS */

  /*     AND THE ODD TERMS IN PLACE, USING THE SYMMETRY RELATIONS TO */
  /*     OBTAIN THE VALUES FOR NEGATIVE SUBSCRIPTS.  THE TRANSFORM OF THE */

  /*     RESULTING SEQUENCE MAY BE SEPARATED BY USING THE FACT THAT THE */
  /*     TRANSFORM OF THE EVEN TERMS IS REAL, WHILE THE PRODCT OF THE */
  /*     TRANSFORM OF THE ODD TERMS AND (COS(PI*T/N) - I*SIN(PI*T/N)) IS */
  /*     IMAGINARY.  THE SCRATCH LOCATION IS REQUIRED BECAUSE THE FORMULA */

  /*     FOR SEPARATING THE TWO TRANSFORMS BREAKS DOWN WHEN T = N/2. */


  /* Parameter adjustments */
  --dim;
  --y;
  --x;

  /* Function Body */
  nover2 = *n / 2;
  if (nover2 << 1 != *n) {
    goto L700;
  }
  twon = (real) (*n << 1);
  twopi = (float)6.2831852;
  d1 = dim[1];
  d2 = dim[2];
  d3 = dim[3];
  d4 = dim[4] - 1;
  d5 = dim[5];

  k0 = (*n - 1) * d2 + 1;
  i__1 = d1;
  i__2 = d3;
  for (k1 = k0; i__2 < 0 ? k1 >= i__1 : k1 <= i__1; k1 += i__2) {
    k2 = k1 + d4;
    i__3 = k2;
    i__4 = d5;
    for (k = k1; i__4 < 0 ? k >= i__3 : k <= i__3; k += i__4) {
      l = k + d2;
      x[l] = x[k];
      /* L100: */
    }
  }
  s = (float)1.;
  nn = *n - 2;
  i__4 = nn;
  for (i = 1; i <= i__4; i += 2) {
    s = -s;
    mn = (*n + 1 - i) * d2;
    k0 = (i - 1) * d2 + 1;
    i__3 = d1;
    i__2 = d3;
    for (k1 = k0; i__2 < 0 ? k1 >= i__3 : k1 <= i__3; k1 += i__2) {
      k2 = k1 + d4;
      i__1 = k2;
      i__5 = d5;
      for (k = k1; i__5 < 0 ? k >= i__1 : k <= i__1; k += i__5) {
        j = k + d2;
        l = k + (d2 << 1);
        m = k + mn;
        x[m] += s * x[j];
        x[k] += x[j];
        x[j] = x[l] - x[j];
        y[k] += y[j];
        y[j] -= y[l];
        /* L150: */
      }
    }
    /* L200: */
  }
  k0 = (*n - 2) * d2 + 1;
  i__4 = d1;
  i__5 = d3;
  for (k1 = k0; i__5 < 0 ? k1 >= i__4 : k1 <= i__4; k1 += i__5) {
    k2 = k1 + d4;
    i__1 = k2;
    i__2 = d5;
    for (k = k1; i__2 < 0 ? k >= i__1 : k <= i__1; k += i__2) {
      l = k + d2;
      x[k] += x[l];
      y[k] += y[l];
      x[l] = -x[l];
      /* L250: */
    }
  }

  /*     REORDER SCRAMBLED FOURIER COEFFICIENTS */

  i__2 = nn;
  for (i = 1; i <= i__2; ++i) {
    k = i;
L300:
    k <<= 1;
    if (k > *n - 1) {
      k = (*n << 1) - 1 - k;
    }
    if (k < i) {
      goto L300;
    }
    if (k == i) {
      goto L400;
    }
    j = (k - i) * d2;
    k0 = i * d2 + 1;
    i__1 = d1;
    i__5 = d3;
    for (k1 = k0; k1 <= i__1; k1 += i__5) {
      k2 = k1 + d4;
      i__4 = k2;
      i__3 = d5;
      for (k = k1;  k <= i__4; k += i__3) {
        l = k + j;
        a = x[k];
        x[k] = x[l];
        x[l] = a;
        a = y[k];
        y[k] = y[l];
        y[l] = a;
        /* L350: */
      }
    }
L400:
    ;
  }

  cmplft_(&x[1], &y[1], n, &dim[1]);

  m = nover2 - 1;
  i__2 = m;
  for (i = 1; i <= i__2; ++i) {
    angle = twopi * (real) i / twon;
    c = (float)cos(angle);
    s = (float)sin(angle);
    k0 = i * d2 + 1;
    fold = true;
    goto L500;

L450:
    c = -c;
    k0 = (*n - i) * d2 + 1;
    fold = false;

L500:
    i__3 = d1;
    i__4 = d3;
    for (k1 = k0; k1 <= i__3; k1 += i__4) {
      k2 = k1 + d4;
      i__5 = k2;
      i__1 = d5;
      for (k = k1; k <= i__5; k += i__1) {
        a = y[k] / c;
        x[k] += s * a;
        y[k] = a;
        /* L550: */
      }
    }
    if (fold) {
      goto L450;
    }
    /* L600: */
  }

  m = nover2 * d2;
  k0 = m + 1;
  i__2 = d1;
  i__1 = d3;
  for (k1 = k0; k1 <= i__2; k1 += i__1) {
    k2 = k1 + d4;
    i__5 = k2;
    i__4 = d5;
    for (k = k1; k <= i__5; k += i__4) {
      j = k - m;
      l = k + m;
      a = x[l] * (float)2.;
      x[k] += a;
      y[k] = a;
      x[l] = x[j];
      y[l] = -y[j];
      /* L650: */
    }
  }

  return 0;

L700:
  fprintf(stderr, "sdiad n odd.  n = %ld\n", *n);
  return (-1);

} /* sdiad_ */

/* Subroutine */
int inv21_(real* x, real* y, integer* n, integer* d) {
  /* System generated locals */
  integer i__1, i__2, i__3, i__4, i__5;

  /* Local variables */
  static real a, b, c;
  static integer i, j, k, l, m;
  static real r, s, c1;
  static integer d1, d2, d3, d4, d5, j1, j2, j3;
  static real s1;
  static integer nover2, kk, ll;
  static real pi;

  /*     INVERTS FOURIER TRANSFORM ALONG A SCREW */
  /*     DIAD. THE RESULT IS SCALED BY N. */


  /* Parameter adjustments */
  --d;
  --y;
  --x;

  /* Function Body */
  pi = (float)3.141593;

  d1 = d[1];
  d2 = d[2];
  d3 = d[3];
  d4 = d[4] - 1;
  d5 = d[5];

  nover2 = *n / 2;
  ll = *n * d2;
  kk = nover2 * d2;
  i__1 = d1;
  i__2 = d3;
  for (j1 = 1; j1 <= i__1; j1 += i__2) {
    j2 = j1 + d4;
    i__3 = j2;
    i__4 = d5;
    for (j = j1; j <= i__3; j += i__4) {
      l = ll + j;
      k = kk + j;
      x[l] = x[j] + x[k];
      x[k] += y[k];
      y[l] = (float)0.;
      /* L100: */
      y[k] = (float)0.;
    }
  }

  c1 = (float)cos(pi / (real) (*n));
  s1 = (float)sin(pi / (real) (*n));
  c = (float)1.;
  s = (float)0.;
  i__4 = nover2;
  for (i = 2; i <= i__4; ++i) {
    kk = (*n + 2 - (i << 1)) * d2;
    ll = (*n + 1 - i) * d2;
    r = c * c1 - s * s1;
    s = c * s1 + s * c1;
    c = r;
    j1 = (i - 1) * d2 + 1;
    i__3 = d1;
    i__2 = d3;
    for (j2 = j1; j2 <= i__3; j2 += i__2) {
      j3 = j2 + d4;
      i__1 = j3;
      i__5 = d5;
      for (j = j2; j <= i__1; j += i__5) {
        l = j + ll;
        k = j + kk;
        x[l] = x[l] + x[j] + x[k];
        x[j] += s * y[j];
        x[k] += s * y[k];
        y[j] = c * y[j];
        /* L150: */
        y[k] = -c * y[k];
      }
    }
    /* L200: */
  }

  cmplft_(&x[1], &y[1], n, &d[1]);

  i__4 = nover2;
  for (i = 1; i <= i__4; ++i) {
    kk = (*n + 1 - (i << 1)) * d2;
    ll = kk + i * d2;
    j1 = (i - 1) * d2 + 1;
    i__5 = d1;
    i__1 = d3;
    for (j2 = j1; j2 <= i__5; j2 += i__1) {
      j3 = j2 + d4;
      i__2 = j3;
      i__3 = d5;
      for (j = j2; j <= i__2; j += i__3) {
        k = j + kk;
        l = j + ll;
        a = x[j] - x[l];
        b = y[j] + y[l];
        x[j] = x[l];
        y[j] = -y[l];
        x[l] = x[k] + a;
        y[l] = y[k] - b;
        x[k] = a;
        /* L250: */
        y[k] = b;
      }
    }
    /* L300: */
  }

  m = *n - 2;
  i__4 = m;
  for (i = 1; i <= i__4; ++i) {
    k = i;
L320:
    j = k;
    k = j / 2;
    if (k << 1 != j) {
      k = *n - 1 - k;
    }
    if ((i__3 = k - i) < 0) {
      goto L320;
    } else if (i__3 == 0) {
      goto L400;
    } else {
      goto L340;
    }
L340:
    kk = (k - i) * d2;
    j1 = i * d2 + 1;
    i__3 = d1;
    i__2 = d3;
    for (j2 = j1; j2 <= i__3; j2 += i__2) {
      j3 = j2 + d4;
      i__1 = j3;
      i__5 = d5;
      for (j = j2; j <= i__1; j += i__5) {
        k = j + kk;
        a = x[k];
        b = y[k];
        x[k] = x[j];
        y[k] = y[j];
        x[j] = a;
        /* L360: */
        y[j] = b;
      }
    }
L400:
    ;
  }

  return 0;

} /* inv21_ */

/* Subroutine */
int cmplft_(real* x, real* y, integer* n, integer* d) {
  /* Subroutine */
  int s_stop();

  /* Local variables */
  static integer pmax;
  static integer psym;
  static logical error;
  static integer unsym[15];
  static integer factor[15], twogrp, sym[15];



  /*     COMPLEX FINITE DISCRETE FOURIER TRANSFORM */
  /*     TRANSFORMS ONE DIMENSION OF MULTI-DIMENSIONAL DATA */
  /*     MODIFIED BY L. F. TEN EYCK FROM A ONE-DIMENSIONAL VERSION WRITTEN
   */
  /*     BY G. T. SANDE, 1969. */

  /*     THIS PROGRAM CALCULATES THE TRANSFORM */
  /*               (X(T) + I*Y(T))*(COS(2*PI*T/N) - I*SIN(2*PI*T/N)) */

  /*     INDEXING -- THE ARRANGEMENT OF THE MULTI-DIMENSIONAL DATA IS */
  /*     SPECIFIED BY THE INTEGER ARRAY D, THE VALUES OF WHICH ARE USED AS
   */
  /*     CONTROL PARAMETERS IN DO LOOPS.  WHEN IT IS DESIRED TO COVER ALL */

  /*     ELEMENTS OF THE DATA FOR WHICH THE SUBSCRIPT BEING TRANSFORMED HAS
   */
  /*     THE VALUE I0, THE FOLLOWING IS USED. */

  /*               I1 = (I0 - 1)*D(2) + 1 */
  /*               DO 100 I2 = I1, D(1), D(3) */
  /*               I3 = I2 + D(4) - 1 */
  /*               DO 100 I = I2, I3, D(5) */
  /*                  . */
  /*                  . */
  /*           100 CONTINUE */

  /*     WITH THIS INDEXING IT IS POSSIBLE TO USE A NUMBER OF ARRANGEMENTS
   */
  /*     OF THE DATA, INCLUDING NORMAL FORTRAN COMPLEX NUMBERS (D(5) = 2) */

  /*     OR SEPARATE STORAGE OF REAL AND IMAGINARY PARTS. */


  /*     P MAX IS THE LARGEST PRIME FACTOR THAT WILL BE TOLERATED BY THIS */

  /*     PROGRAM. */
  /*     TWO GRP IS THE LARGEST POWER OF TWO THAT IS TREATED AS A SPECIAL */

  /*     CASE. */

  /* Parameter adjustments */
  --d;
  --y;
  --x;

  /* Function Body */
  pmax = 19;
  twogrp = 8;

  if (*n <= 1) {
    goto L100;
  }
  srfp_(n, &pmax, &twogrp, factor, sym, &psym, unsym, &error);
  if (error) {
    goto L200;
  }

  mdftkd_(n, factor, &d[1], &x[1], &y[1]);
  diprp_(n, sym, &psym, unsym, &d[1], &x[1], &y[1]);

L100:
  return 0;

L200:
  fprintf(stderr, "invalid number of points for cmpl_ft.  n = %ld\n", *n);
  return (-1);

} /* cmplft_ */

/* Subroutine */
int rsymft_(real* x, integer* n, integer* dim) {

  /* System generated locals */
  integer i__1, i__2, i__3, i__4, i__5;


  /* Local variables */
  static real twon;
  static integer twod2;
  static real a, b, c, d;
  static integer i, j, k, l, m;
  static real angle;
  static integer d1, d2, d3, d4, d5, j0, j1, k0, k1, k2, i0;
  static real twopi;
  static integer i1, i2, nover2, nover4;
  static real co;
  static integer ii, mj, mk, ml, mm;
  static real si;
  static integer nn;


  /*     REAL SYMMETRIC MULTIDIMENSIONAL FOURIER TRANSFORM */
  /*     N MUST BE A MULTIPLE OF 4.  THE TWO UNIQUE ELEMENTS ARE STORED AT
   */
  /*     X(1) AND X(N+1). */

  /*     DECIMATION IN FREQUENCY APPLIED TO A REAL SYMMETRIC SEQUENCE OF */
  /*     LENGTH 2N GIVES A REAL SYMMETRIC SEQUENCE OF LENGTH N, THE */
  /*     TRANSFORM OF WHICH GIVES THE EVEN NUMBERED FOURIER COEFFICIENTS, */

  /*     AND A HERMITIAN SYMMETRIC SEQUENCE OF LENGTH N, THE TRANSFORM OF */

  /*     WHICH GIVES THE ODD NUMBERED FOURIER COEFFICIENTS.  THE SUM OF */
  /*     THE TWO SEQUENCES IS A HERMITIAN SYMMETRIC SEQUENCE OF LENGTH N, */

  /*     WHICH MAY BE STORED IN N/2 COMPLEX LOCATIONS.  THE TRANSFORM OF */
  /*     THIS SEQUENCE IS N REAL NUMBERS REPRESENTING THE TERM BY TERM SUM
   */
  /*     OF THE EVEN AND ODD NUMBERED FOURIER COEFFICIENTS.  THIS SYMMETRIC
   */
  /*     SEQUENCE MAY BE SOLVED IF ANY OF THE FOURIER COEFFICIENTS ARE */
  /*     KNOWN.  FOR THIS PURPOSE X0, WHICH IS SIMPLY THE SUM OF THE */
  /*     ORIGINAL SEQUENCE, IS COMPUTED AND SAVED IN X(N+1). */


  /* Parameter adjustments */
  --dim;
  --x;

  /* Function Body */
  if (*n == 1) {
    goto L1300;
  }
  nover2 = *n / 2;
  nover4 = *n / 4;
  if (nover4 << 2 != *n) {
    goto L1400;
  }
  d1 = dim[1];
  d2 = dim[2];
  d3 = dim[3];
  d4 = dim[4] - 1;
  d5 = dim[5];
  twopi = (float)6.283185;
  twon = (real) (*n << 1);
  twod2 = d2 << 1;

  k0 = *n * d2 + 1;
  i__1 = d1;
  i__2 = d3;
  for (k1 = k0; k1 <= i__1; k1 += i__2) {
    k2 = k1 + d4;
    i__3 = k2;
    i__4 = d5;
    for (k = k1; k <= i__3; k += i__4) {
      x[k] /= (float)2.;
      /* L100: */
    }
  }

  i__4 = nover2;
  for (i = 2; i <= i__4; ++i) {
    angle = twopi * (real) (i - 1) / twon;
    co = (float)cos(angle);
    si = (float)sin(angle);
    k0 = (i - 1) * d2 + 1;
    j0 = (*n + 2 - (i << 1)) * d2;
    j1 = (*n + 1 - i) * d2;
    i__3 = d1;
    i__2 = d3;
    for (k1 = k0; k1 <= i__3; k1 += i__2) {
      k2 = k1 + d4;
      i__1 = k2;
      i__5 = d5;
      for (k = k1; k <= i__1; k += i__5) {
        l = k + j0;

        nn = k + j1;
        a = x[l] + x[k];
        b = x[l] - x[k];
        x[k] = a - b * co;
        x[l] = b * si;
        x[nn] += a;
        /* L200: */
      }
    }
    /* L300: */
  }

  if (nover4 == 1) {
    goto L600;
  }
  j0 = nover4 - 1;
  i__4 = j0;
  for (i = 1; i <= i__4; ++i) {
    k0 = (nover2 + i) * d2 + 1;
    j1 = (nover2 - (i << 1)) * d2;
    i__5 = d1;
    i__1 = d3;
    for (k1 = k0; k1 <= i__5; k1 += i__1) {
      k2 = k1 + d4;
      i__2 = k2;
      i__3 = d5;
      for (k = k1; k <= i__2; k += i__3) {
        l = k + j1;
        a = x[k];
        x[k] = x[l];
        x[l] = a;
        /* L400: */
      }
    }
    /* L500: */
  }

L600:
  j0 = nover2 * d2;
  j1 = *n * d2;
  i__4 = d1;
  i__3 = d3;
  for (k1 = 1; k1 <= i__4; k1 += i__3) {
    k2 = k1 + d4;
    i__2 = k2;
    i__1 = d5;
    for (k = k1; k <= i__2; k += i__1) {
      i = k + j0;
      l = k + j1;
      x[i] *= (float)2.;
      x[l] = x[k] + x[i] + x[l] * (float)2.;
      x[k] *= (float)2.;
      /* L700: */
    }
  }

  k = nover2 * d2 + 1;
  hermft(&x[1], &x[k], &nover2, &dim[1]);

  /*     SOLVE THE EQUATIONS FOR ALL OF THE SEQUENCES */

  i0 = 1 - d2;
  mk = nover2 * d2;
  mj = mk + d2;
  ml = *n * d2 + d2;
  mm = ml;
  i__1 = nover4;
  for (ii = 1; ii <= i__1; ++ii) {
    i0 += d2;
    mj -= twod2;
    ml -= twod2;
    mm -= d2;
    i__2 = d1;
    i__3 = d3;
    for (i1 = i0; i1 <= i__2; i1 += i__3) {
      i2 = i1 + d4;
      i__4 = i2;
      i__5 = d5;
      for (i = i1; i <= i__4; i += i__5) {
        j = i + mj;
        k = i + mk;
        l = i + ml;
        m = i + mm;
        a = x[i] - x[m];
        b = x[l] - a;
        c = x[k] - b;
        d = x[j] - c;
        x[i] = x[m];
        x[j] = a;
        x[k] = b;
        x[l] = c;
        x[m] = d;
        /* L800: */
      }
    }
  }

  /*     THE RESULTS ARE NOW IN A SCRAMBLED DIGIT REVERSED ORDER, I.E. */
  /*     X(1), X(5), X(9), ..., X(10), X(6), X(2), ..., X(3), X(7), X(11),
   */
  /*     ..., X(12), X(8), X(4).  THE FOLLOWING SECTION OF PROGRAM FOLLOWS
   */
  /*     THE PERMUTATION CYCLES AND DOES THE NECESSARY INTERCHANGES. */

  if (nover4 == 1) {
    goto L1300;
  }
  nn = *n - 2;
  i__5 = nn;
  for (i = 1; i <= i__5; ++i) {
    k = i;

L1000:
    k0 = k / 4;
    l = k - (k0 << 2);
    if (l != l / 2 << 1) {
      k0 = nover4 - 1 - k0;
    }
    k = k0 + l * nover4;
    if (k < i) {
      goto L1000;
    }
    if (k == i) {
      goto L1200;
    }

    k0 = i * d2 + 1;
    j0 = (k - i) * d2;
    i__4 = d1;
    i__3 = d3;
    for (k1 = k0; k1 <= i__4; k1 += i__3) {
      k2 = k1 + d4;
      i__2 = k2;
      i__1 = d5;
      for (k = k1; k <= i__2; k += i__1) {
        l = k + j0;
        a = x[k];
        x[k] = x[l];
        x[l] = a;
        /* L1100: */
      }
    }
L1200:
    ;
  }

L1300:
  return 0;

L1400:
  fprintf(stderr, "n not a multiple of 4 in r_sym_ft.  n = %ld\n", *n);
  return (-1);

} /* rsymft_ */

/* Subroutine */
int srfp_(integer* pts, integer* pmax, integer* twogrp, integer* factor, integer* sym, integer* psym, integer* unsym, logical* error) {

  /* System generated locals */
  integer i__1;


  /* Local variables */
  static integer nest, ptwo, f, j, n, p, q, r, jj, pp[14], qq[7];

  /* Fortran I/O blocks */


  /*     SYMMETRIZED REORDERING FACTORING PROGRAMME */



  /* Parameter adjustments */
  --unsym;
  --sym;
  --factor;

  /* Function Body */
  nest = 14;

  n = *pts;
  *psym = 1;
  f = 2;
  p = 0;
  q = 0;
L100:
  if (n <= 1) {
    goto L500;
  }
  i__1 = *pmax;
  for (j = f; j <= i__1; ++j) {
    if (n == n / j * j) {
      goto L300;
    }
    /* L200: */
  }
  goto L1400;
L300:
  if ((p << 1) + q >= nest) {
    goto L1600;
  }
  f = j;
  n /= f;
  if (n == n / f * f) {
    goto L400;
  }
  ++q;
  qq[q - 1] = f;
  goto L100;
L400:
  n /= f;
  ++p;
  pp[p - 1] = f;
  *psym *= f;
  goto L100;

L500:
  r = 1;
  if (q == 0) {
    r = 0;
  }
  if (p < 1) {
    goto L700;
  }
  i__1 = p;
  for (j = 1; j <= i__1; ++j) {
    jj = p + 1 - j;
    sym[j] = pp[jj - 1];
    factor[j] = pp[jj - 1];
    jj = p + q + j;
    factor[jj] = pp[j - 1];
    jj = p + r + j;
    sym[jj] = pp[j - 1];
    /* L600: */
  }
L700:
  if (q < 1) {
    goto L900;
  }
  i__1 = q;
  for (j = 1; j <= i__1; ++j) {
    jj = p + j;
    unsym[j] = qq[j - 1];
    factor[jj] = qq[j - 1];
    /* L800: */
  }
  /* Computing 2nd power */
  i__1 = *psym;
  sym[p + 1] = *pts / (i__1 * i__1);
L900:
  jj = (p << 1) + q;
  factor[jj + 1] = 0;
  ptwo = 1;
  j = 0;
L1000:
  ++j;
  if (factor[j] == 0) {
    goto L1200;
  }
  if (factor[j] != 2) {
    goto L1000;
  }
  ptwo <<= 1;
  factor[j] = 1;
  if (ptwo >= *twogrp) {
    goto L1100;
  }
  if (factor[j + 1] == 2) {
    goto L1000;
  }
L1100:
  factor[j] = ptwo;

  ptwo = 1;
  goto L1000;
L1200:
  if (p == 0) {
    r = 0;
  }
  jj = (p << 1) + r;
  sym[jj + 1] = 0;
  if (q <= 1) {
    q = 0;
  }
  unsym[q + 1] = 0;
  *error = false;

L1300:
  return 0;

L1400:
  fprintf(stderr, "largest factor exceeds %ld.  n = %ld\n", (long int)pmax, (long int)*pts);
  *error = true;
  goto L1300;

L1600:
  fprintf(stderr, "factor count exceeds %ld.  n = %ld\n", nest, *pts);
  *error = true;
  goto L1300;

} /* srfp_ */

/* Subroutine */
int diprp_(integer* pts, integer* sym, integer* psym, integer* unsym, integer* dim, real* x, real* y) {
  /* System generated locals */
  integer i__1, i__2, i__3, i__4, i__5, i__6, i__7, i__8, i__9, i__10,
          i__11, i__12, i__13, i__14, i__15, i__16, i__17, i__18, i__19,
          i__20, i__21, i__22, i__23, i__24, i__25, i__26, i__27, i__28,
          i__29, i__30, i__31;
  static integer equiv_25[14], equiv_26[14];

  /* Local variables */
  static integer mods, nest, size, test, mult, a, b, c, d, e, f, g, h, i, j,
                 k, l, m, n, p;
#define s (equiv_25)
  static real t;
#define u (equiv_26)
  static integer delta, p0, p1, p2, p3, p4, p5;
#define al (equiv_26)
#define bl (equiv_26 + 1)
  static integer dk;
#define cl (equiv_26 + 2)
#define dl (equiv_26 + 3)
#define el (equiv_26 + 4)
#define fl (equiv_26 + 5)
  static integer jj;
#define bs (equiv_25 + 1)
  static integer kk, lk;
#define cs (equiv_25 + 2)
#define ds (equiv_25 + 3)
#define es (equiv_25 + 4)
#define fs (equiv_25 + 5)
#define gs (equiv_25 + 6)
#define hs (equiv_25 + 7)
#define is (equiv_25 + 8)
#define js (equiv_25 + 9)
#define ks (equiv_25 + 10)
#define ls (equiv_25 + 11)
  static integer nt;
#define ms (equiv_25 + 12)
#define ns (equiv_25 + 13)
#define gl (equiv_26 + 6)
#define hl (equiv_26 + 7)
#define il (equiv_26 + 8)
#define jl (equiv_26 + 9)
#define kl (equiv_26 + 10)
#define ll (equiv_26 + 11)
#define ml (equiv_26 + 12)
#define nl (equiv_26 + 13)
  static logical onemod;
  static integer modulo[14], punsym, sep;

  /*     DOUBLE IN PLACE REORDERING PROGRAMME */




  /* Parameter adjustments */
  --y;
  --x;
  --dim;
  --unsym;
  --sym;

  /* Function Body */
  nest = 14;

  nt = dim[1];
  sep = dim[2];
  p2 = dim[3];
  size = dim[4] - 1;
  p4 = dim[5];
  if (sym[1] == 0) {
    goto L500;
  }
  i__1 = nest;
  for (j = 1; j <= i__1; ++j) {
    u[j - 1] = 1;
    s[j - 1] = 1;
    /* L100: */
  }
  n = *pts;
  i__1 = nest;
  for (j = 1; j <= i__1; ++j) {
    if (sym[j] == 0) {
      goto L300;
    }
    jj = nest + 1 - j;
    u[jj - 1] = n;
    s[jj - 1] = n / sym[j];
    n /= sym[j];
    /* L200: */
  }
L300:

  jj = 0;
  i__1 = *al;
  for (a = 1; a <= i__1; ++a) {
    i__2 = *bl;
    i__3 = *bs;
    for (b = a; b <= i__2; b += i__3) {
      i__4 = *cl;
      i__5 = *cs;
      for (c = b; c <= i__4; c += i__5) {
        i__6 = *dl;
        i__7 = *ds;
        for (d = c; d <= i__6; d += i__7) {
          i__8 = *el;
          i__9 = *es;
          for (e = d; e <= i__8; e += i__9) {

            i__10 = *fl;
            i__11 = *fs;
            for (f = e; f <= i__10; f +=
                   i__11) {
              i__12 = *gl;
              i__13 = *gs;
              for (g = f; g <= i__12;
                   g += i__13) {
                i__14 = *hl;
                i__15 = *hs;
                for (h = g; h <=
                     i__14; h += i__15) {
                  i__16 = *il;
                  i__17 = *is;
                  for (i = h; i <=
                       i__16; i += i__17) {
                    i__18 = *jl;
                    i__19 = *js;
                    for (j = i;
                         j <= i__18; j += i__19) {
                      i__20 = *kl;

                      i__21 = *ks;
                      for (k = j; k <= i__20; k +=
                             i__21) {
                        i__22 = *ll;
                        i__23 = *ls;
                        for (l = k; l <= i__22; l
                             += i__23) {
                          i__24 = *ml;
                          i__25 = *ms;
                          for (m = l; m <= i__24;
                               m += i__25) {
                            i__26 = *nl;
                            i__27 = *ns;
                            for (n = m; n <=
                                 i__26; n += i__27) {
                              ++jj;
                              if (jj >= n) {
                                goto L400;
                              }
                              delta = (n - jj) * sep;
                              p1 = (jj - 1) * sep + 1;
                              i__28 = nt;
                              i__29 = p2;
                              for (p0 = p1; p0 <= i__28; p0 += i__29) {
                                p3 = p0 + size;
                                i__30 = p3;
                                i__31 = p4;
                                for (p = p0; p <= i__30; p += i__31) {
                                  p5 = p + delta;
                                  t = x[p];
                                  x[p] = x[p5];
                                  x[p5] = t;
                                  t = y[p];
                                  y[p] = y[p5];
                                  y[p5] = t;
                                  /* L350: */
                                }
                              }
L400:
                              ;
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }

L500:

  if (unsym[1] == 0) {
    goto L1900;
  }
  /* Computing 2nd power */
  i__27 = *psym;
  punsym = *pts / (i__27 * i__27);
  mult = punsym / unsym[1];
  test = (unsym[1] * unsym[2] - 1) * mult * *psym;
  lk = mult;
  dk = mult;
  i__27 = nest;
  for (k = 2; k <= i__27; ++k) {
    if (unsym[k] == 0) {
      goto L700;
    }
    lk *= unsym[k - 1];
    dk /= unsym[k];
    u[k - 1] = (lk - dk) * *psym;
    mods = k;
    /* L600: */
  }
L700:
  onemod = mods < 3;
  if (onemod) {
    goto L900;
  }
  i__27 = mods;
  for (j = 3; j <= i__27; ++j) {
    jj = mods + 3 - j;
    modulo[jj - 1] = u[j - 1];
    /* L800: */
  }
L900:
  modulo[1] = u[1];
  *jl = (punsym - 3) * *psym;
  *ms = punsym * *psym;

  i__27 = *jl;
  i__26 = *psym;
  for (j = *psym; j <= i__27; j += i__26) {
    k = j;

L1000:
    k *= mult;
    if (onemod) {
      goto L1200;
    }
    i__25 = mods;
    for (i = 3; i <= i__25; ++i) {
      k -= k / modulo[i - 1] * modulo[i - 1];
      /* L1100: */
    }
L1200:
    if (k >= test) {
      goto L1300;
    }
    k -= k / modulo[1] * modulo[1];
    goto L1400;
L1300:
    k = k - k / modulo[1] * modulo[1] + modulo[1];
L1400:
    if (k < j) {
      goto L1000;
    }

    if (k == j) {
      goto L1700;
    }
    delta = (k - j) * sep;
    i__25 = *psym;
    for (l = 1; l <= i__25; ++l) {
      i__24 = *pts;
      i__23 = *ms;
      for (m = l; m <= i__24; m += i__23) {
        p1 = (m + j - 1) * sep + 1;
        i__22 = nt;
        i__21 = p2;
        for (p0 = p1; p0 <= i__22; p0 +=
               i__21) {
          p3 = p0 + size;
          i__20 = p3;
          i__19 = p4;
          for (jj = p0; jj <= i__20; jj +=
                 i__19) {
            kk = jj + delta;
            t = x[jj];
            x[jj] = x[kk];
            x[kk] = t;
            t = y[jj];
            y[jj] = y[kk];
            y[kk] = t;
            /* L1500: */
          }
        }
      }
      /* L1600: */
    }
L1700:
    /* L1800: */
    ;
  }

L1900:
  return 0;
} /* diprp_ */

#undef nl
#undef ml
#undef ll
#undef kl
#undef jl
#undef il
#undef hl
#undef gl
#undef ns
#undef ms
#undef ls
#undef ks
#undef js
#undef is
#undef hs

#undef gs
#undef fs
#undef es
#undef ds
#undef cs
#undef bs
#undef fl
#undef el
#undef dl
#undef cl
#undef bl
#undef al
#undef u
#undef s


/* Subroutine */
int mdftkd_(integer* n, integer* factor, integer* dim, real* x, real* y) {

  /* Local variables */
  static integer f, m, p, r, s;


  /*     MULTI-DIMENSIONAL COMPLEX FOURIER TRANSFORM KERNEL DRIVER */



  /* Parameter adjustments */
  --y;
  --x;
  --dim;
  --factor;

  /* Function Body */
  s = dim[2];
  f = 0;
  m = *n;
L100:
  ++f;
  p = factor[f];
  if (p == 0) {
    return 0;
  }
  m /= p;
  r = m * s;
  if (p > 8) {
    goto L700;
  }
  switch ((int)p) {
    case 1:
      goto L100;
    case 2:
      goto L200;
    case 3:
      goto L300;
    case 4:
      goto L400;
    case 5:
      goto L500;
    case 6:
      goto L800;
    case 7:
      goto L700;
    case 8:
      goto L600;
  }
  goto L800;

L200:
  r2cftk_(n, &m, &x[1], &y[1], &x[r + 1], &y[r + 1], &dim[1]);
  goto L100;

L300:
  r3cftk_(n, &m, &x[1], &y[1], &x[r + 1], &y[r + 1], &x[(r << 1) + 1], &y[(
                                                                            r << 1) + 1], &dim[1]);
  goto L100;

L400:
  r4cftk_(n, &m, &x[1], &y[1], &x[r + 1], &y[r + 1], &x[(r << 1) + 1], &y[(
                                                                            r << 1) + 1], &x[r * 3 + 1], &y[r * 3 + 1], &dim[1]);
  goto L100;

L500:
  r5cftk_(n, &m, &x[1], &y[1], &x[r + 1], &y[r + 1], &x[(r << 1) + 1], &y[(
                                                                            r << 1) + 1], &x[r * 3 + 1], &y[r * 3 + 1], &x[(r << 2) + 1], &y[(
                                                                                                                                               r << 2) + 1], &dim[1]);
  goto L100;

L600:
  r8cftk_(n, &m, &x[1], &y[1], &x[r + 1], &y[r + 1], &x[(r << 1) + 1], &y[(
                                                                            r << 1) + 1], &x[r * 3 + 1], &y[r * 3 + 1], &x[(r << 2) + 1], &y[(
                                                                                                                                               r << 2) + 1], &x[r * 5 + 1], &y[r * 5 + 1], &x[r * 6 + 1], &y[r *
                                                                                                                                                                                                             6 + 1], &x[r * 7 + 1], &y[r * 7 + 1], &dim[1]);
  goto L100;

L700:
  rpcftk_(n, &m, &p, &r, &x[1], &y[1], &dim[1]);
  goto L100;

L800:
  fprintf(stderr, "transfer error detected in mdftkd p = %ld\n", p);
  return 0;
} /* mdftkd_ */

/* Subroutine */
int r2cftk_(integer* n, integer* m, real* x0, real* y0, real* x1, real* y1, integer* dim) {
  /* Initialized data */

  static real twopi = (float)6.283185;

  /* System generated locals */
  integer i__1, i__2, i__3, i__4, i__5, i__6, i__7;

  /* Local variables */
  static logical fold;
  static integer size;
  static logical zero;
  static real c;
  static integer j, k, l;
  static real s, angle;
  static integer k0, k1, k2, l1, m2, mover2, kk;
  static real is, iu;
  static integer ns, nt;
  static real rs, ru, fm2;
  static integer mm2, sep;
  static real fjm1;

  /*     RADIX 2 MULTI-DIMENSIONAL COMPLEX FOURIER TRANSFORM KERNEL */


  /* Parameter adjustments */
  --dim;
  --y1;
  --x1;
  --y0;
  --x0;

  /* Function Body */

  nt = dim[1];
  sep = dim[2];
  l1 = dim[3];
  size = dim[4] - 1;
  k2 = dim[5];
  ns = *n * sep;
  m2 = *m << 1;
  fm2 = (real) m2;
  mover2 = *m / 2 + 1;
  mm2 = sep * m2;

  fjm1 = (float)-1.;
  i__1 = mover2;
  for (j = 1; j <= i__1; ++j) {
    fold = j > 1 && j << 1 < *m + 2;
    k0 = (j - 1) * sep + 1;
    fjm1 += (float)1.;
    angle = twopi * fjm1 / fm2;
    zero = angle == (float)0.;
    if (zero) {
      goto L200;
    }
    c = (float)cos(angle);
    s = (float)sin(angle);
    goto L200;
L100:
    fold = false;
    k0 = (*m + 1 - j) * sep + 1;
    c = -c;
L200:

    i__2 = ns;
    i__3 = mm2;
    for (kk = k0; kk <= i__2; kk += i__3) {
      i__4 = nt;
      i__5 = l1;
      for (l = kk; l <= i__4; l += i__5) {
        k1 = l + size;
        i__6 = k1;
        i__7 = k2;
        for (k = l; k <= i__6; k += i__7) {
          rs = x0[k] + x1[k];
          is = y0[k] + y1[k];
          ru = x0[k] - x1[k];
          iu = y0[k] - y1[k];
          x0[k] = rs;
          y0[k] = is;
          if (zero) {
            goto L300;
          }
          x1[k] = ru * c + iu * s;
          y1[k] = iu * c - ru * s;
          goto L400;
L300:
          x1[k] = ru;
          y1[k] = iu;
L400:
          /* L420: */
          ;
        }
        /* L440: */
      }
      /* L500: */
    }
    if (fold) {
      goto L100;
    }
    /* L600: */
  }

  return 0;
} /* r2cftk_ */

/* Subroutine */
int r3cftk_(integer* n, integer* m, real* x0, real* y0, real* x1, real* y1, real* x2, real* y2, integer* dim) {
  /* Initialized data */

  static real twopi = (float)6.283185;
  static real a = (float)-.5;
  static real b = (float).8660254;

  /* System generated locals */
  integer i__1, i__2, i__3, i__4, i__5, i__6, i__7;

  /* Local variables */
  static logical fold;
  static integer size;
  static logical zero;
  static integer j, k, l;
  static real t, angle, c1, c2, i0, i1;
  static integer k0, k1, k2, l1;
  static real i2;
  static integer m3;
  static real r0, s1, s2, r1, r2;
  static integer mover2;
  static real ia, ib, ra, rb;
  static integer kk;
  static real is;
  static integer ns, nt;
  static real rs, fm3;
  static integer mm3, sep;
  static real fjm1;

  /*     RADIX 3 MULTI-DIMENSIONAL COMPLEX FOURIER TRANSFORM KERNEL */


  /* Parameter adjustments */
  --dim;
  --y2;
  --x2;
  --y1;
  --x1;
  --y0;
  --x0;

  /* Function Body */

  nt = dim[1];
  sep = dim[2];
  l1 = dim[3];
  size = dim[4] - 1;
  k2 = dim[5];
  ns = *n * sep;
  m3 = *m * 3;
  fm3 = (real) m3;
  mm3 = sep * m3;
  mover2 = *m / 2 + 1;

  fjm1 = (float)-1.;
  i__1 = mover2;
  for (j = 1; j <= i__1; ++j) {
    fold = j > 1 && j << 1 < *m + 2;
    k0 = (j - 1) * sep + 1;
    fjm1 += (float)1.;
    angle = twopi * fjm1 / fm3;
    zero = angle == (float)0.;
    if (zero) {
      goto L200;
    }
    c1 = (float)cos(angle);
    s1 = (float)sin(angle);
    c2 = c1 * c1 - s1 * s1;
    s2 = s1 * c1 + c1 * s1;
    goto L200;
L100:
    fold = false;
    k0 = (*m + 1 - j) * sep + 1;
    t = c1 * a + s1 * b;
    s1 = c1 * b - s1 * a;
    c1 = t;
    t = c2 * a - s2 * b;
    s2 = -c2 * b - s2 * a;
    c2 = t;
L200:

    i__2 = ns;
    i__3 = mm3;
    for (kk = k0; kk <= i__2; kk += i__3) {
      i__4 = nt;
      i__5 = l1;
      for (l = kk; l <= i__4; l += i__5) {
        k1 = l + size;
        i__6 = k1;
        i__7 = k2;
        for (k = l; k <= i__6; k += i__7) {
          r0 = x0[k];
          i0 = y0[k];
          rs = x1[k] + x2[k];
          is = y1[k] + y2[k];
          x0[k] = r0 + rs;
          y0[k] = i0 + is;
          ra = r0 + rs * a;
          ia = i0 + is * a;
          rb = (x1[k] - x2[k]) * b;
          ib = (y1[k] - y2[k]) * b;
          if (zero) {
            goto L300;
          }
          r1 = ra + ib;
          i1 = ia - rb;
          r2 = ra - ib;
          i2 = ia + rb;
          x1[k] = r1 * c1 + i1 * s1;
          y1[k] = i1 * c1 - r1 * s1;
          x2[k] = r2 * c2 + i2 * s2;
          y2[k] = i2 * c2 - r2 * s2;
          goto L400;
L300:
          x1[k] = ra + ib;
          y1[k] = ia - rb;
          x2[k] = ra - ib;
          y2[k] = ia + rb;
L400:
          /* L420: */
          ;
        }
        /* L440: */
      }
      /* L500: */
    }
    if (fold) {
      goto L100;
    }
    /* L600: */
  }

  return 0;
} /* r3cftk_ */

/* Subroutine */
int r4cftk_(integer* n, integer* m, real* x0, real* y0, real* x1, real* y1, real* x2, real* y2, real* x3, real* y3, integer* dim) {
  /* Initialized data */

  static real twopi = (float)6.283185;

  /* System generated locals */
  integer i__1, i__2, i__3, i__4, i__5, i__6, i__7;

  /* Local variables */
  static logical fold;
  static integer size;
  static logical zero;
  static integer j, k, l;
  static real t, angle, c1, c2, c3, i1;
  static integer k0, k1, k2, l1;
  static real i2, i3;
  static integer m4;
  static real r1, s1, s2, s3, r2, r3;
  static integer mover2, kk, ns, nt;
  static real fm4, is0, is1;
  static integer mm4;
  static real iu0, iu1, rs0, rs1, ru0, ru1;
  static integer sep;
  static real fjm1;

  /*     RADIX 4 MULTI-DIMENSIONAL COMPLEX FOURIER TRANSFORM KERNEL */


  /* Parameter adjustments */
  --dim;
  --y3;
  --x3;
  --y2;
  --x2;
  --y1;
  --x1;
  --y0;
  --x0;

  /* Function Body */

  nt = dim[1];
  sep = dim[2];
  l1 = dim[3];
  size = dim[4] - 1;
  k2 = dim[5];
  ns = *n * sep;
  m4 = *m << 2;
  fm4 = (real) m4;
  mm4 = sep * m4;
  mover2 = *m / 2 + 1;

  fjm1 = (float)-1.;
  i__1 = mover2;
  for (j = 1; j <= i__1; ++j) {
    fold = j > 1 && j << 1 < *m + 2;
    k0 = (j - 1) * sep + 1;
    fjm1 += (float)1.;
    angle = twopi * fjm1 / fm4;
    zero = angle == (float)0.;
    if (zero) {
      goto L200;
    }
    c1 = (float)cos(angle);
    s1 = (float)sin(angle);
    c2 = c1 * c1 - s1 * s1;
    s2 = s1 * c1 + c1 * s1;
    c3 = c2 * c1 - s2 * s1;
    s3 = s2 * c1 + c2 * s1;
    goto L200;
L100:
    fold = false;
    k0 = (*m + 1 - j) * sep + 1;
    t = c1;
    c1 = s1;
    s1 = t;
    c2 = -c2;
    t = c3;
    c3 = -s3;
    s3 = -t;
L200:

    i__2 = ns;
    i__3 = mm4;
    for (kk = k0; kk <= i__2; kk += i__3) {
      i__4 = nt;
      i__5 = l1;
      for (l = kk; l <= i__4; l += i__5) {
        k1 = l + size;
        i__6 = k1;
        i__7 = k2;
        for (k = l; k <= i__6; k += i__7) {
          rs0 = x0[k] + x2[k];
          is0 = y0[k] + y2[k];
          ru0 = x0[k] - x2[k];
          iu0 = y0[k] - y2[k];
          rs1 = x1[k] + x3[k];
          is1 = y1[k] + y3[k];
          ru1 = x1[k] - x3[k];
          iu1 = y1[k] - y3[k];
          x0[k] = rs0 + rs1;
          y0[k] = is0 + is1;
          if (zero) {
            goto L300;
          }
          r1 = ru0 + iu1;
          i1 = iu0 - ru1;
          r2 = rs0 - rs1;
          i2 = is0 - is1;
          r3 = ru0 - iu1;
          i3 = iu0 + ru1;
          x2[k] = r1 * c1 + i1 * s1;
          y2[k] = i1 * c1 - r1 * s1;
          x1[k] = r2 * c2 + i2 * s2;
          y1[k] = i2 * c2 - r2 * s2;
          x3[k] = r3 * c3 + i3 * s3;
          y3[k] = i3 * c3 - r3 * s3;
          goto L400;
L300:
          x2[k] = ru0 + iu1;
          y2[k] = iu0 - ru1;
          x1[k] = rs0 - rs1;
          y1[k] = is0 - is1;
          x3[k] = ru0 - iu1;
          y3[k] = iu0 + ru1;
L400:
          /* L420: */
          ;
        }
        /* L440: */
      }
      /* L500: */
    }
    if (fold) {
      goto L100;
    }
    /* L600: */
  }

  return 0;
} /* r4cftk_ */

/* Subroutine */
int r5cftk_(integer* n, integer* m, real* x0, real* y0, real* x1, real* y1, real* x2, real* y2, real* x3, real* y3, real* x4, real* y4,
            integer* dim) {
  /* Initialized data */

  static real twopi = (float)6.283185;
  static real a1 = (float).30901699;
  static real b1 = (float).95105652;
  static real a2 = (float)-.80901699;
  static real b2 = (float).58778525;

  /* System generated locals */
  integer i__1, i__2, i__3, i__4, i__5, i__6, i__7;

  /* Local variables */
  static logical fold;
  static integer size;
  static logical zero;
  static integer j, k, l;
  static real t, angle, c1, c2, c3, c4, i0;
  static integer k0, k1, k2, l1;
  static real i1, i2, i3;
  static integer m5;
  static real s1, s2, s3, s4, r0, r1, r2, r3, r4, i4;
  static integer mover2, kk, ns, nt;
  static real ia1, ia2, ib1, ib2, ra1, ra2, rb1, rb2, fm5, is1, is2;
  static integer mm5;
  static real iu1, iu2, rs1, rs2, ru1, ru2;
  static integer sep;
  static real fjm1;

  /*     RADIX 5 MULTI-DIMENSIONAL COMPLEX FOURIER TRANSFORM KERNEL */


  /* Parameter adjustments */
  --dim;
  --y4;
  --x4;
  --y3;
  --x3;
  --y2;
  --x2;
  --y1;
  --x1;
  --y0;
  --x0;

  /* Function Body */

  nt = dim[1];
  sep = dim[2];
  l1 = dim[3];
  size = dim[4] - 1;
  k2 = dim[5];
  ns = *n * sep;
  m5 = *m * 5;
  fm5 = (real) m5;
  mm5 = sep * m5;
  mover2 = *m / 2 + 1;

  fjm1 = (float)-1.;
  i__1 = mover2;
  for (j = 1; j <= i__1; ++j) {
    fold = j > 1 && j << 1 < *m + 2;
    k0 = (j - 1) * sep + 1;
    fjm1 += (float)1.;
    angle = twopi * fjm1 / fm5;
    zero = angle == (float)0.;
    if (zero) {
      goto L200;
    }
    c1 = (float)cos(angle);
    s1 = (float)sin(angle);
    c2 = c1 * c1 - s1 * s1;
    s2 = s1 * c1 + c1 * s1;
    c3 = c2 * c1 - s2 * s1;
    s3 = s2 * c1 + c2 * s1;
    c4 = c2 * c2 - s2 * s2;
    s4 = s2 * c2 + c2 * s2;
    goto L200;
L100:
    fold = false;
    k0 = (*m + 1 - j) * sep + 1;
    t = c1 * a1 + s1 * b1;
    s1 = c1 * b1 - s1 * a1;
    c1 = t;
    t = c2 * a2 + s2 * b2;
    s2 = c2 * b2 - s2 * a2;
    c2 = t;
    t = c3 * a2 - s3 * b2;
    s3 = -c3 * b2 - s3 * a2;
    c3 = t;
    t = c4 * a1 - s4 * b1;
    s4 = -c4 * b1 - s4 * a1;
    c4 = t;
L200:

    i__2 = ns;
    i__3 = mm5;
    for (kk = k0; kk <= i__2; kk += i__3) {
      i__4 = nt;
      i__5 = l1;
      for (l = kk; l <= i__4; l += i__5) {
        k1 = l + size;
        i__6 = k1;
        i__7 = k2;
        for (k = l; k <= i__6; k += i__7) {
          r0 = x0[k];
          i0 = y0[k];
          rs1 = x1[k] + x4[k];
          is1 = y1[k] + y4[k];
          ru1 = x1[k] - x4[k];
          iu1 = y1[k] - y4[k];
          rs2 = x2[k] + x3[k];
          is2 = y2[k] + y3[k];
          ru2 = x2[k] - x3[k];
          iu2 = y2[k] - y3[k];
          x0[k] = r0 + rs1 + rs2;
          y0[k] = i0 + is1 + is2;
          ra1 = r0 + rs1 * a1 + rs2 * a2;
          ia1 = i0 + is1 * a1 + is2 * a2;
          ra2 = r0 + rs1 * a2 + rs2 * a1;
          ia2 = i0 + is1 * a2 + is2 * a1;
          rb1 = ru1 * b1 + ru2 * b2;
          ib1 = iu1 * b1 + iu2 * b2;
          rb2 = ru1 * b2 - ru2 * b1;
          ib2 = iu1 * b2 - iu2 * b1;
          if (zero) {
            goto L300;
          }
          r1 = ra1 + ib1;
          i1 = ia1 - rb1;
          r2 = ra2 + ib2;
          i2 = ia2 - rb2;
          r3 = ra2 - ib2;
          i3 = ia2 + rb2;
          r4 = ra1 - ib1;
          i4 = ia1 + rb1;
          x1[k] = r1 * c1 + i1 * s1;
          y1[k] = i1 * c1 - r1 * s1;
          x2[k] = r2 * c2 + i2 * s2;
          y2[k] = i2 * c2 - r2 * s2;
          x3[k] = r3 * c3 + i3 * s3;
          y3[k] = i3 * c3 - r3 * s3;
          x4[k] = r4 * c4 + i4 * s4;
          y4[k] = i4 * c4 - r4 * s4;
          goto L400;
L300:
          x1[k] = ra1 + ib1;
          y1[k] = ia1 - rb1;
          x2[k] = ra2 + ib2;
          y2[k] = ia2 - rb2;
          x3[k] = ra2 - ib2;
          y3[k] = ia2 + rb2;
          x4[k] = ra1 - ib1;
          y4[k] = ia1 + rb1;
L400:
          /* L420: */
          ;
        }
        /* L440: */
      }
      /* L500: */
    }
    if (fold) {
      goto L100;
    }
    /* L600: */
  }

  return 0;
} /* r5cftk_ */

/* Subroutine */
int r8cftk_(integer* n, integer* m, real* x0, real* y0, real* x1, real* y1, real* x2, real* y2, real* x3, real* y3, real* x4, real* y4, real* x5,
            real* y5, real* x6, real* y6, real* x7, real* y7, integer* dim) {
  /* Initialized data */

  static real twopi = (float)6.283185;
  static real e = (float).70710678;

  /* System generated locals */
  integer i__1, i__2, i__3, i__4, i__5, i__6, i__7;

  /* Local variables */
  static logical fold;
  static integer size;
  static logical zero;
  static integer j, k, l;
  static real t, angle, c1, c2, c3, c4, c5, c6, c7;
  static integer k0, k1, k2, l1;
  static real i1, i2, i3, r1, s1;
  static integer m8;
  static real s2, s3, s4, s5, s6, s7, r2, r3, r4, r5, r6, r7, i4, i5, i6,
              i7;
  static integer mover2, kk, ns, nt;
  static real fm8, is0, is1, is2, is3, iu0, iu1;
  static integer mm8;
  static real iu2, iu3, rs0, rs1, rs2, rs3, ru0, ru1, ru2, ru3;
  static integer sep;
  static real fjm1, iss0, iss1, isu0, isu1, ius0, ius1, iuu0, iuu1, rss0,
              rss1, rsu0, rsu1, rus0, rus1, ruu0, ruu1;

  /*     RADIX 8 MULTI-DIMENSIONAL COMPLEX FOURIER TRANSFORM KERNEL */


  /* Parameter adjustments */
  --dim;
  --y7;
  --x7;
  --y6;
  --x6;
  --y5;
  --x5;
  --y4;
  --x4;
  --y3;
  --x3;
  --y2;
  --x2;
  --y1;
  --x1;
  --y0;
  --x0;

  /* Function Body */

  nt = dim[1];
  sep = dim[2];
  l1 = dim[3];
  size = dim[4] - 1;
  k2 = dim[5];
  ns = *n * sep;
  m8 = *m << 3;
  fm8 = (real) m8;
  mm8 = sep * m8;
  mover2 = *m / 2 + 1;

  fjm1 = (float)-1.;
  i__1 = mover2;
  for (j = 1; j <= i__1; ++j) {
    fold = j > 1 && j << 1 < *m + 2;
    k0 = (j - 1) * sep + 1;
    fjm1 += (float)1.;
    angle = twopi * fjm1 / fm8;
    zero = angle == (float)0.;
    if (zero) {
      goto L200;
    }
    c1 = (float)cos(angle);
    s1 = (float)sin(angle);
    c2 = c1 * c1 - s1 * s1;
    s2 = s1 * c1 + c1 * s1;
    c3 = c2 * c1 - s2 * s1;
    s3 = s2 * c1 + c2 * s1;
    c4 = c2 * c2 - s2 * s2;
    s4 = s2 * c2 + c2 * s2;
    c5 = c4 * c1 - s4 * s1;
    s5 = s4 * c1 + c4 * s1;
    c6 = c4 * c2 - s4 * s2;
    s6 = s4 * c2 + c4 * s2;
    c7 = c4 * c3 - s4 * s3;
    s7 = s4 * c3 + c4 * s3;
    goto L200;
L100:
    fold = false;
    k0 = (*m + 1 - j) * sep + 1;
    t = (c1 + s1) * e;
    s1 = (c1 - s1) * e;
    c1 = t;

    t = s2;
    s2 = c2;
    c2 = t;
    t = (-c3 + s3) * e;
    s3 = (c3 + s3) * e;
    c3 = t;
    c4 = -c4;
    t = -(c5 + s5) * e;
    s5 = (-c5 + s5) * e;
    c5 = t;
    t = -s6;
    s6 = -c6;
    c6 = t;
    t = (c7 - s7) * e;
    s7 = -(c7 + s7) * e;
    c7 = t;
L200:

    i__2 = ns;
    i__3 = mm8;
    for (kk = k0; kk <= i__2; kk += i__3) {
      i__4 = nt;
      i__5 = l1;
      for (l = kk; l <= i__4; l += i__5) {
        k1 = l + size;
        i__6 = k1;
        i__7 = k2;
        for (k = l; k <= i__6; k += i__7) {
          rs0 = x0[k] + x4[k];
          is0 = y0[k] + y4[k];
          ru0 = x0[k] - x4[k];
          iu0 = y0[k] - y4[k];
          rs1 = x1[k] + x5[k];
          is1 = y1[k] + y5[k];
          ru1 = x1[k] - x5[k];
          iu1 = y1[k] - y5[k];
          rs2 = x2[k] + x6[k];
          is2 = y2[k] + y6[k];
          ru2 = x2[k] - x6[k];
          iu2 = y2[k] - y6[k];
          rs3 = x3[k] + x7[k];
          is3 = y3[k] + y7[k];
          ru3 = x3[k] - x7[k];
          iu3 = y3[k] - y7[k];
          rss0 = rs0 + rs2;
          iss0 = is0 + is2;
          rsu0 = rs0 - rs2;
          isu0 = is0 - is2;
          rss1 = rs1 + rs3;
          iss1 = is1 + is3;
          rsu1 = rs1 - rs3;
          isu1 = is1 - is3;
          rus0 = ru0 - iu2;
          ius0 = iu0 + ru2;
          ruu0 = ru0 + iu2;
          iuu0 = iu0 - ru2;
          rus1 = ru1 - iu3;
          ius1 = iu1 + ru3;
          ruu1 = ru1 + iu3;
          iuu1 = iu1 - ru3;
          t = (rus1 + ius1) * e;
          ius1 = (ius1 - rus1) * e;
          rus1 = t;
          t = (ruu1 + iuu1) * e;
          iuu1 = (iuu1 - ruu1) * e;
          ruu1 = t;
          x0[k] = rss0 + rss1;
          y0[k] = iss0 + iss1;
          if (zero) {
            goto L300;
          }
          r1 = ruu0 + ruu1;
          i1 = iuu0 + iuu1;
          r2 = rsu0 + isu1;
          i2 = isu0 - rsu1;
          r3 = rus0 + ius1;
          i3 = ius0 - rus1;
          r4 = rss0 - rss1;
          i4 = iss0 - iss1;
          r5 = ruu0 - ruu1;
          i5 = iuu0 - iuu1;
          r6 = rsu0 - isu1;
          i6 = isu0 + rsu1;
          r7 = rus0 - ius1;
          i7 = ius0 + rus1;
          x4[k] = r1 * c1 + i1 * s1;
          y4[k] = i1 * c1 - r1 * s1;
          x2[k] = r2 * c2 + i2 * s2;
          y2[k] = i2 * c2 - r2 * s2;
          x6[k] = r3 * c3 + i3 * s3;
          y6[k] = i3 * c3 - r3 * s3;
          x1[k] = r4 * c4 + i4 * s4;
          y1[k] = i4 * c4 - r4 * s4;
          x5[k] = r5 * c5 + i5 * s5;
          y5[k] = i5 * c5 - r5 * s5;
          x3[k] = r6 * c6 + i6 * s6;
          y3[k] = i6 * c6 - r6 * s6;
          x7[k] = r7 * c7 + i7 * s7;
          y7[k] = i7 * c7 - r7 * s7;
          goto L400;
L300:
          x4[k] = ruu0 + ruu1;
          y4[k] = iuu0 + iuu1;
          x2[k] = rsu0 + isu1;
          y2[k] = isu0 - rsu1;
          x6[k] = rus0 + ius1;
          y6[k] = ius0 - rus1;
          x1[k] = rss0 - rss1;
          y1[k] = iss0 - iss1;
          x5[k] = ruu0 - ruu1;
          y5[k] = iuu0 - iuu1;
          x3[k] = rsu0 - isu1;
          y3[k] = isu0 + rsu1;
          x7[k] = rus0 - ius1;
          y7[k] = ius0 + rus1;
L400:
          /* L420: */
          ;
        }
        /* L440: */
      }
      /* L500: */
    }
    if (fold) {
      goto L100;
    }
    /* L600: */
  }

  return 0;
} /* r8cftk_ */

/* Subroutine */
int rpcftk_(integer* n, integer* m, integer* p, integer* r, real* x, real* y, integer* dim) {
  /* Initialized data */

  static real twopi = (float)6.283185;

  /* System generated locals */
  integer x_dim1, x_offset, y_dim1, y_offset, i__1, i__2, i__3, i__4, i__5,
          i__6, i__7, i__8, i__9;

  /* Local variables */
  static logical fold;
  static integer size;
  static logical zero;
  static real a[18], b[18], c[18];
  static integer j, k, l;
  static real s[18], t;
  static integer u;
  static real angle;
  static integer v, k0, k1, k2, l1;
  static real aa[81] /* was [9][9] */, bb[81] /* was [9][9] */;
  static integer mover2;
  static real ia[9], ib[9], ra[9];
  static integer jj;
  static real rb[9], fp;
  static integer kk;
  static real fu, is;
  static integer mp;
  static real iu;
  static integer pm, pp, ns, nt;
  static real rs, ru, xt, yt, fmp;
  static integer sep, mmp;
  static real fjm1;

  /*     RADIX PRIME MULTI-DIMENSIONAL COMPLEX FOURIER TRANSFORM KERNEL */



  /* Parameter adjustments */
  --dim;
  y_dim1 = *r;
  y_offset = y_dim1 + 1;
  y -= y_offset;
  x_dim1 = *r;
  x_offset = x_dim1 + 1;
  x -= x_offset;

  /* Function Body */

  nt = dim[1];
  sep = dim[2];
  l1 = dim[3];
  size = dim[4] - 1;
  k2 = dim[5];
  ns = *n * sep;
  mover2 = *m / 2 + 1;
  mp = *m * *p;
  fmp = (real) mp;
  mmp = sep * mp;
  pp = *p / 2;
  pm = *p - 1;
  fp = (real) (*p);
  fu = (float)0.;
  i__1 = pp;
  for (u = 1; u <= i__1; ++u) {
    fu += (float)1.;
    angle = twopi * fu / fp;
    jj = *p - u;
    a[u - 1] = (float)cos(angle);
    b[u - 1] = (float)sin(angle);
    a[jj - 1] = a[u - 1];
    b[jj - 1] = -b[u - 1];
    /* L100: */
  }
  i__1 = pp;
  for (u = 1; u <= i__1; ++u) {
    i__2 = pp;
    for (v = 1; v <= i__2; ++v) {
      jj = u * v - u * v / *p * *p;
      aa[v + u * 9 - 10] = a[jj - 1];
      bb[v + u * 9 - 10] = b[jj - 1];
      /* L200: */
    }
    /* L300: */
  }


  fjm1 = (float)-1.;
  i__1 = mover2;
  for (j = 1; j <= i__1; ++j) {
    fold = j > 1 && j << 1 < *m + 2;
    k0 = (j - 1) * sep + 1;
    fjm1 += (float)1.;
    angle = twopi * fjm1 / fmp;
    zero = angle == (float)0.;
    if (zero) {
      goto L700;
    }
    c[0] = (float)cos(angle);
    s[0] = (float)sin(angle);
    i__2 = pm;
    for (u = 2; u <= i__2; ++u) {
      c[u - 1] = c[u - 2] * c[0] - s[u - 2] * s[0];
      s[u - 1] = s[u - 2] * c[0] + c[u - 2] * s[0];
      /* L400: */
    }
    goto L700;
L500:
    fold = false;
    k0 = (*m + 1 - j) * sep + 1;
    i__2 = pm;
    for (u = 1; u <= i__2; ++u) {
      t = c[u - 1] * a[u - 1] + s[u - 1] * b[u - 1];
      s[u - 1] = -s[u - 1] * a[u - 1] + c[u - 1] * b[u - 1];
      c[u - 1] = t;
      /* L600: */
    }
L700:

    i__2 = ns;
    i__3 = mmp;
    for (kk = k0; kk <= i__2; kk += i__3) {
      i__4 = nt;
      i__5 = l1;
      for (l = kk; l <= i__4; l += i__5) {
        k1 = l + size;
        i__6 = k1;
        i__7 = k2;
        for (k = l; k <= i__6; k += i__7) {
          xt = x[k + x_dim1];
          yt = y[k + y_dim1];
          rs = x[k + (x_dim1 << 1)] + x[k + *p * x_dim1];
          is = y[k + (y_dim1 << 1)] + y[k + *p * y_dim1];
          ru = x[k + (x_dim1 << 1)] - x[k + *p * x_dim1];
          iu = y[k + (y_dim1 << 1)] - y[k + *p * y_dim1];
          i__8 = pp;
          for (u = 1; u <= i__8; ++u) {
            ra[u - 1] = xt + rs * aa[u - 1];
            ia[u - 1] = yt + is * aa[u - 1];
            rb[u - 1] = ru * bb[u - 1];
            ib[u - 1] = iu * bb[u - 1];
            /* L800: */
          }
          xt += rs;
          yt += is;
          i__8 = pp;
          for (u = 2; u <= i__8; ++u) {
            jj = *p - u;
            rs = x[k + (u + 1) * x_dim1] + x[k + (jj + 1) *
                                             x_dim1];
            is = y[k + (u + 1) * y_dim1] + y[k + (jj + 1) *
                                             y_dim1];
            ru = x[k + (u + 1) * x_dim1] - x[k + (jj + 1) *
                                             x_dim1];
            iu = y[k + (u + 1) * y_dim1] - y[k + (jj + 1) *
                                             y_dim1];
            xt += rs;
            yt += is;
            i__9 = pp;
            for (v = 1; v <= i__9; ++v) {
              ra[v - 1] += rs * aa[v + u * 9 - 10];
              ia[v - 1] += is * aa[v + u * 9 - 10];
              rb[v - 1] += ru * bb[v + u * 9 - 10];
              ib[v - 1] += iu * bb[v + u * 9 - 10];
              /* L900: */
            }
            /* L1000: */
          }
          x[k + x_dim1] = xt;
          y[k + y_dim1] = yt;
          i__8 = pp;
          for (u = 1; u <= i__8; ++u) {
            jj = *p - u;
            if (zero) {
              goto L1100;
            }
            xt = ra[u - 1] + ib[u - 1];
            yt = ia[u - 1] - rb[u - 1];
            x[k + (u + 1) * x_dim1] = xt * c[u - 1] + yt * s[u -
                                                             1];
            y[k + (u + 1) * y_dim1] = yt * c[u - 1] - xt * s[u -
                                                             1];
            xt = ra[u - 1] - ib[u - 1];
            yt = ia[u - 1] + rb[u - 1];
            x[k + (jj + 1) * x_dim1] = xt * c[jj - 1] + yt * s[jj
                                                               - 1];
            y[k + (jj + 1) * y_dim1] = yt * c[jj - 1] - xt * s[jj
                                                               - 1];
            goto L1200;
L1100:
            x[k + (u + 1) * x_dim1] = ra[u - 1] + ib[u - 1];
            y[k + (u + 1) * y_dim1] = ia[u - 1] - rb[u - 1];
            x[k + (jj + 1) * x_dim1] = ra[u - 1] - ib[u - 1];
            y[k + (jj + 1) * y_dim1] = ia[u - 1] + rb[u - 1];
L1200:
            /* L1300: */
            ;
          }
          /* L1320: */
        }
        /* L1340: */
      }
      /* L1400: */
    }
    if (fold) {
      goto L500;
    }
    /* L1500: */
  }

  return 0;
} /* rpcftk_ */

