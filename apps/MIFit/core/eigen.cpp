#include <cstdio>
#include <cmath>
#include <cstdlib>
#ifndef HUGE
#define HUGE 1e10    /* large positive number */
#endif
typedef double* matrix;
typedef double* vector;
void
eigen(
  matrix a,    /* input matrix (n x n) */
  matrix b,    /* output eigenvectors in columns */
  vector c,    /* output eigenvalues */
  int n) {
  extern void tred2(matrix a, matrix z, vector d, vector e, double tol, int n);
  extern void tql2(matrix z, vector d, vector e, double eps, int n);
  vector e;
  double epsilon, TAU;
  epsilon =  2.6E-17;       /* Set this to precision of a double-precision
                             * floating point value */
  TAU = 1. / HUGE;       /* HUGE is defined in <math.h>, TAU is
                            supposed to be  smallest possible positive value */

  e = (double*)calloc(n, sizeof e[0]);
  if (e == NULL) {
    return;             /* no easy recovery */
  }
  tred2(a, b, c, e, TAU, n);
  tql2(b, c, e, epsilon, n);
  free(e);
}

#define as(i, j) (a[(i)*n + (j)])
#define zs(i, j) (z[(i)*n + (j)])
void tred2(matrix a, matrix z, vector d, vector e, double tol, int n) {
  /*      algorithm from NUMERISCHE MATHEMATIK 11, 181-195 (1968) */
  /*      for Householder tridiagonalization of a, transformation */
  /*      matrix is returned in z, diagonal in d, and subdiagonal */
  /*      in e, tol is machine dependent.                         */
  int i, j, k;
  int L;
  double g, h, f, hh;

  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      zs(i, j) = as(i, j);
    }
  }

  for (i = n-1; i >= 1; i--) {
    L = i - 2;
    f = zs(i, i-1);
    g = 0.0;
    for (k = 0; k <= L; k++) {
      g += (zs(i, k) * zs(i, k));
    }
    h = g + f * f;

    /* skip this transformation if g is too small */
    if (g <= tol) {
      e[i] = f;
      h = 0.0;
    } else {
      L++;
      if (f >= 0.0) {
        g = -sqrt(h);
      } else {          g = sqrt(h);}
      e[i] = g;
      h -= (f*g);
      zs(i, i-1) = f - g;
      f = 0.0;
      for (j = 0; j <= L; j++) {
        zs(j, i) = zs(i, j) / h;
        g = 0.0;
        for (k = 0; k <= j; k++) {
          g += (zs(j, k) * zs(i, k));
        }
        for (k = j+1; k <= L; k++) {
          g += (zs(k, j) * zs(i, k));
        }
        e[j] = g/h;
        f += (g * zs(j, i));
      }
      hh  = f / (h + h);
      for (j = 0; j <= L; j++) {
        f = zs(i, j);
        g = e[j] - hh * f;
        e[j] = g;
        for (k = 0; k <= j; k++) {
          zs(j, k) -= (f*e[k] + g*zs(i, k));
        }
      }
    }
    d[i] = h;
  }
  d[0] = 0.0;
  e[0] = 0.0;
  for (i = 0; i < n; i++) {
    L = i - 1;
    if (d[i] != 0.0) {
      for (j = 0; j <= L; j++) {
        g = 0.0;
        for (k = 0; k <= L; k++) {
          g += (zs(i, k) * zs(k, j));
        }
        for (k = 0; k <= L; k++) {
          zs(k, j) -= (g*zs(k, i));
        }
      }
    }
    d[i] = zs(i, i);
    zs(i, i) = 1.0;
    for (j = 0; j <= L; j++) {
      zs(i, j) = 0.0;
      zs(j, i) = 0.0;
    }
  }
}

#define zs(i, j) (z[(i)*n + (j)])

void tql2(matrix z, vector d, vector e, double eps, int n) {
  /*      Algorithm from NUMERISCHE MATHEMATIK 11, 293-306 (1968)         */
  /*      For QL eigenvalue procedure for tridiagonal matrices.           */
  /*      This procedure provides a virtually orthonormal set of          */
  /*      eigenvectors.  When z is the output of tred2, the eigen-        */
  /*      vectors are those of the original matrix a.  When z is          */
  /*      the identity matrix the eigenvectors are those of the           */
  /*      tridiagonal matrix whose diagonal elements are given in         */
  /*      d and subdiagonal elements in e.                                */
  /*      Eigenvalues are returned in d and eigenvectors in z, by         */
  /*      columns, ordered in increasing order of eigenvalues.            */
  int i, L, j;
  int k, m;
  double h, p, r, b, c, f, g, s;

  for (i = 1; i < n; i++) {
    e[i-1] = e[i];
  }
  e[n-1] = 0.0;
  b = 0.0;
  f = 0.0;
  for (L = 0; L < n; L++) {
    j = 0;
    h = eps * (fabs(d[L]) + fabs(e[L]));
    if (b < h) {
      b = h;
    }

    /* look for small sub-diagonal element */
    for (m = L; fabs(e[m]) > b  && m < n ; m++) {
      ;
    }
    if (m != L) {
      do {
        if (j >= 30) {
          return;                               /* no easy recovery */
        }
        j++;

        /* form shift */
        p = (d[L+1] - d[L]) / (2.0 * e[L]);
        r = sqrt(p*p + 1.0);
        if (p < 0.0) {
          h = p - r;
        } else {          h = p + r;}
        h = d[L] - e[L]/h;
        for (i = L; i < n; i++) {
          d[i] -= h;
        }
        f += h;

        /* QL transformation */
        p = d[m];
        c = 1.0;
        s = 0.0;
        for (i = m-1; i >= L; i--) {
          g = c*e[i];
          h = c*p;
          if (fabs(p) >= fabs(e[i])) {
            c = e[i] / p;
            r = sqrt(c*c + 1.0);
            e[i+1] = s * p * r;
            s = c / r;
            c = 1.0 / r;
          } else {
            c = p / e[i];
            r = sqrt(c*c + 1.0);
            e[i+1] = s * e[i] * r;
            s = 1.0 / r;
            c *= s;
          }
          p = c*d[i] - s*g;
          d[i+1] = h + s * ( c*g + s*d[i]);

          /* form vector */
          for (k = 0; k < n; k++) {
            h = zs(k, i+1);
            zs(k, i+1) = s * zs(k, i) + c*h;
            zs(k, i) = c*zs(k, i) - s*h;
          }
        }
        e[L] = s * p;
        d[L] = c * p;
      } while (fabs(e[L]) > b);
    }
    d[L] += f;
  }

  /* order eigenvalues and eigenvectors */
  for (i = 0; i < n; i++) {
    k = i;
    p = d[i];
    for (j = i+1; j < n; j++) {
      if (d[j] < p) {
        k = j;
        p = d[j];
      }
    }
    if (k != i) {
      d[k] = d[i];
      d[i] = p;
      for (j = 0; j < n; j++) {
        p = zs(j, i);
        zs(j, i) = zs(j, k);
        zs(j, k) = p;
      }
    }
  }
}

