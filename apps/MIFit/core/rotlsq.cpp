#include <cctype>
#include <cstdlib>
#include <cstring>
#include <cassert>

#include "mathlib.h"

/* Finds rotation matrix r and vector t such that       */
/* SUM ( sqr (w* (b - r*a + t))) is minimized, where    */
/* a and b are vectors of m points, weighted by w.      */
/* These transformations are then applied to a to produce c */
/* By Tom Williams. Return of "vector" and "matrix" added by Mike Pique */
//int rotlsqfit(double a[][3], double b[][3], double w[],
//	int m, double matrix[3][3], double vector[3]);
//FILE *fp;
int
rotlsqfit(
  double (* a)[3],      /* input:  m fixed  3d points (double) */
  double (* b)[3],      /* input:  m guide  3d points (double) */
  double* w,                   /* input:  m weights */
  int m,               /* input:  number of points */
  double matrix[3][3], /* output: 3x3 rotation matrix 'r' */
  double vector[3]) {  /* output: 1x3 translation vector 't' (post-rotation) */

  
  double r[3][3];                   /* 3x3 matrix to rotate a onto b */
  double ta[3];                     /* translates origin to a */
  double tb[3];                     /* translates origin to b */
  int i, j;
  double wsum;
  extern int
  rotate(
    double a[][3],             /* input:  m fixed 3d points (double) */
    double b[][3],             /* input:  m guide 3d points (double) */
    double w[],                /* input:  m weights */
    int m,                     /* input:  number of points */
    double r[3][3]);           /* output: 3x3 matrix to rotate a onto b */
  extern int matvec(double*, double*, double*);

  if (m <= 2) { // return identity matrix and vector between
    for (i = 0; i < 3; i++) {
      for (j = 0; j < 3; j++) {
        if (i == j) {
          matrix[i][j] = 1;
        } else {matrix[i][j] = 0;}
      }
    }
    if (m <= 0) {
      return false;
    }
    vector[0] = b[0][0] - a[0][0];
    vector[1] = b[0][1] - a[0][1];
    vector[2] = b[0][2] - a[0][2];
    return true;
  }

#ifdef CDEBUG
  fp = fopen("debug.txt", "w");
  fprintf(fp, "rotlsqfit(%d points:\n", m);
  for (i = 0; i < m; i++) { /* dump arrays   */
    fprintf(fp,
      "pnt %d  weight %.2f fixed(%.1f %.1f %.1f) guide(%.1f %.1f %.1f)\n",
      i, w[i],
      a[i][0], a[i][1],  a[i][2],
      b[i][0], b[i][1],  b[i][2]);
  }
  fflush(fp);
#endif
  /* find center of mass of a and b */
  for (j = 0; j < 3; j++) {
    ta[j] = tb[j] = 0.0;
  }
  wsum = 0.0;

  for (i = 0; i < m; i++) {
    for (j = 0; j < 3; j++) {
      ta[j] += w[i] * a[i][j];
      tb[j] += w[i] * b[i][j];
    }
    wsum += w[i];
  }
  if (wsum == 0) {
    wsum = 1;
  }
  for (j = 0; j < 3; j++) {
    ta[j] = ta[j] / wsum;
    tb[j] = tb[j] / wsum;
  }

  for (i = 0; i < m; i++) {
    for (j = 0; j < 3; j++) {
      a[i][j] = a[i][j] - ta[j];
      b[i][j] = b[i][j] - tb[j];
    }
  }

  if (!rotate(a, b, w, m, r)) {
#ifdef CDEBUG
    dmatdump(r, m, m, "rotlsqfit - rotate failed,  matrix =");
#endif
    return (false);
  }
#ifdef CDEBUG
  fprintf(fp, "rotlsqfit - rotate returned matrix =\n");
  fprintf(fp, " %f %f %f\n %f %f %f\n %f %f %f\n", r[0][0], r[0][1], r[0][2],
    r[1][0], r[1][1], r[1][2],
    r[2][0], r[2][1], r[2][2]);
#endif
  if (matrix != NULL) {
    orthomatrix(r, matrix);
  }
  if (vector != NULL) {
    for (j = 0; j < 3; j++) {
      vector[j] = tb[j] - (r[j][0]*ta[0] + r[j][1]*ta[1] + r[j][2]*ta[2]);
    }
  }

  return true;
}

