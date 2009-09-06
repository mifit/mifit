#include <cstdio> // for the definition of NULL
#include <cmath>
#include "Matrices.h"

#define dtor(x) ( (x) * (3.141592654/180.))
#define X 0
#define Y 1
#define Z 2
#define do3(i) for (i = 0; i < 3; i++)
#define do4(i) for (i = 0; i < 4; i++)

void buildmat(float xdegrees, float ydegrees, float zdegrees, long int mat[3][3]) {
  /* small-angle approximation matrix */
  double sx, cx, sy, cy, sz, cz;
  double msf = MSF;
  if (xdegrees == 0.) {
    sx = 0.; cx = 1.;
  } else {
    sx = sin(dtor(xdegrees));
    cx = cos(dtor(xdegrees));
  }
  if (ydegrees == 0.) {
    sy = 0.; cy = 1.;
  } else {
    sy = sin(dtor(ydegrees));
    cy = cos(dtor(ydegrees));
  }
  if (zdegrees == 0.) {
    sz = 0.; cz = 1.;
  } else {
    sz = sin(dtor(zdegrees));
    cz = cos(dtor(zdegrees));
  }

  mat[X][X] = (long)(cy*cz*msf);
  mat[X][Y] = (long)(sz*msf);
  mat[X][Z] = (long)(sy*msf);

  mat[Y][X] = (long)(-sz*msf);
  mat[Y][Y] = (long)(cx*cz*msf);
  mat[Y][Z] = (long)(-sx*msf);

  mat[Z][X] = (long)(-sy*msf);
  mat[Z][Y] = (long)(sx*cy*msf);
  mat[Z][Z] = (long)(cx*cy*msf);
}

/*
   void buildmat(float xdegrees,float ydegrees,float zdegrees, float mat[3][3])
   {
   // small-angle approximation matrix
    double sx, cx,sy,cy, sz,cz;
    if(xdegrees==0.) {sx=0.; cx=1.;}
    else {
        sx = sin(dtor(xdegrees));
        cx = cos(dtor(xdegrees));
        }
    if(ydegrees==0.) {sy=0.; cy=1.;}
    else {
        sy = sin(dtor(ydegrees));
        cy = cos(dtor(ydegrees));
                }
        if(zdegrees==0.) {sz=0.; cz=1.;}
        else {
                sz = sin(dtor(zdegrees));
                cz = cos(dtor(zdegrees));
                }

    mat[X][X] = (float)(cy*cz);
    mat[X][Y] = (float)(sz);
    mat[X][Z] = (float)(sy);

    mat[Y][X] = (float)(-sz);
    mat[Y][Y] = (float)(cx*cz);
    mat[Y][Z] = (float)(-sx);

    mat[Z][X] = (float)(-sy);
    mat[Z][Y] = (float)(sx*cy);
    mat[Z][Z] = (float)(cx*cy);
   }
 */

void incmatrix(float xinc, float yinc, float zinc, long int oldmat[3][3], long int newmat[3][3]) {
  static long int incmat[3][3];   /* small-angle incremental matrix */

  buildmat(xinc, yinc, zinc, incmat);
  matmul(incmat, oldmat, newmat);
}

void incmatrix(float xinc, float yinc, float zinc, float oldmat[3][3], float newmat[3][3]) {
  static float incmat[3][3];   /* small-angle incremental matrix */


  buildmat(xinc, yinc, zinc, incmat);
  matmul(incmat, oldmat, newmat);
}

void
matmul(long int l[3][3], long int r[3][3], long int prod[3][3]) {
  /* matrix multiply l (left) by r (right) putting result in product.
     It's OK for any of the three to overlap
   */
  static long int tmp[3][3];

  tmp[0][0] = l[0][0]*r[0][0]+l[0][1]*r[1][0]+l[0][2]*r[2][0];
  tmp[0][1] = l[0][0]*r[0][1]+l[0][1]*r[1][1]+l[0][2]*r[2][1];
  tmp[0][2] = l[0][0]*r[0][2]+l[0][1]*r[1][2]+l[0][2]*r[2][2];

  tmp[1][0] = l[1][0]*r[0][0]+l[1][1]*r[1][0]+l[1][2]*r[2][0];
  tmp[1][1] = l[1][0]*r[0][1]+l[1][1]*r[1][1]+l[1][2]*r[2][1];
  tmp[1][2] = l[1][0]*r[0][2]+l[1][1]*r[1][2]+l[1][2]*r[2][2];

  tmp[2][0] = l[2][0]*r[0][0]+l[2][1]*r[1][0]+l[2][2]*r[2][0];
  tmp[2][1] = l[2][0]*r[0][1]+l[2][1]*r[1][1]+l[2][2]*r[2][1];
  tmp[2][2] = l[2][0]*r[0][2]+l[2][1]*r[1][2]+l[2][2]*r[2][2];

  prod[0][0] = tmp[0][0]/MSF;
  prod[0][1] = tmp[0][1]/MSF;
  prod[0][2] = tmp[0][2]/MSF;

  prod[1][0] = tmp[1][0]/MSF;
  prod[1][1] = tmp[1][1]/MSF;
  prod[1][2] = tmp[1][2]/MSF;

  prod[2][0] = tmp[2][0]/MSF;
  prod[2][1] = tmp[2][1]/MSF;
  prod[2][2] = tmp[2][2]/MSF;
}

void
matmul(float l[3][3], float r[3][3], float prod[3][3]) {
  /* matrix multiply l (left) by r (right) putting result in product.
     It's OK for any of the three to overlap
   */
  static float tmp[3][3];

  tmp[0][0] = l[0][0]*r[0][0]+l[0][1]*r[1][0]+l[0][2]*r[2][0];
  tmp[0][1] = l[0][0]*r[0][1]+l[0][1]*r[1][1]+l[0][2]*r[2][1];
  tmp[0][2] = l[0][0]*r[0][2]+l[0][1]*r[1][2]+l[0][2]*r[2][2];

  tmp[1][0] = l[1][0]*r[0][0]+l[1][1]*r[1][0]+l[1][2]*r[2][0];
  tmp[1][1] = l[1][0]*r[0][1]+l[1][1]*r[1][1]+l[1][2]*r[2][1];
  tmp[1][2] = l[1][0]*r[0][2]+l[1][1]*r[1][2]+l[1][2]*r[2][2];

  tmp[2][0] = l[2][0]*r[0][0]+l[2][1]*r[1][0]+l[2][2]*r[2][0];
  tmp[2][1] = l[2][0]*r[0][1]+l[2][1]*r[1][1]+l[2][2]*r[2][1];
  tmp[2][2] = l[2][0]*r[0][2]+l[2][1]*r[1][2]+l[2][2]*r[2][2];

  prod[0][0] = tmp[0][0];
  prod[0][1] = tmp[0][1];
  prod[0][2] = tmp[0][2];

  prod[1][0] = tmp[1][0];
  prod[1][1] = tmp[1][1];
  prod[1][2] = tmp[1][2];

  prod[2][0] = tmp[2][0];
  prod[2][1] = tmp[2][1];
  prod[2][2] = tmp[2][2];
}

void orthomatrix(long int in[3][3], long int out[3][3]) {
  /* set to cleaned-up copy of in. OK if same array as in */
  long int i;
  static long int dotp;
  /* re-orthogonalize and renormalize a 3x3 rotation matrix using
   * Gra[h]m-Schmidt
   * M Pique
   */
  normvect(in[X], out[X]);
  normvect(in[Y], out[Y]);
  /* make Y-row orthogonal to X-row */
  dotp =  dotvect(out[X], out[Y]);
  do3(i) out[Y][i] -= dotp * out[X][i]/MSF;
  normvect(out[Y], out[Y]);      /* renormalize Y-row */
  /*  form cross product Z-row of the two vectors X-row and Y-row */
  cross(out[X], out[Y], out[Z]);
}

void orthomatrix(double in[3][3], double out[3][3]) {
  /* set to cleaned-up copy of in. OK if same array as in */
  long int i;
  static double dotp;
  /* re-orthogonalize and renormalize a 3x3 rotation matrix using
   * Gra[h]m-Schmidt
   * M Pique
   */
  normvect(in[X], out[X]);
  normvect(in[Y], out[Y]);
  /* make Y-row orthogonal to X-row */
  dotp =  dotvect(out[X], out[Y]);
  do3(i) out[Y][i] -= dotp * out[X][i];
  normvect(out[Y], out[Y]);      /* renormalize Y-row */
  /*  form cross product Z-row of the two vectors X-row and Y-row */
  cross(out[X], out[Y], out[Z]);
}

void orthomatrix(float in[3][3], float out[3][3]) {
  /* set to cleaned-up copy of in. OK if same array as in */
  int i;
  static float dotp;
  /* re-orthogonalize and renormalize a 3x3 rotation matrix using
   * Gra[h]m-Schmidt
   * M Pique
   */
  normvect(in[X], out[X]);
  normvect(in[Y], out[Y]);
  /* make Y-row orthogonal to X-row */
  dotp =  dotvect(out[X], out[Y]);
  do3(i) out[Y][i] -= dotp * out[X][i];
  normvect(out[Y], out[Y]);      /* renormalize Y-row */
  /*  form cross product Z-row of the two vectors X-row and Y-row */
  cross(out[X], out[Y], out[Z]);
}

void cross(long a[3], long b[3], long c[3]) {
  /* c := cross product of a and b */
  c[X] = (a[Y] * b[Z] - a[Z] * b[Y])/MSF;
  c[Y] = (a[Z] * b[X] - a[X] * b[Z])/MSF;

  c[Z] = (a[X] * b[Y] - a[Y] * b[X])/MSF;
}

void cross(float a[3], float b[3], float c[3]) {
  /* c := cross product of a and b */
  c[X] = a[Y] * b[Z] - a[Z] * b[Y];
  c[Y] = a[Z] * b[X] - a[X] * b[Z];
  c[Z] = a[X] * b[Y] - a[Y] * b[X];
}

void cross(double a[3], double b[3], double c[3]) {
  /* c := cross product of a and b */
  c[X] = a[Y] * b[Z] - a[Z] * b[Y];
  c[Y] = a[Z] * b[X] - a[X] * b[Z];
  c[Z] = a[X] * b[Y] - a[Y] * b[X];
}

void normvect(long in[3], long out[3]) {
  /* set to renormalized copy of in. OK if same vector as in */
  register long i;
  double normfactor;
  normfactor = 1./sqrt(
                 (double)(in[X]*in[X] + in[Y]*in[Y] + in[Z]*in[Z])/1000000.0);
  do3(i) out[i] = (long)((double)in[i]*normfactor);
}

void normvect(float in[3], float out[3]) {
  /* set to renormalized copy of in. OK if same vector as in */
  int i;
  double normfactor;
  normfactor = 1./sqrt(in[X]*in[X] + in[Y]*in[Y] + in[Z]*in[Z]);
  do3(i) out[i] = (float)((double)in[i]*normfactor);
}

void normvect(double in[3], double out[3]) {
  /* set to renormalized copy of in. OK if same vector as in */
  register int i;
  double normfactor;
  normfactor = 1./sqrt(in[X]*in[X] + in[Y]*in[Y] + in[Z]*in[Z]);
  do3(i) out[i] = (double)in[i]*normfactor;
}

long
dotvect(long in1[3], long in2[3]) {
  return ( (in1[X]*in2[X] + in1[Y]*in2[Y] + in1[Z]*in2[Z])/MSF );
}

double
dotvect(double in1[3], double in2[3]) {
  return ( in1[X]*in2[X] + in1[Y]*in2[Y] + in1[Z]*in2[Z] );
}

float
dotvect(float in1[3], float in2[3]) {
  return ( in1[X]*in2[X] + in1[Y]*in2[Y] + in1[Z]*in2[Z] );
}

void uinv(
  long mat[3][3],
  long imat[3][3]) {
  long det;
  long i, j;
  imat[0][0] = (mat[1][1]*mat[2][2]-mat[2][1]*mat[1][2])/MSF;
  imat[1][0] = (mat[1][2]*mat[2][0]-mat[1][0]*mat[2][2])/MSF;
  imat[2][0] = (mat[1][0]*mat[2][1]-mat[1][1]*mat[2][0])/MSF;
  imat[0][1] = (mat[0][2]*mat[2][1]-mat[0][1]*mat[2][2])/MSF;
  imat[1][1] = (mat[0][0]*mat[2][2]-mat[0][2]*mat[2][0])/MSF;
  imat[2][1] = (mat[0][1]*mat[2][0]-mat[0][0]*mat[2][1])/MSF;
  imat[0][2] = (mat[0][1]*mat[1][2]-mat[0][2]*mat[1][1])/MSF;
  imat[1][2] = (mat[0][2]*mat[1][0]-mat[0][0]*mat[1][2])/MSF;
  imat[2][2] = (mat[0][0]*mat[1][1]-mat[0][1]*mat[1][0])/MSF;
  det = (mat[0][0]*mat[1][1]/MSF*mat[2][2]/MSF+mat[1][0]*mat[2][1]/MSF*mat[0][2]/MSF+
         mat[2][0]*mat[0][1]/MSF*mat[1][2]/MSF-mat[2][2]*mat[0][1]/MSF*mat[1][0]/MSF-
         mat[0][0]*mat[2][1]/MSF*mat[1][2]/MSF-mat[2][0]*mat[1][1]/MSF*mat[0][2]/MSF);
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      imat[j][i] = imat[j][i]*MSF/det;
    }
  }
}

float vectorangle(float v1[3], float v2[3]) {
  /* angle between vectors v1 v2 in radians */
  float d1, d2, a;

  d1 = (float)sqrt(v1[X]*v1[X]+v1[Y]*v1[Y]+v1[Z]*v1[Z]);
  d2 = (float)sqrt(v2[X]*v2[X]+v2[Y]*v2[Y]+v2[Z]*v2[Z]);
  a = (float)acos((v1[X]*v2[X] + v1[Y]*v2[Y] + v1[Z]*v2[Z]) /(d1*d2));
  return (a);
}

