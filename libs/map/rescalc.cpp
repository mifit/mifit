#include "rescalc.h"
#include <cstdio>
#include <cmath>

/* function to calculate the resolution of a reflection
   in Angstroms. The unit cell must be specified on the
   first call and init must be 1 to properly initialize -
   afterwards init is not 1.
 */
float rescalc(int ih, int ik, int il, float A, float B, float C, float ALPHA, float BETA, float GAMMA, int init) {
  return 0.5f/sthol(ih, ik, il, A, B, C, ALPHA, BETA, GAMMA, init);
}

float sthol(int ih, int ik, int il, float A, float B, float C, float ALPHA, float BETA, float GAMMA, int init) {

  double sthol;
  double ah, ak, al, V;
  static double as, bs, cs, cosgs, cosbs, cosas;
  static double cosa, cosb, cosg, sina, sinb, sing;
  static double alpha, beta, gamma;
  /*double t;*/

  if (init == 1) {
    if (A == 0.0 || B == 0.0 || C == 0.0 || ALPHA == 0.0 || BETA == 0.0 || GAMMA == 0.0) {
      return (-1.);
    }

    alpha = ALPHA / 180.0 * 3.1416 ;
    beta = BETA / 180.0 * 3.1416 ;
    gamma = GAMMA / 180.0 * 3.1416 ;
    cosa = cos(alpha);
    cosb = cos(beta);
    cosg = cos(gamma);
    sina = sin(alpha);
    sinb = sin(beta);
    sing = sin(gamma);
    cosas = (cosb*cosg - cosa)/sinb/sing ;
    cosbs = (cosa*cosg - cosb)/sina/sing ;
    cosgs = (cosa*cosb - cosg)/sina/sinb ;
    V = A*B*C * sqrt(1-cosa*cosa-cosb*cosb-cosg*cosg+2.0*cosa*cosg*cosb);
    as = B*C*sina/V;
    bs = A*C*sinb/V;
    cs = A*B*sing/V;
  } /* end of initialize */

  if ((ih == 0) && (ik == 0) && (il == 0)) {
    fprintf(stderr, "Error in rescalc():  ih=ik=il=0.\n");
    return -1.;
  }
  ah = ih ;
  ak = ik ;
  al = il ;
  sthol = 0.5 * sqrt(ah*ah*as*as + ak*ak*bs*bs + al*al*cs*cs
            +2.0*ah*ak*as*bs*cosgs
            +2.0*ah*al*as*cs*cosbs
            +2.0*ak*al*bs*cs*cosas) ;
  return (float)sthol;
}

float Volume(float A, float B, float C, float ALPHA, float BETA, float GAMMA) {

  double V;
  double cosgs, cosbs, cosas;
  double cosa, cosb, cosg, sina, sinb, sing;
  double alpha, beta, gamma;
  /*double t;*/

  alpha = ALPHA / 180.0 * 3.1416 ;
  beta = BETA / 180.0 * 3.1416 ;
  gamma = GAMMA / 180.0 * 3.1416 ;
  cosa = cos(alpha);
  cosb = cos(beta);
  cosg = cos(gamma);
  sina = sin(alpha);
  sinb = sin(beta);
  sing = sin(gamma);
  cosas = (cosb*cosg - cosa)/sinb/sing ;
  cosbs = (cosa*cosg - cosb)/sina/sing ;
  cosgs = (cosa*cosb - cosg)/sina/sinb ;
  V = A*B*C * sqrt(1-cosa*cosa-cosb*cosb-cosg*cosg+2.0*cosa*cosg*cosb);

  return ((float)V);
}

