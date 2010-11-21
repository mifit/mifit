#include <cmath>
#include "Matrices.h"
#include "crystmat.h"

void incmatrix(float xinc, float yinc, float zinc, float oldmat[3][3], float newmat[3][3])
{
    static float incmat[3][3]; /* small-angle incremental matrix */

    buildmat(xinc, yinc, zinc, incmat);
    matmul(incmat, oldmat, newmat);
}

void matmul(float l[3][3], float r[3][3], float prod[3][3])
{
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

void orthomatrix(double in[3][3], double out[3][3])
{
    // re-orthogonalize and renormalize a 3x3 rotation matrix using Gram-Schmidt

    /* set to cleaned-up copy of in. OK if same array as in */
    static double dotp;
    normvect(in[0], out[0]);
    normvect(in[1], out[1]);
    /* make Y-row orthogonal to X-row */
    dotp =  dotvect(out[0], out[1]);
    for (int i = 0; i < 3; ++i) out[1][i] -= dotp * out[0][i];
    normvect(out[1], out[1]);    /* renormalize Y-row */
    /*  form cross product Z-row of the two vectors X-row and Y-row */
    cross(out[0], out[1], out[2]);
}

void orthomatrix(float in[3][3], float out[3][3])
{
    // re-orthogonalize and renormalize a 3x3 rotation matrix using Gram-Schmidt

    /* set to cleaned-up copy of in. OK if same array as in */
    static float dotp;
    normvect(in[0], out[0]);
    normvect(in[1], out[1]);
    /* make Y-row orthogonal to X-row */
    dotp =  dotvect(out[0], out[1]);
    for (int i = 0; i < 3; ++i) out[1][i] -= dotp * out[0][i];
    normvect(out[1], out[1]);    /* renormalize Y-row */
    /*  form cross product Z-row of the two vectors X-row and Y-row */
    cross(out[0], out[1], out[2]);
}


void cross(float a[3], float b[3], float c[3])
{
    /* c := cross product of a and b */
    c[0] = a[1] * b[2] - a[2] * b[1];
    c[1] = a[2] * b[0] - a[0] * b[2];
    c[2] = a[0] * b[1] - a[1] * b[0];
}

void cross(double a[3], double b[3], double c[3])
{
    /* c := cross product of a and b */
    c[0] = a[1] * b[2] - a[2] * b[1];
    c[1] = a[2] * b[0] - a[0] * b[2];
    c[2] = a[0] * b[1] - a[1] * b[0];
}

void normvect(float in[3], float out[3])
{
    /* set to renormalized copy of in. OK if same vector as in */
    double normfactor = 1.0/sqrt(in[0]*in[0] + in[1]*in[1] + in[2]*in[2]);
    for (int i = 0; i < 3; ++i) out[i] = (float)((double)in[i]*normfactor);
}

void normvect(double in[3], double out[3])
{
    /* set to renormalized copy of in. OK if same vector as in */
    double normfactor = 1.0/sqrt(in[0]*in[0] + in[1]*in[1] + in[2]*in[2]);
    for (int i = 0; i < 3; ++i) out[i] = (double)in[i]*normfactor;
}

double dotvect(double in1[3], double in2[3])
{
    return in1[0]*in2[0] + in1[1]*in2[1] + in1[2]*in2[2];
}

float dotvect(float in1[3], float in2[3])
{
    return in1[0]*in2[0] + in1[1]*in2[1] + in1[2]*in2[2];
}

float vectorangle(float v1[3], float v2[3])
{
    /* angle between vectors v1 v2 in radians */
    float d1 = (float)sqrt(v1[0]*v1[0]+v1[1]*v1[1]+v1[2]*v1[2]);
    float d2 = (float)sqrt(v2[0]*v2[0]+v2[1]*v2[1]+v2[2]*v2[2]);
    float a = (float)acos((v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]) /(d1*d2));
    return a;
}

