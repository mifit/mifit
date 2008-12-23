#ifndef CRYSTMAT_H
#define CRYSTMAT_H

void transform(float ut[3][3], float*, float*, float*);
int orthog(float a, float b, float c, float alpha, float beta, float gamma, float ut[3][3]);
void uinv(float mat[3][3], float imat[3][3]);
void uinvd(double mat[3][3], double imat[3][3]);
void initrotate(float, float, float, float, float, float, float, float mat[4][3]);
void rotate(float* x, float* y, float* z, float mat[4][3]);
void rotate(float v[3], float mat[4][3]);
void buildmat(float xdegrees, float ydegrees, float zdegrees, float mat[3][3]);
void xl_rotate(float x, float y, float z, float* xp, float* yp, float* zp, float mat[4][3]);

#endif
