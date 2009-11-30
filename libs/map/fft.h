#ifndef FFT_H
#define FFT_H

// local, private header for maplib

#include "sfcalc.h"
#include "maptypes.h"

typedef struct fcomplex
{
    float re;
    float im;
} fcomplex;

int cmplft_(float*, float *y, long int *n, long int *d);
int hermft(float *x, float *y, long int *n, long int *dim);


//@{
// Puts a reflection's indices into a standard quadrant and at calculates some useful
// values.
// <pre>
//@@ #include <mifit/legacy/Xguicryst.h>
// Example of how to transfer symmetry information from a CMapHeaderBase structure, mhin.
//	integer ih,ik,il;
//	integer iss[864];
//	integer its[288];
//	float eps, wt;
//	integer mult,mk,iflg,iflg2,nsym;
//	for(k=0;k<mhin->nsym;k++){
//		for(j=0;j<3;j++){
//			its[j+k*3]=  ROUND(mhin->symops[j][3][k]*24.0);
//			for(i=0;i<3;i++)
//				iss[i+3*(j+3*k)]= (int)mhin->symops[i][j][k];
//		}
//	}
// </pre>
//@}

int stdrefl_(long int *ih, long int *ik, long int *il, long int *mult, float *eps,
             long int *mk, long int *iflg, long int *iflg2, long int *iss, long int *its,
             long int *nsym);
int MIMapFactor(int ntest, int prime, int even, int inc);
int SigmaA();
float *fft3d(CREFL *refl, int nrefl, CMapHeaderBase *mapheader, int usepsi);
void get_unit_cell();
void cycle(int *x, int *y, int *z);
double psi_(int p, int N, int k);
void read_sf(float *x, int nx, int ny, int nz, float scale, float temp, CREFL *refl, int nrefl, int usepsi);
void expansion_error(int nsymop);
int xpnd(int h1[], int h2[], fcomplex * f1, fcomplex * f2, int n);
void wplane(float *x, int nx, int ny, float *p, int mx, int my, int x0, int y_0, int level);
int nextf(int h[3], fcomplex * f, int i, CREFL * refl, int nrefl);
void cexp(fcomplex *a, float b);
int str_index(char *str1, char *str2);
int my_index(char *str, char ch);
void put_80(char *cptr, FILE *file);
double B_(double x, int j, int N);
int SigmaA_C(CREFL *refl, int nrefl, CMapHeaderBase *mhin);
float sim(float x);
int smi(float *am, long int *nm, long int *n, long int *nfail);
int solv(float *am, float *v, float *dv, float *diag, long int *nm, long int *nv, float *sig);
int matset(float *am, float *dv, long int *nm, long int *nv);

#endif // ifndef FFT_H
