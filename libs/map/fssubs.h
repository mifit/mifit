#ifndef MI_FFSUBS_H
#define MI_FFSUBS_H

#include <stdio.h>
#include <vector>
class CMapHeaderBase;

int fsread(int* rho, int* nx, int* ny, int* nz, float* rhoscale, float* rmsrho, FILE* fp, char* err);
int fssize(FILE* fp);
int fsread_uni(std::vector<int>& rho, int* nx, int* ny, int* nz, float* rhoscale, float* rmsrho, FILE* fp, char* err, int swab);
int fssize_uni(FILE* fp, int* swab);
int fswrite(float* rho, FILE* fp, CMapHeaderBase* mh, int norn, int ncent);
#endif
