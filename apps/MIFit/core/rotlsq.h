#ifndef mifit_legacy_rotlsq_h
#define mifit_legacy_rotlsq_h

int rotlsqfit(double (*a)[3], double (*b)[3], double* w,
              int, double r[3][3], double v[3]);

#endif
