#include "LSQMatrix.h"
#include <stdio.h>

bool LSQMatrix::Save(const char* pathname) {
  FILE* fp;
  if (IsOK()) {
    fp = fopen(pathname, "w");
    if (fp) {
      fprintf(fp, "matrix %0.4f %0.4f %0.4f %0.4f %0.4f %0.4f %0.4f %0.4f %0.4f\n",
        mat[0][0], mat[0][1], mat[0][2],
        mat[1][0], mat[1][1], mat[1][2],
        mat[2][0], mat[2][1], mat[2][2]);
      fprintf(fp, "vector %0.3f %0.3f %0.3f\n", vec[0], vec[1], vec[2]);
      fclose(fp);
      return true;
    }
  }
  return false;
}

bool LSQMatrix::Load(const char* pathname) {
  FILE* fp;
  fp = fopen(pathname, "r");
  if (fp) {
    fscanf(fp, "matrix %f%f%f%f%f%f%f%f%f\n",
      &mat[0][0], &mat[0][1], &mat[0][2],
      &mat[1][0], &mat[1][1], &mat[1][2],
      &mat[2][0], &mat[2][1], &mat[2][2]);
    fscanf(fp, "vector %f%f%f\n", &vec[0], &vec[1], &vec[2]);
    fclose(fp);
    init = true;
    return true;
  }
  return false;
}

float LSQMatrix::Xvalue(float tx, float ty, float tz) {
  float xout = tx;
  if (IsOK()) {
    xout = mat[0][0]*tx+mat[0][1]*ty+mat[0][2]*tz
           + vec[0];
  }
  return xout;
}

float LSQMatrix::Yvalue(float tx, float ty, float tz) {
  float yout = ty;
  if (IsOK()) {
    yout = mat[1][0]*tx+mat[1][1]*ty+mat[1][2]*tz
           + vec[1];
  }
  return yout;
}

float LSQMatrix::Zvalue(float tx, float ty, float tz) {
  float zout = tz;
  if (IsOK()) {
    zout = mat[2][0]*tx+mat[2][1]*ty+mat[2][2]*tz
           + vec[2];
  }
  return zout;
}

void LSQMatrix::SetMatrix(float m[3][3], float v[3]) {
  for (int i = 0; i < 3 ; i++) {
    for (int j = 0; j < 3; j++) {
      mat[i][j] = m[i][j];
    }
  }
  vec[0] = v[0];
  vec[1] = v[1];
  vec[2] = v[2];
  init = true;
}

void LSQMatrix::SetMatrix(double m[3][3], double v[3]) {
  for (int i = 0; i < 3 ; i++) {
    for (int j = 0; j < 3; j++) {
      mat[i][j] = (float)m[i][j];
    }
  }
  vec[0] = (float)v[0];
  vec[1] = (float)v[1];
  vec[2] = (float)v[2];
  init = true;
}

void LSQMatrix::GetMatrix(double m[3][3], double v[3]) {
  for (int i = 0; i < 3 ; i++) {
    for (int j = 0; j < 3; j++) {
      m[i][j] = mat[i][j];
    }
  }
  v[0] = vec[0];
  v[1] = vec[1];
  v[2] = vec[2];
}

