#ifndef MATRICES_H
#define MATRICES_H

// matrix scale factor
//@{
// Matrix scale factor.
//@}
#define MSF 1000L
//@{
// half the matrix scale factor.
//@}
#define MSF2 500L
//@{
// round using the matrix scale factor.
//@}
#define MROUND(x) ( ((x) >= 0) ? ((int)((x+MSF2)/MSF)) : ((int)((x-MSF2)/MSF)))

//@{
// build a matrix for rotating around three axes.
// integer version for fixed-point arithmetic
//@}
void buildmat(float xdegrees, float ydegrees, float zdegrees, long int mat[3][3]);
//@{
// build a matrix for rotating around three axes.
// floating point version
//@}
void buildmat(float xdegrees, float ydegrees, float zdegrees, float mat[3][3]);
//@{
// calculate a matrix to rotate incrementally given the old rotation matrix.
// integer version for fixed-point arithmetic
//@}
void incmatrix(float xinc, float yinc, float zinc, long int oldmat[3][3], long int newmat[3][3]);
//@{
// calculate a matrix to rotate incrementally given the old rotation matrix.
// floating point version
//@}
void incmatrix(float xinc, float yinc, float zinc, float oldmat[3][3], float newmat[3][3]);
//@{
// integer matrix multiplication.
//@}
void matmul(long int l[3][3], long int r[3][3], long int prod[3][3]);
//@{
// floating-point matrix multiplication.
//@}
void matmul(float l[3][3], float r[3][3], float prod[3][3]);
//@{
// orthogonaliz the matrix (make the determinant 1).
//@}
void orthomatrix(long int in[3][3], long int out[3][3]);
//@{
// cross two vectors
// integer version for fixed-point arithmetic
//@}
void cross(long a[3], long b[3], long c[3]);
//@{
// normalize a vector
//@}
void normvect(long in[3], long out[3]);
//@{
// dot product of two vectors
//@}
long dotvect(long in1[3], long in2[3]);
//@{
// orthogonalize a matrix
//@}
void orthomatrix(double in[3][3], double out[3][3]);
//@{
// orthogonalize a matrix
//@}
void orthomatrix(float in[3][3], float out[3][3]);
//@{
// cross product of two vectors
// c is the cross of a and b.
//@}
void cross(double a[3], double b[3], double c[3]);
//@{
// cross product of two vectors
// c is the cross of a and b.
//@}
void cross(float a[3], float b[3], float c[3]);
//@{
// normalize a vector
//@}
void normvect(double in[3], double out[3]);
//@{
// normalize a vector
//@}
void normvect(float in[3], float out[3]);
//@{
// dot product of two vectors
//@}
double dotvect(double in1[3], double in2[3]);
//@{
// dot product of two vectors
//@}
float dotvect(float in1[3], float in2[3]);
//@{
// invert a matrix
//@}
void uinv(long int mat[3][3], long int imat[3][3]);
//@{
// invert a matrix
//@}
void uinv(float mat[3][3], float imat[3][3]);
//@{
// angle between two vectors
//@}
float vectorangle(float v1[3], float v2[3]);


#endif

