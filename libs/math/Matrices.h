#ifndef MATRICES_H
#define MATRICES_H

// matrix scale factor
//@{
// Matrix scale factor.
//@}
#define MSF 1000L

//@{
// calculate a matrix to rotate incrementally given the old rotation matrix.
// floating point version
//@}
void incmatrix(float xinc, float yinc, float zinc, float oldmat[3][3], float newmat[3][3]);

//@{
// floating-point matrix multiplication.
//@}
void matmul(float l[3][3], float r[3][3], float prod[3][3]);

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
// angle between two vectors
//@}
float vectorangle(float v1[3], float v2[3]);


#endif // ifndef MATRICES_H

