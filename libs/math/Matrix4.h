#ifndef mi_math_Matrix4_h
#define mi_math_Matrix4_h

#include <math/Tuple4.h>
#include <math/Vector4.h>
#include <math/IndexOutOfBoundsException.h>
#include <string>

namespace mi {
namespace math {
template<class Type> class Matrix4;
}
}
#include <math/Quaternion.h>

namespace mi {
namespace math {

/**
 * A 4 x 4 matrix of Types.
 */
template<class Type>
class Matrix4 {
public:
  /**
   * The first element of the first row.
   */
  Type m00;

  /**
   * The second element of the first row.
   */
  Type m01;

  /**
   * The third element of the first row.
   */
  Type m02;

  /**
   * The fourth element of the first row.
   */
  Type m03;

  /**
   * The first element of the second row.
   */
  Type m10;

  /**
   * The second element of the second row.
   */
  Type m11;

  /**
   * The third element of the second row.
   */
  Type m12;

  /**
   * The fourth element of the second row.
   */
  Type m13;

  /**
   * The first element of the third row.
   */
  Type m20;

  /**
   * The second element of the third row.
   */
  Type m21;

  /**
   * The third element of the third row.
   */
  Type m22;

  /**
   * The fourth element of the third row.
   */
  Type m23;

  /**
   * The first element of the fourth row.
   */
  Type m30;

  /**
   * The second element of the fourth row.
   */
  Type m31;

  /**
   * The third element of the fourth row.
   */
  Type m32;

  /**
   * The fourth element of the fourth row.
   */
  Type m33;

  /**
   * Constructs and initializes a Matrix4 to all zeros.
   */
  Matrix4();

  /**
   * Constructs and initializes a Matrix4 from the specified 16 values.
   *
   * @param m00 the [0][0] element
   * @param m01 the [0][1] element
   * @param m02 the [0][2] element
   * @param m03 the [0][3] element
   * @param m10 the [1][0] element
   * @param m11 the [1][1] element
   * @param m12 the [1][2] element
   * @param m13 the [1][3] element
   * @param m20 the [2][0] element
   * @param m21 the [2][1] element
   * @param m22 the [2][2] element
   * @param m23 the [2][3] element
   * @param m30 the [3][0] element
   * @param m31 the [3][1] element
   * @param m32 the [3][2] element
   * @param m33 the [3][3] element
   */
  Matrix4(Type m00, Type m01, Type m02, Type m03,
          Type m10, Type m11, Type m12, Type m13,
          Type m20, Type m21, Type m22, Type m23,
          Type m30, Type m31, Type m32, Type m33);

  /**
   * Constructs and initializes a Matrix4 from the specified 16
   * element array.  this.m00=v[0], this.m01=v[1], etc.
   *
   * @param  v the array of length 16 containing in order
   */
  Matrix4(const Type v[]);

  /**
   * Creates a matrix with the value derived from the quaternion.
   *
   * @param q1  q quaternion.
   */
  Matrix4(const Quaternion<Type>& q1);

  /**
   * Sets 16 values
   *
   * @param m00 the [0][0] element
   * @param m01 the [0][1] element
   * @param m02 the [0][2] element
   * @param m03 the [0][3] element
   * @param m10 the [1][0] element
   * @param m11 the [1][1] element
   * @param m12 the [1][2] element
   * @param m13 the [1][3] element
   * @param m20 the [2][0] element
   * @param m21 the [2][1] element
   * @param m22 the [2][2] element
   * @param m23 the [2][3] element
   * @param m30 the [3][0] element
   * @param m31 the [3][1] element
   * @param m32 the [3][2] element
   * @param m33 the [3][3] element
   */
  void set(Type m00, Type m01, Type m02, Type m03,
           Type m10, Type m11, Type m12, Type m13,
           Type m20, Type m21, Type m22, Type m23,
           Type m30, Type m31, Type m32, Type m33);

  /**
   * Sets the values in this Matrix4 equal to the row-major
   * array parameter (ie, the first four elements of the array
   * will be copied into the first row of this matrix, etc.).
   */
  void set(const Type m[]);

  /**
   * Sets the value of this matrix to a copy of the passed
   * matrix m1.
   *
   * @param m1 the matrix to be copied.
   */
  void set(const Matrix4& m1);

  /**
   * Sets the value of this matrix to the value derived from the
   * quaternion.
   *
   * @param q1  q quaternion.
   */
  void set(const Quaternion<Type>& q1);

  /**
   * Sets this Matrix4 to identity.
   */
  void setIdentity();

  /**
   * Sets this matrix to all zeros.
   */
  void setZero();

  /**
   * Returns true if all of the data members of Matrix4 m1 are
   * equal to the corresponding data members in this Matrix4.
   *
   * @param m1 The matrix with which the comparison is made.
   * @return true or false
   */
  bool operator==(const Matrix4& m1) const;

  /**
   * Returns true if the difference of each of the elements of this
   * matrix and the corresponding element of matrix m1 are less than
   * or equal to the epsilon parameter.
   *
   * @param m1  the other matrix.
   * @param epsilonValue  the threshold value.
   */
  bool epsilonEquals(const Matrix4& m1, Type epsilonValue) const;

  /**
   * Adds a scalar to each component of this matrix.
   *
   * @param scalar  the value to add.
   */
  void add(Type scalar);

  /**
   * Adds a scalar to each component of this matrix.
   *
   * @param scalar  the value to add.
   */
  Matrix4& operator+=(Type scalar);


  /**
   * Adds a scalar to each component of the matrix m1 and places
   * the result into this. Matrix m1 is not modified.
   *
   * @param scalar  the value to add.
   * @parm m1  the original matrix values.
   */
  void add(Type scalar, const Matrix4& m1);

  /**
   * Sets the value of this matrix to sum of itself and matrix
   * m1.
   *
   * @param m1  the other matrix.
   */
  void add(const Matrix4& m1);

  /**
   * Sets the value of this matrix to the sum of itself and matrix m1.
   *
   * @param m1  the other matrix.
   */
  Matrix4& operator+=(const Matrix4& m1);

  /**
   * Sets the value of this matrix to the matrix sum of matrices
   * m1 and m2.
   *
   * @param m1  the first matrix.
   * @param m2  the second matrix.
   */
  void add(const Matrix4& m1, const Matrix4& m2);

  /**
   * Returns a new Matrix4 with the value of the sum of this matrix
   * and matrix m1.
   *
   * @param m1  the other matrix.
   * @return  a new Matrix4.
   */
  Matrix4 operator+(const Matrix4& m1) const;

  /**
   * Substracts a scalar from each component of this matrix.
   *
   * @param scalar  the value to subtract.
   */
  void subtract(Type scalar);

  /**
   * Substracts a scalar from each component of this matrix.
   *
   * @param scalar  the value to subtract.
   */
  Matrix4& operator-=(Type scalar);

  /**
   * Subtracts a scalar to each component of the matrix m1 and places
   * the result into this. Matrix m1 is not modified.
   *
   * @param scalar  the value to subtract.
   * @parm m1  the original matrix values.
   */
  void subtract(Type scalar, const Matrix4& m1);

  /**
   * Sets the value of this matrix to the matrix difference of
   * itself and matrix m1.
   *
   * @param m1  the other matrix
   */
  void subtract(const Matrix4& m1);

  /**
   * Sets the value of this matrix to the difference of itself and
   * matrix m1.
   *
   * @param m1 the other matrix.
   */
  Matrix4& operator-=(const Matrix4& m1);

  /**
   * Sets the value of this matrix to the matrix difference of
   * matrices m1 and m2.
   *
   * @param m1 the first matrix
   * @param m2 the second matrix
   */
  void subtract(const Matrix4& m1, const Matrix4& m2);

  /**
   * Returns a new Matrix4 with the value of the difference
   * of this matrix and m1.
   *
   * @param m1  the first matrix.
   * @return  a new Matrix4.
   */
  Matrix4 operator-(const Matrix4& m1) const;

  /**
   * Sets the value of this matrix to its transpose.
   */
  void transpose();

  /**
   * Sets the value of this matrix to the transpose of the
   * argument matrix
   *
   * @param m1 the matrix to be transposed
   */
  void transpose(const Matrix4& m1);

  void getInColumnMajorOrder(Type* matrix) const;

  /**
   * Returns the value at the first row and first column of this matrix.
   */
  Type get00() const;

  /**
   * Returns the value at the first row and second column of this matrix.
   */
  Type get01() const;

  /**
   * Returns the value at the first row and third column of this matrix.
   */
  Type get02() const;

  /**
   * Returns the value at the first row and fourth column of this matrix.
   */
  Type get03() const;

  /**
   * Returns the value at the second row and first column of this matrix.
   */
  Type get10() const;

  /**
   * Returns the value at the second row and second column of this matrix.
   */
  Type get11() const;

  /**
   * Returns the value at the second row and third column of this matrix.
   */
  Type get12() const;

  /**
   * Returns the value at the second row and fourth column of this matrix.
   */
  Type get13() const;

  /**
   * Returns the value at the third row and first column of this matrix.
   */
  Type get20() const;

  /**
   * Returns the value at the third row and second column of this matrix.
   */
  Type get21() const;

  /**
   * Returns the value at the third row and third column of this matrix.
   */
  Type get22() const;

  /**
   * Returns the value at the third row and fourth column of this matrix.
   */
  Type get23() const;

  /**
   * Returns the value at the fourth row and first column of this matrix.
   */
  Type get30() const;

  /**
   * Returns the value at the fourth row and second column of this matrix.
   */
  Type get31() const;

  /**
   * Returns the value at the fourth row and third column of this matrix.
   */
  Type get32() const;

  /**
   * Returns the value at the fourth row and fourth column of this matrix.
   */
  Type get33() const;

  /**
   * Sets the value of the first row and first column of this matrix.
   */
  void set00(const Type& value);

  /**
   * Sets the value of the first row and second column of this matrix.
   */
  void set01(const Type& value);

  /**
   * Sets the value of the first row and third column of this matrix.
   */
  void set02(const Type& value);

  /**
   * Sets the value of the first row and fourth column of this matrix.
   */
  void set03(const Type& value);

  /**
   * Sets the value of the second row and first column of this matrix.
   */
  void set10(const Type& value);

  /**
   * Sets the value of the second row and second column of this matrix.
   */
  void set11(const Type& value);

  /**
   * Sets the value of the second row and third column of this matrix.
   */
  void set12(const Type& value);

  /**
   * Sets the value of the second row and fourth column of this matrix.
   */
  void set13(const Type& value);

  /**
   * Sets the value of the third row and first column of this matrix.
   */
  void set20(const Type& value);

  /**
   * Sets the value of the third row and second column of this matrix.
   */
  void set21(const Type& value);

  /**
   * Sets the value of the third row and third column of this matrix.
   */
  void set22(const Type& value);

  /**
   * Sets the value of the third row and fourth column of this matrix.
   */
  void set23(const Type& value);

  /**
   * Sets the value of the fourth row and first column of this matrix.
   */
  void set30(const Type& value);

  /**
   * Sets the value of the fourth row and second column of this matrix.
   */
  void set31(const Type& value);

  /**
   * Sets the value of the fourth row and third column of this matrix.
   */
  void set32(const Type& value);

  /**
   * Sets the value of the fourth row and fourth column of this matrix.
   */
  void set33(const Type& value);

  /**
   * Sets the specified element of this matrix4 to the value provided.
   *
   * @param row  the row number to be modified (zero indexed)
   * @param column  the column number to be modified (zero indexed)
   * @param value the new value
   */
  void setElement(int row, int column, Type value);

  /**
   * Retrieves the value at the specified row and column of this matrix.
   *
   * @param row  the row number to be retrieved (zero indexed)
   * @param column  the column number to be retrieved (zero indexed)
   * @return the value at the indexed element
   */
  Type getElement(int row, int column) const;

  /**
   * Sets the specified row of this matrix4 to the Vector
   * provided.
   *
   * @param row the row number to be modified (zero indexed)
   * @param v the replacement row
   */
  void setRow(int row, const Vector4<Type>& v);

  /**
   * Copies the matrix values in the specified row into the
   * vector parameter.
   * @param row the matrix row
   * @param v The vector into which the matrix row values will be copied
   */
  void getRow(int row, Vector4<Type>& v) const;

  /**
   * Sets the specified column of this matrix4 to the vector
   * provided.
   *
   * @param column the column number to be modified (zero indexed).
   * @param v the replacement column.
   */
  void setColumn(int column, const Vector4<Type>& v);

  /**
   * Copies the matrix values in the specified column into the
   * vector parameter.
   *
   * @param column the matrix column
   * @param v The vector into which the matrix column values
   *   will be copied
   */
  void getColumn(int column, Vector4<Type>& v) const;

  /**
   * Sets the value of this matrix to the matrix inverse
   * of the passed matrix m1.
   *
   * @param m1 the matrix to be inverted
   */
  void invert(const Matrix4& m1);

  /**
   * Sets the value of this matrix to its inverse.
   */
  void invert();

  /**
   * Computes the determinant of this matrix.
   *
   * @return the determinant of the matrix
   */
  Type determinant() const;

  /**
   * Sets the value of this matrix to a translation matrix by
   * the given coordinates.
   */
  void translation(Type x, Type y, Type z);

  /**
   * Sets the value of this matrix to a translation matrix by
   * the given vector.
   */
  void translation(const Vector3<Type>& vector);

  /**
   * Sets the value of this matrix to a scaling matrix by
   * the given coordinates.
   */
  void scaling(Type x, Type y, Type z);

  /**
   * Sets the value of this matrix to a scaling matrix by
   * the given vector.
   */
  void scaling(const Vector3<Type>& vector);

  /**
   * Sets the value of this matrix to a rotation matrix about
   * the x axis by the passed angle.
   *
   * @param angle the angle to rotate about the X axis in radians
   */
  void rotationX(Type angle);

  /**
   * Sets the value of this matrix to a rotation matrix about
   * the y axis by the passed angle.
   *
   * @param angle the angle to rotate about the Y axis in radians */
  void rotationY(Type angle);

  /**
   * Sets the value of this matrix to a rotation matrix about
   * the z axis by the passed angle.
   *
   * @param angle the angle to rotate about the Z axis in radians
   */
  void rotationZ(Type angle);

  /**
   * Multiplies each element of this matrix by a scalar.
   *
   * @param scalar The scalar multiplier.
   */
  void multiply(Type scalar);

  /**
   * Multiplies each element of matrix m1 by a scalar and places
   * the result into this. Matrix m1 is not modified.
   *
   * @param scalar The scalar multiplier.
   * @param m1 The original matrix.
   */
  void multiply(Type scalar, const Matrix4& m1);

  /**
   * Sets the value of this matrix to the result of multiplying
   * itself with matrix m1.
   *
   * @param m1 the other matrix */
  void multiply(const Matrix4& m1);

  /**
   * Returns a matrix with the value of this matrix multiplying
   * itself with matrix m2.
   *
   * @param m2 the other matrix */
  Matrix4 operator*(const Matrix4& m2) const;

  /**
   * Sets the value of this matrix to the result of multiplying
   * the two argument matrices together.
   *
   * @param m1 the first matrix
   * @param m2 the second matrix
   */
  void multiply(const Matrix4& m1, const Matrix4& m2);

  /**
   * Transform the vector vec using this Matrix4 and place the
   * result into vecOut.
   *
   * @param vec the T precision vector to be transformed
   * @param vecOut the vector into which the transformed values
   *   are placed
   */
  void transform(const Tuple4<Type>& vec, Tuple4<Type>& vecOut) const;

  /**
   * Transform the vector vec using this Matrix4 and place the
   * result back into vec.
   *
   * @param vec the T precision vector to be transformed
   */
  void transform(Tuple4<Type>& vec) const;

  /**
   * Transform the vector vec using this Matrix4 and place the
   * result into vecOut.
   *
   * @param vec the T precision vector to be transformed
   * @param vecOut the vector into which the transformed values
   *   are placed
   */
  void transform(const Tuple3<Type>& vec, Tuple3<Type>& vecOut) const;

  /**
   * Transform the vector vec using this Matrix4 and place the
   * result back into vec.
   *
   * @param vec the T precision vector to be transformed
   */
  void transform(Tuple3<Type>& vec) const;

  /**
   * Negates the value of this matrix: this = -this.
   */
  void negate();

  /**
   * Sets the value of this matrix equal to the negation of of
   * the Matrix4 parameter.
   *
   * @param m1 The source matrix
   */
  void negate(const Matrix4& m1);

  /**
   * Returns a new Matrix4 with the value of the negation of
   * this matrix.
   *
   * @return a new Matrix4.
   */
  Matrix4 operator-() const;

  /**
   * Sets each element to its absolute value.
   */
  void absolute();

  /**
   * Returns a string that contains the values of this Matrix4.
   * @return the String representation
   */
  std::string toString() const;

};

/**
 * Write a matrix to an output stream. This method is designed to
 * be used on a character stream, and may not work for other
 * streams (e.g., binary stream).
 *
 * @param stream  the output stream written to.
 * @param m1  the matrix to write.
 * @return the modified output stream.
 */
template<class Type>
std::ostream& operator<<(std::ostream& stream, const Matrix4<Type>& m1);

/**
 * Adds a scalar to each component of the matrix m1 and places
 * the result into a new matrix. Matrix m1 is not modified.
 *
 * @param scalar  the value to add.
 * @parm m1  the original matrix values.
 */
template<class Type>
Matrix4<Type> operator+(Type scalar, const Matrix4<Type>& m1);

/**
 * Adds a scalar to each component of the matrix m1 and places
 * the result into a new matrix. Matrix m1 is not modified.
 *
 * @param scalar  the value to add.
 * @parm m1  the original matrix values.
 */
template<class Type>
Matrix4<Type> operator+(const Matrix4<Type>& m1, Type scalar);

/**
 * Subtracts a scalar to each component of the matrix m1 and places
 * the result into a new matrix. Matrix m1 is not modified.
 *
 * @param scalar  the value to subtract.
 * @parm m1  the original matrix values.
 */
template<class Type>
Matrix4<Type> operator-(Type scalar, const Matrix4<Type>& m1);

/**
 * Subtracts a scalar to each component of the matrix m1 and places
 * the result into a new matrix. Matrix m1 is not modified.
 *
 * @param scalar  the value to subtract.
 * @parm m1  the original matrix values.
 */
template<class Type>
Matrix4<Type> operator-(const Matrix4<Type>& m1, Type scalar);


template<class Type>
Matrix4<Type>::Matrix4()
  : m00(0.0), m01(0.0), m02(0.0), m03(0.0),
  m10(0.0), m11(0.0), m12(0.0), m13(0.0),
  m20(0.0), m21(0.0), m22(0.0), m23(0.0),
  m30(0.0), m31(0.0), m32(0.0), m33(0.0) {
}

template<class Type>
Matrix4<Type>::Matrix4(Type m00, Type m01, Type m02, Type m03,
                       Type m10, Type m11, Type m12, Type m13,
                       Type m20, Type m21, Type m22, Type m23,
                       Type m30, Type m31, Type m32, Type m33)
  : m00(m00), m01(m01), m02(m02), m03(m03),
  m10(m10), m11(m11), m12(m12), m13(m13),
  m20(m20), m21(m21), m22(m22), m23(m23),
  m30(m30), m31(m31), m32(m32), m33(m33) {
}

template<class Type>
Matrix4<Type>::Matrix4(const Type v[])
  : m00(v[0]), m01(v[1]), m02(v[2]), m03(v[3]),
  m10(v[4]), m11(v[5]), m12(v[6]), m13(v[7]),
  m20(v[8]), m21(v[9]), m22(v[10]), m23(v[11]),
  m30(v[12]), m31(v[13]), m32(v[14]), m33(v[15]) {
}

template<class Type>
void Matrix4<Type>::set(Type m00, Type m01, Type m02, Type m03,
                        Type m10, Type m11, Type m12, Type m13,
                        Type m20, Type m21, Type m22, Type m23,
                        Type m30, Type m31, Type m32, Type m33) {
  this->m00 = m00;
  this->m01 = m01;
  this->m02 = m02;
  this->m03 = m03;

  this->m10 = m10;
  this->m11 = m11;
  this->m12 = m12;
  this->m13 = m13;

  this->m20 = m20;
  this->m21 = m21;
  this->m22 = m22;
  this->m23 = m23;

  this->m30 = m30;
  this->m31 = m31;
  this->m32 = m32;
  this->m33 = m33;
}

template<class Type>
void Matrix4<Type>::set(const Matrix4<Type>& m1) {
  m00 = m1.m00;
  m01 = m1.m01;
  m02 = m1.m02;
  m03 = m1.m03;

  m10 = m1.m10;
  m11 = m1.m11;
  m12 = m1.m12;
  m13 = m1.m13;

  m20 = m1.m20;
  m21 = m1.m21;
  m22 = m1.m22;
  m23 = m1.m23;

  m30 = m1.m30;
  m31 = m1.m31;
  m32 = m1.m32;
  m33 = m1.m33;
}

template<class Type>
Matrix4<Type>::Matrix4(const Quaternion<Type>& q1) {
  Matrix4<Type>::set(q1);
}

template<class Type>
void Matrix4<Type>::set(const Type v[]) {
  m00 = v[0];
  m01 = v[1];
  m02 = v[2];
  m03 = v[3];

  m10 = v[4];
  m11 = v[5];
  m12 = v[6];
  m13 = v[7];

  m20 = v[8];
  m21 = v[9];
  m22 = v[10];
  m23 = v[11];

  m30 = v[12];
  m31 = v[13];
  m32 = v[14];
  m33 = v[15];
}

template<class Type>
void Matrix4<Type>::set(const Quaternion<Type>& q1) {
  Type norm = q1.norm();
  Type s = (Type)0.0;
  if (norm > 0.0) {
    s = (Type)2.0 / norm;
  }
  Type xs = q1.getX()*s;
  Type ys = q1.getY()*s;
  Type zs = q1.getZ()*s;

  Type xx = q1.getX()*xs;
  Type yy = q1.getY()*ys;
  Type zz = q1.getZ()*zs;
  Type xy = q1.getX()*ys;
  Type xz = q1.getX()*zs;
  Type xw = xs*q1.getW();
  Type yz = q1.getY()*zs;
  Type yw = ys*q1.getW();
  Type zw = zs*q1.getW();
  set((Type)1.0-(yy+zz), (xy-zw), (xz+yw), (Type)0.0,
    (xy+zw), (Type)1.0-(xx+zz), (yz-xw), (Type)0.0,
    (xz-yw), (yz+xw), (Type)1.0-(xx+yy), (Type)0.0,
    (Type)0.0, (Type)0.0, (Type)0.0, (Type)1.0);
}

template<class Type>
void Matrix4<Type>::setIdentity() {
  m00 = 1.0;
  m01 = 0.0;
  m02 = 0.0;
  m03 = 0.0;

  m10 = 0.0;
  m11 = 1.0;
  m12 = 0.0;
  m13 = 0.0;

  m20 = 0.0;
  m21 = 0.0;
  m22 = 1.0;
  m23 = 0.0;

  m30 = 0.0;
  m31 = 0.0;
  m32 = 0.0;
  m33 = 1.0;
}

template<class Type>
void Matrix4<Type>::setZero() {
  m00 = 0.0;
  m01 = 0.0;
  m02 = 0.0;
  m03 = 0.0;

  m10 = 0.0;
  m11 = 0.0;
  m12 = 0.0;
  m13 = 0.0;

  m20 = 0.0;
  m21 = 0.0;
  m22 = 0.0;
  m23 = 0.0;

  m30 = 0.0;
  m31 = 0.0;
  m32 = 0.0;
  m33 = 0.0;
}

template<class Type>
bool Matrix4<Type>::operator==(const Matrix4<Type>& m1) const {
  return equals(m00, m1.m00)
         && equals(m01, m1.m01)
         && equals(m02, m1.m02)
         && equals(m03, m1.m03)
         && equals(m10, m1.m10)
         && equals(m11, m1.m11)
         && equals(m12, m1.m12)
         && equals(m13, m1.m13)
         && equals(m20, m1.m20)
         && equals(m21, m1.m21)
         && equals(m22, m1.m22)
         && equals(m23, m1.m23)
         && equals(m30, m1.m30)
         && equals(m31, m1.m31)
         && equals(m32, m1.m32)
         && equals(m33, m1.m33);
}

template<class Type>
void Matrix4<Type>::add(Type scalar) {
  m00 += scalar;
  m01 += scalar;
  m02 += scalar;
  m03 += scalar;

  m10 += scalar;
  m11 += scalar;
  m12 += scalar;
  m13 += scalar;

  m20 += scalar;
  m21 += scalar;
  m22 += scalar;
  m23 += scalar;

  m30 += scalar;
  m31 += scalar;
  m32 += scalar;
  m33 += scalar;
}

template<class Type>
Matrix4<Type>&Matrix4<Type>::operator+=(Type scalar) {
  add(scalar);
  return *this;
}

template<class Type>
void Matrix4<Type>::add(Type scalar, const Matrix4<Type>& m1) {
  set(m1);
  add(scalar);
}

template<class Type>
Matrix4<Type> operator+(Type scalar, const Matrix4<Type>& m1) {
  Matrix4<Type> result(m1);
  result.add(scalar);
  return result;
}

template<class Type>
Matrix4<Type> operator+(const Matrix4<Type>& m1, Type scalar) {
  Matrix4<Type> result(m1);
  result.add(scalar);
  return result;
}

template<class Type>
void Matrix4<Type>::add(const Matrix4<Type>& m1) {
  m00 += m1.m00;
  m01 += m1.m01;
  m02 += m1.m02;
  m03 += m1.m03;

  m10 += m1.m10;
  m11 += m1.m11;
  m12 += m1.m12;
  m13 += m1.m13;

  m20 += m1.m20;
  m21 += m1.m21;
  m22 += m1.m22;
  m23 += m1.m23;

  m30 += m1.m30;
  m31 += m1.m31;
  m32 += m1.m32;
  m33 += m1.m33;
}

template<class Type>
Matrix4<Type>&Matrix4<Type>::operator+=(const Matrix4<Type>& m1) {
  add(m1);
  return *this;
}

template<class Type>
void Matrix4<Type>::add(const Matrix4<Type>& m1, const Matrix4<Type>& m2) {
  set(m1);
  add(m2);
}

template<class Type>
Matrix4<Type> Matrix4<Type>::operator+(const Matrix4<Type>& m1) const {
  Matrix4<Type> result(*this);
  result.add(m1);
  return result;
}

template<class Type>
void Matrix4<Type>::subtract(Type scalar) {
  add(-scalar);
}

template<class Type>
Matrix4<Type>&Matrix4<Type>::operator-=(Type scalar) {
  subtract(scalar);
  return *this;
}

template<class Type>
void Matrix4<Type>::subtract(Type scalar, const Matrix4<Type>& m1) {
  set(m1);
  subtract(scalar);
}

template<class Type>
Matrix4<Type> operator-(Type scalar, const Matrix4<Type>& m1) {
  Matrix4<Type> result(m1);
  result.subtract(scalar);
  return result;
}

template<class Type>
Matrix4<Type> operator-(const Matrix4<Type>& m1, Type scalar) {
  Matrix4<Type> result(m1);
  result.subtract(scalar);
  return result;
}

template<class Type>
void Matrix4<Type>::subtract(const Matrix4<Type>& m1) {
  m00 -= m1.m00;
  m01 -= m1.m01;
  m02 -= m1.m02;
  m03 -= m1.m03;

  m10 -= m1.m10;
  m11 -= m1.m11;
  m12 -= m1.m12;
  m13 -= m1.m13;

  m20 -= m1.m20;
  m21 -= m1.m21;
  m22 -= m1.m22;
  m23 -= m1.m23;

  m30 -= m1.m30;
  m31 -= m1.m31;
  m32 -= m1.m32;
  m33 -= m1.m33;
}

template<class Type>
Matrix4<Type>&Matrix4<Type>::operator-=(const Matrix4<Type>& m1) {
  subtract(m1);
  return *this;
}

template<class Type>
void Matrix4<Type>::subtract(const Matrix4<Type>& m1, const Matrix4<Type>& m2) {
  set(m1);
  subtract(m2);
}

template<class Type>
Matrix4<Type> Matrix4<Type>::operator-(const Matrix4<Type>& m1) const {
  Matrix4<Type> result(*this);
  result.subtract(m1);
  return result;
}

template<class Type>
void Matrix4<Type>::transpose(const Matrix4<Type>& m1) {
  m00 = m1.m00;
  m01 = m1.m10;
  m02 = m1.m20;
  m03 = m1.m30;

  m10 = m1.m01;
  m11 = m1.m11;
  m12 = m1.m21;
  m13 = m1.m31;

  m20 = m1.m02;
  m21 = m1.m12;
  m22 = m1.m22;
  m23 = m1.m32;

  m30 = m1.m03;
  m31 = m1.m13;
  m32 = m1.m23;
  m33 = m1.m33;
}

template<class Type>
void Matrix4<Type>::transpose() {
  transpose(Matrix4<Type>(*this));
}

template<class Type>
void Matrix4<Type>::getInColumnMajorOrder(Type* matrix) const {
  matrix[0] = m00;
  matrix[1] = m10;
  matrix[2] = m20;
  matrix[3] = m30;

  matrix[4] = m01;
  matrix[5] = m11;
  matrix[6] = m21;
  matrix[7] = m31;

  matrix[8] = m02;
  matrix[9] = m12;
  matrix[10] = m22;
  matrix[11] = m32;

  matrix[12] = m03;
  matrix[13] = m13;
  matrix[14] = m23;
  matrix[15] = m33;
}

template<class Type>
Type Matrix4<Type>::get00() const {
  return m00;
}

template<class Type>
Type Matrix4<Type>::get01() const {
  return m01;
}

template<class Type>
Type Matrix4<Type>::get02() const {
  return m02;
}

template<class Type>
Type Matrix4<Type>::get03() const {
  return m03;
}

template<class Type>
Type Matrix4<Type>::get10() const {
  return m10;
}

template<class Type>
Type Matrix4<Type>::get11() const {
  return m11;
}

template<class Type>
Type Matrix4<Type>::get12() const {
  return m12;
}

template<class Type>
Type Matrix4<Type>::get13() const {
  return m13;
}

template<class Type>
Type Matrix4<Type>::get20() const {
  return m20;
}

template<class Type>
Type Matrix4<Type>::get21() const {
  return m21;
}

template<class Type>
Type Matrix4<Type>::get22() const {
  return m22;
}

template<class Type>
Type Matrix4<Type>::get23() const {
  return m23;
}

template<class Type>
Type Matrix4<Type>::get30() const {
  return m30;
}

template<class Type>
Type Matrix4<Type>::get31() const {
  return m31;
}

template<class Type>
Type Matrix4<Type>::get32() const {
  return m32;
}

template<class Type>
Type Matrix4<Type>::get33() const {
  return m33;
}

template<class Type>
void Matrix4<Type>::set00(const Type& value) {
  m00 = value;
}

template<class Type>
void Matrix4<Type>::set01(const Type& value) {
  m01 = value;
}

template<class Type>
void Matrix4<Type>::set02(const Type& value) {
  m02 = value;
}

template<class Type>
void Matrix4<Type>::set03(const Type& value) {
  m03 = value;
}

template<class Type>
void Matrix4<Type>::set10(const Type& value) {
  m10 = value;
}

template<class Type>
void Matrix4<Type>::set11(const Type& value) {
  m11 = value;
}

template<class Type>
void Matrix4<Type>::set12(const Type& value) {
  m12 = value;
}

template<class Type>
void Matrix4<Type>::set13(const Type& value) {
  m13 = value;
}

template<class Type>
void Matrix4<Type>::set20(const Type& value) {
  m20 = value;
}

template<class Type>
void Matrix4<Type>::set21(const Type& value) {
  m21 = value;
}

template<class Type>
void Matrix4<Type>::set22(const Type& value) {
  m22 = value;
}

template<class Type>
void Matrix4<Type>::set23(const Type& value) {
  m23 = value;
}

template<class Type>
void Matrix4<Type>::set30(const Type& value) {
  m30 = value;
}

template<class Type>
void Matrix4<Type>::set31(const Type& value) {
  m31 = value;
}

template<class Type>
void Matrix4<Type>::set32(const Type& value) {
  m32 = value;
}

template<class Type>
void Matrix4<Type>::set33(const Type& value) {
  m33 = value;
}

template<class Type>
Type Matrix4<Type>::getElement(int row, int column) const {
  Type result = 0.0;
  if (row == 0) {
    if (column == 0) {
      result = m00;
    } else if (column == 1) {
      result = m01;
    } else if (column == 2) {
      result = m02;
    } else if (column == 3) {
      result = m03;
    } else {
      throw IndexOutOfBoundsException();
    }
  } else if (row == 1) {
    if (column == 0) {
      result = m10;
    } else if (column == 1) {
      result = m11;
    } else if (column == 2) {
      result = m12;
    } else if (column == 3) {
      result = m13;
    } else {
      throw IndexOutOfBoundsException();
    }
  } else if (row == 2) {
    if (column == 0) {
      result = m20;
    } else if (column == 1) {
      result = m21;
    } else if (column == 2) {
      result = m22;
    } else if (column == 3) {
      result = m23;
    } else {
      throw IndexOutOfBoundsException();
    }
  } else if (row == 3) {
    if (column == 0) {
      result = m30;
    } else if (column == 1) {
      result = m31;
    } else if (column == 2) {
      result = m32;
    } else if (column == 3) {
      result = m33;
    } else {
      throw IndexOutOfBoundsException();
    }
  } else {
    throw IndexOutOfBoundsException();
  }

  return result;
}

template<class Type>
void Matrix4<Type>::getRow(int row, Vector4<Type>& v) const {
  if (row == 0) {
    v.set(m00, m01, m02, m03);
  } else if (row == 1) {
    v.set(m10, m11, m12, m13);
  } else if (row == 2) {
    v.set(m20, m21, m22, m23);
  } else if (row == 3) {
    v.set(m30, m31, m32, m33);
  } else {
    throw IndexOutOfBoundsException();
  }
}

template<class Type>
void Matrix4<Type>::getColumn(int column, Vector4<Type>& v) const {
  if (column == 0) {
    v.set(m00, m10, m20, m30);
  } else if (column == 1) {
    v.set(m01, m11, m21, m31);
  } else if (column == 2) {
    v.set(m02, m12, m22, m32);
  } else if (column == 3) {
    v.set(m03, m13, m23, m33);
  } else {
    throw IndexOutOfBoundsException();
  }
}

template<class Type>
void Matrix4<Type>::setElement(int row, int column, Type value) {
  if (row == 0) {
    if (column == 0) {
      m00 = value;
    } else if (column == 1) {
      m01 = value;
    } else if (column == 2) {
      m02 = value;
    } else if (column == 3) {
      m03 = value;
    } else {
      throw IndexOutOfBoundsException();
    }
  } else if (row == 1) {
    if (column == 0) {
      m10 = value;
    } else if (column == 1) {
      m11 = value;
    } else if (column == 2) {
      m12 = value;
    } else if (column == 3) {
      m13 = value;
    } else {
      throw IndexOutOfBoundsException();
    }
  } else if (row == 2) {
    if (column == 0) {
      m20 = value;
    } else if (column == 1) {
      m21 = value;
    } else if (column == 2) {
      m22 = value;
    } else if (column == 3) {
      m23 = value;
    } else {
      throw IndexOutOfBoundsException();
    }
  } else if (row == 3) {
    if (column == 0) {
      m30 = value;
    } else if (column == 1) {
      m31 = value;
    } else if (column == 2) {
      m32 = value;
    } else if (column == 3) {
      m33 = value;
    } else {
      throw IndexOutOfBoundsException();
    }
  } else {
    throw IndexOutOfBoundsException();
  }
}

template<class Type>
void Matrix4<Type>::setColumn(int column, const Vector4<Type>& v) {
  if (column == 0) {
    m00 = v.getX();
    m10 = v.getY();
    m20 = v.getZ();
    m30 = v.getW();
  } else if (column == 1) {
    m01 = v.getX();
    m11 = v.getY();
    m21 = v.getZ();
    m31 = v.getW();
  } else if (column == 2) {
    m02 = v.getX();
    m12 = v.getY();
    m22 = v.getZ();
    m32 = v.getW();
  } else if (column == 3) {
    m03 = v.getX();
    m13 = v.getY();
    m23 = v.getZ();
    m33 = v.getW();
  } else {
    throw IndexOutOfBoundsException();
  }
}

template<class Type>
void Matrix4<Type>::setRow(int row, const Vector4<Type>& v) {
  if (row == 0) {
    m00 = v.getX();
    m01 = v.getY();
    m02 = v.getZ();
    m03 = v.getW();
  } else if (row == 1) {
    m10 = v.getX();
    m11 = v.getY();
    m12 = v.getZ();
    m13 = v.getW();
  } else if (row == 2) {
    m20 = v.getX();
    m21 = v.getY();
    m22 = v.getZ();
    m23 = v.getW();
  } else if (row == 3) {
    m30 = v.getX();
    m31 = v.getY();
    m32 = v.getZ();
    m33 = v.getW();
  } else {
    throw IndexOutOfBoundsException();
  }
}

template<class Type>
Type Matrix4<Type>::determinant() const {
  return
         (m00*m11 - m01*m10)*(m22*m33 - m23*m32)
         -(m00*m12 - m02*m10)*(m21*m33 - m23*m31)
         +(m00*m13 - m03*m10)*(m21*m32 - m22*m31)
         +(m01*m12 - m02*m11)*(m20*m33 - m23*m30)
         -(m01*m13 - m03*m11)*(m20*m32 - m22*m30)
         +(m02*m13 - m03*m12)*(m20*m31 - m21*m30);
}

template<class Type>
void Matrix4<Type>::invert() {
  Type s = determinant();
  if (s == 0.0) {
    return;
  }
  set(m11*(m22*m33 - m23*m32) + m12*(m23*m31 - m21*m33) + m13*(m21*m32 - m22*m31),
    m21*(m02*m33 - m03*m32) + m22*(m03*m31 - m01*m33) + m23*(m01*m32 - m02*m31),
    m31*(m02*m13 - m03*m12) + m32*(m03*m11 - m01*m13) + m33*(m01*m12 - m02*m11),
    m01*(m13*m22 - m12*m23) + m02*(m11*m23 - m13*m21) + m03*(m12*m21 - m11*m22),

    m12*(m20*m33 - m23*m30) + m13*(m22*m30 - m20*m32) + m10*(m23*m32 - m22*m33),
    m22*(m00*m33 - m03*m30) + m23*(m02*m30 - m00*m32) + m20*(m03*m32 - m02*m33),
    m32*(m00*m13 - m03*m10) + m33*(m02*m10 - m00*m12) + m30*(m03*m12 - m02*m13),
    m02*(m13*m20 - m10*m23) + m03*(m10*m22 - m12*m20) + m00*(m12*m23 - m13*m22),

    m13*(m20*m31 - m21*m30) + m10*(m21*m33 - m23*m31) + m11*(m23*m30 - m20*m33),
    m23*(m00*m31 - m01*m30) + m20*(m01*m33 - m03*m31) + m21*(m03*m30 - m00*m33),
    m33*(m00*m11 - m01*m10) + m30*(m01*m13 - m03*m11) + m31*(m03*m10 - m00*m13),
    m03*(m11*m20 - m10*m21) + m00*(m13*m21 - m11*m23) + m01*(m10*m23 - m13*m20),

    m10*(m22*m31 - m21*m32) + m11*(m20*m32 - m22*m30) + m12*(m21*m30 - m20*m31),
    m20*(m02*m31 - m01*m32) + m21*(m00*m32 - m02*m30) + m22*(m01*m30 - m00*m31),
    m30*(m02*m11 - m01*m12) + m31*(m00*m12 - m02*m10) + m32*(m01*m10 - m00*m11),
    m00*(m11*m22 - m12*m21) + m01*(m12*m20 - m10*m22) + m02*(m10*m21 - m11*m20));

  multiply((Type)1.0/s);
}

template<class Type>
void Matrix4<Type>::invert(const Matrix4<Type>& m1) {
  set(m1);
  invert();
}

template<class Type>
void Matrix4<Type>::multiply(Type scalar) {
  m00 *= scalar;
  m01 *= scalar;
  m02 *= scalar;
  m03 *= scalar;

  m10 *= scalar;
  m11 *= scalar;
  m12 *= scalar;
  m13 *= scalar;

  m20 *= scalar;
  m21 *= scalar;
  m22 *= scalar;
  m23 *= scalar;

  m30 *= scalar;
  m31 *= scalar;
  m32 *= scalar;
  m33 *= scalar;
}

template<class Type>
void Matrix4<Type>::multiply(Type scalar, const Matrix4<Type>& m1) {
  set(m1);
  multiply(scalar);
}

template<class Type>
void Matrix4<Type>::multiply(const Matrix4<Type>& m1) {
  set(m00*m1.m00 + m01*m1.m10 + m02*m1.m20 + m03*m1.m30,
    m00*m1.m01 + m01*m1.m11 + m02*m1.m21 + m03*m1.m31,
    m00*m1.m02 + m01*m1.m12 + m02*m1.m22 + m03*m1.m32,
    m00*m1.m03 + m01*m1.m13 + m02*m1.m23 + m03*m1.m33,

    m10*m1.m00 + m11*m1.m10 + m12*m1.m20 + m13*m1.m30,
    m10*m1.m01 + m11*m1.m11 + m12*m1.m21 + m13*m1.m31,
    m10*m1.m02 + m11*m1.m12 + m12*m1.m22 + m13*m1.m32,
    m10*m1.m03 + m11*m1.m13 + m12*m1.m23 + m13*m1.m33,

    m20*m1.m00 + m21*m1.m10 + m22*m1.m20 + m23*m1.m30,
    m20*m1.m01 + m21*m1.m11 + m22*m1.m21 + m23*m1.m31,
    m20*m1.m02 + m21*m1.m12 + m22*m1.m22 + m23*m1.m32,
    m20*m1.m03 + m21*m1.m13 + m22*m1.m23 + m23*m1.m33,

    m30*m1.m00 + m31*m1.m10 + m32*m1.m20 + m33*m1.m30,
    m30*m1.m01 + m31*m1.m11 + m32*m1.m21 + m33*m1.m31,
    m30*m1.m02 + m31*m1.m12 + m32*m1.m22 + m33*m1.m32,
    m30*m1.m03 + m31*m1.m13 + m32*m1.m23 + m33*m1.m33);
}

template<class Type>
void Matrix4<Type>::multiply(const Matrix4<Type>& m1, const Matrix4<Type>& m2) {
  set(m1);
  multiply(m2);
}

template<class Type>
Matrix4<Type> Matrix4<Type>::operator*(const Matrix4<Type>& m2) const {
  Matrix4<Type> result(*this);
  result.multiply(m2);
  return result;
}

template<class Type>
bool Matrix4<Type>::epsilonEquals(const Matrix4<Type>& m1, Type epsilonValue) const {
  Matrix4<Type> diff(*this);
  diff.subtract(m1);
  diff.absolute();
  bool equalM00 = (diff.m00 < epsilonValue || std::abs(diff.m00 - epsilonValue) <= std::numeric_limits<Type>::epsilon());
  bool equalM01 = (diff.m01 < epsilonValue || std::abs(diff.m01 - epsilonValue) <= std::numeric_limits<Type>::epsilon());
  bool equalM02 = (diff.m02 < epsilonValue || std::abs(diff.m02 - epsilonValue) <= std::numeric_limits<Type>::epsilon());
  bool equalM03 = (diff.m03 < epsilonValue || std::abs(diff.m03 - epsilonValue) <= std::numeric_limits<Type>::epsilon());

  bool equalM10 = (diff.m10 < epsilonValue || std::abs(diff.m10 - epsilonValue) <= std::numeric_limits<Type>::epsilon());
  bool equalM11 = (diff.m11 < epsilonValue || std::abs(diff.m11 - epsilonValue) <= std::numeric_limits<Type>::epsilon());
  bool equalM12 = (diff.m12 < epsilonValue || std::abs(diff.m12 - epsilonValue) <= std::numeric_limits<Type>::epsilon());
  bool equalM13 = (diff.m13 < epsilonValue || std::abs(diff.m13 - epsilonValue) <= std::numeric_limits<Type>::epsilon());

  bool equalM20 = (diff.m20 < epsilonValue || std::abs(diff.m20 - epsilonValue) <= std::numeric_limits<Type>::epsilon());
  bool equalM21 = (diff.m21 < epsilonValue || std::abs(diff.m21 - epsilonValue) <= std::numeric_limits<Type>::epsilon());
  bool equalM22 = (diff.m22 < epsilonValue || std::abs(diff.m22 - epsilonValue) <= std::numeric_limits<Type>::epsilon());
  bool equalM23 = (diff.m23 < epsilonValue || std::abs(diff.m23 - epsilonValue) <= std::numeric_limits<Type>::epsilon());

  bool equalM30 = (diff.m30 < epsilonValue || std::abs(diff.m30 - epsilonValue) <= std::numeric_limits<Type>::epsilon());
  bool equalM31 = (diff.m31 < epsilonValue || std::abs(diff.m31 - epsilonValue) <= std::numeric_limits<Type>::epsilon());
  bool equalM32 = (diff.m32 < epsilonValue || std::abs(diff.m32 - epsilonValue) <= std::numeric_limits<Type>::epsilon());
  bool equalM33 = (diff.m33 < epsilonValue || std::abs(diff.m33 - epsilonValue) <= std::numeric_limits<Type>::epsilon());

  return (equalM00 && equalM01 && equalM02 && equalM03
          && equalM10 && equalM11 && equalM12 && equalM13
          && equalM20 && equalM21 && equalM22 && equalM23
          && equalM30 && equalM31 && equalM32 && equalM33);
}

template<class Type>
void Matrix4<Type>::translation(const Vector3<Type>& vector) {
  m00 = 1.0;
  m01 = 0.0;
  m02 = 0.0;
  m03 = vector.x;

  m10 = 0.0;
  m11 = 1.0;
  m12 = 0.0;
  m13 = vector.y;

  m20 = 0.0;
  m21 = 0.0;
  m22 = 1.0;
  m23 = vector.z;

  m30 = 0.0;
  m31 = 0.0;
  m32 = 0.0;
  m33 = 1.0;
}

template<class Type>
void Matrix4<Type>::translation(Type x, Type y, Type z) {
  m00 = 1.0;
  m01 = 0.0;
  m02 = 0.0;
  m03 = x;

  m10 = 0.0;
  m11 = 1.0;
  m12 = 0.0;
  m13 = y;

  m20 = 0.0;
  m21 = 0.0;
  m22 = 1.0;
  m23 = z;

  m30 = 0.0;
  m31 = 0.0;
  m32 = 0.0;
  m33 = 1.0;
}

template<class Type>
void Matrix4<Type>::scaling(const Vector3<Type>& vector) {
  m00 = vector.x;
  m01 = 0.0;
  m02 = 0.0;
  m03 = 0.0;

  m10 = 0.0;
  m11 = vector.y;
  m12 = 0.0;
  m13 = 0.0;

  m20 = 0.0;
  m21 = 0.0;
  m22 = vector.z;
  m23 = 0.0;

  m30 = 0.0;
  m31 = 0.0;
  m32 = 0.0;
  m33 = 1.0;
}

template<class Type>
void Matrix4<Type>::scaling(Type x, Type y, Type z) {
  m00 = x;
  m01 = 0.0;
  m02 = 0.0;
  m03 = 0.0;

  m10 = 0.0;
  m11 = y;
  m12 = 0.0;
  m13 = 0.0;

  m20 = 0.0;
  m21 = 0.0;
  m22 = z;
  m23 = 0.0;

  m30 = 0.0;
  m31 = 0.0;
  m32 = 0.0;
  m33 = 1.0;
}

template<class Type>
void Matrix4<Type>::rotationX(Type angle) {
  Type c = std::cos(angle);
  Type s = std::sin(angle);
  m00 = 1.0;
  m01 = 0.0;
  m02 = 0.0;
  m03 = 0.0;

  m10 = 0.0;
  m11 = c;
  m12 = -s;
  m13 = 0.0;

  m20 = 0.0;
  m21 = s;
  m22 = c;
  m23 = 0.0;

  m30 = 0.0;
  m31 = 0.0;
  m32 = 0.0;
  m33 = 1.0;
}

template<class Type>
void Matrix4<Type>::rotationY(Type angle) {
  Type c = std::cos(angle);
  Type s = std::sin(angle);
  m00 = c;
  m01 = 0.0;
  m02 = s;
  m03 = 0.0;

  m10 = 0.0;
  m11 = 1.0;
  m12 = 0.0;
  m13 = 0.0;

  m20 = -s;
  m21 = 0.0;
  m22 = c;
  m23 = 0.0;

  m30 = 0.0;
  m31 = 0.0;
  m32 = 0.0;
  m33 = 1.0;

}

template<class Type>
void Matrix4<Type>::rotationZ(Type angle) {
  Type c = std::cos(angle);
  Type s = std::sin(angle);
  m00 = c;
  m01 = -s;
  m02 = 0.0;
  m03 = 0.0;

  m10 = s;
  m11 = c;
  m12 = 0.0;
  m13 = 0.0;

  m20 = 0.0;
  m21 = 0.0;
  m22 = 1.0;
  m23 = 0.0;

  m30 = 0.0;
  m31 = 0.0;
  m32 = 0.0;
  m33 = 1.0;
}

template<class Type>
void Matrix4<Type>::transform(const Tuple4<Type>& vec, Tuple4<Type>& vecOut) const {
  vecOut.set(m00*vec.getX() + m01*vec.getY() + m02*vec.getZ() + m03*vec.getW(),
    m10*vec.getX() + m11*vec.getY() + m12*vec.getZ() + m13*vec.getW(),
    m20*vec.getX() + m21*vec.getY() + m22*vec.getZ() + m23*vec.getW(),
    m30*vec.getX() + m31*vec.getY() + m32*vec.getZ() + m33*vec.getW());
}

template<class Type>
void Matrix4<Type>::transform(Tuple4<Type>& vec) const {
  Type x = vec.getX();
  Type y = vec.getY();
  Type z = vec.getZ();
  Type w = vec.getW();
  vec.set(m00*x + m01*y + m02*z + m03*w,
    m10*x + m11*y + m12*z + m13*w,
    m20*x + m21*y + m22*z + m23*w,
    m30*x + m31*y + m32*z + m33*w);
}

template<class Type>
void Matrix4<Type>::transform(const Tuple3<Type>& vec, Tuple3<Type>& vecOut) const {
  Type w = m30*vec.getX() + m31* vec.getY() + m32* vec.getZ() + m33;
  vecOut.set((m00*vec.getX() + m01*vec.getY() + m02*vec.getZ() + m03) / w,
    (m10*vec.getX() + m11*vec.getY() + m12*vec.getZ() + m13) / w,
    (m20*vec.getX() + m21*vec.getY() + m22*vec.getZ() + m23) / w);
}

template<class Type>
void Matrix4<Type>::transform(Tuple3<Type>& vec) const {
  Type x = vec.getX();
  Type y = vec.getY();
  Type z = vec.getZ();
  Type w = m30*x + m31*y + m32*z + m33;
  vec.set((m00*x + m01*y + m02*z + m03) / w,
    (m10*x + m11*y + m12*z + m13*w) / w,
    (m20*x + m21*y + m22*z + m23*w) / w);
}

template<class Type>
void Matrix4<Type>::negate() {
  m00 = -m00;
  m01 = -m01;
  m02 = -m02;
  m03 = -m03;

  m10 = -m10;
  m11 = -m11;
  m12 = -m12;
  m13 = -m13;

  m20 = -m20;
  m21 = -m21;
  m22 = -m22;
  m23 = -m23;

  m30 = -m30;
  m31 = -m31;
  m32 = -m32;
  m33 = -m33;
}

template<class Type>
void Matrix4<Type>::negate(const Matrix4<Type>& m1) {
  set(m1);
  negate();
}

template<class Type>
Matrix4<Type> Matrix4<Type>::operator-() const {
  return Matrix4<Type>(-m00, -m01, -m02, -m03,
                       -m10, -m11, -m12, -m13,
                       -m20, -m21, -m22, -m23,
                       -m30, -m31, -m32, -m33);
}

template<class Type>
void Matrix4<Type>::absolute() {
  m00 = std::abs(m00);
  m01 = std::abs(m01);
  m02 = std::abs(m02);
  m03 = std::abs(m03);

  m10 = std::abs(m10);
  m11 = std::abs(m11);
  m12 = std::abs(m12);
  m13 = std::abs(m13);

  m20 = std::abs(m20);
  m21 = std::abs(m21);
  m22 = std::abs(m22);
  m23 = std::abs(m23);

  m30 = std::abs(m30);
  m31 = std::abs(m31);
  m32 = std::abs(m32);
  m33 = std::abs(m33);
}

template<class Type>
std::string Matrix4<Type>::toString() const {
  std::string oss;
  oss = "[ [" + m00 + "," + m01 + ","
        + m02 + "," + m03 + "]" + '\n'
        + "  [" + m10 + "," + m11 + ","
        + m12 + "," + m13 + "]" + '\n'
        + "  [" + m20 + "," + m21 + ","
        + m22 + "," + m23 + "]" + '\n'
        + "  [" + m30 + "," + m31 + ","
        + m32 + "," + m33 + "] ]" + std::endl;
  return oss;
}

template<class Type>
std::ostream& operator<<(std::ostream& stream, const Matrix4<Type>& m1) {
  return stream << m1.toString();
}

}
}


#endif
