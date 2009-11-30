#ifndef mi_math_Matrix3_h
#define mi_math_Matrix3_h

#include <cmath>
#include <iostream>
#include <string>
#include <math/Tuple3.h>
#include <math/IndexOutOfBoundsException.h>
#include <math/Vector3.h>

namespace mi
{
    namespace math
    {
        template<class Type> class Matrix3;
    }
}

#include <math/Quaternion.h>

namespace mi
{
    namespace math
    {

/**
 * A 3 x 3 matrix of Types.
 */
        template<class Type>
        class Matrix3
        {
            friend class Matrix3Test;

            /**
             *  Whether to use right or left handed rotation when
             *  converting from a quaternion. Important note: this is only
             *  used for quaternion conversion; it is not a general
             *  setting.
             */
            static const bool useLeftHandedRotation;

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

        public:

            /**
             * Constructs and initializes a Matrix3 to all zeros.
             */
            Matrix3();

            /**
             * Constructs and initializes a Matrix3 from the specified 16 values.
             *
             * @param mm00 the [0][0] element
             * @param mm01 the [0][1] element
             * @param mm02 the [0][2] element
             * @param mm10 the [1][0] element
             * @param mm11 the [1][1] element
             * @param mm12 the [1][2] element
             * @param mm20 the [2][0] element
             * @param mm21 the [2][1] element
             * @param mm22 the [2][2] element
             */
            Matrix3(Type mm00, Type mm01, Type mm02,
                    Type mm10, Type mm11, Type mm12,
                    Type mm20, Type mm21, Type mm22);

            /**
             * Constructs and initializes a Matrix3 from the specified 16
             * element array.  this.m00=v[0], this.m01=v[1], etc.
             *
             * @param  v the array of length 16 containing in order
             */
            Matrix3(const Type v[]);

            /**
             * Sets 9 values
             *
             * @param mm00 the [0][0] element
             * @param mm01 the [0][1] element
             * @param mm02 the [0][2] element
             * @param mm10 the [1][0] element
             * @param mm11 the [1][1] element
             * @param mm12 the [1][2] element
             * @param mm20 the [2][0] element
             * @param mm21 the [2][1] element
             * @param mm22 the [2][2] element
             */
            void set(Type mm00, Type mm01, Type mm02,
                     Type mm10, Type mm11, Type mm12,
                     Type mm20, Type mm21, Type mm22);

            /**
             * Sets the values in this Matrix3 equal to the row-major
             * array parameter (ie, the first three elements of the array
             * will be copied into the first row of this matrix, etc.).
             */
            void set(const Type m[]);

            /**
             * Sets the value of this matrix to a copy of the passed
             * matrix m1.
             *
             * @param m1 the matrix to be copied.
             */
            void set(const Matrix3 &m1);

            /**
             * Sets the value of this matrix to the value derived from the
             * quaternion.
             *
             * @param q1  q quaternion.
             */
            void set(const Quaternion<Type> &q1);

            /**
             * Sets this Matrix3 to identity.
             */
            void setIdentity();

            /**
             * Sets this matrix to all zeros.
             */
            void setZero();

            /**
             * Returns true if all of the data members of Matrix3 m1 are
             * equal to the corresponding data members in this Matrix3.
             *
             * @param m1 The matrix with which the comparison is made.
             * @return true or false
             */
            bool operator==(const Matrix3 &m1) const;

            /**
             * Returns true if the difference of each of the elements of this
             * matrix and the corresponding element of matrix m1 are less than
             * or equal to the epsilon parameter.
             *
             * @param m1  the other matrix.
             * @param epsilon  the threshold value.
             */
            bool epsilonEquals(const Matrix3 &m1, Type epsilon) const;

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
            Matrix3&operator+=(Type scalar);

            /**
             * Adds a scalar to each component of the matrix m1 and places
             * the result into this. Matrix m1 is not modified.
             *
             * @param scalar  the value to add.
             * @parm m1  the original matrix values.
             */
            void add(Type scalar, const Matrix3 &m1);

            /**
             * Sets the value of this matrix to sum of itself and matrix
             * m1.
             *
             * @param m1  the other matrix.
             */
            void add(const Matrix3 &m1);

            /**
             * Sets the value of this matrix to the sum of itself and matrix m1.
             *
             * @param m1  the other matrix.
             */
            Matrix3&operator+=(const Matrix3 &m1);

            /**
             * Sets the value of this matrix to the matrix sum of matrices
             * m1 and m2.
             *
             * @param m1  the first matrix.
             * @param m2  the second matrix.
             */
            void add(const Matrix3 &m1, const Matrix3 &m2);

            /**
             * Returns a new Matrix3 with the value of the sum of this matrix
             * and matrix m1.
             *
             * @param m1  the other matrix.
             * @return  a new Matrix3.
             */
            Matrix3 operator+(const Matrix3 &m1) const;

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
            Matrix3&operator-=(Type scalar);

            /**
             * Subtracts a scalar to each component of the matrix m1 and places
             * the result into this. Matrix m1 is not modified.
             *
             * @param scalar  the value to subtract.
             * @parm m1  the original matrix values.
             */
            void subtract(Type scalar, const Matrix3 &m1);

            /**
             * Sets the value of this matrix to the matrix difference of
             * itself and matrix m1.
             *
             * @param m1  the other matrix
             */
            void subtract(const Matrix3 &m1);

            /**
             * Sets the value of this matrix to the difference of itself and
             * matrix m1.
             *
             * @param m1 the other matrix.
             */
            Matrix3&operator-=(const Matrix3 &m1);

            /**
             * Sets the value of this matrix to the matrix difference of
             * matrices m1 and m2.
             *
             * @param m1 the first matrix
             * @param m2 the second matrix
             */
            void subtract(const Matrix3 &m1, const Matrix3 &m2);

            /**
             * Returns a new Matrix3 with the value of the difference
             * of this matrix and m1.
             *
             * @param m1  the first matrix.
             * @return  a new Matrix3.
             */
            Matrix3 operator-(const Matrix3 &m1) const;

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
            void transpose(const Matrix3 &m1);

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
             * Sets the value of the first row and first column of this matrix.
             */
            void set00(const Type &value);

            /**
             * Sets the value of the first row and second column of this matrix.
             */
            void set01(const Type &value);

            /**
             * Sets the value of the first row and third column of this matrix.
             */
            void set02(const Type &value);

            /**
             * Sets the value of the second row and first column of this matrix.
             */
            void set10(const Type &value);

            /**
             * Sets the value of the second row and second column of this matrix.
             */
            void set11(const Type &value);

            /**
             * Sets the value of the second row and third column of this matrix.
             */
            void set12(const Type &value);

            /**
             * Sets the value of the third row and first column of this matrix.
             */
            void set20(const Type &value);

            /**
             * Sets the value of the third row and second column of this matrix.
             */
            void set21(const Type &value);

            /**
             * Sets the value of the third row and third column of this matrix.
             */
            void set22(const Type &value);

            /**
             * Sets the specified element of this matrix3 to the value provided.
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
             * Sets the specified row of this matrix3 to the Vector
             * provided.
             *
             * @param row the row number to be modified (zero indexed)
             * @param v the replacement row
             */
            void setRow(int row, const Vector3<Type> &v);

            /**
             * Copies the matrix values in the specified row into the
             * vector parameter.
             * @param row the matrix row
             * @param v The vector into which the matrix row values will be copied
             */
            void getRow(int row, Vector3<Type> &v) const;

            /**
             * Sets the specified column of this matrix3 to the vector
             * provided.
             *
             * @param column the column number to be modified (zero indexed).
             * @param v the replacement column.
             */
            void setColumn(int column, const Vector3<Type> &v);

            /**
             * Copies the matrix values in the specified column into the
             * vector parameter.
             *
             * @param column the matrix column
             * @param v The vector into which the matrix column values
             *   will be copied
             */
            void getColumn(int column, Vector3<Type> &v) const;

            /**
             * Sets the value of this matrix to the matrix inverse
             * of the passed matrix m1.
             *
             * @param m1 the matrix to be inverted
             */
            void inverse(const Matrix3 &m1);

            /**
             * Sets the value of this matrix to its inverse.
             */
            void inverse();

            /**
             * Computes the determinant of this matrix.
             *
             * @return the determinant of the matrix
             */
            Type determinant() const;

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
            void multiply(Type scalar, const Matrix3 &m1);

            /**
             * Sets the value of this matrix to the result of multiplying
             * itself with matrix m1.
             *
             * @param m1 the other matrix */
            void multiply(const Matrix3 &m1);

            /**
             * Sets the value of this matrix to the result of multiplying
             * the two argument matrices together.
             *
             * @param m1 the first matrix
             * @param m2 the second matrix
             */
            void multiply(const Matrix3 &m1, const Matrix3 &m2);

            /**
             * Transform the vector vec using this Matrix3 and place the
             * result into vecOut.
             *
             * @param vec the T precision vector to be transformed
             * @param vecOut the vector into which the transformed values
             *   are placed
             */
            void transform(const Tuple3<Type> &vec, Tuple3<Type> &vecOut) const;

            /**
             * Transform the vector vec using this Matrix3 and place the
             * result back into vec.
             *
             * @param vec the T precision vector to be transformed
             */
            void transform(Tuple3<Type> &vec) const;

            /**
             * Negates the value of this matrix: this = -this.
             */
            void negate();

            /**
             * Sets the value of this matrix equal to the negation of of
             * the Matrix3 parameter.
             *
             * @param m1 The source matrix
             */
            void negate(const Matrix3 &m1);

            /**
             * Returns a new Matrix3 with the value of the negation of
             * this matrix.
             *
             * @return a new Matrix3.
             */
            Matrix3 operator-() const;

            /**
             * Sets each element to its absolute value.
             */
            void absolute();

            /**
             * Returns a string that contains the values of this Matrix3.
             *
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
        std::ostream&operator<<(std::ostream &stream, const Matrix3<Type> &m1);

/**
 * Adds a scalar to each component of the matrix m1 and places
 * the result into a new matrix. Matrix m1 is not modified.
 *
 * @param scalar  the value to add.
 * @parm m1  the original matrix values.
 */
        template<class Type>
        Matrix3<Type> operator+(Type scalar, const Matrix3<Type> &m1);

/**
 * Adds a scalar to each component of the matrix m1 and places
 * the result into a new matrix. Matrix m1 is not modified.
 *
 * @param scalar  the value to add.
 * @parm m1  the original matrix values.
 */
        template<class Type>
        Matrix3<Type> operator+(const Matrix3<Type> &m1, Type scalar);

/**
 * Subtracts a scalar to each component of the matrix m1 and places
 * the result into a new matrix. Matrix m1 is not modified.
 *
 * @param scalar  the value to subtract.
 * @parm m1  the original matrix values.
 */
        template<class Type>
        Matrix3<Type> operator-(Type scalar, const Matrix3<Type> &m1);

/**
 * Subtracts a scalar to each component of the matrix m1 and places
 * the result into a new matrix. Matrix m1 is not modified.
 *
 * @param scalar  the value to subtract.
 * @parm m1  the original matrix values.
 */
        template<class Type>
        Matrix3<Type> operator-(const Matrix3<Type> &m1, Type scalar);


        template<class Type>
        const bool Matrix3<Type>::useLeftHandedRotation = false;

        template<class Type>
        Matrix3<Type>::Matrix3()
            : m00(0.0), m01(0.0), m02(0.0),
            m10(0.0), m11(0.0), m12(0.0),
            m20(0.0), m21(0.0), m22(0.0)
        {
        }

        template<class Type>
        Matrix3<Type>::Matrix3(Type mm00, Type mm01, Type mm02,
                               Type mm10, Type mm11, Type mm12,
                               Type mm20, Type mm21, Type mm22)
            : m00(mm00), m01(mm01), m02(mm02),
            m10(mm10), m11(mm11), m12(mm12),
            m20(mm20), m21(mm21), m22(mm22)
        {
        }

        template<class Type>
        Matrix3<Type>::Matrix3(const Type v[])
            : m00(v[0]), m01(v[1]), m02(v[2]),
            m10(v[3]), m11(v[4]), m12(v[5]),
            m20(v[6]), m21(v[7]), m22(v[8])
        {
        }

        template<class Type>
        void Matrix3<Type>::set(Type mm00, Type mm01, Type mm02,
                                Type mm10, Type mm11, Type mm12,
                                Type mm20, Type mm21, Type mm22)
        {
            this->m00 = mm00;
            this->m01 = mm01;
            this->m02 = mm02;

            this->m10 = mm10;
            this->m11 = mm11;
            this->m12 = mm12;

            this->m20 = mm20;
            this->m21 = mm21;
            this->m22 = mm22;
        }

        template<class Type>
        void Matrix3<Type>::set(const Matrix3<Type> &m1)
        {
            m00 = m1.m00;
            m01 = m1.m01;
            m02 = m1.m02;

            m10 = m1.m10;
            m11 = m1.m11;
            m12 = m1.m12;

            m20 = m1.m20;
            m21 = m1.m21;
            m22 = m1.m22;
        }

        template<class Type>
        void Matrix3<Type>::set(const Type v[])
        {
            m00 = v[0];
            m01 = v[1];
            m02 = v[2];

            m10 = v[3];
            m11 = v[5];
            m12 = v[6];

            m20 = v[8];
            m21 = v[9];
            m22 = v[10];
        }

        template<class Type>
        void Matrix3<Type>::set(const Quaternion<Type> &q1)
        {
            Type w;
            if (useLeftHandedRotation)
            {
                // Left handed rotation
                w = -q1.getW();
            }
            else
            {
                // Right handed rotation
                w = q1.getW();
            }

            Type norm = q1.norm();
            double s = 0.0;
            if (norm > 0.0)
            {
                s = 2.0 / norm;
            }
            Type xs = q1.getX()*s;
            Type ys = q1.getY()*s;
            Type zs = q1.getZ()*s;

            Type xx = q1.getX()*xs;
            Type yy = q1.getY()*ys;
            Type zz = q1.getZ()*zs;
            Type xy = q1.getX()*ys;
            Type xz = q1.getX()*zs;
            Type xw = xs*w;
            Type yz = q1.getY()*zs;
            Type yw = ys*w;
            Type zw = zs*w;
            set(1.0-(yy+zz), (xy-zw),     (xz+yw),
                (xy+zw),     1.0-(xx+zz), (yz-xw),
                (xz-yw),     (yz+xw),     1.0-(xx+yy));

            //                  Type w2 = w*w;
            //                  Type x2 = q1.getX()*q1.getX();
            //                  Type y2 = q1.getY()*q1.getY();
            //                  Type z2 = q1.getZ()*q1.getZ();

            //                  Type wx = w*q1.getX();
            //                  Type wy = w*q1.getY();
            //                  Type wz = w*q1.getZ();

            //                  Type xy = q1.getX()*q1.getY();
            //                  Type xz = q1.getX()*q1.getZ();

            //                  Type yz = q1.getY()*q1.getZ();

            //                  set(w2 + x2 - y2 - z2,  2.0 * (xy - wz),    2.0 * (xz + wy),
            //                          2.0 * (xy + wz),    w2 - x2 + y2 - z2,  2.0 * (yz - wx),
            //                          2.0 * (xz - wy),    2.0 * (yz + wx),    w2 - x2 - y2 + z2);

        }

        template<class Type>
        void Matrix3<Type>::setIdentity()
        {
            m00 = 1.0;
            m01 = 0.0;
            m02 = 0.0;

            m10 = 0.0;
            m11 = 1.0;
            m12 = 0.0;

            m20 = 0.0;
            m21 = 0.0;
            m22 = 1.0;
        }

        template<class Type>
        void Matrix3<Type>::setZero()
        {
            m00 = 0.0;
            m01 = 0.0;
            m02 = 0.0;

            m10 = 0.0;
            m11 = 0.0;
            m12 = 0.0;

            m20 = 0.0;
            m21 = 0.0;
            m22 = 0.0;
        }

        template<class Type>
        bool Matrix3<Type>::operator==(const Matrix3<Type> &m1) const
        {
            return equals(m00, m1.m00)
                   && equals(m01, m1.m01)
                   && equals(m02, m1.m02)
                   && equals(m10, m1.m10)
                   && equals(m11, m1.m11)
                   && equals(m12, m1.m12)
                   && equals(m20, m1.m20)
                   && equals(m21, m1.m21)
                   && equals(m22, m1.m22);
        }

        template<class Type>
        void Matrix3<Type>::add(Type scalar)
        {
            m00 += scalar;
            m01 += scalar;
            m02 += scalar;

            m10 += scalar;
            m11 += scalar;
            m12 += scalar;

            m20 += scalar;
            m21 += scalar;
            m22 += scalar;
        }

        template<class Type>
        Matrix3<Type>&Matrix3<Type>::operator+=(Type scalar)
        {
            add(scalar);
            return *this;
        }

        template<class Type>
        void Matrix3<Type>::add(Type scalar, const Matrix3<Type> &m1)
        {
            set(m1);
            add(scalar);
        }

        template<class Type>
        Matrix3<Type> operator+(Type scalar, const Matrix3<Type> &m1)
        {
            Matrix3<Type> result(m1);
            result.add(scalar);
            return result;
        }

        template<class Type>
        Matrix3<Type> operator+(const Matrix3<Type> &m1, Type scalar)
        {
            Matrix3<Type> result(m1);
            result.add(scalar);
            return result;
        }

        template<class Type>
        void Matrix3<Type>::add(const Matrix3<Type> &m1)
        {
            m00 += m1.m00;
            m01 += m1.m01;
            m02 += m1.m02;

            m10 += m1.m10;
            m11 += m1.m11;
            m12 += m1.m12;

            m20 += m1.m20;
            m21 += m1.m21;
            m22 += m1.m22;
        }

        template<class Type>
        Matrix3<Type>&Matrix3<Type>::operator+=(const Matrix3<Type> &m1)
        {
            add(m1);
            return *this;
        }

        template<class Type>
        void Matrix3<Type>::add(const Matrix3<Type> &m1, const Matrix3<Type> &m2)
        {
            set(m1);
            add(m2);
        }

        template<class Type>
        Matrix3<Type> Matrix3<Type>::operator+(const Matrix3<Type> &m1) const
        {
            Matrix3<Type> result(*this);
            result.add(m1);
            return result;
        }

        template<class Type>
        void Matrix3<Type>::subtract(Type scalar)
        {
            add(-scalar);
        }

        template<class Type>
        Matrix3<Type>&Matrix3<Type>::operator-=(Type scalar)
        {
            subtract(scalar);
            return *this;
        }

        template<class Type>
        void Matrix3<Type>::subtract(Type scalar, const Matrix3<Type> &m1)
        {
            set(m1);
            subtract(scalar);
        }

        template<class Type>
        Matrix3<Type> operator-(Type scalar, const Matrix3<Type> &m1)
        {
            Matrix3<Type> result(m1);
            result.subtract(scalar);
            return result;
        }

        template<class Type>
        Matrix3<Type> operator-(const Matrix3<Type> &m1, Type scalar)
        {
            Matrix3<Type> result(m1);
            result.subtract(scalar);
            return result;
        }

        template<class Type>
        void Matrix3<Type>::subtract(const Matrix3<Type> &m1)
        {
            m00 -= m1.m00;
            m01 -= m1.m01;
            m02 -= m1.m02;

            m10 -= m1.m10;
            m11 -= m1.m11;
            m12 -= m1.m12;

            m20 -= m1.m20;
            m21 -= m1.m21;
            m22 -= m1.m22;
        }

        template<class Type>
        Matrix3<Type>&Matrix3<Type>::operator-=(const Matrix3<Type> &m1)
        {
            subtract(m1);
            return *this;
        }

        template<class Type>
        void Matrix3<Type>::subtract(const Matrix3<Type> &m1, const Matrix3<Type> &m2)
        {
            set(m1);
            subtract(m2);
        }

        template<class Type>
        Matrix3<Type> Matrix3<Type>::operator-(const Matrix3<Type> &m1) const
        {
            Matrix3<Type> result(*this);
            result.subtract(m1);
            return result;
        }

        template<class Type>
        void Matrix3<Type>::transpose(const Matrix3<Type> &m1)
        {
            m00 = m1.m00;
            m01 = m1.m10;
            m02 = m1.m20;

            m10 = m1.m01;
            m11 = m1.m11;
            m12 = m1.m21;

            m20 = m1.m02;
            m21 = m1.m12;
            m22 = m1.m22;
        }

        template<class Type>
        void Matrix3<Type>::transpose()
        {
            transpose(Matrix3<Type>(*this));
        }

        template<class Type>
        Type Matrix3<Type>::get00() const
        {
            return m00;
        }

        template<class Type>
        Type Matrix3<Type>::get01() const
        {
            return m01;
        }

        template<class Type>
        Type Matrix3<Type>::get02() const
        {
            return m02;
        }

        template<class Type>
        Type Matrix3<Type>::get10() const
        {
            return m10;
        }

        template<class Type>
        Type Matrix3<Type>::get11() const
        {
            return m11;
        }

        template<class Type>
        Type Matrix3<Type>::get12() const
        {
            return m12;
        }

        template<class Type>
        Type Matrix3<Type>::get20() const
        {
            return m20;
        }

        template<class Type>
        Type Matrix3<Type>::get21() const
        {
            return m21;
        }

        template<class Type>
        Type Matrix3<Type>::get22() const
        {
            return m22;
        }

        template<class Type>
        void Matrix3<Type>::set00(const Type &value)
        {
            m00 = value;
        }

        template<class Type>
        void Matrix3<Type>::set01(const Type &value)
        {
            m01 = value;
        }

        template<class Type>
        void Matrix3<Type>::set02(const Type &value)
        {
            m02 = value;
        }

        template<class Type>
        void Matrix3<Type>::set10(const Type &value)
        {
            m10 = value;
        }

        template<class Type>
        void Matrix3<Type>::set11(const Type &value)
        {
            m11 = value;
        }

        template<class Type>
        void Matrix3<Type>::set12(const Type &value)
        {
            m12 = value;
        }

        template<class Type>
        void Matrix3<Type>::set20(const Type &value)
        {
            m20 = value;
        }

        template<class Type>
        void Matrix3<Type>::set21(const Type &value)
        {
            m21 = value;
        }

        template<class Type>
        void Matrix3<Type>::set22(const Type &value)
        {
            m22 = value;
        }

        template<class Type>
        Type Matrix3<Type>::getElement(int row, int column) const
        {
            Type result = 0.0;
            if (row == 0)
            {
                if (column == 0)
                {
                    result = m00;
                }
                else if (column == 1)
                {
                    result = m01;
                }
                else if (column == 2)
                {
                    result = m02;
                }
                else
                {
                    throw IndexOutOfBoundsException();
                }
            }
            else if (row == 1)
            {
                if (column == 0)
                {
                    result = m10;
                }
                else if (column == 1)
                {
                    result = m11;
                }
                else if (column == 2)
                {
                    result = m12;
                }
                else
                {
                    throw IndexOutOfBoundsException();
                }
            }
            else if (row == 2)
            {
                if (column == 0)
                {
                    result = m20;
                }
                else if (column == 1)
                {
                    result = m21;
                }
                else if (column == 2)
                {
                    result = m22;
                }
                else
                {
                    throw IndexOutOfBoundsException();
                }
            }
            else
            {
                throw IndexOutOfBoundsException();
            }

            return result;
        }

        template<class Type>
        void Matrix3<Type>::getRow(int row, Vector3<Type> &v) const
        {
            if (row == 0)
            {
                v.set(m00, m01, m02);
            }
            else if (row == 1)
            {
                v.set(m10, m11, m12);
            }
            else if (row == 2)
            {
                v.set(m20, m21, m22);
            }
            else
            {
                throw IndexOutOfBoundsException();
            }
        }

        template<class Type>
        void Matrix3<Type>::getColumn(int column, Vector3<Type> &v) const
        {
            if (column == 0)
            {
                v.set(m00, m10, m20);
            }
            else if (column == 1)
            {
                v.set(m01, m11, m21);
            }
            else if (column == 2)
            {
                v.set(m02, m12, m22);
            }
            else
            {
                throw IndexOutOfBoundsException();
            }
        }

        template<class Type>
        void Matrix3<Type>::setElement(int row, int column, Type value)
        {
            if (row == 0)
            {
                if (column == 0)
                {
                    m00 = value;
                }
                else if (column == 1)
                {
                    m01 = value;
                }
                else if (column == 2)
                {
                    m02 = value;
                }
                else
                {
                    throw IndexOutOfBoundsException();
                }
            }
            else if (row == 1)
            {
                if (column == 0)
                {
                    m10 = value;
                }
                else if (column == 1)
                {
                    m11 = value;
                }
                else if (column == 2)
                {
                    m12 = value;
                }
                else
                {
                    throw IndexOutOfBoundsException();
                }
            }
            else if (row == 2)
            {
                if (column == 0)
                {
                    m20 = value;
                }
                else if (column == 1)
                {
                    m21 = value;
                }
                else if (column == 2)
                {
                    m22 = value;
                }
                else
                {
                    throw IndexOutOfBoundsException();
                }
            }
            else
            {
                throw IndexOutOfBoundsException();
            }
        }

        template<class Type>
        void Matrix3<Type>::setColumn(int column, const Vector3<Type> &v)
        {
            if (column == 0)
            {
                m00 = v.getX();
                m10 = v.getY();
                m20 = v.getZ();
            }
            else if (column == 1)
            {
                m01 = v.getX();
                m11 = v.getY();
                m21 = v.getZ();
            }
            else if (column == 2)
            {
                m02 = v.getX();
                m12 = v.getY();
                m22 = v.getZ();
            }
            else
            {
                throw IndexOutOfBoundsException();
            }
        }

        template<class Type>
        void Matrix3<Type>::setRow(int row, const Vector3<Type> &v)
        {
            if (row == 0)
            {
                m00 = v.getX();
                m01 = v.getY();
                m02 = v.getZ();
            }
            else if (row == 1)
            {
                m10 = v.getX();
                m11 = v.getY();
                m12 = v.getZ();
            }
            else if (row == 2)
            {
                m20 = v.getX();
                m21 = v.getY();
                m22 = v.getZ();
            }
            else
            {
                throw IndexOutOfBoundsException();
            }
        }

        template<class Type>
        Type Matrix3<Type>::determinant() const
        {
            return
                m00*(m11*m22 - m21*m12)
                - m01*(m10*m22 - m02*m12)
                + m02*(m10*m21 - m20*m11);
        }

        template<class Type>
        void Matrix3<Type>::inverse()
        {
            Type s = determinant();
            if (s == 0.0)
            {
                return;
            }
            set(m11*m22 - m12*m21, m02*m21 - m01*m22, m01*m12 - m02*m11,
                m12*m20 - m10*m22, m00*m22 - m02*m20, m02*m10 - m00*m12,
                m10*m21 - m11*m20, m01*m20 - m00*m21, m00*m11 - m01*m10);
            multiply(1.0/s);
        }

        template<class Type>
        void Matrix3<Type>::inverse(const Matrix3<Type> &m1)
        {
            set(m1);
            inverse();
        }

        template<class Type>
        void Matrix3<Type>::multiply(Type scalar)
        {
            m00 *= scalar;
            m01 *= scalar;
            m02 *= scalar;

            m10 *= scalar;
            m11 *= scalar;
            m12 *= scalar;

            m20 *= scalar;
            m21 *= scalar;
            m22 *= scalar;
        }

        template<class Type>
        void Matrix3<Type>::multiply(Type scalar, const Matrix3<Type> &m1)
        {
            set(m1);
            multiply(scalar);
        }

        template<class Type>
        void Matrix3<Type>::multiply(const Matrix3<Type> &m1)
        {
            set(m00*m1.m00 + m01*m1.m10 + m02*m1.m20,
                m00*m1.m01 + m01*m1.m11 + m02*m1.m21,
                m00*m1.m02 + m01*m1.m12 + m02*m1.m22,

                m10*m1.m00 + m11*m1.m10 + m12*m1.m20,
                m10*m1.m01 + m11*m1.m11 + m12*m1.m21,
                m10*m1.m02 + m11*m1.m12 + m12*m1.m22,

                m20*m1.m00 + m21*m1.m10 + m22*m1.m20,
                m20*m1.m01 + m21*m1.m11 + m22*m1.m21,
                m20*m1.m02 + m21*m1.m12 + m22*m1.m22);
        }

        template<class Type>
        void Matrix3<Type>::multiply(const Matrix3<Type> &m1, const Matrix3<Type> &m2)
        {
            set(m1);
            multiply(m2);
        }

        template<class Type>
        bool Matrix3<Type>::epsilonEquals(const Matrix3<Type> &m1, Type epsilon) const
        {
            Matrix3<Type> diff(*this);
            diff.subtract(m1);
            diff.absolute();
            bool equalM00 = (diff.m00 < epsilon || std::abs(diff.m00 - epsilon) <= std::numeric_limits<Type>::epsilon());
            bool equalM01 = (diff.m01 < epsilon || std::abs(diff.m01 - epsilon) <= std::numeric_limits<Type>::epsilon());
            bool equalM02 = (diff.m02 < epsilon || std::abs(diff.m02 - epsilon) <= std::numeric_limits<Type>::epsilon());

            bool equalM10 = (diff.m10 < epsilon || std::abs(diff.m10 - epsilon) <= std::numeric_limits<Type>::epsilon());
            bool equalM11 = (diff.m11 < epsilon || std::abs(diff.m11 - epsilon) <= std::numeric_limits<Type>::epsilon());
            bool equalM12 = (diff.m12 < epsilon || std::abs(diff.m12 - epsilon) <= std::numeric_limits<Type>::epsilon());

            bool equalM20 = (diff.m20 < epsilon || std::abs(diff.m20 - epsilon) <= std::numeric_limits<Type>::epsilon());
            bool equalM21 = (diff.m21 < epsilon || std::abs(diff.m21 - epsilon) <= std::numeric_limits<Type>::epsilon());
            bool equalM22 = (diff.m22 < epsilon || std::abs(diff.m22 - epsilon) <= std::numeric_limits<Type>::epsilon());

            return (equalM00 && equalM01 && equalM02
                    && equalM10 && equalM11 && equalM12
                    && equalM20 && equalM21 && equalM22);
        }

        template<class Type>
        void Matrix3<Type>::rotationX(Type angle)
        {
            Type c = std::cos(angle);
            Type s = std::sin(angle);
            m00 = 1.0;
            m01 = 0.0;
            m02 = 0.0;

            m10 = 0.0;
            m11 = c;
            m12 = -s;

            m20 = 0.0;
            m21 = s;
            m22 = c;
        }

        template<class Type>
        void Matrix3<Type>::rotationY(Type angle)
        {
            Type c = std::cos(angle);
            Type s = std::sin(angle);
            m00 = c;
            m01 = 0.0;
            m02 = s;

            m10 = 0.0;
            m11 = 1.0;
            m12 = 0.0;

            m20 = -s;
            m21 = 0.0;
            m22 = c;
        }

        template<class Type>
        void Matrix3<Type>::rotationZ(Type angle)
        {
            Type c = std::cos(angle);
            Type s = std::sin(angle);
            m00 = c;
            m01 = -s;
            m02 = 0.0;

            m10 = s;
            m11 = c;
            m12 = 0.0;

            m20 = 0.0;
            m21 = 0.0;
            m22 = 1.0;
        }

        template<class Type>
        void Matrix3<Type>::transform(const Tuple3<Type> &vec, Tuple3<Type> &vecOut) const
        {
            vecOut.set(m00*vec.x + m01*vec.y + m02*vec.z,
                       m10*vec.x + m11*vec.y + m12*vec.z,
                       m20*vec.x + m21*vec.y + m22*vec.z);
        }

        template<class Type>
        void Matrix3<Type>::transform(Tuple3<Type> &vec) const
        {
            Type x = vec.x;
            Type y = vec.y;
            Type z = vec.z;
            vec.set(m00*x + m01*y + m02*z,
                    m10*x + m11*y + m12*z,
                    m20*x + m21*y + m22*z);
        }

        template<class Type>
        void Matrix3<Type>::negate()
        {
            m00 = -m00;
            m01 = -m01;
            m02 = -m02;

            m10 = -m10;
            m11 = -m11;
            m12 = -m12;

            m20 = -m20;
            m21 = -m21;
            m22 = -m22;
        }

        template<class Type>
        void Matrix3<Type>::negate(const Matrix3<Type> &m1)
        {
            set(m1);
            negate();
        }

        template<class Type>
        Matrix3<Type> Matrix3<Type>::operator-() const
        {
            Matrix3<Type> result(*this);
            result.negate();
            return result;
        }

        template<class Type>
        void Matrix3<Type>::absolute()
        {
            m00 = std::abs(m00);
            m01 = std::abs(m01);
            m02 = std::abs(m02);

            m10 = std::abs(m10);
            m11 = std::abs(m11);
            m12 = std::abs(m12);

            m20 = std::abs(m20);
            m21 = std::abs(m21);
            m22 = std::abs(m22);
        }

        template<class Type>
        std::string Matrix3<Type>::toString() const
        {
            std::string oss;
            oss = "[ [" + m00 + "," + m01 + ","
                  + m02 + "]" + '\n'
                  + "  [" + m10 + "," + m11 + ","
                  + m12 + "]" + '\n'
                  + "  [" + m20 + "," + m21 + ","
                  + m22 + "] ]" + std::endl;
            return oss;
        }

        template<class Type>
        std::ostream&operator<<(std::ostream &stream, const Matrix3<Type> &m1)
        {
            return stream << m1.toString();
        }

    }
}


#endif // ifndef mi_math_Matrix3_h
