#ifndef mi_math_Quaternion_h
#define mi_math_Quaternion_h

#include <math/Tuple4.h>
#include <math/Random.h>

namespace mi
{
    namespace math
    {
        template<class Type> class Quaternion;
    }
}
#include <math/Matrix3.h>
#include <math/Matrix4.h>

namespace mi
{
    namespace math
    {

/**
 * A quaternion.
 */
        template<class Type>
        class Quaternion : public Tuple4<Type>
        {
        public:
            /**
             * Constructs and initializes a Quaternion to (0,0,0,0).
             */
            Quaternion();

            /**
             * Constructs and initializes a Quaternion from x,y,z,w coordinates.
             *
             * @param x  the x coordinate.
             * @param y  the y coordinate.
             * @param z  the z coordinate.
             * @param w  the w coordinate.
             */
            Quaternion(Type x, Type y, Type z, Type w);

            /*
             * Constructs and initializes a Quaternion with the value of a
             * Tuple4.
             *
             * @param source the tuple to be copied.
             */
            Quaternion(const Tuple4<Type> &source);

            /**
             * Creates a quaternion to the specified axis and angle.
             *
             * @param axis  the axis.
             * @param angle  the angle to rotate about the given axis (in radians).
             */
            Quaternion(const Vector3<Type> &axis, Type angle);

            /**
             * Creates a rotation which will rotate the source vector onto the
             * target vector.
             *
             * @param sourceVector the vector from which to rotate
             * @param targetVector the vector to which to rotate
             */
            Quaternion(const Vector3<Type> &sourceVector, const Vector3<Type> &targetVector);


            /**
             * Sets the value of this tuple to the specified x,y,z,w coordinates.
             *
             * @param x  the x coordinate.
             * @param y  the y coordinate.
             * @param z  the z coordinate.
             * @param w  the w coordinate.
             */
            void set(Type x, Type y, Type z, Type w);


            /**
             * Sets to a rotation which will rotate the source vector onto the
             * target vector.
             *
             * @param sourceVector the vector from which to rotate
             * @param targetVector the vector to which to rotate
             */
            void set(const Vector3<Type> &sourceVector, const Vector3<Type> &targetVector);

            /**
             * Sets the value of this quaternion to the specified axis and angle.
             *
             * @param axis  the axis.
             * @param angle  the angle to rotate about the given axis (in radians).
             */
            void set(const Vector3<Type> &axis, Type angle);

            /**
             * Sets the value of this quaternion to the specified axis and angle.
             *
             * @param x  the x coordinate of the axis.
             * @param y  the y coordinate of the axis.
             * @param z  the z coordinate of the axis.
             * @param angle  the angle to rotate about the given axis (in radians).
             */
            void setFromAxisAngle(Type x, Type y, Type z, Type angle);

            /**
             * Sets the value of this quaternion by the specified Euler angles.
             *
             * @param xAngle  rotation about the x axis.
             * @param yAngle  rotation about the y axis.
             * @param zAngle  rotation about the z axis.
             */
            void setFromEulerAngles(Type xAngle, Type yAngle, Type zAngle);

            /**
             * Sets the value of this tuple to the value of another.
             *
             * @param t1  the other tuple.
             */
            void set(const Tuple4<Type> &t1);

            /**
             * Sets the value of this quaternion to the rotational
             * component of the passed matrix.
             */
            void set(const Matrix3<Type> &m1);

            /**
             * Sets the value of this quaternion to the rotational
             * component of the passed matrix.
             */
            void set(const Matrix4<Type> &m1);

            /**
             * Returns the norm of this quaternion. (i.e., the sum of the
             * squares of the components)
             */
            Type norm() const;

            /**
             * Negate the value of of each of this quaternion's x,y,z
             * coordinates in place.
             */
            void conjugate();

            /**
             * Sets the value of this quaternion to the conjugate of
             * quaternion q1.
             */
            void conjugate(const Quaternion &q1);

            /**
             * Sets the value of this quaternion to the quaternion inverse
             * of itself.
             */
            void inverse();

            /**
             * Sets the value of this quaternion to quaternion inverse of
             * quaternion q1.
             */
            void inverse(const Quaternion &q1);

            /**
             * Sets the value of this quaternion to the quaternion product
             * of itself and q1 (this = this * q1).
             */
            void multiply(const Quaternion &q1);

            /**
             * Sets the value of this quaternion to the quaternion product
             * of itself and q1 (this = this * q1).
             */
            Quaternion&operator*=(const Quaternion &q1);

            /**
             * Sets the value of this quaternion to the quaternion product
             * of quaternions q1 and q2 (this = q1 * q2).
             */
            void multiply(const Quaternion &q1, const Quaternion &q2);

            /**
             * Returns a new Quaternion with the value of the quaternion
             * product of quaternions q1 and q2 (q = q1 * q2).
             */
            Quaternion operator*(const Quaternion &q1) const;

            /**
             * Rotate a vector by this quaternion (v2 = this * v1).
             */
            Vector3<Type> rotate(const Vector3<Type> &v1) const;

            /**
             * Rotate a vector by this quaternion (v2 = this * v1).
             */
            Vector3<Type> operator*(const Vector3<Type> &v1) const;

            /**
             * Normalizes the value of this quaternion in place.
             */
            void normalize();

            /**
             * Sets the value of this quaternion to the normalized value
             * of quaternion q1.
             */
            void normalize(const Quaternion &q1);

            /**
             * Performs a great circle interpolation between this
             * quaternion and the quaternion parameter and places the
             * result into this quaternion.
             *
             * @param q1 the other quaternion
             * @param omega the amount of interpolation [0.0, 1.0]
             */
            void interpolate(const Quaternion &q1, Type omega);

            /**
             * Performs a great circle interpolation between quaternion q1
             * and quaternion q2 and places the result into this
             * quaternion.
             *
             * @param q1 the first quaternion
             * @param q2 the second quaternion
             * @param omega the amount of interpolation [0.0, 1.0]
             */
            void interpolate(const Quaternion &q1, const Quaternion &q2, Type omega);

            /**
             *  Generates a random, uniformly distributed, unit quaternion.
             */
            static Quaternion random(Random &rng);

            /**
             *  Generates a unit quaternion with a uniformly random direction,
             *  and an amount of rotation of the specified angle.
             *
             *  @param angle  the amount of rotation about the random axis.
             */
            static Quaternion random(Random &rng, double angle);

        private:

            /**
             * Construct a unit quaternion from rotation matrix. Assumes
             * matrix is used to multiply column vector on the left: vnew
             * = mat vold. Works correctly for right-handed coordinate
             * system and right-handed rotations.  Translation and
             * perspective components ignored.
             */
            void setFromMatrix(Type m00, Type m01, Type m02, Type m03,
                               Type m10, Type m11, Type m12, Type m13,
                               Type m20, Type m21, Type m22, Type m23,
                               Type m30, Type m31, Type m32, Type m33);
        };


        template<class Type>
        Quaternion<Type>::Quaternion()
            : Tuple4<Type>(0.0, 0.0, 0.0, 0.0)
        {
        }

        template<class Type>
        Quaternion<Type>::Quaternion(Type x, Type y, Type z, Type w)
            : Tuple4<Type>(x, y, z, w)
        {
        }

        template<class Type>
        Quaternion<Type>::Quaternion(const Tuple4<Type> &source)
            : Tuple4<Type>(source)
        {
        }

        template<class Type>
        Quaternion<Type>::Quaternion(const Vector3<Type> &axis, Type angle)
        {
            Quaternion<Type>::setFromAxisAngle(axis.x, axis.y, axis.z, angle);
        }

        template<class Type>
        Quaternion<Type>::Quaternion(const Vector3<Type> &sourceVector, const Vector3<Type> &targetVector)
        {
            Quaternion<Type>::set(sourceVector, targetVector);
        }

        template<class Type>
        void Quaternion<Type>::set(Type x, Type y, Type z, Type w)
        {
            Tuple4<Type>::set(x, y, z, w);
        }

        template<class Type>
        void Quaternion<Type>::set(const Vector3<Type> &axis, Type angle)
        {
            Quaternion<Type>::setFromAxisAngle(axis.x, axis.y, axis.z, angle);
        }

        template<class Type>
        void Quaternion<Type>::setFromAxisAngle(Type x, Type y, Type z, Type angle)
        {
            const Type epsilon = (Type) 0.0000001;

            Type length = std::sqrt(x*x + y*y + z*z);
            if (length < epsilon)
            {
                // ~zero length axis, so reset rotation to zero.
                Tuple4<Type>::set(0.0, 0.0, 0.0, 1.0);
                return;
            }

            Type inverseNorm  = (Type)(1.0/length);
            Type cosHalfAngle = (Type) std::cos(0.5*angle);
            Type sinHalfAngle = (Type) std::sin(0.5*angle);

            Tuple4<Type>::set(x * sinHalfAngle * inverseNorm,
                              y * sinHalfAngle * inverseNorm,
                              z * sinHalfAngle * inverseNorm,
                              cosHalfAngle);
        }

        template<class Type>
        void Quaternion<Type>::setFromEulerAngles(Type xAngle, Type yAngle, Type zAngle)
        {
            Quaternion<Type> qx;
            Quaternion<Type> qy;
            Quaternion<Type> qz;

            qx.setFromAxisAngle(1.0, 0.0, 0.0, xAngle);
            qy.setFromAxisAngle(0.0, 1.0, 0.0, yAngle);
            qz.setFromAxisAngle(0.0, 0.0, 1.0, zAngle);

            multiply(qx, qy);
            multiply(qz);
        }

        template<class Type>
        void Quaternion<Type>::set(const Tuple4<Type> &t1)
        {
            Tuple4<Type>::set(t1);
        }

        template<class Type>
        Type Quaternion<Type>::norm() const
        {
            return Tuple4<Type>::x*Tuple4<Type>::x + Tuple4<Type>::y*Tuple4<Type>::y + Tuple4<Type>::z*Tuple4<Type>::z + Tuple4<Type>::w*Tuple4<Type>::w;
        }

        template<class Type>
        void Quaternion<Type>::conjugate()
        {
            Tuple4<Type>::x = -Tuple4<Type>::x;
            Tuple4<Type>::y = -Tuple4<Type>::y;
            Tuple4<Type>::z = -Tuple4<Type>::z;
        }

        template<class Type>
        void Quaternion<Type>::conjugate(const Quaternion<Type> &q1)
        {
            set(q1);
            conjugate();
        }

        template<class Type>
        void Quaternion<Type>::inverse()
        {
            Type n = norm();
            // Divide by zero may occur.
            Tuple4<Type>::x = -Tuple4<Type>::x/n;
            Tuple4<Type>::y = -Tuple4<Type>::y/n;
            Tuple4<Type>::z = -Tuple4<Type>::z/n;
            Tuple4<Type>::w /= n;
        }

        template<class Type>
        void Quaternion<Type>::inverse(const Quaternion<Type> &q1)
        {
            set(q1);
            inverse();
        }

        template<class Type>
        void Quaternion<Type>::normalize()
        {
            Type n = std::sqrt(norm());
            // Divide by zero may occur.
            Tuple4<Type>::x /= n;
            Tuple4<Type>::y /= n;
            Tuple4<Type>::z /= n;
            Tuple4<Type>::w /= n;
        }

        template<class Type>
        void Quaternion<Type>::normalize(const Quaternion<Type> &q1)
        {
            set(q1);
            normalize();
        }

        template<class Type>
        void Quaternion<Type>::multiply(const Quaternion<Type> &q1)
        {
            set(Tuple4<Type>::x*q1.w + Tuple4<Type>::w*q1.x + Tuple4<Type>::y*q1.z - Tuple4<Type>::z*q1.y,
                Tuple4<Type>::y*q1.w + Tuple4<Type>::w*q1.y + Tuple4<Type>::z*q1.x - Tuple4<Type>::x*q1.z,
                Tuple4<Type>::z*q1.w + Tuple4<Type>::w*q1.z + Tuple4<Type>::x*q1.y - Tuple4<Type>::y*q1.x,
                Tuple4<Type>::w*q1.w - Tuple4<Type>::x*q1.x - Tuple4<Type>::y*q1.y - Tuple4<Type>::z*q1.z);
        }

        template<class Type>
        Quaternion<Type>&Quaternion<Type>::operator*=(const Quaternion<Type> &q1)
        {
            multiply(q1);
            return *this;
        }

        template<class Type>
        void Quaternion<Type>::multiply(const Quaternion<Type> &q1, const Quaternion<Type> &q2)
        {
            set(q1);
            multiply(q2);
        }

        template<class Type>
        Quaternion<Type> Quaternion<Type>::operator*(const Quaternion<Type> &q1) const
        {
            Quaternion<Type> result(*this);
            result.multiply(q1);
            return result;
        }

        template<class Type>
        Vector3<Type> Quaternion<Type>::rotate(const Vector3<Type> &v1) const
        {
            return Quaternion<Type>::operator*(v1);
        }

        template<class Type>
        Vector3<Type> Quaternion<Type>::operator*(const Vector3<Type> &v1) const
        {

            Quaternion<Type> copy(*this);
            Quaternion<Type> qv(v1.x, v1.y, v1.z, 0.0);
            copy.multiply(qv);
            Quaternion<Type> conj(*this);
            conj.conjugate();
            copy.multiply(conj);

            return Vector3<Type>(copy.x, copy.y, copy.z);
        }

        template<class Type>
        void Quaternion<Type>::set(const Matrix3<Type> &m1)
        {
            setFromMatrix(m1.get00(), m1.get01(), m1.get02(), 0.0,
                          m1.get10(), m1.get11(), m1.get12(), 0.0,
                          m1.get20(), m1.get21(), m1.get22(), 0.0,
                          0.0, 0.0, 0.0, 1.0);
        }

        template<class Type>
        void Quaternion<Type>::set(const Matrix4<Type> &m1)
        {
            setFromMatrix(m1.get00(), m1.get01(), m1.get02(), m1.get03(),
                          m1.get10(), m1.get11(), m1.get12(), m1.get13(),
                          m1.get20(), m1.get21(), m1.get22(), m1.get23(),
                          m1.get30(), m1.get31(), m1.get32(), m1.get33());
        }

        template<class Type>
        void Quaternion<Type>::setFromMatrix(Type m00, Type m01, Type m02, Type /* m03 */,
                                             Type m10, Type m11, Type m12, Type /* m13 */,
                                             Type m20, Type m21, Type m22, Type /* m23 */,
                                             Type /* m30 */, Type /* m31 */, Type /* m32 */, Type m33)
        {
            // From Ken Shoemake, ftp://ftp.cis.upenn.edu/pub/graphics/shoemake

            /* This algorithm avoids near-zero divides by looking for a large
             * component first w, then x, y, or z. When the trace is greater
             * than zero, |w| is greater than 1/2, which is as small as a
             * largest component can be.  Otherwise, the largest diagonal
             * entry corresponds to the largest of |x|, |y|, or |z|, one of
             * which must be larger than |w|, and at least 1/2.
             */
            Type s;
            Type trace = m00 + m11 + m22;
            if (trace >= 0.0)
            {
                s = (Type) std::sqrt(trace + m33);
                Tuple4<Type>::w = s*(Type)0.5;
                s = (Type)(0.5/s);
                Tuple4<Type>::x = (m21 - m12)*s;
                Tuple4<Type>::y = (m02 - m20)*s;
                Tuple4<Type>::z = (m10 - m01)*s;
            }
            else
            {
                int h = 0;
                if (m11 > m00)
                {
                    if (m22 > m11)
                    {
                        h = 2;
                    }
                    else
                    {
                        h = 1;
                    }
                }
                else
                {
                    if (m22 > m00)
                    {
                        h = 2;
                    }
                }
                switch (h)
                {
                case 0:
                    s = std::sqrt(m00 - (m11 + m22) + m33);
                    Tuple4<Type>::x = (Type)(s*0.5);
                    s = (Type)(0.5/s);
                    Tuple4<Type>::y = (m01 + m10)*s;
                    Tuple4<Type>::z = (m20 + m02)*s;
                    Tuple4<Type>::w = (m21 - m12)*s;
                    break;
                case 1:
                    s = std::sqrt(m11 - (m22 + m00) + m33);
                    Tuple4<Type>::y = (Type)(s*0.5);
                    s = (Type)(0.5/s);
                    Tuple4<Type>::z = (m12 + m21)*s;
                    Tuple4<Type>::x = (m01 + m10)*s;
                    Tuple4<Type>::w = (m02 - m20)*s;
                    break;
                case 2:
                    s = std::sqrt(m22 - (m00 + m11) + m33);
                    Tuple4<Type>::z = (Type)(s*0.5);
                    s = (Type)(0.5/s);
                    Tuple4<Type>::x = (m20 + m02)*s;
                    Tuple4<Type>::y = (m12 + m21)*s;
                    Tuple4<Type>::w = (m10 - m01)*s;
                    break;
                }
            }
            if (m33 != 1.0)
            {
                s = (Type)(1.0/std::sqrt(m33));
                Tuple4<Type>::x *= s;
                Tuple4<Type>::y *= s;
                Tuple4<Type>::z *= s;
                Tuple4<Type>::w *= s;
            }
        }

        template<class Type>
        Quaternion<Type> Quaternion<Type>::random(Random &rng)
        {
            /*
             *  We define the four elements of the quaternion in terms of
             *  four "hyperspherical" coordinates r, alpha, beta, and
             *  gamma.
             *
             *  q0 = r * cos(alpha)
             *  q1 = r * sin(alpha) * cos(beta)
             *  q2 = r * sin(alpha) * sin(beta) * cos(gamma)
             *  q3 = r * sin(alpha) * sin(beta) * sin(gamma)
             *
             *  where 0 < alpha < pi/2
             *        0 < beta < pi
             *        0 < gamma < 2 pi
             *
             *  Note that these choosing r = 1 ensures that q0^2 + q1^2 +
             *  q2^2 + q3^2 = 1
             *
             *  It is straightforward to show that the Jacobian
             *
             *  dq0 dq1 dq2 dq3 = r^3 sin^2(alpha) sin(beta) dr dalpha
             *  dbeta dgamma
             *
             *  We want to generate a uniform distribution in
             *  (q0,q1,q2,q3) subject to the condition that q0^2 + q1^2 +
             *  q2^2 + q3^2 = 1
             *
             *  Choose the new variable (x,y,z) and generate uniform
             *  distributions for each in the ranges:
             *              0  <  x  < 2 pi
             *              -1 <  y  < 1
             *              0  <  z  < pi / 2
             *
             *  Now let
             *      gamma = x
             *      cos(beta) = -y
             *      alpha - sin(2 alpha)/2 = z
             *
             *  Then use the equations for (q0,q1,q2,q3) at the top to
             *  convert to the quaternions
             *
             *
             */

            const double DELTA_ALPHA = 1.0e-6;

            double gamma = rng.nextDouble() * 2.0 * M_PI;

            double beta = std::acos(2.0 * rng.nextDouble() - 1.0);

            double z = rng.nextDouble() * M_PI * 0.5;

            /*
             * Use Newton's method to solve for alpha.
             */

            double alpha = M_PI/4.0;
            double oldAlpha;

            do
            {
                oldAlpha = alpha;

                alpha = oldAlpha - (oldAlpha - std::sin(2.0 * oldAlpha) / 2.0 - z)
                        / (1.0 - std::cos(2.0 * alpha));

            } while (std::abs(alpha - oldAlpha) > DELTA_ALPHA);

            double sinAlpha = std::sin(alpha);
            double sinBeta  = std::sin(beta);

            return Quaternion<Type>(sinAlpha * std::cos(beta),
                                    sinAlpha * sinBeta * std::cos(gamma),
                                    sinAlpha * sinBeta * std::sin(gamma),
                                    std::cos(alpha));
        }

        template<class Type>
        Quaternion<Type> Quaternion<Type>::random(Random &rng, double angle)
        {
            double randomCosTheta = rng.nextDouble() * 2.0 - 1.0;
            double randomPhi      = rng.nextDouble() * 2.0 * M_PI;

            double nx = std::sqrt(1 - randomCosTheta * randomCosTheta) * std::cos(randomPhi);
            double ny = std::sqrt(1 - randomCosTheta * randomCosTheta) * std::sin(randomPhi);
            double nz = randomCosTheta;

            double sinThetaOverTwo = std::sin(angle / 2.0);

            return Quaternion<Type>(nx * sinThetaOverTwo,
                                    ny * sinThetaOverTwo,
                                    nz * sinThetaOverTwo,
                                    std::cos(angle / 2.0));
        }

        template<class Type>
        void Quaternion<Type>::interpolate(const Quaternion &q1, Type omega)
        {
            normalize();
            Type n1 = std::sqrt(q1.norm());

            // zero-div may occur.
            Type x1 = q1.x/n1;
            Type y1 = q1.y/n1;
            Type z1 = q1.z/n1;
            Type w1 = q1.w/n1;

            // t is cosine (dot product)
            Type t = Tuple4<Type>::x*x1 + Tuple4<Type>::y*y1 + Tuple4<Type>::z*z1 + Tuple4<Type>::w*w1;

            // same quaternion (avoid domain error)
            if (1.0 <= std::abs(t))
            {
                return;
            }

            // t is now theta
            t = std::acos(t);

            Type sin_t = std::sin(t);

            // same quaternion (avoid zero-div)
            if (sin_t == 0.0)
            {
                return;
            }

            Type s = std::sin((1.0-omega)*t)/sin_t;
            t = std::sin(omega*t)/sin_t;

            // set values
            Tuple4<Type>::x = s*Tuple4<Type>::x + t*x1;
            Tuple4<Type>::y = s*Tuple4<Type>::y + t*y1;
            Tuple4<Type>::z = s*Tuple4<Type>::z + t*z1;
            Tuple4<Type>::w = s*Tuple4<Type>::w + t*w1;
        }

        template<class Type>
        void Quaternion<Type>::interpolate(const Quaternion &q1, const Quaternion &q2, Type omega)
        {
            set(q1);
            interpolate(q2, omega);
        }

        template<class Type>
        void Quaternion<Type>::set(const Vector3<Type> &sourceVector, const Vector3<Type> &targetVector)
        {

            const Type epsilon = 1.0e-7;
            double dotProdPlus1 = 1.0 + sourceVector.dot(targetVector);
            // Check for degenerate case of full u-turn. Use epsilon for detection
            if (dotProdPlus1 < epsilon)
            {
                // Get an orthogonal vector of the given vector
                // in a plane with maximum vector coordinates.
                // Then use it as quaternion axis with pi angle
                if (std::abs(sourceVector.x) < 0.6)
                {
                    const double norm = std::sqrt(1.0 - sourceVector.x * sourceVector.x);
                    set(0.0, sourceVector.z/norm, -sourceVector.y/norm, 0.0);
                }
                else if (fabs(sourceVector.y) < 0.6)
                {
                    const double norm = std::sqrt(1.0 - sourceVector.y * sourceVector.y);
                    set(-sourceVector.z/norm, 0.0, sourceVector.x/norm, 0.0);
                }
                else
                {
                    const double norm = std::sqrt(1.0 - sourceVector.z * sourceVector.z);
                    set(sourceVector.y/norm, -sourceVector.x/norm, 0.0, 0.0);
                }
            }
            else
            {
                // Find the shortest angle quaternion that transforms normalized vectors
                // into one other.
                const double s = std::sqrt(0.5 * dotProdPlus1);
                Vector3<Type> v;
                v.cross(sourceVector, targetVector);
                v.scale(1.0 / (2.0 * s));
                set(v.x, v.y, v.z, s);
            }
            normalize();
        }

    }
}


#endif // ifndef mi_math_Quaternion_h
