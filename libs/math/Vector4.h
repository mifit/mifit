#ifndef mi_math_Vector4_h
#define mi_math_Vector4_h

#include <math/Vector4.h>
#include <math/Tuple4.h>
#include <cmath>

namespace mi
{
    namespace math
    {

/**
 * A four element vector.
 */
        template<class Type>
        class Vector4 : public Tuple4<Type>
        {
        public:
            /**
             * Constructs and initializes a Vector4 to (0,0,0,0).
             */
            Vector4();

            /**
             * Constructs and initializes a Vector4 from x,y,z coordinates.
             *
             * @param x  the x coordinate.
             * @param y  the y coordinate.
             * @param z  the z coordinate.
             * @param w  the w coordinate.
             */
            Vector4(Type x, Type y, Type z, Type w);

            /*
             * Constructs and initializes a Vector4 with the value of a Tuple4.
             *
             * @param source  the tuple to be copied.
             */
            Vector4(const Tuple4<Type> &source);

            /**
             * Returns the angle in radians between this vector and the
             * vector v1. The return value is constrained to the
             * range [0,pi].
             *
             * @param v1  the other vector
             */
            Type angle(const Vector4 &v1);

            /**
             * Returns the dot product of this vector and vector v1.
             *
             * @param v1  another vector.
             */
            Type dot(const Vector4 &v1) const;

            /**
             * Returns the length of this vector.
             */
            Type length() const;

            /**
             * Returns the squared length of this vector.
             */
            Type lengthSquared() const;

            /**
             * Normalizes this vector.
             */
            void normalize();

            /**
             * Sets the value of this vector to the normalization of
             * vector v1.
             *
             * @param v1  another vector.
             */
            void normalize(const Vector4 &v1);

            /**
             * Sets the length of this vector. Equivalent to normalize() followed by scale(length).
             */
            void setLength(Type length);

        };


        template<class Type>
        Vector4<Type>::Vector4()
            : Tuple4<Type>(0.0, 0.0, 0.0, 0.0)
        {
        }

        template<class Type>
        Vector4<Type>::Vector4(Type x, Type y, Type z, Type w)
            : Tuple4<Type>(x, y, z, w)
        {
        }

        template<class Type>
        Vector4<Type>::Vector4(const Tuple4<Type> &source)
            : Tuple4<Type>(source)
        {
        }

        template<class Type>
        Type Vector4<Type>::angle(const Vector4<Type> &v1)
        {
            return std::acos(dot(v1)/(length() * v1.length()));
        }

        template<class Type>
        Type Vector4<Type>::dot(const Vector4<Type> &v1) const
        {
            return (Tuple4<Type>::x*v1.x + Tuple4<Type>::y*v1.y + Tuple4<Type>::z*v1.z + Tuple4<Type>::w*v1.w);
        }

        template<class Type>
        Type Vector4<Type>::length() const
        {
            return std::sqrt(lengthSquared());
        }

        template<class Type>
        Type Vector4<Type>::lengthSquared() const
        {
            return (Tuple4<Type>::x*Tuple4<Type>::x + Tuple4<Type>::y*Tuple4<Type>::y + Tuple4<Type>::z*Tuple4<Type>::z + Tuple4<Type>::w*Tuple4<Type>::w);
        }

        template<class Type>
        void Vector4<Type>::normalize()
        {
            Type len = length();
            // Divide by zero could occur.

            // Should probably handle with an exception:
            // if (len < DBL_EPSILON) {
            //  throw ArithmeticException("divide by zero");
            // }
            Tuple4<Type>::x /= len;
            Tuple4<Type>::y /= len;
            Tuple4<Type>::z /= len;
            Tuple4<Type>::w /= len;
        }

        template<class Type>
        void Vector4<Type>::setLength(Type length)
        {
            normalize();
            scale(length);
        }

    }
}


#endif // ifndef mi_math_Vector4_h
