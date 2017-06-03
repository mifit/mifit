#ifndef mi_math_Vector3_h
#define mi_math_Vector3_h

#include <math/Tuple3.h>
#include <math/mi_math.h>
#include <cmath>

namespace mi
{
    namespace math
    {

/**
 * A three element vector.
 */
        template<class Type>
        class Vector3 : public Tuple3<Type>
        {
        public:
            /**
             * Constructs and initializes a Vector3 to (0,0,0).
             */
            Vector3();

            /**
             * Constructs and initializes a Vector3 from x,y,z coordinates.
             *
             * @param x  the x coordinate.
             * @param y  the y coordinate.
             * @param z  the z coordinate.
             */
            Vector3(Type x, Type y, Type z);

            /*
             * Constructs and initializes a Vector3 with the value of a Tuple3.
             *
             * @param source  the tuple to be copied.
             */
            Vector3(const Tuple3<Type> &source);

            /**
             * Returns the angle in radians between this vector and the
             * vector v1. The return value is constrained to the
             * range [0,pi].
             *
             * @param v1  the other vector
             */
            Type angle(const Vector3 &v1) const;

            /**
             * Sets this vector to the vector cross product of vectors v1
             * and v2.
             *
             * @param v1  the first vector.
             * @param v2  the second vector.
             */
            void cross(const Vector3 &v1, const Vector3 &v2);

            /**
             * Returns the dot product of this vector and vector v1.
             *
             * @param v1  another vector.
             */
            Type dot(const Vector3 &v1) const;

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
            void normalize(const Vector3 &v1);

            /**
             * Sets the length of this vector. Equivalent to normalize() followed by scale(length).
             */
            void setLength(Type length);


        };


        template<class Type>
        Vector3<Type>::Vector3()
            : Tuple3<Type>(0.0, 0.0, 0.0)
        {
        }

        template<class Type>
        Vector3<Type>::Vector3(Type x, Type y, Type z)
            : Tuple3<Type>(x, y, z)
        {
        }

        template<class Type>
        Vector3<Type>::Vector3(const Tuple3<Type> &source)
            : Tuple3<Type>(source)
        {
        }

        template<class Type>
        Type Vector3<Type>::angle(const Vector3<Type> &v1) const
        {
            return std::acos(dot(v1)/(length() * v1.length()));
        }

        template<class Type>
        void Vector3<Type>::cross(const Vector3<Type> &v1, const Vector3<Type> &v2)
        {
            Vector3<Type>::set(v1.y*v2.z - v1.z*v2.y,
                v1.z*v2.x - v1.x*v2.z,
                v1.x*v2.y - v1.y*v2.x);
        }

        template<class Type>
        Type Vector3<Type>::dot(const Vector3<Type> &v1) const
        {
            return (Tuple3<Type>::x*v1.x + Tuple3<Type>::y*v1.y + Tuple3<Type>::z*v1.z);
        }

        template<class Type>
        Type Vector3<Type>::length() const
        {
            return std::sqrt(lengthSquared());
        }

        template<class Type>
        Type Vector3<Type>::lengthSquared() const
        {
            return (Tuple3<Type>::x*Tuple3<Type>::x + Tuple3<Type>::y*Tuple3<Type>::y + Tuple3<Type>::z*Tuple3<Type>::z);
        }

        template<class Type>
        void Vector3<Type>::normalize()
        {
            Type len = length();
            // Divide by zero could occur.

            // Should probably handle with an exception:
            // if (len < std::numeric_limits<Type>::epsilon()) {
            //  throw ArithmeticException("divide by zero");
            // }
            Tuple3<Type>::x /= len;
            Tuple3<Type>::y /= len;
            Tuple3<Type>::z /= len;
        }

        template<class Type>
        void Vector3<Type>::setLength(Type length)
        {
            normalize();
            scale(length);
        }

    }
}


#endif // ifndef mi_math_Vector3_h
