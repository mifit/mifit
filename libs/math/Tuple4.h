#ifndef mi_math_Tuple4_h
#define mi_math_Tuple4_h

#include <math/IndexOutOfBoundsException.h>
#include <cmath>
#include <limits>
#include <iostream>

#undef absolute
#undef abs
namespace mi
{
    namespace math
    {

/**
 * A four element tuple.
 */
        template<class Type> class Tuple4
        {
        public:
            Type x;
            Type y;
            Type z;
            Type w;

            /**
             * Constructs and initializes a Tuple4 to (0,0,0,0).
             */
            Tuple4();

            /**
             * Constructs and initializes a Tuple4 from x,y,z,w coordinates.
             *
             * @param x  the x coordinate.
             * @param y  the y coordinate.
             * @param z  the z coordinate.
             * @param w  the w coordinate.
             */
            Tuple4(Type x, Type y, Type z, Type w);


            /**
             * Sets the value of this quaternion to the specified x,y,z,w
             * coordinates.
             *
             * @param x  the x coordinate.
             * @param y  the y coordinate.
             * @param z  the z coordinate.
             * @param w the w coordinate.
             */
            void set(Type x, Type y, Type z, Type w);

            /**
             * Sets the value of this quaternion to the value of another.
             *
             * @param t1  the other tuple.
             */
            void set(const Tuple4<Type> &t1);

            /**
             * Returns a reference to the element at the specified index.
             *
             * @param index  the index of the element
             * @throws IndexOutOfBoundsException if (index < 0 || index > 3)
             */
            Type&operator[](int index);

            /**
             * Returns a copy of the element at the specified index.
             *
             * @param index  the index of the element
             * @throws IndexOutOfBoundsException if (index < 0 || index > 3)
             */
            Type operator[](int index) const;

            /**
             * Returns the first element of the tuple.
             */
            Type getX() const;

            /**
             * Returns the second element of the tuple.
             */
            Type getY() const;

            /**
             * Returns the thrid element of the tuple.
             */
            Type getZ() const;

            /**
             * Returns the fourth element of the tuple.
             */
            Type getW() const;

            /**
             * Sets the value of this tuple to the sum of itself and tuple t1.
             *
             * @param t1  the other tuple.
             */
            void add(const Tuple4<Type> &t1);

            /**
             * Sets the value of this tuple to the sum of itself and tuple t1.
             *
             * @param t1  the other tuple.
             */
            Tuple4<Type>&operator+=(const Tuple4<Type> &t1);

            /**
             * Sets the value of this tuple to the sum of tuples t1 and t2.
             *
             * @param t1  the first tuple.
             * @param t2  the second tuple.
             */
            void add(const Tuple4<Type> &t1, const Tuple4<Type> &t2);

            /**
             * Returns a new Tuple4 with the value of the sum of this tuple
             * and tuple t1.
             *
             * @param t1  the other tuple.
             * @return  a new Tuple4.
             */
            Tuple4<Type> operator+(const Tuple4<Type> &t1) const;

            /**
             * Sets the value of this tuple to the difference of itself and
             * tuple t1.
             *
             * @param t1 the other tuple.
             */
            void subtract(const Tuple4<Type> &t1);

            /**
             * Sets the value of this tuple to the difference of itself and
             * tuple t1.
             *
             * @param t1 the other tuple.
             */
            Tuple4<Type>&operator-=(const Tuple4<Type> &t1);

            /**
             * Sets the value of this tuple to the difference of tuples t1 and
             * t2.
             *
             * @param t1 the other tuple.
             */
            void subtract(const Tuple4<Type> &t1, const Tuple4<Type> &t2);

            /**
             * Returns a new Tuple4 with the value of the difference
             * of this tuple and t1.
             *
             * @param t1  the first tuple.
             * @return  a new Tuple4.
             */
            Tuple4<Type> operator-(const Tuple4<Type> &t1) const;

            /**
             * Sets the value of this tuple to the scalar multiplication
             * of itself and the scale factor.
             *
             * @param s  the scale factor.
             */
            void scale(Type s);

            /**
             * Sets the value of this tuple to the scalar multiplicationof
             * itself and the scale factor.
             *
             * @param s  the scale factor.
             */
            Tuple4<Type>&operator*=(Type s);

            /**
             * Sets the value of this tuple to the scalar multiplication
             * of another tuple and the scale factor.
             *
             * @param s  the scale factor.
             * @param t  the other tuple.
             */
            void scale(Type s, const Tuple4<Type> &t);

            /**
             * Negates the value of this tuple in place.
             */
            void negate();

            /**
             * Returns a new Tuple4 with the value of the negation of this tuple.
             *
             * @return  a new Tuple4.
             */
            Tuple4<Type> operator-() const;

            /**
             * Sets this tuple to the negaation of tuple t1.
             *
             * @param t1  the other tuple.
             */
            void negate(const Tuple4<Type> &t1);

            /**
             * Sets each element of this tuple to its absolute value.
             */
            void absolute();

            /**
             * Sets each element of this tuple to the absolute value of tuple t1.
             *
             * @param t1  the other tuple.
             */
            void absolute(const Tuple4<Type> &t1);

            /**
             * Returns true if each of the elements of this tuple are equal to each
             * corresponding element of tuple t1.
             *
             * @param t1  the other tuple.
             */
            bool operator==(const Tuple4<Type> &t1) const;

            /**
             * Returns true if the difference of each of the elements of this
             * tuple and the corresponding element of tuple t1 are less than
             * or equal to the epsilon parameter.
             *
             * @param t1  the other tuple.
             * @param epsilon  the threshold value.
             */
            bool epsilonEquals(const Tuple4<Type> &t1, Type epsilon) const;

        };

/**
 * Write a tuple to an output stream. This method is designed to
 * be used on a character stream, and may not work for other
 * streams (e.g., binary stream).
 *
 * @param stream  the output stream written to.
 * @param t1  the tuple to write.
 * @return the modified output stream.
 */
        template<class Type>
        std::ostream&operator<<(std::ostream &stream, const Tuple4<Type> &t1);

/**
 * Read a tuple to an input stream.
 *
 * @param stream  the input stream read from.
 * @param t1  the tuple to store the read value.
 * @return  the modified input stream.
 */
        template<class Type>
        std::istream&operator>>(std::istream &stream, Tuple4<Type> &t1);


        template<class Type>
        Tuple4<Type>::Tuple4()
            : x(0.0), y(0.0), z(0.0), w(0.0)
        {
        }

        template<class Type>
        Tuple4<Type>::Tuple4(Type x, Type y, Type z, Type w)
            : x(x), y(y), z(z), w(w)
        {
        }

        template<class Type>
        void Tuple4<Type>::set(Type x, Type y, Type z, Type w)
        {
            this->x = x;
            this->y = y;
            this->z = z;
            this->w = w;
        }

        template<class Type>
        void Tuple4<Type>::set(const Tuple4<Type> &t1)
        {
            x = t1.x;
            y = t1.y;
            z = t1.z;
            w = t1.w;
        }

        template<class Type>
        Type& Tuple4<Type>::operator[](int index)
        {
            switch (index)
            {
            case 0:
                return x;
            case 1:
                return y;
            case 2:
                return z;
            case 3:
                return w;
            default:
                throw IndexOutOfBoundsException();
            }
        }

        template<class Type>
        Type Tuple4<Type>::operator[](int index) const
        {
            switch (index)
            {
            case 0:
                return x;
            case 1:
                return y;
            case 2:
                return z;
            case 3:
                return w;
            default:
                throw IndexOutOfBoundsException();
            }
        }

        template<class Type>
        Type Tuple4<Type>::getX() const
        {
            return x;
        }

        template<class Type>
        Type Tuple4<Type>::getY() const
        {
            return y;
        }

        template<class Type>
        Type Tuple4<Type>::getZ() const
        {
            return z;
        }

        template<class Type>
        Type Tuple4<Type>::getW() const
        {
            return w;
        }

        template<class Type>
        void Tuple4<Type>::add(const Tuple4<Type> &t1)
        {
            x += t1.x;
            y += t1.y;
            z += t1.z;
            w += t1.w;
        }

        template<class Type>
        Tuple4<Type>&Tuple4<Type>::operator+=(const Tuple4<Type> &t1)
        {
            add(t1);
            return *this;
        }

        template<class Type>
        void Tuple4<Type>::add(const Tuple4<Type> &t1, const Tuple4<Type> &t2)
        {
            set(t1);
            add(t2);
        }

        template<class Type>
        Tuple4<Type> Tuple4<Type>::operator+(const Tuple4<Type> &t1) const
        {
            Tuple4<Type> result(*this);
            result.add(t1);
            return result;
        }

        template<class Type>
        void Tuple4<Type>::subtract(const Tuple4<Type> &t1)
        {
            x -= t1.x;
            y -= t1.y;
            z -= t1.z;
            w -= t1.w;
        }

        template<class Type>
        Tuple4<Type>&Tuple4<Type>::operator-=(const Tuple4<Type> &t1)
        {
            subtract(t1);
            return *this;
        }

        template<class Type>
        void Tuple4<Type>::subtract(const Tuple4<Type> &t1, const Tuple4<Type> &t2)
        {
            set(t1);
            subtract(t2);
        }

        template<class Type>
        Tuple4<Type> Tuple4<Type>::operator-(const Tuple4<Type> &t1) const
        {
            Tuple4<Type> result(*this);
            result.subtract(t1);
            return result;
        }

        template<class Type>
        void Tuple4<Type>::scale(Type s)
        {
            x *= s;
            y *= s;
            z *= s;
            w *= s;
        }

        template<class Type>
        Tuple4<Type>&Tuple4<Type>::operator*=(Type s)
        {
            scale(s);
            return *this;
        }

        template<class Type>
        void Tuple4<Type>::scale(Type s, const Tuple4<Type> &t)
        {
            set(t);
            scale(s);
        }

        template<class Type>
        void Tuple4<Type>::negate()
        {
            x = -x;
            y = -y;
            z = -z;
            w = -w;
        }

        template<class Type>
        Tuple4<Type> Tuple4<Type>::operator-() const
        {
            Tuple4<Type> result(*this);
            result.negate();
            return result;
        }

        template<class Type>
        void Tuple4<Type>::negate(const Tuple4<Type> &t1)
        {
            set(t1);
            negate();
        }

        template<class Type>
        void Tuple4<Type>::absolute()
        {
            x = std::abs(x);
            y = std::abs(y);
            z = std::abs(z);
            w = std::abs(w);
        }

        template<class Type>
        void Tuple4<Type>::absolute(const Tuple4<Type> &t1)
        {
            set(t1);
            absolute();
        }

        template<class Type>
        bool Tuple4<Type>::operator==(const Tuple4<Type> &t1) const
        {
            return (equals(x, t1.x) && equals(y, t1.y) && equals(z, t1.z) && equals(w, t1.w));
        }

        template<class Type>
        bool Tuple4<Type>::epsilonEquals(const Tuple4<Type> &t1, Type epsilon) const
        {
            Tuple4<Type> diff(*this);
            diff.subtract(t1);
            diff.absolute();
            bool xEqual = (diff.x < epsilon || std::abs(diff.x - epsilon) <= std::numeric_limits<Type>::epsilon());
            bool yEqual = (diff.y < epsilon || std::abs(diff.y - epsilon) <= std::numeric_limits<Type>::epsilon());
            bool zEqual = (diff.z < epsilon || std::abs(diff.z - epsilon) <= std::numeric_limits<Type>::epsilon());
            bool wEqual = (diff.w < epsilon || std::abs(diff.w - epsilon) <= std::numeric_limits<Type>::epsilon());
            return (xEqual && yEqual && zEqual && wEqual);
        }

        template<class Type>
        std::ostream&operator<<(std::ostream &stream, const Tuple4<Type> &t1)
        {
            stream << t1.getX();
            stream << ' ';
            stream << t1.getY();
            stream << ' ';
            stream << t1.getZ();
            stream << ' ';
            stream << t1.getW();
            return stream;
        }

        template<class Type>
        std::istream&operator>>(std::istream &stream, Tuple4<Type> &t1)
        {
            Type x;
            Type y;
            Type z;
            Type w;
            stream >> x;
            stream >> y;
            stream >> z;
            stream >> w;
            t1.set(x, y, z, w);
            return stream;
        }

    }
}


#endif // ifndef mi_math_Tuple4_h
