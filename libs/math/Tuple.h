#ifndef mi_math_Tuple_h
#define mi_math_Tuple_h

#include <iostream>
#include <math/ArrayIndexOutOfBoundsException.h>

namespace mi {
namespace math {

/**
 * A three element tuple.
 */
template<int dimension, class Type>
class Tuple {
protected:
  Type values[dimension];
public:
  /**
   * Constructs and initializes a Tuple to (0,0,0).
   */
  Tuple();

  /**
   * Constructs and initializes a Tuple from x,y,z coordinates.
   *
   * @param x  the x coordinate.
   * @param y  the y coordinate.
   * @param z  the z coordinate.
   */
  Tuple(const Type* valuesArray);

  /**
   * Constructs and initializes a Tuple with the value of another.
   *
   * @param source  the tuple to be copied.
   */
  Tuple(const Tuple& source);

  /**
   * Set the value of this tuple to the value of another.
   *
   * @param source  the tuple to be copied.
   */
  Tuple& operator=(const Tuple& source);

  /**
   * Sets the value of this tuple to the specified x,y coordinates.
   *
   * @param source  the tuple to be copied.
   */
  void set(const Type* valuesArray);

  /**
   * Returns a reference to the element at the specified index.
   *
   * @param index  the index of the element
   * @throws ArrayIndexOutOfBoundsException if (index < 0 || index > dimension)
   */
  double& operator[](int index);

  /**
   * Returns a copy of the element at the specified index.
   *
   * @param index  the index of the element
   * @throws ArrayIndexOutOfBoundsException if (index < 0 || index > dimension)
   */
  double operator[](int index) const;

  /**
   * Sets the value of this tuple to the sum of itself and tuple t1.
   *
   * @param t1  the other tuple.
   */
  Tuple& operator+=(const Tuple& t1);

  /**
   * Returns a new Tuple with the value of the sum of this tuple
   * and tuple t1.
   *
   * @param t1  the other tuple.
   * @return  a new Tuple.
   */
  Tuple operator+(const Tuple& t1) const;

  /**
   * Sets the value of this tuple to the sum of tuples t1 and t2.
   *
   * @param t1  the first tuple.
   * @param t2  the second tuple.
   */
  void add(const Tuple& t1, const Tuple& t2);

  /**
   * Sets the value of this tuple to the difference of itself and
   * tuple t1.
   *
   * @param t1 the other tuple.
   */
  Tuple& operator-=(const Tuple& t1);

  /**
   * Returns a new Tuple with the value of the difference
   * of this tuple and t1.
   *
   * @param t1  the first tuple.
   * @return  a new Tuple.
   */
  Tuple operator-(const Tuple& t1) const;

  /**
   * Sets the value of this tuple to the difference of tuples t1 and
   * t2.
   *
   * @param t1 the other tuple.
   */
  void subtract(const Tuple& t1, const Tuple& t2);

  /**
   * Returns a new Tuple with the value of the negation of this tuple.
   *
   * @return  a new Tuple.
   */
  Tuple operator-() const;

  /**
   * Negates the value of this tuple in place.
   */
  void negate();

  /**
   * Sets this tuple to the negaation of tuple t1.
   *
   * @param t1  the other tuple.
   */
  void negate(const Tuple& t1);

  /**
   * Sets each element of this tuple to its absolute value.
   */
  void absolute();

  /**
   * Sets each element of this tuple to the absolute value of tuple t1.
   *
   * @param t1  the other tuple.
   */
  void absolute(const Tuple& t1);

  /**
   * Returns true if each of the elements of this tuple are equal to each
   * corresponding element of tuple t1.
   *
   * @param t1  the other tuple.
   */
  bool operator==(const Tuple& t1) const;

  /**
   * Returns true if the difference of each of the elements of this
   * tuple and the corresponding element of tuple t1 are less than
   * or equal to the epsilon parameter.
   *
   * @param t1  the other tuple.
   * @param epsilon  the threshold value.
   */
  bool epsilonEquals(const Tuple& t1, double epsilon) const;

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
std::ostream& operator<<(std::ostream& stream, const Tuple& t1);

/**
 * Read a tuple to an input stream.
 *
 * @param stream  the input stream read from.
 * @param t1  the tuple to store the read value.
 * @return  the modified input stream.
 */
std::istream& operator>>(std::istream& stream, Tuple& t1);

}
}


#endif
