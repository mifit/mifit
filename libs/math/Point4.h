#ifndef mi_math_Point4_h
#define mi_math_Point4_h

#include <math/Tuple4.h>
#include <cmath>

namespace mi {
namespace math {

/**
 * A four element point representing Type x,y,z,w coordinates.
 */
template<class Type>
class Point4 : public Tuple4<Type> {
public:
  /**
   * Constructs and initializes a Point4 to (0,0,0,0).
   */
  Point4();

  /**
   * Constructs and initializes a Point4 from x,y,z,w coordinates.
   *
   * @param x  the x coordinate.
   * @param y  the y coordinate.
   * @param z  the z coordinate.
   * @param w  the w coordinate.
   */
  Point4(Type x, Type y, Type z, Type w);

  /**
   * Constructs and initializes a Point4 from a Tuple4.
   */
  Point4(const Tuple4<Type>& t1);

  /**
   * Returns the distance between this point and point p1.
   *
   * @param p1  the other point.
   */
  Type distance(const Point4& p1) const;

  /**
   * Returns the square of the distance between this point and
   * point p1.
   *
   * @param p1  the other point.
   */
  Type distanceSquared(const Point4& p1) const;

  /**
   * Multiplies each of the x,y,z components of the Point4
   * parameter by 1/w, places the projected values into this
   * point, and places a 1 as the w parameter of this point.
   *
   * @param p1 the source Point4.
   */
  void project(const Point4& p1);

};


template<class Type>
Point4<Type>::Point4()
  : Tuple4<Type>(0.0, 0.0, 0.0, 0.0) {
}

template<class Type>
Point4<Type>::Point4(Type x, Type y, Type z, Type w)
  : Tuple4<Type>(x, y, z, w) {
}

template<class Type>
Point4<Type>::Point4(const Tuple4<Type>& t1)
  : Tuple4<Type>(t1) {
}

template<class Type>
Type Point4<Type>::distance(const Point4<Type>& p1) const {
  return std::sqrt(distanceSquared(p1));
}

template<class Type>
Type Point4<Type>::distanceSquared(const Point4<Type>& p1) const {
  Type dx = Tuple4<Type>::x - p1.x;
  Type dy = Tuple4<Type>::y - p1.y;
  Type dz = Tuple4<Type>::z - p1.z;
  Type dw = Tuple4<Type>::w - p1.w;
  return (dx*dx + dy*dy + dz*dz + dw*dw);
}

template<class Type>
void Point4<Type>::project(const Point4<Type>& p1) {
  set(p1.getX()/p1.getW(), p1.getY()/p1.getW(), p1.getZ()/p1.getW(), 1.0);
}

}
}


#endif
