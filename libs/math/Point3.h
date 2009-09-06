#ifndef mi_math_Point3_h
#define mi_math_Point3_h

#include <math/Tuple3.h>
#include <math/Point4.h>

namespace mi {
namespace math {

/**
 * A three element point representing Type x,y,z coordinates.
 */
template<class Type>
class Point3 : public Tuple3<Type> {
public:
  /**
   * Constructs and initializes a Point3 to (0,0,0).
   */
  Point3();

  /**
   * Constructs and initializes a Point3 from x,y,z coordinates.
   *
   * @param x  the x coordinate.
   * @param y  the y coordinate.
   * @param z  the z coordinate.
   */
  Point3(Type x, Type y, Type z);

  /**
   * Constructs and initializes a Point3 from a Tuple3.
   */
  Point3(const Tuple3<Type>& t1);

  /**
   * Returns the distance between this point and point p1.
   *
   * @param p1  the other point.
   */
  Type distance(const Point3& p1) const;

  /**
   * Returns the square of the distance between this point and
   * point p1.
   *
   * @param p1  the other point.
   */
  Type distanceSquared(const Point3& p1) const;

  /**
   * Multiplies each of the x,y,z components of the Point4
   * parameter by 1/w, places the projected values into this
   * point.
   *
   * @param p1 the source Point4.
   */
  void project(const Point4<Type>& p1);

};


template<class Type>
Point3<Type>::Point3()
  : Tuple3<Type>(0.0, 0.0, 0.0) {
}

template<class Type>
Point3<Type>::Point3(Type x, Type y, Type z)
  : Tuple3<Type>(x, y, z) {
}

template<class Type>
Point3<Type>::Point3(const Tuple3<Type>& t1)
  : Tuple3<Type>(t1) {
}

template<class Type>
Type Point3<Type>::distance(const Point3<Type>& p1) const {
  return std::sqrt(distanceSquared(p1));
}

template<class Type>
Type Point3<Type>::distanceSquared(const Point3<Type>& p1) const {
  Type dx = Tuple3<Type>::x - p1.x;
  Type dy = Tuple3<Type>::y - p1.y;
  Type dz = Tuple3<Type>::z - p1.z;
  return (dx*dx + dy*dy + dz*dz);
}

template<class Type>
void Point3<Type>::project(const Point4<Type>& p1) {
  set(p1.getX()/p1.getW(), p1.getY()/p1.getW(), p1.getZ()/p1.getW());
}

}
}


#endif
