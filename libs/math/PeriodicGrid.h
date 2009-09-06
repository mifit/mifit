#ifndef mi_math_PeriodicGrid_h
#define mi_math_PeriodicGrid_h

#include <math/Point3.h>
#include <math/Grid.h>
#include <string>

namespace mi {
namespace math {

/**
 *  Represents a evenly distributed set of points in three-dimensional space.
 *
 *  <p>Arrays of integers are used as indecies into the grid. However,
 *  the length of the arrays are not checked. If an array is not of
 *  the appropriate length, an ArrayIndexOutOfBoundsException may
 *  occur. The appropriate lengths of arrays are specified for each
 *  method.
 */
template<class Type>
class PeriodicGrid : public Grid<Type> {
public:
  /**
   *  Creates a new grid with the specified size, spacing, and origin.
   *
   *  @param size  the xyz sizes (int[3])
   *  @param spacing  the xyz physical spacings (Type[3])
   *  @param origin  the physical origin
   */
  PeriodicGrid(int* size, Type* spacing, Point3<Type>& origin);

  /**
   *  Creates a new grid with the specified physical boundaries and
   *  spacing.  The size of the grid will be calculated to be
   *  (upperBoundary - lowerBoundary) / spacing + 2.
   *
   *  @param lowerBoundary  the lower physical boundary (Type[3])
   *  @param upperBoundary  the upper physical boundary (Type[3])
   *  @param spacing the xyz physical spacings (Type[3])
   */
  PeriodicGrid(Type* lowerBoundary, Type* upperBoundary, Type* spacing);

  /**
   * Destroys this grid.
   */
  virtual ~PeriodicGrid();

  /**
   *  Sets the label for the grid.
   */
  virtual void setName(const std::string& newName);

  /**
   *  Returns the physical point represented by the given index.  If
   *  the index is within the grid, the point returned will lie
   *  within the space represented by the grid.
   *
   *  @param index  an index into the grid (int[3])
   */
  virtual Point3<Type> getPoint(int* index);

  /**
   *  Returns the physical point represented by the given index.  If
   *  the index is within the grid, the point returned will lie
   *  within the space represented by the grid.
   *
   *  @param index an index into the grid (int[3])
   */
  virtual Point3<Type> getPoint(int xIndex, int yIndex, int zIndex);

  /**
   *  Returns the estimated index nearest to the given point. The
   *  index returned may be outside the grid. It is left to the user
   *  to verify the index is in the grid with the
   *  <code>isIndexInGrid</code> method.
   *
   *  @return the estimated index (int[3])
   */
  virtual int* getNearestIndex(Point3<Type>& point);

  virtual void getNearestIndex(Point3<Type>& point, int* index);

  /**
   *  Returns the label for this grid.
   */
  virtual std::string getName();

  /**
   *  Returns the value of the grid at the given index.
   *
   *  @throws ArrayIndexOutOfBoundsException if the index derived
   *  from the point lies outside the grid.
   */
  virtual float getValue(int* index);

  /**
   *  Returns the value of the grid at the index represented by a b c.
   *
   *  @param a the index in the first dimension
   *  @param b the index in the second dimension
   *  @param c the index in the third dimension
   *  @throws ArrayIndexOutOfBoundsException if the index derived
   *  from the point lies outside the grid.
   */
  virtual float getValue(int a, int b, int c);

  /**
   *  Sets the value of the grid at the given index.
   *
   *  @throws ArrayIndexOutOfBoundsException if the index derived
   *  from the point lies outside the grid.
   */
  virtual void setValue(int* index, float value);

  /**
   *  Sets the value of the grid at the given index.
   *
   *  @param a the index in the first dimension
   *  @param b the index in the second dimension
   *  @param c the index in the third dimension
   *  @throws ArrayIndexOutOfBoundsException if the index derived
   *  from the point lies outside the grid.
   */
  virtual void setValue(int a, int b, int c, float value);

  /**
   * Returns whether the given index is in the grid.
   *
   * @return true if the index is contained in the grid, otherwise false.
   */
  virtual bool isIndexInGrid(int* index);

  /**
   * Returns whether the given index is in the grid.
   *
   * @return true if the index is contained in the grid, otherwise false.
   */
  virtual bool isIndexInGrid(int a, int b, int c);

  /**
   * Returns whether the given point is in the grid.
   *
   * @return true if the point is contained in the grid, otherwise false.
   */
  virtual bool isPointInGrid(Point3<Type>& point);

  /**
   *  Linear interpolation of the value of a point which does not lie on a grid point.
   *
   *  Determines the fractional distance is each direction(i.e. x,
   *  y, z). Then determines x value using fractionsl distance, and
   *  determines y and z values similarly. Interpolates xy plane
   *  value using x and y values. Then interpolates xyz value using
   *  xy plane and z value.
   *
   *  @param pt a physical point within the grid space
   *  @return value at the point
   *  @throws ArrayIndexOutOfBoundsException if the point is not in the grid
   */
  virtual float getLinearInterpolatedValue(Point3<Type>& pt);

  /**
   *  Returns the value of the grid nearest the given point.
   *
   *  @param point a physical point within the grid space
   *  @throws ArrayIndexOutOfBoundsException if the index derived
   *  from the point lies outside the grid.
   */
  virtual float getNearestValue(Point3<Type>& point);

  /**
   *  Sets the value of the grid nearest the given point.
   *
   *  @param point a physical point within the grid space
   *  @param value  the value of the point
   *  @ throws ArrayIndexOutOfBoundsException if the index derived
   *  from the point lies outside the grid.
   */
  virtual void setNearestValue(Point3<Type>& point, float value);

protected:
  /**
   * Sets the current index for efficient communication between methods.
   *
   * @param point the point for whom the index is set.
   */
  virtual void setCurrentIndex(Point3<Type>& point);

  /**
   * Applies the boundary conditions to the given index. If an index
   * is above or below its respective grid dimension, the index has the
   * grid size for that dimension added or subtracted from the value until
   * the index is within the grid.
   *
   *  @param index an index into the grid (int[3])
   *  @param index0 an index into the grid
   *  @param index1 an index into the grid
   *  @param index2 an index into the grid
   */
  virtual void applyBoundaryConditions(int& index0, int& index1, int& index2);

  /**
   * Applies the boundary conditions to the given index. If an index
   * is above or below its respective grid dimension, the index has the
   * grid size for that dimension added or subtracted from the value until
   * the index is within the grid.
   *
   *  @param index an index into the grid (int[3])
   */
  virtual void applyBoundaryConditions(int* index);

  /**
   * Applies the boundary conditions to the given point. If an point
   * is above or below its respective grid area, the index has the
   * grid size for that axis added or subtracted from the value until
   * the point is within the grid.
   *
   *  @param point the point to adjust to within grid area.
   */
  virtual void applyBoundaryConditions(Point3<Type>& point);
};

template<class Type>
PeriodicGrid<Type>::PeriodicGrid(int* size, Type* spacing, Point3<Type>& origin)
  : Grid<Type>(size, spacing, origin) {
}

template<class Type>
PeriodicGrid<Type>::PeriodicGrid(Type* lowerBoundary, Type* upperBoundary, Type* spacing)
  : Grid<Type>(lowerBoundary, upperBoundary, spacing) {
}

template<class Type>
PeriodicGrid<Type>::~PeriodicGrid() {
}

template<class Type>
void PeriodicGrid<Type>::setName(const std::string& newName) {
  Grid<Type>::name = newName;
}

template<class Type>
Point3<Type> PeriodicGrid<Type>::getPoint(int* index) {
  return getPoint(index[0], index[1], index[2]);
}

template<class Type>
Point3<Type> PeriodicGrid<Type>::getPoint(int xIndex, int yIndex, int zIndex) {
  return Grid<Type>::getPoint(xIndex, yIndex, zIndex);
}

template<class Type>
int* PeriodicGrid<Type>::getNearestIndex(Point3<Type>& point) {
  int* result = new int[3];
  setCurrentIndex(point);
  result[0] = Grid<Type>::currentIndex0;
  result[1] = Grid<Type>::currentIndex1;
  result[2] = Grid<Type>::currentIndex2;
  return result;
}

template<class Type>
void PeriodicGrid<Type>::getNearestIndex(Point3<Type>& point, int* index) {
  Grid<Type>::getNearestIndex(point, index);
}

template<class Type>
void PeriodicGrid<Type>::setCurrentIndex(Point3<Type>& point) {
  Grid<Type>::setCurrentIndex(point);
}

template<class Type>
std::string PeriodicGrid<Type>::getName() {
  return Grid<Type>::name;
}

template<class Type>
float PeriodicGrid<Type>::getValue(int* index) {
  applyBoundaryConditions(index);
  return Grid<Type>::data[index[0] + Grid<Type>::gridSize0*(index[1] + Grid<Type>::gridSize1*index[2])];
}

template<class Type>
float PeriodicGrid<Type>::getValue(int a, int b, int c) {
  applyBoundaryConditions(a, b, c);
  return Grid<Type>::data[a + Grid < Type > ::gridSize0*(b + Grid < Type > ::gridSize1*c)];
}

template<class Type>
void PeriodicGrid<Type>::setValue(int* index, float value) {
  applyBoundaryConditions(index);
  Grid<Type>::setValue(index, value);
}

template<class Type>
void PeriodicGrid<Type>::setValue(int a, int b, int c, float value) {
  applyBoundaryConditions(a, b, c);
  Grid<Type>::setValue(a, b, c, value);
}

template<class Type>
bool PeriodicGrid<Type>::isIndexInGrid(int* index) {
  return true;
}

template<class Type>
bool PeriodicGrid<Type>::isIndexInGrid(int a, int b, int c) {
  return true;
}

template<class Type>
bool PeriodicGrid<Type>::isPointInGrid(Point3<Type>& point) {
  return true;
}

template<class Type>
float PeriodicGrid<Type>::getLinearInterpolatedValue(Point3<Type>& pt) {
  applyBoundaryConditions(pt);
  return Grid<Type>::getLinearInterpolatedValue(pt);
  /*
     float xYZInterpolation = 0.0f;

     // Coordinates of interpolation point relative to map origin
     relativePoint.set(pt);
     relativePoint.subtract(gridOrigin);

     // Integer origin of eight grid points making up cube around interpolation point
     cubeOrigin[0] = (int)(relativePoint.x / gridSpacing0);
     cubeOrigin[1] = (int)(relativePoint.y / gridSpacing1);
     cubeOrigin[2] = (int)(relativePoint.z / gridSpacing2);
     cout << cubeOrigin[0] << ' ' << cubeOrigin[1] << ' ' << cubeOrigin[2] << endl;

     // Fractional distance in each direction
     fractionalDistance[0] = (float)(relativePoint.x/gridSpacing0 - (double)(cubeOrigin[0]));
     fractionalDistance[1] = (float)(relativePoint.y/gridSpacing1 - (double)(cubeOrigin[1]));
     fractionalDistance[2] = (float)(relativePoint.z/gridSpacing2 - (double)(cubeOrigin[2]));

     float temp = getValue(cubeOrigin);
     float bottomXInterpolation = temp + fractionalDistance[0] * (getValue(cubeOrigin[0] + 1, cubeOrigin[1], cubeOrigin[2]) - temp);
     temp = getValue(cubeOrigin[0], cubeOrigin[1] + 1, cubeOrigin[2]);
     float topXInterpolation = temp + fractionalDistance[0] * (getValue(cubeOrigin[0] + 1, cubeOrigin[1] + 1, cubeOrigin[2]) - temp);
     float xInterpolation = bottomXInterpolation + fractionalDistance[1] * (topXInterpolation - bottomXInterpolation);

     temp = getValue(cubeOrigin[0], cubeOrigin[1], cubeOrigin[2] + 1);
     float bottomXYInterpolation = temp + fractionalDistance[0] * (getValue(cubeOrigin[0] + 1, cubeOrigin[1], cubeOrigin[2] + 1) - temp);

     temp = getValue(cubeOrigin[0], cubeOrigin[1] + 1, cubeOrigin[2] + 1);
     float topXYInterpolation = temp + fractionalDistance[0] * (getValue(cubeOrigin[0] + 1, cubeOrigin[1] + 1, cubeOrigin[2] + 1) - temp);
     float xYInterpolation = bottomXYInterpolation + fractionalDistance[1] * (topXYInterpolation - bottomXYInterpolation);

     xYZInterpolation = xInterpolation + fractionalDistance[2] * (xYInterpolation - xInterpolation);

     return xYZInterpolation;
   */
}

template<class Type>
float PeriodicGrid<Type>::getNearestValue(Point3<Type>& point) {
  setCurrentIndex(point);
  return getValue(Grid<Type>::currentIndex0, Grid<Type>::currentIndex1, Grid<Type>::currentIndex2);
}

template<class Type>
void PeriodicGrid<Type>::setNearestValue(Point3<Type>& point, float value) {
  setCurrentIndex(point);
  setValue(Grid<Type>::currentIndex0, Grid<Type>::currentIndex1, Grid<Type>::currentIndex2, value);
}

template<class Type>
void PeriodicGrid<Type>::applyBoundaryConditions(int* index) {
  applyBoundaryConditions(index[0], index[1], index[2]);
}

template<class Type>
void PeriodicGrid<Type>::applyBoundaryConditions(int& index0, int& index1, int& index2) {
  while (index0 >= Grid<Type>::gridSize0) {
    index0 -= Grid<Type>::gridSize0;
  }
  while (index0 < 0) {
    index0 += Grid<Type>::gridSize0;
  }
  while (index1 >= Grid<Type>::gridSize1) {
    index1 -= Grid<Type>::gridSize1;
  }
  while (index1 < 0) {
    index1 += Grid<Type>::gridSize1;
  }
  while (index2 >= Grid<Type>::gridSize2) {
    index2 -= Grid<Type>::gridSize2;
  }
  while (index2 < 0) {
    index2 += Grid<Type>::gridSize2;
  }
}

template<class Type>
void PeriodicGrid<Type>::applyBoundaryConditions(Point3<Type>& point) {
  Type size0 = Grid<Type>::gridSize0 * Grid<Type>::gridSpacing0;
  Type size1 = Grid<Type>::gridSize1 * Grid<Type>::gridSpacing1;
  Type size2 = Grid<Type>::gridSize2 * Grid<Type>::gridSpacing2;

  while (point.getX() >= Grid<Type>::gridOrigin.getX() + size0) {
    point.setX(point.getX() - size0);
  }
  while (point.getX() < Grid<Type>::gridOrigin.getX()) {
    point.setX(point.getX() + size0);
  }
  while (point.getY() >= Grid<Type>::gridOrigin.getY() + size1) {
    point.setY(point.getY() - size1);
  }
  while (point.getY() < Grid<Type>::gridOrigin.getY()) {
    point.setY(point.getY() + size1);
  }
  while (point.getZ() >= Grid<Type>::gridOrigin.getZ() + size2) {
    point.setZ(point.getZ() - size2);
  }
  while (point.getZ() < Grid<Type>::gridOrigin.getZ()) {
    point.setZ(point.getZ() + size2);
  }
}

}
}

#endif

