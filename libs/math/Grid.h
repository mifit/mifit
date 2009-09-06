#ifndef mi_math_Grid_h
#define mi_math_Grid_h

#include <math/Point3.h>
#include <string>
#include <sstream>

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
class Grid {
protected:

  /**
   *  Stores the values in the grid
   */
  float* data;

  /**
   *  The grid dimensions
   */
  int gridSize0;
  int gridSize1;
  int gridSize2;

  /**
   *  The physical spacing between the integral indecies of the grid
   */
  Type gridSpacing0;
  Type gridSpacing1;
  Type gridSpacing2;

  /**
   *  A label for the grid
   */
  std::string name;

  /**
   *  The physical origin
   */
  Point3<Type> gridOrigin;

  /**
   *  The largest physical point represented
   */
  Point3<Type> gridMaximum;

  /**
   *  Creates a new grid with the specified size, spacing, and origin.
   *
   *  @param size  the xyz sizes (int[3])
   *  @param spacing  the xyz physical spacings (Type[3])
   *  @param origin  the physical origin
   */
  void initialize(int* size, Type* spacing, Point3<Type>& origin);

public:
  /**
   *  Creates a new grid with the specified size, spacing, and origin.
   *
   *  @param size  the xyz sizes (int[3])
   *  @param spacing  the xyz physical spacings (Type[3])
   *  @param origin  the physical origin
   */
  Grid(int* size, Type* spacing, Point3<Type>& origin);

  /**
   *  Creates a new grid with the specified physical boundaries and
   *  spacing.  The size of the grid will be calculated to be
   *  (upperBoundary - lowerBoundary) / spacing + 2.
   *
   *  @param lowerBoundary  the lower physical boundary (Type[3])
   *  @param upperBoundary  the upper physical boundary (Type[3])
   *  @param spacing the xyz physical spacings (Type[3])
   */
  Grid(Type* lowerBoundary, Type* upperBoundary, Type* spacing);

  /**
   * Destroys this grid.
   */
  virtual ~Grid();

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
   *  Returns the number of indecies in each dimension of the grid.
   *
   *  @return the index size (int[3])
   */
  virtual int* getSize();

  /**
   *  Returns the physical origin of the grid space.
   */
  virtual Point3<Type> getOrigin();

  /**
   *  Returns the physical spacing between indecies into the grid space.
   *
   *  @return the physical spacing (Type[3])
   */
  virtual Type* getSpacing();

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
   *  Determines whether the sizes of two grids are equal.
   *
   *  @return true if the sizes of the grids are equal.
   */
  virtual bool isSizeEquals(Grid<Type>& other);

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
   * Stores the current index for efficient communication between methods.
   */
  int currentIndex0;

  /**
   * Stores the current index for efficient communication between methods.
   */
  int currentIndex1;

  /**
   * Stores the current index for efficient communication between methods.
   */
  int currentIndex2;

  /**
   * Sets the current index for efficient communication between methods.
   *
   * @param point the point for whom the index is set.
   */
  virtual void setCurrentIndex(Point3<Type>& point);

  /**
   * Used in linear interpolation.
   */
  Point3<Type> relativePoint;

  /**
   * Used in linear interpolation.
   */
  int cubeOrigin[3];

  /**
   * Used in linear interpolation.
   */
  float fractionalDistance[3];

};

template<class Type>
void Grid<Type>::initialize(int* size, Type* spacing, Point3<Type>& origin) {
  name = "";
  gridSize0 = size[0];
  gridSize1 = size[1];
  gridSize2 = size[2];
  gridSpacing0 = spacing[0];
  gridSpacing1 = spacing[1];
  gridSpacing2 = spacing[2];
  gridOrigin = origin;
  gridMaximum = Point3<Type>(origin.x + gridSize0*gridSpacing0,
                             origin.y + gridSize1*gridSpacing1,
                             origin.z + gridSize2*gridSpacing2);
  data = new float[gridSize0*gridSize1*gridSize2];
  for (int i = 0; i < gridSize0; ++i) {
    for (int j = 0; j < gridSize1; ++j) {
      for (int k = 0; k < gridSize2; ++k) {
        data[i + gridSize0*(j + gridSize1*k)] = 0.0;
      }
    }
  }
}

template<class Type>
Grid<Type>::Grid(int* size, Type* spacing, Point3<Type>& origin) {
  initialize(size, spacing, origin);
}

template<class Type>
Grid<Type>::Grid(Type* lowerBoundary, Type* upperBoundary, Type* spacing) {
  int* size = new int[3];
  size[0] = (int)((upperBoundary[0] - lowerBoundary[0])/spacing[0])+2;
  size[1] = (int)((upperBoundary[1] - lowerBoundary[1])/spacing[1])+2;
  size[2] = (int)((upperBoundary[2] - lowerBoundary[2])/spacing[1])+2;
  Point3<Type> origin(lowerBoundary[0], lowerBoundary[1], lowerBoundary[2]);
  initialize(size, spacing, origin);
  delete[] size;
}

template<class Type>
Grid<Type>::~Grid() {
  delete[] data;
}

template<class Type>
void Grid<Type>::setName(const std::string& newName) {
  name = newName;
}

template<class Type>
Point3<Type> Grid<Type>::getPoint(int* index) {
  return getPoint(index[0], index[1], index[2]);
}

template<class Type>
Point3<Type> Grid<Type>::getPoint(int xIndex, int yIndex, int zIndex) {
  return Point3<Type>(gridOrigin.x + xIndex * gridSpacing0,
                      gridOrigin.y + yIndex * gridSpacing1,
                      gridOrigin.z + zIndex * gridSpacing2);
}

template<class Type>
int* Grid<Type>::getNearestIndex(Point3<Type>& point) {
  int* result = new int[3];
  setCurrentIndex(point);
  result[0] = currentIndex0;
  result[1] = currentIndex1;
  result[2] = currentIndex2;
  return result;
}

template<class Type>
void Grid<Type>::getNearestIndex(Point3<Type>& point, int* index) {
  index[0] = (int) rint((point.x - gridOrigin.x) / gridSpacing0);
  index[1] = (int) rint((point.y - gridOrigin.y) / gridSpacing1);
  index[2] = (int) rint((point.z - gridOrigin.z) / gridSpacing2);
}

template<class Type>
void Grid<Type>::setCurrentIndex(Point3<Type>& point) {
  currentIndex0 = (int) rint((point.x - gridOrigin.x) / gridSpacing0);
  currentIndex1 = (int) rint((point.y - gridOrigin.y) / gridSpacing1);
  currentIndex2 = (int) rint((point.z - gridOrigin.z) / gridSpacing2);
}

template<class Type>
std::string Grid<Type>::getName() {
  return name;
}

template<class Type>
float Grid<Type>::getValue(int* index) {
  return data[index[0] + gridSize0*(index[1] + gridSize1*index[2])];
}

template<class Type>
float Grid<Type>::getValue(int a, int b, int c) {
  return data[a + gridSize0*(b + gridSize1*c)];
}

template<class Type>
void Grid<Type>::setValue(int* index, float value) {
  data[index[0] + gridSize0*(index[1] + gridSize1*index[2])] = value;
}

template<class Type>
void Grid<Type>::setValue(int a, int b, int c, float value) {
  data[a + gridSize0*(b + gridSize1*c)] = value;
}

template<class Type>
int* Grid<Type>::getSize() {
  int* result = new int[3];
  result[0] = gridSize0;
  result[1] = gridSize1;
  result[2] = gridSize2;
  return result;
}

template<class Type>
Point3<Type> Grid<Type>::getOrigin() {
  return gridOrigin;
}

template<class Type>
Type* Grid<Type>::getSpacing() {
  Type* result = new Type[3];
  result[0] = gridSpacing0;
  result[1] = gridSpacing1;
  result[2] = gridSpacing2;
  return result;
}

template<class Type>
bool Grid<Type>::isIndexInGrid(int* index) {
  if ((index[0] >= 0 && index[0] < gridSize0)
      && (index[1] >= 0 && index[1] < gridSize1)
      && (index[2] >= 0 && index[2] < gridSize2)) {
    return true;
  } else {
    return false;
  }
}

template<class Type>
bool Grid<Type>::isIndexInGrid(int a, int b, int c) {
  if ((a >= 0 && a < gridSize0)
      && (b >= 0 && b < gridSize1)
      && (c >= 0 && c < gridSize2)) {
    return true;
  } else {
    return false;
  }
}

template<class Type>
bool Grid<Type>::isPointInGrid(Point3<Type>& point) {
  setCurrentIndex(point);
  return isIndexInGrid(currentIndex0, currentIndex1, currentIndex2);
}

template<class Type>
bool Grid<Type>::isSizeEquals(Grid<Type>& other) {
  int* otherSize = other.getSize();
  if (gridSize0 != otherSize[0]) {
    return false;
  }
  if (gridSize1 != otherSize[1]) {
    return false;
  }
  if (gridSize2 != otherSize[2]) {
    return false;
  }
  return true;
}

template<class Type>
float Grid<Type>::getLinearInterpolatedValue(Point3<Type>& pt) {
  float xYZInterpolation = 0.0f;

  if (isPointInGrid(pt)) {

    // Coordinates of interpolation point relative to map origin
    relativePoint.set(pt);
    relativePoint.subtract(gridOrigin);

    // Integer origin of eight grid points making up cube around interpolation point
    cubeOrigin[0] = (int)(relativePoint.x / gridSpacing0);
    cubeOrigin[1] = (int)(relativePoint.y / gridSpacing1);
    cubeOrigin[2] = (int)(relativePoint.z / gridSpacing2);

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

  } else {
    std::ostringstream message;
    message << "point is not in the grid: (" << pt.getX() << ", "
            << pt.getY() << ", " << pt.getZ() << ")";
    throw IndexOutOfBoundsException(message.str().c_str());
  }

  return xYZInterpolation;
}

template<class Type>
float Grid<Type>::getNearestValue(Point3<Type>& point) {
  setCurrentIndex(point);
  return getValue(currentIndex0, currentIndex1, currentIndex2);
}

template<class Type>
void Grid<Type>::setNearestValue(Point3<Type>& point, float value) {
  setCurrentIndex(point);
  setValue(currentIndex0, currentIndex1, currentIndex2, value);
}

}
}

#endif

