#ifndef mi_opengl_interact_MouseArcBall_h
#define mi_opengl_interact_MouseArcBall_h


#include <math/Point3.h>
#include <math/Quaternion.h>
#include <math/Vector3.h>

namespace mi {
namespace opengl {

class Viewpoint;

namespace interact {

/**
 * MouseArcBall is a mouse interaction object that lets users control the
 * rotation of an object via a mouse using quaterion math.
 */
class MouseArcBall {

  int controlSizeX;

  int controlSizeY;

  float controlAspectX;

  float controlAspectY;

  mi::opengl::Viewpoint* relativeToViewpoint;


  /**
   * Converts window coordinates to sphere coordinates.
   */
  void windowToSphereCoordinates(mi::math::Vector3<float>& vect);

  void setCurrentVector(int x, int y);

protected:

  mi::math::Point3<float> center;

  float radius;

  bool invertRotation;

  mi::math::Quaternion<float> historicalQuat;

  mi::opengl::Viewpoint* viewpoint;

  mi::math::Quaternion<float> currentQuat;

  mi::math::Vector3<float> currentVector;

  mi::math::Vector3<float>* startVector;

  mi::math::Quaternion<float> relativeStartQuat;

  virtual void initializeRotation();

  virtual void applyRotation();

public:

  MouseArcBall(mi::opengl::Viewpoint* viewpoint, mi::math::Point3<float>& center, float radius);

  MouseArcBall(mi::opengl::Viewpoint* viewpoint);

  virtual ~MouseArcBall() {
  }

  void setRelativeToViewpoint(mi::opengl::Viewpoint* relativeToViewpoint);

  void setPlace(const mi::math::Point3<float>& center, float radius);

  mi::math::Point3<float>& getCenter();

  float getRadius();

  void setRadius(float radius);

  bool isInvertRotation();

  void setInvertRotation(bool invertRotation);

  mi::math::Vector3<float>* getStartVector();
  mi::math::Vector3<float>* getCurrentVector();

  void beginRotate(int x, int y);

  void rotate(int x, int y);

  void endRotate();

  void setControlSize(int x, int y);

  void getControlSize(int& x, int& y);

};
}
}
}

#endif
