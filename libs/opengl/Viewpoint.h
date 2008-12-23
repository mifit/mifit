#ifndef mi_opengl_Viewpoint_h
#define mi_opengl_Viewpoint_h

#include <math/Vector3.h>
#include <math/Quaternion.h>

namespace mi {
namespace opengl {

class Viewpoint {

  mi::math::Vector3<float> position;
  mi::math::Quaternion<float> rotation;

protected:

  void transformVector(mi::math::Vector3<float>& vector);

  mi::math::Vector3<float> getTargetVector(float focalLength);

public:

  Viewpoint();
  virtual ~Viewpoint() {
  }

  void set(const Viewpoint& viewpoint);

  mi::math::Vector3<float> getPosition() const;

  void setPosition(const mi::math::Vector3<float>& position);

  void translate(mi::math::Vector3<float>& translation);

  virtual mi::math::Quaternion<float> getRotation() const;

  virtual void setRotation(const mi::math::Quaternion<float>& rotation);

  void rotate(mi::math::Quaternion<float>& rotation);

  void render();

  mi::math::Vector3<float> getTarget(float focalLength);

  mi::math::Vector3<float> getViewVector();

  mi::math::Vector3<float> getRightVector();

  mi::math::Vector3<float> getUpVector();

};
}
}

#endif
