#ifndef mi_opengl_Camera_h
#define mi_opengl_Camera_h

#include <math/Vector3.h>
#include <math/Quaternion.h>
#include <opengl/Viewpoint.h>

namespace mi {
namespace opengl {

class Camera : public Viewpoint {

public:

  Camera();

  mi::math::Vector3<float> getEye();

  void setEye(const mi::math::Vector3<float>& eye);

  void render();

  void lookAt(const mi::math::Vector3<float>& target);

};

}
}

#endif
