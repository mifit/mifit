#ifndef mi_opengl_interact_TargetFeedback_h
#define mi_opengl_interact_TargetFeedback_h

#include <math/Vector3.h>
#include <math/Point3.h>

namespace mi {
namespace opengl {
namespace interact {

class TargetFeedback {

  mi::math::Vector3<float> target;

  mi::math::Point3<float> color;

  float length;

public:

  TargetFeedback(float* color);

  void setLength(float length);

  void setTarget(const mi::math::Vector3<float>& target);

  void render();

};

}
}
}

#endif
