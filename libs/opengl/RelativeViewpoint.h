#ifndef mi_opengl_RelativeViewpoint_h
#define mi_opengl_RelativeViewpoint_h

#include <math/Vector3.h>
#include <math/Quaternion.h>
#include <opengl/Viewpoint.h>

namespace mi {
namespace opengl {

class RelativeViewpoint : public Viewpoint {

protected:

  Viewpoint* reference;

public:

  RelativeViewpoint(Viewpoint* reference);

  virtual mi::math::Quaternion<float> getRotation();

  virtual void setRotation(mi::math::Quaternion<float>& rotation);

};
}
}

#endif
