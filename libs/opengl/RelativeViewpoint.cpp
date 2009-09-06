#include <opengl/RelativeViewpoint.h>

using namespace mi::math;

namespace mi {
namespace opengl {

RelativeViewpoint::RelativeViewpoint(Viewpoint* reference) : Viewpoint(), reference(reference) {
}

Quaternion<float> RelativeViewpoint::getRotation() {
  Quaternion<float> q(reference->getRotation());
  q.multiply(Viewpoint::getRotation());
  return q;
}

void RelativeViewpoint::setRotation(Quaternion<float>& rotation) {
  Quaternion<float> q(reference->getRotation());
  q.inverse();
  q.multiply(rotation);
  Viewpoint::setRotation(q);
}

}
}
