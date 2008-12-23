#include <opengl/ViewportRelativeViewpoint.h>

#include <opengl/Frustum.h>
#include <opengl/Viewport.h>

using namespace mi::math;

namespace mi {
namespace opengl {

ViewportRelativeViewpoint::ViewportRelativeViewpoint(Viewpoint* reference) : RelativeViewpoint(reference) {
}

void ViewportRelativeViewpoint::setPosition(Viewport* viewport, Frustum* frustum, float x, float y) {
  float heightAtTarget = frustum->getFrustumHeight(frustum->getFocalLength());
  float glUnitsPerPixel = heightAtTarget / (float) viewport->getHeight();
  Vector3<float> position(reference->getTarget(frustum->getFocalLength()));
  Vector3<float> rightVector(reference->getRightVector());
  Vector3<float> upVector(reference->getUpVector());
  rightVector.scale(x * glUnitsPerPixel);
  upVector.scale(y * glUnitsPerPixel);
  position.add(rightVector);
  position.add(upVector);
  Viewpoint::setPosition(position);
}

}
}
