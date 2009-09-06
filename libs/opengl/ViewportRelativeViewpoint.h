#ifndef mi_opengl_ViewportRelativeViewpoint_h
#define mi_opengl_ViewportRelativeViewpoint_h

#include <opengl/RelativeViewpoint.h>

namespace mi {
namespace opengl {

class Frustum;
class Viewport;

class ViewportRelativeViewpoint : public RelativeViewpoint {

public:

  ViewportRelativeViewpoint(Viewpoint* reference);

  void setPosition(Viewport* viewport, Frustum* frustum, float x, float y);

};
}
}

#endif
