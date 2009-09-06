#ifndef mi_opengl_interact_MousePicker_h
#define mi_opengl_interact_MousePicker_h

#include <vector>
#include <opengl/OpenGL.h>

namespace mi {
namespace opengl {

class Frustum;
class Renderable;

namespace interact {

class MousePicker {

  int selectBufferLength;
  GLuint* selectBuffer;

public:

  MousePicker();
  ~MousePicker();

  std::vector<GLuint> pick(int x, int y, Frustum* frustum, Renderable* scene);
};
}
}
}

#endif
