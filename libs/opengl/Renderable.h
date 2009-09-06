#ifndef mi_opengl_Renderable_h
#define mi_opengl_Renderable_h

namespace mi {
namespace opengl {

class Renderable {
public:
  virtual ~Renderable() {
  }

  virtual void render() = 0;

};
}
}

#endif
