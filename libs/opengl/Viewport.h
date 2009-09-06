#ifndef mi_opengl_Viewport_h
#define mi_opengl_Viewport_h

namespace mi {
namespace opengl {

class Viewport {

  int width;

  int height;

  int x;

  int y;

public:

  Viewport();

  void set(int x, int y, int width, int height);

  int getX();

  int getY();

  int getWidth();

  int getHeight();

  int* toArray();

  void render();

};

}
}

#endif

