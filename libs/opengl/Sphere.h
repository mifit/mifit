#ifndef mi_opengl_Sphere_h
#define mi_opengl_Sphere_h

namespace mi {
namespace opengl {

class Sphere {

  float* vertexBuffer;

  float* normalBuffer;

  float* textureBuffer;

  int vertexCount;

public:

  Sphere(float radius, int slices, int stacks);

  void render();

};

}
}

#endif
