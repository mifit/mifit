#ifndef mi_opengl_Light_h
#define mi_opengl_Light_h

#include <math/Vector4.h>
#include <opengl/OpenGL.h>


namespace mi {
namespace opengl {

class Light {

  int lightId;

  float* ambientColor;

  float* diffuseColor;

  float* specularColor;

  float* position;

public:

  Light(int lightId = GL_LIGHT1, float* ambientColor = NULL);

  void setPosition(mi::math::Vector4<float>& position);

  void setAmbientColor(float color0, float color1, float color2, float color3);

  void setDiffuseColor(float color0, float color1, float color2, float color3);

  void setSpecularColor(float color0, float color1, float color2, float color3);

  void render();

  mi::math::Vector4<float> getPosition();

  float* getAmbientColor();

  float* getDiffuseColor();

  float* getSpecularColor();

  void glEnable();

  void glDisable();
};

}
}

#endif
