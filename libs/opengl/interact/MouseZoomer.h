#ifndef mi_opengl_interact_MouseZoomer_h
#define mi_opengl_interact_MouseZoomer_h

namespace mi {
namespace opengl {
namespace interact {

class PropertyCommand;

class MouseZoomer {

  PropertyCommand* propertyCommand;

  int previousX;

  int previousY;

  float startValue;

  float currentValue;

  float scaling;

public:

  MouseZoomer(PropertyCommand* propertyCommand, float scaling);

  float getScaling();

  void setScaling(float angle);

  float getStartValue();

  float getCurrentValue();

  void beginZoom(int x, int y);

  void zoom(int x, int y);

  void endZoom();

};

}
}
}

#endif

