#ifndef mi_opengl_interact_FieldOfViewZoomCommand_h
#define mi_opengl_interact_FieldOfViewZoomCommand_h

#include <opengl/interact/PropertyCommand.h>

namespace mi {
namespace opengl {

class Frustum;

namespace interact {

class FieldOfViewZoomCommand : public PropertyCommand {

  Frustum* frustum;

public:

  FieldOfViewZoomCommand(Frustum* frustum);

  virtual float getValue();

  virtual void setValue(float value);
};
}
}
}

#endif
