#ifndef mi_opengl_interact_PropertyCommand_h
#define mi_opengl_interact_PropertyCommand_h

namespace mi {
namespace opengl {
namespace interact {

class PropertyCommand {

public:
  virtual ~PropertyCommand() {
  }

  virtual float getValue() = 0;

  virtual void setValue(float value) = 0;

};

}
}
}

#endif
