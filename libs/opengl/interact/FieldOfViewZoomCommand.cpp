#include <opengl/interact/FieldOfViewZoomCommand.h>

#include <opengl/Frustum.h>

namespace mi
{
namespace opengl
{
namespace interact
{

FieldOfViewZoomCommand::FieldOfViewZoomCommand(Frustum *frustum)
    : frustum(frustum)
{
}

float FieldOfViewZoomCommand::getValue()
{
    return frustum->getFieldOfView();
}

void FieldOfViewZoomCommand::setValue(float value)
{
    frustum->setFieldOfView(value);
}

}
}
}
