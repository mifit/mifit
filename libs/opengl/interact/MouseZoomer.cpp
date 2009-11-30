#include <opengl/interact/MouseZoomer.h>

#include <opengl/interact/PropertyCommand.h>

namespace mi
{
namespace opengl
{
namespace interact
{

MouseZoomer::MouseZoomer(PropertyCommand *propertyCommand, float scaling)
{
    this->propertyCommand = propertyCommand;
    this->scaling = scaling;
}

float MouseZoomer::getScaling()
{
    return scaling;
}

void MouseZoomer::setScaling(float value)
{
    scaling = value;
}

float MouseZoomer::getStartValue()
{
    return startValue;
}

float MouseZoomer::getCurrentValue()
{
    return currentValue;
}

void MouseZoomer::beginZoom(int x, int y)
{
    previousX = x;
    previousY = y;
    startValue = propertyCommand->getValue();
}

void MouseZoomer::endZoom()
{
}

void MouseZoomer::zoom(int, int y)
{
    currentValue = startValue;
    currentValue += (previousY - y) * scaling;
    propertyCommand->setValue(currentValue);
}

}
}
}

