#include <opengl/interact/MouseTranslator.h>

#include <opengl/Axis.h>
#include <opengl/Viewpoint.h>

using namespace mi::math;
using namespace mi::opengl;

namespace mi
{
namespace opengl
{
namespace interact
{

MouseTranslator::MouseTranslator(Viewpoint *viewpoint, Viewpoint *relativeToViewpoint, float scaling)
    : startPosition(NULL),
      currentPosition(new Vector3<float>()),
      axes(new Axis),
      currentAxes(NULL)
{

    this->viewpoint = viewpoint;
    this->relativeToViewpoint = relativeToViewpoint;
    this->scaling = scaling;
}

MouseTranslator::~MouseTranslator()
{
    delete currentPosition;
    delete axes;
}

float MouseTranslator::getScaling()
{
    return scaling;
}

void MouseTranslator::setScaling(float angle)
{
    scaling = angle;
}

Vector3<float>*MouseTranslator::getStartPosition()
{
    return startPosition;
}

Vector3<float>*MouseTranslator::getCurrentPosition()
{
    return currentPosition;
}

void MouseTranslator::setAxes(const Axis &axes)
{
    this->axes->set(axes);
}

void MouseTranslator::beginTranslate(int x, int y)
{
    previousX = x;
    previousY = y;
    startPosition = new Vector3<float>(viewpoint->getPosition());
    Quaternion<float> orientation(relativeToViewpoint->getRotation());
    currentAxes = new Axis();
    currentAxes->set(*axes);
    currentAxes->rotate(orientation);
}

void MouseTranslator::endTranslate()
{
    delete startPosition;
    startPosition = NULL;
    delete currentAxes;
    currentAxes = NULL;
}

void MouseTranslator::translate(int x, int y)
{
    if (startPosition == NULL)
    {
        return;
    }
    currentPosition->set(*startPosition);
    int deltaX = previousX - x;
    int deltaY = y - previousY;
    Vector3<float> xChange(currentAxes->x);
    xChange.scale(deltaX * scaling);
    Vector3<float> yChange(currentAxes->y);
    yChange.scale(deltaY * scaling);
    Vector3<float> translation;
    translation.add(xChange, yChange);
    currentPosition->add(translation);
    viewpoint->setPosition(*currentPosition);
}

}
}
}
